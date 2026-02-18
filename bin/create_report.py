#!/usr/bin/env python3

import json
import argparse
import os
import pandas as pd
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer, PageBreak
from reportlab.platypus.flowables import AnchorFlowable
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_RIGHT
from reportlab.lib import colors
from reportlab.lib.units import inch
from datetime import datetime
from map_taxid import load_taxdump, load_names, load_merged


def get_high_ani_matches(member):
    """Return the pre-computed high-ANI match list embedded in the member dict.

    match_paths.py populates each member's ``high_ani_matches`` field as a list
    of dicts::

        [{"key": "<taxid>", "ani_pct": <float 0-100>}, ...]

    sorted descending by ani_pct.  Returns an empty list when the field is
    absent (e.g. ANI matrix was not enabled during match_paths run).
    """
    return member.get('high_ani_matches') or []


def check_if_any_high_ani_in_dataset(strains):
    """Return True if any strain in the dataset has at least one high-ANI match."""
    for strain in strains:
        if get_high_ani_matches(strain):
            return True
    return False


def check_if_k2_reads_present(strains):
    for strain in strains:
        if strain.get('k2_reads') is not None and strain.get('k2_reads', 0) > 0:
            return True
    return False


def should_include_strain(strain, args):
    category = str(strain.get('microbial_category', 'Unknown'))
    if 'Primary' in category:
        return True
    if "Opportunistic" in category and args.show_opportunistic:
        return True
    if 'Potential' in category and args.show_potentials:
        return True
    if 'Commensal' in category and args.show_commensals:
        return True
    if category == 'Unknown' and args.show_unidentified:
        return True
    return False


def passes_confidence_threshold(strain, threshold):
    if threshold > 1.0:
        threshold = threshold / 100.0
    tass_score = strain.get('tass_score', 0)
    return tass_score >= threshold


def load_json_samples(input_files):
    all_sample_data = []
    for input_file in input_files:
        with open(input_file, 'r') as f:
            data = json.load(f)
            all_sample_data.extend(data)
    return all_sample_data


def organize_data_by_sample(sample_data):
    samples_dict = {}
    for species_group in sample_data:
        sample_name = species_group.get('sample_name', 'Unknown Sample')
        if sample_name not in samples_dict:
            samples_dict[sample_name] = []
        samples_dict[sample_name].append(species_group)
    return samples_dict


def sanitize_bookmark_name(name):
    import re
    sanitized = re.sub(r'[^\w\-]', '_', str(name))
    sanitized = re.sub(r'_+', '_', sanitized)
    sanitized = sanitized.strip('_')
    return sanitized


def collect_all_bookmarks(samples_dict, low_confidence_strains):
    bookmarks = set()
    bookmarks.add('color_key')
    bookmarks.add('column_explanations')
    if low_confidence_strains:
        bookmarks.add('low_confidence')
    bookmarks.add('additional_info')
    # Strain-detail appendix anchors (always registered so links validate)
    bookmarks.add('strain_detail_table')
    for sample_name, species_groups in samples_dict.items():
        sample_bookmark = f"sample_{sanitize_bookmark_name(sample_name)}"
        bookmarks.add(sample_bookmark)
        strain_sample_bm = f"strain_sample_{sanitize_bookmark_name(sample_name)}"
        bookmarks.add(strain_sample_bm)
        for species_group in species_groups:
            species_key = species_group.get('toplevelkey', species_group.get('key', 'unknown'))
            species_bookmark = f"species_{sanitize_bookmark_name(sample_name)}_{species_key}"
            bookmarks.add(species_bookmark)
            strain_row_bm = f"strain_row_{sanitize_bookmark_name(sample_name)}_{species_key}"
            bookmarks.add(strain_row_bm)
    return bookmarks


def create_safe_link(text, bookmark, valid_bookmarks, color="blue"):
    if bookmark in valid_bookmarks:
        return f'<link href="#{bookmark}" color="{color}">{text}</link>'
    else:
        print(f"Warning: Bookmark '{bookmark}' does not exist, skipping link")
        return text


def build_taxid_to_bookmark_map(samples_dict):
    taxid_map = {}
    for sample_name, species_groups in samples_dict.items():
        sample_bookmark_part = sanitize_bookmark_name(sample_name)
        for species_group in species_groups:
            species_key = species_group.get('toplevelkey', species_group.get('key', 'unknown'))
            species_bookmark = f"species_{sample_bookmark_part}_{species_key}"
            taxid_map[str(species_key)] = species_bookmark
            for strain in species_group.get('members', []):
                strain_key = str(strain.get('key', ''))
                if strain_key:
                    taxid_map[strain_key] = species_bookmark
    return taxid_map


def get_sample_stats(species_groups):
    total_alignments = sum(sg.get('numreads', 0) for sg in species_groups)
    primary_pathogen_count = 0
    for sg in species_groups:
        for member in sg.get('members', []):
            if 'Primary' in str(member.get('microbial_category', '')):
                primary_pathogen_count += 1
    return total_alignments, primary_pathogen_count


def has_min_reads(strain, min_reads=1):
    try:
        reads = float(strain.get("numreads", 0) or 0)
    except Exception:
        reads = 0
    return reads >= min_reads


def get_species_group_stats(species_group):
    members = species_group.get('members', [])
    primary_count = sum(1 for m in members if 'Primary' in str(m.get('microbial_category', '')))
    group_tass = species_group.get('tass_score', 0)
    tass_scores = [m.get('tass_score', 0) for m in members]
    max_member_tass = max(tass_scores) if tass_scores else group_tass
    return group_tass, primary_count, max_member_tass


def group_members_by_subkey(members):
    """
    Group member strains by their 'subkey' field.
    Falls back to the strain's own 'key' if 'subkey' is absent.
    Returns: {subkey: [strain, ...sorted by TASS desc]}
    """
    subkey_groups = {}
    for strain in members:
        sk = str(strain.get('subkey', strain.get('key', 'unknown')))
        subkey_groups.setdefault(sk, []).append(strain)
    for sk in subkey_groups:
        subkey_groups[sk].sort(key=lambda s: s.get('tass_score', 0), reverse=True)
    return subkey_groups


def get_category_color(microbial_category, ann_class, alpha=1.0):
    val = str(microbial_category)
    derived = str(ann_class)
    if "Primary" in val and derived == "Direct":
        base_color = colors.lightcoral
    elif "Primary" in val:
        base_color = colors.HexColor('#fab462')
    elif "Commensal" in val:
        base_color = colors.lightgreen
    elif "Opportunistic" in val:
        base_color = colors.HexColor('#ffe6a8')
    elif "Potential" in val:
        base_color = colors.lightblue
    else:
        base_color = colors.white
    if alpha < 1.0:
        return colors.Color(base_color.red, base_color.green, base_color.blue, alpha=alpha)
    return base_color

def _valid_num(x):
    try:
        if x is None:
            return None
        v = float(x)
        if v == 0:
            return None
        return v
    except Exception:
        return None


def strip_species_prefix(strain_name, species_name):
    """
    Remove the leading species name from a strain name so only the
    distinguishing part is shown.
    E.g. strip_species_prefix("Escherichia coli ETEC", "Escherichia coli")
         -> "ETEC"
    Returns the original name if the result would be empty or the prefix
    does not match.
    """
    import re
    if not species_name or not strain_name:
        return strain_name
    # Escape special regex chars in the species name and allow variable whitespace
    pattern = r'^\s*' + re.escape(species_name) + r'\s*'
    result = re.sub(pattern, '', strain_name, count=1, flags=re.IGNORECASE).strip()
    # Strip leading punctuation/separators left behind (e.g. "-", "_", "/")
    result = re.sub(r'^[\-_/,;:\s]+', '', result).strip()
    return result if result else strain_name


def create_combined_sample_table(all_strains, species_group_map, small_style,
                                  show_ani_column, show_k2_column,
                                  taxid_to_bookmark, valid_bookmarks,
                                  sample_total_reads=0, sample_name=None,
                                  available_width=None, use_subkey=True,
                                  show_strains_table=True):
    """
    Create a single table combining all strains from all species groups.

    When use_subkey=True (default):
      - Members sharing the same 'subkey' are collapsed into one row.
      - Left columns show the best-TASS member's stats (TASS, K2, Reads, RPM, Coverage).
      - Right 3 columns contain a nested mini-table listing every individual
        strain with its name, TASS score, and read count.

    When use_subkey=False:
      - Original flat per-strain layout (one row per strain, original behaviour).
    """
    # ── Styles ────────────────────────────────────────────────────────────────
    strain_name_style = ParagraphStyle(
        'StrainName', parent=small_style, fontSize=10, leading=10,
        wordWrap='CJK')
    data_style = ParagraphStyle(
        'DataStyle', parent=small_style, fontSize=8, leading=9,
        wordWrap='CJK')
    ani_style = ParagraphStyle(
        'ANIStyle', parent=small_style, fontSize=6, leading=8,
        wordWrap='CJK')
    group_header_style = ParagraphStyle(
        'GroupHeader', parent=small_style, fontSize=10, leading=11,
        fontName='Helvetica-Bold', wordWrap='CJK')
    group_strain_summary_style = ParagraphStyle(
        'GroupStrainSummary', parent=small_style, fontSize=7, leading=9,
        alignment=TA_RIGHT, fontName='Helvetica-Oblique', wordWrap='CJK')
    mini_style = ParagraphStyle(
        'MiniStyle', parent=small_style, fontSize=7, leading=9,
        wordWrap='CJK')
    mini_header_style = ParagraphStyle(
        'MiniHeaderStyle', parent=small_style, fontSize=7, leading=9,
        fontName='Helvetica-Bold', wordWrap='CJK')

    if available_width is None:
        available_width = 8.5 * inch - 0.02 * 8.5 * inch

    def _has_any_inline_strain_tables():
        if not (use_subkey and show_strains_table):
            return False
        strains_by_group = {}
        for s in all_strains:
            sg = species_group_map[id(s)]
            gk = sg.get('toplevelkey', sg.get('key', 'unknown'))
            strains_by_group.setdefault(gk, []).append(s)
        for gk, group_strains in strains_by_group.items():
            if not group_strains:
                continue
            total_group_members = len(group_strains)
            subkey_groups = {}
            for s in group_strains:
                sk = str(s.get('subkey', s.get('key', 'unknown')))
                subkey_groups.setdefault(sk, []).append(s)
            for sk, members in subkey_groups.items():
                has_multiple_members = len(members) > 1
                has_different_keys = any(str(m.get('key', '')) != str(sk) for m in members)
                sublevel_members = [m for m in members if str(m.get('key', '')) != str(sk)]
                is_name_switch = total_group_members == 1 and has_different_keys
                has_appendix_entry = (
                    (has_multiple_members or has_different_keys)
                    and sublevel_members
                    and not is_name_switch
                )
                if has_appendix_entry:
                    return True
        return False

    show_inline_table = use_subkey and show_strains_table and _has_any_inline_strain_tables()

    # ── Headers ──────────────────────────────────────────────────────────────
    header_style = ParagraphStyle(
        'HeaderStyle', parent=small_style, fontSize=8, leading=9,
        fontName='Helvetica-Bold', textColor=colors.whitesmoke,
        alignment=TA_CENTER, wordWrap='CJK')
    base_headers = ['', Paragraph('Organism', header_style), Paragraph('TASS', header_style)]
    if show_k2_column:
        base_headers.append(Paragraph('K2<br/>Reads', header_style))
    base_headers += [
        Paragraph('Reads', header_style),
        Paragraph('RPM', header_style),
        Paragraph('Cov.', header_style),
    ]
    if show_ani_column:
        base_headers.append(Paragraph('High<br/>ANI', header_style))
    n_base = len(base_headers)

    # Strain detail column in subkey mode (single column holding a nested mini-table).
    # In compact/appendix mode the arrow is embedded inline in the name cell,
    # so no extra column is needed.
    if show_inline_table:
        headers = base_headers + ['']
    else:
        headers = base_headers
    n_total = len(headers)

    # ── Column widths ─────────────────────────────────────────────────────────
    #   Columns: [indicator, Organism, TASS, (K2 Reads), Reads, RPM, Coverage, (High ANI)]
    def _base_col_widths(w):
        if show_k2_column and show_ani_column:
            return [w*0.03, w*0.24, w*0.08, w*0.11, w*0.14, w*0.10, w*0.10, w*0.20]
        elif show_k2_column:
            return [w*0.03, w*0.30, w*0.09, w*0.13, w*0.15, w*0.13, w*0.17]
        elif show_ani_column:
            return [w*0.03, w*0.28, w*0.09, w*0.15, w*0.12, w*0.12, w*0.21]
        else:
            return [w*0.03, w*0.40, w*0.10, w*0.17, w*0.14, w*0.16]

    def build_metrics_mini_table(max_cds, max_mmbert, max_width=None):
        """
        Create a proper table showing CDS and/or mmbert with column headers.
        Returns None if neither metric is present.
        Hides columns when data isn't available.
        max_width constrains the total table width to fit within a container.
        """
        cds_v = _valid_num(max_cds)
        mm_v = _valid_num(max_mmbert)

        # If neither metric is present, return None
        if cds_v is None and mm_v is None:
            return None

        # Build header and data rows based on what's available
        rows = []

        # Determine which columns to include
        has_cds = cds_v is not None
        has_mmbert = mm_v is not None

        # Build header
        header_row = []
        if has_cds:
            header_row.append(Paragraph(f"<b>CDS</b>", mini_header_style))
        if has_mmbert:
            header_row.append(Paragraph(f"<b>mmbert %</b>", mini_header_style))
        rows.append(header_row)

        # Build data row
        data_row = []
        if has_cds:
            data_row.append(Paragraph(f"{cds_v:.0f}", mini_style))
        if has_mmbert:
            mmbert_pct = mm_v * 100  # Convert to percentage
            data_row.append(Paragraph(f"{mmbert_pct:.2f}%", mini_style))
        rows.append(data_row)

        # Calculate column widths, respecting max_width if provided
        desired_cds_w = 0.5 * inch
        desired_mm_w = 0.65 * inch
        col_widths_metrics = []
        if has_cds:
            col_widths_metrics.append(desired_cds_w)
        if has_mmbert:
            col_widths_metrics.append(desired_mm_w)

        # Scale down if total exceeds max_width
        if max_width is not None:
            total_desired = sum(col_widths_metrics)
            if total_desired > max_width:
                scale = max_width / total_desired
                col_widths_metrics = [w * scale for w in col_widths_metrics]

        # Create table with only the necessary columns
        t = Table(rows, colWidths=col_widths_metrics)
        t.setStyle(TableStyle([
            ('FONTSIZE', (0, 0), (-1, -1), 7),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#D5E8FF')),  # Header background
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('TOPPADDING', (0, 0), (-1, -1), 2),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 2),
            ('LEFTPADDING', (0, 0), (-1, -1), 3),
            ('RIGHTPADDING', (0, 0), (-1, -1), 3),
            ('VALIGN', (0, 0), (-1, -1), 'TOP'),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
        ]))
        return t
    if show_inline_table:
        left_w = available_width * 0.62
        right_w = available_width * 0.38
        col_widths = _base_col_widths(left_w) + [right_w]
    else:
        # Compact/appendix mode: no separate right column — arrow is inline in name cell
        right_w = available_width * 0.38  # kept for metrics_max_w calc but not added as column
        col_widths = _base_col_widths(available_width)

    # ── Table data ────────────────────────────────────────────────────────────
    table_data = [headers]
    table_styles = []
    group_row_indices = []

    current_species_key = None
    row_idx = 1  # row 0 is the header
    emitted_subkeys_per_group = {}  # {species_key: set of already-emitted subkeys}

    # Pre-compute strain names per species group for group header summaries
    strains_per_group = {}
    for s in all_strains:
        sg = species_group_map[id(s)]
        gk = sg.get('toplevelkey', sg.get('key', 'unknown'))
        strains_per_group.setdefault(gk, []).append(s.get('name', 'Unknown'))
    group_mmbert_max = {}
    group_dmnd_cds_max = {}
    for s in all_strains:
        sg = species_group_map[id(s)]
        gk = sg.get('toplevelkey', sg.get('key', 'unknown'))

        # Use only the GROUP-level mmbert, not member-level
        group_mm = _valid_num(sg.get('mmbert', None))
        if group_mm is not None:
            group_mmbert_max[gk] = group_mm

        best_dmnd = None
        for c in [
            (sg.get('diamond') or {}).get('cds', None),
            (s.get('diamond') or {}).get('cds', None)
        ]:
            c_v = _valid_num(c)
            if c_v is not None:
                best_dmnd = c_v if best_dmnd is None else max(best_dmnd, c_v)
        if best_dmnd is not None:
            group_dmnd_cds_max[gk] = best_dmnd
    i = 0
    while i < len(all_strains):
        strain = all_strains[i]
        species_group = species_group_map[id(strain)]
        species_key = species_group.get('toplevelkey', species_group.get('key', 'unknown'))

        # ── Species group header row ──────────────────────────────────────────
        if current_species_key != species_key:
            current_species_key = species_key
            emitted_subkeys_per_group[species_key] = set()

            group_name = species_group.get('toplevelname', 'Unknown')
            group_reads = species_group.get('numreads', 0)
            group_k2_reads = species_group.get('k2_reads', 0)

            group_strain_names = strains_per_group.get(species_key, [])
            n_strains = len(group_strain_names)

            group_name_para = Paragraph(
                f'<b><link href="https://www.ncbi.nlm.nih.gov/taxonomy/?term={species_key}" '
                f'color="blue">{group_name}</link></b>',
                group_header_style
            )

            # Group-level metrics: CDS/mmbert (group-level only)
            gm = group_mmbert_max.get(species_key, None)
            dmd_cds = group_dmnd_cds_max.get(species_key, None)
            # Pass right column content width (minus outer cell padding) to cap mini-table
            metrics_max_w = (right_w - 6) if use_subkey else None
            group_metrics_tbl = build_metrics_mini_table(dmd_cds, gm, max_width=metrics_max_w)

            # The group name cell spans columns 1+2
            spanned_name_width = col_widths[1] + col_widths[2]

            # Name cell contains ONLY the name paragraph (metrics moved to right column)
            # Subtract outer cell padding (LEFT=4 + RIGHT=4) so text wraps within bounds
            name_cell_w = spanned_name_width - 8
            group_name_cell = Table(
                [[group_name_para]],
                colWidths=[name_cell_w]
            )
            group_name_cell.setStyle(TableStyle([
                ('LEFTPADDING', (0, 0), (-1, -1), 0),
                ('RIGHTPADDING', (0, 0), (-1, -1), 0),
                ('TOPPADDING', (0, 0), (-1, -1), 0),
                ('BOTTOMPADDING', (0, 0), (-1, -1), 0),
                ('VALIGN', (0, 0), (-1, -1), 'TOP'),
            ]))

            group_row = ['', group_name_cell, '']
            if show_k2_column:
                group_row.append(Paragraph(f'<b>{group_k2_reads:,.0f}</b>', mini_style))
            group_row.append(Paragraph(f'<b>{group_reads:,.0f}</b>', mini_style))
            group_row += ['', '']
            if show_ani_column:
                group_row.append('')

            # Right-side detail column (inline mode only — compact mode has no extra column)
            if show_inline_table:
                # Full inline mode: show CDS/mmbert metrics (organism count is now in the spanned area)
                if group_metrics_tbl is not None:
                    right_cell_rows = [[group_metrics_tbl]]
                    right_cell = Table(right_cell_rows, colWidths=[right_w - 6])
                    right_cell.setStyle(TableStyle([
                        ('LEFTPADDING', (0, 0), (-1, -1), 0),
                        ('RIGHTPADDING', (0, 0), (-1, -1), 0),
                        ('TOPPADDING', (0, 0), (-1, -1), 0),
                        ('BOTTOMPADDING', (0, 0), (-1, -1), 0),
                        ('VALIGN', (0, 0), (-1, -1), 'TOP'),
                        ('ALIGN', (0, 0), (-1, -1), 'RIGHT'),
                    ]))
                    group_row.append(right_cell)
                else:
                    group_row.append('')

            table_data.append(group_row)
            group_row_indices.append(row_idx)

            # Span group name across Organism + TASS columns only
            table_styles.append(('SPAN', (1, row_idx), (2, row_idx)))
            # Span trailing columns: stop before the right detail column when it exists
            rpm_col_idx = 5 if show_k2_column else 4
            span_end = (n_base - 1) if show_inline_table else (n_total - 1)
            if rpm_col_idx <= span_end:
                table_styles.append(('SPAN', (rpm_col_idx, row_idx), (span_end, row_idx)))
            # Always place organism count summary in the RPM slot (first cell of trailing span)
            summary_text = f'<i>{n_strains} detected organism{"s" if n_strains != 1 else ""}</i>'
            group_row[rpm_col_idx] = Paragraph(summary_text, group_strain_summary_style)
            table_styles.append(('BACKGROUND', (0, row_idx), (-1, row_idx), colors.HexColor('#E8E8E8')))
            table_styles.append(('ALIGN', (1, row_idx), (1, row_idx), 'LEFT'))
            table_styles.append(('TOPPADDING', (0, row_idx), (-1, row_idx), 8))
            table_styles.append(('BOTTOMPADDING', (0, row_idx), (-1, row_idx), 8))
            row_idx += 1

        # ── Data rows ─────────────────────────────────────────────────────────
        if use_subkey:
            sk = str(strain.get('subkey', strain.get('key', 'unknown')))

            # Skip if we already emitted this subkey for this group
            if sk in emitted_subkeys_per_group[species_key]:
                i += 1
                continue

            # Collect ALL members in this group sharing this subkey
            subkey_members = [
                s for s in all_strains
                if (species_group_map[id(s)].get(
                        'toplevelkey', species_group_map[id(s)].get('key', 'unknown')
                    ) == species_key
                    and str(s.get('subkey', s.get('key', 'unknown'))) == sk)
            ]
            subkey_members.sort(key=lambda s: s.get('tass_score', 0), reverse=True)
            best = subkey_members[0]
            emitted_subkeys_per_group[species_key].add(sk)

            # Left columns driven by best-TASS member
            microbial_category = best.get('microbial_category', 'Unknown')
            ann_class = best.get('annClass', '')
            is_hc = best.get('high_cons', False)
            indicator_text = '★' if is_hc else ''

            best_key = best.get('key', '')

            # Count total visible members in this species group to decide
            # whether the strain mini-table will be shown or suppressed.
            total_group_members = [
                s for s in all_strains
                if species_group_map[id(s)].get(
                    'toplevelkey', species_group_map[id(s)].get('key', 'unknown')
                ) == species_key
            ]
            single_strain_in_group = len(total_group_members) == 1

            # When there is only one strain (strain table removed), show
            # the member's own name/key instead of the subkey-level name
            # so the row reads e.g. "SARS_CoV-2" rather than the species
            # subkey name which duplicates the group header.
            if single_strain_in_group:
                display_name = best.get('name', 'Unknown')
                display_key = best_key
            elif str(best.get('subkey', best_key)) != str(best_key):
                # In subkey mode, show the species-level subkeyname when key != subkey
                display_name = best.get('subkeyname', best.get('name', 'Unknown'))
                display_key = str(best.get('subkey', best_key))
            else:
                display_name = best.get('name', 'Unknown')
                display_key = best_key
            # Base name HTML (arrow suffix added below once we know if one is needed)
            name_html_base = (
                f'{display_name} '
                f'(<link href="https://www.ncbi.nlm.nih.gov/taxonomy/?term={display_key}" '
                f'color="blue">{display_key}</link>)'
            )

            # ── Determine whether sub-level strains exist ──────────────────────
            # Must be done before building the row so we can embed the arrow in
            # the name cell (compact mode) or build the inline mini-table.
            has_multiple_members = len(subkey_members) > 1
            # Compare each member's key to sk (the subkey), NOT to best.get('key').
            # This correctly detects members that are genuine sub-level strains
            # (key != subkey), including the single-member case where the main
            # table does a name-switch (single_strain_in_group).
            has_different_keys = any(
                str(m.get('key', '')) != str(sk) for m in subkey_members
            )
            subkey_members = [m for m in subkey_members if str(m.get('key', '')) != str(sk)]

            # When the group has exactly one visible strain AND key != subkey,
            # the main table has already switched the display to show the strain
            # directly (name-switch).  There is nothing additional to show in
            # the appendix, and no arrow should be rendered.
            is_name_switch = single_strain_in_group and has_different_keys

            has_appendix_entry = (
                (has_multiple_members or has_different_keys)
                and subkey_members
                and not is_name_switch
            )

            # ── Build name cell with optional inline ↓ arrow (compact mode) ───
            if not show_inline_table and has_appendix_entry:
                strain_row_bm = f"strain_row_{sanitize_bookmark_name(sample_name)}_{species_key}"
                arrow_link = create_safe_link('Pathogenic strain detail \u2193', strain_row_bm,
                                              valid_bookmarks, color='blue')
                name_html = f'{name_html_base} {arrow_link}'
            else:
                name_html = name_html_base

            strain_reads = float(best.get('numreads', 0) or 0)
            pct = (strain_reads / sample_total_reads * 100.0) if sample_total_reads else 0.0
            rpm = best.get('rpm', 0) or 0

            high_ani_text = ''
            if show_ani_column:
                # Use the pre-computed matches embedded by match_paths.py
                matches = get_high_ani_matches(best)
                if matches:
                    ani_links = []
                    for m in matches[:3]:
                        taxid = str(m.get('key', ''))
                        ani_pct = float(m.get('ani_pct', 0))
                        if taxid in taxid_to_bookmark:
                            bm = taxid_to_bookmark[taxid]
                            if bm in valid_bookmarks:
                                ani_links.append(
                                    f'<link href="#{bm}" color="blue">{taxid} ({ani_pct:.1f}%)</link>')
                            else:
                                ani_links.append(f"{taxid} ({ani_pct:.1f}%)")
                        else:
                            ani_links.append(f"{taxid} ({ani_pct:.1f}%)")
                    high_ani_text = Paragraph(", ".join(ani_links), ani_style)
                else:
                    high_ani_text = Paragraph("-", ani_style)

            row = [
                indicator_text,
                Paragraph(name_html, strain_name_style),
                Paragraph(f"{best.get('tass_score', 0)*100:.1f}", data_style),
            ]
            if show_k2_column:
                row.append(Paragraph(f"{best.get('k2_reads', 0):,.0f}", data_style))
            row.append(Paragraph(f"{strain_reads:,.0f} ({pct:.1f}%)", data_style))
            row.append(Paragraph(f"{rpm:,.0f}", data_style))
            row.append(Paragraph(f"{min(100, best.get('coverage', 0)*100):.1f}%", data_style))
            if show_ani_column:
                row.append(high_ani_text)

            # ── Right detail column (inline mini-table, show_strains_table only) ─
            if show_inline_table:
                if has_appendix_entry:
                    mini_rows = [[
                        Paragraph('<b>Additional Strain/Subsp.</b>', mini_header_style),
                        Paragraph('<b>TASS</b>', mini_header_style),
                        Paragraph('<b>Reads Aligned</b>', mini_header_style),
                    ]]
                    for m in subkey_members:
                        m_reads = float(m.get('numreads', 0) or 0)
                        m_pct = (m_reads / sample_total_reads * 100.0) if sample_total_reads else 0.0
                        m_key = m.get('key', '')
                        m_name = m.get('name', 'Unknown')
                        m_star = '★ ' if m.get('high_cons', False) else ''
                        mini_rows.append([
                            Paragraph(
                                f'{m_star}<link href="https://www.ncbi.nlm.nih.gov/taxonomy/?term={m_key}" '
                                f'color="blue">{m_name}</link>',
                                mini_style
                            ),
                            Paragraph(f"{m.get('tass_score', 0)*100:.1f}", mini_style),
                            Paragraph(f"{m_reads:,.0f} ({m_pct:.1f}%)", mini_style),
                        ])
                    mini_avail = right_w - 6
                    mini_col_widths = [mini_avail * 0.50, mini_avail * 0.20, mini_avail * 0.30]
                    mini_tbl = Table(mini_rows, colWidths=mini_col_widths)
                    mini_tbl.setStyle(TableStyle([
                        ('FONTSIZE', (0, 0), (-1, -1), 7),
                        ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
                        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#D5E8FF')),
                        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                        ('TOPPADDING', (0, 0), (-1, -1), 2),
                        ('BOTTOMPADDING', (0, 0), (-1, -1), 2),
                        ('LEFTPADDING', (0, 0), (-1, -1), 3),
                        ('RIGHTPADDING', (0, 0), (-1, -1), 3),
                        ('VALIGN', (0, 0), (-1, -1), 'TOP'),
                        ('ALIGN', (1, 0), (2, -1), 'CENTER'),
                    ]))
                    row.append(mini_tbl)
                else:
                    row.append('')
            # Compact mode: no extra column — arrow already embedded in name_html above

            table_data.append(row)

            print(f"  Row {row_idx}: [sk={sk}] {best.get('name', '')[:40]} "
                  f"({len(subkey_members)} member(s)) - Cat: {microbial_category}")

            ind_color = get_category_color(microbial_category, ann_class, alpha=1.0)
            row_color = get_category_color(microbial_category, ann_class, alpha=0.15)
            table_styles.append(('BACKGROUND', (0, row_idx), (0, row_idx), ind_color))
            table_styles.append(('BACKGROUND', (1, row_idx), (-1, row_idx), row_color))
            # Horizontal separator below each strain row with padding
            table_styles.append(('LINEBELOW', (0, row_idx), (-1, row_idx), 1.5, colors.HexColor('#CCCCCC')))
            table_styles.append(('BOTTOMPADDING', (0, row_idx), (-1, row_idx), 8))
            table_styles.append(('TOPPADDING', (0, row_idx), (-1, row_idx), 6))
            row_idx += 1
            i += 1

        else:
            # ── Original flat per-strain mode ──────────────────────────────────
            microbial_category = strain.get('microbial_category', 'Unknown')
            ann_class = strain.get('annClass', '')
            is_hc = strain.get('high_cons', False)
            indicator_text = '★' if is_hc else ''

            strain_key = strain.get('key', '')
            name_html = (
                f'{strain.get("name", "Unknown")} '
                f'(<link href="https://www.ncbi.nlm.nih.gov/taxonomy/?term={strain_key}" '
                f'color="blue">{strain_key}</link>)'
            )

            high_ani_text = ''
            if show_ani_column:
                # Use the pre-computed matches embedded by match_paths.py
                matches = get_high_ani_matches(strain)
                if matches:
                    ani_links = []
                    for m in matches[:3]:
                        taxid = str(m.get('key', ''))
                        ani_pct = float(m.get('ani_pct', 0))
                        if taxid in taxid_to_bookmark:
                            bm = taxid_to_bookmark[taxid]
                            if bm in valid_bookmarks:
                                ani_links.append(
                                    f'<link href="#{bm}" color="blue">{taxid} ({ani_pct:.1f}%)</link>')
                            else:
                                ani_links.append(f"{taxid} ({ani_pct:.1f}%)")
                        else:
                            ani_links.append(f"{taxid} ({ani_pct:.1f}%)")
                    high_ani_text = Paragraph(", ".join(ani_links), ani_style)
                else:
                    high_ani_text = Paragraph("-", ani_style)

            strain_reads = float(strain.get('numreads', 0) or 0)
            pct = (strain_reads / sample_total_reads * 100.0) if sample_total_reads else 0.0
            rpm = strain.get('rpm', 0) or 0

            row = [
                indicator_text,
                Paragraph(name_html, strain_name_style),
                Paragraph(f"{strain.get('tass_score', 0)*100:.1f}", data_style),
            ]
            if show_k2_column:
                row.append(Paragraph(f"{strain.get('k2_reads', 0):,.0f}", data_style))
            row.append(Paragraph(f"{strain_reads:,.0f} ({pct:.1f}%)", data_style))
            row.append(Paragraph(f"{rpm:,.0f}", data_style))
            row.append(Paragraph(f"{min(100, strain.get('coverage', 0)*100):.1f}%", data_style))
            if show_ani_column:
                row.append(high_ani_text)

            table_data.append(row)

            print(f"  Row {row_idx}: {strain.get('name', 'Unknown')[:40]} "
                  f"- Category: {microbial_category}, Class: {ann_class}")

            ind_color = get_category_color(microbial_category, ann_class, alpha=1.0)
            row_color = get_category_color(microbial_category, ann_class, alpha=0.15)
            table_styles.append(('BACKGROUND', (0, row_idx), (0, row_idx), ind_color))
            table_styles.append(('BACKGROUND', (1, row_idx), (-1, row_idx), row_color))
            # Horizontal separator below each strain row with padding
            table_styles.append(('LINEBELOW', (0, row_idx), (-1, row_idx), 1.5, colors.HexColor('#CCCCCC')))
            table_styles.append(('BOTTOMPADDING', (0, row_idx), (-1, row_idx), 8))
            table_styles.append(('TOPPADDING', (0, row_idx), (-1, row_idx), 6))
            row_idx += 1
            i += 1

    # ── Base table style ──────────────────────────────────────────────────────
    base_ts = [
        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#3498DB')),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
        ('ALIGN', (0, 0), (-1, 0), 'CENTER'),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, 0), 9),
        ('BOTTOMPADDING', (0, 0), (-1, 0), 6),
        ('TOPPADDING', (0, 0), (-1, 0), 6),
        ('LEFTPADDING', (0, 0), (-1, 0), 3),
        ('RIGHTPADDING', (0, 0), (-1, 0), 3),
        ('ALIGN', (0, 1), (0, -1), 'CENTER'),
        ('VALIGN', (0, 1), (0, -1), 'MIDDLE'),
        ('FONTSIZE', (0, 1), (0, -1), 14),
        ('GRID', (0, 0), (-1, -1), 1, colors.black),
        ('VALIGN', (1, 1), (1, -1), 'MIDDLE'),
        ('VALIGN', (2, 1), (-1, -1), 'MIDDLE'),
        ('ALIGN', (2, 1), (-1, -1), 'CENTER'),
        ('ALIGN', (1, 1), (1, -1), 'LEFT'),
        # Organism column (col 1): keep comfortable padding
        ('LEFTPADDING', (1, 1), (1, -1), 4),
        ('RIGHTPADDING', (1, 1), (1, -1), 4),
        ('TOPPADDING', (1, 1), (1, -1), 6),
        ('BOTTOMPADDING', (1, 1), (1, -1), 6),
        # Data columns (col 2 onwards): tight padding to prevent overflow
        ('LEFTPADDING', (2, 1), (n_base - 1, -1), 2),
        ('RIGHTPADDING', (2, 1), (n_base - 1, -1), 2),
        ('TOPPADDING', (2, 1), (n_base - 1, -1), 4),
        ('BOTTOMPADDING', (2, 1), (n_base - 1, -1), 4),
    ]

    if show_inline_table:
        # Right detail column styles + vertical delimiter (inline mode only)
        right_col_ts = [
            ('ALIGN', (n_base, 1), (n_total - 1, -1), 'LEFT'),
            ('VALIGN', (n_base, 1), (n_total - 1, -1), 'TOP'),
            ('LEFTPADDING', (n_base, 1), (n_total - 1, -1), 4),
            ('TOPPADDING', (n_base, 1), (n_total - 1, -1), 2),
            ('RIGHTPADDING', (n_base, 1), (n_total - 1, -1), 2),
            ('BOTTOMPADDING', (n_base, 1), (n_total - 1, -1), 2),
            ('LINEAFTER', (n_base - 1, 0), (n_base - 1, -1), 2, colors.HexColor('#3498DB')),
        ]
        base_ts += right_col_ts

    table = Table(table_data, repeatRows=1, colWidths=col_widths)
    table.setStyle(TableStyle(base_ts + table_styles))
    return table


def create_low_confidence_table(low_confidence_strains, small_style, show_k2_column,
                                 available_width=None):
    strain_name_style = ParagraphStyle(
        'StrainName', parent=small_style, fontSize=8, leading=10)
    data_style = ParagraphStyle(
        'DataStyle', parent=small_style, fontSize=7, leading=9)

    headers = (['Sample', 'Organism', 'TASS', 'K2 Reads', 'Reads']
               if show_k2_column else ['Sample', 'Organism', 'TASS', 'Reads'])
    table_data = [headers]

    for item in low_confidence_strains:
        strain = item['strain']
        sample_name = item['sample_name']

        strain_name_text = strain.get('name', 'Unknown')
        strain_key = strain.get('key', '')
        if strain_key:
            strain_name_text = (
                f'{strain_name_text} '
                f'(<link href="https://www.ncbi.nlm.nih.gov/taxonomy/?term={strain_key}" '
                f'color="blue">{strain_key}</link>)'
            )
        if strain.get('high_cons', False):
            strain_name_text = f"★ {strain_name_text}"

        row = [
            Paragraph(sample_name, strain_name_style),
            Paragraph(strain_name_text, strain_name_style),
            Paragraph(f"{strain.get('tass_score', 0)*100:.1f}", data_style),
        ]
        if show_k2_column:
            row.append(Paragraph(f"{strain.get('k2_reads', 0):,.0f}", data_style))
        row.append(Paragraph(f"{strain.get('numreads', 0):,.0f}", data_style))
        table_data.append(row)

    if available_width is None:
        available_width = 8.5*inch - 0.02*8.5*inch

    if show_k2_column:
        col_widths = [available_width*0.18, available_width*0.50,
                      available_width*0.10, available_width*0.10, available_width*0.12]
    else:
        col_widths = [available_width*0.20, available_width*0.56,
                      available_width*0.11, available_width*0.13]

    table = Table(table_data, repeatRows=1, colWidths=col_widths)
    table_styles = [
        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#3498DB')),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
        ('ALIGN', (0, 0), (-1, 0), 'CENTER'),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, 0), 9),
        ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
        ('GRID', (0, 0), (-1, -1), 1, colors.black),
        ('VALIGN', (0, 1), (-1, -1), 'MIDDLE'),
        ('ALIGN', (2, 1), (-1, -1), 'CENTER'),
        ('ALIGN', (0, 1), (1, -1), 'LEFT'),
        ('LEFTPADDING', (0, 1), (-1, -1), 6),
        ('RIGHTPADDING', (0, 1), (-1, -1), 6),
        ('TOPPADDING', (0, 1), (-1, -1), 6),
        ('BOTTOMPADDING', (0, 1), (-1, -1), 6),
    ]
    for idx, item in enumerate(low_confidence_strains):
        row_color = get_category_color(
            item['strain'].get('microbial_category', 'Unknown'),
            item['strain'].get('annClass', ''),
            alpha=0.25
        )
        table_styles.append(('BACKGROUND', (0, idx+1), (-1, idx+1), row_color))
    table.setStyle(TableStyle(table_styles))
    return table


def create_strain_detail_tables(samples_dict, sorted_groups_by_sample,
                                small_style, show_k2_column, valid_bookmarks,
                                args, available_width):
    """
    Build the strain-detail appendix story elements.

    Returns a list of ReportLab flowables to be appended to the main story.
    One flat table per sample; each row is a single strain member.
    Rows are compact (small font, tight padding) to save space.
    Each species-group block in the table is preceded by an AnchorFlowable
    so that links from the main table can jump directly to it.
    """
    story_items = []

    # ── Section heading ───────────────────────────────────────────────────────
    heading_style = ParagraphStyle(
        'StrainDetailHeading', parent=small_style, fontSize=14,
        textColor=colors.HexColor('#34495E'), spaceBefore=10, spaceAfter=4,
        fontName='Helvetica-Bold')
    note_style = ParagraphStyle(
        'StrainDetailNote', parent=small_style, fontSize=8, leading=10,
        textColor=colors.HexColor('#555555'))
    sample_heading_style = ParagraphStyle(
        'StrainSampleHeading', parent=small_style, fontSize=11,
        fontName='Helvetica-Bold', textColor=colors.HexColor('#2C3E50'),
        spaceBefore=8, spaceAfter=2)
    cell_style = ParagraphStyle(
        'StrainCell', parent=small_style, fontSize=7, leading=8, wordWrap='CJK')
    header_cell_style = ParagraphStyle(
        'StrainHeaderCell', parent=small_style, fontSize=7, leading=8,
        fontName='Helvetica-Bold', textColor=colors.whitesmoke,
        alignment=TA_CENTER, wordWrap='CJK')

    story_items.append(AnchorFlowable('strain_detail_table'))
    story_items.append(Paragraph('<b>Additional Lower-Level Strain Detail</b>', heading_style))
    story_items.append(Paragraph(
        'This appendix lists additional lower-level strains, subspecies, and assemblies '
        'detected beneath each species-level entry in the main table above. '
        'Only species groups with more than one distinct strain-level member are shown here. '
        'Links labeled "Pathogenic strain detail" jump between the main table and this '
        'appendix when pathogenic strain-level entries are present.',
        note_style))
    story_items.append(Spacer(1, 0.06 * inch))

    for sample_name, sorted_groups in sorted_groups_by_sample:
        # Pre-scan: does this sample have any strains with CDS or mmbert?
        sample_has_cds = False
        sample_has_mmbert = False
        sample_strain_rows = []  # list of (sg, strain) pairs to render

        for sg in sorted_groups:
            species_key = sg.get('toplevelkey', sg.get('key', 'unknown'))
            qualifying = [
                m for m in sg.get('members', [])
                if (should_include_strain(m, args)
                    and has_min_reads(m, args.min_reads)
                    and passes_confidence_threshold(m, args.min_conf))
            ]
            # Mirror the main table's right-column visibility logic exactly.
            # A group belongs in the appendix only when it has genuinely
            # additional lower-level strains beyond what is already shown in
            # the main table row.
            #
            # Two cases to exclude:
            #   1. Every member's key == its own subkey — nothing sub-level.
            #   2. Single-strain name-switch: exactly one qualifying member
            #      whose key != subkey.  In this case the main table row has
            #      already switched its display name to that strain, so there
            #      is nothing extra to list in the appendix.
            has_sublevel_members = any(
                str(m.get('key', '')) != str(m.get('subkey', m.get('key', '')))
                for m in qualifying
            )
            if not has_sublevel_members:
                continue
            # Detect single-strain name-switch (mirrors single_strain_in_group
            # logic from create_combined_sample_table): only one visible member
            # in the entire species group and its key differs from its subkey.
            if len(qualifying) == 1 and str(qualifying[0].get('key', '')) != str(qualifying[0].get('subkey', qualifying[0].get('key', ''))):
                continue
            # Only append members that are genuine sub-level strains (key != subkey).
            # Members where key == subkey are already shown as their own species row
            # in the main table and must not appear again in the appendix.
            sublevel = [
                m for m in qualifying
                if str(m.get('key', '')) != str(m.get('subkey', m.get('key', '')))
            ]
            for strain in sublevel:
                sample_strain_rows.append((sg, strain))
                cds_v = _valid_num((strain.get('diamond') or {}).get('cds', None))
                mm_v = _valid_num(strain.get('mmbert', None))
                if cds_v is not None:
                    sample_has_cds = True
                if mm_v is not None:
                    sample_has_mmbert = True

        if not sample_strain_rows:
            continue

        # ── Sample heading with back-link ─────────────────────────────────────
        strain_sample_bm = f"strain_sample_{sanitize_bookmark_name(sample_name)}"
        main_sample_bm = f"sample_{sanitize_bookmark_name(sample_name)}"
        back_link = create_safe_link('\u2191', main_sample_bm, valid_bookmarks, color='blue')
        story_items.append(AnchorFlowable(strain_sample_bm))
        story_items.append(Paragraph(
            f'<b>{sample_name}</b> {back_link}', sample_heading_style))

        # ── Column headers ────────────────────────────────────────────────────
        col_headers = [
            '',  # indicator
            Paragraph('Organism', header_cell_style),
            Paragraph('TASS', header_cell_style),
        ]
        if show_k2_column:
            col_headers.append(Paragraph('K2<br/>Reads', header_cell_style))
        col_headers += [
            Paragraph('Reads', header_cell_style),
            Paragraph('Cov.', header_cell_style),
        ]
        if sample_has_cds:
            col_headers.append(Paragraph('CDS', header_cell_style))
        if sample_has_mmbert:
            col_headers.append(Paragraph('mmbert%', header_cell_style))

        # ── Column widths ─────────────────────────────────────────────────────
        n_extra = (1 if sample_has_cds else 0) + (1 if sample_has_mmbert else 0)
        n_k2 = 1 if show_k2_column else 0
        # Distribute: indicator(2%), organism(rest), TASS(7%), K2(8%), Reads(12%), Cov(7%), CDS(7%), mmbert(8%)
        fixed_pcts = 0.02 + 0.07 + n_k2 * 0.08 + 0.12 + 0.07 + n_extra * 0.075
        org_pct = max(0.25, 1.0 - fixed_pcts)

        w = available_width
        cw = [w * 0.02, w * org_pct, w * 0.07]
        if show_k2_column:
            cw.append(w * 0.08)
        cw += [w * 0.12, w * 0.07]
        if sample_has_cds:
            cw.append(w * 0.075)
        if sample_has_mmbert:
            cw.append(w * 0.075)

        # ── Build table data ──────────────────────────────────────────────────
        table_data = [col_headers]
        table_styles_det = []
        row_idx = 1

        emitted_species = set()

        for sg, strain in sample_strain_rows:
            species_key = sg.get('toplevelkey', sg.get('key', 'unknown'))
            species_name = sg.get('toplevelname', '')

            # Place per-species anchor before the first strain row for that group
            if species_key not in emitted_species:
                emitted_species.add(species_key)
                strain_row_bm = f"strain_row_{sanitize_bookmark_name(sample_name)}_{species_key}"
                story_items.append(AnchorFlowable(strain_row_bm))

            microbial_category = strain.get('microbial_category', 'Unknown')
            ann_class = strain.get('annClass', '')
            is_hc = strain.get('high_cons', False)
            indicator = '★' if is_hc else ''

            strain_key = strain.get('key', '')
            strain_name = strain.get('name', 'Unknown')

            # Back-link ↑ to the species row in the main table
            species_bm = f"species_{sanitize_bookmark_name(sample_name)}_{species_key}"
            up_link = create_safe_link('Pathogenic strain detail \u2191', species_bm,
                                       valid_bookmarks, color='blue')

            name_html = (
                f'<link href="https://www.ncbi.nlm.nih.gov/taxonomy/?term={strain_key}" '
                f'color="blue">{strain_name}</link> {up_link}'
            )

            strain_reads = float(strain.get('numreads', 0) or 0)
            sample_total = max(1, sum(
                sg2.get('numreads', 0)
                for sg2 in samples_dict.get(sample_name, [])
            ))
            pct = strain_reads / sample_total * 100.0

            tass_val = strain.get('tass_score', 0) * 100

            row = [
                indicator,
                Paragraph(name_html, cell_style),
                Paragraph(f'{tass_val:.1f}', cell_style),
            ]
            if show_k2_column:
                row.append(Paragraph(f"{strain.get('k2_reads', 0):,.0f}", cell_style))
            row += [
                Paragraph(f'{strain_reads:,.0f} ({pct:.1f}%)', cell_style),
                Paragraph(f'{min(100, strain.get("coverage", 0) * 100):.1f}%', cell_style),
            ]
            if sample_has_cds:
                cds_v = _valid_num((strain.get('diamond') or {}).get('cds', None))
                row.append(Paragraph(f'{cds_v:.0f}' if cds_v is not None else '—', cell_style))
            if sample_has_mmbert:
                mm_v = _valid_num(strain.get('mmbert', None))
                row.append(Paragraph(
                    f'{mm_v * 100:.2f}%' if mm_v is not None else '—', cell_style))

            table_data.append(row)

            # Row styling
            ind_color = get_category_color(microbial_category, ann_class, alpha=1.0)
            row_color = get_category_color(microbial_category, ann_class, alpha=0.12)
            table_styles_det.append(('BACKGROUND', (0, row_idx), (0, row_idx), ind_color))
            table_styles_det.append(('BACKGROUND', (1, row_idx), (-1, row_idx), row_color))
            table_styles_det.append(('LINEBELOW', (0, row_idx), (-1, row_idx),
                                     0.5, colors.HexColor('#DDDDDD')))
            table_styles_det.append(('TOPPADDING', (0, row_idx), (-1, row_idx), 2))
            table_styles_det.append(('BOTTOMPADDING', (0, row_idx), (-1, row_idx), 2))
            row_idx += 1

        # ── Base table style ──────────────────────────────────────────────────
        base_ts = [
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#3498DB')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 7),
            ('TOPPADDING', (0, 0), (-1, 0), 3),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 3),
            ('LEFTPADDING', (0, 0), (-1, 0), 2),
            ('RIGHTPADDING', (0, 0), (-1, 0), 2),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#AAAAAA')),
            ('ALIGN', (0, 0), (-1, 0), 'CENTER'),
            ('ALIGN', (0, 1), (0, -1), 'CENTER'),
            ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
            ('ALIGN', (1, 1), (1, -1), 'LEFT'),
            ('ALIGN', (2, 1), (-1, -1), 'CENTER'),
            ('LEFTPADDING', (1, 1), (1, -1), 3),
            ('RIGHTPADDING', (1, 1), (1, -1), 3),
            ('LEFTPADDING', (2, 1), (-1, -1), 2),
            ('RIGHTPADDING', (2, 1), (-1, -1), 2),
            ('FONTSIZE', (0, 1), (0, -1), 11),  # indicator column larger for ★
        ]

        det_table = Table(table_data, repeatRows=1, colWidths=cw)
        det_table.setStyle(TableStyle(base_ts + table_styles_det))
        story_items.append(det_table)
        story_items.append(Spacer(1, 0.08 * inch))

    return story_items


def create_pdf_template(output_path, samples_dict, args):
    """
    Create the full PDF report.
    """
    page_width, page_height = letter
    left_margin = page_width * 0.01
    right_margin = page_width * 0.01

    doc = SimpleDocTemplate(
        output_path,
        pagesize=letter,
        leftMargin=left_margin,
        rightMargin=right_margin,
        topMargin=0.5*inch,
        bottomMargin=0.5*inch
    )
    available_width = doc.width
    story = []
    styles = getSampleStyleSheet()
    taxid_to_bookmark = build_taxid_to_bookmark_map(samples_dict)

    legend_text_style = ParagraphStyle('LegendText', parent=styles['Normal'], fontSize=8, leading=10)
    title_style = ParagraphStyle(
        'CustomTitle', parent=styles['Heading1'], fontSize=24,
        textColor=colors.HexColor('#2C3E50'), spaceAfter=20, alignment=TA_CENTER)
    heading_style = ParagraphStyle(
        'CustomHeading', parent=styles['Heading2'], fontSize=16,
        textColor=colors.HexColor('#34495E'), spaceAfter=0, spaceBefore=12)
    indent_style = ParagraphStyle(
        'IndentStyle', parent=styles['Normal'], fontSize=10, leading=10, leftIndent=20)
    small_style = ParagraphStyle('SmallText', parent=styles['Normal'], fontSize=10, leading=10)
    metadata_style = styles['Normal']

    # Collect all strains
    all_strains = []
    for sample_name, species_groups in samples_dict.items():
        for species_group in species_groups:
            all_strains.extend(species_group.get('members', []))

    if not args.enable_matrix:
        show_ani_column = False
    else:
        show_ani_column = check_if_any_high_ani_in_dataset(all_strains)
    show_k2_column = check_if_k2_reads_present(all_strains)
    use_subkey = not args.no_subkey
    if use_subkey:
        # Auto-disable subkey columns if ALL strains have subkey == key
        any_subkey_differs = any(
            str(strain.get('subkey', strain.get('key', ''))) != str(strain.get('key', ''))
            for species_groups in samples_dict.values()
            for sg in species_groups
            for strain in sg.get('members', [])
        )
        if not any_subkey_differs:
            use_subkey = False
            print("Auto-disabled subkey grouping: all subkeys match their keys")

    show_strains_table = getattr(args, 'show_strains_table', False)

    print(f"\nHigh ANI column: {'SHOWN' if show_ani_column else 'HIDDEN'} (pre-computed from match_paths.py)")
    print(f"K2 Reads column: {'SHOWN' if show_k2_column else 'HIDDEN'}")
    print(f"Subkey grouping: {'ENABLED' if use_subkey else 'DISABLED'}")
    print(f"Strain detail table: {'INLINE' if show_strains_table else 'APPENDIX'}")
    print(f"\nFiltering settings:")
    print(f"  Show Potentials: {args.show_potentials}")
    print(f"  Show Commensals: {args.show_commensals}")
    print(f"  Show Opportunistic: {args.show_opportunistic}")
    print(f"  Show Unidentified: {args.show_unidentified}")
    print(f"  Sorting mode: {'Alphabetical' if args.sort_alphabetical else 'TASS Score (descending)'}")
    print(f"  Max members per group: {args.max_members if args.max_members else 'Unlimited'}")
    print(f"  Max TOC groups per sample: {args.max_toc}")
    # Collect low confidence strains
    low_confidence_strains = []
    for sample_name, species_groups in samples_dict.items():
        for species_group in species_groups:
            for strain in species_group.get('members', []):
                if should_include_strain(strain, args) and has_min_reads(strain, args.min_reads):
                    if not passes_confidence_threshold(strain, args.min_conf):
                        if strain.get("high_cons"):
                            low_confidence_strains.append({
                                'strain': strain,
                                'sample_name': sample_name,
                                'species_group': species_group
                            })

    valid_bookmarks = collect_all_bookmarks(samples_dict, low_confidence_strains)
    valid_bookmarks.update(taxid_to_bookmark.values())

    # ── Title ─────────────────────────────────────────────────────────────────
    story.append(Paragraph("Organism Discovery Report", title_style))
    story.append(Spacer(1, 0.05*inch))

    current_time = datetime.now().strftime("%Y-%m-%d %H:%M")
    generation_text = (
        f"This report was generated using "
        f"<link href=\"https://github.com/jhuapl-bio/taxtriage\" color=\"blue\">TaxTriage</link> "
        f"<b>{args.version}</b> on <b>{current_time}</b> and is derived from "
        f"an in <link href=\"https://github.com/jhuapl-bio/taxtriage/blob/main/assets/pathogen_sheet.csv\" "
        f"color=\"blue\">development spreadsheet of human-host pathogens</link>. "
        f"Samples with confidence score below {args.min_conf} were filtered out."
    )
    story.append(Paragraph(generation_text, metadata_style))
    story.append(Spacer(1, 0.02*inch))
    story.append(Paragraph("<b>★</b> = High Consequence Pathogen", small_style))
    story.append(Spacer(1, 0.02*inch))

    # ── Table of Contents ─────────────────────────────────────────────────────
    story.append(Paragraph("Table of Contents", heading_style))
    story.append(Spacer(1, 0.03*inch))
    story.append(Paragraph(
        "<i>Format: Sample Name (Total Alignments # - # Primary Pathogens) → "
        "Category Label ■ (# Primary Strains, Max TASS)</i>", small_style))
    story.append(Spacer(1, 0.15*inch))
    story.append(Paragraph(
        "Click on sample names or species groups to jump to their sections. "
        "Only samples/groups with visible strains are shown here", small_style))
    story.append(Paragraph(
        "The table is organized by samples first, then in order of TASS Score by default or alphabetical if selected. "
        "Each row below a group is attributed to the highest-TASS strain in that group, and the number of strains shown per group is limited in the TOC for readability (see settings).",
        small_style))

    for sample_name in sorted(samples_dict.keys()):
        bookmark_name = f"sample_{sanitize_bookmark_name(sample_name)}"
        species_groups = samples_dict[sample_name]
        total_alignments, primary_count = get_sample_stats(species_groups)

        if args.sort_alphabetical:
            sorted_groups = sorted(species_groups, key=lambda sg: sg.get('toplevelname', 'Unknown'))
        else:
            sorted_groups = sorted(species_groups, key=lambda sg: get_species_group_stats(sg)[0], reverse=True)

        visible_groups = []
        for sg in sorted_groups:
            has_visible = any(
                should_include_strain(s, args)
                and has_min_reads(s, 1)
                and passes_confidence_threshold(s, args.min_conf)
                for s in sg.get('members', [])
            )
            if has_visible:
                visible_groups.append(sg)

        if not visible_groups:
            continue

        link_text = create_safe_link(
            f'{sample_name} ({total_alignments:,} Alignments - {primary_count} Primary Pathogens)',
            bookmark_name, valid_bookmarks)
        story.append(Paragraph(link_text, heading_style))
        story.append(Spacer(1, 0.04*inch))

        toc_groups = visible_groups[:args.max_toc]
        has_more = len(visible_groups) > args.max_toc

        species_links = []
        for sg in toc_groups:
            sp_name = sg.get('toplevelname', 'Unknown')
            sp_key = sg.get('toplevelkey', sg.get('key', 'unknown'))
            sp_bm = f"species_{sanitize_bookmark_name(sample_name)}_{sp_key}"
            g_tass, g_primary, _ = get_species_group_stats(sg)
            g_cat = sg.get('microbial_category', 'Unknown')
            g_ann = sg.get('annClass', '')
            cat_color = get_category_color(g_cat, g_ann, alpha=1.0)
            color_hex = '#{:02x}{:02x}{:02x}'.format(
                int(cat_color.red * 255), int(cat_color.green * 255), int(cat_color.blue * 255))
            species_links.append(
                f'{create_safe_link(sp_name, sp_bm, valid_bookmarks)} '
                f'<font color="{color_hex}">■</font> '
                f'({g_primary}, {g_tass*100:.1f})'
            )

        if has_more:
            first_hidden = visible_groups[args.max_toc]
            h_key = first_hidden.get('toplevelkey', first_hidden.get('key', 'unknown'))
            h_bm = f"species_{sanitize_bookmark_name(sample_name)}_{h_key}"
            species_links.append(create_safe_link(
                f"... ({len(visible_groups) - args.max_toc} more)", h_bm, valid_bookmarks))

        story.append(Paragraph(f"→ {', '.join(species_links)}", indent_style))
        story.append(Spacer(1, 0.02*inch))

    story.append(Spacer(1, 0.05*inch))
    story.append(Paragraph('• ' + create_safe_link('Color Key', 'color_key', valid_bookmarks), styles['Normal']))
    story.append(Paragraph('• ' + create_safe_link('Column Explanations', 'column_explanations', valid_bookmarks), styles['Normal']))
    if low_confidence_strains:
        story.append(Paragraph(
            '• ' + create_safe_link('Low Confidence, High Consequence Detections', 'low_confidence', valid_bookmarks),
            styles['Normal']))
    if not show_strains_table:
        story.append(Paragraph(
            '• ' + create_safe_link('Additional Lower-Level Strain Detail', 'strain_detail_table', valid_bookmarks),
            styles['Normal']))
    story.append(Paragraph('• ' + create_safe_link('Additional Information', 'additional_info', valid_bookmarks), styles['Normal']))
    story.append(Spacer(1, 0.00*inch))

    # ── Per-sample content ────────────────────────────────────────────────────
    # Also build sorted_groups_by_sample for the strain appendix.
    sorted_groups_by_sample = []  # list of (sample_name, sorted_groups)

    for sample_name in sorted(samples_dict.keys()):
        bookmark_name = f"sample_{sanitize_bookmark_name(sample_name)}"
        story.append(AnchorFlowable(bookmark_name))

        sampletype = (samples_dict[sample_name][0].get('sampletype', 'Unspecified Type')
                      if samples_dict[sample_name] else 'Unspecified Type')
        story.append(Paragraph(f"Sample: {sample_name} ({sampletype})", heading_style))
        story.append(Spacer(1, 0.1*inch))

        species_groups = samples_dict[sample_name]
        sample_total_reads = sum(sg.get('numreads', 0) for sg in species_groups)

        if args.sort_alphabetical:
            sorted_groups = sorted(species_groups, key=lambda sg: sg.get('toplevelname', 'Unknown'))
        else:
            sorted_groups = sorted(species_groups, key=lambda sg: get_species_group_stats(sg)[0], reverse=True)

        # Record for appendix
        sorted_groups_by_sample.append((sample_name, sorted_groups))

        for sg in sorted_groups:
            sp_key = sg.get('toplevelkey', sg.get('key', 'unknown'))
            sp_bm = f"species_{sanitize_bookmark_name(sample_name)}_{sp_key}"
            story.append(AnchorFlowable(sp_bm))

        all_sample_strains = []
        species_group_map = {}

        for sg in sorted_groups:
            group_strains = []
            for strain in sg.get('members', []):
                if should_include_strain(strain, args) and has_min_reads(strain, 1):
                    if passes_confidence_threshold(strain, args.min_conf):
                        group_strains.append(strain)
                        species_group_map[id(strain)] = sg
            group_strains.sort(key=lambda s: s.get('tass_score', 0), reverse=True)
            if args.max_members is not None and args.max_members > 0:
                group_strains = group_strains[:args.max_members]
            all_sample_strains.extend(group_strains)

        if all_sample_strains:
            combined_table = create_combined_sample_table(
                all_sample_strains, species_group_map, small_style,
                show_ani_column, show_k2_column,
                taxid_to_bookmark, valid_bookmarks,
                sample_total_reads=sample_total_reads,
                sample_name=sample_name,
                available_width=available_width,
                use_subkey=use_subkey,
                show_strains_table=show_strains_table,
            )
            story.append(combined_table)
        else:
            story.append(Paragraph(
                "<i>No data available for this sample above confidence threshold</i>",
                styles['Italic']))

        if args.max_members is not None and args.max_members > 0:
            story.append(Paragraph(
                f"<i>Showing top {args.max_members} strains per group by TASS score</i>",
                small_style))
        story.append(Spacer(1, 0.1*inch))

    # ── Footer: Color Key ─────────────────────────────────────────────────────
    story.append(Spacer(1, 0.15*inch))
    story.append(AnchorFlowable('color_key'))
    story.append(Paragraph("<b>Color Key:</b>", heading_style))
    url = 'https://github.com/jhuapl-bio/taxtriage/blob/main/assets/pathogen_sheet.csv'
    legend_data = [
        ['', Paragraph('Category', legend_text_style), Paragraph('Description', legend_text_style)],
        ['', Paragraph('Primary Pathogen (Direct)', legend_text_style),
         Paragraph(f'The organism is directly detected according to <link href="{url}" color="blue">taxonomic id</link> or assembly name and is listed as being of importance in your sample type.', legend_text_style)],
        ['', Paragraph('Primary Pathogen', legend_text_style),
         Paragraph('It is directly detected according to taxonomic id or assembly name and is listed as being of importance in a different sample type.', legend_text_style)],
        ['', Paragraph('Commensal', legend_text_style), Paragraph('Normal flora / non-pathogenic', legend_text_style)],
        ['', Paragraph('Opportunistic', legend_text_style), Paragraph('May cause disease in certain conditions', legend_text_style)],
        ['', Paragraph('Potential', legend_text_style), Paragraph('Potential pathogen requiring further investigation', legend_text_style)],
        ['', Paragraph('Unknown', legend_text_style), Paragraph('Classified organism with unknown significance', legend_text_style)],
    ]
    legend_table = Table(legend_data, colWidths=[0.3*inch, 1.6*inch, 3.3*inch])
    legend_table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#3498DB')),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, 0), 9),
        ('BOTTOMPADDING', (0, 0), (-1, 0), 8),
        ('BACKGROUND', (0, 1), (0, 1), colors.lightcoral),
        ('BACKGROUND', (0, 2), (0, 2), colors.HexColor('#fab462')),
        ('BACKGROUND', (0, 3), (0, 3), colors.lightgreen),
        ('BACKGROUND', (0, 4), (0, 4), colors.HexColor('#ffe6a8')),
        ('BACKGROUND', (0, 5), (0, 5), colors.lightblue),
        ('BACKGROUND', (0, 6), (0, 6), colors.white),
        ('GRID', (0, 0), (-1, -1), 1, colors.black),
        ('FONTSIZE', (1, 1), (-1, -1), 8),
        ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
        ('LEFTPADDING', (1, 1), (-1, -1), 6),
    ]))
    story.append(legend_table)
    story.append(Spacer(1, 0.02*inch))

    # ── Footer: Column Explanations ───────────────────────────────────────────
    story.append(AnchorFlowable('column_explanations'))
    story.append(Paragraph("<b>Column Explanations:</b>", heading_style))
    story.append(Spacer(1, 0.03*inch))
    explanations = [
        "• <b>Specimen ID (Taxonomic ID #):</b> The unique identifier for the sample including its name and taxonomic ID. The taxonomic ID is a link to the NCBI Taxonomy Browser for that organism.",
        "• <b>Detected Organism:</b> The organism detected in the sample, which could be a bacterium, virus, fungus, or parasite.",
        "• <b>TASS Score:</b> A metric between 0 and 100 that reflects the confidence of the organism's detection, with 100 being the highest value.",
        "• <b>K2 Reads:</b> The number of reads classified by Kraken2, a tool for taxonomic classification of sequencing data.",
        "• <b># Reads Aligned:</b> The number of reads from the sequencing data that align to the organism's genome. (%) refers to all alignments for that species across the entire sample.",
        "• <b>RPM:</b> Reads Per Million (RPM). This normalized metric allows for comparison of abundance across samples and organisms of different sizes.",
        "• <b>Cov.:</b> The coverage of the organism's genome by aligned reads, expressed as a percentage.",
        "• <b>Strain Detail Columns (subkey mode):</b> When subkey grouping is enabled (default), the rightmost section shows each individual strain within the subkey group — its name, TASS score, and read count. The main row values are taken from the highest-confidence (best-TASS) member of that subkey group. Pass <b>--no_subkey</b> to disable this and revert to one row per strain.",
        "• <b>High ANI Matches:</b> When the High ANI column is shown, this column lists any strains in the dataset that have an Average Nucleotide Identity (ANI) above the specified threshold with the given strain. Each match includes a link to its section in the report if it is present, or just the taxonomic ID and ANI percentage if not present.",
    ]
    for e in explanations:
        story.append(Paragraph(e, metadata_style))
        story.append(Spacer(1, 0.02*inch))
    story.append(Spacer(1, 0.1*inch))

    # ── Footer: Low Confidence ────────────────────────────────────────────────
    if low_confidence_strains:
        story.append(AnchorFlowable('low_confidence'))
        story.append(Paragraph("<b>Low Confidence, High Consequence Detections:</b>", heading_style))
        story.append(Paragraph(
            f"The following strains were detected but fell below the confidence threshold "
            f"of {args.min_conf} and are listed here for reference.", metadata_style))
        story.append(Spacer(1, 0.05*inch))
        story.append(create_low_confidence_table(
            low_confidence_strains, small_style, show_k2_column, available_width))
        story.append(Spacer(1, 0.1*inch))

    # ── Strain Detail Appendix (when inline table is off) ────────────────────
    if not show_strains_table:
        appendix_items = create_strain_detail_tables(
            samples_dict, sorted_groups_by_sample,
            small_style, show_k2_column, valid_bookmarks,
            args, available_width,
        )
        story.extend(appendix_items)

    # ── Footer: Additional Information ───────────────────────────────────────
    story.append(AnchorFlowable('additional_info'))
    story.append(Paragraph("<b>Additional Information:</b>", heading_style))
    story.append(Spacer(1, 0.01*inch))
    url2 = "https://github.com/jhuapl-bio/taxtriage/blob/main/docs/usage.md#confidence-scoring"
    story.append(Paragraph(
        f'Please visit our <a href="{url2}"><b><font color="blue">DOCUMENTATION PAGE</font></b></a> '
        f'for more information on how confidence is calculated.', metadata_style))
    story.append(Spacer(1, 0.01*inch))
    story.append(Paragraph(
        "The following information highlights the description for the color combinations "
        "for each organism class in the annotated table(s).", metadata_style))
    story.append(Spacer(1, 0.01*inch))
    story.append(Paragraph(
        "Please see the relevant Discovery Analysis txt file for low confidence, "
        "high consequence annotations that were not present in the pdf.", metadata_style))
    story.append(Spacer(1, 0.01*inch))
    story.append(Paragraph(
        "Read amounts are represented as the <b>total number of aligned reads</b> of sufficient "
        "mapping quality <b>(% aligned for all reads in sample)</b>.", metadata_style))
    story.append(Spacer(1, 0.01*inch))
    story.append(Paragraph(
        "If there are questions or issues with your report, please open an issue on GitHub as a "
        "discussion <link href=\"https://github.com/jhuapl-bio/taxtriage/discussions\" color=\"blue\">"
        "here</link>. Issues should be tracked/submitted at "
        "<link href=\"https://github.com/jhuapl-bio/taxtriage/issues\" color=\"blue\">this link</link>.",
        metadata_style))

    doc.build(story)
    print(f"\nPDF created successfully: {output_path}")


def create_tabular_output(output_path, samples_dict, args):
    """
    Create a tabular output file (CSV/TSV/TXT/XLSX) with strain-level data.
    Includes ALL strains (not filtered). A 'Subkey' column is included.
    """
    file_ext = os.path.splitext(output_path)[1].lower()

    headers = [
        'Index', 'index', 'Detected Organism', 'Specimen ID', 'Sample Type',
        '% Reads', '# Reads Aligned', '% Aligned Reads', 'Coverage',
        'HHS Percentile', 'IsAnnotated', 'AnnClass', 'Microbial Category',
        'High Consequence', 'Taxonomic ID #', 'Status', 'Gini Coefficient',
        'Mean BaseQ', 'Mean MapQ', 'Mean Depth', 'isSpecies',
        'Pathogenic Subsp/Strains', 'K2 Reads', 'RPKM', 'RPM',
        'Parent K2 Reads', 'MapQ Score', 'Disparity Score', 'Minhash Score',
        'Diamond Identity', 'K2 Disparity Score', 'Siblings score',
        'Breadth Weight Score', 'TASS Score', 'MicrobeRT Probability',
        'MicrobeRT Model', 'Reads Aligned', 'Group', 'Subkey'
    ]

    all_rows = []
    global_index = 0

    for sample_name in sorted(samples_dict.keys()):
        species_groups = samples_dict[sample_name]

        if args.sort_alphabetical:
            sorted_groups = sorted(species_groups, key=lambda sg: sg.get('toplevelname', 'Unknown'))
        else:
            sorted_groups = sorted(
                species_groups, key=lambda sg: get_species_group_stats(sg)[0], reverse=True)

        sample_total_reads = max(1, sum(sg.get('numreads', 0) for sg in species_groups))

        for sg in sorted_groups:
            group_key = sg.get('toplevelkey', sg.get('key', ''))
            sample_type = sg.get('sampletype', 'unknown')
            strains = sorted(sg.get('members', []),
                             key=lambda s: s.get('tass_score', 0), reverse=True)

            for local_idx, strain in enumerate(strains):
                if not has_min_reads(strain, 1):
                    continue
                strain_reads = float(strain.get('numreads', 0) or 0)
                pct_reads = strain_reads / sample_total_reads * 100.0

                all_rows.append([
                    global_index, local_idx,
                    strain.get('name', 'Unknown'), sample_name, sample_type,
                    f"{pct_reads:.4f}", int(strain_reads), f"{pct_reads:.4f}",
                    f"{(min(1, strain.get('coverage', 0)) or 0)*100:.0f}%",
                    '100.0',
                    'Yes' if strain.get('isAnnotated', True) else 'No',
                    strain.get('annClass', ''),
                    strain.get('microbial_category', 'Unknown'),
                    'True' if strain.get('high_cons', False) else 'False',
                    strain.get('key', ''), strain.get('status', ''),
                    f"{(strain.get('gini_coefficient', 0) or 0):.2f}",
                    f"{(strain.get('meanbaseq', 0) or 0):.2f}",
                    f"{(strain.get('meanmapq', 0) or 0):.2f}",
                    f"{(strain.get('meandepth', 0) or 0):.1f}",
                    'True' if strain.get('isSpecies', False) else 'False',
                    '',
                    int(strain.get('k2_reads', 0) or 0),
                    strain.get("rpkm", 0) or 0,
                    strain.get("rpm", 0) or 0,
                    int(strain.get('parent_k2_reads', 0) or 0),
                    f"{(strain.get('mapq_score', 0) or 0):.2f}",
                    f"{(strain.get('disparity', 0) or 0):.2f}",
                    f"{(strain.get('minhash_reduction', 0) or 0):.2f}",
                    f"{(strain.get('diamond_identity', 0) or 0):.1f}",
                    f"{(strain.get('k2_disparity_score', 0) or 0):.1f}",
                    f"{(strain.get('siblings_score', 0) or 0):.1f}",
                    f"{(strain.get('breadth_log_score', 0) or 0):.2f}",
                    int((strain.get('tass_score', 0) or 0) * 100),
                    f"{(strain.get('mmbert', 0) or 0):.4f}",
                    strain.get('mmbert_model', '') or '',
                    int(strain_reads),
                    group_key,
                    strain.get('subkey', strain.get('key', ''))  # NEW: subkey column
                ])
                global_index += 1

    df = pd.DataFrame(all_rows, columns=headers)

    if file_ext == '.csv':
        df.to_csv(output_path, index=False)
        print(f"CSV output created: {output_path}")
    elif file_ext in ['.tsv', '.txt']:
        df.to_csv(output_path, sep='\t', index=False)
        print(f"TSV output created: {output_path}")
    elif file_ext in ['.xlsx', '.xls']:
        df.to_excel(output_path, index=False, engine='openpyxl')
        print(f"Excel output created: {output_path}")
    else:
        df.to_csv(output_path, sep='\t', index=False)
        print(f"TSV output created (default): {output_path}")

    print(f"  Total strains in output: {len(all_rows)}")
    print(f"  Note: Output includes ALL strains (not filtered by category or confidence)")


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate PDF report from pathogen discovery data",
        epilog="Example: python create_report.py -i confidences.json -o report.pdf",
    )
    parser.add_argument(
        "-i", "--input", metavar="INPUT", required=True, nargs="+", default=[],
        help="Base pathogen discovery table file(s), JSON format only. Can specify more than one",
    )
    parser.add_argument(
        "-d", "--distributions", metavar="DISTRIBUTIONS", required=False,
        help="TSV file that contains all the distribution information for body sites and organisms",
    )
    parser.add_argument(
        "-a", "--abundance_col", metavar="ABU", required=False, default='% Aligned Reads',
        help="Name of abundance column, default is abundance",
    )
    parser.add_argument(
        "-r", "--min_reads", metavar="READS", required=False, default=1, type=int,
        help="Minimum number of reads required to consider an organism for reporting. Default is 1.",
    )
    parser.add_argument(
        "-c", "--min_conf", metavar="MINCONF", required=False, default=0.3, type=float,
        help="Value that must be met for a table to report an organism due to confidence column.",
    )
    parser.add_argument(
        "-x", "--id_col", metavar="IDCOL", required=False, default="Detected Organism",
        help="Name of id column, default is id",
    )
    parser.add_argument(
        "-p", "--percentile", metavar="PERCENTILE", required=False, type=float, default=0.75,
        help="Only show organisms that are in the top percentile of healthy subjects expected abu",
    )
    parser.add_argument(
        "-v", "--version", metavar="VERSION", required=False, default='Local Build',
        help="What version of TaxTriage is in use",
    )
    parser.add_argument('--sorttass', action="store_true", required=False,
                        help="Sort by TASS score if available")
    parser.add_argument('--sort_alphabetical', action="store_true", required=False,
                        help="Sort groups alphabetically instead of by TASS score")
    parser.add_argument('--enable_matrix', action="store_true", required=False,
                        help="Enable matrix view if available")
    parser.add_argument(
        "--max_members", metavar="MAX_MEMBERS", required=False, type=int, default=None,
        help="Maximum number of top strains (by TASS) to show per group. Default: show all",
    )
    parser.add_argument(
        "--max_toc", metavar="MAX_TOC", required=False, type=int, default=4,
        help="Maximum number of species groups to show in TOC per sample. Default: 4",
    )
    parser.add_argument("--show_commensals", action="store_true", required=False,
                        help="Show the commensals table")
    parser.add_argument("--show_unidentified", action="store_true", required=False,
                        help="Show the all organisms now listed as commensal or pathogen")
    parser.add_argument("--show_potentials", action="store_true", required=False,
                        help="Show the potentials table")
    parser.add_argument("--show_opportunistic", action="store_true", required=False,
                        help="Show the opportunistic table")
    parser.add_argument(
        "-m", "--missing_samples", metavar="MISSING", required=False, default=None, nargs="+",
        help="Missing samples if any",
    )
    parser.add_argument(
        "-s", "--sitecol", metavar="SCOL", required=False, default='Sample Type',
        help="Name of site column, default is body_site",
    )
    parser.add_argument(
        "-t", "--type", metavar="TYPE", required=False, default='Detected Organism',
        help="What type of data is being processed. Options: 'Taxonomic ID #' or 'Detected Organism'.",
        choices=['Taxonomic ID #', 'Detected Organism'],
    )
    parser.add_argument(
        "--taxdump", metavar="TAXDUMP", required=False, default=None,
        help="Merge the entries on a specific rank args.rank, importing files from nodes.dmp, names.dmp and potentially merged.dmp",
    )
    parser.add_argument(
        "--rank", metavar="RANK", required=False, default="genus",
        help='IF merging with taxdump, what rank to merge on',
    )
    parser.add_argument(
        "-o", "--output", metavar="OUTPUT", required=True, type=str,
        help="Path of output file (pdf)",
    )
    parser.add_argument(
        "-u", "--output_txt", metavar="OUTPUT_TXT", required=False, type=str,
        help="Path of tabular output file. Format determined by extension: .csv, .tsv, .txt (TSV), or .xlsx",
    )
    # ── NEW: subkey grouping control ──────────────────────────────────────────
    parser.add_argument(
        "--no_subkey",
        action="store_true",
        default=False,
        help=(
            "Disable subkey grouping (default: enabled). "
            "By default, members sharing the same 'subkey' value are collapsed into a "
            "single row. The row's left columns show the best-TASS member's metrics "
            "(TASS, K2 Reads, Reads, RPM, Coverage); the right 3 columns contain a "
            "nested mini-table listing every individual strain with its name, TASS score, "
            "and read count. Pass --no_subkey to revert to the original flat "
            "one-row-per-strain layout."
        ),
    )
    # ── Strain detail table toggle ────────────────────────────────────────────
    parser.add_argument(
        "--show_strains_table",
        action="store_true",
        default=False,
        help=(
            "Show the full strain detail mini-table inline in the right column of the "
            "main report table (default: off). When off (default), strain details are "
            "moved to a separate appendix table placed after the Low Confidence section, "
            "and the main table right column shows a compact strain-count link instead."
        ),
    )

    return parser.parse_args()


def main():
    args = parse_args()

    taxdump_dict = {}
    names_map = {}
    merged_tax_data = {}

    if args.taxdump:
        if os.path.exists(f"{args.taxdump}/nodes.dmp"):
            taxdump_dict = load_taxdump(f"{args.taxdump}/nodes.dmp")
            print(f"Loaded nodes.dmp: {len(taxdump_dict)} entries")
        if os.path.exists(f"{args.taxdump}/names.dmp"):
            names_map = load_names(f"{args.taxdump}/names.dmp")
            print(f"Loaded names.dmp: {len(names_map)} entries")
        if os.path.exists(f"{args.taxdump}/merged.dmp"):
            merged_tax_data = load_merged(f"{args.taxdump}/merged.dmp")
            print(f"Loaded merged.dmp: {len(merged_tax_data)} entries")

    sample_data = load_json_samples(args.input)
    print(f"Loaded {len(sample_data)} species groups from JSON file(s)")

    samples_dict = organize_data_by_sample(sample_data)
    print(f"Found {len(samples_dict)} unique sample(s)")

    print(f"\nConfiguration:")
    print(f"  Output PDF: {args.output}")
    if args.output_txt:
        print(f"  Output TXT: {args.output_txt}")
    print(f"  Min Confidence: {args.min_conf}")
    print(f"  Show Potentials: {args.show_potentials}")
    print(f"  Show Unidentified: {args.show_unidentified}")
    print(f"  Show Commensals: {args.show_commensals}")
    print(f"  Subkey grouping: {'DISABLED' if args.no_subkey else 'ENABLED'}")
    print(f"  High ANI matches: read from 'high_ani_matches' field in JSON (set by match_paths.py)")

    create_pdf_template(args.output, samples_dict, args)

    if args.output_txt:
        create_tabular_output(args.output_txt, samples_dict, args)

    return taxdump_dict, names_map, merged_tax_data, sample_data


if __name__ == "__main__":
    taxdump_dict, names_map, merged_tax_data, sample_data = main()
