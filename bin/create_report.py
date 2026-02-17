#!/usr/bin/env python3

import json
import argparse
import os
import pandas as pd
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer, PageBreak
from reportlab.platypus.flowables import AnchorFlowable
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_LEFT, TA_CENTER
from reportlab.lib import colors
from reportlab.lib.units import inch
from datetime import datetime
from map_taxid import load_taxdump, load_names, load_merged


def load_ani_matrix(ani_matrix_files):
    import pandas as pd
    if not ani_matrix_files:
        return {}
    ani_data = {}
    for ani_file in ani_matrix_files:
        if not os.path.exists(ani_file):
            print(f"Warning: ANI matrix file not found: {ani_file}")
            continue
        df = pd.read_csv(ani_file, index_col=0)
        print(f"Loaded ANI matrix from {ani_file}: {df.shape[0]} x {df.shape[1]}")
        for taxid1 in df.index:
            if str(taxid1) not in ani_data:
                ani_data[str(taxid1)] = {}
            for taxid2 in df.columns:
                ani_value = df.loc[taxid1, taxid2]
                if pd.notna(ani_value):
                    ani_data[str(taxid1)][str(taxid2)] = float(ani_value)
    return ani_data


def get_high_ani_matches(strain_key, ani_data, ani_threshold, all_strain_keys):
    matches = []
    strain_key_str = str(strain_key)
    if strain_key_str not in ani_data:
        return matches
    for other_key, ani_value in ani_data[strain_key_str].items():
        if other_key == strain_key_str:
            continue
        if ani_value >= ani_threshold:
            matches.append((other_key, ani_value))
    matches.sort(key=lambda x: x[1], reverse=True)
    return matches


def check_if_any_high_ani_in_dataset(strains, ani_data, ani_threshold):
    for strain in strains:
        strain_key = strain.get('key', '')
        matches = get_high_ani_matches(strain_key, ani_data, ani_threshold, set())
        if matches:
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
    for sample_name, species_groups in samples_dict.items():
        sample_bookmark = f"sample_{sanitize_bookmark_name(sample_name)}"
        bookmarks.add(sample_bookmark)
        for species_group in species_groups:
            species_key = species_group.get('toplevelkey', species_group.get('key', 'unknown'))
            species_bookmark = f"species_{sanitize_bookmark_name(sample_name)}_{species_key}"
            bookmarks.add(species_bookmark)
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


def create_combined_sample_table(all_strains, species_group_map, small_style, ani_data,
                                  ani_threshold, show_ani_column, show_k2_column,
                                  taxid_to_bookmark, valid_bookmarks,
                                  sample_total_reads=0, sample_name=None,
                                  available_width=None, use_subkey=True):
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
        'StrainName', parent=small_style, fontSize=10, leading=10)
    data_style = ParagraphStyle(
        'DataStyle', parent=small_style, fontSize=8, leading=9)
    ani_style = ParagraphStyle(
        'ANIStyle', parent=small_style, fontSize=6, leading=8)
    group_header_style = ParagraphStyle(
        'GroupHeader', parent=small_style, fontSize=10, leading=11,
        fontName='Helvetica-Bold')
    mini_style = ParagraphStyle(
        'MiniStyle', parent=small_style, fontSize=7, leading=9)
    mini_header_style = ParagraphStyle(
        'MiniHeaderStyle', parent=small_style, fontSize=7, leading=9,
        fontName='Helvetica-Bold')

    if available_width is None:
        available_width = 8.5 * inch - 0.02 * 8.5 * inch

    # ── Headers ──────────────────────────────────────────────────────────────
    base_headers = ['', 'Organism', 'TASS']
    if show_k2_column:
        base_headers.append('K2 Reads')
    base_headers += ['Reads', 'RPM', 'Coverage']
    if show_ani_column:
        base_headers.append('High ANI')
    n_base = len(base_headers)

    # Strain detail column in subkey mode (single column holding a nested mini-table)
    if use_subkey:
        headers = base_headers + ['']
    else:
        headers = base_headers
    n_total = len(headers)

    # ── Column widths ─────────────────────────────────────────────────────────
    def _base_col_widths(w):
        if show_k2_column and show_ani_column:
            return [w*0.03, w*0.28, w*0.08, w*0.10, w*0.12, w*0.10, w*0.09, w*0.20]
        elif show_k2_column:
            return [w*0.03, w*0.33, w*0.09, w*0.13, w*0.14, w*0.12, w*0.16]
        elif show_ani_column:
            return [w*0.03, w*0.32, w*0.09, w*0.14, w*0.11, w*0.11, w*0.20]
        else:
            return [w*0.03, w*0.42, w*0.10, w*0.17, w*0.13, w*0.15]

    if use_subkey:
        left_w = available_width * 0.60
        right_w = available_width * 0.40
        col_widths = _base_col_widths(left_w) + [right_w]
    else:
        col_widths = _base_col_widths(available_width)

    # ── Table data ────────────────────────────────────────────────────────────
    table_data = [headers]
    table_styles = []
    group_row_indices = []

    current_species_key = None
    row_idx = 1  # row 0 is the header
    emitted_subkeys_per_group = {}  # {species_key: set of already-emitted subkeys}

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

            group_name_para = Paragraph(
                f'<b><link href="https://www.ncbi.nlm.nih.gov/taxonomy/?term={species_key}" '
                f'color="blue">{group_name}</link></b>',
                group_header_style
            )

            group_row = ['', group_name_para, '']
            if show_k2_column:
                group_row.append(Paragraph(f'<b>{group_k2_reads:,.0f}</b>', group_header_style))
            group_row.append(Paragraph(f'<b>{group_reads:,.0f}</b>', group_header_style))
            group_row += ['', '']
            if show_ani_column:
                group_row.append('')
            if use_subkey:
                group_row += ['']

            table_data.append(group_row)
            group_row_indices.append(row_idx)

            # Span group name across Organism + TASS columns only
            table_styles.append(('SPAN', (1, row_idx), (2, row_idx)))
            # Span the empty trailing columns (RPM onward) together
            rpm_col_idx = 5 if show_k2_column else 4
            if rpm_col_idx < n_total:
                table_styles.append(('SPAN', (rpm_col_idx, row_idx), (n_total - 1, row_idx)))
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
            # In subkey mode, show the species-level subkeyname when key != subkey
            if str(best.get('subkey', best_key)) != str(best_key):
                display_name = best.get('subkeyname', best.get('name', 'Unknown'))
                display_key = str(best.get('subkey', best_key))
            else:
                display_name = best.get('name', 'Unknown')
                display_key = best_key
            name_html = (
                f'{display_name} '
                f'(<link href="https://www.ncbi.nlm.nih.gov/taxonomy/?term={display_key}" '
                f'color="blue">{display_key}</link>)'
            )

            strain_reads = float(best.get('numreads', 0) or 0)
            pct = (strain_reads / sample_total_reads * 100.0) if sample_total_reads else 0.0
            rpm = best.get('rpm', 0) or 0

            high_ani_text = ''
            if show_ani_column:
                matches = get_high_ani_matches(best_key, ani_data, ani_threshold, set())
                if matches:
                    ani_links = []
                    for taxid, ani_value in matches[:3]:
                        ani_pct = ani_value * 100
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

            # Right side: always render a nested mini-table for all subkey members
            mini_rows = [[
                Paragraph('<b>Strain</b>', mini_header_style),
                Paragraph('<b>TASS</b>', mini_header_style),
                Paragraph('<b>Reads</b>', mini_header_style),
            ]]
            for m in subkey_members:
                m_reads = float(m.get('numreads', 0) or 0)
                m_pct = (m_reads / sample_total_reads * 100.0) if sample_total_reads else 0.0
                m_key = m.get('key', '')
                m_star = '★ ' if m.get('high_cons', False) else ''
                mini_rows.append([
                    Paragraph(
                        f'{m_star}<link href="https://www.ncbi.nlm.nih.gov/taxonomy/?term={m_key}" '
                        f'color="blue">{m.get("name", "Unknown")}</link>',
                        mini_style
                    ),
                    Paragraph(f"{m.get('tass_score', 0)*100:.1f}", mini_style),
                    Paragraph(f"{m_reads:,.0f} ({m_pct:.1f}%)", mini_style),
                ])

            mini_tbl = Table(mini_rows)
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
                matches = get_high_ani_matches(strain_key, ani_data, ani_threshold, set())
                if matches:
                    ani_links = []
                    for taxid, ani_value in matches[:3]:
                        ani_pct = ani_value * 100
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
        ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
        ('ALIGN', (0, 1), (0, -1), 'CENTER'),
        ('VALIGN', (0, 1), (0, -1), 'MIDDLE'),
        ('FONTSIZE', (0, 1), (0, -1), 14),
        ('GRID', (0, 0), (-1, -1), 1, colors.black),
        ('VALIGN', (1, 1), (1, -1), 'MIDDLE'),
        ('VALIGN', (2, 1), (-1, -1), 'MIDDLE'),
        ('ALIGN', (2, 1), (-1, -1), 'CENTER'),
        ('ALIGN', (1, 1), (1, -1), 'LEFT'),
        ('LEFTPADDING', (1, 1), (1, -1), 6),
        ('RIGHTPADDING', (1, 1), (1, -1), 6),
        ('TOPPADDING', (1, 1), (1, -1), 6),
        ('BOTTOMPADDING', (1, 1), (1, -1), 6),
    ]

    if use_subkey:
        base_ts += [
            ('ALIGN', (n_base, 1), (n_total - 1, -1), 'LEFT'),
            ('VALIGN', (n_base, 1), (n_total - 1, -1), 'TOP'),
            ('LEFTPADDING', (n_base, 1), (n_total - 1, -1), 4),
            ('TOPPADDING', (n_base, 1), (n_total - 1, -1), 2),
            ('RIGHTPADDING', (n_base, 1), (n_total - 1, -1), 2),
            ('BOTTOMPADDING', (n_base, 1), (n_total - 1, -1), 2),
            # Vertical delimiter between species columns and strain detail columns
            ('LINEAFTER', (n_base - 1, 0), (n_base - 1, -1), 2, colors.HexColor('#3498DB')),
        ]

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


def create_pdf_template(output_path, samples_dict, ani_data, args):
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

    show_ani_column = check_if_any_high_ani_in_dataset(all_strains, ani_data, args.ani_threshold)
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

    print(f"\nHigh ANI column: {'SHOWN' if show_ani_column else 'HIDDEN'} (threshold: {args.ani_threshold})")
    print(f"K2 Reads column: {'SHOWN' if show_k2_column else 'HIDDEN'}")
    print(f"Subkey grouping: {'ENABLED' if use_subkey else 'DISABLED'}")
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
    story.append(Paragraph('• ' + create_safe_link('Additional Information', 'additional_info', valid_bookmarks), styles['Normal']))
    story.append(Spacer(1, 0.00*inch))

    # ── Per-sample content ────────────────────────────────────────────────────
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
                all_sample_strains, species_group_map, small_style, ani_data,
                args.ani_threshold, show_ani_column, show_k2_column,
                taxid_to_bookmark, valid_bookmarks,
                sample_total_reads=sample_total_reads,
                sample_name=sample_name,
                available_width=available_width,
                use_subkey=use_subkey
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
        "• <b>Specimen ID (Type):</b> The unique identifier for the sample, including the type of specimen (e.g., blood, tissue).",
        "• <b>Detected Organism:</b> The organism detected in the sample, which could be a bacterium, virus, fungus, or parasite.",
        "• <b>Microbial Category:</b> The classification of the organism, indicating whether it is primary, opportunistic, commensal, or potential.",
        "• <b># Reads Aligned:</b> The number of reads from the sequencing data that align to the organism's genome. (%) refers to all alignments for that species across the entire sample.",
        "• <b>RPM:</b> Reads Per Million (RPM). This normalized metric allows for comparison of abundance across samples and organisms of different sizes.",
        "• <b>TASS Score:</b> A metric between 0 and 100 that reflects the confidence of the organism's detection, with 100 being the highest value.",
        "• <b>Taxonomic ID #:</b> The taxid for the organism according to NCBI Taxonomy. The parenthesis (if present) is the group it belongs to, usually the genus.",
        "• <b>K2 Reads:</b> The number of reads classified by Kraken2, a tool for taxonomic classification of sequencing data.",
        "• <b>Strain Detail Columns (subkey mode):</b> When subkey grouping is enabled (default), the rightmost section shows each individual strain within the subkey group — its name, TASS score, and read count. The main row values are taken from the highest-confidence (best-TASS) member of that subkey group. Pass <b>--no_subkey</b> to disable this and revert to one row per strain.",
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
        "--ani_threshold", metavar="ANI_THRESHOLD", required=False, type=float, default=0.90,
        help="Threshold for ANI to consider 'High ANI' in the report",
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
        "--ani_matrix", metavar="ANI_MATRIX", required=False, type=str, nargs="+",
        help="Path to organism ANI matrix CSV file from Sourmash signatures. Optional",
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

    ani_data = {}
    if args.ani_matrix:
        ani_data = load_ani_matrix(args.ani_matrix)
        print(f"Loaded ANI data for {len(ani_data)} taxa")

    sample_data = load_json_samples(args.input)
    print(f"Loaded {len(sample_data)} species groups from JSON file(s)")

    samples_dict = organize_data_by_sample(sample_data)
    print(f"Found {len(samples_dict)} unique sample(s)")

    print(f"\nConfiguration:")
    print(f"  Output PDF: {args.output}")
    if args.output_txt:
        print(f"  Output TXT: {args.output_txt}")
    if args.ani_matrix:
        print(f"  ANI Matrix: {args.ani_matrix}")
        print(f"  ANI Threshold: {args.ani_threshold}")
    print(f"  Min Confidence: {args.min_conf}")
    print(f"  Show Potentials: {args.show_potentials}")
    print(f"  Show Unidentified: {args.show_unidentified}")
    print(f"  Show Commensals: {args.show_commensals}")
    print(f"  Subkey grouping: {'DISABLED' if args.no_subkey else 'ENABLED'}")

    create_pdf_template(args.output, samples_dict, ani_data, args)

    if args.output_txt:
        create_tabular_output(args.output_txt, samples_dict, args)

    return taxdump_dict, names_map, merged_tax_data, sample_data


if __name__ == "__main__":
    taxdump_dict, names_map, merged_tax_data, sample_data = main()
