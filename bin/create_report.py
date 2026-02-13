#!/usr/bin/env python3

import json
import argparse
import os
import pandas as pd
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer, PageBreak
from reportlab.platypus.flowables import AnchorFlowable  # FIX: proper named destinations
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_LEFT, TA_CENTER
from reportlab.lib import colors
from reportlab.lib.units import inch
from datetime import datetime
from map_taxid import load_taxdump, load_names, load_merged


def load_ani_matrix(ani_matrix_files):
    """
    Load ANI matrix files and return a dictionary for quick lookups.

    Args:
        ani_matrix_files: List of CSV file paths containing ANI data

    Returns:
        dict: {taxid1: {taxid2: ani_value, ...}, ...}
    """
    import pandas as pd

    if not ani_matrix_files:
        return {}

    ani_data = {}

    for ani_file in ani_matrix_files:
        if not os.path.exists(ani_file):
            print(f"Warning: ANI matrix file not found: {ani_file}")
            continue

        # Read the CSV file
        df = pd.read_csv(ani_file, index_col=0)
        print(f"Loaded ANI matrix from {ani_file}: {df.shape[0]} x {df.shape[1]}")

        # Convert to dictionary format for easy lookup
        for taxid1 in df.index:
            if str(taxid1) not in ani_data:
                ani_data[str(taxid1)] = {}

            for taxid2 in df.columns:
                ani_value = df.loc[taxid1, taxid2]
                if pd.notna(ani_value):
                    ani_data[str(taxid1)][str(taxid2)] = float(ani_value)

    return ani_data


def get_high_ani_matches(strain_key, ani_data, ani_threshold, all_strain_keys):
    """
    Find all strains with ANI >= threshold for a given strain.

    Args:
        strain_key: The taxid/key for the strain to check
        ani_data: ANI dictionary from load_ani_matrix
        ani_threshold: Minimum ANI value to report
        all_strain_keys: Set of all strain keys in the dataset

    Returns:
        list: List of (taxid, ani_value) tuples for matches
    """
    matches = []
    strain_key_str = str(strain_key)

    if strain_key_str not in ani_data:
        return matches

    for other_key, ani_value in ani_data[strain_key_str].items():
        # Don't report self-matches
        if other_key == strain_key_str:
            continue

        # Only report if >= threshold
        if ani_value >= ani_threshold:
            matches.append((other_key, ani_value))

    # Sort by ANI value (highest first)
    matches.sort(key=lambda x: x[1], reverse=True)

    return matches


def check_if_any_high_ani_in_dataset(strains, ani_data, ani_threshold):
    """
    Check if ANY strain in the dataset has high ANI matches.

    Args:
        strains: List of all strains across all species groups
        ani_data: ANI dictionary
        ani_threshold: Minimum ANI value

    Returns:
        bool: True if at least one strain has high ANI matches
    """
    for strain in strains:
        strain_key = strain.get('key', '')
        matches = get_high_ani_matches(strain_key, ani_data, ani_threshold, set())
        if matches:
            return True

    return False


def check_if_k2_reads_present(strains):
    """
    Check if any strain has k2_reads data.

    Args:
        strains: List of all strains

    Returns:
        bool: True if at least one strain has k2_reads
    """
    for strain in strains:
        if strain.get('k2_reads') is not None and strain.get('k2_reads', 0) > 0:
            return True
    return False


def should_include_strain(strain, args):
    """
    Determine if a strain should be included based on its microbial category.

    Args:
        strain: Strain dictionary
        args: Command line arguments

    Returns:
        bool: True if strain should be included
    """
    category = str(strain.get('microbial_category', 'Unknown'))

    # Always include Primary Pathogens and Opportunistic
    if 'Primary' in category:
        return True
    if "Opportunistic" in category and args.show_opportunistic:
        return True
    # Include based on flags
    if 'Potential' in category and args.show_potentials:
        return True
    if 'Commensal' in category and args.show_commensals:
        return True
    if category == 'Unknown' and args.show_unidentified:
        return True

    return False


def passes_confidence_threshold(strain, threshold):
    """
    Check if a strain passes the confidence threshold.

    Args:
        strain: Strain dictionary
        threshold: Minimum confidence value (0-100 scale)

    Returns:
        bool: True if strain passes threshold
    """
    # Convert threshold to 0-1 scale if needed
    if threshold > 1.0:
        threshold = threshold / 100.0

    tass_score = strain.get('tass_score', 0)
    return tass_score >= threshold


def load_json_samples(input_files):
    """
    Load and merge multiple JSON files containing sample data.

    Args:
        input_files: List of JSON file paths

    Returns:
        list: All sample data from all JSON files
    """
    all_sample_data = []
    for input_file in input_files:
        with open(input_file, 'r') as f:
            data = json.load(f)
            all_sample_data.extend(data)

    return all_sample_data


def organize_data_by_sample(sample_data):
    """
    Organize the loaded JSON data by sample name.

    - "species_groups" = top-level organisms (first level of dict)
    - "strains" = members within each species group (where tables are made)

    Args:
        sample_data: List of dictionaries from JSON

    Returns:
        dict: {sample_name: [list of species_groups]}
    """
    samples_dict = {}

    for species_group in sample_data:
        # Get sample name
        sample_name = species_group.get('sample_name', 'Unknown Sample')

        if sample_name not in samples_dict:
            samples_dict[sample_name] = []

        samples_dict[sample_name].append(species_group)

    return samples_dict


def sanitize_bookmark_name(name):
    """
    Sanitize a string to be used as a bookmark/anchor name.

    Args:
        name: String to sanitize

    Returns:
        str: Sanitized bookmark name
    """
    # Replace spaces and special characters with underscores
    import re
    sanitized = re.sub(r'[^\w\-]', '_', str(name))
    # Remove consecutive underscores
    sanitized = re.sub(r'_+', '_', sanitized)
    # Remove leading/trailing underscores
    sanitized = sanitized.strip('_')
    return sanitized


def collect_all_bookmarks(samples_dict, low_confidence_strains):
    """
    Collect all bookmark names that will be created in the PDF.
    This is used to validate links before creating them.

    Args:
        samples_dict: Dictionary organized by sample name
        low_confidence_strains: List of low confidence strains

    Returns:
        set: Set of all bookmark names that will exist
    """
    bookmarks = set()

    # Add reference section bookmarks
    bookmarks.add('color_key')
    bookmarks.add('column_explanations')

    # Only add low_confidence bookmark if there are low confidence strains
    if low_confidence_strains:
        bookmarks.add('low_confidence')

    bookmarks.add('additional_info')

    # Add sample and species bookmarks
    for sample_name, species_groups in samples_dict.items():
        sample_bookmark = f"sample_{sanitize_bookmark_name(sample_name)}"
        bookmarks.add(sample_bookmark)

        for species_group in species_groups:
            species_key = species_group.get('toplevelkey', species_group.get('key', 'unknown'))
            species_bookmark = f"species_{sanitize_bookmark_name(sample_name)}_{species_key}"
            bookmarks.add(species_bookmark)

    return bookmarks


def create_safe_link(text, bookmark, valid_bookmarks, color="blue"):
    """
    Create a safe link that only links if the bookmark exists.

    Args:
        text: Display text
        bookmark: Bookmark name to link to
        valid_bookmarks: Set of valid bookmark names
        color: Link color

    Returns:
        str: HTML string with link if bookmark exists, plain text otherwise
    """
    if bookmark in valid_bookmarks:
        return f'<link href="#{bookmark}" color="{color}">{text}</link>'
    else:
        # Return plain text if bookmark doesn't exist
        print(f"Warning: Bookmark '{bookmark}' does not exist, skipping link")
        return text


def build_taxid_to_bookmark_map(samples_dict):
    """
    Build a mapping of taxid (key) to bookmark names for internal linking.

    Args:
        samples_dict: Dictionary organized by sample name

    Returns:
        dict: {taxid: bookmark_name}
    """
    taxid_map = {}

    for sample_name, species_groups in samples_dict.items():
        sample_bookmark_part = sanitize_bookmark_name(sample_name)

        for species_group in species_groups:
            species_key = species_group.get('toplevelkey', species_group.get('key', 'unknown'))

            # Create bookmark for this species group
            species_bookmark = f"species_{sample_bookmark_part}_{species_key}"

            # Map the top-level key
            taxid_map[str(species_key)] = species_bookmark

            # Also map all member strain keys to the same bookmark (they're in the same table)
            for strain in species_group.get('members', []):
                strain_key = str(strain.get('key', ''))
                if strain_key:
                    taxid_map[strain_key] = species_bookmark

    return taxid_map


def get_sample_stats(species_groups):
    """
    Calculate total alignments and primary pathogen count for a sample.

    Args:
        species_groups: List of species groups for a sample

    Returns:
        tuple: (total_alignments, primary_pathogen_count)
    """
    total_alignments = sum(sg.get('numreads', 0) for sg in species_groups)

    primary_pathogen_count = 0
    for sg in species_groups:
        for member in sg.get('members', []):
            if 'Primary' in str(member.get('microbial_category', '')):
                primary_pathogen_count += 1

    return total_alignments, primary_pathogen_count
def has_min_reads(strain, min_reads=1):
    """
    Return True if the strain has at least min_reads aligned reads.
    Treat missing/None as 0.
    """
    try:
        reads = float(strain.get("numreads", 0) or 0)
    except Exception:
        reads = 0
    return reads >= min_reads

def get_species_group_stats(species_group):
    """
    Calculate stats for a species group.

    Args:
        species_group: Species group dictionary

    Returns:
        tuple: (group_tass_score, primary_pathogen_count, max_member_tass_score)
    """
    members = species_group.get('members', [])

    primary_count = sum(1 for m in members if 'Primary' in str(m.get('microbial_category', '')))

    # Get the group's own TASS score
    group_tass = species_group.get('tass_score', 0)

    # Get max TASS score from members for reference
    tass_scores = [m.get('tass_score', 0) for m in members]
    max_member_tass = max(tass_scores) if tass_scores else group_tass

    return group_tass, primary_count, max_member_tass


def create_pdf_template(output_path, samples_dict, ani_data, args):
    """
    Create a PDF report with:
    - Table of contents with clickable sample names and species groups
    - Tables for each strain within each species group
    - Optional High ANI column if any matches exist
    - PDF bookmarks for navigation

    Args:
        output_path: Path to save PDF
        samples_dict: Dictionary organized by sample name
        ani_data: ANI matrix dictionary
        args: Command line arguments
    """
    # Create the PDF document with narrow margins (1% on left/right)
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

    # Build taxid to bookmark mapping for ANI links
    taxid_to_bookmark = build_taxid_to_bookmark_map(samples_dict)

    # Custom styles
    legend_text_style = ParagraphStyle(
        'LegendText',
        parent=styles['Normal'],
        fontSize=8,
        leading=10,
    )
    title_style = ParagraphStyle(
        'CustomTitle',
        parent=styles['Heading1'],
        fontSize=24,
        textColor=colors.HexColor('#2C3E50'),
        spaceAfter=20,
        alignment=TA_CENTER
    )

    # NOTE: anchor_style removed — we now use AnchorFlowable() instead of
    # Paragraph('<a name="..."/>', anchor_style). AnchorFlowable creates proper
    # PDF named destinations that <link href="#..."> can target correctly.

    heading_style = ParagraphStyle(
        'CustomHeading',
        parent=styles['Heading2'],
        fontSize=16,
        textColor=colors.HexColor('#34495E'),
        spaceAfter=0,
        spaceBefore=12
    )
    indent_style = ParagraphStyle(
        'SmallText',
        parent=styles['Normal'],
        fontSize=10,
        leading=10,
        leftIndent=20,
    )

    small_style = ParagraphStyle(
        'SmallText',
        parent=styles['Normal'],
        fontSize=10,
        leading=10
    )

    # Collect all strains across all samples to check for ANY high ANI matches
    all_strains = []
    for sample_name, species_groups in samples_dict.items():
        for species_group in species_groups:
            all_strains.extend(species_group.get('members', []))

    # Check if we should show the High ANI column
    show_ani_column = check_if_any_high_ani_in_dataset(all_strains, ani_data, args.ani_threshold)
    print(f"\nHigh ANI column will be {'SHOWN' if show_ani_column else 'HIDDEN'} (threshold: {args.ani_threshold})")

    # Check if we should show K2 Reads column
    show_k2_column = check_if_k2_reads_present(all_strains)
    print(f"K2 Reads column will be {'SHOWN' if show_k2_column else 'HIDDEN'}")

    # Filter strains based on category flags
    print(f"\nFiltering settings:")
    print(f"  Show Potentials: {args.show_potentials}")
    print(f"  Show Commensals: {args.show_commensals}")
    print(f"  Show Opportunistic: {args.show_opportunistic}")
    print(f"  Show Unidentified: {args.show_unidentified}")
    print(f"  Sorting mode: {'Alphabetical' if args.sort_alphabetical else 'TASS Score (descending)'}")
    print(f"  Max members per group: {args.max_members if args.max_members else 'Unlimited'}")
    print(f"  Max TOC groups per sample: {args.max_toc}")

    # Collect low confidence strains across all samples for later display
    low_confidence_strains = []

    # First pass: collect low confidence strains
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

    # Collect all valid bookmarks that will be created (AFTER determining low_confidence_strains)
    valid_bookmarks = collect_all_bookmarks(samples_dict, low_confidence_strains)

    # Add taxid bookmarks to valid set
    valid_bookmarks.update(taxid_to_bookmark.values())

    # Title
    story.append(Paragraph("Organism Discovery Report", title_style))
    story.append(Spacer(1, 0.05*inch))

    # Report generation info
    metadata_style = styles['Normal']
    current_time = datetime.now().strftime("%Y-%m-%d %H:%M")
    generation_text = (f"This report was generated using <link href=\"https://github.com/jhuapl-bio/taxtriage\" color=\"blue\">TaxTriage</link> <b>{args.version}</b> on <b>{current_time}</b> and is derived from "
                      f"an in <link href=\"https://github.com/jhuapl-bio/taxtriage/blob/main/assets/pathogen_sheet.csv\" color=\"blue\">development spreadsheet of human-host pathogens</link>. Samples with confidence score below "
                      f"{args.min_conf} were filtered out.")
    story.append(Paragraph(generation_text, metadata_style))
    story.append(Spacer(1, 0.02*inch))

    story.append(Paragraph("<b>★</b> = High Consequence Pathogen", small_style))
    story.append(Spacer(1, 0.02*inch))

    # Table of Contents with clickable links
    story.append(Paragraph("Table of Contents", heading_style))
    story.append(Spacer(1, 0.03*inch))

    # Add explanation for TOC format
    toc_explanation = (
        "<i>Format: Sample Name (Total Alignments # - # Primary Pathogens) → "
        "Category Label ■ (# Primary Strains, Max TASS)</i>"
    )
    story.append(Paragraph(toc_explanation, small_style))
    expl = "Click on sample names or species groups to jump to their sections. Only samples/groups with visible strains are shown here"
    story.append(Spacer(1, 0.15*inch))
    story.append(Paragraph(expl, small_style))

    # Sort samples alphabetically
    for sample_name in sorted(samples_dict.keys()):
        # Create a bookmark/anchor name for this sample
        bookmark_name = f"sample_{sanitize_bookmark_name(sample_name)}"

        # Calculate sample stats
        species_groups = samples_dict[sample_name]
        total_alignments, primary_count = get_sample_stats(species_groups)

        # Sort species groups by TASS score (descending) or alphabetically
        if args.sort_alphabetical:
            sorted_groups = sorted(
                species_groups,
                key=lambda sg: sg.get('toplevelname', 'Unknown')
            )
        else:
            # Sort by group TASS score (descending)
            sorted_groups = sorted(
                species_groups,
                key=lambda sg: get_species_group_stats(sg)[0],  # group_tass
                reverse=True
            )

        # Filter species groups that have at least one visible member
        visible_groups = []
        for species_group in sorted_groups:
            # Check if this group has any members that would be shown
            has_visible_members = False
            for strain in species_group.get('members', []):
                if (should_include_strain(strain, args)
                        and has_min_reads(strain, 1)
                        and passes_confidence_threshold(strain, args.min_conf)):
                    has_visible_members = True
                    break

            if has_visible_members:
                visible_groups.append(species_group)

        # Only show sample in TOC if it has visible groups
        if not visible_groups:
            continue

        # Create clickable link to the sample section with stats
        link_text = create_safe_link(
            f'{sample_name} ({total_alignments:,} Alignments - {primary_count} Primary Pathogens)',
            bookmark_name,
            valid_bookmarks
        )
        story.append(Paragraph(f"{link_text}", heading_style))
        story.append(Spacer(1, 0.04*inch))

        # Apply max_toc limit to visible groups in TOC
        toc_display_groups = visible_groups[:args.max_toc]
        has_more_groups = len(visible_groups) > args.max_toc

        species_links = []
        for species_group in toc_display_groups:
            species_name = species_group.get('toplevelname', 'Unknown')
            species_key = species_group.get('toplevelkey', species_group.get('key', 'unknown'))
            species_bookmark = f"species_{sanitize_bookmark_name(sample_name)}_{species_key}"

            # Get species group stats and color
            group_tass, primary_count_sg, max_member_tass = get_species_group_stats(species_group)
            group_category = species_group.get('microbial_category', 'Unknown')
            group_ann_class = species_group.get('annClass', '')

            # Get color for the square
            category_color = get_category_color(group_category, group_ann_class, alpha=1.0)
            color_hex = '#{:02x}{:02x}{:02x}'.format(
                int(category_color.red * 255),
                int(category_color.green * 255),
                int(category_color.blue * 255)
            )

            species_link = (
                f'{create_safe_link(species_name, species_bookmark, valid_bookmarks)} '
                f'<font color="{color_hex}">■</font> '
                f'({primary_count_sg}, {group_tass*100:.1f})'
            )
            species_links.append(species_link)

        # If there are more groups, add "..." with link to first hidden group
        if has_more_groups:
            first_hidden = visible_groups[args.max_toc]
            hidden_species_key = first_hidden.get('toplevelkey', first_hidden.get('key', 'unknown'))
            hidden_bookmark = f"species_{sanitize_bookmark_name(sample_name)}_{hidden_species_key}"
            more_link = create_safe_link(
                f"... ({len(visible_groups) - args.max_toc} more)",
                hidden_bookmark,
                valid_bookmarks
            )
            species_links.append(more_link)

        species_list = ", ".join(species_links)
        story.append(Paragraph(f"→ {species_list}", indent_style))
        story.append(Spacer(1, 0.02*inch))

    # Add links to reference sections at bottom of TOC
    story.append(Spacer(1, 0.05*inch))
    story.append(Paragraph(
        '• ' + create_safe_link('Color Key', 'color_key', valid_bookmarks),
        styles['Normal']
    ))
    story.append(Paragraph(
        '• ' + create_safe_link('Column Explanations', 'column_explanations', valid_bookmarks),
        styles['Normal']
    ))

    # Only add link to low confidence section if it exists
    if low_confidence_strains:
        story.append(Paragraph(
            '• ' + create_safe_link('Low Confidence, High Consequence Detections', 'low_confidence', valid_bookmarks),
            styles['Normal']
        ))

    story.append(Paragraph(
        '• ' + create_safe_link('Additional Information', 'additional_info', valid_bookmarks),
        styles['Normal']
    ))

    story.append(Spacer(1, 0.00*inch))

    # Generate content for each sample - tables appear right after TOC
    for sample_name in sorted(samples_dict.keys()):
        bookmark_name = f"sample_{sanitize_bookmark_name(sample_name)}"

        # FIX: Use AnchorFlowable instead of Paragraph with <a name> tag.
        # AnchorFlowable creates a proper PDF named destination at the correct
        # page position, enabling <link href="#bookmark"> to navigate accurately.
        story.append(AnchorFlowable(bookmark_name))

        # Sample header (smaller, appears right above its table)
        # get the sampletype
        if len(samples_dict[sample_name]) > 0:
            sampletype = samples_dict[sample_name][0].get('sampletype', 'Unspecified Type')
        else:
            sampletype = 'Unspecified Type'
        story.append(Paragraph(f"Sample: {sample_name} ({sampletype})", heading_style))
        # add spacer
        story.append(Spacer(1, 0.1*inch))
        species_groups = samples_dict[sample_name]

        # Calculate sample total reads
        sample_total_reads = sum(sg.get('numreads', 0) for sg in species_groups)
        if sample_total_reads <= 0:
            sample_total_reads = 0

        # Sort species groups by TASS score (descending) or alphabetically
        if args.sort_alphabetical:
            sorted_groups = sorted(
                species_groups,
                key=lambda sg: sg.get('toplevelname', 'Unknown')
            )
        else:
            # Sort by group TASS score (descending)
            sorted_groups = sorted(
                species_groups,
                key=lambda sg: get_species_group_stats(sg)[0],  # group_tass
                reverse=True
            )

        # FIX: Add AnchorFlowable for each species group BEFORE the combined table.
        # Previously these anchors were scattered before the table was built; now we
        # place them sequentially so each gets the correct page/position.
        for species_group in sorted_groups:
            species_key = species_group.get('toplevelkey', species_group.get('key', 'unknown'))
            species_bookmark = f"species_{sanitize_bookmark_name(sample_name)}_{species_key}"
            story.append(AnchorFlowable(species_bookmark))

        # Collect all strains from all species groups for this sample
        all_sample_strains = []
        species_group_map = {}  # Map each strain to its species group

        for species_group in sorted_groups:
            # Get strains for this group
            group_strains = []
            for strain in species_group.get('members', []):
                # Check category filter first
                if should_include_strain(strain, args) and has_min_reads(strain, 1):
                    if passes_confidence_threshold(strain, args.min_conf):
                        group_strains.append(strain)
                        species_group_map[id(strain)] = species_group

            # Sort strains within this group by TASS score (descending)
            group_strains.sort(key=lambda s: s.get('tass_score', 0), reverse=True)
            # Apply max_members filter if specified
            if args.max_members is not None and args.max_members > 0:
                group_strains = group_strains[:args.max_members]

            # Add sorted strains to the overall list
            all_sample_strains.extend(group_strains)

        if all_sample_strains:
            # Create one combined table for all strains in this sample
            combined_table = create_combined_sample_table(
                all_sample_strains, species_group_map, small_style, ani_data,
                args.ani_threshold, show_ani_column, show_k2_column,
                taxid_to_bookmark, valid_bookmarks,
                sample_total_reads=sample_total_reads,
                sample_name=sample_name,
                available_width=available_width
            )
            story.append(combined_table)
        else:
            story.append(Paragraph("<i>No data available for this sample above confidence threshold</i>", styles['Italic']))
        if args.max_members is not None and args.max_members > 0:
            story.append(Paragraph(f"<i>Showing top {args.max_members} strains per group by TASS score</i>", small_style))
        story.append(Spacer(1, 0.1*inch))

    # Add footer documentation section
    story.append(Spacer(1, 0.15*inch))

    # FIX: AnchorFlowable for Color Key section
    story.append(AnchorFlowable('color_key'))
    story.append(Paragraph("<b>Color Key:</b>", heading_style))
    url = 'https://github.com/jhuapl-bio/taxtriage/blob/main/assets/pathogen_sheet.csv'

    # Create color legend table
    legend_data = [
        ['', Paragraph('Category', legend_text_style), Paragraph('Description', legend_text_style)],
        ['', Paragraph('Primary Pathogen (Direct)', legend_text_style),
        Paragraph(f'The organism is directly detected according to <link href="{url}" color="blue">taxonomic id</link> or assembly name and is listed as being of importance in your sample type.', legend_text_style)],
        ['', Paragraph('Primary Pathogen', legend_text_style),
        Paragraph('It is directly detected according to taxonomic id or assembly name and is listed as being of importance in a different sample type.', legend_text_style)],
        ['', Paragraph('Commensal', legend_text_style),
        Paragraph('Normal flora / non-pathogenic', legend_text_style)],
        ['', Paragraph('Opportunistic', legend_text_style),
        Paragraph('May cause disease in certain conditions', legend_text_style)],
        ['', Paragraph('Potential', legend_text_style),
        Paragraph('Potential pathogen requiring further investigation', legend_text_style)],
        ['', Paragraph('Unknown', legend_text_style),
        Paragraph('Classified organism with unknown significance', legend_text_style)],
    ]


    legend_table = Table(legend_data, colWidths=[0.3*inch, 1.6*inch, 3.3*inch])
    legend_table.setStyle(TableStyle([
        # Header
        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#3498DB')),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, 0), 9),
        ('BOTTOMPADDING', (0, 0), (-1, 0), 8),

        # Color indicators
        ('BACKGROUND', (0, 1), (0, 1), colors.lightcoral),
        ('BACKGROUND', (0, 2), (0, 2), colors.HexColor('#fab462')),
        ('BACKGROUND', (0, 3), (0, 3), colors.lightgreen),
        ('BACKGROUND', (0, 4), (0, 4), colors.HexColor('#ffe6a8')),
        ('BACKGROUND', (0, 5), (0, 5), colors.lightblue),
        ('BACKGROUND', (0, 6), (0, 6), colors.white),
        # Grid and alignment
        ('GRID', (0, 0), (-1, -1), 1, colors.black),
        ('FONTSIZE', (1, 1), (-1, -1), 8),
        ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
        ('LEFTPADDING', (1, 1), (-1, -1), 6),
    ]))

    story.append(legend_table)
    story.append(Spacer(1, 0.02*inch))

    # FIX: AnchorFlowable for Column Explanations section
    story.append(AnchorFlowable('column_explanations'))
    story.append(Paragraph("<b>Column Explanations:</b>", heading_style))
    story.append(Spacer(1, 0.03*inch))

    column_explanations = [
        "• <b>Specimen ID (Type):</b> The unique identifier for the sample, including the type of specimen (e.g., blood, tissue).",
        "• <b>Detected Organism:</b> The organism detected in the sample, which could be a bacterium, virus, fungus, or parasite.",
        "• <b>Microbial Category:</b> The classification of the organism, indicating whether it is primary, opportunistic, commensal, or potential.",
        "• <b># Reads Aligned:</b> The number of reads from the sequencing data that align to the organism's genome, indicating its presence. (%) refers to all alignments (more than 1 alignment per read can take place) for that species across the entire sample. The format is (total % of aligned reads in sample).",
        "• <b>RPM:</b> Reads Per Million (RPM). This normalized metric allows for comparison of abundance across samples and organisms of different sizes.",
        "• <b>TASS Score:</b> A metric between 0 and 100 that reflects the confidence of the organism's detection, with 100 being the highest value.",
        "• <b>Taxonomic ID #:</b> The taxid for the organism according to NCBI Taxonomy, which provides a unique identifier for each species. The parenthesis (if present) is the group it belongs to, usually the genus.",
        "• <b>Pathogenic Subsp/Strains:</b> Indicates specific pathogenic subspecies, serotypes, or strains, if detected in the sample. (%) indicates the percent of all aligned reads belonging to that strain.",
        "• <b>K2 Reads:</b> The number of reads classified by Kraken2, a tool for taxonomic classification of sequencing data.",
    ]

    for explanation in column_explanations:
        story.append(Paragraph(explanation, metadata_style))
        story.append(Spacer(1, 0.02*inch))

    story.append(Spacer(1, 0.1*inch))

    # FIX: AnchorFlowable for Low Confidence section (only if it exists)
    if low_confidence_strains:
        story.append(AnchorFlowable('low_confidence'))
        story.append(Paragraph("<b>Low Confidence, High Consequence Detections:</b>", heading_style))
        explanation_text = (
            f"The following strains were detected but fell below the confidence threshold "
            f"of {args.min_conf} and are listed here for reference."
        )
        story.append(Paragraph(explanation_text, metadata_style))
        story.append(Spacer(1, 0.05*inch))

        low_conf_table = create_low_confidence_table(
            low_confidence_strains, small_style, show_k2_column, available_width
        )
        story.append(low_conf_table)
        story.append(Spacer(1, 0.1*inch))

    # FIX: AnchorFlowable for Additional Information section
    story.append(AnchorFlowable('additional_info'))
    story.append(Paragraph("<b>Additional Information:</b>", heading_style))
    story.append(Spacer(1, 0.01*inch))

    url = "https://github.com/jhuapl-bio/taxtriage/blob/main/docs/usage.md#confidence-scoring"
    text_with_url = f'Please visit our <a href="{url}"><b><font color="blue">DOCUMENTATION PAGE</font></b></a> for more information on how confidence is calculated.'
    story.append(Paragraph(text_with_url, metadata_style))
    story.append(Spacer(1, 0.01*inch))

    story.append(Paragraph("The following information highlights the description for the color combinations for each organism class in the annotated table(s).", metadata_style))
    story.append(Spacer(1, 0.01*inch))

    story.append(Paragraph("Please see the relevant Discovery Analysis txt file for low confidence, high consequence annotations that were not present in the pdf.", metadata_style))
    story.append(Spacer(1, 0.01*inch))

    story.append(Paragraph("Read amounts are represented as the <b>total number of aligned reads</b> of sufficient mapping quality <b>(% aligned for all reads in sample)</b>.", metadata_style))
    story.append(Spacer(1, 0.01*inch))
    story.append(Paragraph("If there are questions or issues with your report, please open an issue on GitHub as a discussion <link href=\"https://github.com/jhuapl-bio/taxtriage/discussions\" color=\"blue\">here</link>. Issues should be tracked/submitted at <link href=\"https://github.com/jhuapl-bio/taxtriage/issues\" color=\"blue\">this link</link>.", metadata_style))

    # FIX: Build the PDF directly — no PyPDF2 post-processing needed.
    # The AnchorFlowable calls above embed proper named destinations in the PDF
    # during doc.build(), so all <link href="#..."> internal links work correctly.
    # The old PyPDF2 block was broken because it hardcoded all bookmarks to page 0.
    doc.build(story)

    print(f"\nPDF created successfully: {output_path}")


def get_category_color(microbial_category, ann_class, alpha=1.0):
    """
    Get the color for a row based on microbial category and annotation class.

    Args:
        microbial_category: String like "Primary Pathogen", "Commensal", etc.
        ann_class: String like "Direct" or other classification
        alpha: Transparency level (0.0 to 1.0)

    Returns:
        Color: ReportLab color object with specified alpha
    """
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

    # Apply alpha if not fully opaque
    if alpha < 1.0:
        # Create a new color with alpha
        return colors.Color(base_color.red, base_color.green, base_color.blue, alpha=alpha)

    return base_color


def create_combined_sample_table(all_strains, species_group_map, small_style, ani_data, ani_threshold, show_ani_column, show_k2_column, taxid_to_bookmark, valid_bookmarks, sample_total_reads=0, sample_name=None, available_width=None):
    """
    Create a single table combining all strains from all species groups.
    Each species group gets its own header row, followed by its strains.

    Args:
        all_strains: List of all strain dictionaries from all species groups
        species_group_map: Map of strain id() to its parent species group
        small_style: ParagraphStyle for small text
        ani_data: ANI matrix dictionary
        ani_threshold: Minimum ANI value to report
        show_ani_column: Whether to include the High ANI column
        show_k2_column: Whether to include the K2 Reads column
        taxid_to_bookmark: Mapping of taxid to bookmark names for internal links
        valid_bookmarks: Set of valid bookmark names for link validation
        sample_total_reads: Total reads for the sample
        sample_name: Name of the sample
        available_width: Available width for the table

    Returns:
        Table: Formatted ReportLab table object
    """
    # Create styles
    strain_name_style = ParagraphStyle(
        'StrainName',
        parent=small_style,
        fontSize=10,
        leading=10
    )

    data_style = ParagraphStyle(
        'DataStyle',
        parent=small_style,
        fontSize=8,
        leading=9
    )

    ani_style = ParagraphStyle(
        'ANIStyle',
        parent=small_style,
        fontSize=6,
        leading=8
    )

    group_header_style = ParagraphStyle(
        'GroupHeader',
        parent=small_style,
        fontSize=10,
        leading=11,
        fontName='Helvetica-Bold'
    )

    # Table headers - build dynamically based on what columns to show
    # Order: Star, Strain Name, TASS, K2 Reads (if shown), Reads, RPKM/RPM, Coverage, High ANI (if shown)
    headers = ['', 'Strain Name', 'TASS', 'K2 Reads', 'Reads', 'RPM', 'Coverage', 'High ANI'] if show_k2_column and show_ani_column else (
        ['', 'Strain Name', 'TASS', 'K2 Reads', 'Reads', 'RPM', 'Coverage'] if show_k2_column else (
            ['', 'Strain Name', 'TASS', 'Reads', 'RPM', 'Coverage', 'High ANI'] if show_ani_column else
            ['', 'Strain Name', 'TASS', 'Reads', 'RPM', 'Coverage']
        )
    )

    table_data = [headers]

    # Track which species groups we've seen
    current_species_group = None
    group_row_indices = []  # Track row indices that are group headers

    row_idx = 1  # Start after header

    for strain in all_strains:
        species_group = species_group_map[id(strain)]
        species_key = species_group.get('toplevelkey', species_group.get('key', 'unknown'))

        # Check if we're starting a new species group
        if current_species_group != species_key:
            # Start new group - add group header row
            current_species_group = species_key

            group_name = species_group.get('toplevelname', 'Unknown')
            group_reads = species_group.get('numreads', 0)
            group_k2_reads = species_group.get('k2_reads', 0)
            group_category = species_group.get('microbial_category', 'Unknown')
            group_ann_class = species_group.get('annClass', '')

            # Get color for the group row background
            category_color = get_category_color(group_category, group_ann_class, alpha=1.0)

            # Create group name with link
            group_name_para = Paragraph(
                f'<b><link href="https://www.ncbi.nlm.nih.gov/taxonomy/?term={species_key}" color="blue">{group_name}</link></b>',
                group_header_style
            )

            # Build group header row
            # Order: Star, Strain Name, TASS, K2 Reads (if shown), Reads, RPKM/RPM, Coverage, High ANI (if shown)
            group_row = [
                '',  # Star column - empty
                group_name_para,  # Group name
                '',  # TASS - blank for groups
            ]

            if show_k2_column:
                group_row.append(Paragraph(f'<b>{group_k2_reads:,.0f}</b>', group_header_style))  # K2 Reads

            group_row.append(Paragraph(f'<b>{group_reads:,.0f}</b>', group_header_style))  # Reads aligned
            group_row.append('')  # RPKM/RPM - blank for groups
            group_row.append('')  # Coverage - blank

            if show_ani_column:
                group_row.append('')  # High ANI - blank

            table_data.append(group_row)
            group_row_indices.append(row_idx)
            row_idx += 1

        # Build strain data
        microbial_category = strain.get('microbial_category', 'Unknown')
        ann_class = strain.get('annClass', '')
        is_high_consequence = strain.get('high_cons', False)

        strain_name_text = strain.get('name', 'Unknown')
        strain_key = strain.get('key', '')

        strain_name_text = (
            f'{strain_name_text} '
            f'(<link href="https://www.ncbi.nlm.nih.gov/taxonomy/?term={strain_key}" '
            f'color="blue">{strain_key}</link>)'
        )

        strain_name = Paragraph(strain_name_text, strain_name_style)
        indicator_text = '★' if is_high_consequence else ''

        # Get High ANI matches if showing column
        high_ani_text = ''
        if show_ani_column:
            matches = get_high_ani_matches(strain_key, ani_data, ani_threshold, set())
            if matches:
                ani_links = []
                for taxid, ani_value in matches[:3]:
                    ani_percent = ani_value * 100
                    if taxid in taxid_to_bookmark:
                        bookmark = taxid_to_bookmark[taxid]
                        # Use safe link that checks if bookmark exists
                        if bookmark in valid_bookmarks:
                            link = f'<link href="#{bookmark}" color="blue">{taxid} ({ani_percent:.1f}%)</link>'
                        else:
                            # Bookmark doesn't exist, just show text
                            link = f"{taxid} ({ani_percent:.1f}%)"
                        ani_links.append(link)
                    else:
                        ani_links.append(f"{taxid} ({ani_percent:.1f}%)")
                high_ani_text = Paragraph(", ".join(ani_links), ani_style)
            else:
                high_ani_text = Paragraph("-", ani_style)

        strain_reads = float(strain.get('numreads', 0) or 0)
        pct = (strain_reads / sample_total_reads * 100.0) if sample_total_reads else 0.0
        reads_text = f"{strain_reads:,.0f} ({pct:.1f}%)"
        reads_paragraph = Paragraph(reads_text, data_style)

        # Get RPKM and RPM values
        rpkm = strain.get('rpkm', 0) or 0
        rpm = strain.get('rpm', 0) or 0
        rpkm_rpm_text = f"{rpm:,.0f}"
        rpkm_rpm_paragraph = Paragraph(rpkm_rpm_text, data_style)

        # Build row - Order: Star, Strain Name, TASS, K2 Reads (if shown), Reads, RPKM/RPM, Coverage, High ANI (if shown)
        row = [
            indicator_text,
            strain_name,
            Paragraph(f"{strain.get('tass_score', 0)*100:.1f}", data_style),  # TASS
        ]

        if show_k2_column:
            row.append(Paragraph(f"{strain.get('k2_reads', 0):,.0f}", data_style))  # K2 Reads

        row.append(reads_paragraph)  # Reads
        row.append(rpkm_rpm_paragraph)  # RPKM (RPM)
        row.append(Paragraph(f"{min(100,strain.get('coverage', 0)*100):.1f}%", data_style))  # Coverage

        if show_ani_column:
            row.append(high_ani_text)  # High ANI

        table_data.append(row)
        row_idx += 1

    # Create table with ADAPTIVE column widths based on available_width
    if available_width is None:
        # Fallback to letter size calculation
        available_width = 8.5*inch - 0.02*8.5*inch  # letter width minus 1% margins on each side

    # Calculate proportional widths based on available space
    # Order: Star, Strain Name, TASS, K2 Reads (if shown), Reads, RPKM/RPM, Coverage, High ANI (if shown)
    if show_k2_column and show_ani_column:
        col_widths = [
            available_width * 0.03,  # Star
            available_width * 0.30,  # Strain Name
            available_width * 0.08,  # TASS
            available_width * 0.10,  # K2 Reads
            available_width * 0.11,  # Reads
            available_width * 0.10,  # RPKM (RPM)
            available_width * 0.09,  # Coverage
            available_width * 0.19   # High ANI
        ]
    elif show_k2_column:
        col_widths = [
            available_width * 0.03,  # Star
            available_width * 0.36,  # Strain Name
            available_width * 0.09,  # TASS
            available_width * 0.13,  # K2 Reads
            available_width * 0.13,  # Reads
            available_width * 0.12,  # RPKM (RPM)
            available_width * 0.14   # Coverage
        ]
    elif show_ani_column:
        col_widths = [
            available_width * 0.03,  # Star
            available_width * 0.34,  # Strain Name
            available_width * 0.09,  # TASS
            available_width * 0.13,  # Reads
            available_width * 0.11,  # RPKM (RPM)
            available_width * 0.10,  # Coverage
            available_width * 0.20   # High ANI
        ]
    else:
        col_widths = [
            available_width * 0.03,  # Star
            available_width * 0.43,  # Strain Name
            available_width * 0.10,  # TASS
            available_width * 0.15,  # Reads
            available_width * 0.13,  # RPKM (RPM)
            available_width * 0.16   # Coverage
        ]

    table = Table(table_data, repeatRows=1, colWidths=col_widths)

    # Build table style
    table_styles = [
        # Header styling
        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#3498DB')),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
        ('ALIGN', (0, 0), (-1, 0), 'CENTER'),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, 0), 9),
        ('BOTTOMPADDING', (0, 0), (-1, 0), 12),

        # Indicator column styling (star column)
        ('ALIGN', (0, 1), (0, -1), 'CENTER'),
        ('VALIGN', (0, 1), (0, -1), 'MIDDLE'),
        ('FONTSIZE', (0, 1), (0, -1), 14),

        # Grid
        ('GRID', (0, 0), (-1, -1), 1, colors.black),
        ('VALIGN', (1, 1), (1, -1), 'MIDDLE'),
        ('VALIGN', (2, 1), (-1, -1), 'MIDDLE'),

        # General body styling
        ('ALIGN', (2, 1), (-1, -1), 'CENTER'),  # All data columns centered
        ('ALIGN', (1, 1), (1, -1), 'LEFT'),  # Strain name left-aligned

        # Padding for strain name column
        ('LEFTPADDING', (1, 1), (1, -1), 6),
        ('RIGHTPADDING', (1, 1), (1, -1), 6),
        ('TOPPADDING', (1, 1), (1, -1), 6),
        ('BOTTOMPADDING', (1, 1), (1, -1), 6),
    ]

    # Add group header row styling
    for group_row_idx in group_row_indices:
        table_styles.append(('BACKGROUND', (0, group_row_idx), (-1, group_row_idx), colors.HexColor('#E8E8E8')))
        table_styles.append(('ALIGN', (1, group_row_idx), (1, group_row_idx), 'LEFT'))
        table_styles.append(('ALIGN', (2, group_row_idx), (-1, group_row_idx), 'CENTER'))
        table_styles.append(('TOPPADDING', (0, group_row_idx), (-1, group_row_idx), 8))
        table_styles.append(('BOTTOMPADDING', (0, group_row_idx), (-1, group_row_idx), 8))

    # Add row-specific colors for strain rows (skip group header rows)
    row_idx = 1
    for strain in all_strains:
        species_group = species_group_map[id(strain)]

        # Skip to next row if this is the start of a new group (group header row)
        if row_idx in group_row_indices:
            row_idx += 1

        microbial_category = strain.get('microbial_category', 'Unknown')
        ann_class = strain.get('annClass', '')

        # Debug output
        print(f"  Row {row_idx}: {strain.get('name', 'Unknown')[:40]} - Category: {microbial_category}, Class: {ann_class}")

        indicator_color = get_category_color(microbial_category, ann_class, alpha=1.0)
        table_styles.append(('BACKGROUND', (0, row_idx), (0, row_idx), indicator_color))

        row_color = get_category_color(microbial_category, ann_class, alpha=0.15)
        table_styles.append(('BACKGROUND', (1, row_idx), (-1, row_idx), row_color))

        row_idx += 1

    table.setStyle(TableStyle(table_styles))
    return table


def create_low_confidence_table(low_confidence_strains, small_style, show_k2_column, available_width=None):
    """
    Create a simple table for low confidence detections.

    Args:
        low_confidence_strains: List of dicts with strain, sample_name, species_group
        small_style: ParagraphStyle for small text
        show_k2_column: Whether to include K2 Reads column
        available_width: Available width for the table

    Returns:
        Table: Formatted ReportLab table object
    """
    # Create styles
    strain_name_style = ParagraphStyle(
        'StrainName',
        parent=small_style,
        fontSize=8,
        leading=10
    )

    data_style = ParagraphStyle(
        'DataStyle',
        parent=small_style,
        fontSize=7,
        leading=9
    )

    # Table headers - build dynamically based on what columns to show
    # Order: Sample, Strain Name, TASS, K2 Reads (if shown), Reads
    headers = ['Sample', 'Strain Name', 'TASS', 'K2 Reads', 'Reads'] if show_k2_column else ['Sample', 'Strain Name', 'TASS', 'Reads']

    table_data = [headers]

    # Add strain rows
    for item in low_confidence_strains:
        strain = item['strain']
        sample_name = item['sample_name']

        microbial_category = strain.get('microbial_category', 'Unknown')
        ann_class = strain.get('annClass', '')
        is_high_consequence = strain.get('high_cons', False)

        strain_name_text = strain.get('name', 'Unknown')
        strain_key = strain.get('key', '')
        if strain_key:
            strain_name_text = (
                f'{strain_name_text} '
                f'(<link href="https://www.ncbi.nlm.nih.gov/taxonomy/?term={strain_key}" '
                f'color="blue">{strain_key}</link>)'
            )
        # Add star if high consequence
        if is_high_consequence:
            strain_name_text = f"★ {strain_name_text}"

        strain_name = Paragraph(strain_name_text, strain_name_style)

        # Build row - Order: Sample, Strain Name, TASS, K2 Reads (if shown), Reads
        row = [
            Paragraph(sample_name, strain_name_style),
            strain_name,
            Paragraph(f"{strain.get('tass_score', 0)*100:.1f}", data_style),  # TASS
        ]

        if show_k2_column:
            row.append(Paragraph(f"{strain.get('k2_reads', 0):,.0f}", data_style))  # K2 Reads

        row.append(Paragraph(f"{strain.get('numreads', 0):,.0f}", data_style))  # Reads

        table_data.append(row)

    # Create table with ADAPTIVE column widths based on available_width
    if available_width is None:
        available_width = 8.5*inch - 0.02*8.5*inch

    if show_k2_column:
        col_widths = [
            available_width * 0.18,  # Sample
            available_width * 0.50,  # Strain Name
            available_width * 0.10,  # TASS
            available_width * 0.10,  # K2 Reads
            available_width * 0.12   # Reads
        ]
    else:
        col_widths = [
            available_width * 0.20,  # Sample
            available_width * 0.56,  # Strain Name
            available_width * 0.11,  # TASS
            available_width * 0.13   # Reads
        ]

    table = Table(table_data, repeatRows=1, colWidths=col_widths)

    # Build table style
    table_styles = [
        # Header styling
        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#3498DB')),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
        ('ALIGN', (0, 0), (-1, 0), 'CENTER'),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, 0), 9),
        ('BOTTOMPADDING', (0, 0), (-1, 0), 12),

        # Grid
        ('GRID', (0, 0), (-1, -1), 1, colors.black),
        ('VALIGN', (0, 1), (-1, -1), 'MIDDLE'),

        # General body styling
        ('ALIGN', (2, 1), (-1, -1), 'CENTER'),
        ('ALIGN', (0, 1), (1, -1), 'LEFT'),

        # Padding
        ('LEFTPADDING', (0, 1), (-1, -1), 6),
        ('RIGHTPADDING', (0, 1), (-1, -1), 6),
        ('TOPPADDING', (0, 1), (-1, -1), 6),
        ('BOTTOMPADDING', (0, 1), (-1, -1), 6),
    ]

    # Add row-specific colors
    for idx, item in enumerate(low_confidence_strains):
        row_idx = idx + 1
        strain = item['strain']
        microbial_category = strain.get('microbial_category', 'Unknown')
        ann_class = strain.get('annClass', '')

        row_color = get_category_color(microbial_category, ann_class, alpha=0.25)
        table_styles.append(('BACKGROUND', (0, row_idx), (-1, row_idx), row_color))

    table.setStyle(TableStyle(table_styles))
    return table


def create_tabular_output(output_path, samples_dict, args):
    """
    Create a tabular output file (CSV/TSV/TXT/XLSX) with strain-level data.
    This includes ALL strains, not just those that pass filters (unlike the PDF).

    Args:
        output_path: Path to save the output file
        samples_dict: Dictionary organized by sample name
        args: Command line arguments
    """
    # Determine file format from extension
    file_ext = os.path.splitext(output_path)[1].lower()

    # Define column headers matching the report.txt format
    headers = [
        'Index',
        'index',
        'Detected Organism',
        'Specimen ID',
        'Sample Type',
        '% Reads',
        '# Reads Aligned',
        '% Aligned Reads',
        'Coverage',
        'HHS Percentile',
        'IsAnnotated',
        'AnnClass',
        'Microbial Category',
        'High Consequence',
        'Taxonomic ID #',
        'Status',
        'Gini Coefficient',
        'Mean BaseQ',
        'Mean MapQ',
        'Mean Depth',
        'isSpecies',
        'Pathogenic Subsp/Strains',
        'K2 Reads',
        'RPKM',
        'RPM',
        'Parent K2 Reads',
        'MapQ Score',
        'Disparity Score',
        'Minhash Score',
        'Diamond Identity',
        'K2 Disparity Score',
        'Siblings score',
        'Breadth Weight Score',
        'TASS Score',
        'MicrobeRT Probability',
        'MicrobeRT Model',
        'Reads Aligned',
        'Group'
    ]

    # Collect all strain data (ALL strains, not filtered)
    all_rows = []
    global_index = 0

    for sample_name in sorted(samples_dict.keys()):
        species_groups = samples_dict[sample_name]

        # Sort species groups
        if args.sort_alphabetical:
            sorted_groups = sorted(species_groups, key=lambda sg: sg.get('toplevelname', 'Unknown'))
        else:
            sorted_groups = sorted(species_groups, key=lambda sg: get_species_group_stats(sg)[0], reverse=True)

        # Calculate sample total reads for percentages
        sample_total_reads = sum(sg.get('numreads', 0) for sg in species_groups)
        if sample_total_reads <= 0:
            sample_total_reads = 1  # Avoid division by zero

        for species_group in sorted_groups:
            group_key = species_group.get('toplevelkey', species_group.get('key', ''))

            # Get ALL strains from this group (no filtering for output file)
            strains = species_group.get('members', [])

            # Sort by TASS score
            strains.sort(key=lambda s: s.get('tass_score', 0), reverse=True)

            # Get top-level sample_type from species_group
            sample_type = species_group.get('sampletype', 'unknown')

            # Create row for each strain
            for local_idx, strain in enumerate(strains):
                if not has_min_reads(strain, 1):
                    continue
                strain_reads = float(strain.get('numreads', 0) or 0)
                pct_reads = (strain_reads / sample_total_reads * 100.0) if sample_total_reads else 0.0

                # Get organism name
                organism_name = strain.get('name', 'Unknown')

                # Build the row with corrected field mappings
                row = [
                    global_index,  # Index
                    local_idx,  # index (within group)
                    organism_name,  # Detected Organism
                    sample_name,  # Specimen ID
                    sample_type,  # Sample Type
                    f"{pct_reads:.4f}",  # % Reads
                    int(strain_reads),  # # Reads Aligned
                    f"{pct_reads:.4f}",  # % Aligned Reads
                    f"{(min(1, strain.get('coverage', 0)) or 0)*100:.0f}%",  # Coverage
                    '100.0',  # HHS Percentile (placeholder)
                    'Yes' if strain.get('isAnnotated', True) else 'No',  # IsAnnotated
                    strain.get('annClass', ''),  # AnnClass
                    strain.get('microbial_category', 'Unknown'),  # Microbial Category
                    'True' if strain.get('high_cons', False) else 'False',  # High Consequence
                    strain.get('key', ''),  # Taxonomic ID #
                    strain.get('status', ''),  # Status
                    f"{(strain.get('gini_coefficient', 0) or 0):.2f}",  # Gini Coefficient
                    f"{(strain.get('meanbaseq', 0) or 0):.2f}",  # Mean BaseQ
                    f"{(strain.get('meanmapq', 0) or 0):.2f}",  # Mean MapQ
                    f"{(strain.get('meandepth', 0) or 0):.1f}",  # Mean Depth
                    'True' if strain.get('isSpecies', False) else 'False',  # isSpecies
                    '',  # Pathogenic Subsp/Strains
                    int(strain.get('k2_reads', 0) or 0),  # K2 Reads
                    strain.get("rpkm", 0) or 0,  # RPKM
                    strain.get("rpm", 0) or 0,  # RPM
                    int(strain.get('parent_k2_reads', 0) or 0),  # Parent K2 Reads
                    f"{(strain.get('mapq_score', 0) or 0):.2f}",  # MapQ Score
                    f"{(strain.get('disparity', 0) or 0):.2f}",  # Disparity Score
                    f"{(strain.get('minhash_reduction', 0) or 0):.2f}",  # Minhash Score
                    f"{(strain.get('diamond_identity', 0) or 0):.1f}",  # Diamond Identity
                    f"{(strain.get('k2_disparity_score', 0) or 0):.1f}",  # K2 Disparity Score
                    f"{(strain.get('siblings_score', 0) or 0):.1f}",  # Siblings score
                    f"{(strain.get('breadth_log_score', 0) or 0):.2f}",  # Breadth Weight Score
                    int((strain.get('tass_score', 0) or 0) * 100),  # TASS Score (as integer 0-100)
                    f"{(strain.get('mmbert', 0) or 0):.4f}",  # MicrobeRT Probability
                    strain.get('mmbert_model', '') or '',  # MicrobeRT Model
                    int(strain_reads),  # Reads Aligned
                    group_key  # Group
                ]

                all_rows.append(row)
                global_index += 1

    # Create DataFrame
    df = pd.DataFrame(all_rows, columns=headers)
    # Save based on file extension
    if file_ext in ['.csv']:
        df.to_csv(output_path, index=False)
        print(f"CSV output created: {output_path}")
    elif file_ext in ['.tsv', '.txt']:
        df.to_csv(output_path, sep='\t', index=False)
        print(f"TSV output created: {output_path}")
    elif file_ext in ['.xlsx', '.xls']:
        df.to_excel(output_path, index=False, engine='openpyxl')
        print(f"Excel output created: {output_path}")
    else:
        # Default to TSV
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
        "-i",
        "--input",
        metavar="INPUT",
        required=True,
        nargs="+",
        default=[],
        help="Base pathogen discovery table file(s), JSON format only. Can specify more than one",
    )
    parser.add_argument(
        "-d",
        "--distributions",
        metavar="DISTRIBUTIONS",
        required=False,
        help="TSV file that contains all the distribution information for body sites and organisms",
    )
    parser.add_argument("-a", "--abundance_col", metavar="ABU", required=False, default='% Aligned Reads',
                        help="Name of abundance column, default is abundance")
    parser.add_argument("-r", "--min_reads", metavar="READS", required=False, default=1, type=int,
                        help="Minimum number of reads required to consider an organism for reporting. Default is 'Min Reads' column in input file.")
    parser.add_argument("-c", "--min_conf", metavar="MINCONF", required=False, default=0.3, type=float,
                        help="Value that must be met for a table to report an organism due to confidence column.")
    parser.add_argument("-x", "--id_col", metavar="IDCOL", required=False, default="Detected Organism",
                        help="Name of id column, default is id")
    parser.add_argument("-p", "--percentile", metavar="PERCENTILE", required=False, type=float, default=0.75,
                        help="Only show organisms that are in the top percentile of healthy subjects expected abu")
    parser.add_argument("-v", "--version", metavar="VERSION", required=False, default='Local Build',
                        help="What version of TaxTriage is in use")
    parser.add_argument('--sorttass', action="store_true", required=False,
                        help="Sort by TASS score if available")
    parser.add_argument('--sort_alphabetical', action="store_true", required=False,
                        help="Sort groups alphabetically instead of by TASS score")
    parser.add_argument("--max_members", metavar="MAX_MEMBERS", required=False, type=int, default=None,
                        help="Maximum number of top strains (by TASS) to show per group. Default: show all")
    parser.add_argument("--max_toc", metavar="MAX_TOC", required=False, type=int, default=4,
                        help="Maximum number of species groups to show in TOC per sample. Default: 4")
    parser.add_argument("--show_commensals", action="store_true", required=False,
                        help="Show the commensals table")
    parser.add_argument("--show_unidentified", action="store_true", required=False,
                        help="Show the all organisms now listed as commensal or pathogen")
    parser.add_argument("--show_potentials", action="store_true", required=False,
                        help="Show the potentials table")
    parser.add_argument("--show_opportunistic", action="store_true", required=False,
                        help="Show the opportunistic table")
    parser.add_argument("-m", "--missing_samples", metavar="MISSING", required=False, default=None,
                        help="Missing samples if any", nargs="+")
    parser.add_argument("-s", "--sitecol", metavar="SCOL", required=False, default='Sample Type',
                        help="Name of site column, default is body_site")
    parser.add_argument("--ani_threshold", metavar="ANI_THRESHOLD", required=False, type=float, default=0.90,
                        help="Threshold for ANI to consider 'High ANI' in the report")
    parser.add_argument("-t", "--type", metavar="TYPE", required=False, default='Detected Organism',
                        help="What type of data is being processed. Options: 'Taxonomic ID #' or 'Detected Organism'.",
                        choices=['Taxonomic ID #', 'Detected Organism'])
    parser.add_argument("--taxdump", metavar="TAXDUMP", required=False, default=None,
                        help="Merge the entries on a specific rank args.rank, importing files from nodes.dmp, names.dmp and potentially merged.dmp")
    parser.add_argument("--rank", metavar="TAXDUMP", required=False, default="genus",
                        help='IF merging with taxdump, what rank to merge on')
    parser.add_argument(
        "-o",
        "--output",
        metavar="OUTPUT",
        required=True,
        type=str,
        help="Path of output file (pdf)",
    )
    parser.add_argument(
        "--ani_matrix",
        metavar="ANI_MATRIX",
        required=False,
        type=str,
        nargs="+",
        help="Path to organism ANI matrix CSV file from Sourmash signatures. Optional",
    )
    parser.add_argument(
        "-u",
        "--output_txt",
        metavar="OUTPUT_TXT",
        required=False,
        type=str,
        help="Path of tabular output file. Format determined by extension: .csv, .tsv, .txt (TSV), or .xlsx",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    # Initialize empty dictionaries for taxonomy data
    taxdump_dict = {}
    names_map = {}
    merged_tax_data = {}

    # Load taxdump files if provided
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

    # Load ANI matrix if provided
    ani_data = {}
    if args.ani_matrix:
        ani_data = load_ani_matrix(args.ani_matrix)
        print(f"Loaded ANI data for {len(ani_data)} taxa")

    # Load JSON sample data
    sample_data = load_json_samples(args.input)
    print(f"Loaded {len(sample_data)} species groups from JSON file(s)")

    # Organize data by sample
    samples_dict = organize_data_by_sample(sample_data)
    print(f"Found {len(samples_dict)} unique sample(s)")

    # Show configuration
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

    # Create the PDF report
    create_pdf_template(args.output, samples_dict, ani_data, args)

    # Create tabular output if requested
    if args.output_txt:
        create_tabular_output(args.output_txt, samples_dict, args)

    return taxdump_dict, names_map, merged_tax_data, sample_data


if __name__ == "__main__":
    taxdump_dict, names_map, merged_tax_data, sample_data = main()
