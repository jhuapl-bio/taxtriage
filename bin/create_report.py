#!/usr/bin/env python3

import json
import argparse
import os
import pandas as pd
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer, PageBreak
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
    if 'Primary' in category or 'Opportunistic' in category:
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


def collect_all_bookmarks(samples_dict):
    """
    Collect all bookmark names that will be created in the PDF.
    This is used to validate links before creating them.

    Args:
        samples_dict: Dictionary organized by sample name

    Returns:
        set: Set of all bookmark names that will exist
    """
    bookmarks = set()

    # Add reference section bookmarks
    bookmarks.add('color_key')
    bookmarks.add('column_explanations')
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

    Args:
        output_path: Path to save PDF
        samples_dict: Dictionary organized by sample name
        ani_data: ANI matrix dictionary
        args: Command line arguments
    """
    # Create the PDF document with narrow margins (2.5% on left/right)
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

    # Collect all valid bookmarks that will be created
    valid_bookmarks = collect_all_bookmarks(samples_dict)

    # Add taxid bookmarks to valid set
    valid_bookmarks.update(taxid_to_bookmark.values())

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
    anchor_style = ParagraphStyle(
        "Anchor",
        parent=styles["Normal"],
        fontSize=1,
        leading=1,
        spaceBefore=0,
        spaceAfter=0,
    )

    heading_style = ParagraphStyle(
        'CustomHeading',
        parent=styles['Heading2'],
        fontSize=16,
        textColor=colors.HexColor('#34495E'),
        spaceAfter=8
    )

    small_style = ParagraphStyle(
        'SmallText',
        parent=styles['Normal'],
        fontSize=8,
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
    print(f"  Show Unidentified: {args.show_unidentified}")
    print(f"  Sorting mode: {'Alphabetical' if args.sort_alphabetical else 'TASS Score (descending)'}")
    print(f"  Max members per group: {args.max_members if args.max_members else 'Unlimited'}")

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
        "Category Label ■ (# Primary Strains, Max TASS%)</i>"
    )
    story.append(Paragraph(toc_explanation, small_style))
    story.append(Spacer(1, 0.05*inch))
    sample_total_reads = sum(sg.get('numreads', 0) for sg in species_groups)

    if sample_total_reads <= 0:
        sample_total_reads = 0


    # Sort samples alphabetically
    for sample_name in sorted(samples_dict.keys()):
        # Create a bookmark/anchor name for this sample
        bookmark_name = f"sample_{sanitize_bookmark_name(sample_name)}"

        # Calculate sample stats
        species_groups = samples_dict[sample_name]
        total_alignments, primary_count = get_sample_stats(species_groups)

        # Create clickable link to the sample section with stats
        link_text = create_safe_link(
            f'{sample_name} ({total_alignments:,}# - {primary_count}*)',
            bookmark_name,
            valid_bookmarks
        )
        story.append(Paragraph(f"• {link_text}", styles['Normal']))

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

        species_links = []
        for species_group in sorted_groups:
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
                f'({primary_count_sg}, {group_tass*100:.1f}%)'
            )
            species_links.append(species_link)

        species_list = ", ".join(species_links)
        story.append(Paragraph(f"  → {species_list}", small_style))
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
    story.append(Paragraph(
        '• ' + create_safe_link('Low Confidence, High Consequence Detections', 'low_confidence', valid_bookmarks),
        styles['Normal']
    ))
    story.append(Paragraph(
        '• ' + create_safe_link('Additional Information', 'additional_info', valid_bookmarks),
        styles['Normal']
    ))

    story.append(Spacer(1, 0.00*inch))

    # Collect low confidence strains across all samples for later display
    low_confidence_strains = []

    # Generate content for each sample - tables appear right after TOC
    for sample_name in sorted(samples_dict.keys()):
        bookmark_name = f"sample_{sanitize_bookmark_name(sample_name)}"

        # Add bookmark/anchor for this sample
        story.append(Paragraph(f'<a name="{bookmark_name}"/>', anchor_style))

        # Sample header (smaller, appears right above its table)
        story.append(Paragraph(f"Sample: {sample_name}", heading_style))
        # NO spacer here - table appears immediately after header

        species_groups = samples_dict[sample_name]

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

        # Add bookmarks for each species group (for TOC links to work)
        for species_group in sorted_groups:
            species_key = species_group.get('toplevelkey', species_group.get('key', 'unknown'))
            species_bookmark = f"species_{sanitize_bookmark_name(sample_name)}_{species_key}"
            story.append(Paragraph(f'<a name="{species_bookmark}"/>', anchor_style))

        # Collect all strains from all species groups for this sample
        all_sample_strains = []
        species_group_map = {}  # Map each strain to its species group

        for species_group in sorted_groups:
            # Get strains for this group
            group_strains = []
            for strain in species_group.get('members', []):
                # Check category filter first
                if should_include_strain(strain, args):
                    # Then check confidence threshold
                    if passes_confidence_threshold(strain, args.min_conf):
                        group_strains.append(strain)
                        species_group_map[id(strain)] = species_group
                    else:
                        # Collect low confidence strains for later
                        if strain.get("high_cons"):
                            low_confidence_strains.append({
                                'strain': strain,
                                'sample_name': sample_name,
                                'species_group': species_group
                            })

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
                sample_name=sample_name
            )
            story.append(combined_table)
        else:
            story.append(Paragraph("<i>No data available for this sample above confidence threshold</i>", styles['Italic']))
        if args.max_members is not None and args.max_members > 0:
            story.append(Paragraph(f"<i>Showing top {args.max_members} strains per group by TASS score</i>", small_style))
        story.append(Spacer(1, 0.1*inch))

    # Add footer documentation section
    story.append(Spacer(1, 0.15*inch))

    # Add anchor for Color Key
    story.append(Paragraph('<a name="color_key"/>', anchor_style))
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

    # Add anchor for Column Explanations
    story.append(Paragraph('<a name="column_explanations"/>', anchor_style))
    story.append(Paragraph("<b>Column Explanations:</b>", heading_style))
    story.append(Spacer(1, 0.03*inch))

    column_explanations = [
        "• <b>Specimen ID (Type):</b> The unique identifier for the sample, including the type of specimen (e.g., blood, tissue).",
        "• <b>Detected Organism:</b> The organism detected in the sample, which could be a bacterium, virus, fungus, or parasite.",
        "• <b>Microbial Category:</b> The classification of the organism, indicating whether it is primary, opportunistic, commensal, or potential.",
        "• <b># Reads Aligned:</b> The number of reads from the sequencing data that align to the organism's genome, indicating its presence. (%) refers to all alignments (more than 1 alignment per read can take place) for that species across the entire sample. The format is (total % of aligned reads in sample).",
        "• <b>TASS Score:</b> A metric between 0 and 100 that reflects the confidence of the organism's detection, with 100 being the highest confidence.",
        "• <b>Taxonomic ID #:</b> The taxid for the organism according to NCBI Taxonomy, which provides a unique identifier for each species. The parenthesis (if present) is the group it belongs to, usually the genus.",
        "• <b>Pathogenic Subsp/Strains:</b> Indicates specific pathogenic subspecies, serotypes, or strains, if detected in the sample. (%) indicates the percent of all aligned reads belonging to that strain.",
        "• <b>K2 Reads:</b> The number of reads classified by Kraken2, a tool for taxonomic classification of sequencing data.",
        # "• <b>HMP Percentile:</b> What percentile the abundance falls under relative to the given sample type based on HMP NCBI taxonomy classification information."
    ]

    for explanation in column_explanations:
        story.append(Paragraph(explanation, metadata_style))
        story.append(Spacer(1, 0.02*inch))

    story.append(Spacer(1, 0.1*inch))

    # Add anchor and table for Low Confidence Detections
    if low_confidence_strains:
        story.append(Paragraph('<a name="low_confidence"/>', anchor_style))
        story.append(Paragraph("<b>Low Confidence, High Consequence Detections:</b>", heading_style))
        story.append(Spacer(1, 0.03*inch))

        explanation_text = (
            f"The following strains were detected but fell below the confidence threshold "
            f"of {args.min_conf} and are listed here for reference."
        )
        story.append(Paragraph(explanation_text, metadata_style))
        story.append(Spacer(1, 0.05*inch))

        # Create low confidence table
        low_conf_table = create_low_confidence_table(
            low_confidence_strains, small_style, show_k2_column
        )
        story.append(low_conf_table)
        story.append(Spacer(1, 0.1*inch))

    # Add anchor for Additional Information
    story.append(Paragraph('<a name="additional_info"/>', anchor_style))
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
    # Build the PDF
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


def create_combined_sample_table(all_strains, species_group_map, small_style, ani_data, ani_threshold, show_ani_column, show_k2_column, taxid_to_bookmark, valid_bookmarks, sample_total_reads=0, sample_name=None):
    """
    Create a single table combining all strains from all species groups.
    Each species group's strains appear consecutively with group info on the left.

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

    ani_style = ParagraphStyle(
        'ANIStyle',
        parent=small_style,
        fontSize=6,
        leading=8
    )

    # Table headers - build dynamically based on what columns to show
    headers = ['Group Info', '', 'Strain Name', 'Reads', 'Coverage', 'TASS', 'K2 Reads', 'High ANI'] if show_k2_column and show_ani_column else (
        ['Group Info', '', 'Strain Name', 'Reads', 'Coverage', 'TASS', 'K2 Reads'] if show_k2_column else (
            ['Group Info', '', 'Strain Name', 'Reads', 'Coverage', 'TASS', 'High ANI'] if show_ani_column else
            ['Group Info', '', 'Strain Name', 'Reads', 'Coverage', 'TASS']
        )
    )

    table_data = [headers]

    # Track which species groups we've seen and their row spans
    current_species_group = None
    group_start_row = None
    span_info = []  # List of (start_row, end_row, group_summary)

    row_idx = 1  # Start after header

    for strain in all_strains:
        species_group = species_group_map[id(strain)]
        species_key = species_group.get('toplevelkey', species_group.get('key', 'unknown'))

        # Check if we're starting a new species group
        if current_species_group != species_key:
            # Save span info for previous group
            if current_species_group is not None:
                span_info.append((group_start_row, row_idx - 1))

            # Start new group
            current_species_group = species_key
            group_start_row = row_idx

            # Create group summary with colored square
            group_name = species_group.get('toplevelname', 'Unknown')
            group_tass = species_group.get('tass_score', 0)
            group_reads = species_group.get('numreads', 0)
            group_category = species_group.get('microbial_category', 'Unknown')
            group_ann_class = species_group.get('annClass', '')

            # Get color for the square
            category_color = get_category_color(group_category, group_ann_class, alpha=1.0)
            color_hex = '#{:02x}{:02x}{:02x}'.format(
                int(category_color.red * 255),
                int(category_color.green * 255),
                int(category_color.blue * 255)
            )

            group_summary = Paragraph(
                f"<link color=\"blue\"  href=\"https://www.ncbi.nlm.nih.gov/taxonomy/?term=Escherichia\">{group_name}</link> "
                f'<font color="{color_hex}">■</font><br/>'
                f"TASS: {group_tass*100:.1f}%<br/>"
                f"# Alignments: {group_reads:,.0f}",
                small_style
            )
        else:
            group_summary = ''  # Will be spanned from above

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
            )# Add star if high consequence

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
        # Build row - add columns conditionally
        row = [
            group_summary,
            indicator_text,
            strain_name,
            reads_paragraph,
            Paragraph(f"{min(100,strain.get('coverage', 0)*100):.1f}%", data_style),
            Paragraph(f"{strain.get('tass_score', 0)*100:.1f}", data_style),
        ]

        if show_k2_column:
            row.append(Paragraph(f"{strain.get('k2_reads', 0):,.0f}", data_style))

        if show_ani_column:
            row.append(high_ani_text)

        table_data.append(row)
        row_idx += 1

    # Save final group span
    if current_species_group is not None:
        span_info.append((group_start_row, row_idx - 1))

    # Create table with dynamic column widths - optimized to prevent overflow
    # Calculate available width (letter size minus margins)
    page_width = 50.5 * inch

    margins = page_width * 0.09  # 2.5% on each side

    if show_k2_column and show_ani_column:
        # All columns: need to fit 8 columns
        col_widths = [1.0*inch, 0.2*inch, 2.3*inch, 0.6*inch, 0.6*inch, 0.5*inch, 0.6*inch, 1.2*inch]
    elif show_k2_column:
        # 7 columns: no ANI
        col_widths = [1.1*inch, 0.2*inch, 2.8*inch, 0.7*inch, 0.7*inch, 0.5*inch, 0.7*inch]
    elif show_ani_column:
        # 7 columns: no K2
        col_widths = [1.1*inch, 0.2*inch, 2.5*inch, 0.7*inch, 0.7*inch, 0.5*inch, 1.3*inch]
    else:
        # 6 columns: base case
        col_widths = [1.2*inch, 0.2*inch, 3.2*inch, 0.8*inch, 0.8*inch, 0.6*inch]

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

        # Group info column styling
        ('VALIGN', (0, 1), (0, -1), 'TOP'),
        ('ALIGN', (0, 1), (0, -1), 'LEFT'),
        ('FONTSIZE', (0, 1), (0, -1), 8),
        ('LEFTPADDING', (0, 1), (0, -1), 6),
        ('RIGHTPADDING', (0, 1), (0, -1), 6),

        # Indicator column styling (star column)
        ('ALIGN', (1, 1), (1, -1), 'CENTER'),
        ('VALIGN', (1, 1), (1, -1), 'MIDDLE'),
        ('FONTSIZE', (1, 1), (1, -1), 14),

        # Grid
        ('GRID', (0, 0), (-1, -1), 1, colors.black),
        ('VALIGN', (2, 1), (2, -1), 'MIDDLE'),
        ('VALIGN', (3, 1), (-1, -1), 'MIDDLE'),

        # General body styling
        ('ALIGN', (3, 1), (-1, -1), 'CENTER'),
        ('ALIGN', (2, 1), (2, -1), 'LEFT'),

        # Padding for strain name column
        ('LEFTPADDING', (2, 1), (2, -1), 6),
        ('RIGHTPADDING', (2, 1), (2, -1), 6),
        ('TOPPADDING', (2, 1), (2, -1), 6),
        ('BOTTOMPADDING', (2, 1), (2, -1), 6),
    ]

    # Add spans for each species group
    for start_row, end_row in span_info:
        if start_row <= end_row:
            table_styles.append(('SPAN', (0, start_row), (0, end_row)))

    # Add row-specific colors
    row_idx = 1
    for strain in all_strains:
        microbial_category = strain.get('microbial_category', 'Unknown')
        ann_class = strain.get('annClass', '')

        # Debug output
        print(f"  Row {row_idx}: {strain.get('name', 'Unknown')[:40]} - Category: {microbial_category}, Class: {ann_class}")

        indicator_color = get_category_color(microbial_category, ann_class, alpha=1.0)
        table_styles.append(('BACKGROUND', (1, row_idx), (1, row_idx), indicator_color))

        row_color = get_category_color(microbial_category, ann_class, alpha=0.15)
        table_styles.append(('BACKGROUND', (2, row_idx), (-1, row_idx), row_color))

        row_idx += 1

    table.setStyle(TableStyle(table_styles))
    return table


def create_low_confidence_table(low_confidence_strains, small_style, show_k2_column):
    """
    Create a simple table for low confidence detections.

    Args:
        low_confidence_strains: List of dicts with strain, sample_name, species_group
        small_style: ParagraphStyle for small text
        show_k2_column: Whether to include K2 Reads column

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
    headers = ['Sample', 'Strain Name', 'Reads', 'TASS']
    if show_k2_column:
        headers.insert(3, 'K2 Reads')

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

        # Build row - add columns conditionally
        row = [
            Paragraph(sample_name, strain_name_style),
            strain_name,
            Paragraph(f"{strain.get('numreads', 0):,.0f}", data_style),
        ]

        if show_k2_column:
            row.append(Paragraph(f"{strain.get('k2_reads', 0):,.0f}", data_style))

        row.append(Paragraph(f"{strain.get('tass_score', 0)*100:.1f}", data_style))

        table_data.append(row)

    # Create table with dynamic column widths
    if show_k2_column:
        col_widths = [1.2*inch, 3.8*inch, 0.8*inch, 0.7*inch, 0.6*inch]
    else:
        col_widths = [1.3*inch, 4.2*inch, 0.9*inch, 0.7*inch]

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
        ('ALIGN', (2, 1), (-1, -1), 'CENTER'),  # Numbers centered
        ('ALIGN', (0, 1), (1, -1), 'LEFT'),  # Names left-aligned

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
    parser.add_argument("--show_commensals", action="store_true", required=False,
                        help="Show the commensals table")
    parser.add_argument("--show_unidentified", action="store_true", required=False,
                        help="Show the all organisms now listed as commensal or pathogen")
    parser.add_argument("--show_potentials", action="store_true", required=False,
                        help="Show the potentials table")
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
        help="Path of output file (txt)",
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

    return taxdump_dict, names_map, merged_tax_data, sample_data


if __name__ == "__main__":
    taxdump_dict, names_map, merged_tax_data, sample_data = main()
