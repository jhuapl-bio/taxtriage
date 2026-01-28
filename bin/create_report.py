#!/usr/bin/env python3

##############################################################################################
# Copyright 2022 The Johns Hopkins University Applied Physics Laboratory LLC
# All rights reserved.
# Permission is hereby granted, free of charge, to any person obtaining a copy of this
# software and associated documentation files (the "Software"), to deal in the Software
# without restriction, including without limitation the rights to use, copy, modify,
# merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
# PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
# OR OTHER DEALINGS IN THE SOFTWARE.
#

import os
import pandas as pd
import numpy as np
from pandas.api.types import CategoricalDtype

import matplotlib.pyplot as plt
from reportlab.lib import colors
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter, landscape
from distributions import import_distributions, body_site_map, make_vplot
from reportlab.graphics.shapes import Line

from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer, Image, ListFlowable, ListItem
from reportlab.platypus.flowables import HRFlowable, Flowable

from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch

from datetime import datetime
from io import StringIO
from reportlab.lib.colors import Color
from map_taxid import load_taxdump, load_names, load_merged, get_root

import argparse
def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv",
    )
    parser.add_argument(
        "-i",
        "--input",
        metavar="INPUT",
        required=True,
        nargs="+",
        default=[],
        help="Base pathogen discovery table file(s), TSV format only. Can specify more than one",
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
    parser.add_argument("--show_commensals", action="store_true", required=False,
                        help="Show the commensals table")
    parser.add_argument("--show_unidentified",   action="store_true", required=False,
                        help="Show the all organisms now listed as commensal or pathogen")
    parser.add_argument("--show_potentials",  action="store_true", required=False,
                        help="Show the potentials table")
    parser.add_argument("-m", "--missing_samples", metavar="MISSING", required=False, default=None,
                        help="Missing samples if any", nargs="+")
    parser.add_argument("-s", "--sitecol", metavar="SCOL", required=False, default='Sample Type',
                        help="Name of site column, default is body_site")
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
        "-u",
        "--output_txt",
        metavar="OUTPUT_TXT",
        required=False,
        type=str,
        help="Path of output file (txt)",
    )

    return parser.parse_args(argv)




# Function to adjust font size based on text length
def adjust_font_size(text, max_length=10, default_font_size=10, min_font_size=6):
    if len(text) > max_length:
        # Calculate new font size (simple linear reduction, could be improved)
        new_size = max(default_font_size - (len(text) - max_length) // 5, min_font_size)
        return f'<font size={new_size}>{text}</font>'
    else:
        return f'<font size={default_font_size}>{text}</font>'

# Ensure cell content that is empty is displayed as a blank space in the PDF
def format_cell_content(cell):
    # Convert NaN or None to an empty string
    if pd.isna(cell):
        return ""
    else:
        # Adjust font size based on content length
        return adjust_font_size(str(cell))

def with_alpha(color, alpha):
    """
    Return a new Color with the same RGB as `color` but different alpha.
    Works for reportlab Color objects.
    """
    return colors.Color(color.red, color.green, color.blue, alpha=alpha)

def import_data(inputfiles ):
    # Load your TSV data into a DataFrame
    # tsv_data = """
    # Name\Specimen ID\tSpecimen Type\t% Reads\t% Aligned Reads\t# Reads Aligned\tIsAnnotated\tPathogenic Sites\tType\tTaxonomic ID #\tStatus\tGini Coefficient\tMean BaseQ\tMean MapQ\tMean Coverage\tMean Depth\tAnnClass\tisSpecies\tPathogenic Subsp/Strains
    # Escherichia coli\tNasal Swab\tnasal\t100.0\t1.28\t22\tYes\tblood, urine\tCommensal\t562\testablished\t0.054\t14.0\t26.0\t2.59\t0.028\tDirect\tFalse\tCommensal Listing
    # Salmonella enterica\tNasal Swab\tnasal\t24.13\t0.49\t7\tYes\tstool\tPrimary\t28901\testablished\t0.061\t14.0\t16.3\t1.63\t0.016\tDirect\tFalse\t
    # Staphylococcus aureus\tNasal Swab\tnasal\t96.86\t52.39\t896\tYes\tblood\tCommensal\t1280\testablished\t0.63\t14.0\t59.0\t57.66\t0.85\tDirect\tFalse\tCommensal Listing
    # Pseudomonas aeruginosa\tNasal Swab\tnasal\t0.32\t0.21\t3\tYes\tabscess, blood, urine\tPrimary\t287\testablished\t0.024\t14.0\t59.7\t1.86\t0.018\tDirect\tFalse\t
    # Listeria monocytogenes\tNasal Swab\tnasal\t17.44\t11.52\t197\tYes\tblood, gut\tPrimary\t1639\testablished\t.87\t14.0\t52.4\t29.28\t0.35\tDirect\tFalse\t
    # Bacillus subtilis\tNasal Swab\tnasal\t2.75\t1.87\t32\tYes\t\tCommensal\t1423\testablished\t0.055\t14.0\t18.1\t1.75\t0.021\tDirect\tFalse\tCommensal Listing
    # Limosilactobacillus fermentum\tNasal Swab\tnasal\t15.88\t12.74\t218\tNo\t\tUnknown\t1613\t\t0.65\t14.0\t59.72\t66.56\t9.37\tNone\tTrue\t
    # Enterococcus faecalis\tNasal Swab\tnasal\t6.44\t5.55\t95\tYes\tblood, urine\tOpportunistic\t1351\testablished\t0.53\t14.0\t60.0\t53.42\t2.19\tDirect\tFalse\t
    # Saccharomyces cerevisiae\tNasal Swab\tnasal\t13.81\t13.81\t236\tNo\t\tPotential\t4932\t0.17\t13.97\t59.42\t13.099\t0.14\tDerived\tFalse\tS288C (0.1%)
    # Staphylococcus aureus\tStool\tgut\t100.0\t63.62\t1427\tYes\tblood\tCommensal\t1280\testablished\t0.083\t40.0\t15.33\t4.57\t0.047\tDirect\tFalse\tCommensal Listing
    # Neisseria gonorrhoeae\tStool\tgut\t36.15\t36.023\t808\tYes\tblood, urine\tPrimary\t485\testablished\t0.775\t40.0\t41.4\t4.55\t0.046\tDirect\tTrue\tTUM19854 (0.4%)
    # Metabacillus litoralis\tStool\tgut\t0.35\t0.35\t8\tNo\t\tUnknown\t152268\t\t0.0078\t40.0\t17.87\t0.61\t0.0061\tNone\tTrue\t
    # """.strip()
    # df = pd.read_csv(StringIO(tsv_data), sep='\t')

    # # # Simulating additional data
    # np.random.seed(42)
    # df['Gini Coefficient'] = np.random.uniform(0, 1, df.shape[0])
    # df['MeanBaseQ'] = np.random.uniform(20, 40, df.shape[0])
    # df['MeanMapQ'] = np.random.uniform(30, 60, df.shape[0])
    # df['Breadth of Coverage'] = np.random.uniform(50, 100, df.shape[0])
    # df['Depth of Coverage'] = np.random.uniform(10, 100, df.shape[0])
    if len(inputfiles) ==0:
        raise ValueError("No input files given")
    dfs = []
    for inputfile in inputfiles:
        df = pd.read_csv(inputfile, sep='\t')
        dfs.append(df)
    df = pd.concat(dfs)

    df['% Reads'] = df['% Reads'].apply(lambda x: float(x) if not pd.isna(x) else 0)
    # set # Reads Aligned as int
    df['# Reads Aligned'] = df['# Reads Aligned'].apply(lambda x: int(x) if not pd.isna(x) else 0)
    # set TASS Score as a float
    # sort the dataframe by the Sample THEN the # Reads
    # df = df.sort_values(by=["Specimen ID",  "TASS Score", "Microbial Category"], ascending=[True, False, True])
    # trim all of NAme column  of whitespace either side
    df["Detected Organism"] = df["Detected Organism"].str.strip()
    dictnames = {
        11250: "human respiratory syncytial virus B",
        12814: "human respiratory syncytial virus A",
    }
    # df['Organism'] = df["Detected Organism"]
    # set if putative to it with *  in Detected organism using lambda x
    df['Detected Organism'] = df.apply(lambda x: f'{x["Detected Organism"]}*' if x['Status'] == 'putative' else x["Detected Organism"], axis=1)
    df['Detected Organism'] = df.apply(lambda x: f'{x["Detected Organism"]}°' if x['AnnClass'] == 'Derived' else x["Detected Organism"], axis=1)
    df['Detected Organism'] = df.apply(lambda x: f'{x["Detected Organism"]}' if x['High Consequence'] == True else x["Detected Organism"], axis=1)
    df["Detected Organism"] = df[["Detected Organism", 'Taxonomic ID #']].apply(lambda x: dictnames[x['Taxonomic ID #']] if x['Taxonomic ID #'] in dictnames else x["Detected Organism"], axis=1)
    # replace all NaN with ""
    df = df.fillna("")
    return df
def extract_reads(value):
    """
    Extract leading integer read count from strings like:
    '142 (63.62%)' or '8 (0.35% - 0.12%)'
    """
    if pd.isna(value):
        return 0
    s = str(value).strip()
    if s == "":
        return 0
    try:
        # take everything before first space or '('
        for sep in [' ', '(']:
            if sep in s:
                s = s.split(sep)[0]
                break
        return int(float(s))
    except Exception:
        return 0

def prepare_three_layer_table(df, sample_col, group_col, plot_dict, names_map=None, include_headers=True, columns=None):
    """
    Build a table with three levels:
      - Sample header (spans full width)
      - Group header (spans full width) inside each sample
      - Member rows (one per df row)

    Ensures every table row has the same number of columns (important when some rows
    include a plot image cell and others do not).
    Returns: data (list of rows), style_cmds (list of TableStyle commands)
    """
    data = []
    style_cmds = []

    # Work on a copy so we don't mutate caller's df
    df_proc = df.copy()

    if 'K2 Reads' in df_proc.columns:
        df_proc['K2 Reads'] = df_proc['K2 Reads'].apply(lambda x: int(x) if not pd.isna(x) and x != "" else 0)

    if columns is None:
        columns = [c for c in df_proc.columns.values if c not in (sample_col, group_col, 'plot')]

    # Determine whether we need an explicit plot column (so all rows have same width)
    has_plot_column = len(plot_dict.keys()) > 0

    # number of columns in member area: len(columns) + (1 plot col if present)
    member_cols_count = len(columns) + (1 if has_plot_column else 0)
    # total columns including the left "sample/group" placeholder
    num_cols = 1 + member_cols_count

    # small empty Flowable to use for blank cells (avoids plain strings)
    empty_cell = Paragraph('', small_font_style)
    empty_cell_normal = Paragraph('', styles['Normal'])  # for header-like blanks if needed

    # Header (column titles) - single row above everything
    if include_headers:
        header_cells = [Paragraph('', styles['Normal'])]
        for col in columns:
            header_cells.append(Paragraph(f'<b>{col}</b>', styles['Normal']))
        if has_plot_column:
            header_cells.append(Paragraph(f'<b>Plot</b>', styles['Normal']))
        # Ensure header row length equals num_cols (use Flowables)
        while len(header_cells) < num_cols:
            header_cells.append(empty_cell)
        data.append(header_cells)
        style_cmds.append(('BACKGROUND', (0,0), (-1,0), colors.slategray))
        style_cmds.append(('TEXTCOLOR', (0,0), (-1,0), colors.whitesmoke))
        style_cmds.append(('ALIGN', (0,0), (-1,0), 'CENTER'))
        start_row = 1
    else:
        start_row = 0

    # helper to find reads column
    possible_reads_cols = ['Reads Aligned', '# Reads Aligned', '# Reads Aligned to Sample', 'K2 Reads', 'Quant']
    reads_col = next((c for c in possible_reads_cols if c in df_proc.columns), None)

    # Build ordered list of samples (preserve df order)
    sample_groups = []
    for sample_key, sample_df in df_proc.groupby(sample_col, sort=False):
        sample_groups.append((sample_key, sample_df))

    current_row = start_row
    for s_idx, (sample_key, sample_df) in enumerate(sample_groups):
        # sample-level stats
        sample_label = sample_key if sample_key and str(sample_key).strip() != "" else "Unknown Sample"
        if reads_col:
            # use extract_reads helper to robustly parse numeric-like strings
            sample_reads = int(sample_df[reads_col].apply(extract_reads).sum())
        else:
            sample_reads = 0
        sample_count = sample_df.shape[0]

        # add sample header (span full width)
        sample_header_text = f"<b>{sample_label}</b><br/>{sample_count} hits — {sample_reads:,} alignments"
        sample_para = Paragraph(sample_header_text, styles['Normal'])
        sample_row = [sample_para] + [empty_cell_normal] * (num_cols - 1)
        # ensure length (all Flowables)
        while len(sample_row) < num_cols:
            sample_row.append(empty_cell_normal)
        data.append(sample_row)
        style_cmds.append(('SPAN', (0, current_row), (num_cols-1, current_row)))
        # Neutral header background
        style_cmds.append(('BACKGROUND', (0, current_row), (num_cols-1, current_row), colors.grey))
        style_cmds.append(('TEXTCOLOR', (0, current_row), (num_cols-1, current_row), colors.black))
        style_cmds.append(('ALIGN', (0, current_row), (num_cols-1, current_row), 'LEFT'))
        style_cmds.append(('VALIGN', (0, current_row), (num_cols-1, current_row), 'TOP'))
        style_cmds.append(('LEFTPADDING', (0, current_row), (num_cols-1, current_row), 6))
        style_cmds.append(('RIGHTPADDING', (0, current_row), (num_cols-1, current_row), 6))


        current_row += 1

        # Build ordered list of groups for this sample
        group_list = []
        for group_key, group_df in sample_df.groupby(group_col, sort=False):
            group_list.append((group_key, group_df))

        for g_idx, (group_key, group_df) in enumerate(group_list):
            group_taxid = str(group_key) if group_key is not None else ""
            if names_map and group_taxid in names_map:
                group_label = f"{names_map[group_taxid]} ({group_taxid})"
            else:
                group_label = group_taxid if group_taxid and group_taxid.lower() != "unknown" else group_taxid

            # group-level stats
            if reads_col:
                group_reads = int(group_df[reads_col].apply(extract_reads).sum())
            else:
                group_reads = 0
            group_count = group_df.shape[0]
            high_con_sequence_count = group_df[group_df['High Consequence'] == True].shape[0]

            group_header_text = f"<b>{group_label}</b><br/>{group_count} hits ({high_con_sequence_count} High Consequence) — {group_reads:,} alignments"
            group_para = Paragraph(group_header_text, styles['Normal'])
            group_row = [group_para] + [empty_cell_normal] * (num_cols - 1)
            while len(group_row) < num_cols:
                group_row.append(empty_cell_normal)
            data.append(group_row)

            # style for group header
            style_cmds.append(('SPAN', (0, current_row), (num_cols-1, current_row)))
            # Neutral group header (span full width) so it doesn't get mistaken for coloring the data rows
            style_cmds.append(('BACKGROUND', (0, current_row), (num_cols-1, current_row), colors.lightgrey))
            style_cmds.append(('TEXTCOLOR', (0, current_row), (num_cols-1, current_row), colors.black))
            style_cmds.append(('ALIGN', (0, current_row), (num_cols-1, current_row), 'CENTER'))
            style_cmds.append(('VALIGN', (0, current_row), (num_cols-1, current_row), 'MIDDLE'))
            style_cmds.append(('LEFTPADDING', (0, current_row), (num_cols-1, current_row), 6))
            style_cmds.append(('RIGHTPADDING', (0, current_row), (num_cols-1, current_row), 6))


            current_row += 1

            # add member rows under this group
            star_style = ParagraphStyle(
                name='StarStyle',
                parent=styles['Normal'],
                alignment=1,          # 1 = CENTER
                fontSize=10,
                leading=10
            )
            for idx, row in group_df.iterrows():
                # LEFT column: show a star if High Consequence else keep blank, keep it center aligned
                if row.get('High Consequence', False) in (True, 'True', 'true', 1):
                    left_cell = Paragraph('★', star_style)
                else:
                    left_cell = empty_cell

                member_cells = [left_cell]

                # Fill data columns as Paragraph flowables
                for col in columns:
                    cell_val = row.get(col, "")
                    member_cells.append(Paragraph(format_cell_content(cell_val), small_font_style))

                # append plot image cell or empty Flowable to ensure constant width
                if has_plot_column:
                    plot_key = (row.get('Detected Organism'), row.get('Specimen Type'))
                    if plot_key in plot_dict:
                        plot_image = Image(plot_dict[plot_key])
                        plot_image.drawHeight = 0.4 * inch
                        plot_image.drawWidth = 0.9 * inch
                        member_cells.append(plot_image)
                        plot_is_empty = False
                    else:
                        member_cells.append(empty_cell)
                        plot_is_empty = True
                else:
                    plot_is_empty = False  # no plot column at all

                # Ensure member_cells length == num_cols (all Flowables)
                while len(member_cells) < num_cols:
                    member_cells.append(empty_cell)

                data.append(member_cells)

                # determine color for the member data cells (based on Microbial Category / AnnClass)
                val = row.get('Microbial Category', "")
                derived = row.get('AnnClass', "")

                if "Primary" in str(val) and derived == "Direct":
                    color = colors.lightcoral
                elif "Primary" in str(val):
                    color = colors.HexColor('#fab462')
                elif "Commensal" in str(val):
                    color = colors.lightgreen
                elif "Opportunistic" in str(val):
                    color = colors.HexColor('#ffe6a8')
                elif "Potential" in str(val):
                    color = colors.lightblue
                else:
                    color = colors.white

                # Determine indices:
                # left spacer column index = 0
                # data columns indices = 1 .. (num_cols-1)  if no plot column
                # if plot column exists, plot column index = num_cols-1, so data columns end at num_cols-2
                if has_plot_column:
                    data_start_col = 1
                    data_end_col = num_cols - 2
                    plot_col_index = num_cols - 1
                else:
                    data_start_col = 1
                    data_end_col = num_cols - 1
                    plot_col_index = None
                # Determine indices
                data_start_col = 1
                data_end_col = num_cols - 1 if plot_col_index is None else num_cols - 2


                # Semi-transparent version for data cells
                faded_color = with_alpha(color, 0.22)
                # Apply background color ONLY to the data columns (not left spacer, not plot cell)
                # if data_end_col >= data_start_col:
                #     style_cmds.append(('BACKGROUND', (data_start_col, current_row), (data_end_col, current_row), color))
                style_cmds.append(('BACKGROUND', (0, current_row), (0, current_row), color))  # left spacer cell colored very lightly
                # Data cells: same color but faded
                if data_end_col >= data_start_col:
                    style_cmds.append(('BACKGROUND', (data_start_col, current_row), (data_end_col, current_row), faded_color))

                # center the star (if present) in that left cell
                style_cmds.append(('ALIGN', (0, current_row), (0, current_row), 'CENTER'))
                style_cmds.append(('VALIGN', (0, current_row), (0, current_row), 'MIDDLE'))
                # If plot column exists and plot cell is empty, color just the plot cell subtly
                if plot_col_index is not None and plot_is_empty:
                    style_cmds.append(('BACKGROUND', (plot_col_index, current_row), (plot_col_index, current_row), colors.HexColor("#f2f2f2b8")))
                # Center data cells (but do NOT color them)
                if data_end_col >= data_start_col:
                    style_cmds.append(('ALIGN', (data_start_col, current_row), (data_end_col, current_row), 'CENTER'))
                    style_cmds.append(('VALIGN', (data_start_col, current_row), (data_end_col, current_row), 'TOP'))

                # Keep entire-row vertical alignment and subtle divider
                style_cmds.append(('VALIGN', (0, current_row), (-1, current_row), 'TOP'))
                style_cmds.append(('LINEBELOW', (0, current_row), (-1, current_row), 0.25, colors.lightgrey))

                current_row += 1

            # Only add a spacer after the group if it's not the last group in the sample
            if g_idx < (len(group_list) - 1):
                spacer_cells = [empty_cell] * num_cols
                data.append(spacer_cells)
                current_row += 1

        # Only add a sample spacer if it's not the last sample overall
        if s_idx < (len(sample_groups) - 1):
            big_spacer = [empty_cell] * num_cols
            data.append(big_spacer)
            current_row += 1

    return data, style_cmds


def split_df(df_full):
    # Filter DataFrame for IsAnnotated == 'Yes' and 'No'
    # df_ = df_full[df_full['IsAnnotated'] == 'Yes'].copy()
    # append (taxid) from taxid column to Detected Organism

    df_yes = df_full[~df_full['Microbial Category'].isin([ 'Unknown', 'N/A', np.nan, "Commensal", "Potential" ] ) ].copy()
    df_opp = df_full[df_full['Microbial Category'].isin([ "Potential"])].copy()
    df_comm = df_full[df_full['Microbial Category'].isin(['Commensal'])].copy()
    df_unidentified = df_full[(df_full['Microbial Category'].isin([ 'Unknown', 'N/A', np.nan, ""] ))].copy()
    df_yes.reset_index(drop=True, inplace=True)
    df_opp.reset_index(drop=True, inplace=True)
    return df_yes, df_opp, df_comm, df_unidentified


def safe_get(row, column_name):
    """
    Safely access a value from a Pandas row (Series) and handle missing or empty values.

    Parameters:
    - row: Pandas Series representing a row.
    - column_name: The column name to access.

    Returns:
    - The value if present, stripped of leading/trailing spaces if it's a string.
    - None if the column is missing or the value is empty.
    """
    # Check if the column exists in the row
    if column_name in row.index:
        val = row[column_name]
        # Handle string values by stripping spaces
        if isinstance(val, str):
            val = val.strip()
            return val if val else None  # Return None if the string is empty after stripping
        return val  # Return the value directly if it's not a string
    return None  # Return None if the column is not present


# Custom styles for title and subtitle with right alignment
left_align_style = ParagraphStyle(
    name='leftAlign',
    parent=getSampleStyleSheet()['Normal'],
    alignment=0,  # 2 is for right alignment
    # fontSize=12,
    spaceAfter=10,
)
subtitle_style = ParagraphStyle(
    name='leftAlignSubtitle',
    parent=getSampleStyleSheet()['Normal'],
    alignment=0,  # Right align
    # fontSize=12,
    spaceAfter=10,
)

title_style = ParagraphStyle(
    name='leftAlignTitle',
    parent=getSampleStyleSheet()['Title'],
    alignment=0,  # Right align
    # fontSize=14,
    spaceAfter=10,
)
styles = getSampleStyleSheet()
small_font_style = ParagraphStyle(name='SmallFont', parent=styles['Normal'], fontSize=2)
normal_style = styles['Normal']

def return_table_style(df=pd.DataFrame(), color_pathogen=False):
    """
    Basic table style: header/grid/valign. Actual body coloring is applied by the
    table builder (prepare_three_layer_table) because it knows the exact row indices.
    """
    table_style = TableStyle([
        ('BACKGROUND', (0,0), (-1,0), colors.gray),
        ('TEXTCOLOR', (0,0), (-1,0), colors.whitesmoke),
        ('GRID', (0,0), (-1,-1), 1, colors.black),
        ('VALIGN', (0,0), (-1,-1), 'TOP'),
        ('LEFTPADDING', (0,0), (-1,-1), 4),
        ('RIGHTPADDING', (0,0), (-1,-1), 4),
    ])
    return table_style

def draw_vertical_line(canvas, doc):
        """
        Draw a vertical line 5% from the left of the page, starting 10% down from the top
        and ending 10% up from the bottom.
        """
        page_width, page_height = letter
        line_x = 0.05 * page_width
        start_y = 0.1 * page_height
        end_y = page_height - (0.1 * page_height)
        canvas.saveState()
        canvas.setStrokeColor(colors.black)
        canvas.setLineWidth(1)
        canvas.line(line_x, start_y, line_x, end_y)
        canvas.restoreState()
def make_table(data, table_style=None, col_widths=None):
    """
    Create a ReportLab Table with optional explicit column widths.
    """
    if col_widths:
        table = Table(data, repeatRows=1, colWidths=col_widths)
    else:
        table = Table(data, repeatRows=1)

    if table_style:
        table.setStyle(table_style)
    return table



def create_report(
    output,
    df_identified,
    df_potentials,
    df_unidentified,
    df_commensals,
    df_high_cons_low_conf,
    plotbuffer,
    version=None,
    missing_samples=None,
    min_conf=None,
    names_map = {}
):

    # PDF file setup
    pdf_file = output
    doc = SimpleDocTemplate(pdf_file, pagesize=landscape(letter))
    # Placeholder values for version and date
    #### Section to style things up a bit
    # Set the left margin to 10% of the width of a letter size page (8.5 inches)
    # Set custom margins based on percentage of the page size
    left_margin = 0.1 * letter[0]  # 10% of the width of a letter page (landscape width)
    right_margin = left_margin / 5  # 1/5th of the left margin for right margin
    top_margin = bottom_margin = 0.1 * letter[1]  # 10% of the height of a letter page

    # Modify the doc setup to include custom margins and landscape mode
    doc = SimpleDocTemplate(
        pdf_file,
        pagesize=landscape(letter),  # Ensure the document is landscape
        leftMargin=left_margin,
        rightMargin=right_margin,
        topMargin=top_margin,
        bottomMargin=bottom_margin
    )
    # version = "1.3.2"  # Example version
    if not version:
        version = "Local Build"
    # get datetime of year-mont-day hour:min
    date = datetime.now().strftime("%Y-%m-%d %H:%M")  # Current date
    # sort df_identified by TASS Score
    # filter out so only Class is PAthogen
    # df_identified = df_identified.sort_values(by=['Specimen ID', '# Reads Aligned'], ascending=[False, True])
    # df_potentials = df_potentials.sort_values(by=['Specimen ID', '# Reads Aligned'], ascending=[False, False])
    df_identified_paths = df_identified
    df_identified_others = df_commensals
    # df_identified_others = df_identified[df_identified['Class'] != 'Pathogen']
    # df_unidentified = df_unidentified.sort_values(by=['Specimen ID', '# Reads Aligned'], ascending=[False, False])
    elements = []

    ##########################################################################################
    ##### Section to make the Top Table - all annotated commensal or otherwise
    title = Paragraph("Organism Discovery Report", title_style)
    # Add the title and subtitle
    if missing_samples:
        samples_missing = f"Missing Samples due to filtering, failed step(s), or non-presence of organisms: <b>{', '.join(missing_samples)}</b><br></br> "
    else:
        samples_missing = ""
    subtitle = Paragraph(f"This report was generated using TaxTriage \
        <b>{version}</b> on <b>{date}</b> and is derived from an in development spreadsheet of \
        human-host pathogens.<br></br><br></br> \
        {samples_missing}\
        Samples with confidence score below <b>{str(min_conf)}</b> were filtered out. \
        ",
        subtitle_style
    )
    elements.append(title)
    elements.append(subtitle)
    extra_text = Paragraph(f"★ denotes high consequence pathogens")
    elements.append(extra_text)
    elements.append(Spacer(1, 12))
    if not df_identified_paths.empty:
        columns_yes = df_identified_paths.columns.values
        # print only rows in df_identified with Gini Coeff above 0.2
        columns_yes = [
            "Detected Organism",
            "# Reads Aligned",
            "TASS Score",
            "Taxonomic ID #", "Pathogenic Subsp/Strains",
            "Coverage",
            "K2 Reads"
        ]
        # check if all K2 reads column are 0 or nan
        if df_identified_paths['K2 Reads'].sum() == 0:
            columns_yes = columns_yes[:-1]
        if "MicrobeRT Probability" in df_identified_paths.columns.values:
            columns_yes.insert(4, "MicrobeRT Probability")
        # Now, call prepare_data_with_headers for both tables without manually preparing headers

        page_width, _ = landscape(letter)

        # total usable width (match your document margins)
        usable_width = page_width * 0.88  # ~matches your left/right margins

        # Make left spacer column narrow
        left_col_width = 0.04 * usable_width   # ~4% of table width

        # Remaining width split across real data columns
        num_data_cols = len(columns_yes) + (1 if len(plotbuffer) > 0 else 0)
        remaining_width = usable_width - left_col_width
        data_col_width = remaining_width / num_data_cols

        col_widths = [left_col_width] + [data_col_width] * num_data_cols

        data_yes, style_cmds = prepare_three_layer_table(
            df_identified_paths,
            sample_col='Specimen ID (Type)',
            group_col='Group',
            plot_dict=plotbuffer,
            names_map=names_map,
            include_headers=True,
            columns=columns_yes
        )
        table_style = return_table_style(df_identified_paths, color_pathogen=False)
        for cmd in style_cmds:
            table_style.add(*cmd)
        table = make_table(data_yes, table_style=table_style, col_widths=col_widths)


        elements.append(table)
        elements.append(Spacer(1, 12))  # Space between tables
    else:
        # Make a large message saying that no pathogen was found above min confidence
        title = Paragraph("No Primary Pathogens Found Above TASS Score Threshold", title_style)
        subtitle = Paragraph(f"No pathogens were found above the minimum confidence score of <b>{str(min_conf)}</b> in the given samples. <br></br>Consider lowering the confidence score OR looking at the <b>all.organisms.report.txt</b> file or <b>MultiQC</b> report HTML file for more information.", subtitle_style)
        elements.append(title)
        elements.append(subtitle)

    ##########################################################################################

    if min_conf:
        subtitle = Paragraph(f"See the <b>`report/all.organisms.report.txt`</b> file for a full list of everything identified", subtitle_style)
        elements.append(subtitle)

    # Adding regular text

    styles = getSampleStyleSheet()
    # Adding subtext (you can adjust the style to make it look like subtext)
    subtext_style = styles["BodyText"]
    subtext_style.fontSize = 10  # Smaller font size for subtext
    subtext_style.leading = 12
    subtext_para = Paragraph("Organisms marked with * are putative and have relatively lower references listing their annotations as a pathogen in the given sample types. Classifications of pathogens are described as:", subtext_style)
    elements.append(subtext_para)
    elements.append(Spacer(1, 12))  # Space between tables
    bullet_list_items = [
        "This is an unannotated organism or commensal organism for the given site.",
        "Primary Pathogen annotated in sample type(s) listed.",
        "Primary Pathogen annotated in sample type(s) other than your listed one.",
        "Opportunistic Pathogen",
        "Potential Pathogen",
    ]

    # Create a list of bullet items with specified colors
    bullet_colors = [colors.lightgreen,   colors.coral, '#fab462' , '#ffe6a8', colors.lightblue  ]
    style = styles['Normal']

    # Create custom ListItems with colored bullets
    custom_list_items = [
        ListItem(Paragraph(item, style), bulletBorder=colors.black,  bulletColor=bullet_colors[idx], bulletFontSize=9 )
        for idx, item in enumerate(bullet_list_items)
    ]

    # Create the ListFlowable
    bullet_list = ListFlowable(
        custom_list_items,
        start="square",

        bulletType='bullet'  # '1' for numbered list
    )
    # add horizontal line in reportlab
    elements.append(bullet_list)
    elements.append(Spacer(1, 12))  # Space between tables
    # Generate the explanation paragraph and append it to the elements
    # After adding your table to the elements
    subtext_style = styles['Normal']

    # Create an HRFlowable for the horizontal line
    horizontal_line = HRFlowable(width="100%", thickness=1, color=colors.black, spaceBefore=12, spaceAfter=12)

    # Create explanatory bullets with better formatting
    bullet_list_items = [
        "Primary: Exposure to the agent generally results in a diseased state in both immunocompromised and immunocompetent individuals.",
        "Opportunistic: Exposure to the agent causes a diseased state under certain conditions, including immunocompromised status, wound infections, and nosocomial infections.",
        "Commensal: Organisms typically found in the human microbiota.",
        "Potential: Organisms that have been associated with disease states but are not extensively studied.",
        "°: Indicates a pathogenic subspecies/serotype/strain/etc with the same name as the species listed, just different taxids."
    ]

    bullet_list = ListFlowable(
        [ListItem(Paragraph(item, subtext_style)) for item in bullet_list_items],
        bulletType='bullet',
        start='circle',
        # size of bullet
        bulletFontSize=9
    )
    elements.append(bullet_list)

    elements.append(horizontal_line)

    elements.append(Spacer(1, 12))

    # Define the column explanations
    column_explanations = [
        "Specimen ID (Type): The unique identifier for the sample, including the type of specimen (e.g., blood, tissue).",
        "Detected Organism: The organism detected in the sample, which could be a bacterium, virus, fungus, or parasite.",
        "Microbial Category: The classification of the organism, indicating whether it is primary, opportunistic, commensal, or potential.",
        "# Reads Aligned: The number of reads from the sequencing data that align to the organism's genome, indicating its presence. (%) refers to all alignments (more than 1 alignment per read can take place) for that species across the entire sample. The format is (total % of aligned reads in sample).",
        "TASS Score: A metric between 0 and 1 that reflects the confidence of the organism's detection, with 1 being the highest confidence.",
        "Taxonomic ID #: The taxid for the organism according to NCBI Taxonomy, which provides a unique identifier for each species. The parenthesis (if present) is the group it belongs to, usually the genus.",
        "Pathogenic Subsp/Strains: Indicates specific pathogenic subspecies, serotypes, or strains, if detected in the sample. (%) indicates the percent of all aligned reads belonging to that strain.",
        "K2 Reads: The number of reads classified by Kraken2, a tool for taxonomic classification of sequencing data."
        "HMP Plot: What percentile the abundance falls under relative to the given sample type based on Healthy Human Subject NCBI taxonomy classification information"
    ]

    # Create bullet points for each column explanation
    bullet_list_items = [
        Paragraph(item, subtext_style) for item in column_explanations
    ]

    # Create the bullet list
    bullet_list = ListFlowable(
        [ListItem(item) for item in bullet_list_items],
        bulletType='bullet',
        start='circle',
        bulletFontSize=9
    )

    # Append the bullet list to the elements
    elements.append(bullet_list)

    # Add the horizontal line to the elements
    elements.append(horizontal_line)

    paragraph_style = styles['Normal']

    # Create a paragraph with a clickable URL
    url = "https://github.com/jhuapl-bio/taxtriage/blob/main/docs/usage.md#confidence-scoring"
    text_with_url = f'Please visit our <a href="{url}"><b><font color="blue">DOCUMENTATION PAGE</font></b></a> for more information on how confidence is calculated.'

    # Add the paragraph with the URL
    elements.append(Paragraph(text_with_url, paragraph_style))



    subtext_para = Paragraph("The following information highlights the description for the color combinations for each organism class in the annotated table(s)", subtext_style)
    elements.append(subtext_para)
    elements.append(Spacer(1, 12))
    subtext_para = Paragraph("Please see the relevant Discovery Analysis txt file for low confidence annotations that were not present in the pdf", subtext_style)
    elements.append(subtext_para)
    elements.append(Spacer(1, 12))

    subtext_para = Paragraph("Read amounts are represented as the <b>total number of aligned reads</b> of sufficient mapping quality <b>(% aligned for all reads in sample)</b>", subtext_style)
    elements.append(subtext_para)

    elements.append(horizontal_line)
    ##########################################################################################
    if not df_high_cons_low_conf.empty:
        # print only rows in df_identified with Gini Coeff above 0.2
        columns_yes = [
                       "Detected Organism",
                       'TASS Score',
                       "# Reads Aligned", "Taxonomic ID #", "Coverage",
                       "K2 Reads"]
        # check if all K2 reads column are 0 or nan
        if df_identified_paths['K2 Reads'].sum() == 0:
            columns_yes = columns_yes[:-1]
        if "MicrobeRT Probability" in df_high_cons_low_conf.columns.values:
            columns_yes.insert(4, "MicrobeRT Probability")
        # if all of Group is Unknown, then remove it from list

        # Now, call prepare_data_with_headers for both tables without manually preparing headers
        # call new grouped builder — pass names_map (received into create_report)
        data_hc, style_cmds = prepare_three_layer_table(
            df_high_cons_low_conf,
            sample_col='Specimen ID (Type)',
            group_col='Group',
            plot_dict={},
            names_map=names_map,
            include_headers=True,
            columns=columns_yes
        )




        # if data shape is >=1 then append, otherwise make text saying it is empty
        if df_high_cons_low_conf.shape[0] >= 1:
            table_style = return_table_style(df_high_cons_low_conf, color_pathogen=False)
            for cmd in style_cmds:
                table_style.add(*cmd)
            table = make_table(data_hc, table_style=table_style, col_widths=col_widths)
            # # Add the title and subtitle
            title = Paragraph("SUPPLEMENTARY: High Consequence Low Confidence", title_style)
            subtitle = Paragraph(f"These were identified as high consequence pathogens but with low confidence. The below list of microorganisms represent pathogens of heightened concern, to which reads mapped.  The confidence metrics did not meet criteria set forth to be included in the above table; however, the potential presence of these organisms should be considered for biosafety, follow-up diagnostic testing (if clinical presentation warrants), and situational awareness purposes.", subtitle_style)
            elements.append(title)
            elements.append(subtitle)
            elements.append(Spacer(1, 12))
            elements.append(table)
        else:
            elements.append(Paragraph("No High Consequence Low Confidence Organisms", subtext_style))

        elements.append(Spacer(1, 12))

    ##########################################################################################
    #### Table on opportunistic pathogens
    if not df_potentials.empty:
        columns_opp = ["Detected Organism",
                       "# Reads Aligned",
                       "TASS Score", "Taxonomic ID #",
                       "Pathogenic Subsp/Strains", "Coverage", "K2 Reads"
                       ]
        if df_potentials['K2 Reads'].sum() == 0:
            columns_opp = columns_opp[:-1]
        if "MicrobeRT Probability" in df_potentials.columns.values:
            columns_opp.insert(4, "MicrobeRT Probability")
        data_pot, style_cmds = prepare_three_layer_table(
            df_potentials,
            sample_col='Specimen ID (Type)',
            group_col='Group',
            plot_dict=plotbuffer,  # if you don't want plots here
            names_map=names_map,
            include_headers=True,
            columns=columns_yes
        )
        table_style = return_table_style(df_potentials, color_pathogen=True)
        for cmd in style_cmds:
            table_style.add(*cmd)
        table = make_table(data_pot, table_style=table_style, col_widths=col_widths)


        # Add the title and subtitle
        Title = Paragraph("Low Potential Pathogens", title_style)
        elements.append(Title)
        elements.append(Spacer(1, 12))
        elements.append(table)
        elements.append(Spacer(1, 12))  # Space between tables

    # ################################################################################################
    ### Table on commensals
    if not df_identified_others.empty:
        columns_yes = df_identified_others.columns.values
        # print only rows in df_identified with Gini Coeff above 0.2
        columns_yes = [
                       "Detected Organism",
                       "# Reads Aligned", "TASS Score", "Taxonomic ID #", "Coverage",
                        "K2 Reads"]
        # check if all K2 reads column are 0 or nan
        if df_identified_paths['K2 Reads'].sum() == 0:
            columns_yes = columns_yes[:-1]
        if "MicrobeRT Probability" in df_identified_others.columns.values:
            columns_yes.insert(4, "MicrobeRT Probability")
        # if all of Group is Unknown, then remove it from list
        # Now, call prepare_data_with_headers for both tables without manually preparing headers
        data_comm, style_cmds = prepare_three_layer_table(
            df_identified_others,
            sample_col='Specimen ID (Type)',
            group_col='Group',
            plot_dict=plotbuffer,
            names_map=names_map,
            include_headers=True,
            columns=columns_yes
        )
        table_style = return_table_style(df_identified_others, color_pathogen=False)
        for cmd in style_cmds:
            table_style.add(*cmd)
        table = make_table(data_comm, table_style=table_style, col_widths=col_widths)

        title = Paragraph("Commensals", title_style)
        subtitle = Paragraph(f"These were identified & were listed as a commensal directly", subtitle_style)
        elements.append(title)
        elements.append(subtitle)
        elements.append(Spacer(1, 12))  # Space between tables

        elements.append(table)
        elements.append(Spacer(1, 12))  # Space between tables


    elements.append(Spacer(1, 12))
    if not df_unidentified.empty:
    #     ##########################################################################################
        ### Section to Make the "Unannotated" Table
        second_title = "Unannotated Organisms"
        second_subtitle = "The following table displays the unannotated organisms and their alignment statistics. Be aware that this is the exhaustive list of all organisms (species only) contained within the samples that had atleast one read aligned"
        elements.append(Paragraph(second_title, title_style))
        elements.append(Paragraph(second_subtitle, subtitle_style))

        columns_no = ['Detected Organism','# Reads Aligned', "TASS Score", "Coverage", "K2 Reads"]
        data_no, style_cmds = prepare_three_layer_table(
            df_unidentified,
            sample_col='Specimen ID (Type)',
            group_col='Group',
            plot_dict=plotbuffer,
            names_map=names_map,
            include_headers=True,
            columns=columns_no
        )
        table_style = return_table_style(df_unidentified, color_pathogen=False)
        for cmd in style_cmds:
            table_style.add(*cmd)
        table_no = make_table(data_no, table_style=table_style, col_widths=col_widths)
        elements.append(table_no)
    elements.append(Spacer(1, 12))  # Space between tables

    ##########################################################################################
    ##### Add equations HERE - Work in progress ##################################################

    ##########################################################################################
    #####  Build the PDF
    # Adjust the build method to include the draw_vertical_line function
    doc.build(elements, onFirstPage=draw_vertical_line, onLaterPages=draw_vertical_line)

    print(f"PDF generated: {pdf_file}")


def main():
    args = parse_args()
    taxdump_dict, names_map, merged_data = {}, {}, {}
    if args.taxdump:
        # load the taxdump file
        if os.path.exists(f"{args.taxdump}/nodes.dmp"):
            taxdump_dict = load_taxdump(f"{args.taxdump}/nodes.dmp")
        if os.path.exists(f"{args.taxdump}/names.dmp"):
            names_map = load_names(f"{args.taxdump}/names.dmp")
        if os.path.exists(f"{args.taxdump}/merged.dmp"):
            merged_data = load_merged(f"{args.taxdump}/merged.dmp")
    df_full = import_data(args.input)
    # Set tass score as a flaot
    # fill High Consequence with False if it is NaN
    df_full['High Consequence'].fillna(False, inplace=True)
    df_full['Reads Aligned'] = df_full['# Reads Aligned'].apply(lambda x: int(x) if not pd.isna(x) else 0)
    df_full['Coverage'] = df_full['Coverage'].apply(lambda x: f"{100*x:.0f}%" if not pd.isna(x) else 0)
    # fill empty string with False for high consequence
    df_full['High Consequence'] = df_full['High Consequence'].apply(lambda x: False if x == "" else x)
    df_full['TASS Score'] = df_full['TASS Score'].apply(lambda x: f"{100*x:.0f}" if not pd.isna(x) else 0)
    df_full["Group"] = df_full["Taxonomic ID #"].apply(lambda x: get_root(x, args.rank, taxdump_dict))
    # for Detected Organism, add (Group) if it is not null or not Unknown
    # df_full["Taxonomic ID #"] = df_full.apply(lambda x: f"{x['Taxonomic ID #']}, axis=1)
    # 1) create the rank
    # If your column is actually booleans + NaN, you can do:
    # fill all None for Group as "Unknown"
    df_full['Group'].fillna('Others', inplace=True)
    # 2) sort by rank desc, then TASS Score desc
    df_full.sort_values(
        by=['High Consequence', 'TASS Score'],
        ascending=[False, False],
        inplace=True
    )
    # if MicrobeRT Probability is in df_full columns then set to float and 2 decimals

    if "MicrobeRT Probability" in df_full.columns:
        # Replace empty strings with NaN
        # if all is empty strings or empty rows in total then drop it
        if df_full["MicrobeRT Probability"].replace("", np.nan).isna().all():
            df_full.drop(columns=["MicrobeRT Probability"], inplace=True)
        else:
            df_full["MicrobeRT Probability"].replace("", np.nan, inplace=True)

            # Convert to float where possible
            df_full["MicrobeRT Probability"] = pd.to_numeric(df_full["MicrobeRT Probability"], errors="coerce")

            # Format to 2 decimals, keep NaN as is
            df_full["MicrobeRT Probability"] = df_full["MicrobeRT Probability"].apply(
                lambda x: f"{x:.2f}" if not pd.isna(x) else np.nan
            )
            # reformat nan to empty string
            df_full["MicrobeRT Probability"] = df_full["MicrobeRT Probability"].replace(np.nan, "")
    if args.output_txt:
        # write out the data to a txt file
        # sort df_full on "TASS Score"
        df_full = df_full.reset_index(drop=False)
        df_full.to_csv(args.output_txt, sep="\t", index=True, index_label="Index")
    # change column "id" in avbundance_data to "tax_id" if args.type is "Detected Organism"
    # df_full = df_full.rename(columns={args.id_col: args.type})
    df_full = df_full.rename(columns={args.sitecol: 'body_site'})
    df_full = df_full.rename(columns={args.abundance_col: 'abundance'})
    df_full = df_full.dropna(subset=[args.type])
    # df_identified = df_identified[[args.type, 'body_site', 'abundance']]
    # convert all body_site with map
    df_full['body_site'] = df_full['body_site'].map(lambda x: body_site_map(x) )
    # Sort on # Reads aligned
    # make new column that is # of reads aligned to sample (% reads in sample) string format
    def quantval(x):
        reads_aligned_formatted = f"{x['# Reads Aligned']:,}"
        if x['abundance'] == x['% Reads']:
            # format # reads aligned to have commas
            return f"{reads_aligned_formatted} ({x['abundance']:.2f}%)"
        else:
            return f"{reads_aligned_formatted} ({x['abundance']:.2f}% - {x['% Reads']:.2f}%)"
    df_full['Quant'] = df_full.apply(lambda x: quantval(x), axis=1)
    # append MicrobeRT Probability to TASS Score if it is present
    # add body site to Sample col with ()
    if args.percentile:
        # filter where value is above percentile only
        df_full = df_full[df_full['HHS Percentile'] >= args.percentile]
    def make_sample(x):
        if not x['body_site']:
            return ""
        else:
            return f"{x['Specimen ID']} ({x['body_site']})"
    df_full['Specimen ID (Type)'] = df_full.apply(lambda x: make_sample(x), axis=1)
    # remove anything in name where "phage" is in it
    df_full = df_full[~df_full[args.type].str.contains("phage")]
    # group on sampletype and get sum of abundance col
    # get the sum of abundance for each sample

    plotbuffer = dict()
    # check if all of the 'Sample Type' is Unknown only, if so then set vairbale
    isUnknownAll = df_full['body_site'].str.contains("Unknown").all()

    mapped_colids = {
        "Detected Organism": "name",
        "Taxonomic ID #": "tax_id",
    }
    if not isUnknownAll and args.distributions and os.path.exists(args.distributions):
        stats_dict, site_counts = import_distributions(
            args.distributions,
            mapped_colids.get(args.type, args.type),
            []
        )
        # # for each of the entries, figure out the top 10 norm_abundance for each site and the corresponding name and tax_id from stats_dict
        # sites_tops = {}
        # body_sites = df_full['body_site'].unique().tolist()
        # for entry in stats_dict:
        #     taxid, body_site = entry
        #     if stats_dict[entry].get('rank') != "species":
        #         continue
        #     if body_site not in sites_tops:
        #         sites_tops[body_site] = []
        #     mean_abundance = stats_dict[entry]['mean']
        #     sites_tops[body_site].append((mean_abundance, taxid))

        # for site, top_list in sites_tops.items():
        #     # sort top_list by mean_abundance desc
        #     top_list_sorted = sorted(top_list, key=lambda x: x[0], reverse=True)
        #     # take top 10
        #     sites_tops[site] = top_list_sorted[:10]
        # # print the top 10 for each site
        # for site, top_list in sites_tops.items():
        #     print(f"Top 10 organisms for site {site}:")
        #     for mean_abundance, taxid in top_list:
        #         print(f"  TaxID: {taxid}, Mean Abundance: {mean_abundance}")
        # exit()

        for _, row in df_full.iterrows():
            # taxid, body_site, stats, args, result_df
            # if taxid and body site not in stats dict then make it empty or 0
            # if (row['tax_id'], row['body_site']) in dists or row['tax_id'] not in taxidsonly or len(body_sites)== 0:
            if (row[args.type], row['body_site']) not in stats_dict:
                stats = {
                    'mean': 0,
                    'std': 0,
                    "min_abundance":0,
                    "max_abundance": 0,
                    "variance": 0,
                    "body_site": "Unknown",
                    "tax_id": "Unknown",
                    "Detected Organism": "Unknown",
                    "rank": "Unknown",
                    "abundances": [],
                    "TASS Score": 0,
                    "gini_coefficient": 0,
                }
            else:
                stats = stats_dict[(row[args.type], row['body_site'])]
            buffer = make_vplot(
                row[args.type],
                stats,
                args.type if args.type == "Detected Organism" else "tax_id",
                df_full,
                percentile_column="HHS Percentile"
            )
            plotbuffer[(row[args.type], row['body_site'])] = buffer
    df_high_cons_low_conf = pd.DataFrame()
    if args.min_conf and args.min_conf > 0:
        df_high_cons_low_conf = df_full[(df_full['High Consequence'] == True) & (df_full['TASS Score'].astype(float) < args.min_conf)]
        df_full = df_full[df_full['TASS Score'].astype(float) >= args.min_conf]
    df_full['TASS Score'] = df_full['TASS Score'].apply(lambda x: f"{x}" if not pd.isna(x) else 0)
    print(f"Size of of full list of organisms: {df_full.shape[0]}")
    print(f"Size of of low confidence high consequence pathogens: {df_high_cons_low_conf.shape[0]}")
    # lambda x add the % Reads column to name column
    df_identified, df_potentials, df_commensal, df_unidentified= split_df(df_full)
    remap_headers = {
        "name": "Detected Organism",
        'taxid': "Taxonomic ID #",
        "# Reads Aligned": "# Reads Aligned to Sample",
        "body_site": "Specimen Type",
        "abundance": "% of Aligned",
        # "Pathogenic Sites": "Locations",
        "% Reads": "% Reads of Organism",
        "Microbial Category": "Microbial Category",
        'Quant': "# Reads Aligned",
        "Gini Coefficient": "Gini Coeff",
    }
    df_identified= df_identified.rename(columns=remap_headers)

    if args.show_unidentified:
        df_unidentified= df_unidentified.rename(columns=remap_headers)
    else:
        df_unidentified = pd.DataFrame()

    if args.show_commensals:
        df_commensal = df_commensal.rename(columns=remap_headers)
    else:
        df_commensal = pd.DataFrame()

    if args.show_potentials:
        df_potentials = df_potentials.rename(columns=remap_headers)
    else:
        df_potentials = pd.DataFrame()
    version = args.version
    missing_samples = args.missing_samples
    create_report(
        args.output,
        df_identified,
        df_potentials,
        df_unidentified,
        df_commensal,
        df_high_cons_low_conf,
        plotbuffer,
        version,
        missing_samples,
        args.min_conf,
        names_map=names_map
    )






if __name__ == "__main__":
    main()
