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
import matplotlib.pyplot as plt
from reportlab.lib import colors
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter, landscape
from distributions import import_distributions, make_vplot, body_site_map
from reportlab.graphics.shapes import Line

from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer, Image, ListFlowable, ListItem
from reportlab.platypus.flowables import HRFlowable, Flowable

from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch

from datetime import datetime
from io import StringIO
from reportlab.lib.colors import Color

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
        help="Base pathogen discovery table file, TSV format",
    )
    parser.add_argument(
        "-distributions",
        "--distributions",
        metavar="DISTRIBUTIONS",
        required=False,
        help="TSV file that contains all the distribution information for body sites and organisms",
    )
    parser.add_argument("-a", "--abundance_col", metavar="ABU", required=False, default='% Aligned Reads',
                        help="Name of abundance column, default is abundance")
    parser.add_argument("-c", "--min_conf", metavar="MINCONF", required=False, default=0.5, type=float,
                        help="Value that must be met for a table to report an organism due to confidence column.")
    parser.add_argument("-x", "--id_col", metavar="IDCOL", required=False, default="Detected Organism",
                        help="Name of id column, default is id")
    parser.add_argument("-v", "--version", metavar="VERSION", required=False, default='Local Build',
                        help="What version of TaxTriage is in use")
    parser.add_argument("-s", "--sitecol", metavar="SCOL", required=False, default='Sample Type',
                        help="Name of site column, default is body_site")
    parser.add_argument("-t", "--type", metavar="TYPE", required=False, default='name',
                        help="What type of data is being processed. Options: 'tax_id' or 'name'.",
                        choices=['tax_id', 'name'])
    parser.add_argument(
        "-o",
        "--output",
        metavar="OUTPUT",
        required=True,
        type=str,
        help="Path of output file",
    )

    return parser.parse_args(argv)




# Function to adjust font size based on text length
def adjust_font_size(text, max_length=8, default_font_size=8, min_font_size=6):
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



def import_data(inputfile ):
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

    df = pd.read_csv(inputfile, sep='\t')

    # set % Reads aligned as float
    df['% Reads'] = df['% Reads'].apply(lambda x: float(x) if not pd.isna(x) else 0)
    # set # Reads Aligned as int
    df['# Reads Aligned'] = df['# Reads Aligned'].apply(lambda x: int(x) if not pd.isna(x) else 0)

    # sort the dataframe by the Sample THEN the # Reads
    df = df.sort_values(by=["Specimen ID", "Microbial Category",  "# Reads Aligned"], ascending=[True, False, True])
    # trim all of NAme column  of whitespace either side
    df["Detected Organism"] = df["Detected Organism"].str.strip()
    dictnames = {
        11250: "human respiratory syncytial virus B",
        12814: "human respiratory syncytial virus A",
    }
    # df['Organism'] = df["Detected Organism"]
    for row_idx, row in df.iterrows():
        if row['Status'] == 'putative':
            #  update index of row to change organism name to bold
            df.at[row_idx, "Detected Organism"] = f'{row["Detected Organism"]}*'
        # change the Name column if the mapnames for taxid is in the dict
        df["Detected Organism"] = df[["Detected Organism", 'Taxonomic ID #']].apply(lambda x: dictnames[x['Taxonomic ID #']] if x['Taxonomic ID #'] in dictnames else x["Detected Organism"], axis=1)
    # replace all NaN with ""
    df = df.fillna("")
    return df

def split_df(df_full):
    # Filter DataFrame for IsAnnotated == 'Yes' and 'No'
    # df_ = df_full[df_full['IsAnnotated'] == 'Yes'].copy()
    # append (taxid) from taxid column to Detected Organism

    df_yes = df_full[~df_full['Microbial Category'].isin([ 'Unknown', 'N/A', np.nan, "Commensal", "Potential" ] ) ].copy()
    df_opp = df_full[df_full['Microbial Category'].isin([  "Potential"])].copy()
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

def prepare_data_with_headers(df, plot_dict, include_headers=True, columns=None):
    data = []
    # convert k2 reads to int
    df['K2 Reads'] = df['K2 Reads'].apply(lambda x: int(x) if not pd.isna(x) else 0)
    if not columns:
        columns = df.columns.values[:-1]  # Assuming last column is for plots which should not be included in text headers
    if include_headers:
        headers = [Paragraph('<b>{}</b>'.format(col), styles['Normal']) for col in columns]
        if len(plot_dict.keys()) > 0:
            headers.append(Paragraph('<b>Percentile of Healthy Subject (HHS)</b>', styles['Normal']))  # Plot column header
        data.append(headers)
    for index, row in df.iterrows():
        row_data = [Paragraph(format_cell_content(str(cell)), small_font_style ) for cell in row[columns][:]]  # Exclude plot data
        # Insert the plot image
        if len(plot_dict.keys()) > 0:
            plot_key = (row['Detected Organism'], row['Specimen Type'])
            if plot_key in plot_dict:
                plot_image = Image(plot_dict[plot_key])
                plot_image.drawHeight = 0.5 * inch  # Height of the image
                plot_image.drawWidth = 1* inch  # Width of the image, adjusted from your figsize
                row_data.append(plot_image)
        data.append(row_data)
    return data

def return_table_style(df, color_pathogen=False):
    # Start with the basic style
    table_style = TableStyle([
        ('BACKGROUND', (0,0), (-1,0), colors.gray), # header
        ('TEXTCOLOR', (0,0), (-1,0), colors.whitesmoke),
        ('GRID', (0,0), (-1,-1), 1, colors.black),
        ('VALIGN', (0,0), (-1,-1), 'TOP'),

    ])
    if color_pathogen:
        # Placeholder for cells to color (row_index, col_index) format
        cells_to_color = []
        colorindexcol = 2

        sampleindx = df.columns.get_loc('Microbial Category')
        # Example post-processing to mark cells
        for row_idx, row in df.iterrows():
            val = row['Microbial Category']
            sites =  row['Locations']
            # if nan then set to empty string
            if pd.isna(sites):
                sites = ""
            # Get Sample Type value from row
            sampletype = row['Specimen Type']
            if val != "Commensal" and sampletype in sites:
                color = 'lightgreen'
            elif val != "Commensal" and row.AnnClass == 'Derived':
                color = 'lightblue'
            elif val != "Commensal":
                color = "papayawhip"
            elif val == "Commensal" and row.AnnClass == "Derived":
                color = 'lightblue'
            else:
                color = "white"
            # Ensure indices are within the table's dimensions
            style_command = ('BACKGROUND', (colorindexcol, row_idx+1), (colorindexcol, row_idx+1), color)  # Or lightorange based on condition
            table_style.add(*style_command)
    else:
        table_style.add(*('BACKGROUND', (0,1), (-1,-1), colors.white))
    return table_style

def make_table(data, table_style=None):
    # Set table style to the return value of the function

    # Set style
    # Applying the custom style to the title and subtitle
    # Table configuration
    table = Table(data, repeatRows=1)
    # Apply this style to your tables
    table.setStyle(table_style)
    return table
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


def create_report(
    output,
    df_identified,
    df_opportunistic,
    df_unidentified,
    df_commensals,
    plotbuffer,
    version=None
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
    # sort df_identified by Confidence Metric (0-1)
    # filter out so only Class is PAthogen
    # df_identified = df_identified.sort_values(by=['Specimen ID', '# Reads Aligned'], ascending=[False, True])
    # df_opportunistic = df_opportunistic.sort_values(by=['Specimen ID', '# Reads Aligned'], ascending=[False, False])
    df_identified_paths = df_identified
    df_identified_others = df_commensals
    # df_identified_others = df_identified[df_identified['Class'] != 'Pathogen']
    # df_unidentified = df_unidentified.sort_values(by=['Specimen ID', '# Reads Aligned'], ascending=[False, False])
    elements = []
    ##########################################################################################
    ##########################################################################################
    ##### Section to make the Top Table - all annotated commensal or otherwise
    if not df_identified_paths.empty:
        columns_yes = df_identified_paths.columns.values
        # print only rows in df_identified with Gini Coeff above 0.2
        columns_yes = [
            "Specimen ID (Type)",
            "Detected Organism",
            "Microbial Category",
            "# Reads Aligned",
            "Confidence Metric (0-1)",
            "Taxonomic ID #", "Pathogenic Subsp/Strains",
            "K2 Reads"
            ]
        # check if all K2 reads column are 0 or nan
        if df_identified_paths['K2 Reads'].sum() == 0:
            columns_yes = columns_yes[:-1]
        # Now, call prepare_data_with_headers for both tables without manually preparing headers
        data_yes = prepare_data_with_headers(df_identified_paths, plotbuffer, include_headers=True, columns=columns_yes)
        table_style = return_table_style(df_identified_paths, color_pathogen=True)
        table = make_table(
            data_yes,
            table_style=table_style
        )
        # Add the title and subtitle
        title = Paragraph("Organism Discovery Report", title_style)
        subtitle = Paragraph(f"This report was generated using TaxTriage <b>{version}</b> on <b>{date}</b> and is derived from an in development spreadsheet of human-host pathogens.", subtitle_style)
        elements.append(title)
        elements.append(subtitle)
        elements.append(Spacer(1, 12))
        elements.append(table)
        elements.append(Spacer(1, 12))  # Space between tables
    # Adding regular text

    styles = getSampleStyleSheet()
    # Adding subtext (you can adjust the style to make it look like subtext)
    subtext_style = styles["BodyText"]
    subtext_style.fontSize = 10  # Smaller font size for subtext
    subtext_style.leading = 12
    subtext_para = Paragraph("Organisms marked with * are putative and have relatively lower references listing their annotations as a pathogen in the given sample types. Classifications of pathogens are described as:", subtext_style)
    elements.append(subtext_para)
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
        "≡: Indicates a pathogenic subspecies/serotype/strain/etc with the same name as the species listed, just different taxids."
    ]

    bullet_list = ListFlowable(
        [ListItem(Paragraph(item, subtext_style)) for item in bullet_list_items],
        bulletType='bullet',
        start='circle'
    )
    elements.append(bullet_list)

    elements.append(horizontal_line)

    elements.append(Spacer(1, 12))

    # Define the column explanations
    column_explanations = [
        "Specimen ID (Type): The unique identifier for the sample, including the type of specimen (e.g., blood, tissue).",
        "Detected Organism: The organism detected in the sample, which could be a bacterium, virus, fungus, or parasite.",
        "Microbial Category: The classification of the organism, indicating whether it is primary, opportunistic, commensal, or potential.",
        "# Reads Aligned: The number of reads from the sequencing data that align to the organism's genome, indicating its presence. (%) refers to all alignments (more than 1 alignment per read can take place) for that species across the entire sample.",
        "Confidence Metric (0-1): A metric between 0 and 1 that reflects the confidence of the organism's detection, with 1 being the highest confidence.",
        "Taxonomic ID #: The taxid for the organism according to NCBI Taxonomy, which provides a unique identifier for each species.",
        "Pathogenic Subsp/Strains: Indicates specific pathogenic subspecies, serotypes, or strains, if detected in the sample. (%) indicates the percent of all aligned reads belonging to that strain.",
        "K2 Reads: The number of reads classified by Kraken2, a tool for taxonomic classification of sequencing data."
    ]

    # Create bullet points for each column explanation
    bullet_list_items = [
        Paragraph(item, subtext_style) for item in column_explanations
    ]

    # Create the bullet list
    bullet_list = ListFlowable(
        [ListItem(item) for item in bullet_list_items],
        bulletType='bullet',
        start='circle'
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

    bullet_list_items = [
        "Green/White: Direct match for the taxid/organism name with your sample type from the database.",
        "Blue: Derived Pathogenicity from any listed pathogenic strains of a given organism.",
        "Beige: Pathogens annotated in sample type(s) other than your listed one.",
    ]

    # Create a list of bullet items with specified colors
    bullet_colors = [colors.lightgreen, colors.lightblue, colors.papayawhip,  ]
    style = styles['Normal']

    # Create custom ListItems with colored bullets
    custom_list_items = [
        ListItem(Paragraph(item, style), bulletBorder=colors.black,  bulletColor=bullet_colors[idx], )
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

    subtext_para = Paragraph("Read amounts are represented as the <b>total number of aligned reads</b> of sufficient mapping quality <b>(% aligned for all reads in sample)</b>", subtext_style)
    elements.append(subtext_para)

    elements.append(horizontal_line)

    ##########################################################################################
    #### Table on opportunistic pathogens
    if not df_opportunistic.empty:
        columns_opp = ["Specimen ID (Type)", "Detected Organism",
                       "Microbial Category", "# Reads Aligned",
                       "Confidence Metric (0-1)", "Taxonomic ID #",
                       "Pathogenic Subsp/Strains", "K2 Reads"
                       ]
        if df_opportunistic['K2 Reads'].sum() == 0:
            columns_opp = columns_opp[:-1]
        data_opp = prepare_data_with_headers(df_opportunistic, plotbuffer, include_headers=True, columns=columns_opp)
        table_style = return_table_style(df_opportunistic, color_pathogen=True)
        table = make_table(
            data_opp,
            table_style=table_style
        )
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
        columns_yes = ["Specimen ID (Type)",
                       "Detected Organism", "Microbial Category",
                       "# Reads Aligned", "Confidence Metric (0-1)", "Taxonomic ID #",
                       "K2 Reads"]
        # check if all K2 reads column are 0 or nan
        # if df_identified_paths['K2 Reads'].sum() == 0:
        #     columns_yes = columns_yes[:-1]
        # Now, call prepare_data_with_headers for both tables without manually preparing headers
        data_yes = prepare_data_with_headers(df_identified_others, plotbuffer, include_headers=True, columns=columns_yes)
        table_style = return_table_style(df_identified_others, color_pathogen=False)
        table = make_table(
            data_yes,
            table_style=table_style
        )
        title = Paragraph("Commensals", title_style)
        subtitle = Paragraph(f"These were identified & were listed as a commensal directly", subtitle_style)
        elements.append(title)
        elements.append(subtitle)
        elements.append(Spacer(1, 12))  # Space between tables

        elements.append(table)
        elements.append(Spacer(1, 12))  # Space between tables
    # # Adding regular text

    elements.append(Spacer(1, 12))
    if not df_unidentified.empty:
        ##########################################################################################
        ### Section to Make the "Unannotated" Table
        second_title = "Unannotated Organisms"
        second_subtitle = "The following table displays the unannotated organisms and their alignment statistics. Be aware that this is the exhaustive list of all organisms (species only) contained within the samples that had atleast one read aligned"
        elements.append(Paragraph(second_title, title_style))
        elements.append(Paragraph(second_subtitle, subtitle_style))

        columns_no = ['Specimen ID (Type)', 'Detected Organism','# Reads Aligned', "Confidence Metric (0-1)", "K2 Reads" ]
        data_no = prepare_data_with_headers(df_unidentified, plotbuffer, include_headers=True, columns=columns_no)
        if df_unidentified['K2 Reads'].sum() == 0:
            columns_no = columns_no[:-1]
        table_style = return_table_style(df_unidentified, color_pathogen=False)
        table_no = make_table(
            data_no,
            table_style=table_style
        )
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
    df_full = import_data(args.input)


    # change column "id" in avbundance_data to "tax_id" if args.type is "Detected Organism"
    df_full = df_full.rename(columns={args.id_col: args.type})
    df_full = df_full.rename(columns={args.sitecol: 'body_site'})
    df_full = df_full.rename(columns={args.abundance_col: 'abundance'})
    df_full = df_full.dropna(subset=[args.type])
    # df_identified = df_identified[[args.type, 'body_site', 'abundance']]
    # convert all body_site with map
    df_full['body_site'] = df_full['body_site'].map(lambda x: body_site_map(x) )
    # Sort on # Reads aligned
    df_full = df_full.sort_values(by=["# Reads Aligned"], ascending=False)
    # make new column that is # of reads aligned to sample (% reads in sample) string format
    df_full['Quant'] = df_full.apply(lambda x: f"{x['# Reads Aligned']} ({x['abundance']:.2f}%)", axis=1)
    # add body sit to Sample col with ()
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
    if args.distributions and os.path.exists(args.distributions):
        stats_dict, site_counts = import_distributions(
            args.distributions,
            args.type,
            []
        )

        for index, row in df_full.iterrows():
            # taxid, body_site, stats, args, result_df
            # if taxid and body site not in stats dict then make it empty or 0
            taxidsonly = [key[0] for key in stats_dict.keys()]
            bodysites = [key[1] for key in stats_dict.keys()]
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
            rank = stats['rank']
            buffer = make_vplot(
                row[args.type],
                stats,
                args.type,
                df_full
            )
            plotbuffer[(row[args.type], row['body_site'])] = buffer
    # convert all locations nan to "Unknown"
    df_full['Pathogenic Sites'] = df_full['Pathogenic Sites'].fillna("Unknown")
    if args.min_conf and args.min_conf > 0:
        df_full = df_full[df_full['TASS Score'].astype(float) >= args.min_conf]
    # df_full['name'] = df_full['name'] + " (" + df_full['Taxonomic ID #'].astype(str) + ")"
    df_full['Confidence Metric (0-1)'] = df_full['TASS Score'].apply(lambda x: f"{x:.2f}" if not pd.isna(x) else 0)
    print(f"Size of of full list of organisms: {df_full.shape[0]}")
    df_identified, df_opportunistic, df_commensal, df_unidentified= split_df(df_full)
    remap_headers = {
        "name": "Detected Organism",
        'taxid': "Taxonomic ID #",
        "# Reads Aligned": "# Reads Aligned to Sample",
        "body_site": "Specimen Type",
        "abundance": "% of Aligned",
        "Pathogenic Sites": "Locations",
        "% Reads": "% Reads of Organism",
        "Microbial Category": "Microbial Category",
        'Quant': "# Reads Aligned",
        "Gini Coefficient": "Gini Coeff",
    }
    df_identified= df_identified.rename(columns=remap_headers)
    df_unidentified= df_unidentified.rename(columns=remap_headers)
    df_commensal = df_commensal.rename(columns=remap_headers)
    df_opportunistic = df_opportunistic.rename(columns=remap_headers)
    version = args.version
    create_report(
        args.output,
        df_identified,
        df_opportunistic,
        df_unidentified,
        df_commensal,
        plotbuffer,
        version,
    )






if __name__ == "__main__":
    main()

