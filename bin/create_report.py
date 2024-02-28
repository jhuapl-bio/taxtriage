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


import pandas as pd
import numpy as np
from reportlab.lib import colors
from reportlab.pdfgen import canvas

from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer
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
        "-o",
        "--output",
        metavar="OUTPUT",
        required=True,
        type=str,
        help="Path of output file",
    )

    return parser.parse_args(argv)




# Function to adjust font size based on text length
def adjust_font_size(text, max_length=40, default_font_size=10, min_font_size=6):
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
    tsv_data = """
    Name\t% Aligned\t% Total Reads\t# Aligned\tIsAnnotated\tSites\tType\tTaxid\tStatus
    Staphylococcus aureus\t0.6\t0.6\t32\tYes\tblood\tCommensal\t1280\testablished
    Klebsiella pneumoniae\t46.22\t46.22\t2447\tYes\t"abscess, blood , skin, urine"\tCommensal\t573\testablished
    Dickeya fangzhongdai\t2.87\t2.87\t152\tNo\t\t\t1778540\tN/A
    Pediococcus acidilactici\t15.94\t15.94\t844\tNo\t\t\t1254\tN/A
    Neisseria gonorrhoeae\t15.17\t15.17\t803\tYes\t"blood, urine"\tPathogen\t485\testablished
    Staphylococcus aureus subsp. aureus NCTC 8325\t1.23\t1.23\t65\tNo\t\t\t93061\tN/A
    Metabacillus litoralis\t0.08\t0.08\t4\tNo\t\t\t152268\tN/A
    Fluviibacter phosphoraccumulans\t12.18\t12.18\t645\tNo\t\t\t1751046\tN/A
    Diaphorobacter ruginosibacter\t5.7\t5.7\t302\tNo\t\t\t1715720\tN/A
    """.strip()
    df = pd.read_csv(inputfile, sep='\t')
    # df = pd.read_csv(StringIO(tsv_data), sep='\t')
    remap_headers = {
        "Name": "Organism",
        "# Aligned": "# Aligned to Sample",
        "Sites": "Locations",
        "% Total Reads": "% Reads in Sample",
        "Type": "Class",
    }
    df = df.rename(columns=remap_headers)
    # sort the dataframe by the Sample THEN the # Reads
    df = df.sort_values(by=[ "Class", "Sample", "# Aligned to Sample"], ascending=[False, True, False])
    for row_idx, row in df.iterrows():
        if row['Status'] == 'putative':
            #  update index of row to change organism name to bold
            df.at[row_idx, 'Organism'] = f'{row.Organism}*'

    # Filter DataFrame for IsAnnotated == 'Yes' and 'No'
    df_yes = df[df['IsAnnotated'] == 'Yes'].copy()
    df_no = df[df['IsAnnotated'] == 'No'].copy()
    # reset index
    df_yes.reset_index(drop=True, inplace=True)
    df_no.reset_index(drop=True, inplace=True)
    return df_yes, df_no




# Custom styles for title and subtitle with right alignment
left_align_style = ParagraphStyle(
    name='leftAlign',
    parent=getSampleStyleSheet()['Normal'],
    alignment=0,  # 2 is for right alignment
    fontSize=12,
    spaceAfter=10,
)
subtitle_style = ParagraphStyle(
    name='leftAlignSubtitle',
    parent=getSampleStyleSheet()['Normal'],
    alignment=0,  # Right align
    fontSize=12,
    spaceAfter=10,
)

title_style = ParagraphStyle(
    name='leftAlignTitle',
    parent=getSampleStyleSheet()['Title'],
    alignment=0,  # Right align
    fontSize=14,
    spaceAfter=10,
)

def prepare_data_with_headers(df, include_headers=True, color_status=False, columns=None):
    styles = getSampleStyleSheet()
    header_style = styles['Normal']  # Customize this style for headers
    data = []
    if not columns:
        columns = df.columns.values
    # Prepare header row
    if include_headers:
        headers = [Paragraph('<b>{}</b>'.format(col), header_style) for col in columns]
        data.append(headers)

    # Prepare data rows with conditional coloring for the "Sample" cell
    for index, row in df.iterrows():
        row_data = []
        for col in columns:
            cell = row[col]
            text = format_cell_content(cell)
            paragraph = Paragraph(text, styles['Normal'])
            row_data.append(paragraph)
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
        colorindexcol = 3

        sampleindx = df.columns.get_loc('Sample Type')
        # Example post-processing to mark cells
        print(df.columns.values)
        for row_idx, row in enumerate(df.itertuples(index=False)):
            val = row.Class
            status = row.Status
            sites = row.Locations
            # if nan then set to empty string
            if pd.isna(sites):
                sites = ""
            # Get Sample Type value from row
            sampletype = row[sampleindx]
            if val == 'Pathogen' and sampletype in sites:
                color = 'lightgreen'
            elif val == 'Pathogen' :
                color = 'lightyellow'
            else:
                color = 'white'
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


def create_report(output, df_yes, df_no):

    # PDF file setup
    pdf_file = output
    doc = SimpleDocTemplate(pdf_file, pagesize=letter)
    # Placeholder values for version and date
    #### Section to style things up a bit
    # Set the left margin to 10% of the width of a letter size page (8.5 inches)
    left_margin = 0.1 * letter[0]
    right_margin = left_margin / 5
    top_margin = bottom_margin = 0.1 * letter[1]


    # Modify the doc setup to include custom margins
    doc = SimpleDocTemplate(
        pdf_file,
        pagesize=letter,
        leftMargin=left_margin,
        rightMargin=right_margin,
        topMargin=top_margin,
        bottomMargin=bottom_margin
    )
    version = "1.3.2"  # Example version
    # get datetime of year-mont-day hour:min
    date = datetime.now().strftime("%Y-%m-%d %H:%M")  # Current date

    elements = []
    ##########################################################################################
    ##### Section to make the Top Table - all annotated commensal or otherwise
    if not df_yes.empty:
        columns_yes = df_yes.columns.values
        columns_yes = ["Sample", "Sample Type", "Organism", "Class", "% Reads in Sample", "# Aligned to Sample", "Locations"]
        # Now, call prepare_data_with_headers for both tables without manually preparing headers
        data_yes = prepare_data_with_headers(df_yes, columns=columns_yes, include_headers=True, color_status=True)
        table_style = return_table_style(df_yes, color_pathogen=True)
        table = make_table(
            data_yes,
            table_style=table_style
        )
        # Add the title and subtitle
        title = Paragraph("Pathogen Discovery Analysis", title_style)
        subtitle = Paragraph(f"This report was generated using TaxTriage {version} on {date} and is derived from an in development spreadsheet of human-host pathogens. It will likely change performance as a result of rapid development practices.", subtitle_style)
        elements = [title, subtitle, Spacer(1, 12)]
        elements.append(table)
        elements.append(Spacer(1, 12))  # Space between tables
    # Adding regular text

    styles = getSampleStyleSheet()

    # Adding subtext (you can adjust the style to make it look like subtext)
    subtext_style = styles["BodyText"]
    subtext_style.fontSize = 10  # Smaller font size for subtext
    subtext_style.leading = 12
    subtext_para = Paragraph("Organisms marked with * are putative and have relatively lower references listing their annotations as a pathogen in the given sample types", subtext_style)
    elements.append(subtext_para)
    elements.append(Spacer(1, 12))
    subtext_para = Paragraph("Light yellow cells represent pathogens annotated in sample type(s) other than your listed one. Green represents a match with your sample type", subtext_style)
    elements.append(subtext_para)
    elements.append(Spacer(1, 12))
    if not df_no.empty:
        ##########################################################################################
        ### Section to Make the "Unannotated" Table
        second_title = "Unannotated Organisms"
        second_subtitle = "The following table displays the unannotated organisms and their alignment statistics. Be aware that this is the exhaustive list of all organisms contained within the samples that had atleast one read aligned"
        elements.append(Paragraph(second_title, title_style))
        elements.append(Paragraph(second_subtitle, subtitle_style))
        columns_no = ['Sample',  "Sample Type", 'Organism', '% Reads in Sample', '# Aligned to Sample' ]
        data_no = prepare_data_with_headers(df_no, include_headers=True, color_status=False, columns=columns_no)
        table_style = return_table_style(df_no, color_pathogen=False)
        table_no = make_table(
            data_no,
            table_style=table_style
        )
        elements.append(table_no)
    elements.append(Spacer(1, 12))  # Space between tables
    ##########################################################################################
    #####  Build the PDF
    # Adjust the build method to include the draw_vertical_line function
    doc.build(elements, onFirstPage=draw_vertical_line, onLaterPages=draw_vertical_line)

    print(f"PDF generated: {pdf_file}")


def main():
    args = parse_args()
    df_yes, df_no = import_data(args.input)
    create_report(args.output, df_yes, df_no)
if __name__ == "__main__":
    main()

