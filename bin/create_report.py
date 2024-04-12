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
from reportlab.lib import colors
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter
from distributions import import_distributions, make_vplot

from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer, Image
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
    parser.add_argument("-a", "--abundance_col", metavar="ABU", required=False, default='% Total Reads',
                        help="Name of abundance column, default is abundance")
    parser.add_argument("-x", "--id_col", metavar="IDCOL", required=False, default='Name',
                        help="Name of id column, default is id")
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
    # tsv_data = """
    # Name\tSample\tSample Type\t% Aligned\t% Total Reads\t# Aligned\tIsAnnotated\tSites\tType\tTaxid\tStatus\tGini Coefficient\tMeanBaseQ\tMeanMapQ\tBreadth of Coverage\tDepth of Coverage
    # Staphylococcus aureus\tshortreads\tstool\t0.0008\t0.0008\t2\tYes\tstool\tCommensal\t1280\testablished\t0.4
    # Klebsiella pneumoniae\tshortreads\tstool\t0.002\t0.002\t20\tYes\t"abscess, stool, skin, urine"\tCommensal\t573\testablished\t0.23
    # Dickeya fangzhongdai\tshortreads\tstool\t0.0002\t0.0002\t2\tNo\t\t\t1778540\tN/A\t0.95
    # Pediococcus acidilactici\tlongreads\toral\t0.0005\t0.0005\t5\tNo\t\t\t1254\tN/A\t0.9
    # Neisseria gonorrhoeae\tlongreads\toral\t0.025\t0.025\t120\tYes\t"blood, oral, stool, urine"\tPathogen\t485\testablished\t0.02
    # Escherichia coli\tshortreads\tstool\t0.01\t0.01\t100\tNo\t\t\t93061\tN/A\t0.48
    # Metabacillus litoralis\tshortreads\tstool\t0.08\t0.08\t800\tNo\t\t\t152268\tN/A\t0.80
    # Fluviibacter phosphoraccumulans\tlongreads\toral\t0.0005\t0.0005\t5\tNo\t\t\t1751046\tN/A\t0.96
    # Diaphorobacter ruginosibacter\tlongreads\toral\t0.00003\t0.00003\t1\tNo\t\t\t1715720\tN/A\t0.97
    # """.strip()
    # df = pd.read_csv(StringIO(tsv_data), sep='\t')
    # # Simulating additional data
    # np.random.seed(42)
    # # df['Gini Coefficient'] = np.random.uniform(0, 1, df.shape[0])
    # df['MeanBaseQ'] = np.random.uniform(20, 40, df.shape[0])
    # df['MeanMapQ'] = np.random.uniform(30, 60, df.shape[0])
    # df['Breadth of Coverage'] = np.random.uniform(50, 100, df.shape[0])
    # df['Depth of Coverage'] = np.random.uniform(10, 100, df.shape[0])

    df = pd.read_csv(inputfile, sep='\t')


    # sort the dataframe by the Sample THEN the # Reads
    df = df.sort_values(by=[ "Type", "Sample", "# Aligned"], ascending=[False, True, False])
    # trim all of NAme column  of whitespace either side
    df['Name'] = df['Name'].str.strip()

    dictnames = {
        11250: "human respiratory syncytial virus B",
        12814: "human respiratory syncytial virus A",
    }

    # df['Organism'] = df['Name']
    for row_idx, row in df.iterrows():
        if row['Status'] == 'putative':
            #  update index of row to change organism name to bold
            df.at[row_idx, 'Name'] = f'{row.Name}*'
        # change the Name column if the mapnames for taxid is in the dict
        df['Name'] = df[['Name', 'Taxid']].apply(lambda x: dictnames[x['Taxid']] if x['Taxid'] in dictnames else x['Name'], axis=1)
    return df

def split_df(df_full):
    # Filter DataFrame for IsAnnotated == 'Yes' and 'No'
    df_yes = df_full[df_full['IsAnnotated'] == 'Yes'].copy()
    df_no = df_full[df_full['IsAnnotated'] == 'No'].copy()
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

def prepare_data_with_headers(df, plot_dict, include_headers=True, columns=None):
    styles = getSampleStyleSheet()
    data = []
    if not columns:
        columns = df.columns.values[:-1]  # Assuming last column is for plots which should not be included in text headers
    if include_headers:
        headers = [Paragraph('<b>{}</b>'.format(col), styles['Normal']) for col in columns]
        if len(plot_dict.keys()) > 0:
            headers.append(Paragraph('<b>Percentile of Healthy Subject</b>', styles['Normal']))  # Plot column header
        data.append(headers)
    for index, row in df.iterrows():
        row_data = [Paragraph(format_cell_content(str(cell)), styles['Normal']) for cell in row[columns][:]]  # Exclude plot data
        # Insert the plot image
        if len(plot_dict.keys()) > 0:
            plot_key = (row['Organism'], row['Sample Type'])
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
        colorindexcol = 3

        sampleindx = df.columns.get_loc('Sample Type')
        # Example post-processing to mark cells
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


def create_report(
    output,
    df_identified,
    df_unidentified,
    plotbuffer
):

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
    # sort df_identified by Alignment Conf
    df_identified = df_identified.sort_values(by=['Alignment Conf'], ascending=False)
    # filter out so only Class is PAthogen
    df_identified['Class']
    df_identified_paths = df_identified[df_identified['Class'] == 'Pathogen']
    df_identified_others = df_identified[df_identified['Class'] != 'Pathogen']

    elements = []
    ##########################################################################################
    ##### Section to make the Top Table - all annotated commensal or otherwise
    if not df_identified_paths.empty:
        columns_yes = df_identified_paths.columns.values
        # print only rows in df_identified with Gini Coeff above 0.2
        columns_yes = ["Sample", "Sample Type", "Organism", "Class", "% Reads in Sample", "# Aligned to Sample", "Alignment Conf", "Locations"]
        # Now, call prepare_data_with_headers for both tables without manually preparing headers
        data_yes = prepare_data_with_headers(df_identified_paths, plotbuffer, include_headers=True, columns=columns_yes)
        table_style = return_table_style(df_identified_paths, color_pathogen=True)
        table = make_table(
            data_yes,
            table_style=table_style
        )
        # Add the title and subtitle
        title = Paragraph("Organism Discovery Analysis", title_style)
        subtitle = Paragraph(f"This report was generated using TaxTriage {version} on {date} and is derived from an in development spreadsheet of human-host pathogens. It will likely change performance as a result of rapid development practices.", subtitle_style)
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
    subtext_para = Paragraph("Organisms marked with * are putative and have relatively lower references listing their annotations as a pathogen in the given sample types", subtext_style)
    elements.append(subtext_para)
    elements.append(Spacer(1, 12))
    subtext_para = Paragraph("Light yellow cells represent pathogens annotated in sample type(s) other than your listed one. Green represents a match with your sample type", subtext_style)
    elements.append(subtext_para)

    if not df_identified_others.empty:
        columns_yes = df_identified_others.columns.values
        # print only rows in df_identified with Gini Coeff above 0.2
        columns_yes = ["Sample", "Sample Type", "Organism", "Class", "% Reads in Sample", "# Aligned to Sample", "Alignment Conf", "Locations"]
        # Now, call prepare_data_with_headers for both tables without manually preparing headers
        data_yes_others = prepare_data_with_headers(df_identified_others, plotbuffer, include_headers=True, columns=columns_yes)
        table_style = return_table_style(df_identified_others, color_pathogen=True)
        table = make_table(
            data_yes_others,
            table_style=table_style
        )
        # Add the title and subtitle
        title = Paragraph("Other Organisms Annotated", title_style)
        subtitle = Paragraph(f"All Organisms were identified but were not listed as a pathogen", subtitle_style)
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
        second_subtitle = "The following table displays the unannotated organisms and their alignment statistics. Be aware that this is the exhaustive list of all organisms contained within the samples that had atleast one read aligned"
        elements.append(Paragraph(second_title, title_style))
        elements.append(Paragraph(second_subtitle, subtitle_style))

        columns_no = ['Sample',  "Sample Type", 'Organism', '% Reads in Sample', '# Aligned to Sample', "Alignment Conf" ]
        data_no = prepare_data_with_headers(df_unidentified, plotbuffer, include_headers=True, columns=columns_no)
        table_style = return_table_style(df_unidentified, color_pathogen=False)
        table_no = make_table(
            data_no,
            table_style=table_style
        )
        elements.append(table_no)
    elements.append(Spacer(1, 12))  # Space between tables
    ##########################################################################################
    #####  Build the PDF
    print(len(elements))
    # Adjust the build method to include the draw_vertical_line function
    doc.build(elements, onFirstPage=draw_vertical_line, onLaterPages=draw_vertical_line)

    print(f"PDF generated: {pdf_file}")



def main():
    args = parse_args()
    df_full = import_data(args.input)

    # change column "id" in avbundance_data to "tax_id" if args.type is "name"
    df_full = df_full.rename(columns={args.id_col: args.type})
    df_full = df_full.rename(columns={args.sitecol: 'body_site'})
    df_full = df_full.rename(columns={args.abundance_col: 'abundance'})
    df_full = df_full.dropna(subset=[args.type])
    # df_identified = df_identified[[args.type, 'body_site', 'abundance']]
    # Map several names for common groups for body_site
    body_site_map = {
        "gut": "stool",
        "nose": "nasal",
        "vagina": "vaginal",
        "teeth": "oral",
        "resp": "nasal",
        "abscess": "skin",
        "absscess": "skin",
        "sputum": "oral",
        "mouth": "oral",
        "urinary tract": "urine",
        "ear": "skin",
        "lung": "nasal",
        "eye": "eye",
        "sinus": "nasal",
        "urogenital": "skin",
        "cornea": "eye",
    }
    # convert all body_site with map
    df_full['body_site'] = df_full['body_site'].map(lambda x: body_site_map[x] if x in body_site_map else x)
    plotbuffer = dict()
    if args.distributions and os.path.exists(args.distributions):
        stats_dict = import_distributions(
            args.distributions,
            args.type,
            []
        )

        for index, row in df_full.iterrows():
            # taxid, body_site, stats, args, result_df
            # if taxid and body site not in stats dict then make it empty or 0
            if (row[args.type], row['body_site']) not in stats_dict:
                stats = {
                    'mean': 0,
                    'std': 0,
                    "min_abundance":0,
                    "max_abundance": 0,
                    "variance": 0,
                    "body_site": "Unknown",
                    "tax_id": "Unknown",
                    "name": "Unknown",
                    "rank": "Unknown",
                    "abundances": [],
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
    df_full['Sites'] = df_full['Sites'].fillna("Unknown")

    df_full['Alignment Conf'] = df_full['Gini Coefficient'].apply(lambda x: f"{1-x:.2f}" if not pd.isna(x) else 0)
    print(f"Size of of full list of organisms: {df_full.shape[0]}")
    df_identified, df_unidentified = split_df(df_full)
    remap_headers = {
        "Name": "Organism",
        "name": "Organism",
        "# Aligned": "# Aligned to Sample",
        "body_site": "Sample Type",
        "abundance": "% Reads in Sample",
        "Sites": "Locations",
        "% Total Reads": "% Reads in Sample",
        "Type": "Class",
        "Gini Coefficient": "Gini Coeff",
    }
    df_identified= df_identified.rename(columns=remap_headers)
    df_unidentified= df_unidentified.rename(columns=remap_headers)
    create_report(
        args.output,
        df_identified,
        df_unidentified,
        plotbuffer
    )






if __name__ == "__main__":
    main()

