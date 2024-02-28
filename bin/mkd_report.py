import csv
import subprocess
import sys
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
import pypandoc

def csv_to_markdown(csv_file_path):
    """Convert CSV file to Markdown table."""
    with open(csv_file_path, 'r', newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        headers = next(reader)
        markdown_table = '| ' + ' | '.join(headers) + ' |\n'
        markdown_table += '|-' * len(headers) + '|\n'
        for row in reader:
            markdown_table += '| ' + ' | '.join(row) + ' |\n'
        return markdown_table

def convert_to_pdf(markdown_data, output_pdf_path):
    """Convert markdown data to PDF using pypandoc."""
    output = pypandoc.convert_text(markdown_data, 'pdf', format='md', outputfile=output_pdf_path, extra_args=['-V', 'geometry:margin=1in'])
    return output  # Not needed unless you want to capture the output

def main(csv_file_path, output_pdf_path):
    markdown_data = csv_to_markdown(csv_file_path)
    convert_to_pdf(markdown_data, output_pdf_path)
    print(f"PDF generated: {output_pdf_path}")


if __name__ == "__main__":
    args = parse_args()

    csv_file_path = args.input
    output_pdf_path = args.output
    main(csv_file_path, output_pdf_path)

