import csv

from collections import defaultdict

from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfbase import pdfmetrics

metadata = "/Users/merribb1/Documents/Projects/APHL/taxtriage/test/metadata.csv"
data_filename = "/Users/merribb1/Documents/Projects/APHL/taxtriage/test/data.csv"


def extract_metadata_from_csv(filename):
    """Extract metadata from a CSV file."""
    metadata = {}
    with open(filename, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            key, value = row
            metadata[key] = value
    return metadata


def render_metadata(c, metadata, start_x, start_y, width_limit):
    """Render metadata dynamically based on available space."""
    font_name = "Helvetica"
    font_size = 12
    bold_font_name = "Helvetica-Bold"

    x = start_x
    y = start_y

    for key, value in metadata.items():
        # Render the key in bold
        key_text = f"{key}: "
        c.setFont(bold_font_name, font_size)  # Set to bold for the key
        key_width = pdfmetrics.stringWidth(key_text, bold_font_name, font_size)
        c.drawString(x, y, key_text)

        # Render the value in regular font
        value_text = f"{value}   "
        c.setFont(font_name, font_size)  # Set to regular for the value
        value_width = pdfmetrics.stringWidth(value_text, font_name, font_size)
        c.drawString(x + key_width, y, value_text)

        # Update position for next key-value pair
        x += key_width + value_width
        if x > width_limit:
            x = start_x
            y -= 0.5 * inch

    return y  # Return the final y position


def generate_pdf_report(data_filename, metadata_filename, output_filename="report.pdf"):
    # Extract metadata
    metadata = extract_metadata_from_csv(metadata_filename)

    # Create a new canvas
    c = canvas.Canvas(output_filename, pagesize=letter)
    width, height = letter

    # Render the Study Name centered at the top as a title
    study_name = metadata.get("Study Name", "")
    c.setFont("Helvetica-Bold", 14)  # Larger font for the title
    title_width = pdfmetrics.stringWidth(study_name, "Helvetica-Bold", 14)
    c.drawString((width - title_width) / 2, height - inch, study_name)
    # Remove study name from metadata to prevent it from being displayed again
    del metadata["Study Name"]

    # Render the remaining metadata
    y_pos = render_metadata(c, metadata, inch, height - 2*inch, width-inch)

    # Draw a line to separate metadata from data
    y_pos -= 0.5 * inch
    c.line(inch, y_pos, width-inch, y_pos)
    y_pos -= 0.5*inch

    # Read data from the data CSV and group by sample name
    data = defaultdict(list)
    with open(data_filename, 'r') as f:
        reader = csv.DictReader(f)  # Using DictReader to handle headers
        for row in reader:
            data[row["Sample"]].append(row)

    # Render grouped data onto the PDF
    for sample, rows in data.items():
        c.setFont("Helvetica", 12)  # Set font for sample name
        c.drawString(inch, y_pos, f"Sample: {sample}")
        y_pos -= 0.5*inch
        for row in rows:
            target, present, abu = row["Target"], int(
                row["Present"]), row['Abundance']
            sentence = f"This sample {'does' if present else 'does not'} contain {target} at abundance {abu}%."
            c.drawString(inch, y_pos, sentence)
            y_pos -= 0.5*inch

        # Add spacing between samples
        y_pos -= 0.5*inch

    # Save the PDF
    c.save()


# Generate the report
generate_pdf_report(data_filename=data_filename, metadata_filename=metadata)
