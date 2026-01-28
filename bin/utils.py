import math
import os
from collections import defaultdict

def load_matchfile(mapfile_path: str,
                   accession_col: str = "Accession",
                   taxid_col: str = "TaxID",
                   desc_col: str = "Description"):
    """
    Load tab/CSV mapfile with columns like Accession, TaxID, Description.
    Returns:
      accession_to_taxid: dict accession -> taxid (str)
      taxid_to_desc: dict taxid -> description (str)
      taxid_to_accessions: dict taxid -> set(accession)
    If file missing or empty, returns empty dicts.
    """
    import csv
    accession_to_taxid = {}
    taxid_to_desc = {}
    taxid_to_accessions = defaultdict(set)

    if not mapfile_path or not os.path.exists(mapfile_path):
        return accession_to_taxid, taxid_to_desc, taxid_to_accessions

    # try to autodetect delimiter (tab or csv)
    with open(mapfile_path, newline='') as fh:
        sample = fh.read(8192)
        fh.seek(0)
        # choose delimiter
        delim = '\t' if '\t' in sample and sample.count('\t') >= sample.count(',') else ','
        reader = csv.DictReader(fh, delimiter=delim)
        # allow alternative column names (case-insensitive)
        header_map = {h.lower(): h for h in reader.fieldnames} if reader.fieldnames else {}

        # find actual columns to use
        acc_col = header_map.get(accession_col.lower(), accession_col) if header_map else accession_col
        tax_col = header_map.get(taxid_col.lower(), taxid_col) if header_map else taxid_col
        desc_col_use = header_map.get(desc_col.lower(), desc_col) if header_map else desc_col

        for row in reader:
            acc = row.get(acc_col) or row.get(accession_col) or None
            tax = row.get(tax_col) or row.get(taxid_col) or None
            desc = row.get(desc_col_use) or row.get(desc_col) or ""
            if not acc or not tax:
                # skip incomplete lines
                continue
            acc = acc.strip()
            tax = str(tax).strip()
            accession_to_taxid[acc] = tax
            taxid_to_desc.setdefault(tax, desc)
            taxid_to_accessions[tax].add(acc)

    return accession_to_taxid, taxid_to_desc, taxid_to_accessions

# Function to apply weights and format the result (assuming `format_non_zero_decimals` is defined)
def apply_weight(value, weight):
    try:
        # Convert the value to float before multiplying by the weight
        value = float(value)
        return value * weight
    except (TypeError, ValueError):
        # If value cannot be converted to a float, return 0 or handle as needed
        return 0

def logarithmic_weight(breadth, min_breadth=1e-3):
    """
    Returns a weight between 0 and 1 for a given breadth value.
    - When breadth == min_breadth, the weight is 1.
    - When breadth == 1, the weight is 0.
    - Values in between are mapped logarithmically.

    Parameters:
    breadth (float): The raw breadth value (expected to be between 0 and 1).
    min_breadth (float): The minimum expected breadth (should be > 0 to avoid log(0)).
    """
    # Clamp breadth to the [min_breadth, 1] range
    breadth = max(min_breadth, min(breadth, 1.0))

    # Normalize the log value
    normalized = (math.log(breadth) - math.log(min_breadth)) / (0 - math.log(min_breadth))
    weight = 1 - normalized
    return weight


def format_non_zero_decimals(number):
    # Convert number from the scientific notation to the tenth digit
    number = "{:.10f}".format(number).rstrip('0').rstrip('.')
    # Convert the number to a string
    num_str = str(number)
    if '.' not in num_str:
        # If there's no decimal point, return the number as is
        return number
    else:
        # Split into integer and decimal parts
        integer_part, decimal_part = num_str.split('.')
        non_zero_decimals = ''.join([d for d in decimal_part if d != '0'][:2])
        # Count leading zeros in the decimal part
        leading_zeros = len(decimal_part) - len(decimal_part.lstrip('0'))
        # Construct the new number with two significant decimal digits
        formatted_number = f"{integer_part}.{('0' * leading_zeros) + non_zero_decimals}"
        # Convert back to float, then to string to remove trailing zeros
        return str(float(formatted_number))

