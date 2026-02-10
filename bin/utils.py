import math
import os
from collections import defaultdict
import csv


def normalize_mapq(mapq_score, max_mapq=60, min_mapq=0):
        # Normalize the MAPQ score to be between 0 and 1
        probability = 1-(10 ** (-mapq_score / 10))
        # convert the min and max mapq
        # Normalize the MAPQ score to be between 0 and 1 based on the min and max
        return probability  # Ensure values stay between 0 and 1


def taxid_to_rank(tid: str, taxdump_dict: dict, target_rank: str = "species", max_hops: int = 2000):
    """
    Walk up taxonomy until we hit target_rank. Returns normalized taxid string or None.
    """
    if tid is None:
        return None
    tid = str(tid).strip()
    if not tid or tid == "0":
        return None

    cur = tid
    for _ in range(max_hops):
        node = taxdump_dict.get(cur)
        if not node:
            return cur  # best effort: return what we have
        if node.get("rank") == target_rank:
            return cur
        parent = node.get("parent") or node.get("parent_taxid") or node.get("parent_id")
        if not parent or str(parent) == cur or str(parent) == "0":
            return cur
        cur = str(parent)
    return cur

def calculate_var(read_counts):
    """
    Manually calculate variance for a list of read counts.
    read_counts: A list of aligned reads for each organism
    """
    n = len(read_counts)
    if n == 0:
        return 0  # Avoid division by zero if no organisms

    # Calculate the mean of the reads
    mean_reads = sum(read_counts) / n

    # the squared differences
    squared_diffs = [(x - mean_reads) ** 2 for x in read_counts]

    # variance
    variance = sum(squared_diffs) / n
    return variance


def load_matchfile(mapfile_path: str,
                   accession_col: int = 0,
                   taxid_col: int = 4,
                   desc_col: int = 2,
                   has_header: bool = True):
    """
    Load tab/CSV mapfile using column *indexes* (0-based).

    Parameters
    ----------
    accession_col : int
        Column index for accession (e.g. NC_XXXX)
    taxid_col : int
        Column index for taxid
    desc_col : int
        Column index for description / organism name
    has_header : bool
        Skip first row if True

    Returns
    -------
    accession_to_taxid : dict
        accession -> taxid (str)
    taxid_to_desc : dict
        taxid -> description (str)
    taxid_to_accessions : dict
        taxid -> set(accession)
    """

    accession_to_taxid = {}
    taxid_to_desc = {}
    taxid_to_accessions = defaultdict(set)

    if not mapfile_path or not os.path.exists(mapfile_path):
        return accession_to_taxid, taxid_to_desc, taxid_to_accessions

    # autodetect delimiter
    with open(mapfile_path, newline='') as fh:
        sample = fh.read(8192)
        fh.seek(0)
        delim = '\t' if '\t' in sample and sample.count('\t') >= sample.count(',') else ','
        reader = csv.reader(fh, delimiter=delim)

        if has_header:
            next(reader, None)

        for row in reader:
            # guard against short rows
            if len(row) <= max(accession_col, taxid_col, desc_col):
                continue

            acc = row[accession_col].strip()
            tax = str(row[taxid_col]).strip()
            desc = row[desc_col].strip() if desc_col is not None else ""

            if not acc or not tax:
                continue

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
    weight = 1 - normalized ** 4
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

