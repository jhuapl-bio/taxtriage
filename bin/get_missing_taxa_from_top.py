import requests
import os
def read_file_1(file_path):
    """Reads the first file and extracts the 5th column."""
    with open(file_path, 'r') as file:
        return {line.split()[4] for line in file}

def read_file_2(file_path):
    """Reads the second file and extracts the 2nd column (delimited by '|')."""
    with open(file_path, 'r') as file:
        for line in file:
            f = line.split('\t')[0]
            # split on "|" and if index 1 is present then yield it
            if len(f.split("|")) > 1:
                yield f.split("|")[1]


def read_file_3(file_path):
    """Reads the third file and creates a mapping from ID to URL."""
    id_to_url = {}
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.split()
            # if lengths of columns is greater than 18 then yield the 5th and 19th column
            if len(parts) > 20:
                ftp_site = parts[20]
                fullid = os.path.basename(ftp_site) + '_genomic.fna.gz'
                ftp_site = ftp_site+'/'+fullid
                # print(ftp_site)
                id_to_url[parts[4]] = fullid   # Adjust indices if necessary

    return id_to_url

def check_url(url):
    """Checks if a URL is active."""
    try:
        response = requests.head(url, allow_redirects=True)
        return response.status_code == 200
    except requests.RequestException:
        return False



def calculate_presence_and_percentage(file_1_data, file_2_data):
    """Calculates the numbers present in file 1 but not in file 2 and the percentage."""
    unique_file_1 = set(file_1_data)
    unique_file_2 = set(file_2_data)

    not_in_file_2 = unique_file_1 - unique_file_2
    percentage_present = (len(unique_file_1) - len(not_in_file_2)) / len(unique_file_1) * 100

    return not_in_file_2, percentage_present

# Replace 'file1.txt' and 'file2.txt' with your actual file paths
file_1_data = read_file_1('/Users/merribb1/Downloads/14_Sample2_DNA.top_report.tsv')
file_2_data = read_file_2('/Users/merribb1/Downloads/27_Sample2_DNAf.txt')
id_to_url = read_file_3('/Users/merribb1/Downloads/assembly_summary_refseq.txt')



not_in_file_2, percentage_present = calculate_presence_and_percentage(file_1_data, file_2_data)

# Checking URLs for IDs not in file 2
valid_urls = 0
dead_urls = 0
for id in not_in_file_2:
    if id in id_to_url:
        if check_url(id_to_url[id]):
            valid_urls += 1
        else:
            dead_urls += 1

print(f"Percentage of numbers in File 1 present in File 2: {percentage_present}%")
print(f"Valid URLs: {valid_urls}")
print(f"Dead URLs: {dead_urls}")
