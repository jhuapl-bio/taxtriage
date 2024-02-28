
#!/usr/bin/env python
import gzip
import os
import urllib.request

def process_assembly_summary(filename):
    processed_taxids = {}
    output_lines = []

    with open(filename, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue  # Skip header lines

            cols = line.strip().split('\t')
            taxid = cols[5]
            refseq_category = cols[4]
            assembly_level = cols[11]
            accession = cols[0]

            # Check if this taxid has already been processed
            if taxid in processed_taxids:
                continue

            # Check for reference, representative or complete genome
            if refseq_category in ['reference genome', 'representative genome'] or assembly_level == 'Complete':
                processed_taxids[taxid] = True
                output_lines.append(line)
            else:
                # If none of the conditions are met, store the line to add later if the taxid is not found with higher priority
                if taxid not in processed_taxids:
                    processed_taxids[taxid] = line

    # Add the first seen accession for taxids that weren't matched with higher priority
    for taxid, line in processed_taxids.items():
        if line is not True:
            output_lines.append(line)

    return output_lines

def download_files(output_lines, download_dir):
    downloaded_files = []
    for line in output_lines:
        cols = line.strip().split('\t')
        ftp_url = cols[19]
        taxid = cols[5]
        if ftp_url:
            download_url = f"{ftp_url}/{os.path.basename(ftp_url)}_genomic.fna.gz"
            local_filename = os.path.join(download_dir, f"{taxid}.fna.gz")

            print(f"Downloading {download_url} to {local_filename}")
            try:
                urllib.request.urlretrieve(download_url, local_filename)
                downloaded_files.append((local_filename, taxid))
            except Exception as e:
                print(f"Error downloading {download_url}: {e}")

    return downloaded_files

def append_taxid_to_fasta(downloaded_files):
    for file_path, taxid in downloaded_files:
        modified_content = []

        with gzip.open(file_path, 'rt') as file:
            for line in file:
                if line.startswith('>'):
                    modified_content.append(f">{taxid}_{line[1:]}")
                else:
                    modified_content.append(line)

        modified_file_path = file_path.replace('.gz', '')
        with open(modified_file_path, 'w') as modified_file:
            modified_file.writelines(modified_content)

input_file="/Users/merribb1/Downloads/assembly_summary_refseq.txt"
output_file="/Users/merribb1/Downloads/unique_taxids_output.txt"
download_dir="/Users/merribb1/Downloads/genomes/ncbi"

# Ensure download directory exists
os.makedirs(download_dir, exist_ok=True)

output_lines = process_assembly_summary(input_file)
with open(output_file, 'w') as out_file:
    out_file.writelines(output_lines)

downloaded_files = download_files(output_lines, download_dir)
append_taxid_to_fasta(downloaded_files)
