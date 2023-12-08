import pandas as pd
import argparse
import os
import matplotlib.pyplot as plt

# Load the TSV file into a DataFrame
parser = argparse.ArgumentParser()
parser.add_argument("-d", metavar="depth_path",
                    required = True,
                    help="Depth file for plots")
parser.add_argument("-i", metavar="file_path",
                    required = True,
                    help="Input Stats sheet for coverage")
args = parser.parse_args()

file_path = args.i
depth_file_path = args.d
data = pd.read_csv(file_path, sep='\t')

# Extract taxid from the #rname column by splitting on '|'
data['taxid'] = data['#rname'].apply(lambda x: x.split('|')[1])

# Calculate the full length for each row as endpos - startpos
data['full_length'] = data['endpos'] - data['startpos']

# Calculate the weighted mean depth and coverage for each taxid
# Weighted by the full length of the row
grouped_data = data.groupby('taxid').apply(
    lambda x: pd.Series({
        'weighted_meandepth': (x['meandepth'] * x['full_length']).sum() / x['full_length'].sum(),
        'weighted_coverage': (x['coverage'] * x['full_length']).sum() / x['full_length'].sum(),
        'total_length': x['full_length'].sum(),
        'total_covbases': x['covbases'].sum()
    })
)
filename_with_extension = os.path.basename(file_path)
samplename = os.path.splitext(filename_with_extension)[0]


# Calculate the percentage of covered bases
grouped_data['percentage_covbases'] = (grouped_data['total_covbases'] / grouped_data['total_length']) * 100

# Reset index to make 'taxid' a column again
grouped_data = grouped_data.reset_index()

# Display the resulting DataFrame
grouped_data[['taxid', 'weighted_meandepth', 'weighted_coverage', 'percentage_covbases']]
depth_data = pd.read_csv(depth_file_path, sep='\t')
depth_data = pd.read_csv(depth_file_path, sep='\t', header=None, names=['rname', 'position', 'depth'])

# Extract the taxid and the reference accession from the rname column
depth_data[['discard', 'taxid', 'accession', 'ref']] = depth_data['rname'].str.split('|', expand=True)

# Now we can drop the 'rname' and 'discard' columns as they are no longer needed
depth_data.drop(['rname', 'discard'], axis=1, inplace=True)

# Get unique reference accessions within each taxid
unique_accessions = depth_data.groupby('taxid')['accession'].unique()

unique_accessions.head()



# Check the unique taxids in the depth data
unique_taxids = depth_data['taxid'].unique()
unique_taxids, len(unique_taxids)
# Function to plot and save individual SVG plots for each taxid
def plot_and_save_accessions_for_taxid(taxid, accessions, depth_data):
    # Determine the number of unique accessions for the taxid
    num_accessions = len(accessions)

    # Create a figure with a subplot for each accession
    fig, axes = plt.subplots(nrows=num_accessions, ncols=1, figsize=(10, 3 * num_accessions), sharex=True)

    # If there's only one accession, wrap axes in a list for consistent indexing
    if num_accessions == 1:
        axes = [axes]

    # Plot the depth data for each accession
    for i, accession in enumerate(accessions):
        # Filter the depth data for the current accession
        accession_data = depth_data[depth_data['accession'] == accession]

        # Plot with vertical lines on the corresponding subplot
        axes[i].vlines(accession_data['position'], ymin=0, ymax=accession_data['depth'], color='blue')
        axes[i].set_title(f'Accession {accession}')
        axes[i].set_ylabel('Depth')

    # Set common labels
    plt.xlabel('Position')
    fig.suptitle(f'Depth Profiles for TaxID {taxid}', fontsize=16)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Adjust the layout to make room for the suptitle

    # Save the plot as an SVG file
    svg_filename = os.path.join(os.path.dirname(args.d), f"{samplename}_{taxid}_depth_profile_taxid.svg")
    plt.savefig(svg_filename, format='svg' )
    return svg_filename

# Create and save a plot for each unique taxid
test_taxid = unique_accessions.index[0]
test_accessions = unique_accessions[test_taxid]
svg_test_filename = plot_and_save_accessions_for_taxid(test_taxid, test_accessions, depth_data)
svg_test_filename
