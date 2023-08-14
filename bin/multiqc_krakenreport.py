from multiqc.modules.base_module import BaseMultiqcModule
import logging
import pandas as pd
import os

# Initialise the logger
log = logging.getLogger('multiqc')

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Kraken2', anchor='kraken2',
        href="http://www.ccb.jhu.edu/software/kraken2/",
        info="is a system for assigning taxonomic labels to short DNA sequences, usually obtained through metagenomic studies.")

        # Define the file paths
        self.file_paths = self.find_log_files('kraken2')

        # Define the rank codes
        self.rank_codes = ["S", "G", "F", "O", "C", "P", "K", "D"]
        print(self.file_paths)
        # Process each file
        for file_path in self.file_paths:
            self.parse_kraken2_report(file_path)

    def parse_kraken2_report(self, file_path):
        # Read the file into a DataFrame
        df = pd.read_csv(file_path, sep="\t", header=None)

        # Rename the columns for clarity
        df.columns = ["percent_abundance", "reads_clade", "reads_taxon", "rank", "tax_id", "name"]

        # Filter the DataFrame to only include the specified rank codes
        df = df[df["rank"].isin(self.rank_codes)]

        # Sort the DataFrame by percent abundance, in descending order
        df = df.sort_values("percent_abundance", ascending=False)

        # Get the top 2 organisms for each rank
        top_organisms_df = df.groupby("rank").head(2)

        # Pivot the DataFrame to get the desired format
        pivot_df = top_organisms_df.pivot(index="rank", columns="name", values="percent_abundance")

        # Get the file name
        file_name = os.path.basename(file_path).split(".")[0]

        # Add the file name as a column
        pivot_df["file"] = file_name

        # Add to the general stats table
        self.general_stats_addcols(pivot_df, file_name)

        # Generate the report section
        self.add_section(
            name='Top Organisms',
            anchor='kraken2-top-organisms',
            plot=pivot_df
        )
