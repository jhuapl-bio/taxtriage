import ncbi_genome_download as ngd
import io
import contextlib




import logging
import re

# Setting up logging to capture output
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()

# List of species or strains
species_list = [
    "Staphylococcus aureus subsp. aureus NCTC 8325",
    "Klebsiella pneumoniae",
    "Metabacillus litoralis",
    "Pediococcus acidilactici",
    "Neisseria gonorrhoeae",
    "Fluviibacter phosphoraccumulans",
    "Diaphorobacter ruginosibacter",
    "Dickeya fangzhongdai"
]

def query_species(species_list, portion=1):
    species_string = ','.join(species_list)

    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer):
        if portion == 1:
            ngd.download(
                genera=species_string,
                flat_output=True,
                refseq_categories='reference',
                assembly_levels='complete',
                dry_run=True
            )
        elif portion == 2:
            ngd.download(
                genera=species_string,
                flat_output=True,
                assembly_levels='complete',
                dry_run=True
            )
        else:
            ngd.download(
                genera=species_string,
                flat_output=True,
                dry_run=True
            )

    output = buffer.getvalue()
    return output.splitlines()

def get_first_hit(entries, species):
    try:
        for entry in entries:
            if species.lower() in entry.lower():
                return entry
        return None
    except Exception as e:
        print(f"Error during processing for {species}: {e}")
        return None


def prioritize_queries(species_list):
    results = {}

    try:
        # First priority
        first_priority_hits = query_species(species_list, 1)
        results = {species: get_first_hit(first_priority_hits, species) for species in species_list}

        # Second priority for missing species
        missing_species = [s for s, hit in results.items() if hit is None]
        if missing_species:
            second_priority_hits = query_species(missing_species, 2)
            for species in missing_species:
                results[species] = get_first_hit(second_priority_hits, species)

        # Final priority for remaining missing species
        still_missing = [s for s, hit in results.items() if hit is None]
        if still_missing:
            final_hits = query_species(still_missing,3)
            for species in still_missing:
                results[species] = get_first_hit(final_hits, species)

    except Exception as e:
        print(f"Error during prioritized queries: {e}")
        results = {species: None for species in species_list}

    return results

# Perform prioritized queries
species_results = prioritize_queries(species_list)

# Output the results
for species, result in species_results.items():
    print(f"{species}: {result}")
