import re
def determine_priority_assembly(line):
    cols = line.strip().split("\t")
    if len(cols) > 7:
        name = cols[7]
        level = cols[4]
        if level == "representative genome" :
            return 0
        elif cols[4] == "reference genome"  :
            return 1
        elif cols[11] == "Complete Genome"  :
            return 2
        else:
            return 3
def format_description(id, description):
    delimiters = "[ ]"  # Split on underscore or comma
    description = description.split(',', 1)[0]
    linesplit = re.split(delimiters, description)
    id = id.replace(">", "")
    # remove all after first comma

    if len(linesplit) > 1:
        return id, " ".join(linesplit[1:])
    else:
        return id, description
def extract_strain_from_description(description):
    delimiters = "[ ]"  # Split on underscore or comma
    # regex match for strain, get value after space until next space
def parse_fasta_header(header):
    parts = header.split()
    organism_parts = []
    strain = None
    for part in parts:
        if part.lower() == "strain":
            strain_index = parts.index(part) + 1
            if strain_index < len(parts):
                strain = parts[strain_index]
            break
        organism_parts.append(part)
    organism_name = " ".join(organism_parts)
    return organism_name, strain
