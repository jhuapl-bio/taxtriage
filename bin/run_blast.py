# to run
# conda install biopython
# run_blast.py /Users/ernluaw1/Desktop/testt/draft_assembly_read1.fasta /Users/ernluaw1/Desktop/testt
import os
from Bio.Blast.Applications import NcbiblastnCommandline
import sys

# assembly path + fasta file
assembly = str(sys.argv[1])

# directory path for output
output_directory = str(sys.argv[2])

# remove path and file extension and maintain base file name
out_file_base = assembly.split('/')[-1].split('.')[0]

# output file
out = '{0}/{1}.txt'.format(output_directory, out_file_base)
os.makedirs(output_directory, exist_ok=True)

# output format for blast - read_id, accession, evalue, bitscore, organism
outfmt = '6 qseqid sseqid evalue bitscore stitle'

# set up blast
cline = NcbiblastnCommandline(query=assembly, db="nt", out=out, outfmt=outfmt,
    evalue=0.001, max_target_seqs=1, remote=True, ungapped=True)

print("running blast")
# run blast
stdout, stderr = cline()

# format output file with header
with open(out) as f:
    mod = f.read()

# header
header = "read_id\taccession\tevalue\tbitscore\torganism\n"

# add header to top of the blast file
new_mod = header + mod

# overwrites original file with header
with open(out, 'w') as f:
    f.write(new_mod)
