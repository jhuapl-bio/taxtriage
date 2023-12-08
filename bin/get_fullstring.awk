#!/usr/bin/gawk -f

# **********************************************************************
# Copyright (C) 2019 Johns Hopkins University Applied Physics Laboratory
#
# All Rights Reserved.
# For any other permission, please contact the Legal Office at JHU/APL.
# **********************************************************************


# INPUT 1: a file containing all of the taxids, one on each line

# INPUT 2: a joined nodes.dmp / names.dmp file (joined.dmp)
#   created from the NCBI taxonomy names.dmp / nodes.dmp files via:
#     join -j 1 -t $'\t' -o 0 2.2 1.2 2.3 \
#       <(grep "scientific name" names.dmp | cut -f1,3 | sort -k1,1) \
#       <(cut -f1,3,5 nodes.dmp | sort -k1,1) | sort -n > joined.dmp

# OUTPUT: the full taxonomic string for each taxid,
#  delimited by the "|" character
#  with each entry in the format: taxid;name(level)

#-------------------------------------------------

BEGIN {
	FS="\t"
}

#-------------------------------------------------

{
	if(NR == FNR) {
		LIST[NR] = $1;
	} else {
		parent[$1] = $2;
		name[$1] = $3;
		level[$1] = $4;
	}
}

#-------------------------------------------------

END {
	PROCINFO["sorted_in"] = "@ind_num_asc";
	for(INDEX in LIST) {
		TAXID = LIST[INDEX];
		out = sprintf("%s;%s(%s)", TAXID, name[TAXID], level[TAXID]);
		if(TAXID > 1) {
			PARENT = parent[TAXID];
			out = sprintf("%s;%s(%s)|%s", PARENT, name[PARENT], level[PARENT], out);
			while(PARENT > 1) {
				PARENT = parent[PARENT];
				out = sprintf("%s;%s(%s)|%s", PARENT, name[PARENT], level[PARENT], out);
			}
		}
		printf("%s\n", out);
	}
}
