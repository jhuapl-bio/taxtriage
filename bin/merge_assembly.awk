#!/usr/bin/awk -f

# Check if the line starts with ">"
BEGIN {FS=OFS="\t"}
NR==FNR {map[$1]=$2; next}
{
    if ($2 == "") {
        print $0, $3  # If assembly is empty, duplicate the third column
    } else {
        print $0, (map[$2] != "" ? map[$2] : $3)  # Otherwise, use mapped value if present, else default to third column
    }
}
