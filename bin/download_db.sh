#!/usr/bin/env bash

db=$1
url=$2
checksum=$3
# db_dir=$4

set -e
# if [ -r $db_dir/$db ]; then
#     echo "File found!"
#     python3 $db_dir/../bin/checksum.py $db_dir/$db/taxo.k2d $checksum
#     exit 0
# fi
echo "Downloading $db from $url $checksum"
wget $url -O ${db}.tar.gz
# curl -o ${db}.tar.gz $url

mkdir -p $db
tar -xvf ${db}.tar.gz -C $db
rm ${db}.tar.gz

if [[  ${db} == 'flukraken2' ]]; then
    mv -f ${db}/${db}/* ${db}/
fi
if [[  ${db} == 'minikraken2' ]]; then
    mv -f ${db}/minikraken2_v2_8GB_201904_UPDATE/* ${db}/
fi

# python3 $db_dir/../bin/checksum.py $db_dir/$db/taxo.k2d $checksum
