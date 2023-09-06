db=$1
db_dir=$2
url=$3
checksum=$4

set -e
if [ -r $db_dir/$db ]; then
    echo "File found!"
    python3 $db_dir/../bin/checksum.py $db_dir/$db/taxo.k2d $checksum
    exit 0
fi
echo "Downloading $db from $url, $checksum, $db_dir"
wget $url -O $db_dir/${db}.tar.gz
mkdir -p $db_dir/$db
tar -xvzf $db_dir/${db}.tar.gz -C $db_dir/$db
rm $db_dir/${db}.tar.gz
python3 $db_dir/../bin/checksum.py $db_dir/$db/taxo.k2d $checksum
