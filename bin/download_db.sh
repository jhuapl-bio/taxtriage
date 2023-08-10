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

wget $url -O $db_dir/${db}.tar.gz
tar -xvzf $db_dir/${db}.tar.gz
rm $db_dir/${db}.tar.gz
mv -f $db* $db_dir/
python3 $db_dir/../bin/checksum.py $db_dir/$db/taxo.k2d $checksum
