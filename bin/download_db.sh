db=$1
db_dir=$2
url=$3
checksum=$4

if [ -r $db_dir/$db ]; then
    echo "File found!"
    exit 0
fi

wget $url -O $db_dir/${db}.tar.gz
python $PWD/bin/checksum.py $db_dir/${db}.tar.gz $checksum
tar -xvzf $db_dir/${db}.tar.gz
rm $db_dir/${db}.tar.gz
mv -f $db* $db_dir/
