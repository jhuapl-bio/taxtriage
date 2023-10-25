import hashlib
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("file", help="File for which you want to create a hash")
parser.add_argument(
    "hashvalue", help="Hash value to match against the calculated hash")
args = parser.parse_args()

file = args.file
hashval = args.hashvalue

algo = hashlib.sha256()

with open(file, "rb") as f:
    for byte_block in iter(lambda: f.read(4096), b""):
        algo.update(byte_block)
    digest = algo.hexdigest()
    print(digest)
exit(0) if digest == hashval else exit(
    "Hash values do not match. Data may be corrupted!")
