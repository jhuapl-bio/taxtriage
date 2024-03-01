import argparse
import os
import json
import glob
import re
import pandas as pd
parser = argparse.ArgumentParser()
parser.add_argument("-i",
                    required = True,
                    help="Input Samplesheet filename")
parser.add_argument(
    "-t",
    type=str,
    required=False,
    default=None,
    help="Where is the top hits report? It is a directory of files with top_report.tsv with samplename in the front"
)
parser.add_argument(
    "-c",
    type=str,
    required=True,
    help="Where is the configuration file? it is a json that contains the publish directories and filenames for each test"
)
parser.add_argument(
    "-w",
    type=str,
    required=True,
    help="Where is the workding directory where your outputs are? "
)
# This script works at the species level to calls tests on each species/taxid associated with a samplesheet list of samples.
# it reqres a samplesheet and a taxtriage working outdir with identifiable structures according to the -c argument (json).
# see assets/pass.json for an example format of the json file on the simulated test data used in taxtriage's test profile



args = parser.parse_args()

file = args.i
# tests = args.t



def val_test(config):
    match = config['match']
    sample = config['sample']
    regex = config['regex']
    directory = config['directory']
    value = config['value']
    test = config['test']
    # if regex is true then check the -w directory for the file using the match. should return a file
    if regex:
        filename = os.path.join(args.w, directory, sample+match)
        # convert filename to regex format
        # use glob to find the file

        filenames = glob.glob(filename)
        # if the file is not found then return false
        if len(filenames) == 0:
            print("File not found, "+filename)
            return False
        else:
            filename = filenames[0]
    else:
        filename = os.path.join(args.w, directory, match)
    # # read the file into a dataframe
    df = pd.read_csv(filename, sep="\t")
    # file all df on rank == S and abundance greater or equal to value
    df = df[df['rank'].str.contains("^S$")]
    df = df[df['abundance'] >= value]

    # top_hits = os.path.join(top, sample + ".top_report.tsv")
    # # read top_hits into a dataframe
    # df_top = pd.read_csv(top_hits, sep="\t")
    # # sort the dataframe by the abundance column and only select those with rank contains S and an optional number like S, S1, S4, etc.
    # df_top = df_top[df_top['rank'].str.contains("^S$")]
    # # sort on abundance
    # df_top = df_top.sort_values(by=['abundance'], ascending=False)
    # # For each sample loop through the taxids and check each of the tests defined as functions. Tests are defined in the config dict imported



def __main__():
    top = args.t
    working_dir = args.w
    config = args.c

    if not top :
        print("No top hits file provided, opting for default location at working directory / top / samplename.tsv")
        top = os.path.join(working_dir, "top")
    df = pd.read_csv(file, sep=",")
    samplenames = df['sample']
    data = df.to_dict(orient='records')
    # import the config -c argument as a dict
    with open(config) as f:
        config = json.load(f)
    for i in data:
        sample = i['sample']
        conf_one = next(test for test in config['tests'] if test['id'] == 'abundance')
        conf_one['sample'] = sample
        val_test(conf_one)

    return
__main__()
