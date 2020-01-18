#!/usr/bin/env python3

import argparse
import pandas as pd

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--repRID',help="The replicate RID.",required=True)
    parser.add_argument('-m', '--metaFile',help="The metadata file to extract.",required=True)
    parser.add_argument('-p', '--parameter',help="The parameter to extract.",required=True)
    args = parser.parse_args()
    return args


def main():
    args = get_args()
    metaFile = pd.read_csv(args.metaFile,sep="\t",header=None)
    if (args.parameter == "repRID"):
        if (len(metaFile.Replicate_RID.unique()) > 1):
            #ERROR
        if not (metaFile.Replicate_RID == arg$repRID):
            #ERROR
        if (len(fileFile[fileFile["File_Type"] == "FastQ"].RID) > 2):
            #ERROR


if __name__ == '__main__':
    main()