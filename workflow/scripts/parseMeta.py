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
    metaFile = pd.read_csv(args.metaFile,sep=",",header=0)
    if (args.parameter == "repRID"):
        if (len(metaFile.Replicate_RID.unique()) > 1):
            print("There are multiple replicate RID's in the metadata: " + " ".join(metaFile.Replicate_RID.unique()))
            exit(1)
        if not (metaFile.Replicate_RID.unique() == args.repRID):
            print("Replicate RID in metadata does not match run parameters: " + metaFile.Replicate_RID.unique() + " vs " + args.repRID)
            exit(1)
        else:
            print(metaFile["Replicate_RID"].unique())
        if (len(metaFile[metaFile["File_Type"] == "FastQ"]) > 2):
            print("There are more then 2 fastq's in the metadata: " + " ".join(metaFile[metaFile["File_Type"] == "FastQ"].RID))
            exit(1)


if __name__ == '__main__':
    main()