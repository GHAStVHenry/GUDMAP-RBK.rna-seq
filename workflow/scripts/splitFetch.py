#!/usr/bin/env python3

import argparse
import pandas as pd
import os

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fetchFile',help="The fetch file from bdgap.zip.",required=True)
    args = parser.parse_args()
    return args

def main():
    args = get_args()
    fetchFile = pd.read_csv(args.fetchFile+"/fetch.txt",sep="\t",header=None)
    fileFile = pd.read_csv(args.fetchFile+"/data/File.csv",sep=",",header=0)
    replicateRID = fileFile.Replicate_RID.unique()
    fetchArray = {i:fileFile.URI[(fileFile.Replicate_RID == i) & (fileFile.File_Type == "FastQ")] for i in replicateRID}
    for i in replicateRID:
        if not os.path.exists(i):
            os.mkdir("Replicate_"+i)
        fetchFile[fetchFile[0].str.contains('|'.join(fetchArray[i]))].to_csv("Replicate_"+i+"/fetch.txt",sep="\t",header=False,index=False)

if __name__ == '__main__':
    main()