#!/usr/bin/env python3

import argparse
import pandas as pd

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fetchFile',help="The fetch file from bdgap.zip.",required=True)
    args = parser.parse_args()
    return args

def main():
    args = get_args()
    fetchFile = pd.read_csv(args.fetchFile+"/fetch.txt",sep="\t",header=None)
    fetchFile_filtered = fetchFile[fetchFile[2].str[-9:]==".fastq.gz"]
    fetchFile_filtered.to_csv(args.fetchFile+"/fetch.txt",sep="\t",header=False,index=False)

if __name__ == '__main__':
    main()