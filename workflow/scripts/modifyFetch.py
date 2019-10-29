#!/usr/bin/env python3

import argparse
import pandas as pd
import re

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fetchFile',help="The fetch file from bdgap.zip.",required=True)
    args = parser.parse_args()
    return args

def main():
    args = get_args()
    fetchFile = pd.read_csv(args.fetchFile+"/fetch.txt",sep="\t",header=None)
    fileFile = pd.read_csv(args.fetchFile+"/data/File.csv",sep=",",header=0)
    fileFile_filtered = fileFile[fileFile["File_Type"]=="FastQ"]
    fetchFile_filtered = fetchFile[fetchFile[2].str[-9:]==".fastq.gz"]
    fetchFile_filtered_renamed = fetchFile_filtered
    for i in fileFile_filtered["File_Name"]:
        fetchFile_filtered_renamed[2][fetchFile_filtered_renamed[2].str.contains(i,regex=False)] = fetchFile_filtered_renamed[2][fetchFile_filtered_renamed[2].str.contains(i,regex=False)].values[0].replace(re.sub("\.R.\.fastq\.gz","",i),fileFile_filtered["Replicate_RID"][fileFile_filtered["File_Name"]==i].values[0])
    fetchFile_filtered_renamed.to_csv(args.fetchFile+"/fetch.txt",sep="\t",header=False,index=False)

if __name__ == '__main__':
    main()