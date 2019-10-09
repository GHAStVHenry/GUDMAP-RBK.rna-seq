import argparse
import pandas as pd

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fetchFile',help="The fetch file from bdgap.zip.",required=True)
    args = parser.parse_args()
    return args

def main():
    args = get_args()
    fetch = pd.read_csv(args.fetchFile+"/fetch.txt",sep="\t",header=None)
    fetch_filtered = fetch[fetch[2].str[-9:]==".fastq.gz"]
    fetch_filtered.to_csv(args.fetchFile+"/fetch.txt",sep="\t",header=False,index=False)

if __name__ == '__main__':
    main()