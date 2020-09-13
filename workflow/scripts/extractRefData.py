#!/usr/bin/env python3

import argparse
import pandas as pd
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--returnParam',help="The parameter to return (URL or MD5).",required=True)
    args = parser.parse_args()
    return args

def main():
    args = get_args()
    refQuery=pd.read_json("refQuery.json")
    if refQuery["File_URL"].count() == 1:
        if args.returnParam == "URL":
            print(refQuery["File_URL"].values[0])
        elif args.returnParam == "fName":
            print(refQuery["File_Name"].values[0])
        elif args.returnParam == "MD5":
            print(refQuery["File_MD5"].values[0])
    else:
        raise Exception("Multple references found: \n%s" %
            refQuery["RID"])

if __name__ == '__main__':
    main()
