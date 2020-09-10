#!/usr/bin/env python3

import argparse
import pandas as pd
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

def main():
    refQuery=pd.read_json("refQuery.json")
    if refQuery["File_URL"].count() == 1:
       print(refQuery["File_URL"].values[0])
    else:
        raise Exception("Multple references found: \n%s" %
            refQuery["RID"])

if __name__ == '__main__':
    main()
