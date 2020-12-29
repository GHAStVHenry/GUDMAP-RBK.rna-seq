#!/usr/bin/env python3

import argparse
import pandas as pd
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--studyRID',
                        help="The study RID.", required=True)
    args = parser.parse_args()
    return args


def main():
    args = get_args()
    studyRID = pd.read_json(args.studyRID+"_studyRID.json")
    if studyRID["RID"].count() > 0:
        studyRID["RID"].to_csv(
            args.studyRID+"_studyRID.csv", header=False, index=False)
    else:
        raise Exception("No associated replicates found: %s" %
                        studyRID)


if __name__ == '__main__':
    main()
