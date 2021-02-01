#!/usr/bin/env python3

import argparse
import pandas as pd
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-r', '--repRID', help="The replicate RID.", required=True)
    parser.add_argument('-m', '--metaFile',
                        help="The metadata file to extract.", required=True)
    parser.add_argument('-p', '--parameter',
                        help="The parameter to extract.", required=True)
    args = parser.parse_args()
    return args


def main():
    args = get_args()
    metaFile = pd.read_csv(args.metaFile, sep=",", header=0)

    # Check replicate RID metadata from 'File.csv'
    if (args.parameter == "repRID"):
        if (len(metaFile.Replicate_RID.unique()) > 1):
            print("There are multiple replicate RID's in the metadata: " +
                  " ".join(metaFile.Replicate_RID.unique()))
            exit(1)
        if not (metaFile.Replicate_RID.unique() == args.repRID):
            print("Replicate RID in metadata does not match run parameters: " +
                  metaFile.Replicate_RID.unique() + " vs " + args.repRID)
            exit(1)
        else:
            rep = metaFile["Replicate_RID"].unique()[0]
            print(rep)

    # Check experiment RID metadata from 'Experiment.csv'
    if (args.parameter == "expRID"):
        if (len(metaFile.Experiment_RID.unique()) > 1):
            print("There are multiple experiment RID's in the metadata: " +
                  " ".join(metaFile.Experiment_RID.unique()))
            exit(1)
        else:
            exp = metaFile["Experiment_RID"].unique()[0]
            print(exp)

    # Check study RID metadata from 'Experiment.csv'
    if (args.parameter == "studyRID"):
        if (len(metaFile.Study_RID.unique()) > 1):
            print("There are multiple study RID's in the metadata: " +
                  " ".join(metaFile.Study_RID.unique()))
            exit(1)
        else:
            study = metaFile["Study_RID"].unique()[0]
            print(study)

    # Get endedness metadata from 'Experiment Settings.csv'
    if (args.parameter == "endsMeta"):
        endsMeta = metaFile.Paired_End.unique()[0]
        print(endsMeta)

    # Get strandedness metadata from 'Experiment Settings.csv'
    if (args.parameter == "stranded"):
        stranded = metaFile.Strandedness.unique()[0]
        print(stranded)

    # Get spike-in metadata from 'Experiment Settings.csv'
    if (args.parameter == "spike"):
        spike = metaFile.Used_Spike_Ins.unique()[0]
        print(spike)

    # Get species metadata from 'Experiment.csv'
    if (args.parameter == "species"):
        species = metaFile.Species.unique()[0]
        print(species)

    # Get read length metadata from 'Experiment Settings.csv'
    if (args.parameter == "readLength"):
        readLength = metaFile.Read_Length.unique()
        print(str(readLength).strip('[]'))


if __name__ == '__main__':
    main()
