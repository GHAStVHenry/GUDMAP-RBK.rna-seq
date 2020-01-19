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
            rep=metaFile["Replicate_RID"].unique()[0]
            print(rep)
        if (len(metaFile[metaFile["File_Type"] == "FastQ"]) > 2):
            print("There are more then 2 fastq's in the metadata: " + " ".join(metaFile[metaFile["File_Type"] == "FastQ"].RID))
            exit(1)
    if (args.parameter == "ends"):
        if (metaFile.Paired_End.unique() == "Single End"):
            ends = "se"
        elif (metaFile.Paired_End.unique() == "Paired End"):
            ends = "pe"
        else:
            ends = "uk"
        print(ends)
    if (args.parameter == "endsManual"):
        if (len(metaFile[metaFile["File_Type"] == "FastQ"].all()) == 1):
            endsManual = "se"
        elif (len(metaFile[metaFile["File_Type"] == "FastQ"].all()) == 2):
            endsManual = "pe"
        print(endsManual)
    if (args.parameter == "stranded"):
        if (metaFile.Has_Strand_Specific_Information.unique() == "yes"):
            stranded = "stranded"
        elif (metaFile.Has_Strand_Specific_Information.unique() == "no"):
            stranded = "unstranded"
        else:
            print("Stranded metadata not match expected options: " + metaFile.Has_Strand_Specific_Information.unique())
            exit(1)
        print(stranded)
    if (args.parameter == "spike"):
        if (metaFile.Used_Spike_Ins.unique() == "yes"):
            spike = "yes"
        elif (metaFile.Used_Spike_Ins.unique() == "no"):
            spike = "no"
        else:
            print("Spike-ins metadata not match expected options: " + metaFile.Used_Spike_Ins.unique())
            exit(1)
        print(spike)
    if (args.parameter == "specie"):
        if (metaFile.Species.unique() == "Mus musculus"):
            specie = "Mus musculus"
        elif (metaFile.Species.unique() == "Homo sapiens"):
            specie = "Homo sapiens"
        else:
            print("Species metadata not match expected options: " + metaFile.Species.unique())
            exit(1)
        print(specie)

if __name__ == '__main__':
    main()