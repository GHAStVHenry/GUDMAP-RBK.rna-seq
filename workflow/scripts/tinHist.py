#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--repRID',help="The replicate RID.",required=True)
    args = parser.parse_args()
    return args

def main():
    args = get_args()
    tin = pd.read_csv(args.repRID + '.sorted.deduped.tin.xls',sep="\t",header=0)
    
    hist = pd.cut(tin['TIN'],bins=pd.interval_range(start=0,freq=10,end=100,closed='right')).value_counts(sort=False)
    labels = ["{0} - {1}".format(i, i + 9) for i in range(1, 100, 10)]
    #labels[0] = '0 - 10'
    binned = tin.assign(Bins=lambda x: pd.cut(tin['TIN'],range(0,105,10),labels=labels,include_lowest=False,right=True))
    binned['chrom'] = binned['chrom'] = binned['chrom'].replace('chr1','chr01')
    binned['chrom'] = binned['chrom'].replace('chr2','chr02')
    binned['chrom'] = binned['chrom'].replace('chr3','chr03')
    binned['chrom'] = binned['chrom'].replace('chr4','chr04')
    binned['chrom'] = binned['chrom'].replace('chr5','chr05')
    binned['chrom'] = binned['chrom'].replace('chr6','chr06')
    binned['chrom'] = binned['chrom'].replace('chr7','chr07')
    binned['chrom'] = binned['chrom'].replace('chr8','chr08')
    binned['chrom'] = binned['chrom'].replace('chr9','chr09')
    hist = pd.pivot_table(binned, values='geneID', index = 'Bins', columns = 'chrom', aggfunc=np.size)
    hist['TOTAL'] = hist.sum(axis=1)
    hist = hist[['TOTAL'] + [ i for i in hist.columns if i != 'TOTAL']]
    hist = hist.T.fillna(0.0).astype(int)
    #hist = hist.apply(lambda x: x/x.sum()*100, axis=1)
    hist.to_csv(args.repRID + '.tin.hist.tsv',sep='\t')
    medFile = open(args.repRID + '.tin.med.csv',"w")
    medFile.write(str(round(tin['TIN'][(tin['TIN']!=0)].median())))
    medFile.close()

if __name__ == '__main__':
    main()
