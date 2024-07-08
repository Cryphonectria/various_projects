#!/usr/bin/env python3
import re
import os
import glob
import sys
import pandas as pd
from functools import reduce
import argparse
from argparse import RawTextHelpFormatter

#------------- parser args
parser = argparse.ArgumentParser(
    formatter_class=RawTextHelpFormatter,
    description='report snippy results \n',
    usage='\n \n python3 %(prog)s -d /path/to/cluster -o summary_analysis.tsv')
parser.add_argument('-d', '--path', type=str, help='path to cluster directory')
parser.add_argument('-o', '--output', type=str, help='output file')
args = parser.parse_args()


# ------------------- functions
def GetSampleID(path):
    IDs = [name for name in os.listdir(path) if os.path.isdir(os.path.join(path, name))]
    return IDs

def LogFileLines(path):
    SNPinfo = []
    strings = ["SNP", "outdir"]
    for filename in glob.glob(path+"/*/snps.log"):
        with open(filename) as f:
            info = f.readlines()
            SNPs = [s for s in info if any(info in s for info in strings)]
            SNPs = [re.sub('.*--outdir | --R1.*|\n|Converted | to TAB format.','', word) for word in SNPs]
            SNPinfo.append(SNPs)

    return pd.DataFrame(SNPinfo[0:])


def GetSnpTable(path):
    filenames=glob.glob(path+"/*/snps.tab")
    list_of_dfs=[pd.read_csv(filename, sep="\t") for filename in filenames]
    list_of_dfs_edit=[]
    for dataframe, filename in zip(list_of_dfs, filenames):
        dataframe = dataframe.fillna("unknown")
        if dataframe.empty:
            dataframe["CHROM"]=[0]
        samplename=re.sub("\\/snps.tab","", filename)
        samplename=re.sub(".*\\/|","", samplename)
        dataframe["EFFECT"] = dataframe["EFFECT"].str.replace(" .*", "", regex=True)
        dataframe["EFFECT"] = dataframe["EFFECT"].str.replace("_", " ", regex=False)
        dataframe = dataframe.rename(columns={"TYPE": "Type"})
        dataframe["sample"] = samplename
        list_of_dfs_edit.append(dataframe)
    
    return pd.concat(list_of_dfs_edit)

def GetSNPInfo(path, column):
    df = GetSnpTable(path)
    df=df.groupby(['sample', column], as_index = False).size()
    return df

def summarize_results(directory, column):
    df = GetSNPInfo(directory, column)
    df_dict = df.set_index('sample') \
           .groupby(level=0) \
           .apply(lambda x: x.set_index(column).T.to_dict('records')[0]) \
           .to_dict()
    summary_df = pd.DataFrame.from_dict(df_dict, orient='index')
    summary_df['sample'] = summary_df.index
    
    return summary_df


# ----------------------- code
directory=args.path

types = summarize_results(directory, "Type")
ftypes = summarize_results(directory, "FTYPE")
effect = summarize_results(directory, "EFFECT")

snps = LogFileLines(directory)
snps.columns = ["sample", "SNPs"]
snps["SNPs"] = snps["SNPs"].str.replace(" SNPs", "", regex=False)

# merge all dfs
all_list = [snps, types, ftypes, effect]
all_df = reduce(lambda  left,right: pd.merge(left,right,on='sample', how='outer'), all_list).fillna(str(0)).transpose()
all_df, all_df.columns = all_df[1:] , all_df.iloc[0]

# remove decimal .0
all_df = all_df.astype('int64')

# rename index
all_df.index = all_df.index.str.replace("unknown_x", "unknown_variant_region", regex=False)
all_df.index = all_df.index.str.replace("unknown_y", "unknown_variant_effect", regex=False)


all_df.to_csv(args.output, sep='\t', encoding='utf-8', index = True)
