#!/usr/bin/env python3
"""
Author:  Lea Stauber
Purpose: Summary statistics for RBK amplicon reads
Version: 0.1
"""

import argparse
import re
import os
from dataclasses import dataclass
import glob
import subprocess
from Bio import SeqIO
import gzip
import matplotlib.pyplot as plt
import numpy as np
import statistics as stat
import pandas as pd

# --------------------------------------------------
def get_args():
    """Get the command-line arguments"""

    parser = argparse.ArgumentParser(description='Fastq summary stats')
    parser.add_argument('-f', '--file', default='sample-id.txt', help='RunID and BC info', type=str)
    return parser.parse_args()


# --------------------------------------------------

def get_runname(sample_file):
    """ Get runname from sample_id.txt file """
    with open(sample_file) as f:
        run_name = f.readline().strip('\n')
        run_name = re.sub('RunName:', '', run_name)
        
        return run_name

# --------------------------------------------------
def read_lengths(fastq):
    """get read lengths as list"""
    try:
        sizes=[]
        with gzip.open(fastq[0], "rt") as fq: # fastq[0] --> returns path as string and not as list
            for record in SeqIO.parse(fq, "fastq"):
                size = len(record)
                sizes.append(size)
        return sizes

    except IndexError as exc:
        print(f"Error: Concatinated fastq file not found!\n"
                f"Returned: {exc}")
    except PermissionError as e:
        print(f"Error: cannot read file from barcode directory\n"
            f"Returned: {e}")

# --------------------------------------------------

def read_lengths_hist(sizes, object):
    """plot read length distribution"""
    try:
        if os.path.exists(object.dir_path):
            y, x, _ = plt.hist(sizes, bins=30, color="grey", alpha = 0.7)
            plt.title(
                "%i sequences\nLengths %ibp to %ibp" % (len(sizes), min(sizes), max(sizes))
            )
            plt.axvline(x=object.amplicon_length, color = "red", ls='--', lw=2)
            plt.text(object.amplicon_length-60, y.max()/2, str(object.amplicon_length)+"bp", rotation=90)
            plt.xlabel("Sequence length (bp)")
            plt.ylabel("Count")
            plt.savefig(os.path.join(object.dir_path, object.BC + '.png'), dpi=200)
            plt.close()
    
    except (TypeError, AttributeError) as exc:
        print(f"Error: barcode directory not found, no data for plotting available!\n"
                f"Returned: {exc}")
    except PermissionError as e:
        print(f"Error: cannot write to barcode directory\n"
            f"Returned: {e}")
    except UnboundLocalError:
        print(f"Warning: gene name and expected amplicon size not specified in database, cannot create histogram for " + object.BC)

# --------------------------------------------------

def summary_statistics(sizes, object):
    """create read length summary statistics as table"""
    try:
        if os.path.exists(object.dir_path):
            summary_stats = {\
                "BC":object.BC, \
                "N_reads": [len(sizes)], \
                "min_length": [min(sizes)], \
                "max_length": [max(sizes)], \
                "mean_length": [int(round(stat.mean(sizes),0))]}
            summary_df = pd.DataFrame(summary_stats)
            summary_df.to_csv(os.path.join(object.dir_path, object.BC + '_stats.txt'), sep="\t", index=False)

    except PermissionError as exc:
        print(f"Error: cannot write to barcode directory \n"
            f"Refurned: {exc}")

# --------------------------------------------------

def concat_summarystats():
    """concatinate summary statistics tables from summary_statistics()"""
    wd=os.getcwd()
    df_list = []
    bad_files = []
    filenames=glob.glob(wd + '/BC*/*_stats.txt')
    for f in filenames:
        try:
            df_list.append(pd.read_csv(f, sep="\t"))
        except:
            bad_files.append(f)
    df = pd.concat(df_list)
    df.to_csv(wd + "/summary_stats.txt", sep="\t", index=False)


# --------------------------------------------------
# --------------------------------------------------
@dataclass
class BarcodeAttributes:
    BC: str
    IDIS: str
    species: str
    gene: str
    file: str ## maybe replace later with run name directly
    
    @property
    def amplicon_length(self):
        """specify amplicon length by gene name"""
        if  "VP1" in self.gene:
            length = 400
        if  "VP2" in self.gene:
            length = 550
        if  "S-Region" in self.gene:
            length = 1000
        return length


    @property
    def dir_path(self):
        """make directory paths for each BC in working directory"""
        wd = os.getcwd()
        path = os.path.join(wd, self.BC)
        return path

    @property
    def get_raw_fastq(self):
        """Get list of raw fastq files for BC"""
        rn = get_runname(self.file)
        idx = re.sub("BC", "", self.BC)
        BC_tmp=glob.glob("/storage/tmp" +"/*" + rn + ".tmp/*/*" + idx + "/*.fastq.gz") ## remove hard coded tmp path
        return BC_tmp

    @property
    def concat_raw_fastq(self):
        try:
            if not os.path.exists(self.dir_path):
                os.mkdir(self.dir_path)

            subprocess.run(f'cat {" ".join(self.get_raw_fastq)} > {"".join(self.dir_path) + "/" + self.BC + ".fastq.gz"}',\
                shell=True, check=True, timeout=60)
                
        except subprocess.CalledProcessError as exc:
            print(
                f"Concatinate fastq files failed. Check if input fastqs and output directory exists\n"
                f"Returned {exc.returncode}\n{exc}"
                )




def main():
    args = get_args()
    try:
        with open(args.file) as f:
            for row in f.readlines()[2:]:
                lines = row.strip().split(':')
                BC = BarcodeAttributes(lines[0], lines[1], lines[2], lines[3], args.file)

                if not os.path.exists(BC.dir_path):
                    BC.concat_raw_fastq
                    fastq = glob.glob(BC.dir_path + "/*.fastq.gz")
                    sizes = read_lengths(fastq)
                    summary_statistics(sizes, BC)
                    read_lengths_hist(sizes, BC)
                    
        concat_summarystats()

    except ValueError:
        print("ERROR: cannot parse input file")



if __name__ == '__main__':
    main()

