#!/usr/bin/env python3
"""
Author:  Lea Stauber
Purpose: check in list.mapped.taxo files how many reads of a specific species (e.g contaminations) are present
Version: 0.1
"""

import argparse
import re
import os
import glob
import pandas as pd

# -------------------------------------------------

def get_args():
    """Get the command-line arguments"""

    parser = argparse.ArgumentParser(description='get number of mapped reads for a species')
    parser.add_argument('-dir', '--directory', help='Path to LorcanX output', type=str)
    parser.add_argument('-s', '--species', help='species name, e.g. Streptococcus_agalactiae ', type=str)
    parser.add_argument('-out', '--output', help='output file', type=str)
    return parser.parse_args()

# -------------------------------------------------

def find_taxo_files(directory):
    """
    Find all files with the name 'list.mapped.taxo'
    """
    taxo_files = glob.glob(f"{directory}/**/01_mapping_reads_to_refDB/list.mapped.taxo", recursive=True)
    return taxo_files

# -------------------------------------------------




def find_species(taxofiles, species):
    """
    fetch species per BC
    """
    result_dict = {}
    for file_path in taxofiles:
        with open(file_path, 'r') as f:
            our_lines = f.readlines()
            BC_match = re.search(r'BC\d{2}', file_path)
            
            if BC_match:
                BC = BC_match.group()
            else:
                BC = "UnknownBC"
            
            lines_with_species = []
            for l in our_lines:
                if species in l:
                    cleaned_line = re.sub(r"=.*", "", l.strip())  # Remove '=.*' pattern from the line
                    lines_with_species.append(cleaned_line)
            
            # Store the current BC and line in the result dictionary
            result_dict[BC] = lines_with_species
    
    return result_dict

# -------------------------------------------------

def dict_to_df(dict):
    df = pd.DataFrame.from_dict(dict, orient='index')
    df.reset_index(inplace=True)
    df = pd.melt(df, id_vars=['index'], var_name='BC', value_name='species')
    df = df.drop(columns=['BC'])
    df.rename(columns={'index': 'BC'}, inplace=True)
    df = df.dropna(subset=['species'])
    # separate species column
    df[['reads', 'mapped_species']] = df['species'].str.split(" ", expand=True)
    df.drop(columns=['species'], inplace=True)
    # Transform column to numeric
    df["reads"] = pd.to_numeric(df["reads"], errors='coerce')
    df.sort_values(by=['BC', 'reads'], inplace=True, ascending=[True, False])
    return df


def main():
    args = get_args()
    files=find_taxo_files(args.directory)
    species=find_species(files, args.species)
    df=dict_to_df(species)
    df.to_csv(args.output, index=False, sep="\t")



if __name__ == '__main__':
    main()