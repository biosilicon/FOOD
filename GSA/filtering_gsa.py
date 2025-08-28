#!/usr/bin/env python3
"""
This script filters GSA RunInfo CSV files based on specific criteria.

Usage:
    python3 filtering_gsa.py -i <input_file> -o <output_file>

Arguments:
    -i, --input_file   Path to the GSA RunInfo CSV file
    -o, --output_file  Path to save the filtered CSV file
"""
import argparse
import subprocess
import pandas as pd

def filter_gsa_runinfo(input_file, output_file):
    '''
    Filter GSA RunInfo CSV file based on specific criteria and collect biosample info.
    
    Parameters:
    input_file (str): Path to the GSA RunInfo CSV file
    output_file (str): Path to save the filtered CSV file

    Returns:
    None
    '''
    def taxon_filter(df):
        # extract taxon id and get the lineage
        df_taxid = df["TaxID"]
        df_taxid = df_taxid.drop_duplicates()
        df_taxid = df_taxid.astype(int)
        df_taxid.to_csv("taxid.tmp", index=False, header=False)
        cmd = "taxonkit lineage taxid.tmp| tee lineage.tmp"
        subprocess.run(cmd, shell=True)
        df_lineage = pd.read_csv("lineage.tmp", sep="\t", header=None, index_col=None)
        df_lineage.columns = ["TaxID", "Lineage"]
        lineage_dict = dict(zip(df_lineage["TaxID"], df_lineage["Lineage"]))
        df["Lineage"] = df["TaxID"].map(lineage_dict)
        # filter rows where Lineage contains 'Metazoa'
        df = df[df["Lineage"].str.contains("Metazoa", case=False, na=False)].copy()
        df.loc[:, "If_chor"] = df["Lineage"].apply(lambda x: "Y" if "Chordata" in x else "N")
        df.loc[:, "If_vert"] = df["Lineage"].apply(lambda x: "Y" if "Vertebrata" in x else "N")
        df.loc[:, "If_mamm"] = df["Lineage"].apply(lambda x: "Y" if "Mammalia" in x else "N")
        # df = df.drop(columns=["Lineage"])
        cmd = "rm -f taxid.tmp lineage.tmp"
        subprocess.run(cmd, shell=True)
        return df

    df = pd.read_csv(input_file, sep=",", header=0, index_col=None, low_memory=False)
    df = taxon_filter(df)
    df.to_csv(output_file, sep=',', index=False, header=True, float_format='%.0f')
    print(f"Filtered data saved to {output_file}")


def main():
    parser = argparse.ArgumentParser(description="Filter GSA RunInfo CSV file")
    parser.add_argument("-i", "--input_file", help="Path to GSA RunInfo CSV file")
    parser.add_argument("-o", "--output_file", help="Path to save the filtered CSV file")
    
    args = parser.parse_args()
    
    filter_gsa_runinfo(args.input_file, args.output_file)

if __name__ == "__main__":
    main()