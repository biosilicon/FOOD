#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
This script is used to fetch and process GSA data. The processing is quite similar to that of SRA data,
including filtering, integrating biosample information, add refseq assembly information and estimate the depth.
'''
import os
import argparse
import subprocess
import logging
import requests
from bs4 import BeautifulSoup
import pandas as pd

class gsa_data_fetch:
    def __init__(self, input_file, output_dir, refinfo, overwrite=False):
        self.input_file = input_file
        self.output_dir = output_dir
        self.refinfo = refinfo
        self.overwrite = overwrite
        self.fetch_data()
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    
    def prefiltering(self) -> pd.DataFrame:
        '''
        Filter GSA RunInfo CSV file based on specific criteria and collect biosample info.
        '''
        df = pd.read_csv(self.input_file, sep=",", header=0, index_col=None, low_memory=False)
        if "Organization" in df.columns:
            df = df.drop(columns=["Organization"])
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
        cmd = "rm -f taxid.tmp lineage.tmp"
        subprocess.run(cmd, shell=True)
        # df.to_csv(f'{self.output_dir}/gsa_runinfo_taxon_filtered.csv', sep=',', index=False, header=True, float_format='%.0f')
        # logging.info(f"Filtered data saved to {self.output_dir}/gsa_runinfo_taxon_filtered.csv")
        return df

    def integrate_biosample(self, df: pd.DataFrame):
        '''
        Integrate the information by Biosample ID and add sample information and refseq assembly information
        '''

        def rearrange(df: pd.DataFrame) -> pd.DataFrame:
            grouped = df.groupby('BioSample')
            new_rows = []
            for name, group in grouped:
                row_count = len(group)
                bioproject = group['BioProject'].dropna().unique().tolist()
                taxid = group['TaxID'].dropna().unique().tolist()
                taxon = group['ScientificName'].dropna().unique().tolist()
                lineage = group['Lineage'].dropna().unique().tolist()
                if_chor = group['If_chor'].dropna().unique().tolist()
                if_vert = group['If_vert'].dropna().unique().tolist()
                if_mamm = group['If_mamm'].dropna().unique().tolist()
                data_size = group['FileSize'].sum()

                # Check for consistency
                if (len(bioproject) > 1 or 
                    len(taxid) > 1 or
                    len(taxon) > 1 or
                    len(lineage) > 1 or 
                    len(if_chor) > 1 or 
                    len(if_vert) > 1 or 
                    len(if_mamm) > 1 ):
                    logging.info(f"Multiple info involved in BioSample group {name}")

                new_row = {
                    'BioSample': name,
                    'Sex': '',
                    'Tissue': '',
                    'ScientificName': '/'.join(map(str, taxon)),
                    'TaxID': '/'.join(map(str, taxid)),
                    'Lineage': '/'.join(map(str, lineage)),
                    'If_chor': '/'.join(map(str, if_chor)),
                    'If_vert': '/'.join(map(str, if_vert)),
                    'If_mamm': '/'.join(map(str, if_mamm)),
                    'RunCount': row_count,
                    'DataSize_GB': data_size/1024/1024/1024,
                    'BioProject': '/'.join(map(str, bioproject))
                }
                new_rows.append(new_row)

            new_df = pd.DataFrame(new_rows)
            return new_df

        def add_sampleinfo(df):
            for index, row in df.iterrows():
                biosample = row['BioSample']
                logging.info(f"Fetching data for {biosample}")
                url = f"https://ngdc.cncb.ac.cn/biosample/browse/{biosample}"
                try:
                    response = requests.get(url)
                    response.raise_for_status()
                    soup = BeautifulSoup(response.text, 'html.parser')
                    sex = soup.find('th', string=lambda text: text.strip() == 'Sex')
                    tissue = soup.find('th', string=lambda text: text.strip() == 'Tissue')
                    if sex:
                        sex_td = sex.find_next('td')
                        df.at[index, 'Sex'] = sex_td.text.strip() if sex_td else ''
                    if tissue:
                        tissue_td = tissue.find_next('td')
                        df.at[index, 'Tissue'] = tissue_td.text.strip() if tissue_td else ''
                except requests.RequestException as e:
                    logging.warning(f"Error fetching data from {url}: {e}")
                    continue
            return df

        def add_ref_info(df):
            df_info = pd.read_csv(self.refinfo, sep='\t', header=1, index_col=None, low_memory=False)
            # Keep only the taxon id that is in our GSA data
            taxid = df['TaxID'].unique().tolist()
            taxid = [int(x) for x in taxid]
            df_info = df_info[df_info['taxid'].isin(taxid)]
            # Only keep the latest refseq
            df_info = df_info[df_info['version_status'] == 'latest']
            # Do not use the partial genome
            df_info = df_info[df_info['genome_rep'] == 'Full']
            
            for taxon in taxid:
                df_taxon = df_info[df_info['taxid'] == taxon].copy()
                if df_taxon.shape[0] > 0:
                    df_taxon = df_taxon.sort_values(by='refseq_category', ascending=False)
                    if df_taxon['refseq_category'].values[0] == 'reference genome':
                        df_taxon = df_taxon.head(1)
                        df.loc[df['TaxID'] == str(taxon), 'RefSeq'] = 'Y'
                    else:
                        df_taxon['assembly_level'] = pd.Categorical(df_taxon['assembly_level'], categories=['Complete Genome', 'Chromosome', 'Scaffold', 'Contig'], ordered=True)
                        df_taxon = df_taxon.sort_values(by='assembly_level', ascending=True)
                        df_taxon = df_taxon.head(1)
                        df.loc[df['TaxID'] == str(taxon), 'RefSeq'] = 'N'
                    df.loc[df['TaxID'] == str(taxon), 'Assembly'] = str(df_taxon['#assembly_accession'].values[0])
                    df.loc[df['TaxID'] == str(taxon), 'Assembly_level'] = str(df_taxon['assembly_level'].values[0])
                    df.loc[df['TaxID'] == str(taxon), 'FTP_path'] = str(df_taxon['ftp_path'].values[0])
                    df.loc[df['TaxID'] == str(taxon), 'Genome_size'] = df_taxon['genome_size'].values[0]
                else:
                    logging.warning(f'{taxon} has no refseq, please check')
            for index, row in df.iterrows():
                try:
                    genome_size = row.get('Genome_size', None)
                    data_size_gb = row.get('DataSize_GB', None)
                    if not pd.isna(genome_size) and not pd.isna(data_size_gb):
                        depth = (float(data_size_gb) * 1024*1024*1024) / int(genome_size)
                        df.at[index, 'Depth'] = round(depth, 2)
                    else:
                        df.at[index, 'Depth'] = ''
                except Exception as e:
                    logging.warning(f"Error calculating depth for index {index}: {e}")
                    df.at[index, 'Depth'] = ''
            return df
            

        df_new = rearrange(df)
        df_new = add_sampleinfo(df_new)
        df_new = add_ref_info(df_new)
        df_new.to_csv(f'{self.output_dir}/gsa_biosample.tsv', sep='\t', index=False, header=True, float_format='%.2f')
        df = df.merge(
            df_new[['BioSample', 'Sex', 'Tissue']],
            on='BioSample',
            how='left'
        )
        df['FTP_path'] = df['FTP_path'].str.replace(
            'ftp://download.big.ac.cn', 'https://download.cncb.ac.cn', regex=False
        )
            
        df.to_csv(f'{self.output_dir}/gsa_runinfo.tsv', sep='\t', header=True, index=False, float_format='%.2f')
        logging.info('Done!')

    def fetch_data(self):
        # Placeholder for data fetching logic
        if not os.path.exists(f'{self.output_dir}/gsa_biosample.tsv') or self.overwrite:
            logging.info("Annotating taxon info and prefiltering the RunInfo data...")
            df = self.prefiltering()
            logging.info("Integrating biosample information and adding refseq assembly information...")
            self.integrate_biosample(df)
        else:
            logging.info(f"Prefiltered file already exists. Use --overwrite to re-generate.")
def main():
    parser = argparse.ArgumentParser(description="Process GSA data.")
    parser.add_argument('-i', '--input', required=True, help='Input GSA RunInfo CSV file')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    parser.add_argument('-r', '--refinfo', required=True, help='RefSeq assembly summary file')
    parser.add_argument('-f', '--overwrite', action='store_true', help='Overwrite existing files')
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    if not os.path.exists(args.output):
        os.makedirs(args.output)

    gsa_data_fetch(
        input_file=args.input,
        output_dir=args.output,
        refinfo=args.refinfo,
        overwrite=args.overwrite
    )
if __name__ == "__main__":
    main()
