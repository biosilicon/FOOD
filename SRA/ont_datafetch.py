#!/usr/bin/env python3
# ont_datafetch.py
# Description: Data fetching utilities for the FindMe pipeline.
import shutil
import os
import time
import logging
import subprocess
import pandas as pd
import xml.etree.ElementTree as ET

class ONTDataFetch:
    def __init__(self, output_dir:str, overwrite:bool=False):
        self.output_dir = output_dir
        self.overwrite = overwrite
        os.makedirs(output_dir, exist_ok=True)
        os.chdir(output_dir)
        self.check_environment()
    
    def check_environment(self) -> None:
        flag = False
        for command in ['esearch', 'efetch', 'taxonkit', 'tee', 'esummary']:
            if not shutil.which(command):
                logging.error(f"Required command '{command}' is not found in the system PATH.")
                flag = True
        if flag:
            raise EnvironmentError("One or more required commands are not available in the system PATH.")
        else:
            logging.info("All required commands are available. Checking taxonomy database...")
        required_files = ['names.dmp', 'nodes.dmp', 'delnodes.dmp', 'merged.dmp']
        for required_file in required_files:
            if not os.path.exists(f'~/.taxonkit/{required_file}'):
                logging.info('Missing taxonkit database, downloading...')
                os.makedirs('~/.taxonkit', exist_ok=True)
                cmd = 'cd ~/.taxonkit && wget https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz && tar -zxvf taxdump.tar.gz && cd -'
                result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
                if result.returncode == 0:
                    logging.info("Taxonomy database installed successfully.")
                    break
                else:
                    raise RuntimeError(f"Failed to install taxonomy database, please install it mannually following https://bioinf.shenwei.me/taxonkit/#dataset")
        logging.info("Environment check passed.")

    def fetchsra(self) -> None:
        # Check if output file exists and handle overwrite
        if os.path.exists(f'{self.output_dir}/sra_runinfo.txt') and not self.overwrite:
            logging.info("SRA runinfo file already exists and overwrite is set to False. Skipping fetch.")
            return
        # Ensure output directory exists
        try:
            # Fetch SRA data using esearch and efetch
            # We need ONT data, but not RNA data, and exclude certain organisms to accelerate the search
            # If you are adding new runs, please change PDAT to a more recent date
            cmd = f'esearch -db sra -query "Oxford Nanopore [Platform] NOT rna data [Filter] NOT coronavirus [Organism] NOT metagenome [Organism] NOT Escherichia coli [Organism] AND 2000/01/01:3000/12/31[PDAT]" | efetch -format runinfo > {self.output_dir}/sra_runinfo.txt'
            # self.step_info.update({'fetchsra': True if os.system(cmd)==0 else False}) # Use subprocess for better error handling
            result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
            if result.returncode == 0:
                logging.info("SRA data fetched successfully.")
            else:
                raise RuntimeError(f"Error fetching SRA data: {result.stderr}")
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Failed to execute SRA fetch command: {e.stderr}") from e
        except Exception as e:
            raise RuntimeError(f"Unexpected error during SRA fetch: {str(e)}") from e
    
    def filtersra(self) -> pd.DataFrame:
        '''
        Filter SRA RunInfo CSV file based on specific criteria and collect biosample info.

        Functions:
        pre_filter(df): Filter the dataframe with specific criteria
        taxon_filter(df): Annotate and filter the dataframe with lineage information
        extract_biosample_info(biosample, xml_data): Extract the sex and tissue information from the XML data of a biosample
        get_biosample_info(df): Collect biosample info of each run including sexual and tissue

        Returns:
        dataframe: Filtered dataframe with additional biosample information
        '''
        def pre_filter(df:pd.DataFrame) -> pd.DataFrame:
            # Filter with data availability
            df = df[df["bases"] > 0]
            # Filter with library type
            df = df[df["LibraryStrategy"] == "WGS"]
            df = df[df["LibrarySelection"].str.lower().isin(["random", "unspecified", "other"])]
            df = df[df["LibrarySource"] == "GENOMIC"]
            # Filter with disease and tumor
            df = df[~df["Disease"].notna()]
            df = df[df["Tumor"] == "no"]
            # Keep only the columns of interest
            head = ['Run','ReleaseDate','LoadDate','spots','bases','avgLength','size_MB',
                    'download_path','Experiment','LibraryStrategy','LibrarySelection',
                    'Platform','Model','SRAStudy','BioProject','Study_Pubmed_id','ProjectID',
                    'Sample','BioSample','SampleType','TaxID','ScientificName',
                    'SampleName','Submission']
            df = df[head]
            return df
        
        def taxon_filter(df:pd.DataFrame) -> pd.DataFrame:
            # Extract taxon id and get the lineage
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
            # Keep rows where Lineage contains 'Metazoa'
            df = df[df["Lineage"].str.contains("Metazoa", case=False, na=False)].copy()
            df.loc[:, "If_chor"] = df["Lineage"].apply(lambda x: "Y" if "Chordata" in x else "N")
            df.loc[:, "If_vert"] = df["Lineage"].apply(lambda x: "Y" if "Vertebrata" in x else "N")
            df.loc[:, "If_mamm"] = df["Lineage"].apply(lambda x: "Y" if "Mammalia" in x else "N")
            cmd = "rm -f taxid.tmp lineage.tmp"
            subprocess.run(cmd, shell=True)
            return df
        
        def extract_biosample_info(biosample:str, xml_data:str) -> tuple:
            '''
            Deal with the XML data of a biosample returned from esummary
            and extract the sex and tissue information.

            Parameters:
            xml_data (str): XML data of a biosample returned from esummary

            Returns:
            sex, tissue (str) if found, otherwise None
            '''
            try:
                root = ET.fromstring(xml_data)
            except ET.ParseError:
                logging.warning(f"Failed to parse XML data of biosample {biosample}")
                return None, None
            sex = None
            tissue = None
            biosample = root.find(".//BioSample")
            if (biosample is not None):
                attributes = biosample.find("Attributes")
                if (attributes is not None):
                    for attribute in attributes.findall("Attribute"):
                        if (attribute.get("attribute_name") == "sex" or
                            attribute.get("harmonized_name") == "sex" or
                            attribute.get("display_name") == "sex"):
                            sex = attribute.text
                        elif (attribute.get("attribute_name") == "tissue" or
                            attribute.get("harmonized_name") == "tissue" or
                            attribute.get("display_name") == "tissue"):
                            tissue = attribute.text
            return sex, tissue

        def get_biosample_info(df):
            biosample_list = df["BioSample"].tolist()
            biosample_list = list(set(biosample_list))
            biosample_list.sort()
            sample_sex = {}
            sample_tissue = {}
            count = 1
            for biosample in biosample_list:
                logging.info(f'{biosample}........... ({count}/{len(biosample_list)})')
                start_time = time.time()
                cmd = f"esummary -db biosample -id {biosample}"
                biosample_info = subprocess.check_output(cmd, shell=True).decode("utf-8")
                sex, tissue = extract_biosample_info(biosample, biosample_info)
                sample_sex[biosample] = sex
                sample_tissue[biosample] = tissue
                end_time = time.time()
                # sleep for 0.5 seconds to avoid being blocked by NCBI
                if (end_time - start_time) < 0.5:
                    time.sleep(0.5 - (end_time - start_time))
                count += 1
            df['Sex'] = df['BioSample'].map(sample_sex)
            df['Tissue'] = df['BioSample'].map(sample_tissue)
            # Mannually remove some unwanted runs
            df = df[~df['Tissue'].str.contains('culture|lympho', case=False, na=False)]
            df = df[~df['BioProject'].isin(['PRJEB66174', 'PRJNA740254', 'PRJNA1146547'])]
            return df
        df = pd.read_csv(f'{self.output_dir}/sra_runinfo.txt', sep=",", header=0, index_col=None, low_memory=False)
        df = pre_filter(df)
        df = taxon_filter(df)
        df = get_biosample_info(df)
        df.to_csv(f'{self.output_dir}/sra_runinfo_filtered.tsv', sep='\t', index=False, header=True, float_format='%.0f')
        logging.info(f"Filtered data saved to {self.output_dir}/sra_runinfo_filtered.tsv")
        return df

    def addrawformat(self, df:pd.DataFrame) -> pd.DataFrame:
        '''
        Get all raw format of each run from xml file
        '''
        def esearch(runid:str) -> str:
            cmd = "esearch -db sra -query " + runid + " | efetch -format xml"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            return result.stdout
        
        def parse_xml(run_info: str):
            # parse xml file to get all raw format of the run
            root = ET.fromstring(run_info)
            semantic_names = []
            for sra_file in root.findall(".//SRAFile[@supertype='Original']"):
                semantic_name = sra_file.get('semantic_name')
                if semantic_name:
                    semantic_names.append(semantic_name)
            return '/'.join(list(set(semantic_names)))
        
        logging.info("Adding raw data format information to each run...")
        run_list = df['Run'].tolist()
        for run_id in run_list:
            try:
                run_info = esearch(run_id)
                raw_format = parse_xml(run_info)
                df.loc[df['Run'] == run_id, 'Raw_Data_Format'] = raw_format
            except:
                logging.warning(f"Warning: No raw format found for {run_id}")
                df.loc[df['Run'] == run_id, 'Raw_Data_Format'] = None
        df.to_csv(f'{self.output_dir}/sra_runinfo_filtered_with_rawformat.tsv', sep='\t', index=False, header=True, float_format='%.0f')
        logging.info(f"Data with raw format saved to {self.output_dir}/sra_runinfo_filtered_with_rawformat.tsv")
        return df
        
    def assignref(self) -> pd.DataFrame:
        genbank_url = 'https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt'
        # Download GenBank assembly summary
        if not os.path.exists(f'{self.output_dir}/assembly_summary_genbank.txt') or self.overwrite:
            try:
                df = pd.read_csv(genbank_url, sep='\t', header=1, index_col=None, low_memory=False)
                df.to_csv(f'{self.output_dir}/assembly_summary_genbank.txt', sep='\t', header=True, index=False)
                logging.info("GenBank assembly summary downloaded successfully.")
            except Exception as e:
                raise RuntimeError(f"Failed to download GenBank assembly summary: {str(e)}") from e
        else:
            logging.info("GenBank assembly summary already exists and overwrite is set to False. Loading existing file.")
            df = pd.read_csv(f'{self.output_dir}/assembly_summary_genbank.txt', sep='\t', header=0, index_col=None, low_memory=False)
        return df
    
    def integrate(self, df:pd.DataFrame) -> pd.DataFrame:
        '''
        Integrate the information by Biosample ID
        '''
        def rearrange(df):
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
                biosample_sex = group['Sex'].dropna().unique().tolist()
                biosample_tissue = group['Tissue'].dropna().unique().tolist()
                data_size = group['size_MB'].sum()
                base_count = group['bases'].sum()
                raw_data_format = group['Raw_Data_Format'].dropna().unique().tolist()

                # Check for consistency
                if (len(bioproject) > 1 or 
                    len(taxid) > 1 or
                    len(taxon) > 1 or
                    len(lineage) > 1 or 
                    len(if_chor) > 1 or 
                    len(if_vert) > 1 or 
                    len(if_mamm) > 1 or 
                    len(biosample_sex) > 1 or 
                    len(biosample_tissue) > 1):
                    logging.warning(f"Multiple info involved in BioSample group {name}")

                new_row = {
                    'BioSample': name,
                    'Sex': '/'.join(map(str, biosample_sex)),
                    'Tissue': '/'.join(map(str, biosample_tissue)),
                    'ScientificName': '/'.join(map(str, taxon)),
                    'TaxID': '/'.join(map(str, taxid)),
                    'Lineage': '/'.join(map(str, lineage)),
                    'If_chor': '/'.join(map(str, if_chor)),
                    'If_vert': '/'.join(map(str, if_vert)),
                    'If_mamm': '/'.join(map(str, if_mamm)),
                    'RunCount': row_count,
                    'DataSize_GB': data_size/1024,
                    'BaseCount_GB': base_count/1024/1024/1024,
                    'BioProject': '/'.join(map(str, bioproject)),
                    'Raw_Data_Format': '/'.join(map(str, raw_data_format))
                }
                new_rows.append(new_row)

            new_df = pd.DataFrame(new_rows)
            return new_df

        def addrefinfo(df:pd.DataFrame, df_info:pd.DataFrame) -> pd.DataFrame:
            # Keep only the taxon id that is in our SRA data
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
                    df.loc[df['TaxID'] == str(taxon), 'GenomeSize'] = str(df_taxon['genome_size'].values[0])
                else:
                    logging.warning(f'{taxon} has no refseq, please check')
            # Calculate the depth of each sample by using the genome size and base count
            df['GenomeSize'] = pd.to_numeric(df['GenomeSize'], errors='coerce')
            df['BaseCount_GB'] = pd.to_numeric(df['BaseCount_GB'], errors='coerce')
            df['Depth'] = df['BaseCount_GB']*1024*1024*1024 / df['GenomeSize']
            df['GenomeSize'] = df['GenomeSize'].astype('Int64')
            return df

        df = rearrange(df)
        df_info = self.assignref()
        df = addrefinfo(df, df_info)
        df.to_csv(f'{self.output_dir}/sra_runinfo_integrated.tsv', sep='\t', index=False, header=True, float_format='%.2f')
        logging.info(f"Integrated data saved to {self.output_dir}/sra_runinfo_integrated.tsv")

    def process_data(self):
        self.fetchsra()
        df = self.filtersra()
        df = self.addrawformat(df)
        self.integrate(df)
        logging.info("Data processing completed.")
        