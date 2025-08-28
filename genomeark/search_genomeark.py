#!/usr/bin/env python3
"""
This script searches the GenomeArk website and retrieves species information.

Usage:
python3 search_genomeark.py -d <metadata_dir> -o <output_file>

Arguments:
-d/--dir: Dir path of genomeark metadata(.yaml) files, can be downloaded from github: https://github.com/genomeark/genomeark-metadata/tree/main/species
-o/--output_file: Path to save the result file

"""
import re
import yaml
import argparse
import pandas as pd
import requests
from bs4 import BeautifulSoup

def search_genomeark(metadata_dir, output_file):
    def retrieve_metadata(metadata_dir, species_name):
        taxid = common_name = genome_size = biosample = sex = individual_name = ''
        try:
            with open (f"{metadata_dir}/{species_name}.yaml", 'r') as file:
                metadata = yaml.safe_load(file)
                taxid = metadata.get('species', {}).get('taxon_id', '')
                common_name = metadata.get('species', {}).get('common_name', '')
                genome_size = metadata.get('species', {}).get('genome_size', '')
                individuals = metadata.get('species', {}).get('individuals', [{}])
                if individuals == ['']:
                    individuals = [{}]
                individual_name = '/'.join([str(individual.get('short_name', '')) for individual in individuals if individual.get('short_name') is not None])
                biosample = '/'.join([str(individual.get('biosample_id', '')) for individual in individuals if individual.get('biosample_id') is not None])
                sex = '/'.join([str(individual.get('sex', '')) for individual in individuals if individual.get('sex') is not None])
        except FileNotFoundError:
            print(f"{species_name} seems have no metadata")
        return taxid, common_name, genome_size, individual_name, biosample, sex

    def check_ont(url):
        response = requests.get(url)
        if response.status_code == 200:
            soup = BeautifulSoup(response.content, 'html.parser')
            table = soup.find('table', class_='raw-data-table')
            withont = 'N'
            onts3 = ont42bp = ''
            if table:
                for row in table.find_all('tr', class_='raw-data-toprow'):
                    cells = row.find_all('td')
                    if cells[0].get('class') == ['label']:
                        label_text = cells[0].text.strip().upper()
                        

                        if 'NANOPORE' in label_text or 'ONT' in label_text:
                            withont = 'Y'
                            links = row.find_all('a')
                            onts3 = links[0]['href'] if len(links) > 0 else ''
                            ont42bp = links[1]['href'] if len(links) > 1 else ''
        return withont, onts3, ont42bp

    url = 'https://www.genomeark.org/genomeark-all/'
    response = requests.get(url)

    if response.status_code == 200:
        df = pd.DataFrame(columns=['Species', 
                                   'CommonName', 
                                   'TaxID', 
                                   'GenomeSize', 
                                   'Platform',
                                   'WithONT', 
                                   'BioSample', 
                                   'IndividualName', 
                                   'Sex', 
                                   'AssemblyInfo', 
                                   'SpeciesPage', 
                                   'ONTAccessS3', 
                                   'ONTAccess42BP'])
        soup = BeautifulSoup(response.content, 'html.parser')

        species_list = soup.find('div', class_='species-list')

        for species in species_list.find_all('div', class_='species-listing'):
            
            name = species.find('div', class_='species-list-info-name').find('i').text.strip().replace(' ', '_')
            taxid, common_name, genome_size, individual_name, biosample, sex = retrieve_metadata(metadata_dir, name)
            sub_page = species.find('a')['href']
            species_page = 'https://www.genomeark.org' + sub_page if sub_page.startswith('/') else sub_page
            info_body = species.find('div', class_='species-list-info-body').text.strip()
            reads_match = re.search(r'Reads:\s*(.*?)\s*Assembly:', info_body, re.DOTALL)
            assembly_match = re.search(r'Assembly:\s*(.*)', info_body, re.DOTALL)
            if reads_match:
                reads_info = re.sub(r'\s+', '', reads_match.group(1)).replace(':::', '/')
            else:
                reads_info = None
            withont = 'N'
            onts3 = ont42bp = ''
            if ('ONT' in reads_info.upper() or 'NANOPORE' in reads_info.upper()):
                withont, onts3, ont42bp =  check_ont(species_page)
            if assembly_match:
                assembly_info = re.sub(r'\s+', '', assembly_match.group(1)).replace(':::', '/')
            else:
                assembly_info = None
            
            df = df._append({'Species': name,
                            'CommonName': common_name,
                            'TaxID': taxid,
                            'GenomeSize': genome_size,
                            'Platform': reads_info,
                            'WithONT': withont,
                            'BioSample': biosample,
                            'IndividualName': individual_name,
                            'Sex': sex,
                            'AssemblyInfo': assembly_info,
                            'SpeciesPage': species_page,
                            'ONTAccessS3': onts3,
                            'ONTAccess42BP': ont42bp}, ignore_index=True)
        df = df.sort_values(by='WithONT', ascending=False)
        df.to_csv(output_file, sep=',', index=False, header=True)
        print(f"Result saved to {output_file}")

    else:
        print(f"Failed to retrieve the webpage. Status code: {response.status_code}")


def main():
    parser = argparse.ArgumentParser(description="Search the GenomeArk website and retrieve species information")
    parser.add_argument("-d", "--dir", help="Dir path of genomeark metadata")
    parser.add_argument("-o", "--output_file", help="Path to save the result file")
    args = parser.parse_args()

    search_genomeark(args.dir, args.output_file)  

if __name__ == "__main__":
    main()