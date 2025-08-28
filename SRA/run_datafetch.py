#!/usr/bin/env python3
import argparse
import logging
from ont_datafetch import ONTDataFetch

def get_args():
    parser = argparse.ArgumentParser(description="Fetch and process ONT sequencing data from NCBI SRA.",
                                     usage="python run_datafetch.py -o <output_dir> [--overwrite]")
    parser.add_argument('-o', '--output_dir', type=str, required=True, help='Output directory to save fetched data.')
    parser.add_argument('-f', '--overwrite', action='store_true', help='Overwrite existing processed data files if they exist.')

    args = parser.parse_args()
    return args

def main():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    args = get_args()
    data_fetcher = ONTDataFetch(output_dir=args.output_dir, overwrite=args.overwrite)
    data_fetcher.process_data()

if __name__ == "__main__":
    main()
