#!/usr/bin/env python3
import argparse
import os
import json
from concurrent.futures import ThreadPoolExecutor

def get_args():
    parser = argparse.ArgumentParser(description='Predicting a list of runs using Longbow')
    parser.add_argument('-i', '--input', type=str, required=True, help='Input file')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output dir')
    parser.add_argument('-f', '--overwrite', action='store_true', help='Overwrite the fastq')
    return parser.parse_args()

def pred(input_file, output_dir, overwrite=False):
    def download(run_id, output_dir):
        cmd = f'fastq-dump \
            {run_id} -N 1 -X 10000 -O {output_dir}/{run_id}'
        if not os.path.exists(f'{output_dir}/{run_id}/{run_id}.fastq') or overwrite:
            os.system(cmd)
        return True if os.path.exists(f'{output_dir}/{run_id}/{run_id}.fastq') else False

    def run_longbow(run_id, output_dir):
        cmd = f'longbow -i {output_dir}/{run_id}/{run_id}.fastq -o {output_dir}/{run_id}/{run_id}.json -t 12'
        os.system(cmd)
        return True if os.path.exists(f'{output_dir}/{run_id}/{run_id}.json') else False

    def data_collection(output_dir, run_ids, failed_run_ids):
        success_run_ids = [run_id for run_id in run_ids if run_id not in failed_run_ids]
        with open (f'{output_dir}/summary.tsv', 'w') as outf:
            outf.write('Run_ID\tFlowcell-Software-Version-Mode\tConfidenceLevel\n')
            for run_id in success_run_ids:
                with open(f'{output_dir}/{run_id}/{run_id}.json', 'r') as f:
                    data = json.load(f)
                    flowcell = data.get("Flowcell", "UNKNOWN")
                    software = data.get("Software", "UNKNOWN")
                    version = data.get("Version", "UNKNOWN")
                    mode = data.get("Mode", "UNKNOWN")
                    conf = data.get("Confidence level", "UNKNOWN")
                    if conf is None:
                        conf = '-'
                    combined_value = f"{flowcell}-{software}-{version}-{mode}"
                    outf.write(f'{run_id}\t{combined_value}\t{conf}\n')
            for run_id in failed_run_ids:
                outf.write(f'{run_id}\tFAILED\tFAILED\n')

    def processing(run_id, output_dir):
        nonlocal count
        count += 1
        print(f'Processing {run_id} in {output_dir} ...... {count}/{len(run_ids)}', flush=True)
        flag = download(run_id, output_dir)
        if flag:
            flag = run_longbow(run_id,  output_dir)
            if not flag:
                print(f'Failed to run longbow on {run_id}', flush=True)
                faild_run_ids.append(run_id)
        else:
            print(f'Failed to download {run_id}', flush=True)
            faild_run_ids.append(run_id)
        
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    with open (input_file, 'r') as f:
        run_ids = f.readlines()
    run_ids = [run_id.strip() for run_id in run_ids]
    faild_run_ids = []
    count = 0
    with ThreadPoolExecutor(max_workers=28) as executor:
        executor.map(processing, run_ids, [output_dir]*len(run_ids))

    data_collection(output_dir, run_ids, faild_run_ids)
        
def main():
    args = get_args()
    pred(args.input, args.output, args.overwrite)

if __name__ == '__main__':
    main()