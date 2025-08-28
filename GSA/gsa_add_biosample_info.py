#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Add Biosample information to runinfo file.
'''
import sys
import pandas as pd



def main():
    runinfo_file = sys.argv[1]
    biosampleinfo_file = sys.argv[2]
    output_file = sys.argv[3]
    df_runinfo = pd.read_csv(runinfo_file, sep=',', header=0, index_col=None)
    df_biosampleinfo = pd.read_csv(biosampleinfo_file, sep='\t', header=0, index_col=None)
    df_runinfo = df_runinfo.merge(
        df_biosampleinfo[['BioSample', 'Sex', 'Tissue']],
        on='BioSample',
        how='left'
    )
            
    df_runinfo.to_csv(output_file, sep='\t', header=True, index=False)
    print('Done!')



if __name__ == "__main__":
    main()