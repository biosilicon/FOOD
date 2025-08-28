#!/bin/bash
# We may need to check those files had not loaded by NCBI but on EBI, especially those with pod5 format.
awk -F, -v OFS="\t" '$4 == 0' sra_runinfo.txt > spots0.txt
ids=$(cut -d, -f 1 spots0.txt)
for acc in $ids; do
    curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${acc}&result=read_run&fields=run_accession,base_count,study_accession,sample_accession,scientific_name,tax_id,fastq_ftp,sra_ftp,submitted_format,submitted_bytes,submitted_ftp,submitted_md5&format=tsv" | tail -n +2 >> spots0.tsv
    sleep 0.5
done
