**Run SRA part first to install the environment**

Raw data fetching due to GSA does not have a tool like edirect:

go to https://ngdc.cncb.ac.cn/gsa/search/ search with:

((("GENOMIC"[source]) AND "WGS"[strategy]) AND "RANDOM"[selection]) NOT "sra"[fileType]

choose Platform "OXFORD_NANOPORE" on the left side of the webpage

Total Items: 1039 Experiments with 1230 Runs up to 20250218

click "send to" button -> choose file -> select RunInfo at Format -> click Create Files 
to download RunInfo.csv
and change several rows and columns because of confusing ',' using 
