#!/bin/bash
set -ex

rm -f dbsnp_index.html dbsnp_files.txt
rm -f snpsWithCodingChangings.tsv

wget -O dbsnp_index.html ftp://ftp.ncbi.nih.gov/snp/latest_release/JSON/

less dbsnp_index.html | grep -oP "ftp://[^\"]*.bz2" | grep chr | sort -u > dbsnp_files.txt

cat dbsnp_files.txt | xargs -I FILE echo "python findProteinMutationsInDBSNP.py --dbsnpFile <(curl -s --output - FILE | bzip2 -dc -) --proteinMappings proteinMappings.tsv" | bash > snpsWithCodingChangings.tsv

rm -f dbsnp_index.html dbsnp_files.txt
