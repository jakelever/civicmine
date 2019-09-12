#!/bin/bash
set -ex

rm -f cgi_biomarkers_latest.zip catalog_of_cancer_genes_latest.zip cgi_cancer_acronyms.tsv cgi_cancer_genes_upon_mutations_or_CNAs.tsv cgi_biomarkers_per_variant.tsv
rm -f cgi_predictive.tsv cgi_predisposing.tsv

wget -O cgi_biomarkers_latest.zip https://www.cancergenomeinterpreter.org/data/cgi_biomarkers_latest.zip
wget -O catalog_of_cancer_genes_latest.zip https://www.cancergenomeinterpreter.org/data/catalog_of_cancer_genes_latest.zip

unzip -p catalog_of_cancer_genes_latest.zip cancer_acronyms.tsv > cgi_cancer_acronyms.tsv
unzip -p catalog_of_cancer_genes_latest.zip cancer_genes_upon_mutations_or_CNAs.tsv > cgi_cancer_genes_upon_mutations_or_CNAs.tsv
unzip -p cgi_biomarkers_latest.zip cgi_biomarkers_per_variant.tsv > cgi_biomarkers_per_variant.tsv

rm -f cgi_biomarkers_latest.zip catalog_of_cancer_genes_latest.zip

python makeCancerGenomeInterpreterComparable.py --cgiBiomarkers cgi_biomarkers_per_variant.tsv --cancers terms_cancers.tsv --drugs terms_drugs.tsv --cgiCancerAcronyms cgi_cancer_acronyms.tsv --cgiCancerGenes cgi_cancer_genes_upon_mutations_or_CNAs.tsv --outPredictive cgi_predictive.tsv --outPredisposing cgi_predisposing.tsv

rm -f cgi_cancer_acronyms.tsv cgi_cancer_genes_upon_mutations_or_CNAs.tsv cgi_biomarkers_per_variant.tsv

