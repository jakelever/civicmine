#!/bin/bash

rm -f oncoKB_allActionableVariants.txt
wget -O oncoKB_allActionableVariants.txt https://oncokb.org/api/v1/utils/allActionableVariants.txt

python makeOncoKBComparable.py --oncokbActionable oncoKB_allActionableVariants.txt --cancers terms_cancers.tsv --drugs terms_drugs.tsv --outFile oncoKB_predictive.tsv

rm -f oncoKB_allActionableVariants.txt

