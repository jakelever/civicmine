#!/bin/bash
set -ex

less CosmicFusionExport.tsv.gz | cut -f 12 -d $'\t' | perl -pe 's/{[^}]*}:[^{]*\_/\t/' | cut -f 1 -d '{' | sed -e '/^\s*$/d' | sort -u > cosmic_fusions.tsv

less CosmicMutantExport.tsv.gz | cut -f 1,19,30 -d $'\t' | sed -e 's/\tp\./\t/g' | grep -F "Confirmed somatic variant" | sort -u | cut -f 1,2 -d $'\t' > cosmic_proteinMutations.tsv
