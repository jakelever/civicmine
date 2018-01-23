#!/bin/bash
set -euxo pipefail

dir=$1

for f in `grep -l drug $dir/*.a2`
do
	cat $f | grep drug | cut -f 2- -d $'\t' | sed -e 's/ drug:T[0-9]*//' > toRemove
	grep -vFf toRemove $f > fixed
	
	mv fixed $f
	rm toRemove
done

