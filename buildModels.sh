#!/bin/bash
set -eux pipefail

if [ ! -d models ]; then
	echo "Extracting annotation data and build Kindred models"
	mkdir -p dataForModels models
	tar xzf data/training.tar.gz -C dataForModels/ --strip-components 1
	tar xzf data/testing.tar.gz -C dataForModels/ --strip-components 1
	python buildModels.py --inTrain dataForModels --outDir models/
fi

