# CIViCmine

<p>
<a href="https://travis-ci.org/jakelever/civicmine">
  <img src="https://travis-ci.org/jakelever/civicmine.svg?branch=master" />
</a>
<a href="http://bionlp.bcgsc.ca/civicmine/">
  <img src="https://img.shields.io/badge/data-viewer-blue.svg" />
</a>
</p>

This is a text mining project to assist curation in the [Clinical Interpretation of Variants in Cancer (CIViC) database](https://civicdb.org). CIViC catalogs information on diagnostic, predictive, predisposing and prognostic variants in cancer. This project aims to text mine this data from PubMed and Pubmed Central Open Access subset

In brief: it finds sentences that describe genes, variants, cancer types and optionally drugs. We then annotated many of them to create training data. This is then used to build a classifier using the [Kindred package](https://github.com/jakelever/kindred). This is then applied to all of PubMed and PubMed Central Open Access to extract sentences and structured information. This is then collated for download and viewing through a [web viewer](http://bionlp.bcgsc.ca/civicmine/).

## Dependencies

The project requires [Kindred](https://github.com/jakelever/kindred) and [PubRunner](https://github.com/jakelever/pubrunner) for execution. Kindred manages the parsing and classification. PubRunner manages the corpora download, format conversions, cluster execution and results upload to Zenodo. These can both be installed through [pip](https://pypi.org/). All other dependencies are managed by PubRunner.

## Inputs

The project uses wordlists from the [BioWordlists](https://github.com/jakelever/biowordlists) project for cancers, genes, drugs, variants and conflicting terms. The corpora used are PubMed abstracts and full-text papers from the PubMed Central Open Access Subset.

We have annotated 1500 sentences that comprise the civicmine_corpus which are used for training a Kindred classifier.

## Outputs

The output of CIViCmine is outlined below:

## Running It

A full run of CIViCmine will take some time as the large corpora need to be downloaded and parsing takes a long time. A cluster is strongly recommended for this task. PubRunner uses Snakemake for execution and can therefore be used on most clusters.

A test run (which is what the Travis-CI test does) can be executed with:
```
pubrunner --test .
```

The full run can be executed with:
```
pubrunner .
```

## Contributors



## Citing

A paper is in the works.
