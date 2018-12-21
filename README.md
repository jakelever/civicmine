# CIViCmine

<p>
<a href="https://travis-ci.org/jakelever/civicmine">
  <img src="https://travis-ci.org/jakelever/civicmine.svg?branch=master" />
</a>
<a href="http://bionlp.bcgsc.ca/civicmine/">
   <img src="https://img.shields.io/badge/data-viewer-9e42f4.svg" />
</a>
<a href="https://doi.org/10.5281/zenodo.1472826">
   <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.1472826.svg" />
</a>
<a href="https://doi.org/10.1101/500686">
   <img src="https://img.shields.io/badge/bioRxiv-preprint-67baea.svg" />
</a>
</p>

This is a text mining project to assist curation in the [Clinical Interpretation of Variants in Cancer (CIViC) database](https://civicdb.org). CIViC catalogs information on diagnostic, predictive, predisposing and prognostic variants in cancer. This project aims to text mine this data from PubMed and Pubmed Central Open Access subset

In brief: it finds sentences that describe genes, variants, cancer types and optionally drugs. We then annotated many of them to create training data. This is then used to build a classifier using the [Kindred package](https://github.com/jakelever/kindred). This is then applied to all of PubMed and PubMed Central Open Access to extract sentences and structured information. This is then collated for download and viewing through a [web viewer](http://bionlp.bcgsc.ca/civicmine/).

## System Requirements

This is a Python3 project which has been tested on Centos 6/7 but should work on other Linux operating systems and MacOS. An individual process of this can be run on a laptop or desktop computer. But in order to process all of the literature (PubMed, etc), this should really be run on a cluster or server-like machine. A cluster that uses Slurm or the SunGrid engine (SGE) are supported. Each node needs only 4 GBs on RAM.

This project relies on text mining using [Kindred](https://github.com/jakelever/kindred) and resource management with [PubRunner](https://github.com/jakelever/pubrunner). These can be installed through pip.

## Installation Guide

You can clone this repo using Git or download the [ZIP file](https://github.com/jakelever/civicmine/archive/master.zip) of it.

```
git clone https://github.com/jakelever/civicmine.git
```

The dependencies can be installed with the command below. Remember to install the English language model for Spacy.

```
pip install kindred pubrunner
python -m spacy download en
```

Installation should take a maximum of 15 minutes (mostly due to the Spacy and language models installation).

## Instructions for Use

A full run of CIViCmine will take some time as the large corpora need to be downloaded and parsing takes a long time. A cluster is strongly recommended for this task. PubRunner uses Snakemake for execution and can therefore be used on most clusters.

A test run (which is what the Travis-CI test does) can be executed with:
```
pubrunner --test .
```

The full run can be executed with:
```
pubrunner .
```

## Inputs

The project uses wordlists from the [BioWordlists](https://github.com/jakelever/biowordlists) project for cancers, genes, drugs, variants and conflicting terms. The corpora used are PubMed abstracts and full-text papers from the PubMed Central Open Access Subset. While processing the full-text articles, subsection headers (e.g. Results) are also extracted to make it easier to locate where statements are made. The list of possible headers extracted can be found [here](https://github.com/jakelever/pubrunner/blob/master/subsectionHeaders.md). The actual extraction is managed by PubRunner.

We have annotated 1500 sentences that comprise the civicmine_corpus which are used for training a Kindred classifier.

## Outputs

The three output files of CIViCmine is outlined below. You likely want **civicmine\_collated.tsv** if you just want the list of cancer biomarkers. If you want the supporting sentences, look at **civicmine\_sentences.tsv**. You can use the *matching\_id* column to connect the two files. If you want to dig further and are okay with a higher false positive rate, look at **civicmine\_unfiltered.tsv**.

**civicmine\_collated.tsv:** This contains the cancer biomarkers with citation counts supporting them. It contains the normalized cancer and gene names along with IDs for HUGO, Entrez Gene and the Disease Ontology.

**civicmine\_sentences.tsv:** This contains the supporting sentences for the cancer biomarker in the collated file. Each row is a single supporting sentence for one cancer biomarker. This file contains information on the source publication (e.g. journal, publication date, etc), the actual sentence and the cancer biomarker extracted.

**civicmine\_unfiltered.tsv:** This is the raw output of the applyModelsToSentences.py script across all of PubMed, Pubmed Central Open Access and PubMed Central Author Manuscript Collection. It contains every predicted relation with a prediction score above 0.5. So this may contain many false positives. Each row contain information on the publication (e.g. journal, publication date, etc) along with the sentence and the specific cancer biomarker extracted (with HUGO, Entrez Gene and Disease Ontology IDs). This file is further processed to create the other two.

## Shiny App

The code in [shiny/](https://github.com/jakelever/civicmine/tree/master/shiny) is the Shiny code used for the [web viewer](http://bionlp.bcgsc.ca/civicmine/). If it is helpful, please use the code for your own projects. The list of dependencies is found at the top of the [app.R](https://github.com/jakelever/civicmine/blob/master/shiny/app.R) file.

## Paper

The code to generate all the figures and text for the paper can be found in [paper/](https://github.com/jakelever/civicmine/tree/master/paper). This may be useful for generating an up-to-date version of the plots for a newer version of CIViCmine.

## Citing

A paper is in the works.

## Issues

If you encounter any problems, please [file an issue](https://github.com/jakelever/civicmine/issues) along with a detailed description.

