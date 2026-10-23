This describes the output files for the [CIViCmine](https://github.com/jakelever/civicmine) project. These files are loaded directly by the [CIViCmine viewer](https://civicmine.ai4biomed.org). The code for this viewer is available in the CIViCmine Github repo if you want to run it independently. Each file is a gzipped tab-delimited file with a header, no comments and no quoting.

You likely want **civicmine\_collated.tsv.gz** if you just want the list of cancer biomarkers. If you want the supporting sentences, look at **civicmine\_sentences.tsv.gz**. You can use the *matching\_id* column to connect the two files. If you want to dig further and are okay with a higher false positive rate, look at **civicmine\_unfiltered.tsv.gz**.

**civicmine\_collated.tsv.gz:** This contains the cancer biomarkers with citation counts supporting them. It contains the normalized cancer and gene names along with IDs for HUGO, Entrez Gene and the Disease Ontology.

**civicmine\_sentences.tsv.gz:** This contains the supporting sentences for the cancer biomarker in the collated file. Each row is a single supporting sentence for one cancer biomarker. This file contains information on the source publication (e.g. journal, publication date, etc), the actual sentence and the cancer biomarker extracted.

**civicmine\_unfiltered.tsv.gz:** This is the raw output of the applyModelsToSentences.py script across all of PubMed and the PubMed Central articles available through the PMC Cloud Service (prepared with [BioText](https://github.com/jakelever/biotext)). It contains every predicted relation with a prediction score above 0.5. So this may contain many false positives. Each row contain information on the publication (e.g. journal, publication date, etc) along with the sentence and the specific cancer biomarker extracted (with HUGO, Entrez Gene and Disease Ontology IDs). This file is further processed to create the other two.

