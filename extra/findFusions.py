import sys
import itertools
import kindred
import pickle
import argparse
import codecs
import time
import re
import string
from collections import defaultdict,Counter
import json

def now():
	return time.strftime("%Y-%m-%d %H:%M:%S")

def findFusions(biocFile,genesFile,wordlistPickle,outFile):
	print("%s : start" % now())

	with open(wordlistPickle,'rb') as f:
		termLookup = pickle.load(f)

	hugo2Name = {}
	with open(genesFile) as f:
		for line in f:
			hugo_gene_id,gene_name,synoyms,entrez_gene_id = line.strip('\n').split('\t')
			hugo2Name[hugo_gene_id] = gene_name

	print("%s : processing..." % now())
	parser = kindred.Parser(model='en_core_sci_sm')
	ner = kindred.EntityRecognizer(lookup=termLookup,detectFusionGenes=True,detectMicroRNA=True,acronymDetectionForAmbiguity=True,mergeTerms=True,detectVariants=True)
	with open(outFile,'w') as outF:
		for corpusno,corpus in enumerate(kindred.iterLoad('biocxml',biocFile)):
			parser.parse(corpus)
			ner.annotate(corpus)

			for doc in corpus.documents:
				pmid = ''
				if 'pmid' in doc.metadata:
					pmid = doc.metadata['pmid']
				#fusions = [ e for  e in doc.entities if e.entityType == 'Gene' ]
				for e in doc.entities:
					if e.entityType == 'gene' and e.externalID.startswith('combo|'):
						gene_ids = e.externalID.split('|')[1:]

						if len(gene_ids) != 2:
							continue

						if any('&' in gene_id for gene_id in gene_ids):
							continue

						for gene_id in gene_ids:
							assert gene_id in hugo2Name, 'Unable to find HUGO gene name for ID: %s (text=%s)' % (gene_id, e.text)

						gene_names = [ hugo2Name[gene_id] for gene_id in gene_ids ]

						assert len(gene_names) == 2

						outData = [pmid, e.text] + gene_ids + gene_names
						outF.write("\t".join(outData) + "\n")
						#print(outData)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Finds relations in Pubmed file')
	parser.add_argument('--biocFile',required=True,help='BioC XML file to use')
	parser.add_argument('--genes',required=True,help='Genes file')
	parser.add_argument('--wordlistPickle',required=True)
	parser.add_argument('--outFile',required=True)

	args = parser.parse_args()

	findFusions(args.biocFile,args.genes,args.wordlistPickle,args.outFile)

