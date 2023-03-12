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

def filterCorpus(corpus,filterTerms):
	filtered = kindred.Corpus()
	for doc in corpus.documents:
		termsFound = any( ft in doc.text.lower() for ft in filterTerms )
		if termsFound:
			filtered.addDocument(doc)
	return filtered

dashCharacters = ["-", "\u00ad", "\u2010", "\u2011", "\u2012", "\u2013", "\u2014", "\u2043", "\u2053"]
def fixDashes(text):
	if any (dc in text for dc in dashCharacters):
		for dc in dashCharacters:
			text = text.replace(dc,'-')
	return text

def tidyWhitespace(text):
	return re.sub(r'\s+', ' ', text)

def cleanCorpus(corpus):
	for doc in corpus.documents:
		if doc.text:
			doc.text = tidyWhitespace(fixDashes(doc.text))
		if doc.metadata['title']:
			doc.metadata['title'] = tidyWhitespace(fixDashes(doc.metadata['title']))

# Deal with table data stored in tab-delimited form
def splitTabbedCorpus(corpus):
	new_corpus = kindred.Corpus()
	for doc in corpus.documents:
		for block in doc.text.split('\t'):
			block = block.strip()
			if block:
				new_doc = kindred.Document(block)
				new_doc.metadata = doc.metadata
				new_corpus.addDocument(new_doc)

	return new_corpus

def parseAndFindEntities(biocFile,filterTermsFile,wordlistPickle,variantStopwordsFile,outSentencesFilename,verbose=False):
	if verbose:
		print("%s : start" % now())

	with open(wordlistPickle,'rb') as f:
		termLookup = pickle.load(f)

	filterTerms = None
	if filterTermsFile:
		with open(filterTermsFile,'r') as f:
			filterTerms = [ line.strip().lower() for line in f ]

	with open(variantStopwordsFile) as f:
		variantStopwords = [ line.strip() for line in f ]

	timers = Counter()

	outSentences = []

	currentID = None
	duplicateCheck = set()

	if verbose:
		print("%s : processing..." % now())
	parser = kindred.Parser(model='en_core_sci_sm')
	ner = kindred.EntityRecognizer(lookup=termLookup,detectFusionGenes=True,detectMicroRNA=True,acronymDetectionForAmbiguity=True,mergeTerms=True,detectVariants=True,variantStopwords=variantStopwords)
	for corpusno,corpus in enumerate(kindred.iterLoad('biocxml',biocFile)):
		if filterTerms:
			startTime = time.time()
			corpus = filterCorpus(corpus,filterTerms)
			timers['filter'] += time.time() - startTime

		corpus = splitTabbedCorpus(corpus)

		cleanCorpus(corpus)

		startTime = time.time()
		parser.parse(corpus)
		timers['parser'] += time.time() - startTime
		#print("%s : parsed" % now())

		# Filter extremely long sentences
		for doc in corpus.documents:
			doc.sentences = [ s for s in doc.sentences if len(s.text) < 500 ]

		startTime = time.time()
		ner.annotate(corpus)
		timers['ner'] += time.time() - startTime
		#print("%s : ner" % now())

		startTime = time.time()

		for doc in corpus.documents:

			# Reset the duplicate check set for each new PMID
			if doc.metadata['id'] != currentID:
				currentID = doc.metadata['id']
				duplicateCheck = set()

			for sentence in doc.sentences:
				sentenceTextLower = sentence.text.lower()

				# Trim extremely long sentences
				if len(sentenceTextLower) > 1000:
					continue

				if filterTerms:
					containsFilterTerm = any( ft in sentenceTextLower for ft in filterTerms)
					if not containsFilterTerm:
						continue

				entityTypesInSentence = set([ entity.entityType for entity,tokenIndices in sentence.entityAnnotations ])
				foundCancer = 'cancer' in entityTypesInSentence
				foundGene = 'gene' in entityTypesInSentence
				foundVariant = 'variant' in entityTypesInSentence

				if foundCancer and foundGene and foundVariant:
					sentenceText = sentence.text.strip(string.whitespace + ',')

					if not sentenceText in duplicateCheck:
						tmpData = dict(doc.metadata)
						tmpData['sentence'] = sentenceText
						outSentences.append(tmpData)
						duplicateCheck.add(sentenceText)

		timers['entitiesAdded'] += time.time() - startTime

		#print("%s : entities added" % now())
		sys.stdout.flush()

	with open(outSentencesFilename,'w') as f:
		json.dump(outSentences,f,indent=2)

	if verbose:
		print("%s : done" % now())

		for section,sectiontime in timers.items():
			print("%s\t%f" % (section,sectiontime))
		print("%s\t%f" % ("parseAndFindEntities total", sum(timers.values())))

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Finds relations in Pubmed file')
	parser.add_argument('--biocFile',required=True,help='BioC XML file to use')
	parser.add_argument('--filterTerms',required=False,default=None,type=str,help='File with terms to filter sentences')
	parser.add_argument('--wordlistPickle',required=True,type=str,help='Pickled wordlist of gene/cancer/drug terms')
	parser.add_argument('--variantStopwords',required=True,type=str,help='File of variants to skip')
	parser.add_argument('--outSentencesFilename',required=True,type=str,help='Output file')
	parser.add_argument('--verbose', action='store_true', help='Whether to print out information about run')

	args = parser.parse_args()

	parseAndFindEntities(args.biocFile,args.filterTerms,args.wordlistPickle,args.variantStopwords,args.outSentencesFilename,verbose=args.verbose)

