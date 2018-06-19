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

def getNormalizedTerm(text,externalID,IDToTerm):
	normalizedTerms = [ IDToTerm[eid] for eid in externalID.split(';') ]
	normalizedTerms = sorted(list(set(normalizedTerms)))

	normalizedTermsLower = [ st.lower() for st in normalizedTerms ]
	textLower = text.lower()

	if textLower in normalizedTermsLower:
		index = normalizedTermsLower.index(textLower)
		normalizedTerm = normalizedTerms[index]
	else:
		normalizedTerm = ";".join(normalizedTerms)
	return normalizedTerm

def normalizeMIRName(externalID):
	assert externalID.startswith('mirna|'), "Unexpected ID: %s" % externalID
	normalizedName = externalID[4:]

	search = re.search('mirna\|\D*(?P<id>\d+[A-Za-z]*)',externalID)
	if search:
		mirID = search.groupdict()['id']
		if not mirID is None:
			normalizedName = "miR-%s" % mirID

	return normalizedName

headers = None
def applyFinalFilter(row):
	global headers

	# Filter out incorrect output with some rules
	if headers is None:
		with open('header.tsv') as f:
			headers = f.read().strip().split('\t')
	assert len(row) == len(headers), "Number of columns in output data (%d) doesn't  match with header count (%d)" % (len(row),len(headers))

	row = { h:v for h,v in zip(headers,row) }

	# Check for the number of semicolons (suggesting a list)
	if row['sentence'].count(';') > 5:
		return False

	if row['section'] == 'back':
		return False

	# Filter some erroneous variant associations for very common substitutions
	expectedVariantAssociations = {'R399Q': 'XRCC1', 'R194W': 'XRCC1', 'T241M': 'XRCC3', 'V600E': 'BRAF', 'T790M': 'EGFR', 'L858R': 'EGFR'}
	if row['variant_normalized'] == 'substitution':
		substitution = row['variant_id'].split('|')[1]
		if substitution in expectedVariantAssociations and not row['gene_normalized'] == expectedVariantAssociations[substitution]:
			return False

	return True

def civicmine(sentenceFile,modelFilenames,filterTerms,wordlistPickle,genes,cancerTypes,drugs,variants,outData):
	print("%s : start" % now())

	models = {}
	assert isinstance(modelFilenames,list)
	for modelFilename in modelFilenames:
		with open(modelFilename,'rb') as f:
			models[modelFilename] = pickle.load(f)

	IDToTerm = {}
	HugoToEntrez = defaultdict(lambda : "N/A")
	with codecs.open(genes,'r','utf-8') as f:
		for line in f:
			geneid,singleterm,_,entrez_geneid = line.strip().split('\t')
			IDToTerm[geneid] = singleterm
			HugoToEntrez[geneid] = entrez_geneid

	with codecs.open(cancerTypes,'r','utf-8') as f:
		for line in f:
			cancerid,singleterm,_ = line.strip().split('\t')
			IDToTerm[cancerid] = singleterm

	with codecs.open(drugs,'r','utf-8') as f:
		for line in f:
			drugid,singleterm,_ = line.strip().split('\t')
			IDToTerm[drugid] = singleterm

	with codecs.open(variants,'r','utf-8') as f:
		for line in f:
			variantid,singleterm,_ = line.strip().split('\t')
			IDToTerm[variantid] = singleterm

	with codecs.open(filterTerms,'r','utf-8') as f:
		filterTerms = [ line.strip().lower() for line in f ]

	with open(wordlistPickle,'rb') as f:
		termLookup = pickle.load(f)
	
	# Truncate the output file
	with codecs.open(outData,'w','utf-8') as outF:
		pass

	timers = Counter()

	print("%s : loading..." % now())
	with open(sentenceFile) as f:
		sentenceData = json.load(f)

	corpus = kindred.Corpus()
	for sentence in sentenceData:
		metadata = dict(sentence)
		del metadata["sentence"]
		doc = kindred.Document(sentence["sentence"],metadata=metadata)
		corpus.addDocument(doc)

	print("%s : loaded..." % now())
	startTime = time.time()
	parser = kindred.Parser()
	parser.parse(corpus)
	timers['parser'] += time.time() - startTime
	print("%s : parsed" % now())

	startTime = time.time()
	ner = kindred.EntityRecognizer(lookup=termLookup,detectVariants=True,detectFusionGenes=True,detectMicroRNA=True,acronymDetectionForAmbiguity=True,mergeTerms=True,removePathways=True)
	ner.annotate(corpus)
	timers['ner'] += time.time() - startTime
	print("%s : ner" % now())

	with codecs.open(outData,'a','utf-8') as outF:
		startTime = time.time()
		for modelname,model in models.items():
			model.predict(corpus)
		timers['predicted'] += time.time() - startTime

		print("%s : predicted" % now())

		startTime = time.time()

		for doc in corpus.documents:
			# Skip if no relations are found
			if len(doc.relations) == 0:
				continue

			# Skip if there isn't an associated PMID
			if not doc.metadata["pmid"]:
				continue

			entity_to_sentence = {}
			for sentence in doc.sentences:
				for entity,tokenIndices in sentence.entityAnnotations:
					entity_to_sentence[entity] = sentence

			geneID2Variant = defaultdict(list)
			for relation in doc.relations:
				# We're only dealing with Variants in this loop
				if relation.relationType != 'AssociatedVariant':
					continue

				typeToEntity = {}
				for entity in relation.entities:
					typeToEntity[entity.entityType] = entity

				geneID = typeToEntity['gene'].entityID
			
				v = typeToEntity['variant']
				prob = relation.probability
				#variant = (typeToEntity['variant'].externalID,typeToEntity['variant'].text,)
				geneID2Variant[geneID].append((v,prob))

			for relation in doc.relations:
				# IgnoreVariant as we deal with them seperately (above)
				if relation.relationType == 'AssociatedVariant':
					continue

				#print(relation)

				sentence = entity_to_sentence[relation.entities[0]]
				sentenceTextLower = sentence.text.lower()

				hasFilterTerm = any( filterTerm in sentenceTextLower for filterTerm in filterTerms )
				if not hasFilterTerm:
					continue
				#words = [ t.word for t in sentence.tokens ]
				#text = " ".join(words)

				sentenceStart = sentence.tokens[0].startPos

				relType = relation.relationType
				entityData = {'cancer':['' for _ in range(5)], 'drug':['' for _ in range(5)], 'gene':['' for _ in range(5)] }
				geneID = None
				for entity in relation.entities:
					if entity.entityType == 'gene':
						assert geneID is None, 'Relation should only contain a single gene'
						geneID = entity.entityID

					if entity.externalID.startswith('combo'):
						externalIDsplit = entity.externalID.split('|')
						normalizedTerms = [ getNormalizedTerm("",st.replace('&',';'),IDToTerm) for st in externalIDsplit[1:] ]
						normalizedTerm = "|".join(normalizedTerms)
					elif entity.externalID.startswith('mirna|'):
						normalizedTerm = normalizeMIRName(entity.externalID)
					else:
						normalizedTerm = getNormalizedTerm(entity.text,entity.externalID,IDToTerm)

					assert len(entity.position) == 1, "Expecting entities that are contigious and have only one start and end position within the text"
					startPos,endPos = entity.position[0]

					tmp = []
					tmp.append(entity.externalID)
					if entity.entityType == 'gene':
						tmp.append(HugoToEntrez[entity.externalID])

					tmp.append(entity.text)
					tmp.append(normalizedTerm)
					tmp.append(startPos - sentenceStart)
					tmp.append(endPos - sentenceStart)

					entityData[entity.entityType] = tmp

				associatedVariants = []
				if geneID in geneID2Variant:
					for entity,prob in list(geneID2Variant[geneID]):
						startPos,endPos = entity.position[0]
						if entity.externalID.startswith('substitution|'):
							normalizedTerm = 'substitution'
						else:
							normalizedTerm = getNormalizedTerm(entity.text,entity.externalID,IDToTerm)
						tmp = []
						tmp.append(entity.externalID)
						tmp.append(entity.text)
						tmp.append(normalizedTerm)
						tmp.append(startPos - sentenceStart)
						tmp.append(endPos - sentenceStart)
						tmp.append(prob)
						associatedVariants.append(tmp)
				else:
					blank = ['' for _ in range(6)]
					associatedVariants.append(blank)
					

				for associatedVariant in associatedVariants:
					m = doc.metadata
					if not 'subsection' in m:
						m['subsection'] = None

					prob = relation.probability
					combinedEntityData = entityData['cancer'] + entityData['gene'] + entityData['drug'] + associatedVariant
					outData = [m['pmid'],m['title'],m['journal'],m['year'],m['month'],m['day'],m['section'],m['subsection'],relType,prob] + combinedEntityData + [sentence.text]
					if applyFinalFilter(outData):
						outLine = "\t".join(map(str,outData))
						outF.write(outLine+"\n")

		timers['output'] += time.time() - startTime

		print("%s : output" % now())

	sys.stdout.flush()

	print("%s : done" % now())
	
	for section,sectiontime in timers.items():
		print("%s\t%f" % (section,sectiontime))
	print("%s\t%f" % ("applyModelsToSentences total", sum(timers.values())))


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Finds relations in Pubmed file')
	parser.add_argument('--sentenceFile',required=True,help='BioC XML file to use')
	parser.add_argument('--models',required=True)
	parser.add_argument('--filterTerms',required=True)
	parser.add_argument('--wordlistPickle',required=True)
	parser.add_argument('--genes',required=True)
	parser.add_argument('--cancerTypes',required=True)
	parser.add_argument('--drugs',required=True)
	parser.add_argument('--variants',required=True)
	parser.add_argument('--outData',required=True)

	args = parser.parse_args()

	civicmine(args.sentenceFile,args.models.split(','),args.filterTerms,args.wordlistPickle,args.genes,args.cancerTypes,args.drugs,args.variants,args.outData)
