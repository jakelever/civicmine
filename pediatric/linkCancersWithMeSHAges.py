import kindred
import argparse
import re
from collections import defaultdict,Counter
import json
from tqdm import tqdm


def main():
	parser = argparse.ArgumentParser(description='Pull MeSH age groups associated with cancer types')
	parser.add_argument('--biocFile',required=True,help='BioC XML file to use')
	parser.add_argument('--cancers',required=True,type=str,help='Cancer list')
	parser.add_argument('--prettyPrint',action='store_true',help='Whether to pretty print the output JSON file')
	parser.add_argument('--outFile',required=True,type=str,help='Output JSON file')

	args = parser.parse_args()

	term_lookup = defaultdict(set)
	cancer_id_to_name = {}
	with open(args.cancers) as f:
		for line in f:
			cancer_id, cancer_name, synonyms = line.strip('\n').split('\t')

			cancer_id_to_name[cancer_id] = cancer_name
			for synonym in synonyms.split('|'):
				term_lookup[synonym].add(('cancer',cancer_id))

	parser = kindred.Parser(model='en_core_sci_sm')
	ner = kindred.EntityRecognizer(lookup=term_lookup)

	age_terms = set(['Pediatrics','Infant','Infant, Newborn','Child','Child, Preschool','Adolescent','Young Adult','Adult','Middle Aged','Aged, 80 and over','Frail Elderly'])

	mesh_ages = {}
	cancer_mentions = defaultdict(set)

	for corpus in kindred.iterLoad('biocxml',args.biocFile):
		corpus.documents = [ doc for doc in corpus.documents if doc.metadata['pmid'] and doc.metadata['pmid'] != 'None' ]

		for doc in corpus.documents:

			pmid = doc.metadata['pmid']
			found_age_terms = []
			if doc.metadata['meshHeadings']:
				mesh_headings = [ heading.strip() for heading in doc.metadata['meshHeadings'].split('\t') ]
				mesh_headings = [ heading.split('~') for heading in mesh_headings if heading ]

				for descriptor_qualifiers in mesh_headings:
					descriptor = descriptor_qualifiers[0].split('|')
					assert len(descriptor) == 4, "Expected four pipe-delimited columns. Got: %s" % descriptor_qualifiers[0]
					_, mesh_id, isMajorYN, name = descriptor
					if name in age_terms:
						found_age_terms.append(name)

				found_age_terms = sorted(set(found_age_terms))

			mesh_ages[pmid] = found_age_terms

		corpus.documents = [ doc for doc in corpus.documents if len(mesh_ages[doc.metadata['pmid']]) > 0 ]

		parser.parse(corpus)

		ner.annotate(corpus)

		for doc in corpus.documents:

			pmid = doc.metadata['pmid']
			found_cancer_ids = sorted(set( e.externalID for e in doc.entities if e.entityType == 'cancer' ))

			if not found_cancer_ids:
				continue

			found_cancer_names = sorted( cancer_id_to_name[cancer_id] for cancer_id in found_cancer_ids )

			cancer_mentions[pmid].update(found_cancer_names)

	output = {}
	for pmid in cancer_mentions:
		output[pmid] = { 'ages':sorted(mesh_ages[pmid]), 'cancer_names':sorted(cancer_mentions[pmid]) }

	with open(args.outFile,'w') as outF:
		if args.prettyPrint:
			json.dump(output,outF,indent=2,sort_keys=True)
		else:
			json.dump(output,outF)

	print("Done.")

if __name__ == '__main__':
	main()

