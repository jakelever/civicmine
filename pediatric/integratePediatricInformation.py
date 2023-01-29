import argparse
import csv
from collections import defaultdict,Counter
import xml.etree.cElementTree as etree
import io
import re
import unicodedata
import calendar
import json
import hashlib
from tqdm import tqdm
from Bio import Entrez
Entrez.email = 'jake.lever@glasgow.ac.uk'

def chunks(lst, n):
	"""Yield successive n-sized chunks from lst."""
	for i in range(0, len(lst), n):
		yield lst[i:i + n]

def getMeSHTerms(requested_pmids):
	document_mesh = defaultdict(set)
	for pmid_chunk in tqdm(list(chunks(requested_pmids, 500))):
		#print("  Progress:", len(document_mesh))
		handle = Entrez.efetch(db='pubmed', id=pmid_chunk, rettype="gb", retmode="xml")
		xml_data = handle.read().decode('utf-8')
		for event, elem in etree.iterparse(io.StringIO(xml_data), events=("start", "end", "start-ns", "end-ns")):
			if event == "end" and elem.tag == "PubmedArticle":  # MedlineCitation'):
				pmid_field = elem.find("./MedlineCitation/PMID")
				assert pmid_field is not None
				pmid = int(pmid_field.text)

				mesh_elems = elem.findall("./MedlineCitation/MeshHeadingList/MeshHeading")
				for mesh_elem in mesh_elems:
					descriptor_elem = mesh_elem.find("./DescriptorName")
					mesh_id = descriptor_elem.attrib["UI"]
					major_topic_yn = descriptor_elem.attrib["MajorTopicYN"]
					descriptor_name = descriptor_elem.text
					
					document_mesh[pmid].add(descriptor_name)

					
		
	return document_mesh

def main():
	parser = argparse.ArgumentParser(description='Take in a CIViCmine output and integrate in information about the pediatric status')
	parser.add_argument('--inSentences',required=True,type=str,help='Input sentences file')
	parser.add_argument('--cancers',required=True,type=str,help='Cancer list')
	parser.add_argument('--meshAges',required=True,type=str,help='MeSH data on different cancer types')
	parser.add_argument('--outSentences',required=True,type=str,help='Output sentences file')
	parser.add_argument('--outCollated',required=True,type=str,help='Output collated file')
	args = parser.parse_args()
	
	mesh_pediatric_groups = set(['Pediatrics','Infant','Infant, Newborn','Child','Child, Preschool','Adolescent'])
	mesh_adult_groups = set(['Adult','Aged','Middle Aged','Young Adult','Aged, 80 and over','Frail Elderly'])
	mesh_age_groups = mesh_pediatric_groups.union(mesh_adult_groups)
	
	pediatric_keywords = [ 'pediatric', 'paediatric', 'childhood', 'infantile', 'juvenile', 'teenage', 'adolescent' ]
	
	print("Loading cancer types...")
	cancer_synonyms = defaultdict(list)
	with open(args.cancers,encoding='utf-8') as f:
		for line in f:
			cancer_id,mainterm,synonyms = line.strip().split('\t')
			
			synonyms = [ s.lower().strip() for s in synonyms.split('|') ]
			synonyms = sorted(set(synonyms))
			for synonym in synonyms:
				cancer_synonyms[synonym].append( (cancer_id,mainterm) )
				
	#cancer_synonyms['ependymoma'].append( ('DOID:7497','brain ependymoma') )
	#cancer_synonyms['infratentorial ependymoma'].append( ('DOID:7497','brain ependymoma') )
	
	print("Loading MeSH data on cancer types...")
	with open(args.meshAges) as f:
		mesh_ages = json.load(f)
		
	#cancer_age_spread = {}
	pediatric_cancers, adult_cancers = set(), set()
	for name,counts in mesh_ages.items():
		pediatric_count = sum( count for age_group,count in counts.items() if age_group in mesh_pediatric_groups )
		adult_count = sum( count for age_group,count in counts.items() if age_group in mesh_adult_groups )

		total_count = pediatric_count + adult_count
		pediatric_percentage = pediatric_count / total_count

		#cancer_age_spread[name] = {'total':total_count, 'pediatric_percentage':pediatric_count / total_count}

		if total_count > 100 and pediatric_count > 0.7:
			pediatric_cancers.add(name)
		elif total_count > 100 and pediatric_count < 0.3:
			adult_cancers.add(name)

	print("Identified %d cancer types that are predominantly pediatric" % len(pediatric_cancers))
	print("Identified %d cancer types that are predominantly adult" % len(adult_cancers))
	
	print("Loading CIViCmine sentences and collated files...")
	with open(args.inSentences, encoding='utf8') as inF:
		reader = csv.DictReader(inF, delimiter="\t")
		sentences = [ row for row in reader ]
		
	pmids = sorted(set([ int(s['pmid']) for s in sentences]))
	print("Gathering MeSH terms for %d documents" % len(pmids))
	mesh_for_documents = getMeSHTerms(pmids)
			
	print("Processing sentences to highlight pediatric information...")
	all_paper_counts = defaultdict(set)
	ped_paper_counts = defaultdict(set)
			
	journal_keywords = ['pediatric','paediatric','child']
		
	collatedKeyFields = 'evidencetype,gene_hugo_id,gene_entrez_id,gene_normalized,cancer_id,cancer_normalized,drug_id,drug_normalized,variant_group,variant_withsub'.split(',')
		
	tidied_sentences = []
	collated_by_matching_id = {}
	# Do stuff
	for s in sentences:
		if s['gene_text'].lower() == 'osteosarcoma':
			continue
	
		pmid = int(s['pmid'])
				
		is_pediatric_journal = any( k in s['journal'].lower() for k in journal_keywords)
		
		specific_pediatric_cancer_mention = any( k in s['cancer_text'].lower() for k in pediatric_keywords )

		is_pediatric_cancer = s['cancer_normalized'] in pediatric_cancers
		is_adult_cancer = s['cancer_normalized'] in adult_cancers
		is_pediatric_not_adult_cancer = is_pediatric_cancer and not is_adult_cancer
		
		is_pediatric_paper = any( k in mesh_for_documents[pmid] for k in mesh_pediatric_groups )
		is_adult_paper = any( k in mesh_for_documents[pmid] for k in mesh_adult_groups )
		is_pediatric_not_adult_paper = is_pediatric_paper and not is_adult_paper
		
		is_pediatric = specific_pediatric_cancer_mention or is_pediatric_not_adult_cancer or is_pediatric_not_adult_paper
		#is_pediatric = specific_pediatric_cancer_mention # is_pediatric_not_adult_paper
		#is_pediatric = is_pediatric_not_adult_paper
		
		s['is_pediatric'] = is_pediatric
		
		
		s['variant_group'] = s['variant_normalized']
		s['variant_withsub'] = s['variant_normalized']
		if s['variant_normalized'] == 'substitution':
			substitution = s['variant_id'].split('|')[1]
			s['variant_withsub'] = '%s (substitution)' % substitution

		collatedKey = tuple( [ s[k] for k in collatedKeyFields ] )
		matching_id = hashlib.md5("|".join(list(collatedKey)).encode('utf-8')).hexdigest()
		collated_by_matching_id[matching_id] = collatedKey
		s['matching_id'] = matching_id
		
		if is_pediatric:
			ped_paper_counts[matching_id].add(pmid)
		all_paper_counts[matching_id].add(pmid)
		
		tidied_sentences.append(s)
		
		
	sentences = tidied_sentences
		
	print("Calculating paper counts for new collated file...")
	collated = []
	for matching_id,collatedKey in collated_by_matching_id.items():
		c = {}
		c['matching_id'] = matching_id
		
		c.update( { k:v for k,v in zip(collatedKeyFields,collatedKey) } )
		
		c['citation_count'] = len(all_paper_counts[matching_id])
		c['ped_citation_count'] = len(ped_paper_counts[matching_id]) if matching_id in ped_paper_counts else 0
		collated.append(c)
		
	print("Saving new sentence and collated files with pediatric data...")
	with open(args.outSentences, 'w', newline='', encoding='utf8') as outF:
		fieldnames = sentences[0].keys()
		writer = csv.DictWriter(outF, fieldnames=fieldnames, delimiter="\t")

		writer.writeheader()
		for s in sentences:
			writer.writerow(s)
			
	with open(args.outCollated, 'w', newline='', encoding='utf8') as outF:
		fieldnames = collated[0].keys()
		writer = csv.DictWriter(outF, fieldnames=fieldnames, delimiter="\t")

		writer.writeheader()
		for s in collated:
			writer.writerow(s)
			
	print("%d of %d associations were flagged as pediatric" % (len(ped_paper_counts),len(all_paper_counts)))
	print("Done.")

if __name__ == '__main__':
	main()
	
