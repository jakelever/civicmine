import argparse
from collections import defaultdict

if __name__ == '__main__':
	parser = argparse.ArgumentParser('Make CancerGenomeInterpreter data comparable to CIViCmine')
	parser.add_argument('--cgiCancerAcronyms',required=True,type=str,help='CGI Cancer Acronyms file')
	parser.add_argument('--cgiCancerGenes',required=True,type=str,help='Combined cancer_genes_upon_mutations_or_CNAs.tsv and cancer_genes_upon_trans.tsv')
	parser.add_argument('--cgiBiomarkers',required=True,type=str,help='CGI Biomakers file')
	parser.add_argument('--cancers',required=True,type=str,help='Cancer terms file')
	parser.add_argument('--drugs',required=True,type=str,help='Drugs terms file')
	parser.add_argument('--outPredictive',required=True,type=str,help='Output file for predictive associations to compare')
	parser.add_argument('--outPredisposing',required=True,type=str,help='Output file for predisposing associations to compare')
	args = parser.parse_args()

	cgiCancerAcronyms = {}
	with open(args.cgiCancerAcronyms) as f:
		header = f.readline()
		for line in f:
			acronym,term = line.strip('\n').split('\t')
			cgiCancerAcronyms[acronym] = term

	cancerLookup = defaultdict(list)
	with open(args.cancers) as inF:
		for line in inF:
			id,single_term,synonyms = line.strip('\n').split('\t')
			for synonym in synonyms.split('|'):
				cancerLookup[synonym.lower()].append(id)
	cancerLookup = dict(cancerLookup)

	drugLookup = defaultdict(list)
	with open(args.drugs) as inF:
		for line in inF:
			id,single_term,synonyms = line.strip('\n').split('\t')
			for synonym in synonyms.split('|'):
				drugLookup[synonym.lower()].append(id)
	drugLookup = dict(drugLookup)

	cancerLookup['cancer'] = ['*']
	cancerLookup['solid tumors'] = ['*']
	cancerLookup['hematological malignancy'] = ['DOID:2531']
	cancerLookup['malignant rhaboid tumor'] = ['DOID:3672']

	cgiCancerAcronyms['DT-PR'] = 'gastrointestinal system cancer'


	# gene    alteration      cancer_acronym  source
	with open(args.cgiCancerGenes) as inF, open(args.outPredisposing,'w') as outF:
		headers = inF.readline().strip('\n').split('\t')
		for lineno,line in enumerate(inF):
			row = { h:v for h,v in zip(headers,line.strip('\n').split('\t')) }
			if not 'predisposing' in row['source'].lower():
				continue

			gene = row['gene']

			# ANother header
			if gene == 'gene':
				continue

			alteration = row['alteration']

			cancer_acronym = row['cancer_acronym']
			if not cancer_acronym in cgiCancerAcronyms:
				print("Predisposing Cancer Acronym Warning: %s" % cancer_acronym)
				continue

			cancer = cgiCancerAcronyms[cancer_acronym]

			cancers = cancer.replace('predisposing','').strip('-– ')

			for cancer in cancers.split('/'):

				if not cancer in cancerLookup and cancer+' cancer' in cancerLookup:
					cancer = cancer+' cancer'
				if not cancer in cancerLookup and cancer+' tumor' in cancerLookup:
					cancer = cancer+' tumor'

				if not cancer in cancerLookup:
					print("Predisposing Cancer Warning: %s" % cancer)
					continue

				cancer_ids = cancerLookup[cancer]

				for cancer_id in cancer_ids:
					outData = [gene, alteration, cancer_id, cancer]
					outF.write("\t".join(outData) + "\n")
			
		

	cancerLookup['any cancer type'] = ['*']
	cancerLookup['all tumors'] = ['*']
	cancerLookup['all solid tumors'] = ['*']
	cancerLookup['b-lymphoblastic leukemia/lymphoma'] = ['DOID:0060592','DOID:7061']
	cancerLookup['anaplastic thyroid cancer'] = ['DOID:1781']
	cancerLookup['erdheim-chester disease'] = ['DOID:4329']
	cancerLookup['peritoneal serous carcinoma'] = ['DOID:4901']
	cancerLookup['esophagogastric cancer'] = ['DOID:5041','DOID:10534']

	cancerLookup['chronic eosinophilic leukemia, nos'] = ['DOID:0080367']
	cancerLookup['myelodysplastic/myeloproliferative neoplasms'] = ['DOID:4972']
	cancerLookup['dedifferentiated liposarcoma'] = ['DOID:3382']
	cancerLookup['uterine serous carcinoma/uterine papillary serous carcinoma'] = ['DOID:5750','DOID:5747']
	cancerLookup['histiocytosis'] = ['?']
	cancerLookup['mastocytosis'] = ['DOID:350']
	cancerLookup['low-grade serous ovarian cancer'] = ['DOID:0050933']
	cancerLookup['skin cancer, non-melanoma'] = ['DOID:4159']
	cancerLookup['ewing sarcoma of soft tissue'] = ['DOID:3369']
	cancerLookup['myelodysplastic syndromes'] = ['DOID:0050908']

	cancerLookup['adrenal adenoma'] = ['DOID:656']
	cancerLookup['anaplastic oligodendroglioma'] = ['DOID:3181']
	cancerLookup['billiary tract'] = ['DOID:4607']
	cancerLookup['cervix squamous cell'] = ['DOID:3744']
	cancerLookup['endometrium'] = ['DOID:1380']
	cancerLookup['eosinophilic chronic leukemia'] = ['DOID:0080367']
	cancerLookup['erdheim-chester histiocytosis'] = ['DOID:4329']
	cancerLookup['female germ cell tumor'] = ['DOID:2994']
	cancerLookup['gastrointestinal stromal'] = ['DOID:9253']
	cancerLookup['giant cell astrocytoma'] = ['DOID:5077']
	cancerLookup['head an neck'] = ['DOID:11934']
	cancerLookup['head an neck squamous'] = ['DOID:5520']
	cancerLookup['hematologic malignancies'] = ['DOID:2531']
	cancerLookup['hyper eosinophilic advanced snydrome'] = ['DOID:999']
	cancerLookup['inflammatory myofibroblastic'] = ['DOID:0050905']
	cancerLookup['lagerhans cell histiocytosis'] = ['DOID:2571']
	cancerLookup['lung squamous cell'] = ['DOID:3907']
	cancerLookup['lymphangioleiomyomatosis'] = ['DOID:3319']
	cancerLookup['malignant peripheral nerve sheat tumor'] = ['DOID:5940']
	cancerLookup['myelodisplasic proliferative syndrome'] = ['DOID:4972']
	cancerLookup['myelodisplasic syndrome'] = ['DOID:0050908']
	cancerLookup['pediatric glioma'] = ['DOID:6383']
	cancerLookup['renal angiomyolipoma'] = ['DOID:3314']
	cancerLookup['salivary glands'] = ['DOID:8850']
	cancerLookup['sarcoma'] = ['DOID:1115']
	cancerLookup['schwannoma'] = ['DOID:3196']
	cancerLookup['solid tumors'] = ['*']
	cancerLookup['systemic mastocytosis'] = ['DOID:349']
	cancerLookup['urinary tract carcinoma'] = ['?']


	#drugLookup['brigatinib'] = ['WIKIDATA:NOID']
	#drugLookup['lorlatinib'] = ['Q27285820']
	#drugLookup['encorafenib'] = ['Q15409405']
	#drugLookup['binimetinib'] = ['Q19903515']
	#drugLookup['niraparib'] = ['Q25326660']
	#drugLookup['dacomitinib'] = ['Q17130597']
	#drugLookup['neratinib'] = ['Q6995920']
	#drugLookup['erdafitinib'] = ['Q27077213']
	#drugLookup['gilteritinib'] = ['Q27077802']
	#drugLookup['midostaurin'] = ['Q6842945']
	#drugLookup['high dose chemotherapy'] = ['CUSTOM0002']
	#drugLookup['ivosidenib'] = ['Q27895417']
	#drugLookup['enasidenib'] = ['Q27077182']
	#drugLookup['larotrectinib'] = ['Q27081513']
	#drugLookup['alpelisib'] = ['Q27074391']
	#drugLookup['talazoparib'] = ['Q25100990']
	#drugLookup['abemaciclib'] = ['Q23901483']
	#drugLookup['poziotinib'] = ['Q27088426']
	#drugLookup['carboplatin-taxol regimen'] = ['Q415588','Q423762']
	#drugLookup['asciminib'] = ['Q27074535']
	#drugLookup['azd5363'] = []
	#drugLookup['azd9496'] = ['Q27461938']
	#drugLookup['azd4547'] = []
	#drugLookup['bgj398'] = ['Q27075200']
	#drugLookup['debio1347'] = ['Q27285013']
	#drugLookup['crenolanib'] = ['Q5184160']
	#drugLookup['tazemetostat'] = ['Q27088941']
	#drugLookup['quizartinib'] = ['Q7272714']
	#drugLookup['tipifarnib'] = ['Q7808830']
	#drugLookup['avapritinib'] = []
	#drugLookup['milademetan tosylate'] = []
	#drugLookup['ro5045337'] = ['Q27287118']
	#drugLookup['capmatinib'] = ['Q27075685']
	#drugLookup['ribociclib'] = ['Q27088552']
	#drugLookup['iodine i 131-6-beta-iodomethyl-19-norcholesterol'] = []
	#drugLookup['selumetinib'] = ['Q7448840']
	#drugLookup['entrectinib'] = ['Q25323953']
	#drugLookup['buparlisib'] = ['Q25100534']
	#drugLookup['serabelisib'] = ['Q27078071']
	#drugLookup['copanlisib'] = ['Q19903876']
	#drugLookup['gdc-0077'] = []
	#drugLookup['taselisib'] = ['Q27088940']
	#drugLookup['loxo-292'] = []
	#drugLookup['blu-667'] = []
	#drugLookup['plx8394'] = ['Q27088419']
	#drugLookup['cemiplimab'] = ['Q55606818']
	#drugLookup['tk216'] = []
	#drugLookup['gsk2636771'] = ['Q27077889']
	#drugLookup['azd8186'] = ['Q27074798']
	#drugLookup['h3b-8800'] = []

	# Isoform RefSeq  Entrez Gene ID  Hugo Symbol     Alteration      Protein Change  Cancer Type     Level   Drugs(s)        PMIDs for drug  Abstracts for drug
	with open(args.cgiBiomarkers) as inF, open(args.outPredictive,'w') as outF:
		headers = inF.readline().strip('\n').split('\t')
		for lineno,line in enumerate(inF):
			row = { h:v for h,v in zip(headers,line.strip('\n').split('\t')) }
			gene_name = row['Gene']
			alteration = row['Alteration']

			if ';' in alteration:
				alteration = 'combo'
			elif alteration.endswith(':amp'):
				alteration = 'amplification'
			elif '__' in alteration:
			 	alteration = 'fusion'
			else:
			 	assert ':' in alteration, ': not found in %s' % alteration
			 	alteration = alteration.split(':')[1].strip()
			 	if alteration == '.':
			 		alteration = 'mutation'
			 	elif alteration == 'del':
			 		alteration = 'deletion'
			 	elif alteration == 'over':
			 		alteration = 'overexpression'
			 	elif ',' in alteration:
			 		alteration = alteration.split(',')[0]

			cancers = row['Primary Tumor type full name'].lower()
			drugs = row['Drug full name'].lower()

			drugs = drugs.split('(')[0].strip()
			drugs = drugs.replace('novel ','')
			drugs = drugs.replace('liposomal ','')

			drugs = drugs.replace('entrictinib','entrectinib')
			drugs = drugs.replace('flourouracil','fluorouracil')
			drugs = drugs.replace('fluvestrant','fulvestrant')
			drugs = drugs.replace('tensirolimus','temsirolimus')
			drugs = drugs.replace('mytomycin c','mitomycin c')
			drugs = drugs.replace('platinum agent','platinum-based chemotherapy')
			drugs = drugs.replace('chemotherapys','chemotherapy')

			#print(lineno, cancer, drug)

			for drug in drugs.split('+'):
				drug = drug.strip()
				for cancer in cancers.split(';'):
					cancer = cancer.strip()

					if not cancer in cancerLookup and cancer+' cancer' in cancerLookup:
						cancer = cancer+' cancer'

					if not drug in drugLookup:
						print("Drug Warning: %s" % drug)
						drug_ids = ['?']
						#continue
					else:
						drug_ids = drugLookup[drug]
						
					if not cancer in cancerLookup:
						print("Cancer Warning: %s" % drug)
						#continue
						cancer_ids = ['?']
					else:
						cancer_ids = cancerLookup[cancer]


					for cancer_id in cancer_ids:
						for drug_id in drug_ids:
							outData = [gene_name,alteration,cancer_id,cancer,drug_id,drug]
							outF.write("\t".join(outData) + "\n")
					#assert cancer in cancerLookup, ""
					#if cancer in cancerLookup and 
					#cancer_id = cancerLookup[cancer]
					#drug_id = drugLookup[drug]
					#print('CANCER',cancer,cancer in cancerLookup)
					#print('DRUG',drug,drug in drugLookup)

					#if not drug in drugLookup:
					#	print(drug)
					#if not cancer in cancerLookup:
					#	print("cancerLookup['%s'] = ['']" % cancer)



