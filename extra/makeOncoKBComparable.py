import argparse
from collections import defaultdict

if __name__ == '__main__':
	parser = argparse.ArgumentParser('Make OncoKB data comparable to CIViCmine')
	parser.add_argument('--oncokbActionable',required=True,type=str,help='Actionable variants OncoKB file')
	parser.add_argument('--cancers',required=True,type=str,help='Cancer terms file')
	parser.add_argument('--drugs',required=True,type=str,help='Drugs terms file')
	parser.add_argument('--outFile',required=True,type=str,help='Output file')
	args = parser.parse_args()

	cancerLookup = defaultdict(list)
	cancerList = {}
	with open(args.cancers) as inF:
		for line in inF:
			id,single_term,synonyms = line.strip('\n').split('\t')
			cancerList[id] = single_term
			for synonym in synonyms.split('|'):
				cancerLookup[synonym.lower()].append(id)
	cancerLookup = dict(cancerLookup)

	drugLookup = defaultdict(list)
	drugList = {}
	with open(args.drugs) as inF:
		for line in inF:
			id,single_term,synonyms = line.strip('\n').split('\t')
			drugList[id] = single_term
			for synonym in synonyms.split('|'):
				drugLookup[synonym.lower()].append(id)
	drugLookup = dict(drugLookup)

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

	for term,ids in cancerLookup.items():
		for id in ids:
			if not id in cancerList:
				cancerList[id] = term
	for term,ids in drugLookup.items():
		for id in ids:
			if not id in drugList:
				drugList[id] = term

	# Isoform RefSeq  Entrez Gene ID  Hugo Symbol     Alteration      Protein Change  Cancer Type     Level   Drugs(s)        PMIDs for drug  Abstracts for drug
	with open(args.oncokbActionable) as inF, open(args.outFile,'w') as outF:
		headers = inF.readline().strip('\n').split('\t')
		for lineno,line in enumerate(inF):
			row = { h:v for h,v in zip(headers,line.strip('\n').split('\t')) }
			entrez_gene_id = row['Entrez Gene ID']
			gene_name = row['Hugo Symbol']

			variant = row['Alteration']

			variant = variant.replace('Amplification','amplification')
			variant = variant.replace('Fusions','fusion')
			variant = variant.replace('Oncogenic Mutations','mutation')
			variant = variant.replace('Truncating Mutations','truncation')
			variant = variant.replace('Wildtype','wildtype')

			cancer = row['Cancer Type'].lower()
			drugs = row['Drugs(s)'].lower()

			drugs = drugs.replace('high dose chemotherapy','chemotherapy')

			print(lineno, cancer, drugs)

			#assert cancer in cancerLookup, ""
			cancer_ids = cancerLookup[cancer]

			drug_ids = ['?']
			for drug in drugs.replace('+',',').split(','):
				drug = drug.strip()
				if drug in drugLookup:
					drug_ids = drugLookup[drug.strip()]
				else:
				 	print("WARNING. Couldn't find %s" % drug)

			for cancer_id in cancer_ids:
				for drug_id in drug_ids:
					outData = [entrez_gene_id,gene_name,variant,cancer_id,cancer,drug_id,drugs]
					outF.write("\t".join(outData) + "\n")



