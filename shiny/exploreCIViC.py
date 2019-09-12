import csv

if __name__ == '__main__':
	drugs = set()
	with open('nightly-ClinicalEvidenceSummaries.tsv') as f:
		reader = csv.DictReader(f,delimiter='\t')
		for row in reader:
			if row['drugs'] != '' and not ',' in row['drugs']:
				drugs.add(row['drugs'].lower())

	drugs = sorted(list(drugs))

	with open('nightly-ClinicalEvidenceSummaries.tsv') as f:
		reader = csv.DictReader(f,delimiter='\t')
		for row in reader:
			#if row['evidence_type'] == 'Predictive' and row['drugs'] == '':
			#	print(row['evidence_id'],row['gene'],row['variant'],row['evidence_type'],row['drugs'])
			#if row['evidence_type'] != 'Predictive' and row['drugs'] != '':
			#	print(row['evidence_id'],row['gene'],row['variant'],row['evidence_type'],row['drugs'])
			if row['evidence_type'] == 'Predictive' and row['drugs'] != '':
				evidenceStatement = row['evidence_statement'].lower()
				foundDrugs = [ d for d in drugs if d in evidenceStatement ]
				if not (foundDrugs == [] or row['drugs'].lower() in foundDrugs):
			
			#and not row['drugs'].lower() in row['evidence_statement'].lower():
					print(row['evidence_id'],row['gene'],row['variant'],row['evidence_type'],row['drugs'],foundDrugs)

