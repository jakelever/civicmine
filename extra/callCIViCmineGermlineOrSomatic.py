import argparse

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Call subtitutions in CIViCmine as likely germline or somatic using COSMIC and dbSNP')
	parser.add_argument('--dbSNP',required=True,type=str,help='Pre-processed dbSNP file')
	parser.add_argument('--COSMIC',required=True,type=str,help='Pre-processed COSMIC file')
	parser.add_argument('--civicmine',required=True,type=str,help='CIViCmine collated file')
	parser.add_argument('--outFile',required=True,type=str,help='Output file')
	args = parser.parse_args()

	dbsnp = set()
	# rs991   NP_001091970    56243   KIAA1217        A1247A  0.9999960123457775
	with open(args.dbSNP) as f:
		for line in f:
			snpid,protein_accession,gene_id,gene,variant,majorAlleleFreq = line.strip('\n').split('\t')

			if majorAlleleFreq == 'Unknown':
				continue

			majorAlleleFreq = float(majorAlleleFreq)
			if majorAlleleFreq > 0.99:
				continue # Less than 1% prevalance

			dbsnp.add ( (gene_id,variant) )
			#if len(dbsnp) > 10000:
			#	break

	print("Loaded dbSNP")

	cosmic = set()
	with open(args.COSMIC) as f:
		for line in f:
			gene,variant = line.strip('\n').split('\t')
			cosmic.add( (gene,variant) )

			#if len(cosmic) > 10000:
			#	break
	print("Loaded COSMIC")

	output = set()
	# matching_id     evidencetype    gene_hugo_id    gene_entrez_id  gene_normalized cancer_id       cancer_normalized       drug_id drug_normalized variant_group	variant_withsub  citation_count
	with open(args.civicmine) as inF:
		header = inF.readline().strip('\n').split('\t')
		for line in inF:
			row = { h:v for h,v in zip(header,line.strip('\n').split('\t')) }

			if row['variant_group'] != 'substitution':
				continue

			evidencetype = row['evidencetype']

			# Skip gene fusions / multi-genes
			gene_hugo_id = row['gene_hugo_id']
			if gene_hugo_id.startswith('combo'):
				continue

			gene_entrez_id = row['gene_entrez_id']
			gene = row['gene_normalized']
			variant = row['variant_withsub'].replace('(substitution)','').strip()

			inDBSNP = (gene_entrez_id,variant) in dbsnp
			inCOSMIC = (gene,variant) in cosmic

			outData = [ evidencetype, gene_entrez_id, gene, variant, inDBSNP, inCOSMIC ]

			output.add( tuple(outData) )

	with open(args.outFile,'w') as outF:
		for outData in sorted(output):
			outF.write("\t".join(map(str,outData)) + "\n")
