import argparse
import json

if __name__ == '__main__':
	parser = argparse.ArgumentParser('Find protein coding changes in dbSNP and catalog with gene info')
	parser.add_argument('--dbsnpFile',required=True,type=str,help='dbSNP file')
	parser.add_argument('--proteinMappings',required=True,type=str,help='File with mappings from protein accessions to gene IDs')
	args = parser.parse_args()

	proteinAccessionToGeneInfo = {}
	with open(args.proteinMappings) as f:
		for line in f:
			proteinAccession,geneName,geneID = line.strip('\n').split('\t')
			assert not proteinAccession in proteinAccessionToGeneInfo, 'Duplicate protein accession: %s' % proteinAccession
			proteinAccessionToGeneInfo[proteinAccession] = (geneName,geneID)

	with open(args.dbsnpFile) as f:
		for line in f:
			data = json.loads(line)
			dbsnpid = 'rs' + data['refsnp_id']

			codingChanges = set()

			alleleDiff = {}
			for p in data["primary_snapshot_data"]["placements_with_allele"]:
				for i,a in enumerate(p["alleles"]):
					hgvs = a["hgvs"]
					if hgvs.startswith('NC_'):
						deleted = a['allele']['spdi']['deleted_sequence']
						inserted = a['allele']['spdi']['inserted_sequence']
						alleleDiff[i] = (deleted != inserted)


			for p in data["primary_snapshot_data"]["placements_with_allele"]:
				for i,a in enumerate(p["alleles"]):
					if not alleleDiff[i]:
						continue

					hgvs = a["hgvs"]
					if hgvs.startswith('NP_'):
						proteinAccession = hgvs.split('.')[0]
						codingChange = hgvs.split(':')[1]

						if proteinAccession in proteinAccessionToGeneInfo and codingChange.startswith('p.') and 'spdi' in a['allele']:
							deleted = a['allele']['spdi']['deleted_sequence']
							inserted = a['allele']['spdi']['inserted_sequence']
							position = int(a['allele']['spdi']['position'])

							if len(deleted) == 1 and len(inserted) == 1:
								singleLetterCodingChange = "%s%d%s" % (deleted,position+1,inserted)
								geneName,geneID = proteinAccessionToGeneInfo[proteinAccession]

								codingChanges.add( (proteinAccession,geneID,geneName,singleLetterCodingChange) )

			if len(codingChanges) == 0:
				continue

			majorAlleleFrequency = 'Unknown'
			bestTotalCount = -1
			for a in data['primary_snapshot_data']["allele_annotations"]:
				for f in a["frequency"]:
					study_name = f['study_name']
					deleted = f['observation']['deleted_sequence']
					inserted = f['observation']['inserted_sequence']
					allele_count = f['allele_count']
					total_count = f['total_count']

					if deleted == inserted and total_count > bestTotalCount:
						bestTotalCount = total_count
						majorAlleleFrequency = allele_count / float(total_count)

			#assert majorAlleleFrequency is not None, "Error finding allele frequency for %s" % dbsnpid

			codingChanges = sorted(codingChanges)
			for proteinAccession,geneID,geneName,codingChange in codingChanges:
				if not '=' in codingChange:
					outData = [dbsnpid,proteinAccession,geneID,geneName,codingChange,str(majorAlleleFrequency)]
					print("\t".join(outData))


