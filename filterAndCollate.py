import argparse
import csv
import os
import hashlib
from collections import defaultdict

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Filter CIViCmine for more conservative predictions')
	parser.add_argument('--inData',required=True,type=str,help='Input directory with TSV files to be filtered')
	parser.add_argument('--outUnfiltered',required=True,type=str,help='Output unfiltered data')
	parser.add_argument('--outCollated',required=True,type=str,help='Output collated and filtered data')
	parser.add_argument('--outSentences',required=True,type=str,help='Output filtered sentences that match collated data')
	args = parser.parse_args()

	assert os.path.isdir(args.inData)

	thresholds = {}
	thresholds['AssociatedVariant'] = 0.7
	thresholds['Diagnostic'] = 0.63
	thresholds['Predictive' ] = 0.93
	thresholds['Predisposing'] = 0.86
	thresholds['Prognostic'] = 0.65

	collated = defaultdict(set)
	collatedMatchingID = {}

	collatedKeyFields = 'evidencetype,gene_hugo_id,gene_entrez_id,gene_normalized,cancer_id,cancer_normalized,drug_id,drug_normalized,variant_group,variant_withsub'

	inputFiles = sorted( [ os.path.join(args.inData,f) for f in os.listdir(args.inData) if f.endswith('.tsv') ] )

	pubmedInputFiles = [ f for f in os.listdir(args.inData) if f.startswith('pubmed') and f.endswith('.tsv') ]
	pmcInputFiles = [ f for f in os.listdir(args.inData) if f.startswith('pmc') and f.endswith('.tsv') ]

	inputFiles = sorted(pmcInputFiles,reverse=True) + sorted(pubmedInputFiles,reverse=True)

	inputFilesHeader = None

	pmidsAlreadySeen = set()

	recordCount,filteredRecordCount = 0,0
	with open(args.outUnfiltered,'w') as outUnfiltered, open(args.outSentences,'w') as outSentences:
		for inputFile in inputFiles:
			with open(os.path.join(args.inData,inputFile)) as inF:
				
				headers = inF.readline().strip('\n').split('\t')
				if inputFilesHeader is None:
					inputFilesHeader = headers

					outUnfiltered.write("matching_id\t" + "\t".join(headers) + '\n')
					outSentences.write("matching_id\t" + "\t".join(headers) + '\n')
				else:
					assert inputFilesHeader == headers, "Headers don't match expected in file %s" % inputFile

				pmidsInThisFile = set()

				for i,line in enumerate(inF):
					row = line.strip('\n').split('\t')

					assert len(row) == len(headers), "Got %d columns, expected %d in row %d, file %s" % (len(row),len(headers),i+1,inputFile)
					r = { h:v for h,v in zip(headers,row) }

					evidencetype = r['evidencetype']
					evidencetype_prob = float(r['evidencetype_prob'])
					threshold = thresholds[evidencetype]
					recordCount += 1

					pmid = r['pmid']
					if pmid in pmidsAlreadySeen:
						continue
					pmidsInThisFile.add(pmid)
					
					r['variant_group'] = r['variant_normalized']

					r['variant_withsub'] = r['variant_normalized']
					if r['variant_normalized'] == 'substitution':
						substitution = r['variant_id'].split('|')[1]
						r['variant_withsub'] = '%s (substitution)' % substitution

					collatedKey = tuple( [ r[k] for k in collatedKeyFields.split(',') ] )
					matchingID = hashlib.md5("|".join(list(collatedKey)).encode('utf-8')).hexdigest()

					outUnfiltered.write(matchingID + "\t" + "\t".join(r[h] for h in headers) + "\n")

					keepIt = evidencetype_prob > threshold and pmid != 'None'
					if keepIt:
						collated[collatedKey].add(pmid)

						# Make a field using the key data that can be used to match between tables
						collatedMatchingID[collatedKey] = matchingID

						outSentences.write(matchingID + "\t" + "\t".join(r[h] for h in headers) + "\n")
						filteredRecordCount += 1

				pmidsAlreadySeen.update(pmidsInThisFile)


	with open(args.outCollated,'w') as outF:
		headers = 'matching_id,%s,citation_count' % collatedKeyFields
		headerCount = len(headers.split(','))
		outF.write(headers.replace(',','\t') + '\n')

		collatedCounts = [ (len(pmids),key) for key,pmids in collated.items() ]
		collatedCounts = sorted(collatedCounts,reverse=True)
		for citation_count,collatedKey in collatedCounts:

			matchingID = collatedMatchingID[collatedKey]

			outData = [matchingID] + list(collatedKey) + [str(citation_count)]
			assert len(outData) == headerCount

			outLine = "\t".join(outData)
			outF.write(outLine + "\n")

	print("%d records filtered to %d sentences and collated to %d gene/cancer associations" % (recordCount, filteredRecordCount, len(collated)))
	print("Written to %s and %s" % (args.outSentences, args.outCollated))

