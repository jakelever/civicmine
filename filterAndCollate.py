import sys
import argparse
from collections import Counter, defaultdict
import hashlib

if __name__ == '__main__':
	assert sys.version_info[0] == 3, "Must use Python 3"

	parser = argparse.ArgumentParser(description='Filter Civicmine for more conservative predictions, collate the results and prepare the sentences for viewing')
	parser.add_argument('--inUnfiltered',required=True,type=str,help='Civicmine TSV to filter')
	parser.add_argument('--outSentences',required=True,type=str,help='Output filtered Civicmine')
	parser.add_argument('--outCollated',required=True,type=str,help='Output filtered Civicmine')
	args = parser.parse_args()

	thresholds = {}
	thresholds['AssociatedVariant'] = 0.7
	thresholds['Diagnostic'] = 0.63
	thresholds['Predictive' ] = 0.93
	thresholds['Predisposing'] = 0.86
	thresholds['Prognostic'] = 0.65

	collated = defaultdict(set)
	collatedMatchingID = {}
	
	sentenceOrdering = set()

	sentenceData = {}
	geneLocs = defaultdict(set)
	cancerLocs = defaultdict(set)
	drugLocs = defaultdict(set)
	variantLocs = defaultdict(set)

	sentenceKeyFields = 'pmid,title,journal,journal_short,year,month,day,section,subsection,evidencetype,cancer_id,cancer_normalized,gene_hugo_id,gene_entrez_id,gene_normalized,drug_id,drug_normalized,variant_id,variant_normalized,sentence'

	readCount,writeCount = 0,0
	with open(args.inUnfiltered,'r') as inF, open(args.outSentences,'w') as outF:
		headers = inF.readline().strip('\n').split('\t')

		#columnsToRemove = [ i for i,h for enumerate(headers) if '_start' in h or '_end' in h ]

		#newheaders = ['matchingid'] + [ h for i,h in enumerate(headers) if not i in columnsToRemove ] + ['formatted_sentence']
		newheaders = ['matching_id'] + headers + ['formatted_sentence']
		outF.write("\t".join(newheaders) + "\n")

		variantCols = [ h for i,h in enumerate(headers) if h.startswith('variant_') ]

		for line in inF:
			r = { h:v for h,v in zip(headers,line.strip('\n').split('\t')) }

			evidencetype_prob = float(r['evidencetype_prob'])
			evidencetype = r['evidencetype']

			readCount += 1
			if evidencetype_prob < thresholds[evidencetype]:
				continue

			pmid = r['pmid']
			gene_hugo_id = r['gene_hugo_id']
			gene_entrez_id = r['gene_entrez_id']
			gene_normalized = r['gene_normalized']
			cancer_id = r['cancer_id']
			cancer_normalized = r['cancer_normalized']
			drug_id = r['drug_id']
			drug_normalized = r['drug_normalized']
			variant_id = r['variant_id']
			variant_normalized = r['variant_normalized']


			ambigious = any( ';' in v for v in [gene_hugo_id,gene_entrez_id,gene_normalized,cancer_id,cancer_normalized,drug_id,drug_normalized,variant_id,variant_normalized] )

			if ambigious:
				# Ambigious terms so skip
				continue

			if pmid == 'None':
				# Uncitable sentence so skip
				continue

			# Blank the variant data if the associated probability is below our threshold
			if r['variant_prob'] != '':
				variant_prob = float(r['variant_prob'])
				if variant_prob <= thresholds['AssociatedVariant']:
					for h in variantCols:
						r[h] = ''
					variant_id = r['variant_id']
					variant_normalized = r['variant_normalized']

			variant_group = variant_normalized
		
			variant_withsub = variant_normalized
			if variant_normalized == 'substitution':
				substitution = variant_id.split('|')[1]
				variant_withsub = '%s (substitution)' % substitution

			r['journal_short'] = r['journal']
			if len(r['journal_short']) > 50:
				r['journal_short'] = r['journal_short'][:50] + '...'

			collatedKey = (evidencetype,gene_hugo_id,gene_entrez_id,gene_normalized,cancer_id,cancer_normalized,drug_id,drug_normalized,variant_group,variant_withsub)
			
			matchingID = hashlib.md5("|".join(list(collatedKey)).encode('utf-8')).hexdigest()
			collatedMatchingID[collatedKey] = matchingID

			collated[collatedKey].add(pmid)

			allHeaders = 'pmid,title,journal,journal_short,year,month,day,section,subsection,evidencetype,evidencetype_prob,cancer_id,cancer_text,cancer_normalized,cancer_start,cancer_end,gene_hugo_id,gene_entrez_id,gene_text,gene_normalized,gene_start,gene_end,drug_id,drug_text,drug_normalized,drug_start,drug_end,variant_id,variant_text,variant_normalized,variant_start,variant_end,variant_prob,sentence'

			sentenceKey = tuple( [ r[k] for k in sentenceKeyFields.split(',') ] )
			sentenceData[sentenceKey] = matchingID
			sentenceOrdering.add( (r['year'],sentenceKey) )

			geneLocs[sentenceKey].add((int(r['gene_start']),int(r['gene_end'])))
			cancerLocs[sentenceKey].add((int(r['cancer_start']),int(r['cancer_end'])))
			if r['drug_start'] and r['drug_end']:
				drugLocs[sentenceKey].add((int(r['drug_start']),int(r['drug_end'])))
			if r['variant_start'] and r['variant_end']:
				variantLocs[sentenceKey].add((int(r['variant_start']),int(r['variant_end'])))



	with open(args.outCollated,'w') as outF:
		headers = [ 'matching_id', 'evidencetype', 'gene_hugo_id','gene_entrez_id','gene_normalized','cancer_id','cancer_normalized','drug_id','drug_normalized','variant_group','variant_withsub', 'citation_count' ]
		outF.write("\t".join(headers) + "\n")

		collatedCounts = [ (len(pmids),key) for key,pmids in collated.items() ]
		collatedCounts = sorted(collatedCounts,reverse=True)
		for citationCount,key in collatedCounts:
			assert len(key) == (len(headers)-2)
			matchingID = collatedMatchingID[key]
			outF.write("%s\t%s\t%d\n" % (matchingID,"\t".join(list(key)),citationCount))



	with open(args.outSentences,'w') as outF:
		headers = 'matching_id,pmid,title,journal,journal_short,year,month,day,section,subsection,evidencetype,cancer_id,cancer_normalized,gene_hugo_id,gene_entrez_id,gene_normalized,drug_id,drug_normalized,variant_id,variant_normalized,sentence,formatted_sentence'
		headerCount = len(headers.split(','))
		outF.write(headers.replace(',','\t') + '\n')

		sentenceOrdering = sorted(list(sentenceOrdering),reverse=True)

		for _,sentenceKey in sentenceOrdering:
			matchingID = sentenceData[sentenceKey]
			s = { h:v for h,v in zip(sentenceKeyFields.split(','), list(sentenceKey)) }

			sentence = s['sentence']

			charByChar = list(sentence)
			for locDict in [geneLocs,cancerLocs,drugLocs,variantLocs]:
				for start,end in locDict[sentenceKey]:
					charByChar[start] = '<b>' + charByChar[start]
					charByChar[end-1] += '</b>'

			formattedSentence = "".join(charByChar)

			outData = [ matchingID] + list(sentenceKey) + [ formattedSentence ]
			assert len(outData) == headerCount

			outLine = "\t".join(outData)
			outF.write(outLine + "\n")

	print("%d of %d rows written out to %s" % (len(sentenceData),readCount,args.outSentences))
	print("%d biomarkers written out to %s" % (len(collatedCounts),args.outCollated))

