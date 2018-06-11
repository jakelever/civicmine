import argparse

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Filter Civicmine for more conservative predictions')
	parser.add_argument('--inCivicmine',required=True,type=str,help='Civicmine TSV to filter')
	parser.add_argument('--outFiltered',required=True,type=str,help='Output filtered Civicmine')
	args = parser.parse_args()

	thresholds = {'Driver':0.80, 'Oncogene': 0.76, 'Tumor_Suppressor': 0.92}

	thresholds = {}
	thresholds['AssociatedVariant'] = 0.7
	thresholds['Diagnostic'] = 0.63
	thresholds['Predictive' ] = 0.93
	thresholds['Predisposing'] = 0.86
	thresholds['Prognostic'] = 0.65

	readCount,writeCount = 0,0
	with open(args.inCivicmine,'r') as inF, open(args.outFiltered,'w') as outF:
		headers = inF.readline().strip('\n').split('\t')
		outF.write("\t".join(headers) + "\n")

		variantCols = [ i for i,h in enumerate(headers) if h.startswith('variant_') ]

		seen = set()

		for line in inF:
			row = line.strip('\n').split('\t')
			rowDict = { h:v for h,v in zip(headers,row) }
			evidencetype_prob = float(rowDict['evidencetype_prob'])
			evidencetype = rowDict['evidencetype']
			readCount += 1

			if evidencetype_prob > thresholds[evidencetype]:
			
				# Blank the variant data if the associated probability is below our threshold
				if rowDict['variant_prob'] != '':
					variant_prob = float(rowDict['variant_prob'])
					if variant_prob <= thresholds['AssociatedVariant']:
						for c in variantCols:
							row[c] = ''

				rowTuple = tuple(row)
				if not rowTuple in seen:
					outF.write("\t".join(row) + "\n")
					writeCount += 1
					seen.add(rowTuple)

	print("%d of %d rows written out to %s" % (writeCount,readCount,args.outFiltered))

