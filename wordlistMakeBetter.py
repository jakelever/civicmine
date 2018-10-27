
import argparse
import sys
import os
from collections import Counter,defaultdict

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Help identify conflicting terms in a set of wordlists')
	parser.add_argument('--knowledgebase',required=True,type=str,help='Knowledge base to use for weighting conflicts')
	args = parser.parse_args()

	termtypes = ['cancers','genes','drugs','variants']
	#termtypes = ['proteins']



	for termtype in termtypes:
		lookup = defaultdict(set)
		singleterms = {}
		synonymLookup = {}

		deletions = defaultdict(set)

		filename = "deletions_%s.tsv" % termtype
		if os.path.isfile(filename):
			with open(filename) as f:
				for line in f:
					split = line.strip('\n').split('\t')
					termid,singleterm,synonyms = split[:3]
					deletions[termid].update( [ s.lower() for s in synonyms.split('|') ] )

		allTermIDs = []
		filename = "terms_%s.tsv" % termtype
		with open(filename) as f:
			for line in f:
				split = line.strip('\n').split('\t')
				termid,singleterm,synonyms = split[:3]

				singleterms[termid] = singleterm
				synonyms = [ s.lower() for s in synonyms.split('|') ]
				synonyms = [ s for s in synonyms if not s in deletions[termid] ]

				allTermIDs.append(termid)

				synonymLookup[termid] = synonyms
				for s in synonyms:
					lookup[s].add(termid)



		conflicting = set([ term for term,termids in lookup.items() if len(termids) > 1 ])

		seen = Counter()
		with open(args.knowledgebase) as f:
		#with open('civicmine_unfiltered.tsv') as f:
			moo = "matchingid      evidencetype    gene_hugo_id    gene_entrez_id  gene_normalized cancer_id       cancer_normalized       drug_id drug_normalized variant_normalized      citation_count"
			headers = f.readline().strip('\n').split('\t')
			#for line in f:
				#rowDict = { h:v for h,v in zip(headers,line.strip('\n').split('\t')) }
				#seen[rowDict['cancer_text']] += 1
				#seen[rowDict['gene_text']] += 1
				#seen[rowDict['drug_text']] += 1
				#seen[rowDict['variant_text']] += 1
				#for v in line.strip('\n').split('\t'):
				#	seen[v] += 1
			#	vals = line.strip('\n').split('\t')
			#	seen += Counter(vals)

			print("before")
			seen = Counter(f.read().lower().replace('\t','\n').split('\n'))
			print("after")
			

		toSolve = sorted([ (count,term) for term,count in seen.items() if term in conflicting ],reverse=True)
		seen = None

		for j,(count,term) in enumerate(toSolve):
			conflicting_termids = sorted(lookup[term])

			print()
			print('#'*40 + " (%d/%d)" % (j+1,len(toSolve)))
			print("Clash: %s (%d)" % (term,count))
			print()
			possibleAnswers = {}
			for i,termid in enumerate(conflicting_termids):
				singleterm = singleterms[termid]
				synonyms = synonymLookup[termid]
				print("%d: %s [%s]" % (i,singleterm,termid))
				print("   %s" % str(synonyms))
				print()
				possibleAnswers[singleterm.lower()] = i


			if term in possibleAnswers:
				response = possibleAnswers[term]
			else:
				response = None
				selections = list(map(str,range(len(conflicting_termids))))
				allowed = set(['x','s'] + selections + allTermIDs)
				while not response in allowed:
					response = input('Which one to keep? (x to skip, s to add to stopwords or ID to add to another term) ')

			if response == 'x':
				continue
			elif response == 's':
				print("Adding to %s stoplist" % termtype)
				filename = 'stopwords_%s.txt' % termtype
				with open(filename,'a') as outF:
					outF.write("%s\n" % term)

				continue
			elif response in allTermIDs:
				termid = response
				response = -1

				singleterm = singleterms[termid]
				print("Adding to %s %s [%s]" % (termtype, singleterm, termid))
				filename = 'additions_%s.tsv' % termtype
				with open(filename,'a') as outF:
					outF.write("%s\t%s\t%s\n" % (termid,singleterm,term))

			response = int(response)

			#toDelete = "|".join(sorted(intersection))
			for i,termid in enumerate(conflicting_termids):
				if i == response:
					continue
				singleterm = singleterms[termid]
				filename = 'deletions_%s.tsv' % termtype
				with open(filename,'a') as outF:
					outF.write("%s\t%s\t%s\n" % (termid,singleterm,term))

		#break
