import argparse
import re

def main():
	parser = argparse.ArgumentParser('Extend cancer wordlist with syndromes and extra terms')
	parser.add_argument('--old',required=True,type=str,help='Base version of cancer list')
	parser.add_argument('--syndromes',required=True,type=str,help='Syndrome list')
	parser.add_argument('--extra',required=True,type=str,help='Extra list of terms')
	parser.add_argument('--outFile',required=True,type=str,help='Output term list')
	args = parser.parse_args()

	id_to_main = {}
	id_to_synonyms = {}

	existing_synonyms = set()

	print("Adding syndromes...")
	with open(args.syndromes) as f:
		for i,line in enumerate(f):
			term_id = "SYNDROME_%04d" % (i+1)
			terms = line.strip().split('|')

			main = terms[0]
			synonyms = [ s.lower().strip() for s in terms ]
			synonyms += [ s.replace('syndrome','syndromes') for s in synonyms ]
			synonyms += [ s.replace('-',' ') for s in synonyms ]
			synonyms += [ s.replace('-','- ') for s in synonyms ]
			synonyms += [ s.replace('-',' - ') for s in synonyms ]
			synonyms = [ re.sub(r'\ +', ' ', s) for s in synonyms ]
			synonyms = sorted(set(synonyms))

			id_to_main[term_id] = main
			id_to_synonyms[term_id] = synonyms

			existing_synonyms.update(synonyms)
	print(f"Got {len(id_to_main)} cancer types\n")

	print("Adding custom list of extra cancer terms...")
	with open(args.extra) as f:
		for line in f:
			term_id, main, synonyms = line.strip().split('\t')
			synonyms = synonyms.split('|')

			assert not term_id in id_to_main
			for s in synonyms:
				assert not s in existing_synonyms
			existing_synonyms.update(synonyms)

			id_to_main[term_id] = main
			id_to_synonyms[term_id] = synonyms
	print(f"Got {len(id_to_main)} cancer types\n")

	print("Loading prior list of cancer types...")
	with open(args.old) as f:
		for line in f:
			term_id, main, synonyms = line.strip().split('\t')
			synonyms = synonyms.split('|')

			synonyms = [ s for s in synonyms if not s in existing_synonyms ]

			if synonyms:
				id_to_main[term_id] = main
				id_to_synonyms[term_id] = synonyms
	print(f"Got {len(id_to_main)} cancer types\n")

	age_terms = ['adult','juvenile','pediatric','paediatric','childhood','adolescent','infantile']
	age_regexes = [ re.compile(r'\b%s\b' % t, re.I) for t in age_terms ]

	print("Filtering out cancer terms with age terms in their main name")
	id_to_main = { term_id:main for term_id,main in id_to_main.items() if not any( regex.search(main) for regex in age_regexes ) }
	print(f"Got {len(id_to_main)} cancer types\n")

	print("Adding age terms in to every synonym")
	for term_id in id_to_main:
		synonyms = id_to_synonyms[term_id]
		synonyms += [ f"{at} {s}" for s in synonyms for at in age_terms ]
		id_to_synonyms[term_id] = sorted(synonyms)

	print("Outputting...")
	with open(args.outFile,'w') as outF:
		for term_id in id_to_main:
			main = id_to_main[term_id]
			synonyms = "|".join(id_to_synonyms[term_id])

			outF.write(f"{term_id}\t{main}\t{synonyms}\n")

	print("Done")

if __name__ == '__main__':
	main()

