import argparse
import json
from tqdm import tqdm
import os
from collections import defaultdict,Counter
import itertools

def main():
	parser = argparse.ArgumentParser(description='Load up MeSH cancer age data and aggregate it')
	parser.add_argument('--meshDataDir',required=True,type=str,help='Directory containing MeSH data on different cancer types')
	parser.add_argument('--outFile',required=True,type=str,help='Output file')
	args = parser.parse_args()
	
	input_files = sorted( f for f in os.listdir(args.meshDataDir) if f.endswith('.json') )

	print("Loading paper data...")
	data_by_paper = {}
	for input_file in tqdm(input_files):
		with open(os.path.join(args.meshDataDir,input_file)) as f:
			tmp_data = json.load(f)
			data_by_paper.update(tmp_data)

	print(f"Loaded data from {len(data_by_paper)} papers")

	cancer_age_counts = defaultdict(Counter)

	print("Aggregating cancer age counts...")
	for pmid in data_by_paper:
		ages = data_by_paper[pmid]['ages']
		cancers = data_by_paper[pmid]['cancer_names']

		for age,cancer in itertools.product(ages,cancers):
			cancer_age_counts[cancer][age] += 1

	with open(args.outFile,'w') as f:
		json.dump(cancer_age_counts,f,indent=2,sort_keys=True)
	print("Done.")

if __name__ == '__main__':
	main()
	
