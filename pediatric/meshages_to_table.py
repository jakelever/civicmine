import argparse
import json
import pandas as pd

def main():
	parser = argparse.ArgumentParser('Convert a JSON file with MeSH age data into a nice table')
	parser.add_argument('--inJSON',required=True,type=str,help='Input JSON file')
	parser.add_argument('--outTSV',required=True,type=str,help='Output TSV file')
	args = parser.parse_args()

	with open(args.inJSON) as f:
		data = json.load(f)

	reformatted = []
	for cancer_type,agedata in data.items():
		row = {}
		row['cancer'] = cancer_type
		row.update(agedata)
		reformatted.append(row)

	df = pd.DataFrame(reformatted)
	df.to_csv(args.outTSV, sep="\t", index=False)
	print(f"Saved {df.shape[0]} rows to {args.outTSV}")

if __name__ == '__main__':
	main()
