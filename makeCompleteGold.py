
import os
import sys

mainDir = 'trainingFixed'
variantDir = 'gold'
outDir = 'completegold'

def loadA2(filename):
	data = []
	with open(filename) as f:
		for line in f:
			relationWithoutID = line.strip().split('\t')[1]

			if "drug" in relationWithoutID:
				relationWithoutID = relationWithoutID.replace('Predictive/Prognostic','Predictive')
			else:
				relationWithoutID = relationWithoutID.replace('Predictive/Prognostic','Prognostic')
			relationWithoutID = relationWithoutID.replace('omicevent','variant')
			relationWithoutID = relationWithoutID.replace('Yes','AssociatedVariant')

			data.append(relationWithoutID)

	return data

def mergeA2(firstA2,secondA2,merged):
	firstRelations = loadA2(firstA2)
	secondRelations = loadA2(secondA2)
	with open(merged,'w') as f:
		for i,r in enumerate(firstRelations + secondRelations):
			line = "R%d\t%s\n" % (i+1,r)
			f.write(line)

for filename in sorted(os.listdir(variantDir)):
	if not filename.endswith('.a2'):
		continue
	print(filename)

	mainA2 = os.path.join(mainDir,filename.replace('.0','0'))
	variantA2 = os.path.join(variantDir,filename)
	assert os.path.isfile(mainA2)
	assert os.path.isfile(variantA2)

	outA2 = os.path.join(outDir,filename)
	mergeA2(mainA2,variantA2,outA2)
	#break
