import kindred
import argparse
import pickle
import os

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Build and save a classifier')
	parser.add_argument('--inTrain',type=str,required=True,help='Directory with corpus in standoff format')
	parser.add_argument('--outDir',type=str,required=True,help='Directory to store output model files')

	args = parser.parse_args()

	relationInfo = []
	relationInfo.append(('AssociatedVariant',0.6,('gene','variant')))
	relationInfo.append(('Diagnostic',0.7,('cancer','gene')))
	relationInfo.append(('Predictive',0.92,('cancer','drug','gene')))
	relationInfo.append(('Prognostic',0.7,('cancer','gene')))
	relationInfo.append(('Predisposing',0.96,('cancer','gene')))

	for relationType,threshold,entityTypes in relationInfo:
		print("Building %s model" % relationType)
		print("  Loading training")
		trainCorpus = kindred.loadDir(dataFormat='standoff',directory=args.inTrain)
	
		for doc in trainCorpus.documents:
			doc.relations = [ r for r in doc.relations if r.relationType == relationType ]
			doc.relations = [ r for r in doc.relations if len(r.entityIDs) == len(entityTypes) ]

		print("  Doing training")
		classifier = kindred.RelationClassifier(classifierType='LogisticRegression',threshold=threshold,entityCount=len(entityTypes),acceptedEntityTypes=[entityTypes])
		classifier.train(trainCorpus)

		print("  Saving classifer")
		outModel = os.path.join(args.outDir,"%s.model" % relationType)
		with open(outModel,'wb') as f:
			pickle.dump(classifier,f)

		print("  Output done!")

