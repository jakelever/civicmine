import kindred
import argparse
import pickle
import os

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Build and save a classifier')
	parser.add_argument('--inTrain',type=str,required=True)
	#parser.add_argument('--outModel_Driver',type=str,required=True)
	#parser.add_argument('--outModel_Oncogene',type=str,required=True)
	#parser.add_argument('--outModel_TumorSuppressor',type=str,required=True)

	args = parser.parse_args()

	thresholds = {'Driver':0.89, 'Oncogene': 0.94, 'Tumor_Suppressor': 0.97}

	relationInfo = []
	relationInfo.append(('Diagnostic','Diagnostic',('cancer','gene')))
	relationInfo.append(('Predictive/Prognostic','Prognostic',('cancer','gene')))
	relationInfo.append(('Predictive/Prognostic','Predictive',('cancer','drug','gene')))
	relationInfo.append(('Predisposing','Predisposing',('cancer','gene')))

	#for relationType,outModel in zip(['Driver','Oncogene','Tumor_Suppressor'], [args.outModel_Driver,args.outModel_Oncogene,args.outModel_TumorSuppressor] ):
	#for relationType in ['Diagnostic','Predictive/Prognostic','Predisposing']:
	for relationType,replacementType,entityTypes in relationInfo:
		print("Building %s model" % relationType)
		print("  Loading training")
		goldDir = 'gold'
		trainCorpus = kindred.loadDir(dataFormat='standoff',directory=args.inTrain)
	
		for doc in trainCorpus.documents:
			doc.relations = [ r for r in doc.relations if r.relationType == relationType ]
			doc.relations = [ r for r in doc.relations if len(r.entityIDs) == len(entityTypes) ]
			for r in doc.relations:
				r.relationType = replacementType

		print("  Doing training")
		features = "entityTypes,unigramsBetweenEntities,bigrams,dependencyPathEdges,dependencyPathEdgesNearEntities".split(',')
		threshold = 0.7 #thresholds[relationType]
		classifier = kindred.RelationClassifier(classifierType='LogisticRegression',threshold=threshold,features=features,entityCount=len(entityTypes),acceptedEntityTypes=[entityTypes])
		classifier.train(trainCorpus)

		print("  Saving classifer")
		outModel = os.path.join('models',"%s.model" % replacementType)
		with open(outModel,'wb') as f:
			pickle.dump(classifier,f)

		print("  Output done!")

