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

	threshold = 0.60

	relationInfo = []
	relationInfo.append(('Yes','AssociatedOmicEvent',('gene','omicevent')))

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
		classifier = kindred.RelationClassifier(classifierType='LogisticRegression',threshold=threshold,features=features,entityCount=len(entityTypes),acceptedEntityTypes=[entityTypes])
		classifier.train(trainCorpus)

		print("  Saving classifer")
		outModel = os.path.join('models',"%s.model" % replacementType)
		with open(outModel,'wb') as f:
			pickle.dump(classifier,f)

		print("  Output done!")

