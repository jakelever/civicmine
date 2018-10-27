import argparse
import sys
import kindred
import random
import numpy as np
import sklearn.metrics
import os
import random

if __name__ == '__main__':
	entityTypes = {}
	entityTypes['AssociatedVariant'] = ('gene','variant')
	entityTypes['Diagnostic'] = ('cancer','gene')
	entityTypes['Predictive'] = ('cancer','drug','gene')
	entityTypes['Prognostic'] = ('cancer','gene')
	entityTypes['Predisposing'] = ('cancer','gene')

	reltypes = ",".join(sorted(list(entityTypes.keys())))

	parser = argparse.ArgumentParser(description='Calculate classifier metrics for a particulate relation type and data set size')
	parser.add_argument('--inDir',type=str,required=True,help='Directory containing corpus to load (in standoff format)')
	parser.add_argument('--reltype',type=str,required=True,help='Relation type to analyze. Must be one of %s' % reltypes)
	parser.add_argument('--iterations',type=int,required=False,default=1,help='Number of iterations to run')
	parser.add_argument('--outCurve',type=str,required=True,help='File to output curve data to')
	args = parser.parse_args()

	with open(args.outCurve,'w') as outF:
		outF.write("%s\t%s\t%s\n" % ('threshold','precision','recall'))

		for _ in range(args.iterations):
			corpus = kindred.loadDir('standoff',args.inDir)
			train,gold = corpus.split(0.75)
			
			# Trim back to relation type of choice
			for doc in train.documents:
				doc.relations = [ r for r in doc.relations if r.relationType == args.reltype ]
			for doc in gold.documents:
				doc.relations = [ r for r in doc.relations if r.relationType == args.reltype ]

			entityType = entityTypes[args.reltype]
			entityCount = len(entityType)
			classifier = kindred.RelationClassifier(classifierType='LogisticRegression',threshold=-100.0,entityCount=entityCount,acceptedEntityTypes=[entityType])

			classifier.train(train)

			predictions = gold.clone()
			predictions.removeRelations()

			classifier.predict(predictions)
			
			scores = []
			for goldDoc,testDoc in zip(gold.documents,predictions.documents):
				goldTuples = { (r.relationType,tuple(r.entityIDs)) for r in goldDoc.relations }
				testTuples = { (r.relationType,tuple(r.entityIDs)):r.probability for r in testDoc.relations }

				tmpScores = [ (99999.9,True) for k in goldTuples if not k in testTuples ]
				tmpScores += [ (prob,k in goldTuples) for k,prob in testTuples.items() ]

				scores += tmpScores

			scores = sorted(scores,reverse=True)

			for threshold in [-0.1] + list(np.arange(0,1,0.01)) + [1.0]:
				TP = sum( 1 for score,trueOrFalse in scores if score > threshold and trueOrFalse == True )
				FP = sum( 1 for score,trueOrFalse in scores if score > threshold and trueOrFalse == False )
				FN = sum( 1 for score,trueOrFalse in scores if score <= threshold and trueOrFalse == True )

				if (TP+FP) != 0 and (TP+FN) != 0:
					precision = TP / float(TP+FP)
					recall = TP / float(TP+FN)

					outF.write("%f\t%f\t%f\n" % (threshold,precision,recall))
					#outF.write("TP:%d\tFP:%d\tFN:%d\t%f\t%f\t%f\n" % (TP,FP,FN,threshold,precision,recall))
					print(TP,FP,FN,threshold,precision,recall)
					sys.stdout.flush()

