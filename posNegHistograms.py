#AML
#graphing histograms of positive and negative raw APARENET values


import pyensembl
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 
import sklearn.metrics as metrics
from Bio import SeqIO
import random
import time



def extractPredictionValuesSeparateArrays(name, negatives, positives, pasType = ""):
	#returns negative values and positive values
	dummyNegatives = np.array([])
	dummyPositives = np.array([])
	predName = "chr" + name
	forward, reverse = openForwardReverse("../../aparentGenomeTesting/chromosomePredictions50/", predName)
	flippedReverse = np.flip(reverse)
	for index, row in negatives.iterrows():
		startInt = int(row['start'])
		endInt = int(row['end'])
		if row['strand'] == "+":
			if pasType == "":
				#forward strand
				#negatives are 0-indexed still
				sliceVals = forward[startInt:endInt+ 1]
				dummyNegatives = np.concatenate((dummyNegatives,sliceVals))
			else:
				if row['type'] == pasType:
					#forward strand
					#negatives are 0-indexed still
					sliceVals = forward[startInt:endInt+ 1]
		elif row['strand'] == "-":
			if pasType == "":
				#reverse strand
				sliceVals = flippedReverse[startInt: endInt + 1]
				dummyNegatives = np.concatenate((dummyNegatives,sliceVals))
			else:
				if row['type'] == pasType:
					#reverse strand
					sliceVals = flippedReverse[startInt: endInt + 1]
					dummyNegatives = np.concatenate((dummyNegatives,sliceVals))
		else:
			print ("ERROR!  Strand not listed for negative example")
	for index,row in positives.iterrows():
		startInt = int(row['start'])
		endInt = int(row['end'])
		if row['strand'] == "+":
			if pasType == "":
				#forward strand
				#negatives are 0-indexed still
				sliceVals = forward[startInt:endInt+ 1]
				dummyPositives = np.concatenate((dummyPositives,sliceVals))
			else:
				if row['type'] == pasType:
					#forward strand
					#negatives are 0-indexed still
					sliceVals = forward[startInt:endInt+ 1]
					dummyPositives = np.concatenate((dummyPositives,sliceVals))
		elif row['strand'] == "-":
			if pasType == "":
				#reverse strand
				flippedReverse = flippedReverse[startInt: endInt + 1]
				dummyPositives = np.concatenate((dummyPositives,sliceVals))
			else:
				if row['type'] == pasType:
					#reverse strand
					flippedReverse = flippedReverse[startInt: endInt + 1]
					dummyPositives = np.concatenate((dummyPositives,sliceVals))
		else:
			print ("ERROR!  Strand not listed for negative example")
	return dummyNegatives, dummyPositives
	

def extractAllPositiveAndNegativePredictionValues(names, stem, b, s, pasType = ""):
	posVals = np.array([])
	negVals = np.array([])
	for name in names:
		fileName = stem + "chro" + name + "_NegSpaces" + str(b) + "_shifted" + str(s) + "Nts"
		negatives = openBalancedNegatives(fileName)
		positives = openBalancedPositives(fileName)
		dummyNegatives, dummyPositives = extractPredictionValuesSeparateArrays(name, negatives, positives, pasType)
		posVals = np.concatenate((posVals, dummyPositives))
		negVals = np.concatenate((negVals, dummyNegatives))
	return negVals, posVals

def openBalancedNegatives(name):
	negativesName = name + "BalancedNegatives.csv"
	#colNames = ["seqName", "start", "end", "clusterID", "strand", "type", "side"]
	return pd.read_csv( negativesName, dtype = {"seqName": str}) 

def openBalancedPositives(name):
	positivesName = name + "BalancedPositives.csv"
	#colNames = ["seqName", "start", "end", "clusterID", "strand", "type", "side"]
	return pd.read_csv( positivesName, dtype = {"seqName": str}) 
	
def openForwardReverse(stem, name):
	totalNameFor = stem + name + ".npy"
	print (totalNameFor)
	forward = np.load(totalNameFor)
	reverse = np.load(stem + name + "RC.npy")
	return forward, reverse

def placePeaksWithTolerance(peaks, clusters, tolerance, sign, lenSeq):
	clusterRanges = clusters.keys()
	for peak in peaks:
		if sign == "-":
			peak = flipSequenceIndex(peak, lenSeq)
		placed = False
		for rng in clusterRanges:
			if rng != 'Unplaced':
				lower = rng[0] - tolerance
				upper = rng[1] + tolerance
				if peak >= lower and peak <= upper: #if the peak is in [start, end]
					clusters[rng].append(peak)
					placed = True
					break
		if not placed: #wasn't placed
			clusters['Unplaced'].append(peak)
	return clusters


#using the peak-cluster dictionaries 
#need to count TP, FP, FN (placed peaks, unplaced peaks, empty clusters)
def fillConfMatrix(dictForward, dictRC):
	countTP = 0 #peaks in true cluster
	countFP = 0 #peak in false cluster 
	countFN = 0 #no peak in true cluster
	countTN = 0 #no peak in false cluster
	for key in dictForward:
		if key != 'Unplaced':
			inCluster = len(dictForward[key])
			if inCluster != 0:
				countTP += inCluster
				#print (key, " contains: ", inCluster)
			else:
				countFN += 1
		else: #unplaced peaks
			countFP += len(dictForward[key])
	for key in dictRC:
		if key != 'Unplaced':
			inCluster = len(dictRC[key])
			if inCluster != 0:
				countTP += inCluster
				#print (key, " contains: ", inCluster)
			else:
				countFN += 1
		else: #unplaced peaks
			countFP += len(dictRC[key])
	return countTP, countFP, countFN
	
def calculatePrecisionRecall(tp, fp, fn):
	precision = tp / (tp + fp)
	recall = tp / (tp + fn)
	if tp == 0:
		return precision, recall, None
	else: 
		f1 =  (2 * (precision * recall))/(precision + recall)
	return precision, recall, f1
	

def openTrueValuesForType(name, pasType):
	
	
	clustersForward = {}
	clustersRC = {} #key of (Start, End) will use to track clusters with no peaks for the FN  
	if pasType == "All":
		for index, row in currentTrueVals.iterrows():
			if row['strand'] == "+": #forward strand
				clustersForward[(row['start'], row['end'])] = []
			else: #negative strand, RC cluster
				clustersRC[(row['start'], row['end'])] = []
		clustersForward['Unplaced'] = []
		clustersRC['Unplaced'] = []
	else:
		
		maskType = currentTrueVals["type"] == pasType
		maskedTrue = currentTrueVals[maskType]
		for index, row in maskedTrue.iterrows():
			if row['strand'] == "+": #forward strand
				clustersForward[(row['start'], row['end'])] = []
			else: #negative strand, RC cluster
				clustersRC[(row['start'], row['end'])] = []
		clustersForward['Unplaced'] = []
		clustersRC['Unplaced'] = []
	#print (clustersForward)	
	return clustersForward, clustersRC


stem = "./datasets/"
'''
negatives = openBalancedNegatives(stem)
positives = openBalancedPositives(stem)
negativeVals, positiveVals = extractPredictionValuesSeparateArrays("21", negatives, positives)
'''
names = ["1","2","3","4","5","6","7","8","9","10","11","12","13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "Y"]
negativeVals, positiveVals = extractAllPositiveAndNegativePredictionValues(names, stem, 1, 50, "IN")
overallMax = max([max(negativeVals), max(positiveVals)])
overallMin = min([min(negativeVals), min(positiveVals)])
bins = np.linspace(overallMin, overallMax, 1000)
plt.hist(negativeVals, bins, alpha = 0.5, label = 'negative values')
plt.hist(positiveVals, bins, alpha = 0.5, label = 'positive values')
plt.legend()
plt.show()
