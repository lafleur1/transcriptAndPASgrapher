#Alyssa La Fleur
#2/27/20
#Comparing APARENT and DeepPASTA Output for Completed Data
import pyensembl
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 
import sklearn.metrics as metrics
from Bio import SeqIO
import random
import time
import math
from scipy.signal import find_peaks
import sys
from csv import writer
from ast import literal_eval



#graph the predicted values in the positive/negative clusters from the balanced datasets
#not using maximum peaks, just the values predicted in those cluster ranges
#inputs: trueValsArray (the boolean labels for the predicted values in the positive and negative clusters) (numpy array), preds (the predicted values in the positive/negative clusters) (numpy array)
#outputs: fpr (false positive rates from roc_curve),tpr (true positive rates from roc_curve),thresholds (thresholds from scipy roc_curve), auc_score (auc score for the ROC curve), prec (precisions for PR curve), rec (recalls for PR curve), thresholdsPR (thresholds from precision_recall_curve), auprc_score (average precision score)
def computeAndGraphAllROCs(trueValsArray, preds):
	#set up true value labels
	print ("Total true labelled positions: ", np.sum(trueValsArray))
	######
	#Compute ROC curve and ROC AUC for each stride length
	fpr, tpr, thresholds = metrics.roc_curve(trueValsArray, preds)
	prec, rec, thresholdsPR = metrics.precision_recall_curve(trueValsArray, preds)
	auc_score = metrics.roc_auc_score(trueValsArray,preds)
	auprc_score = metrics.average_precision_score(trueValsArray, preds)
	posDistFPR, posDistTPR, posDistThresh, x  = maxPosDist45DegreeLine(fpr,tpr,thresholds)
	print ("Best threshold by 45 degree line distance: ", posDistThresh)
	equalFPR, equalTPR, equalThresh, x  = findSpecifictySensitivityEqualityPoint(fpr,tpr,thresholds)
	print ("Best threshold by specificity Sensisitivity Equality: ", equalThresh)
	closeFPR, closeTPR, closeThresh, x  = minDistanceTopLeftCorner(fpr,tpr,thresholds)
	print ("Best threshold by closest to top left corner: ", closeThresh)
	return fpr,tpr,thresholds,auc_score, prec, rec, thresholdsPR, auprc_score


#creates ROC and PR curves for genome, using default values unless changed 
#TODO: add X to the total graph creation once I find those numpy values
def makeGraphsAverageCleavageValues(b = 1, s = 50, typePas = ""):
	chromosomes = ["1","2","3","4","5","6","7","8","9","10","11","12","13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "Y"]
	values, bools = extractAllPredictionValues(chromosomes, "./datasets/", b, s, typePas)
	fpr,tpr,thresholds,auc_score, prec, rec, thresholdsPR, auprc_score = computeAndGraphAllROCs(bools, values)
	print ("AUC Score: ", auc_score)
	print ("Avg Precision Score: ", auprc_score)
	plt.plot(fpr, tpr)
	plt.title(typePas + " ROC")
	plt.xlabel("FPR")
	plt.ylabel("TPR")
	plt.show()
	plt.plot(prec, rec)
	plt.xlabel("Recall")
	plt.ylabel("Precision")
	plt.title(typePas + " PR")
	plt.show()


#https://stackoverflow.com/questions/32742976/how-to-read-a-column-of-csv-as-dtype-list-using-pandas 
def readDeepPASTACSV(name, b, s):
	fileName = "./predictions/chro" + name + "_NegSpaces" + str(b) + "_shifted" + str(s) + "Nts"
	return pd.read_csv(fileName + "BalancedPositives.csv", converters = {"DeepPASTAPredictions": literal_eval}), pd.read_csv(fileName + "BalancedNegatives.csv", converters = {"DeepPASTAPredictions": literal_eval})	
	
#list of chromosomes 15-Y (when all are finished)
def openDeepPASTAForGraphing(dataframe, pasType = "", usingStrand = ""):
	dummyPlusStrand = np.array([])
	dummyMinusStrand = np.array([])
	numberNeg1 = [0,0]
	lengths = [[],[],[],[]]
	if pasType != "":
		mask = dataframe['type'] == pasType
		dataframe = dataframe[mask]
	if usingStrand != "":
		mask = dataframe['strand'] == usingStrand
		dataframe = dataframe[mask]
	for index, row in dataframe.iterrows():
		if row['strand'] == "+":
			#print (row['DeepPASTAPredictions'])
			if row['DeepPASTAPredictions'] != [-1]:
				sliceVals = np.asarray(row['DeepPASTAPredictions'])
				dummyPlusStrand = np.concatenate((dummyPlusStrand,sliceVals)) 
				lengths[0].append(float(row['end']) - float(row['start']))
			else:
				#print ("FOUND -1")
				numberNeg1[0] += 1
				lengths[1].append(float(row['end']) - float(row['start']))
		elif row['strand'] == "-":
			#print (row['DeepPASTAPredictions'])
			if row['DeepPASTAPredictions'] != [-1]:
				sliceVals = np.asarray(row['DeepPASTAPredictions'])
				dummyMinusStrand = np.concatenate((dummyMinusStrand,sliceVals)) 
				lengths[2].append(float(row['end']) - float(row['start']))
			else:
				#print ("FOUND -1")
				numberNeg1[1] += 1
				lengths[3].append(float(row['end']) - float(row['start']))
	return dummyPlusStrand, dummyMinusStrand #, numberNeg1, lengths


def extractPredictionValuesSeparateArrays(negatives, positives, pasType = "", usingStrand = ""):
	positiveForward, positiveReverse = openDeepPASTAForGraphing(positives, pasType, usingStrand)
	negativeForward, negativeReverse =  openDeepPASTAForGraphing(negatives, pasType, usingStrand)
	if usingStrand == "+":
		#return only forward strand positive and negative values
		return negativeForward, positiveForward
	elif usingStrand == "-":
		#return only reverse strand positive and negative values
		return negativeReverse, positiveReverse
	else:
		#return both forward/reverse strand positive and negative values
		return np.concatenate((negativeForward, negativeReverse)), np.concatenate((positiveForward, positiveReverse))


def extractAllPositiveAndNegativePredictionValues(names, pasType = "", usingStrand = ""):
	posVals = np.array([])
	negVals = np.array([])
	for name in names:
		print ("ON: ", name)
		positives, negatives = readDeepPASTACSV(name,1 , 50)
		dummyNegatives, dummyPositives = extractPredictionValuesSeparateArrays(negatives, positives, pasType, usingStrand)
		posVals = np.concatenate((posVals, dummyPositives))
		negVals = np.concatenate((negVals, dummyNegatives))
	return negVals, posVals


def createBoxAndWhiskerForAllPASTypes(chromosomes):
	pasTypes = ['', 'IN', 'TE', 'IG', 'AI', 'EX', 'DS', 'AE', 'AU']
	toPlot = []
	tick_labels = []
	for typePAS in pasTypes:
		print ("ON PAS TYPE: ", typePAS)
		stem = "./datasets/"
		negVals, posVals = extractAllPositiveAndNegativePredictionValues(chromosomes, pasType = typePAS)
		toPlot.append(negVals)
		toPlot.append(posVals)
		tick_labels.append("Negative " + typePAS)
		tick_labels.append("Positive " + typePAS)
	fig = plt.figure(1, figsize=(9, 6))
	ax = fig.add_subplot(111)
	bp = ax.boxplot(toPlot)
	ax.set_xticklabels(tick_labels)
	plt.title("Distribution of DeepPASTA Prediction Values")
	plt.show()

#currently finished with:
chromosomes = ["15", "16", "18", "20", "21", "22", "Y"]
print ("TEST ON Y")
createBoxAndWhiskerForAllPASTypes(chromosomes)



negatives, positives = extractAllPositiveAndNegativePredictionValues(chromosomes, pasType = "", usingStrand = "")

#box and whisker
fig = plt.figure(1, figsize=(9, 6))
ax = fig.add_subplot(111)
bp = ax.boxplot([negatives, positives])
ax.set_xticklabels(['negatives', 'positives'])
plt.show()


'''
#histograms
overallMax = max([max(negatives), max(positives)])
overallMin = min([min(negatives), min(positives)])
bins = np.linspace(overallMin, overallMax, 100)
plt.hist(negatives, bins, alpha = 0.5, label = 'Negatives')
plt.hist(positives, bins, alpha = 0.5, label = 'Positives')
plt.legend()
plt.show()

'''



'''
dummyPlusStrand, dummyMinusStrand, numberFail, failedLength = openDeepPASTAForGraphing(balPosY)
#print (dummyPlusStrand)
print (balPosY.shape)
print ("Failed on +, failed on -: ", numberFail)
print ("Failed lengths: ")
fig = plt.figure(1, figsize=(9, 6))
ax = fig.add_subplot(111)
bp = ax.boxplot(failedLength)
ax.set_xticklabels(['Passed +', 'Failed +', 'Passed -', 'Failed -'])
plt.show()


overallMax = max([max(failedLength[0]), max(failedLength[1]), max(failedLength[2]), max(failedLength[3])])
overallMin = min([min(failedLength[0]), min(failedLength[1]), min(failedLength[2]), min(failedLength[3])])
bins = np.linspace(overallMin, overallMax, 100)
plt.hist(failedLength[0], bins, alpha = 0.5, label = 'Passed +')
plt.hist(failedLength[1], bins, alpha = 0.5, label = 'Failed +')
plt.hist(failedLength[2], bins, alpha = 0.5, label = 'Passed -')
plt.hist(failedLength[3], bins, alpha = 0.5, label = 'Failed -')
plt.legend()
plt.show()

dummyPlusStrand, dummyMinusStrand, numberFail, failedLength = openDeepPASTAForGraphing(balNegY)
#print (dummyPlusStrand)
print (balPosY.shape)
print ("Failed on +, failed on -: ", numberFail)
print ("Failed lengths: ")
fig = plt.figure(1, figsize=(9, 6))
ax = fig.add_subplot(111)
bp = ax.boxplot(failedLength)
ax.set_xticklabels(['Passed +', 'Failed +', 'Passed -', 'Failed -'])
plt.show()

overallMax = max([max(failedLength[0]), max(failedLength[1]), max(failedLength[2]), max(failedLength[3])])
overallMin = min([min(failedLength[0]), min(failedLength[1]), min(failedLength[2]), min(failedLength[3])])
bins = np.linspace(overallMin, overallMax, 100)
plt.hist(failedLength[0], bins, alpha = 0.5, label = 'Passed +')
plt.hist(failedLength[1], bins, alpha = 0.5, label = 'Failed +')
plt.hist(failedLength[2], bins, alpha = 0.5, label = 'Passed -')
plt.hist(failedLength[3], bins, alpha = 0.5, label = 'Failed -')
plt.legend()
plt.show()
'''

