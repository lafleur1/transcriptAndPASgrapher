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


################################# ROC Curve Threshold stuff #############################################

def maxPosDist45DegreeLine(fpr,tpr, threshs):
	#fpr is x, tpr is y
	#find the fpr,tpr,and threshold of the ROC point with the maximum positive distance from the 45 degree line
	maxPosDist = -1
	maxIndex = 0
	for i in range(0,len(fpr)):
		if tpr[i] >=fpr[i]: #above or on the 45 degree line
			currPosDist = tpr[i] - fpr[i] #since 45 degree line is y=x
			if currPosDist >= maxPosDist:
				maxPosDist = currPosDist 
				maxIndex = i
	if maxPosDist == -1:
		return None
	else:
		return fpr[maxIndex], tpr[maxIndex], threshs[maxIndex], maxPosDist


def findSpecifictySensitivityEqualityPoint(fpr,tpr,threshs):
	#find the prediction closest to where sensitivity=specificity
	minDiff = math.sqrt(2) #maximum possible distance for the unit cube of the ROC curve
	minIndex = 0
	for i in range(0,len(fpr)):
		if fpr[i] != 0.0 and tpr[i] != 0.0: #will always choose (0,0) if not blocked from doing so
			se = tpr[i]
			sp = 1 - fpr[i]
			currDiff = math.fabs(se-sp)
			if currDiff < minDiff:
				minDiff = currDiff
				minIndex = i
	if minDiff != math.sqrt(2):
		return fpr[minIndex], tpr[minIndex], threshs[minIndex], minDiff
	else:
		return None
    
def minDistanceTopLeftCorner(fpr,tpr,threshs):
	#find the prediction closest to (1,1)
	minDist = math.sqrt(2) #maximum possible distance for the unit cube of the ROC curve
	minIndex = 0
	for i in range(0,len(fpr)):
		currDist = math.sqrt((fpr[i])**2 + (1-tpr[i])**2)
		#print (currDist)
		if currDist < minDist:
			minDist = currDist
			minIndex = i
	if minDist != math.sqrt(2):
		return fpr[minIndex], tpr[minIndex], threshs[minIndex], minDist
	else:
		return None
########################################################################################################################################


#################################   DeepPASTA Stuff ###################################################

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


def extractPredictionValuesSeparateArraysDeepPASTA(negatives, positives, pasType = "", usingStrand = ""):
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


def extractAllPositiveAndNegativePredictionValuesDeepPASTA(names, pasType = "", usingStrand = ""):
	posVals = np.array([])
	negVals = np.array([])
	for name in names:
		print ("ON: ", name)
		positives, negatives = readDeepPASTACSV(name,1 , 50)
		dummyNegatives, dummyPositives = extractPredictionValuesSeparateArraysDeepPASTA(negatives, positives, pasType, usingStrand)
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
		negVals, posVals = extractAllPositiveAndNegativePredictionValuesDeepPASTA(chromosomes, pasType = typePAS)
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





##############################################################################################################################################################################

#################################################################  APARENT ##########################################################################################


#opens forward and reverse complement APARENT predictions for a chromosome
#input: stem (file path the predictions are stored in), name (chromosome name)
#output: forward strand predictions, reverse complement strand predictions
def openForwardReverse(stem, name):
	totalNameFor = stem + name + ".npy"
	print (totalNameFor)
	forward = np.load(totalNameFor)
	reverse = np.load(stem + name + "RC.npy")
	return forward, reverse

def openBalancedNegatives(name):
	negativesName = name + "BalancedNegatives.csv"
	return pd.read_csv( negativesName, dtype = {"seqName": str}) 

def openBalancedPositives(name):
	positivesName = name + "BalancedPositives.csv"
	return pd.read_csv( positivesName, dtype = {"seqName": str}) 

#for extracting the values from a prediction numpy for a given chromosome
#input: name (chromosome name to open values for and extract true/false cluster range values from), negatives (dataset to use to get values from numpy array), positives (dataset to use to get values from numpy array), pasType (does all pasTypes if not set.) 
#outputs: dummybools, dummy #return all positve/negative values in one array, boolean labels for those values in other array
def extractPredictionValuesAPARENT(name, negatives, positives, pasType = ""):
	#making a sklearn AUC and AUPRC plot to look at different threshold values
	dummy = np.array([])
	dummybools = np.array([])
	predName = "chr" + name
	forward, reverse = openForwardReverse("../../../aparentGenomeTesting/chromosomePredictions50/", predName) #open APARENT prediction values for the current chromosome
	flippedReverse = np.flip(reverse) #flip reverse strand so that predictions are in the same order as the forward 5'->3' direction since all the clusters are indexed off of the forward 5'->3' in polyAsite2.0
	#extract negative values
	for index, row in negatives.iterrows():
		startInt = int(row['start'])
		endInt = int(row['end'])
		if row['strand'] == "+":
			if pasType == "":
				#forward strand
				#negatives are 0-indexed still
				sliceVals = forward[startInt:endInt+ 1]
				dummy = np.concatenate((dummy,sliceVals))
				dummybools = np.concatenate((dummybools, np.zeros(sliceVals.size)))
			else:
				if row['type'] == pasType:
					#forward strand
					#negatives are 0-indexed still
					sliceVals = forward[startInt:endInt+ 1]
					dummy = np.concatenate((dummy,sliceVals))
					dummybools = np.concatenate((dummybools, np.zeros(sliceVals.size)))
		elif row['strand'] == "-":
			if pasType == "":
				#reverse strand
				sliceVals = flippedReverse[startInt: endInt + 1]
				dummy = np.concatenate((dummy,sliceVals))
				dummybools = np.concatenate((dummybools, np.zeros(sliceVals.size)))
			else:
				if row['type'] == pasType:
					#reverse strand
					sliceVals = flippedReverse[startInt: endInt + 1]
					dummy = np.concatenate((dummy,sliceVals))
					dummybools = np.concatenate((dummybools, np.zeros(sliceVals.size)))
		else:
			print ("ERROR!  Strand not listed for negative example")
	#extract positive values
	for index,row in positives.iterrows():
		startInt = int(row['start'])
		endInt = int(row['end'])
		if row['strand'] == "+":
			if pasType == "":
				#forward strand
				#negatives are 0-indexed still
				sliceVals = forward[startInt:endInt+ 1]
				dummy = np.concatenate((dummy,sliceVals))
				dummybools = np.concatenate((dummybools, np.ones(sliceVals.size)))
			else:	
				if row['type'] == pasType:
					#forward strand
					#negatives are 0-indexed still
					sliceVals = forward[startInt:endInt+ 1]
					dummy = np.concatenate((dummy,sliceVals))
					dummybools = np.concatenate((dummybools, np.ones(sliceVals.size)))
		elif row['strand'] == "-":
			if pasType == "":
				#reverse strand
				flippedReverse = flippedReverse[startInt: endInt + 1]
				dummy = np.concatenate((dummy,sliceVals))
				dummybools = np.concatenate((dummybools, np.ones(sliceVals.size)))
			else:
				if row['type'] == pasType:
					#reverse strand
					flippedReverse = flippedReverse[startInt: endInt + 1]
					dummy = np.concatenate((dummy,sliceVals))
					dummybools = np.concatenate((dummybools, np.ones(sliceVals.size)))
		else:
			print ("ERROR!  Strand not listed for negative example")
	return dummybools, dummy #return all positve/negative values in one array, boolean labels for those values in other array

#use above function to extract values for multiple chromosomes
#inputs: list of chromosome names to extract values for, stem (file path the balanced datasets are stored in), b (spacing value of how many copies of the cluster should fit between the negative cluster and the next true cluster), s (shift value, how far each positive cluster should be shifted to make the negative cluster), pasType (if values should only be extracted for one pasType)
#output: aray of all extracted values, array of all extracted boolean labels for use with scipy's AUC functions
def extractAllPredictionValuesAPARENT(names, stem, b, s, pasType = ""):
	values = np.array([])
	bools = np.array([])
	for name in names:
		fileName = stem + "chro" + name + "_NegSpaces" + str(b) + "_shifted" + str(s) + "Nts"
		negatives = openBalancedNegatives(fileName)
		positives = openBalancedPositives(fileName)
		boolsCurrent, valsCurrent = extractPredictionValuesAPARENT(name, negatives, positives, pasType) #use above function to open the predicted APARENT values for the chromosome clusters 
		values = np.concatenate((values, valsCurrent))
		bools = np.concatenate((bools, boolsCurrent))
	return values, bools
	



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
	values, bools = extractAllPredictionValues(chromosomes, "../datasets/", b, s, typePas)
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

############################################################################################

#Comparing ROC & PR Curves for both models across prediction types 

def compareAPARENTandDeepPASTA(chromosomes):
	#for all chromosomes that are finished
	pasTypes = ['', 'IN', 'TE', 'IG', 'AI', 'EX', 'DS', 'AE', 'AU']
	deepPASTAValues = []
	deepPASTABooleans = []
	aparentValues = []
	aparentBooleans = []
	for pType in pasTypes:
		#for each type, pull the prediction values for all chromsomes
		values, bools = extractAllPredictionValuesAPARENT(chromosomes, "../datasets/", b = 1, s = 50, pasType = pType )
		dvalues, dbools = negVals, posVals = extractAllPositiveAndNegativePredictionValuesDeepPASTA(chromosomes, pasType = pType)
		#APARENT 
		print ("PASTYPE: ", pType)
		fprA,tprA,thresholdsA,auc_scoreA, precA, recA, thresholdsPRA, auprc_scoreA = computeAndGraphAllROCs(bools, values)
		print ("	APARRENT AUC Score: ", auc_scoreA)
		print ("	APARENT Avg Precision Score: ", auprc_scoreA)
		fprD,tprD,thresholdsD,auc_scoreD, precD, recD, thresholdsPRD, auprc_scoreD = computeAndGraphAllROCs(dbools, dvalues)
		print ("	DeepPASTA AUC Score: ", auc_scoreD)
		print ("	DeepPASTA Avg Precision Score: ", auprc_scoreD)
		plt.figure()
		plt.plot(fprA, tprA, label = "APARENT")
		plt.plot(fprD, tprD, label = "DeepPASTA")
		plt.title(typePas + " ROC")
		plt.xlabel("FPR")
		plt.ylabel("TPR")
		plt.legend()
		plt.savefig(pType + "ROC_v1.png")
		plt.close()
		plt.figure()
		plt.plot(precA, recA, label = "APARENT")
		plt.plot(precD, recD, label = "DeepPASTA")
		plt.xlabel("Recall")
		plt.ylabel("Precision")
		plt.legend()
		plt.title(typePas + " PR")
		plt.savefig(pType + "PR_v1.png")
		plt.close()
	
		




#currently finished with:
chromosomes = ["15", "16", "17", "18", "19", "20", "21", "22", "Y"]
#print ("TEST ON Y")
#createBoxAndWhiskerForAllPASTypes(chromosomes)

compareAPARENTandDeepPASTA(chromosomes)

'''
#histograms
overallMax = max([max(negatives), max(positives)])
overallMin = min([min(negatives), min(positives)])
bins = np.linspace(overallMin, overallMax, 100)
plt.hist(negatives, bins, alpha = 0.5, label = 'Negatives')
plt.hist(positives, bins, alpha = 0.5, label = 'Positives')
plt.legend()
plt.show()


#box and whisker
fig = plt.figure(1, figsize=(9, 6))
ax = fig.add_subplot(111)
bp = ax.boxplot([negatives, positives])
ax.set_xticklabels(['negatives', 'positives'])
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

