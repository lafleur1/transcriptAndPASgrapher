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
from scipy.signal import find_peaks

def find_peaks_ChromosomeVersion(avgPreds, peak_min_height, peak_min_distance, peak_prominence):
	#forwardPeaks = find_peaks_ChromosomeVersion(forward, minh, dist, (0.01, None)) 
	peaks, _ = find_peaks(avgPreds, height=peak_min_height, distance=peak_min_distance, prominence=peak_prominence) 
	return peaks

def extractValueOfOneDataframe(name, dataframe, pasType = "", usingStrand = ""):
	#returns values separated by strand or not for one dataframe
	dummyPlusStrand = np.array([])
	dummyMinusStrand = np.array([])
	predName = "chr" + name
	forward, reverse = openForwardReverse("../../aparentGenomeTesting/chromosomePredictions50/", predName)
	flippedReverse = np.flip(reverse)
	if pasType != "":
		mask = dataframe['type'] == pasType
		dataframe = dataframe[mask]
	if usingStrand != "":
		mask = dataframe['strand'] == usingStrand
		dataframe = dataframe[mask]
	for index, row in dataframe.iterrows():
		startInt = int(row['start'])
		endInt = int(row['end'])
		if row['strand'] == "+":
			#add all on positive strand  to forward
			sliceVals = forward[startInt:endInt+ 1]
			dummyPlusStrand = np.concatenate((dummyPlusStrand,sliceVals)) 
		elif row['strand'] == "-":
			#add all on negative strand to reverse 
			sliceVals = flippedReverse[startInt:endInt+ 1]
			dummyMinusStrand = np.concatenate((dummyMinusStrand,sliceVals))
	return dummyPlusStrand, dummyMinusStrand

def extractPredictionValuesSeparateArrays(name, negatives, positives, pasType = "", usingStrand = ""):
	#returns negative values and positive values
	#predName = "chr" + name
	#forward, reverse = openForwardReverse("../../aparentGenomeTesting/chromosomePredictions50/", predName)
	#flippedReverse = np.flip(reverse)
	positiveForward, positiveReverse = extractValueOfOneDataframe(name, positives, pasType, usingStrand)
	negativeForward, negativeReverse =  extractValueOfOneDataframe(name, negatives, pasType, usingStrand)
	if usingStrand == "+":
		#return only forward strand positive and negative values
		return negativeForward, positiveForward
	elif usingStrand == "-":
		#return only reverse strand positive and negative values
		return negativeReverse, positiveReverse
	else:
		#return both forward/reverse strand positive and negative values
		return np.concatenate((negativeForward, negativeReverse)), np.concatenate((positiveForward, positiveReverse))
	

def extractAllPositiveAndNegativePredictionValues(names, stem, b, s, pasType = "", usingStrand = ""):
	posVals = np.array([])
	negVals = np.array([])
	for name in names:
		fileName = stem + "chro" + name + "_NegSpaces" + str(b) + "_shifted" + str(s) + "Nts"
		print ("ON: ", name)
		negatives = openBalancedNegatives(fileName)
		positives = openBalancedPositives(fileName)
		dummyNegatives, dummyPositives = extractPredictionValuesSeparateArrays(name, negatives, positives, pasType, usingStrand)
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

def createHistogramForPASType(typePAS = ""):
	stem = "./datasets/"
	names = ["1","2","3","4","5","6","7","8","9","10","11","12","13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "Y"]
	print ("Forward strand: ")
	negativeValsForward, positiveValsForward = extractAllPositiveAndNegativePredictionValues(names, stem, 1, 50, typePAS, "+")
	print ("Reverse strand:")
	negativeValsReverse, positiveValsReverse = extractAllPositiveAndNegativePredictionValues(names, stem, 1, 50, typePAS, "-")
	print ("Number positive on forward strand: ", positiveValsForward.size, " Number negative on forward strand: ", negativeValsForward.size)
	print ("Number positive on reverse strand: ", positiveValsReverse.size, " Number negative on reverse strand: ", negativeValsReverse.size)
	overallMax = max([max(negativeValsForward), max(positiveValsForward), max(negativeValsReverse), max(positiveValsReverse)])
	overallMin = min([min(negativeValsForward), min(positiveValsForward), min(negativeValsReverse), min(positiveValsReverse)])

	bins = np.linspace(overallMin, overallMax, 1000)
	plt.hist(negativeValsForward, bins, alpha = 0.5, label = 'negative values + strand')
	plt.hist(positiveValsForward, bins, alpha = 0.5, label = 'positive values + strand')
	plt.hist(negativeValsReverse, bins, alpha = 0.5, label = 'negative values - strand')
	plt.hist(positiveValsReverse, bins, alpha = 0.5, label = 'positive values - strand')
	plt.title("Distribution of " + typePAS + " Average Cleavage Values")
	plt.legend()
	plt.show()

createHistogramForPASType("TE")

