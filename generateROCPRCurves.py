 #AML 2/18/20
 #opening CSVs and generating curves
 
 
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
from mpl_toolkits.mplot3d import Axes3D


#confusion matrix class to use when calculating the ROC and PR curves for APARNET
class peakConfusionMatrix(object):
	def __init__(self, countTP, countFP, countFN, countTN , threshold, unplacedPredictions = 0, totalPredictions = 0):
		self.truePositives = countTP
		self.falsePositives = countFP
		self.falseNegatives = countFN
		self.trueNegatives = countTN
		self.threshold = threshold
		self.unplacedPredictions = unplacedPredictions #number of unplaced predictions across the genomes(to use with dictionary based methods)
		self.totalPredictionsMade = totalPredictions #number of peaks at the threshold level across the genome (to use with dictionary based methods)
	
	def truePositiveRate(self):
	#AKA sensitivity, how good the model is at predicting the positive class when the actual outcome is positive
		if not (self.truePositives + self.falseNegatives) == 0:
			return self.truePositives/(self.truePositives + self.falseNegatives)
		else:
			print ("ERROR!  Dividing by 0 in TPR calculation")
			return -1
	
	def falsePositiveRate(self):
	#how often positive is predicted when actual outcome is false
		if not (self.falsePositives + self.trueNegatives) == 0 :
			return self.falsePositives/(self.falsePositives + self.trueNegatives)
		else:
			print ("ERROR! Dividing by 0 in FPR")
			return None
	
	def specificity(self):
		if not  (self.trueNegatives + self.falsePositives) == 0:
			return self.trueNegatives/ (self.trueNegatives + self.falsePositives)
		else:
			print ("ERROR! Dividing by 0 in specificity calc")
			return None
	
	def precision(self):
	#how good model is at predicting positive class
		if not (self.truePositives + self.falsePositives) == 0:
			return self.truePositives/(self.truePositives + self.falsePositives)
		else:
			print ("ERROR! Dividing by 0 in precision calc")
			return None
	
	def recall(self):
		#Same as sensitvity AKA TPR
		return self.truePositiveRate()
	
	#return fraction of peaks not falling in a true/false cluster across the entire genome
	def fractionUnplaced(self):
		if self.totalPredictionsMade != 0:
			return self.unplacedPredictions/self.totalPredictionsMade
		else:
			return 0 

#format of saved CSV:
'''
columns_dict = {'name': [], 'pasType': [], 'bufferVal':[], 'spacing':[], 'minh':[], 'distance':[], 'tolerance':[], 'TruePositives':[], 'FalsePositives':[], 'FalseNegatives':[], 'TrueNegatives':[]}
'''


def openCMDataOneChromosome(location, name):
	total_name = location + name + "_ConfusionMatrices.csv"
	return pd.read_csv( total_name)
	
	
def makeCMSForOneChromosome(location, name, pasType = "TE", bufferVal = 1, spacing = 50, distance = 1, tolerance = 0):
	#for the given values, generates the ROC and PR curves
	chromosomeCMValues = openCMDataOneChromosome(location, name)
	#make Confusion Matrices using peakConfusionMatrix from extractingNegatives.py
	confusionMatrices = []
	for index, row in chromosomeCMValues.iterrows():
		if row['pasType'] == pasType and row['bufferVal'] == bufferVal and row['spacing'] == spacing and row['distance'] == distance and row['tolerance'] == tolerance:
			#def __init__(self, countTP, countFP, countFN, countTN , unplacedPredictions = 0, totalPredictions = 0):
			confusionMatrices.append(peakConfusionMatrix(row['TruePositives'], row['FalsePositives'], row['FalseNegatives'], row['TrueNegatives'], row['minh']))
	#print ("number CMs selected: ", len(confusionMatrices))
	return confusionMatrices
	
def graphAllChrROC(location, names, pasType = "TE", bufferVal = 1, spacing = 50, distance = 1, tolerance = 0):
	#plotting all ROC curves on same curve
	for name in names: 
		confusionMatrices = makeCMSForOneChromosome(location, name, pasType, bufferVal, spacing, distance, tolerance)
		#recalls = [c.recall() for c in confusionMatrices] #x axis on Precision Recall Curve
		##precisions = [c.precision() for c in confusionMatrices] #y axis on Precision Recall Curve
		falsePositiveRates = [c.falsePositiveRate() for c in confusionMatrices] #x axis ROC curve
		truePositiveRates = [c.truePositiveRate() for c in confusionMatrices] #y axis ROC Curve
		plt.plot(falsePositiveRates, truePositiveRates, 'o--')	
	#plt.plot(falsePositiveRates, truePositiveRates, 'o--')
	plt.plot([0,1.0],[0,1.0], '--')
	plt.title("ROC Curve")
	plt.xlabel("FPR")
	plt.ylabel("TPR")
	plt.xlim(0,1.0)
	plt.ylim(0,1.0)
	#print (thresholds)
	plt.show()
	
def graphAllChrPR(location, names, pasType = "TE", bufferVal = 1, spacing = 50, distance = 1, tolerance = 0):
	#plotting all ROC curves on same curve
	for name in names: 
		confusionMatrices = makeCMSForOneChromosome(location, name, pasType, bufferVal, spacing, distance, tolerance)
		recalls = [c.recall() for c in confusionMatrices] #x axis on Precision Recall Curve
		precisions = [c.precision() for c in confusionMatrices] #y axis on Precision Recall Curve
		plt.plot(recalls, precisions, 'o--')
	#plot PR Curve
	plt.plot([0,1.0],[1.0,0], '--')
	plt.title("PR Curve")
	plt.xlabel("Recall")
	plt.ylabel("Precision")
	plt.ylim(0,1.0)
	plt.xlim(0,1.0)
	plt.show()
	

def generateGraphs(names, location, pasType = "TE", bufferVal = 1, spacing = 50, distance = 1, tolerance = 0):
	cmVals  = {} #key will be threshold values
	for name in names:
		#for each chromsome open stored CM information 
		chromosomeCMValues = openCMDataOneChromosome(location, name)
		foundMinHeights = []
		#make Confusion Matrices using peakConfusionMatrix from extractingNegatives.py
		for index, row in chromosomeCMValues.iterrows():
			if row['pasType'] == pasType and row['bufferVal'] == bufferVal and row['spacing'] == spacing and row['distance'] == distance and row['tolerance'] == tolerance:
				keyVal = round(row['minh'], 8)
				if keyVal not in cmVals:
					cmVals[keyVal] = [row['TruePositives'], row['FalsePositives'], row['FalseNegatives'], row['TrueNegatives'], row['minh']]
					foundMinHeights.append(keyVal)
				else:
					if keyVal not in foundMinHeights: 
						#update values 
						cmVals[keyVal][0] += row['TruePositives']
						cmVals[keyVal][1] += row['FalsePositives']
						cmVals[keyVal][2] += row['FalseNegatives']
						cmVals[keyVal][3] += row['TrueNegatives']
						foundMinHeights.append(keyVal)
	#now change all threshold values into confusion matrices
	confusionMatrices = []
	for key in cmVals:
		confusionMatrices.append(peakConfusionMatrix(cmVals[key][0], cmVals[key][1], cmVals[key][2], cmVals[key][3], key))
	recalls = [c.recall() for c in confusionMatrices] #x axis on Precision Recall Curve
	precisions = [c.precision() for c in confusionMatrices] #y axis on Precision Recall Curve
	falsePositiveRates = [c.falsePositiveRate() for c in confusionMatrices] #x axis ROC curve
	truePositiveRates = [c.truePositiveRate() for c in confusionMatrices] #y axis ROC Curve
	thresholds = [c.threshold for c in confusionMatrices]
	#graph ROC curve 
	fprs, tprs = zip(*sorted(zip(falsePositiveRates, truePositiveRates)))
	plt.plot(fprs, tprs, 'o--')
	plt.plot([0,1.0],[0,1.0], '--')
	plt.title("ROC Curve")
	plt.xlabel("FPR")
	plt.ylabel("TPR")
	plt.xlim(0,1.0)
	plt.ylim(0,1.0)
	plt.show()
	#plot PR Curve
	rs, ps = zip(*sorted(zip(recalls, precisions)))
	plt.plot(rs, ps, 'o--')
	plt.plot([0,1.0],[1.0,0], '--')
	plt.title("PR Curve")
	plt.xlabel("Recall")
	plt.ylabel("Precision")
	plt.ylim(0,1.0)
	plt.xlim(0,1.0)
	plt.show()
	#making 3d plot with thresholds/ROC values and 3d plot with thresholds/PR values
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.scatter(falsePositiveRates, truePositiveRates, thresholds)
	ax.set_xlabel('FPR')
	ax.set_ylabel('TPR')
	ax.set_zlabel('Thresholds')
	ax.set_xlim3d(0,1.0)
	ax.set_ylim3d(0,1.0)
	ax.set_zlim3d(0,1.0)
	plt.show()
					
def generateGraphCMS(names, location, pasType = "TE", bufferVal = 1, spacing = 50, distance = 1, tolerance = 0):
	cmVals  = {} #key will be threshold values
	for name in names:
		#for each chromsome open stored CM information 
		chromosomeCMValues = openCMDataOneChromosome(location, name)
		foundMinHeights = []
		#make Confusion Matrices using peakConfusionMatrix from extractingNegatives.py
		for index, row in chromosomeCMValues.iterrows():
			if row['pasType'] == pasType and row['bufferVal'] == bufferVal and row['spacing'] == spacing and row['distance'] == distance and row['tolerance'] == tolerance:
				keyVal = round(row['minh'], 8)
				if keyVal not in cmVals:
					cmVals[keyVal] = [row['TruePositives'], row['FalsePositives'], row['FalseNegatives'], row['TrueNegatives'], row['minh']]
					foundMinHeights.append(keyVal)
				else:
					if keyVal not in foundMinHeights: 
						#update values 
						cmVals[keyVal][0] += row['TruePositives']
						cmVals[keyVal][1] += row['FalsePositives']
						cmVals[keyVal][2] += row['FalseNegatives']
						cmVals[keyVal][3] += row['TrueNegatives']
						foundMinHeights.append(keyVal)
	#now change all threshold values into confusion matrices
	confusionMatrices = []
	for key in cmVals:
		confusionMatrices.append(peakConfusionMatrix(cmVals[key][0], cmVals[key][1], cmVals[key][2], cmVals[key][3], key))
	return confusionMatrices
	

def graphROCMultipleDistances(names, location, pasType = "TE", bufferVal = 1, spacing = 50, distances = [1,25,50], tolerance = 0):
	for d in distances:
		confusionMatrices = generateGraphCMS(names, location, pasType, bufferVal, spacing, d, tolerance)
		falsePositiveRates = [c.falsePositiveRate() for c in confusionMatrices] #x axis ROC curve
		truePositiveRates = [c.truePositiveRate() for c in confusionMatrices] #y axis ROC Curve
		#graph ROC curve 
		fprs, tprs = zip(*sorted(zip(falsePositiveRates, truePositiveRates)))
		plt.plot(fprs, tprs, 'o--', label = 'distance ' + str(d))
	plt.plot([0,1.0],[0,1.0], '--')
	plt.title(pasType + " ROC")
	plt.xlabel("FPR")
	plt.ylabel("TPR")
	plt.legend()
	plt.xlim(0,1.0)
	plt.ylim(0,1.0)
	#print (thresholds)
	plt.show()
	
	
def graphPRMultipleDistances(names, location, pasType = "TE", bufferVal = 1, spacing = 50, distances = [1,25,50], tolerance = 0):
	for d in distances:
		confusionMatrices = generateGraphCMS(names, location, pasType, bufferVal, spacing, d, tolerance)
		recalls = [c.recall() for c in confusionMatrices] #x axis on Precision Recall Curve
		precisions = [c.precision() for c in confusionMatrices] #y axis on Precision Recall Curve
		rs, ps = zip(*sorted(zip(recalls, precisions)))
		plt.plot(rs, ps, 'o--', label = 'distance ' + str(d))
	plt.plot([0,1.0],[1.0,0], '--')
	plt.title(pasType + " PR Curve")
	plt.xlabel("Recall")
	plt.ylabel("Precision")
	plt.legend()
	plt.ylim(0,1.0)
	plt.xlim(0,1.0)
	plt.show()

	
def graphThresholdsMultipleDistances(names, location, pasType = "TE", bufferVal = 1, spacing = 50, distances = [1,25,50], tolerance = 0):
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	for d in distances:
		confusionMatrices = generateGraphCMS(names, location, pasType, bufferVal, spacing, d, tolerance)
		falsePositiveRates = [c.falsePositiveRate() for c in confusionMatrices] #x axis ROC curve
		truePositiveRates = [c.truePositiveRate() for c in confusionMatrices] #y axis ROC Curve
		thresholds = [c.threshold for c in confusionMatrices]
		#making 3d plot with thresholds/ROC values and 3d plot with thresholds/PR values
		ax.scatter(falsePositiveRates, truePositiveRates, thresholds, label = 'distance ' + str(d))
	ax.set_xlabel('FPR')
	ax.set_ylabel('TPR')
	ax.set_zlabel('Thresholds')
	ax.set_xlim3d(0,1.0)
	ax.set_ylim3d(0,1.0)
	ax.set_zlim3d(0,1.0)
	plt.legend()
	plt.show()

def condenseCSV(names, location):
	#condenses all values with same pasType/bufferVal/spacing/minh/distance/tolerance together into one row
	compiledData = pd.DataFrame({'name': [], 'pasType': [], 'bufferVal':[], 'spacing':[], 'minh':[], 'distance':[], 'tolerance':[], 'TruePositives':[], 'FalsePositives':[], 'FalseNegatives':[], 'TrueNegatives':[]})
	for name in names:
		resultsForName = openCMDataOneChromosome(location, name)
		#index, row in balancedPos.iterrows()
		for i, row in resultsForName.iterrows():	
			#filteredCompileData = compiledData.index[compiledData['pasType'] == row['pasType'], compiledData['bufferVal'] == row['bufferVal'], compiledData['spacing'] == row['spacing'], compiledData['minh'] == row['minh'], compiledData['distance'] == row['distance'], compiledData['tolerance'] == row['tolerance']]
			if not filteredCompileData: #that combination is not in the combined dataframe yet
				print ("needs to be added")
			else: #it is in the combine dataframe and should be updated
				print ("is in the dataframe")
				



location = "./ConfusionMatrices/"
names = ["1", "Y","2", "3","4", "5", "6","7","8","9","10","11","12","13", "14", "15", "16", "17", "18", "19", "20", "21", "22"]
'''
#TE pastype
graphROCMultipleDistances(names, location, pasType = "TE", bufferVal = 1, spacing = 50, distances = [1,25,50], tolerance = 0)
graphPRMultipleDistances(names, location, pasType = "TE", bufferVal = 1, spacing = 50, distances = [1,25,50], tolerance = 0)
graphThresholdsMultipleDistances(names, location, pasType = "TE", bufferVal = 1, spacing = 50, distances = [1,25,50], tolerance = 0)
'''
#IN pastype
graphROCMultipleDistances(names, location, pasType = "IN", bufferVal = 1, spacing = 50, distances = [1,25,50], tolerance = 0)
graphPRMultipleDistances(names, location, pasType = "IN", bufferVal = 1, spacing = 50, distances = [1,25,50], tolerance = 0)
graphThresholdsMultipleDistances(names, location, pasType = "IN", bufferVal = 1, spacing = 50, distances = [1,25,50], tolerance = 0)
