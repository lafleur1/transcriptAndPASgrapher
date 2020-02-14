#extracting sequences near the true values in same regions for true negative datasets for APA model evaluation
#APARENT: 206 nt, overlapping
#DeepPASTA: 200 nts, polyA site in middle
#
#using same method as Leung et al (see  Inference of the human polyadenylation code supplementary material):
 #same padding requirement (i.e. that spacing between two PAS signals must allow four negative regions 
 #shift left and right by 50 nts but still within the region boundaries 
 
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

#Note: Ensembl is 1-based, polyAsite is 0-based

#rules used for polyAsite2.0 cluster classification: 
'''
if PAS in exon at the end of any transcript:
	TE
elif PAS in any exon in a transcript:
	EX
elif PAS in any intron
	IN
elif PAS within 1000 ds of any transcript terminal exon
	DS
elif antisense to an exon
	AE
elif antisense to an intron 
	AI
elif 1000 nt antisense upstream of a start site ??? (there are some errors with this category (potentially)
	AU #some errors here, possibly.  Will need to check this category to see if they are faulty or not
	#may have off by one error???
else:
	IG
	
'''
#PolyASite2.0 created with version 96 of Ensembl release
ensembl = pyensembl.EnsemblRelease(release = '96')


######################################################################MAKING NEGATIVE DATASETS #######################################################################

#opens human cluster data for given chromosome.  Can optionally filter by PAS type
def openPASClustersForChromosome(name, pasType = 'All'):
	#opening all the true values from PolyASite2.0 
	colnames = ["seqName",  "start" , "end",  "clusterID",  "avgTPM",  "strand",   "percentSupporting",   "protocolsSupporting",  "avgTPM2",   "type",   "upstreamClusters"]
	pas_stuff =pd.read_csv('atlas.clusters.hg38.2-0.bed',delimiter='\t', names = colnames, dtype = {"seqName": str}) 
	trueValBoolMask = pas_stuff['seqName'] == name
	currentTrueVals = pas_stuff[trueValBoolMask] #filtered true vals
	if pasType == "All":
		return currentTrueVals
	else:
		maskType = currentTrueVals["type"] == pasType
		maskedTrue = currentTrueVals[maskType]
		return maskedTrue


#creates false dataset points per true dataset point
#using the DeepPASTA/established APA DL method for creating negative datapoints by shifting the true valye +/- 50 nts and declaring that a 'false' APA site
#fails the pos signal if there are "N" in the true region, or in the region used to predict the true pas region [start - 205, end + 205]
#inputs: current row in the polyaSite2.0 dataframe, all the other PAS clusters in the polyAsite2.0 dataframe, all current negatives, current balanced negatives, the value used for spacing (how many of the cluster must fit between shifted groups), how much the cluster is shifted to create the false cluster position, the fasta sequence for the current chromosome (to look for N's in the input sequence for the network), strictNFiltering (true if sequences should be removed if there are N's in the input sequence used in APARENT), strictFilteringValue (used to find range to check for N's))
#outputs:  leftPassedAll (if the left negative passed and shoutd be added to the overall dataset), rightPassedAll (if the right negative passed and should be added to the overall dataset), leftPassedBalanced (if the left negative passed and can be added to the balanced dataset), rightPassedBalanced (if the right negative passed and can be added to the balanced dataset), leftStart, leftEnd, rightStart, rightEnd (start/end values for the clusters)
def shiftLeftAndRightNegativesFilterPositives(pasRow, allPAS, currentAllNegs, currentBalancedNegs, spacingValue, shiftValue, fastaSeq, strictNFiltering = True, strictFilteringValue = 205):
	#following Leung et al's true negative rules, count it as a negative if it can be shifted 50 nts and still fit four multiples of the negative region between it and the next polyA signal
	#filtering positive signals for sequences with N's
	
	sliceStart = pasRow['start'] - strictFilteringValue
	if sliceStart < 0:
		sliceStart = 0 #if the PAS occurs at the very beginning of the sequence correct for this
	sliceEnd = pasRow['end'] + strictFilteringValue
	if sliceEnd > len(fastaSeq) -1:
		sliceEnd = len(fastaSeq) -1 
	filteringTrueSequence = fastaSeq[sliceStart:sliceEnd + 1] #remove string portion to test for true
	countN = filteringTrueSequence.count("N")
	if countN == 0:  
		width = pasRow['end'] - pasRow['start']
		checkRange = (width * spacingValue) + shiftValue
		#  SET UP
		leftStart = pasRow['start'] - shiftValue
		leftEnd = pasRow['end'] - shiftValue
		leftSliceStart = leftStart - strictFilteringValue
		leftSliceEnd = leftEnd + strictFilteringValue
		if leftSliceStart < 0 :
			leftSliceStart = 0
		if leftSliceEnd > len(fastaSeq) - 1:
			leftSliceEnd = len(fastaSeq) - 1
		leftSlice = fastaSeq[leftSliceStart: leftSliceEnd + 1]
		leftNCount = leftSlice.count("N")
		leftPassedAll = False
		leftPassedBalanced = False
		
		rightStart = pasRow['start'] + shiftValue
		rightEnd = pasRow['end'] + shiftValue
		rightSliceStart = leftStart - strictFilteringValue
		rightSliceEnd = leftEnd + strictFilteringValue
		if rightSliceStart < 0 :
			leftSliceStart = 0
		if rightSliceEnd > len(fastaSeq) - 1:
			leftSliceEnd = len(fastaSeq) - 1
		rightSlice = fastaSeq[rightSliceStart:rightSliceEnd + 1]
		rightNCount = rightSlice.count("N")
		rightPassedAll = False
		rightPassedBalanced = False
		### LEFT
		if leftNCount == 0: 
			#check if a true positive will occur in this range
			leftCheck = pasRow['start'] - checkRange
			filterVals1 = (allPAS['strand'] == pasRow['strand']) & (allPAS['seqName'] == pasRow['seqName']) & (allPAS['end'] >= leftCheck) & (allPAS['start'] >= leftCheck) & (allPAS['start'] < pasRow['start'])
			filterVals2 = (allPAS['strand'] == pasRow['strand']) & (allPAS['seqName'] == pasRow['seqName']) & (allPAS['end'] >= leftCheck) & (allPAS['start'] <= leftCheck)
			otherPAS1Left = allPAS[filterVals1]
			otherPAS2Left = allPAS[filterVals2]
			sumInForbiddenRange = otherPAS1Left.shape[0] + otherPAS2Left.shape[0]
			if sumInForbiddenRange == 0: 
				negativeAllFilter = ((currentAllNegs['start'] <= leftStart) & (currentAllNegs['end'] >= leftStart)) | ((currentAllNegs['start'] <= leftEnd) & (currentAllNegs['end'] >= leftEnd)) | ((currentAllNegs['start'] >= leftStart) & (currentAllNegs['start'] <= leftEnd))				
				negativeBalancedFilter = ((currentBalancedNegs['start'] <= leftStart) & (currentBalancedNegs['end'] >= leftStart)) | ((currentBalancedNegs['start'] <= leftEnd) & (currentBalancedNegs['end'] >= leftEnd)) | ((currentBalancedNegs['start'] >= leftStart) & (currentBalancedNegs['start'] <= leftEnd))	
				negativeAll = currentAllNegs[negativeAllFilter]
				negativeBalanced = currentBalancedNegs[negativeBalancedFilter]
				if negativeAll.shape[0] == 0:
					leftPassedAll = True
				if negativeBalanced.shape[0] == 0:
					leftPassedBalanced = True
		else:
			print ("Sequence rejected on left negative: ", leftSlice)
			leftStart = -1
			leftEnd = -1
		#####
		##### RIGHT 
		if rightNCount ==0:
			#check if the negative range contains/overlaps with an existing positive
			rightCheck = pasRow['end'] + checkRange
			filterVals1 = (allPAS['strand'] == pasRow['strand']) & (allPAS['seqName'] == pasRow['seqName']) & (allPAS['start'] <= rightCheck) & (allPAS['end'] <= rightCheck) & (allPAS['start'] >= pasRow['end']) 
			filterVals2 = (allPAS['strand'] == pasRow['strand']) & (allPAS['seqName'] == pasRow['seqName']) & (allPAS['start'] <= rightCheck) & (allPAS['end'] >= rightCheck)
			otherPAS1Right = allPAS[filterVals1]
			otherPAS2Right = allPAS[filterVals2]
			sumInForbiddenRangeRight = otherPAS1Right.shape[0] + otherPAS2Right.shape[0]
			if sumInForbiddenRangeRight == 0:
				negativeAllFilter = ((currentAllNegs['start'] <= leftStart) & (currentAllNegs['end'] >= leftStart)) | ((currentAllNegs['start'] <= leftEnd) & (currentAllNegs['end'] >= leftEnd)) | ((currentAllNegs['start'] >= leftStart) & (currentAllNegs['start'] <= leftEnd))				
				negativeBalancedFilter = ((currentBalancedNegs['start'] <= leftStart) & (currentBalancedNegs['end'] >= leftStart)) | ((currentBalancedNegs['start'] <= leftEnd) & (currentBalancedNegs['end'] >= leftEnd)) | ((currentBalancedNegs['start'] >= leftStart) & (currentBalancedNegs['start'] <= leftEnd))	
				negativeAll = currentAllNegs[negativeAllFilter]
				negativeBalanced = currentBalancedNegs[negativeBalancedFilter]
				if negativeAll.shape[0] == 0:
					rightPassedAll = True
				if negativeBalanced.shape[0] == 0:
					rightPassedBalanced = True
		else:
			print ("Sequence rejected on right negative: ", rightSlice)
			rightStart = -1
			rightEnd = -1
		return leftPassedAll, rightPassedAll, leftPassedBalanced, rightPassedBalanced, leftStart, leftEnd, rightStart, rightEnd
	else:
		print ("Sequence rejected on positive region: ", filteringTrueSequence)
		return False, False, False, False, -2, -2, -2, -2
	


#create the overall negative dataset and the balanced dataset for one chromosome
#inputs: trueVals (polyAsite2.0 cluster dataframe for the chromosome), spacingValue (how many copies of the cluster to space between the created negative and any other positive clusters), shiftValue (how many nts to shift the positive cluster to make a negative cluster), fileName (name to save the datasets under), fastaSeq (overall chromosome sequence to use to remove sequences with N's)
#outputs: saves a report file under ./reports/fileNameDatasetsReport.txt, saves datasets under ./datasets/fileNameBalancedPositives.csv, ./datasets/fileNameAllNegatives.csv,
#./datasets/fileNameBalancedNegatives.csv
def createNegativeDataSet(trueVals, spacingValue, shiftValue, fileName, fastaSeq):
	#create dictionary for negative valeus left and right 
	#["seqName",  "start" , "end",  "clusterID",  "avgTPM",  "strand",   "percentSupporting",   "protocolsSupporting",  "avgTPM2",   "type",   "upstreamClusters"]
	blankDict = {"seqName":[], "start":[], "end":[], "clusterID":[], "strand":[], "type":[], "side":[]}
	balancedDict = {"seqName":[], "start":[], "end":[], "clusterID":[], "strand":[], "type":[], "side":[]}
	droppedTypeList = []
	copyTrueValues = trueVals.copy(deep = True) #make deep copy of dataframe of the true pas values to filter
	random.seed()
	#setup the negative datasets
	negativesAllDF = pd.DataFrame(blankDict)
	negativesBalancedDF = pd.DataFrame(balancedDict)
	for index, row in trueVals.iterrows():
		leftPassedAll, rightPassedAll, leftPassedBalanced, rightPassedBalanced, leftStart, leftEnd, rightStart, rightEnd = shiftLeftAndRightNegativesFilterPositives(row, trueVals, negativesAllDF, negativesBalancedDF, spacingValue, shiftValue, fastaSeq) #determine which datasets to add the left/right negative clusters to
		leftRow = {"seqName":[], "start":[], "end":[], "clusterID":[], "strand":[], "type":[], "side":[]}
		rightRow = {"seqName":[], "start":[], "end":[], "clusterID":[], "strand":[], "type":[], "side":[]}
		#fill left row
		leftRow['seqName'].append(row['seqName'])
		leftRow['start'].append(leftStart)
		leftRow['end'].append(leftEnd)
		newID = row['clusterID'] + "_LeftNegative"
		leftRow['clusterID'].append(newID)
		leftRow['strand'].append(row['strand'])
		leftRow['type'].append(row['type'])
		leftRow['side'].append("Left")
		#fill right row
		rightRow['seqName'].append(row['seqName'])
		rightRow['start'].append(rightStart)
		rightRow['end'].append(rightEnd)
		newID = row['clusterID'] + "_RightNegative"
		rightRow['clusterID'].append(newID)
		rightRow['strand'].append(row['strand'])
		rightRow['type'].append(row['type'])
		rightRow['side'].append("Right")
		#print if the left/right clusters fail due to N's in the true cluster input ranges
		if leftStart == leftEnd == -1: 
			print ("Failed due to N's in left sequence")
		if rightStart == rightEnd == -1:
			print ("Failed due to N's in right sequence")
		if leftStart == leftEnd == rightStart == rightEnd == -2:
			print ("Failed due to N's in true positive sequence")
		#prepare rows to be added to the overall datasets
		leftRowDF = pd.DataFrame(leftRow)
		rightRowDF = pd.DataFrame(rightRow)
		if leftPassedAll and rightPassedAll:
			#add both to All dataset
			negativesAllDF = negativesAllDF.append(leftRowDF)
			negativesAllDF = negativesAllDF.append(rightRowDF)
		elif leftPassedAll and not rightPassedAll:
			#add left to dataset
			negativesAllDF = negativesAllDF.append(leftRowDF)
		elif not leftPassedAll and rightPassedAll:
			#add right to dataset
			negativesAllDF = negativesAllDF.append(rightRowDF)
		#balanced dataset 
		if leftPassedBalanced and rightPassedBalanced:
			#if both passed, flip coin and add one to balanced dataset
			randInt = random.randint(0,1)
			if randInt == 0: 
				#add left to balanced
				negativesBalancedDF = negativesBalancedDF.append(leftRowDF)
			else:
				#add right to balanced
				negativesBalancedDF = negativesBalancedDF.append(rightRowDF)
		elif leftPassedBalanced and not rightPassedBalanced:
			#if only left, add left to balanced dataset
			negativesBalancedDF = negativesBalancedDF.append(leftRowDF)
		elif not leftPassedBalanced and rightPassedBalanced:
			#if only right, add right to balanced dataset	
			negativesBalancedDF = negativesBalancedDF.append(rightRowDF)
		if not leftPassedBalanced and not rightPassedBalanced: 
			#delete the row from the copy of the true values
			copyTrueValues.drop(index = index, inplace = True)
			print ("DROPPED: ", row['clusterID'], " ", row['type'])
			droppedTypeList.append(row['type'])
	###
	#create a report of the types dropped for each chromosome
	report = open("./reports/" + fileName + "DatasetsReport.txt", "w")
	pasTypes = ['All', 'IN', 'TE', 'IG', 'AI', 'EX', 'DS', 'AE', 'AU'] 
	report.write("Positives dropped: " + "\n")
	for t in pasTypes:
		print ("dropped ", droppedTypeList.count(t), " of type ", t)
		report.write("dropped " + str(droppedTypeList.count(t)) + " of type " + t + "\n")
	report.close()
	####
	#save the datasets with to_csv for each chromosome 
	filteredTrueName = "./datasets/" + fileName + "BalancedPositives.csv"
	allNegatives = "./datasets/" +  fileName + "AllNegatives.csv"
	balancedNegatives ="./datasets/" +  fileName + "BalancedNegatives.csv"
	copyTrueValues.to_csv(filteredTrueName)
	negativesAllDF.to_csv(allNegatives)
	negativesBalancedDF.to_csv(balancedNegatives)
	

#opens forward and reverse complement APARENT predictions for a chromosome
#input: stem (file path the predictions are stored in), name (chromosome name)
#output: forward strand predictions, reverse complement strand predictions
def openForwardReverse(stem, name):
	totalNameFor = stem + name + ".npy"
	print (totalNameFor)
	forward = np.load(totalNameFor)
	reverse = np.load(stem + name + "RC.npy")
	return forward, reverse


#create balanced positive/negative datasets for one chromosome 
#input: chroName (name of the chromosome the dataset is curently being generated for), bufferRegions (list of spacings to make datasets for), spacings (how much to shift each positive cluster to make a negative cluster), fastaPath (path where all fasta files are stored on the computer, uses these to remove any sequences with N's in the APARENT input range around a cluster)
def createBalancedDatasets(chroName, bufferRegions, spacings, fastaPath):
	#open fasta
	totalName = fastaPath + "chr" + chroName + ".fasta"
	chroSeq = SeqIO.read(totalName, "fasta")
	print ("sequence length: ", len(chroSeq.seq))
	currentTime = time.time()
	print ("___________________________________")
	print ("FASTA OPENED: ")
	print (chroSeq)
	print ("___________________________________")
	#open PAS
	openPAS = openPASClustersForChromosome(chroName, pasType= "All") #making a csv which can be filtered for individual PAS types as needed
	for s in spacings:
		for b in bufferRegions:
			print ("___________________________________")
			currentTime = time.time()
			fileName = "chro" + chroName + "_NegSpaces" + str(b) + "_shifted" + str(s) + "Nts"
			print ("generating balanced datasets for: ")
			print (fileName)
			createNegativeDataSet(openPAS, b, s, fileName, chroSeq.seq)
			print ("___________________________________")
			print ("elapsed time: ", time.time() - currentTime)
			print ("___________________________________")
			print (" ")
			
#makes the balanced datasets for the entire genome.  
#RUN THIS TO CREATE ALL BALANCED DATASETS FROM THE POLYASITE2.0 DATASET	
def makeNegativeDatasetsForGenome():
	fastaPath = "../../aparentGenomeTesting/fastas/"
	bufferRegions = [1,4]
	spacing = [ 50]
	chromosomes = ["Y", "1","2","3","4","5","6","7","8","9","10","11","12","13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"]
	for c in chromosomes:
		createBalancedDatasets(c, bufferRegions, spacing, fastaPath)


def makeNegativeDatasetsForOneChromosome(chroName):
	fastaPath = "../../aparentGenomeTesting/fastas/"
	bufferRegions = [1,4]
	spacing = [ 50, 75, 100]
	createBalancedDatasets([chroName], bufferRegions, spacing, fastaPath)

###############################################################################################################################################################################

############################################################# APARENT ASSESSMENT #################################################################################

############################# Scipy metrics ROC and PR Curves ##########################
#docs for the curve functions at: 
#https://scikit-learn.org/stable/modules/generated/sklearn.metrics.roc_curve.html 
#https://scikit-learn.org/stable/modules/generated/sklearn.metrics.precision_recall_curve.html 
#scores:
#https://scikit-learn.org/stable/modules/generated/sklearn.metrics.average_precision_score.html#sklearn.metrics.average_precision_score 
#https://scikit-learn.org/stable/modules/generated/sklearn.metrics.roc_auc_score.html#sklearn.metrics.roc_auc_score

#functions to open created datasets:
#################
def openAllNegatives(name):
	negName = name + "AllNegatives.csv"
	return pd.read_csv( negativesName,  dtype = {"seqName": str}) 

def openBalancedNegatives(name):
	negativesName = name + "BalancedNegatives.csv"
	return pd.read_csv( negativesName, dtype = {"seqName": str}) 

def openBalancedPositives(name):
	positivesName = name + "BalancedPositives.csv"
	return pd.read_csv( positivesName, dtype = {"seqName": str}) 
###############	

#for extracting the values from a prediction numpy for a given chromosome
#input: name (chromosome name to open values for and extract true/false cluster range values from), negatives (dataset to use to get values from numpy array), positives (dataset to use to get values from numpy array), pasType (does all pasTypes if not set.) 
#outputs: dummybools, dummy #return all positve/negative values in one array, boolean labels for those values in other array
def extractPredictionValues(name, negatives, positives, pasType = ""):
	#making a sklearn AUC and AUPRC plot to look at different threshold values
	dummy = np.array([])
	dummybools = np.array([])
	predName = "chr" + name
	forward, reverse = openForwardReverse("../../aparentGenomeTesting/chromosomePredictions50/", predName) #open APARENT prediction values for the current chromosome
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
def extractAllPredictionValues(names, stem, b, s, pasType = ""):
	values = np.array([])
	bools = np.array([])
	for name in names:
		fileName = stem + "chro" + name + "_NegSpaces" + str(b) + "_shifted" + str(s) + "Nts"
		negatives = openBalancedNegatives(fileName)
		positives = openBalancedPositives(fileName)
		boolsCurrent, valsCurrent = extractPredictionValues(name, negatives, positives, pasType) #use above function to open the predicted APARENT values for the chromosome clusters 
		values = np.concatenate((values, valsCurrent))
		bools = np.concatenate((bools, boolsCurrent))
	return values, bools
	

#extract the values for the positive/negative clusters separately to look at distributions
#inputs: name (chromsome name to extract values for), negatives (balanced dataset clusters to extract predictions for), positives (balance dataset clusters dataframe to extract predictions for), pasType (defaults to all the types)
#outputs: dummyNegatives (all the negative prediction values), dummyPositives (all the positive prediction values)
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
					sliceVals = forward[startInt:endInt+ 1]
					dummyNegatives = np.concatenate((dummyNegatives,sliceVals))
		elif row['strand'] == "-":
			if pasType == "":
				#reverse strand
				sliceVals = flippedReverse[startInt: endInt + 1]
				dummyNegatives = np.concatenate((dummyNegatives,sliceVals))
			else:
				if row['type'] == pasType:
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

#functions to find the best thresholds to use for a classifier for a AUC ROC plot
############################################################3
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
#############################################3

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
	plt.title("AUC All But X IN")
	plt.show()
	plt.plot(prec, rec)
	plt.title("PR Curve All but X IN")
	plt.show()
	

#########################################################################################################################################


################## CREATING CONFUSION MATRICS FOR PEAKS IN THE DATA #####################################################

'''
Since we cannot use the scipy roc_curve or precision_recall_curve functions, we will have to manually extract peaks and calculate the false/true positives and negatives using the following criteria:
TP:  If a true polyASite2.0 cluster has a predicted value peak in it
FP: if a false cluster has a predicted value peak in it
TN: if a false cluster does not have a predicted value peak in it
FN: if a true cluster does not have a predicted value peak in it
'''

#confusion matrix class to use when calculating the ROC and PR curves for APARNET
class peakConfusionMatrix(object):
	def __init__(self, countTP, countFP, countFN, countTN , unplacedPredictions = 0, totalPredictions = 0):
		self.truePositives = countTP
		self.falsePositives = countFP
		self.falseNegatives = countFN
		self.trueNegatives = countTN
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
			return 0
	
	def recall(self):
		#Same as sensitvity AKA TPR
		return self.truePositiveRate()
	
	#return fraction of peaks not falling in a true/false cluster across the entire genome
	def fractionUnplaced(self):
		if self.totalPredictionsMade != 0:
			return self.unplacedPredictions/self.totalPredictionsMade
		else:
			return 0 

		
#use to correct the indexes for the RC predictions so that they can be matched with the true cluster labels (since those are indexed according to the forward strand 5-> 3 direction)
def flipSequenceIndex(index, lenSeq):
	b = lenSeq - 1
	return -1 * index + b


def find_peaks_ChromosomeVersion(avgPreds, peak_min_height, peak_min_distance, peak_prominence = None):
	peaks, _ = find_peaks(avgPreds, height=peak_min_height, distance=peak_min_distance, prominence=peak_prominence) 
	return peaks


################################################################# DICTIONARY BASED METHODS #######################################################


#creates forward and reverse strand dictionaries for the balanced dataset given to it
#inputs: toSeparate (dataframe of true or false clusters to break into dictionary of clusters), pasType (pasType to place in dictionary)
#outputs: clustersForward (dictionary of clusters on + strand), clustersRC (dictionary of clusters on - strand) #key is [start, end] for the cluster, contents is empty list
def openBalancedValuesForType(toSeparate, pasType = ""):
	clustersForward = {}
	clustersRC = {} #key of (Start, End) will use to track clusters with no peaks for the FN  
	if pasType == "":
		for index, row in toSeparate.iterrows():
			if row['strand'] == "+": #forward strand
				clustersForward[(row['start'], row['end'])] = []
			else: #negative strand, RC cluster
				clustersRC[(row['start'], row['end'])] = []
		clustersForward['Unplaced'] = []
		clustersRC['Unplaced'] = []
	else:
		#dont bother with unplaced peaks if not doing all pasTypes if not looking for all types
		maskType = toSeparate["type"] == pasType
		maskedTrue = toSeparate[maskType]
		for index, row in maskedTrue.iterrows():
			if row['strand'] == "+": #forward strand
				clustersForward[(row['start'], row['end'])] = []
			else: #negative strand, RC cluster
				clustersRC[(row['start'], row['end'])] = []	
	return clustersForward, clustersRC

#place peaks from APARENT into the cluster dictionaries
#input: peaks (list of peaks to place in the dictionaries), clustersNegative (negative clusters dictionary), clustersPositive (positive clusters dictionary), tolerance (tolerance to place a peak in the cluster), sign (+/- strand ), lenSeq (lenght of the sequence peaks were extracted from )
#output: clustersNegative (filled negative dictionary), clustersPositive (filled positive dictionary), numberUnplacedPeaks (number peaks not placed in a dictionary)
def placePeaksWithTolerance(peaks, clustersNegative, clustersPositive, tolerance, sign, lenSeq):
	clusterRangesNegative = clustersNegative.keys()
	clusterRangesPositive = clustersPositive.keys()
	numberUnplacedPeaks = 0
	for peak in peaks:
		if sign == "-":
			peak = flipSequenceIndex(peak, lenSeq)
		placed = False
		for rng in clusterRangesPositive:
			if rng != 'Unplaced':
				lower = rng[0] - tolerance
				upper = rng[1] + tolerance
				if peak >= lower and peak <= upper: #if the peak is in [start, end]
					clustersPositive[rng].append(peak)
					placed = True
					break
		if not placed:
			for rng in clusterRangesNegative:
				if rng != 'Unplaced':
					lower = rng[0] - tolerance
					upper = rng[1] + tolerance
					if peak >= lower and peak <= upper: #if the peak is in [start, end]
						clustersNegative[rng].append(peak)
						placed = True
						break			
		if not placed: #wasn't placed
			numberUnplacedPeaks += 1
	return clustersNegative, clustersPositive, numberUnplacedPeaks


#using the peak-cluster dictionaries fill a confidence matrix for a given peak threshold
#need to count TP, FP, FN (placed peaks, unplaced peaks, empty clusters)
#inputs: dictForwardPlus (filled forward strand true cluster dictionary), dictForwardNegative (filled negative cluster forward strand dictionary), 
#dictRCPlus (filled reverse complement true cluster dictionary), dictRCNegative (filled reverse complement false cluster dictionary), countOncePerCluster = True
#outputs: countTP (count of true positives), countFP (count of false positives), countFN (count of false negatives), countTN (count of true negatives)
def fillConfMatrix(dictForwardPlus, dictForwardNegative, dictRCPlus, dictRCNegative, countOncePerCluster = True):
	countTP = 0 #peaks in true cluster
	countFP = 0 #peak in false cluster 
	countFN = 0 #no peak in true cluster
	countTN = 0 #no peak in false cluster
	for key in dictForwardPlus:
		inCluster = len(dictForwardPlus[key])
		if inCluster != 0: #peak in a positive cluster
			countTP += 1
		else: #empty true cluster
			countFN += 1
	for key in dictRCPlus:
		inCluster = len(dictRCPlus[key])
		if inCluster != 0:
			countTP += 1 #peak in a positive cluster
		else:
			countFN += 1 #empty true cluster
	for key in dictForwardNegative:
		inCluster = len(dictForwardNegative[key])
		if inCluster != 0:
			#peak in false cluster
			countFP += 1
		else:
			#no peak in false cluster
			countTN += 1
	for key in dictRCNegative:
		inCluster = len(dictRCNegative[key])
		if inCluster != 0:
			#peak in false cluster
			countFP += 1
		else:
			#no peak in false cluster
			countTN += 1
	return countTP, countFP, countFN, countTN
			
#fill a confidence matrix for one chromosome
#input: name (chromosome name), stem (file path to the balanced dataset), b (spacing used to create dataset), s (shifted value used to create negative dataset), minh (minimum peak height for find_peaks), dist (distance between peaks for find_peaks), tolerance (tolerance for placing peaks into the clusters), pasType (pasType to create confidence matrix for) 
#output: countTP, countFP, countFN, countTN, number unplaced peaks, total peaks
def buildConfidenceMatrixOneChro(name, stem, b, s, minh, dist, tolerance, pasType = ""):
	fileName = stem + "chro" + name + "_NegSpaces" + str(b) + "_shifted" + str(s) + "Nts"
	print ("Opening: ", fileName)
	balancedPositives = openBalancedPositives(fileName)
	balancedNegatives = openBalancedNegatives(fileName)
	print ("Size balanced datasets: ", balancedPositives.shape[0])
	clustersFPositve, clustersRCPositive = openBalancedValuesForType(balancedPositives, pasType)
	clustersFNegative, clustersRCNegative = openBalancedValuesForType(balancedNegatives, pasType)
	#open peaks for the forward and reverse strand predictions
	predName = "chr" + name
	forward, reverse = openForwardReverse("../../aparentGenomeTesting/chromosomePredictions50/", predName)
	forwardPeaks = find_peaks_ChromosomeVersion(forward, minh, dist, (0.01, None)) 
	reversePeaks = find_peaks_ChromosomeVersion(reverse, minh, dist, (0.01, None)) 
	print ("Peaks found: ", len(forwardPeaks), " ", len(reversePeaks))
	#place peaks
	clustersFNegative, clustersFPositve, numberUnplacedForward = placePeaksWithTolerance(forwardPeaks, clustersFNegative, clustersFPositve, tolerance, "+", forward.shape[0])
	print ("Placed + strand peaks")
	clustersRCNegative, clustersRCPositive, numberUnplacedReverse = placePeaksWithTolerance(reversePeaks, clustersRCNegative, clustersRCPositive, tolerance, "-", forward.shape[0])
	print ("Placed - strand peaks")
	countTP, countFP, countFN, countTN = fillConfMatrix(clustersFPositve, clustersFNegative, clustersRCPositive, clustersRCNegative)
	return countTP, countFP, countFN, countTN, numberUnplacedForward + numberUnplacedReverse, len(forwardPeaks) + len(reversePeaks)
	

#build confidence matrices for multiple chromosomes 
#input: names (list of chromosome names), stem (file path), bufferVal (buffer value to use to make confidence matrix), spacing (how much clusters are shifted to make negative datasets), minh (minimum height to use for find_peaks), dist (distance between peaks for find_peaks), tolerance (tolerance to use for placing peaks in clusters), pasType (pasType to calculate confidence matrix for)
#output: peakConfusionMatrix object
def buildConfidenceMatrixMultipleChromosomes(names, stem, bufferVal, spacing, minh, dist, tolerance, pasType = ""):
	countTPOverall = 0
	countFPOverall = 0
	countFNOverall = 0
	countTNOverall = 0
	countUnplaced = 0
	totalPeaks = 0
	for name in names: #for each chromosome add values to overall count for making confidence matrix
		print ("On chromosome: ", name)
		countTP, countFP, countFN, countTN, unplaced, numberPeaks = buildConfidenceMatrixOneChro(name, stem, bufferVal, spacing, minh, dist, tolerance, pasType)
		print ("TP:", "FP: ", "FN: ", "TN: ", "unplaced: ", "total: ")
		countTPOverall += countTP
		countFPOverall += countFP
		countFNOverall += countFN
		countTNOverall += countTN
		countUnplaced += unplaced
		totalPeaks += numberPeaks
		print ("True Positives: ", countTPOverall, "False Positives: ", countFPOverall, "False Negatives: ", countFNOverall, "True Negatives: ", countTNOverall, "Unplaced: ",countUnplaced,  "Total Peaks: ", totalPeaks)
	return peakConfusionMatrix(countTPOverall, countFPOverall, countFNOverall, countTNOverall, countUnplaced, totalPeaks)
	
	

def buildConfusionMatricesForGraphing(names, stem, bufferVal, spacing, minhs, dist, tolerance, rocTitle, prTitle, pasType = ""):
	confusionMatrices = []
	for minh in minhs:
		print ("ON: ", minh)
		confusionMatrices.append(buildConfidenceMatrixMultipleChromosomes(names, stem, bufferVal, spacing, minh, dist, tolerance, pasType))
	recalls = [c.recall() for c in confusionMatrices] #x axis on Precision Recall Curve
	precisions = [c.precision() for c in confusionMatrices] #y axis on Precision Recall Curve
	falsePositiveRates = [c.falsePositiveRate() for c in confusionMatrices] #x axis ROC curve
	truePositiveRates = [c.truePositiveRate() for c in confusionMatrices] #y axis ROC Curve
	fractionPeaksUnplaced = [c.fractionUnplaced() for c in confusionMatrices] 
	#graph ROC curve 
	#add (0,0) and (1,0) points if they are not present?
	print ("FPRS", falsePositiveRates)
	print ("TPRS", truePositiveRates)
	plt.plot(falsePositiveRates, truePositiveRates)
	plt.title("ROC Curve")
	plt.xlabel("FPR")
	plt.ylabel("TPR")
	plt.xlim(0,1.0)
	plt.ylim(0,1.0)
	plt.show()
	#plot PR Curve
	print ("recalls", recalls)
	print ("precisions: ", precisions)
	plt.plot(recalls, precisions)
	plt.title("PR Curve")
	plt.xlabel("Recall")
	plt.ylabel("Precision")
	plt.ylim(0,1.0)
	plt.xlim(0,1.0)
	plt.show()
	#plot ratio of unplaced to minH value
	plt.plot(minhs, fractionPeaksUnplaced)
	plt.title("Fraction of peaks unplaced vs Minimum Height of Peak")
	plt.xlabel("Peak Min Height")
	plt.ylabel("Fraction of total peaks unplaced")
	plt.ylim(0,1.0)
	plt.xlim(0,1.0)
	plt.show()

	
######################################################## Dataframe based methods #################################

#creates forward and reverse strand dictionaries for the balanced dataset given to it
def openBalancedValuesForTypeDataFrameVersion(toSeparate, pasType = ""):
	#add another column of all 0's
	toSeparate['peaksInRange'] = toSeparate.shape[0] * [0] #add row at the end of the dataframe with 0's to increment when a peak falls into it
	if pasType != "All":
		maskType = toSeparate["type"] == pasType
		maskedTrue = toSeparate[maskType]
		return maskedTrue
	else:
		return toSeparate

def placePeaksWithToleranceDataFrameVersion(peaks, balancedPos, balancedNegs, tolerance, sign, lenSeq):
	numberUnplacedPeaks = 0
	for peak in peaks:
		if sign == "-":
			peak = flipSequenceIndex(peak, lenSeq)
		balancedPos.loc[(balancedPos['strand'] == sign) & (((balancedPos['start'] <= peak) & (balancedPos['end'] >= peak)) 
			| ((balancedPos['start'] - tolerance <= peak) & (balancedPos['start'] >= peak)) 
			| ((balancedPos['end'] <= peak) & (balancedPos['end'] + tolerance >= peak))), "peaksInRange"] += 1
		balancedNegs.loc[(balancedNegs['strand'] == sign) & (((balancedNegs['start'] <= peak) & (balancedNegs['end'] >= peak)) 
			| ((balancedNegs['start'] - tolerance <= peak) & (balancedNegs['start'] >= peak)) 
			| ((balancedNegs['end'] <= peak) & (balancedNegs['end'] + tolerance >= peak))), "peaksInRange"] += 1
	return balancedPos, balancedNegs 



def fillConfMatrixDataFrameVersion( balancedPos, balancedNegs, totalPeaks = 0, countOncePerCluster = True):
	countTP = (balancedPos['peaksInRange'] > 0).sum() #peaks in true cluster
	countFP = (balancedNegs['peaksInRange'] > 0).sum() #peak in false cluster 
	countFN = (balancedPos['peaksInRange'] == 0).sum() #no peak in true cluster
	countTN = (balancedNegs['peaksInRange'] == 0).sum() #no peak in false cluster
	#print ("POS VALS: ")
	#print (balancedPos[(balancedPos['strand'] == '-') & (balancedPos['peaksInRange'] != 0)])
	#print (balancedPos[(balancedPos['strand'] == '+') & (balancedPos['peaksInRange'] != 0)])
	#print ("NEG VALS: ")
	#print (balancedNegs[(balancedNegs['strand'] == '-') & (balancedNegs['peaksInRange'] != 0)])
	#print (balancedNegs[(balancedNegs['strand'] == '+') & (balancedNegs['peaksInRange'] != 0)])
	if totalPeaks != 0:
		totalPeaksPlaced = balancedPos['peaksInRange'].sum() + balancedNegs['peaksInRange'].sum()
		totalPeaksUnplaced = totalPeaks - totalPeaksPlaced
		return countTP, countFP, countFN, countTN, totalPeaksUnplaced
	else: #no total peaks are given...piecewise computation
		return countTP, countFP, countFN, countTN, -1

def buildConfidenceMatrixOneChroDataFrameVersion(name, stem, b, s, minh, dist, tolerance, pasType = ""):
	fileName = stem + "chro" + name + "_NegSpaces" + str(b) + "_shifted" + str(s) + "Nts"
	print ("Opening: ", fileName)
	balancedPositives = openBalancedPositives(fileName)
	balancedNegatives = openBalancedNegatives(fileName)
	print ("Size balanced datasets: ", balancedPositives.shape[0])
	balancedPosWithCounts = openBalancedValuesForTypeDataFrameVersion(balancedPositives, pasType) 
	balancedNegsWithCounts = openBalancedValuesForTypeDataFrameVersion(balancedNegatives, pasType)
	print ("Size balanced w/ new col: ", balancedPosWithCounts.shape)
	#open peaks for the forward and reverse strand predictions
	predName = "chr" + name
	forward, reverse = openForwardReverse("../../aparentGenomeTesting/chromosomePredictions50/", predName)
	forwardPeaks = find_peaks_ChromosomeVersion(forward, minh, dist, (0.01, None)) 
	reversePeaks = find_peaks_ChromosomeVersion(reverse, minh, dist, (0.01, None)) 
	print ("Peaks found: ", len(forwardPeaks), " ", len(reversePeaks))
	#place peaks
	balancedPos, balancedNegs  = placePeaksWithToleranceDataFrameVersion(forwardPeaks, balancedPosWithCounts, balancedNegsWithCounts, tolerance, "+", forward.shape[0])
	print ("Placed + strand peaks")
	balancedPos, balancedNegs  = placePeaksWithToleranceDataFrameVersion(reversePeaks, balancedPos, balancedNegs, tolerance, "-", forward.shape[0])
	print ("Placed - strand peaks")
	totalPeaks = len(forwardPeaks) + len(reversePeaks)
	countTP, countFP, countFN, countTN, totalUnplaced = fillConfMatrixDataFrameVersion( balancedPos, balancedNegs, totalPeaks) 
	return countTP, countFP, countFN, countTN, totalUnplaced, totalPeaks


def buildConfidenceMatrixMultipleChromosomesDataFrameVersion(names, stem, bufferVal, spacing, minh, dist, tolerance, pasType = ""):
	countTPOverall = 0
	countFPOverall = 0
	countFNOverall = 0
	countTNOverall = 0
	countUnplaced = 0
	totalPeaks = 0
	for name in names:
		print ("On chromosome: ", name)
		countTP, countFP, countFN, countTN, unplaced, numberPeaks = buildConfidenceMatrixOneChroDataFrameVersion(name, stem, bufferVal, spacing, minh, dist, tolerance, pasType)
		countTPOverall += countTP
		countFPOverall += countFP
		countFNOverall += countFN
		countTNOverall += countTN
		countUnplaced += unplaced
		totalPeaks += numberPeaks
		print ("True Positives: ", countTPOverall, "False Positives: ", countFPOverall, "False Negatives: ", countFNOverall, "True Negatives: ", countTNOverall, "Unplaced: ",countUnplaced,  "Total Peaks: ", totalPeaks)
	return peakConfusionMatrix(countTPOverall, countFPOverall, countFNOverall, countTNOverall, countUnplaced, totalPeaks)	
	
def buildConfusionMatricesForGraphingDataFrameVersion(names, stem, bufferVal, spacing, minhs, dist, tolerance, rocTitle, prTitle, pasType = ""):
	confusionMatrices = []
	for minh in minhs:
		print ("ON: ", minh)
		confusionMatrices.append(buildConfidenceMatrixMultipleChromosomesDataFrameVersion(names, stem, bufferVal, spacing, minh, dist, tolerance, pasType))
	recalls = [c.recall() for c in confusionMatrices] #x axis on Precision Recall Curve
	precisions = [c.precision() for c in confusionMatrices] #y axis on Precision Recall Curve
	falsePositiveRates = [c.falsePositiveRate() for c in confusionMatrices] #x axis ROC curve
	truePositiveRates = [c.truePositiveRate() for c in confusionMatrices] #y axis ROC Curve
	fractionPeaksUnplaced = [c.fractionUnplaced() for c in confusionMatrices] 
	#graph ROC curve 
	#add (0,0) and (1,0) points if they are not present?
	print ("FPRS", falsePositiveRates)
	print ("TPRS", truePositiveRates)
	plt.plot(falsePositiveRates, truePositiveRates)
	plt.title("ROC Curve")
	plt.xlabel("FPR")
	plt.ylabel("TPR")
	plt.xlim(0,1.0)
	plt.ylim(0,1.0)
	plt.show()
	#plot PR Curve
	print ("recalls", recalls)
	print ("precisions: ", precisions)
	plt.plot(recalls, precisions)
	plt.title("PR Curve")
	plt.xlabel("Recall")
	plt.ylabel("Precision")
	plt.ylim(0,1.0)
	plt.xlim(0,1.0)
	plt.show()
	#plot ratio of unplaced to minH value
	plt.plot(minhs, fractionPeaksUnplaced)
	plt.title("Fraction of peaks unplaced vs Minimum Height of Peak")
	plt.xlabel("Peak Min Height")
	plt.ylabel("Fraction of total peaks unplaced")
	plt.ylim(0,1.0)
	plt.xlim(0,1.0)
	plt.show()
	
########################################################  Predicting peaks on a portion of the APARNET numpy arrays ####################################


def piecewisePeakPlacement(predictionsForward, predictionsRC, balancedPos, balancedNegs, tolerance, minh, dist):
	#using scipy find peaks on a slice around the balanced pos/neg regions to speed up process, since find_peaks only works on nearby values
	#newID = row['clusterID'] + "_LeftNegative"
	#newID = row['clusterID'] + "_RightNegative"
	for index, row in balancedPos.iterrows():
		#for each row, slice the prediction numpy around it (1000 nts)
		leftId = row['clusterID'] + "_LeftNegative"
		rightId = row['clusterID'] + '_RightNegative'
		if row['strand'] == "+":
			if row['start'] - 500 < 0:
				start = 0 
			else:
				start = row['start'] - 500
			if row['end'] + 500 > predictionsForward.size -1:
				end = predictionsForward.size-1
			else:
				end = row['end'] + 500
			sliceAround = predictionsForward[start:end]
			peaksInSlice = find_peaks_ChromosomeVersion(sliceAround, minh, dist)
			#print ("start: ", start,  " end: ", end)
			#print (peaksInSlice)
			#print ("number of peaks in slice: ", len(peaksInSlice), "minh: ", minh, "distance: ", dist)
			peaksInSliceCorrected = [x + start for x in peaksInSlice] #fix indexing issues that will have occured
			#print (peaksInSliceCorrected)
			x = [(peak >= row['start'] - tolerance) and (peak <= row['end'] + tolerance) for peak in peaksInSliceCorrected]
			#print (x)
			#input()
			if sum(x) > 0:
				#print ("INDEX: ", index)
				balancedPos.loc[index, "peaksInRange"] += 1
			#since negatives are either +/- 50 nts from this position, they are also contained in this slice
			#check if left in negatives
			searchLeftNegative = balancedNegs.index[balancedNegs['clusterID'] == leftId]
			if not searchLeftNegative.empty:
				#used left negative for balanced dataset
				#print ("left start: ", balancedNegs.start[searchLeftNegative].values[0])
				x = [(peak >= balancedNegs.start[searchLeftNegative].values[0] - tolerance) 
					and (peak <= balancedNegs.end[searchLeftNegative].values[0] + tolerance) for peak in peaksInSliceCorrected]
				if sum(x) > 0:
					balancedNegs.loc[searchLeftNegative, 'peaksInRange'] += 1
			#search if used right negative in balanced dataset and increment peaks which fall in that range
			searchRightNegative = balancedNegs.index[balancedNegs['clusterID'] == rightId]
			if not searchRightNegative.empty:
				x = [(peak >= balancedNegs.start[searchRightNegative].values[0] - tolerance) 
					and (peak <= balancedNegs.end[searchRightNegative].values[0] + tolerance) for peak in peaksInSliceCorrected]
				if sum(x) > 0:
					balancedNegs.loc[searchRightNegative, 'peaksInRange'] += 1
		else:
			#correct indexes
			rcIndexEnd = (predictionsForward.size -1) - row['start']
			#print ("total size: ", predictionsForward.size)
			#print ("row start: ", row['start'], " corr end: ", rcIndexEnd)
			if rcIndexEnd + 500 > predictionsForward.size - 1:
				end = predictionsForward.size -1
			else:
				end = rcIndexEnd + 500
			rcIndexStart = (predictionsForward.size -1) - row['end']
			#print ("row ens: ", row['end'], " corr start: ", rcIndexStart)
			if rcIndexStart - 500 < 0 :
				start = 0
			else:
				start = rcIndexStart - 500
			rcSlice = predictionsRC[start:end]
			peaksInSlice = find_peaks_ChromosomeVersion(rcSlice, minh, dist)
			peaksInSliceCorrected = [x + rcIndexStart for x in peaksInSlice] #correct index with RC indexing
			x = [(peak >= rcIndexStart - tolerance) and (peak <= rcIndexEnd + tolerance) for peak in peaksInSliceCorrected]
			#since negatives are either +/- 50 nts from this position, they are also contained in this slice
			if sum(x) > 0:
				balancedPos.loc[index, "peaksInRange"] += 1
			#check negatives in the same slice
			searchLeftNegative = balancedNegs.index[balancedNegs['clusterID'] == leftId]
			if not searchLeftNegative.empty:
				#used left negative for balanced dataset
				correctedStart = (predictionsForward.size - 1) - balancedNegs.end[searchLeftNegative].values[0]
				correctedEnd = (predictionsForward.size - 1) - balancedNegs.start[searchLeftNegative].values[0]
				x = [(peak >= correctedStart - tolerance) 
					and (peak <= correctedEnd + tolerance) for peak in peaksInSliceCorrected]
				if sum(x) > 0:
					balancedNegs.loc[searchLeftNegative, 'peaksInRange'] += 1
			#search if used right negative in balanced dataset and increment peaks which fall in that range
			searchRightNegative = balancedNegs.index[balancedNegs['clusterID'] == rightId]
			if not searchRightNegative.empty:
				correctedStart = (predictionsForward.size - 1) - balancedNegs.end[searchRightNegative].values[0]
				correctedEnd = (predictionsForward.size - 1) - balancedNegs.start[searchRightNegative].values[0]
				x = [(peak >= correctedStart - tolerance) 
					and (peak <= correctedEnd + tolerance) for peak in peaksInSliceCorrected]
				if sum(x) > 0:
					balancedNegs.loc[searchRightNegative, 'peaksInRange'] += 1
	return balancedPos, balancedNegs
			
def buildConfidenceMatrixOneChroPiecewise(predictionLocation, name, stem, b, s, minh, dist, tolerance, pasType = ""):
	fileName = stem + "chro" + name + "_NegSpaces" + str(b) + "_shifted" + str(s) + "Nts"
	print ("Opening: ", fileName)
	balancedPositives = openBalancedPositives(fileName)
	balancedNegatives = openBalancedNegatives(fileName)
	print ("BP Size: ", balancedPositives.shape)
	print ("BN Size: ", balancedNegatives.shape)
	#add peaksInRange column to count peaks in the cluster/empty clusters
	balancedPosWithCounts = openBalancedValuesForTypeDataFrameVersion(balancedPositives, pasType) 
	balancedNegsWithCounts = openBalancedValuesForTypeDataFrameVersion(balancedNegatives, pasType)
	print ("BP Size: ", balancedPosWithCounts.shape)
	print ("BN Size: ", balancedNegsWithCounts.shape)
	#print (balancedPosWithCounts)
	#open peaks for the forward and reverse strand predictions
	predName = "chr" + name
	forward, reverse = openForwardReverse(predictionLocation, predName)
	print ("Total predictions length: ", forward.size)
	balancedPos, balancedNegs = piecewisePeakPlacement(forward, reverse, balancedPosWithCounts, balancedNegsWithCounts, tolerance, minh, dist)
	countTP, countFP, countFN, countTN, dummy = fillConfMatrixDataFrameVersion( balancedPos, balancedNegs) #should still be able to score these the same way, adding default of totalPeaks = [] if not provided, will return -1 for the dummy var
	return countTP, countFP, countFN, countTN
	


def buildConfidenceMatrixMultipleChromosomesPiecewise(predictionLocation, names, stem, bufferVal, spacing, minh, dist, tolerance, pasType = ""):
	countTPOverall = 0
	countFPOverall = 0
	countFNOverall = 0
	countTNOverall = 0
	for name in names:
		print ("--------------------------------------------------")
		print ("On chromosome: ", name)
		countTP, countFP, countFN, countTN = buildConfidenceMatrixOneChroPiecewise(predictionLocation, name, stem, bufferVal, spacing, minh, dist, tolerance, pasType)
		countTPOverall += countTP
		countFPOverall += countFP
		countFNOverall += countFN
		countTNOverall += countTN
		print ("True Positives: ", countTPOverall, "False Positives: ", countFPOverall, "False Negatives: ", countFNOverall, "True Negatives: ", countTNOverall)
		print ("--------------------------------------------------")
	return peakConfusionMatrix(countTPOverall, countFPOverall, countFNOverall, countTNOverall)					


	

def buildConfusionMatricesForGraphingPiecewise(predictionLocation, names, stem, bufferVal, spacing, minhs, dist, tolerance, rocTitle, prTitle, pasType = ""):
	confusionMatrices = []
	for minh in minhs:
		print ("ON: ", minh)
		confusionMatrices.append(buildConfidenceMatrixMultipleChromosomesPiecewise(predictionLocation, names, stem, bufferVal, spacing, minh, dist, tolerance, pasType))
	recalls = [c.recall() for c in confusionMatrices] #x axis on Precision Recall Curve
	precisions = [c.precision() for c in confusionMatrices] #y axis on Precision Recall Curve
	falsePositiveRates = [c.falsePositiveRate() for c in confusionMatrices] #x axis ROC curve
	truePositiveRates = [c.truePositiveRate() for c in confusionMatrices] #y axis ROC Curve 
	#graph ROC curve 
	print ("FPRS", falsePositiveRates)
	print ("TPRS", truePositiveRates)
	plt.plot(falsePositiveRates, truePositiveRates)
	plt.title("ROC Curve")
	plt.xlabel("FPR")
	plt.ylabel("TPR")
	plt.xlim(0,1.0)
	plt.ylim(0,1.0)
	plt.show()
	#plot PR Curve
	print ("recalls", recalls)
	print ("precisions: ", precisions)
	plt.plot(recalls, precisions)
	plt.title("PR Curve")
	plt.xlabel("Recall")
	plt.ylabel("Precision")
	plt.ylim(0,1.0)
	plt.xlim(0,1.0)
	plt.show()
	

def extractMaximumValues(names, predictionLocation):
	maxVal = float('-inf')
	for name in names:
		predName = "chr" + name
		forward, reverse = openForwardReverse(predictionLocation, predName)
		maxVal = max(maxVal, max(forward), max(reverse))
	return maxVal


def openBestThresholds():
	f_open = open("bestThresholds.txt", "r")
	ts = [float(x.strip()) for x in f_open.readlines()]
	return ts	


def estimateThresholds(name, predictionLocation):
	predName = "chr" + name
	forward, reverse = openForwardReverse(predictionLocation, predName)
	both = np.concatenate((forward, reverse))
	print ("forward: ", forward.shape, " reverse complement: ", reverse.shape)
	print ("combined: ", both.shape)
	thresholds = []
	for i in range(1,101):
		#figuring out best thresholds to use overall
		print ("ON %: ", i)
		pfound = np.percentile(both, i)
		thresholds.append(pfound)
		print ("Threshold found: ", pfound)
	print (thresholds)
	f = open("bestThresholds.txt", "w")
	for t in thresholds:
		f.write(str(t) + '\n')
	f.close()
	print ("saved thresholds")
	

#making version to create confusion matrices for a single chromosome across mult. minhs
#keeping bufferValue, spacing, tolerance fixed at 1, 50, 0 for now
#varying spacings used in find_peaks, thresholds used in find_peaks, pasType being evalauted
#saved output as csv files (becuase pandas dataframe to csv is easy)
def singleChromosomeCMSpanning(cmLocation, predictionLocation, names, datasetsLocation, bufferValue, spacing, minhs, distances, tolerance, pasTypes):
	for name in names:
		columns_dict = {'name': [], 'pasType': [], 'bufferVal':[], 'spacing':[], 'minh':[], 'distance':[], 'tolerance':[], 'TruePositives':[], 'FalsePositives':[], 'FalseNegatives':[], 'TrueNegatives':[]}
		for dist in distances: 
			for minh in minhs:
				for ptype in pasTypes:	
					print ("h: ", minh, " distance: ", dist, "pasType: ", ptype)
					countTP, countFP, countFN, countTN = buildConfidenceMatrixOneChroPiecewise(predictionLocation,
						name,
						datasetsLocation,
						bufferValue,
						spacing,
						minh,
						dist,
						tolerance,
						pasType = ptype)
					#nameStorage = name + "_bufferVal_"+ str(bufferValue) + "_spacing_" + str(spacing) #used to specify which dataset is being used
					columns_dict['name'].append(name)
					columns_dict['pasType'].append(ptype)
					columns_dict['bufferVal'].append(bufferValue)
					columns_dict['spacing'].append(spacing)
					columns_dict['minh'].append(minh)
					columns_dict['distance'].append(dist)
					columns_dict['tolerance'].append(tolerance)
					columns_dict['TruePositives'].append(countTP)
					columns_dict['FalsePositives'].append(countFP)
					columns_dict['FalseNegatives'].append(countFN)
					columns_dict['TrueNegatives'].append(countTN)
					print ("True Positives: ", countTP, "False Positives: ", countFP, "False Negatives: ", countFN, "True Negatives: ", countTN)
					print (" ")
		toDF = pd.DataFrame(columns_dict)
		#save csv for chromosome
		saveName = cmLocation + name + "_ConfusionMatrices.csv"
		toDF.to_csv(saveName)
						


names = ["22"]
bufferVals = 1
spacing = 50
distances = [1, 25, 50]
tolerances = 0
pasTypes = ['All', 'IN', 'TE', 'IG', 'AI', 'EX', 'DS', 'AE', 'AU']
cutoffs = openBestThresholds() #made using percentiles 1-100 of chr1 prediction values 
cutoffs = [0.0,0.05,0.1,0.2]
# local laptop locations
chromsomeLocationsLocal ="../../aparentGenomeTesting/chromosomePredictions50/"
fastaLocationLocal = "../../aparentGenomeTesting/fastas/"
cMLocation = "./ConfusionMatrices/"
#for use on research servers
fastaLocationBicycle =  "../fastas/" 
predictionLocationBicycle = "../chromosomePredictions50/"
datasetsLocation = "./datasets/"
#buildConfusionMatricesForGraphingPiecewise(chromsomeLocationsLocal, names, "./datasets/", 1, 50, cutoffs, 1, 20, "testingChrYROC", "restingChrYPrecision", "TE")
singleChromosomeCMSpanning(cMLocation, predictionLocationBicycle, names, datasetsLocation, bufferVals, spacing, cutoffs, distances, tolerances, pasTypes)
		
		
	
	


