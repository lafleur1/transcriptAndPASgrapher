import re, sys, math, random, os, time
import numpy as np
import signal, subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Alphabet import IUPAC
import sklearn.metrics as metrics
import pandas as pd

from keras.preprocessing import sequence
from keras.models import Model
from keras.layers import Bidirectional, Input, concatenate, add
from keras.layers.core import Dense, Dropout, Flatten
from keras.layers.convolutional import Conv1D
from keras.layers.pooling import MaxPooling1D, AveragePooling1D
from keras.layers.recurrent import LSTM
from keras.callbacks import ModelCheckpoint, EarlyStopping
from keras.layers.advanced_activations import PReLU
from sklearn.metrics import average_precision_score, roc_auc_score
import tensorflow as tf
tf.get_logger().setLevel('INFO')
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2' 
os.environ["KMP_WARNINGS"] = "FALSE"

import tensorflow as tf
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)
#functions from DeepPASTA shape_asign_per_nucleotide.py 

def replaceSupportingCharacter(finalShape):
	for i in range(len(finalShape)):
		if finalShape[i] == '<':
			finalShape[i] = 'L'
		elif finalShape[i] == '>':
			finalShape[i] = 'R'
	return finalShape

def updateForIAndM(finalShape, pattern):
	firstIndex = -1
	for i in range(len(finalShape)):
		if finalShape[i] == '(' and finalShape[i+1] != '(':
			firstIndex = (i+1)
		elif (finalShape[i] == '.' or finalShape[i] == '>') and finalShape[i+1] == ')' and firstIndex != -1:
			lastIndex = i
			modified = False
			if len(re.findall(pattern, ''.join(finalShape)[firstIndex:lastIndex])) == 1:
				l = firstIndex
				while (l <= lastIndex):
					if finalShape[l] == '<':
						finalShape[l] = 'L'
					elif finalShape[l] == '>':
						finalShape[l] = 'R'
					elif finalShape[l] == '.':
						finalShape[l] = 'I'
					modified = True
					l = l + 1
			elif len(re.findall(pattern, ''.join(finalShape)[firstIndex:lastIndex])) > 1:
				l = firstIndex
				while (l <= lastIndex):
					if finalShape[l] == '<':
						finalShape[l] = 'L'
					elif finalShape[l] == '>':
						finalShape[l] = 'R'
					elif finalShape[l] == '.':
						finalShape[l] = 'M'
					modified = True
					l = l + 1
			if modified == True:
				l = firstIndex - 1
				m = lastIndex + 1
				while (finalShape[l] == '(' and finalShape[m] == ')'):
					finalShape[l] = '<'
					finalShape[m] = '>'
					l = l - 1
					m = m + 1
			firstIndex = -1
		elif finalShape[i] == '(' and finalShape[i+1] == '(':
			firstIndex = -1
	return finalShape

def assignUnpairedPosition(shape):
	for i in range(len(shape)):
		if shape[i] == '.':
			shape[i] = 'U'
	return shape

def shapeAtEachPosition(abstractShape):
	beforeOrAfterAnyParanthesis = False
	finalShape = ['*'] * len(abstractShape)
	lastParanthesisIndex = 0
	for i in range(len(abstractShape)):
		if abstractShape[i] == ')':
			lastParanthesisIndex = i
		finalShape[i] = abstractShape[i]
	lastbracket = ''
	firstIndex = 0
	lastIndex = 0
	for i in range(len(abstractShape)):
		if beforeOrAfterAnyParanthesis == False and abstractShape[i] == '(':
			beforeOrAfterAnyParanthesis = True
		if beforeOrAfterAnyParanthesis == True and i > lastParanthesisIndex:
			beforeOrAfterAnyParanthesis = False
			
		if abstractShape[i] == "." and beforeOrAfterAnyParanthesis == False:
			finalShape[i] = 'E'
		elif (abstractShape[i] == '(' or abstractShape[i] == ')') and abstractShape[i+1] == '.':
			lastbracket = abstractShape[i]
			firstIndex = (i+1)
		elif abstractShape[i] == '.' and (abstractShape[i+1] == ')' or abstractShape[i+1] == '('):
			lastIndex = i
			if lastbracket == '(' and abstractShape[i+1] == ')':
				l = firstIndex
				while (l <= lastIndex):
					finalShape[l] = 'H'
					l = l+1
				l = firstIndex - 1
				m = lastIndex + 1
				while (abstractShape[l] == '(' and abstractShape[m] == ')'):
					finalShape[l] = '<'
					finalShape[m] = '>'
					l = l - 1
					m = m + 1
	
			lastbracket = ''
	
	finalShape = updateForIAndM(finalShape, '<+\w+>+')
	count = 0
	while ('(' in ''.join(finalShape)) or (')' in ''.join(finalShape)):
		newFinalShape = updateForIAndM(finalShape, '<+\w*>+')
		if newFinalShape == finalShape:
			if count < 30:
				count = count + 1
			else:
				break
		finalShape = newFinalShape
	finalShape = replaceSupportingCharacter(finalShape)	
	if '.' in ''.join(finalShape):
		finalShape = assignUnpairedPosition(finalShape)

	return finalShape


#From DeepPASTA_polyA_site_prediction_testing.py 
#
def oneHotEncodingForSeq(rawSeqList):
	if len(rawSeqList) != 0:
		encodedSeq = np.zeros((len(rawSeqList), len(rawSeqList[0]), 5))
		for i in range(len(rawSeqList)):
			sequence = rawSeqList[i]
			j = 0
			for s in sequence:
				if s == 'A' or s == 'a':
					encodedSeq[i][j] = [1,0,0,0,0]
				elif s == 'T' or s == 't':
					encodedSeq[i][j] = [0,1,0,0,0]
				elif s == 'C' or s == 'c':
					encodedSeq[i][j] = [0,0,1,0,0]
				elif s == 'G' or s == 'g':
					encodedSeq[i][j] = [0,0,0,1,0]
				elif s == 'N' or s == 'n':
					encodedSeq[i][j] = [0,0,0,0,1]
				else:
					#print>>sys.stderr, 'ERROR: Unwanted nucleotide: ' + s
					print ('ERROR: Unwanted nucleotide: ' + s)
				j = j + 1
		return encodedSeq
	else:
		return 0

def oneHotEncodingForSS(rawStructureList):
	if len(rawStructureList) != 0:
		encodedStructure = np.zeros((len(rawStructureList), len(rawStructureList[0]), 7))
		for i in range(len(rawStructureList)):
			structure = rawStructureList[i]
			j = 0
			for s in structure:
				if s == 'U':
					encodedStructure[i][j] = [1,0,0,0,0,0,0]
				elif s == 'E':
					encodedStructure[i][j] = [0,1,0,0,0,0,0]
				elif s == 'L':
					encodedStructure[i][j] = [0,0,1,0,0,0,0]
				elif s == 'R':
					encodedStructure[i][j] = [0,0,0,1,0,0,0]
				elif s == 'H':
					encodedStructure[i][j] = [0,0,0,0,1,0,0]
				elif s == 'M':
					encodedStructure[i][j] = [0,0,0,0,0,1,0]
				elif s == 'I':
					encodedStructure[i][j] = [0,0,0,0,0,0,1]
				else:
					#print>>sys.stderr, 'Warning: Unwanted character ' + s
					print ('Warning: Unwanted character ' + s)
				j = j + 1
		return encodedStructure
	else:
		return 0


def sequenceModel(seqInput):
	seqCov = Conv1D(filters=512,
	kernel_size=8,
	padding = "valid",
	input_shape =(200, 5),
	activation="relu",
	strides=1)(seqInput) 

	seqPool = MaxPooling1D(pool_size = 3, strides = 3)(seqCov)
	seqDout1 = Dropout(rate = 0.7)(seqPool)
	seqBiLstm = Bidirectional(LSTM(units = 128, return_sequences = True))(seqDout1)
	seqDout2 = Dropout(rate = 0.7)(seqBiLstm)
	seqFlat = Flatten()(seqDout2)
	seqDen2 = Dense(256, kernel_initializer='glorot_uniform', activation = 'relu')(seqFlat)
	seqDout4 = Dropout(rate = 0.7)(seqDen2)

	return seqDout4

def structureSubModel(ssInput):
	ssConv = Conv1D(filters=256,
	kernel_size=12,
	padding = "valid",
	activation="relu",
	strides=1)(ssInput)
	ssPool = AveragePooling1D(pool_size = 5, strides = 5)(ssConv)
	ssDout1 = Dropout(rate=0.7)(ssPool)
	seqBiLstm = Bidirectional(LSTM(units = 128, return_sequences = True))(ssDout1)
	seqDout2 = Dropout(rate = 0.7)(seqBiLstm)
	ssFlat = Flatten()(seqDout2)
	ssDen1 = Dense(256, kernel_initializer='glorot_uniform', activation = 'relu')(ssFlat)
	ssDout2 = Dropout(rate=0.7)(ssDen1)

	return ssDout2
	
#RNA shapes: https://academic.oup.com/bioinformatics/article/22/4/500/184565 

#FIRST PY: COMBINING_SUBSTRUCTURE ######################################################################

#combines first portion of 1-100 

def processRNAShapesOutput(inputFile):	
	firstSubStructure = []
	secondSubStructure = []
	outputList = []
	current = 'none'
	if inputFile != "":
		for line in open(inputFile):
			if ">" in line: #first line in file (fasta name)
				###
				if current == 'second': #enters if there is another fasta, not important for me
					for first in firstSubStructure:
						for second in secondSubStructure:
							outputList.append(first+second)
					current = 'none'
				###
				outputList.append(line[:-1])
				firstSubStructure = []
				secondSubStructure = []
			elif ('101' in line) and ('200' in line): 
				current = 'second'
			elif ('1' in line) and ('100' in line): #second line in file
				current = 'first' #set that it is in 1-100 area
			elif ("(" in line or ")" in line or "." in line):
				if current == 'first':
					firstSubStructure.append(line[:-1]) #removing the end line character between lines
				else:
					secondSubStructure.append(line)
	if current == 'second':
		for first in firstSubStructure:
			for second in secondSubStructure:
				outputList.append(first+second)
		current = 'none'
	passedOn = outputList[0:4]
	if len(passedOn) != 4:
		print (passedOn)
		print (outputList)
	count = 0
	outputList3 = []
	lastShape = ''
	for line in passedOn:
		if ">" in line:
			outputList3.append(line)
			count = 0
		elif "(" in line or ")" in line or "." in line:
			lastShape = ''.join(shapeAtEachPosition(line))
			outputList3.append(lastShape[:-1])
			count = count + 1
	return outputList3

################################################
#condensed DeepPASTA scripts 
################################################

def compileModel():
	# Building deep learning model
	training_net = []
	structureSeqLength = 200
	# deep learning sub-model for sequence
	seqInput = Input(shape = (200, 5))
	seqModel = sequenceModel(seqInput)
	training_net.append(seqModel)

	# deep learning sub-model for structure
	ss_training_net = []
	ssInput1 = Input(shape = (structureSeqLength, 7))
	ssInput2 = Input(shape = (structureSeqLength, 7))
	ssInput3 = Input(shape = (structureSeqLength, 7))

	ss_training_net.append(structureSubModel(ssInput1))
	ss_training_net.append(structureSubModel(ssInput2))
	ss_training_net.append(structureSubModel(ssInput3))

	ss_merged_model = add(ss_training_net)
	ss_den1 = Dense(256, kernel_initializer = 'glorot_uniform', activation = 'relu')(ss_merged_model)
	ss_dout1 = Dropout(rate = 0.7)(ss_den1)
	training_net.append(ss_dout1)
	merged_model = concatenate(training_net)

	den1 = Dense(256, kernel_initializer = 'glorot_uniform', activation = 'relu')(merged_model)
	dout1 = Dropout(rate = 0.7)(den1)

	den2 = Dense(128, kernel_initializer = 'glorot_uniform', activation = 'relu')(dout1)
	dout2 = Dropout(rate = 0.7)(den2)
	den3 = Dense(64, activation = 'relu')(dout2)
	den4 = Dense(1, activation = 'sigmoid')(den3)
	model = Model(inputs = [seqInput, ssInput1, ssInput2, ssInput3], outputs = den4)

	model.compile(loss='binary_crossentropy', optimizer='nadam', metrics=['accuracy'])
	model.load_weights('DeepPASTA_polyA_site_learned.hdf5')
	return model

#input of 200 nt sequence 
#input file name from RNAShapes
#input of loaded model 
#output of model score 
def runModel(sequence,inputFile, model ):
	outSS = processRNAShapesOutput(inputFile) #list of three 
	if len(outSS) != 4:
		print ("ERROR!")
		print (outSS)
	encodedTestingSeq = oneHotEncodingForSeq([sequence])
	encodedTestingStructure1 = oneHotEncodingForSS([outSS[1]])
	encodedTestingStructure2 = oneHotEncodingForSS([outSS[2]])
	encodedTestingStructure3 = oneHotEncodingForSS([outSS[3]])
	testingData = []
	testingData.append(encodedTestingSeq)
	testingData.append(encodedTestingStructure1)
	testingData.append(encodedTestingStructure2)
	testingData.append(encodedTestingStructure3)
	testresult1 = model.predict(testingData, batch_size = 2042, verbose = 0)
	return testresult1




'''
#predicting and storing DeepPASTA for each Postiive/False Cluster 
def piecewisePeakPlacement(forwardStrand, reverseCompStrand, balancedDataset):
	#using scipy find peaks on a slice around the balanced pos/neg regions to speed up process, since find_peaks only works on nearby values
	#newID = row['clusterID'] + "_LeftNegative"
	#newID = row['clusterID'] + "_RightNegative"
	lenStrand = len(forwardStrand)
	for index, row in balancedDataset.iterrows():
		#for each row, slice the prediction numpy around it (1000 nts)
		if row['strand'] == "+":
			if row['start'] - 500 < 0:
				start = 0 
			else:
				start = row['start'] - 500
			if row['end'] + 500 > lenStrand -1:
				end = lenStrand-1
			else:
				end = row['end'] + 500
			predictionString = prepTargetAreaForPrediction(start, end, forwardStrand)
			
			
			
		else:
			#correct indexes
			rcIndexEnd = (lenStrand -1) - row['start']
			if rcIndexEnd + 500 > lenStrand - 1:
				end = lenStrand -1
			else:
				end = rcIndexEnd + 500
			rcIndexStart = (lenStrand -1) - row['end']
			#print ("row ens: ", row['end'], " corr start: ", rcIndexStart)
			if rcIndexStart - 500 < 0 :
				start = 0
			else:
				start = rcIndexStart - 500
			rcSlice = reverseCompStrand[start:end]
'''	

def runPASTA(sequence, loadedModel):
	#make fasta from sequence slice
	outputs = []
	current_output = 0
	for i in range(0,len(sequence) - (200-1)):
		strSeq = str(sequence[i:i+200])
		dname = "1"
		f = open("dummy.fa", "w")
		f.write(">" + dname + "\n" + strSeq + "\n")
		f.close()
		# running RNAshapes
		comm1 = "./RNAshapes -f dummy.fa -s -c 5 -t 1 -w 100 -W 100 -O 'D{%s\n}' > rnaShapesOutput.txt"
		p = subprocess.Popen(comm1, shell = True, stdout = subprocess.PIPE)
		p.wait()
		out = runModel(strSeq,"rnaShapesOutput.txt", loadedModel)
		outputs.append(out[0][0])
	return outputs
	
def prepTargetAreaForPrediction(targetStart, targetEnd, chrS):
	return str(chrS[targetStart-99:targetEnd + 101])

def getTargetNucleotides(strPredictedOn):
	return strPredictedOn[99:len(strPredictedOn)-99-1]

	
def openAllNegatives(name):
	negName = name + "AllNegatives.csv"
	return pd.read_csv( negName,  dtype = {"seqName": str}) 

def openBalancedNegatives(name):
	negativesName = name + "BalancedNegatives.csv"
	return pd.read_csv( negativesName, dtype = {"seqName": str}) 

def openBalancedPositives(name):
	positivesName = name + "BalancedPositives.csv"
	return pd.read_csv( positivesName, dtype = {"seqName": str}) 
	
	

#predicting and storing DeepPASTA for each Postiive/False Cluster 
def predictInClusters(forwardStrand, reverseCompStrand, balancedDataset, loadedModel):
	lenStrand = len(forwardStrand)
	for index, row in balancedDataset.iterrows():
		#print (row)
		#for each row, slice the prediction numpy around it (1000 nts)
		if row['strand'] == "+":
			predictionString = prepTargetAreaForPrediction(row['start'], row['end'], forwardStrand)
			try:
				deepPASTAPreds = runPASTA(predictionString,loadedModel)
		#		print ("Index: ", index)
		#		print ("Prediction for: ", getTargetNucleotides(predictionString))
		#		print ("Predicted values: ", deepPASTAPreds)
				balancedDataset.at[index, "DeepPASTAPredictions"] = deepPASTAPreds
			except:
				print ("FAILURE ON: ", row['clusterID'])
						
		else:
			#print (row['start'], row['end'])
			rcIndexEnd = (lenStrand -1) - row['start']
			rcIndexStart = (lenStrand -1) - row['end']
		#	print (rcIndexStart, rcIndexEnd)
			rcPredictionString = prepTargetAreaForPrediction(rcIndexStart, rcIndexEnd, reverseCompStrand)
			#print (rcPredictionString)
			try:
				deepPASTAPreds = runPASTA(rcPredictionString,loadedModel)
		#		print ("Index: ", index)
		#		print ("Prediction for: ", getTargetNucleotides(rcPredictionString))
		#		print ("Predicted values: ", deepPASTAPreds)
				balancedDataset.at[index, "DeepPASTAPredictions"] = deepPASTAPreds
			except:
				print ("FAILURE ON: ", row['clusterID'])
	return balancedDataset
	



def runDeepPASTAClusterPredictions(loadedModel, fastaLocation, name, stem, b, s):
	fileName = stem + "chro" + name + "_NegSpaces" + str(b) + "_shifted" + str(s) + "Nts"
	print ("Opening: ", fileName)
	balancedPositives = openBalancedPositives(fileName)
	balancedNegatives = openBalancedNegatives(fileName)
	#adding column to both dataframes
	#toSeparate['DeepPASTAPredictions'] = toSeparate.shape[0] * [0]
	balancedPositives['DeepPASTAPredictions'] = balancedPositives.shape[0] * [[-1]]
	balancedNegatives['DeepPASTAPredictions'] = balancedNegatives.shape[0] * [[-1]]
	#print (balancedPositives)
	print ("BP Size: ", balancedPositives.shape)
	print ("BN Size: ", balancedNegatives.shape)
	#open fastas for predictions
	totalName = fastaLocation + "chr" + name + ".fasta"
	seqObj = SeqIO.read(totalName, "fasta")
	forward_seq = str(seqObj.seq)
	reverseComp_seq = str(seqObj.reverse_complement().seq)
	print ("Forward len: ", len(forward_seq))
	print ("Reverse len: ", len(reverseComp_seq))
	#predicting positive values
	balPos = predictInClusters(forward_seq, reverseComp_seq, balancedPositives, loadedModel)
	print ("saving positives")
	balPos.to_csv(fileName + "BalancedPositives.csv")
	#predicting negatives values
	balNeg = predictInClusters(forward_seq, reverseComp_seq, balancedNegatives, loadedModel)
	print ("saving negatives")
	balNeg.to_csv(fileName + "BalancedNegatives.csv")
	

#load DeepPASTA model
loadedModel = compileModel()

#names for chromosomes
chromosomes = ["15", "16", "17", "18", "19", "20", "21", "22", "Y"]
fastaLocationLocal = "../../../aparentGenomeTesting/fastas/"

names = []
if len(sys.argv) != 1:
	for arg in sys.argv[1:]:
		names.append(arg)
else:
	names = ["22"]

runDeepPASTAClusterPredictions(loadedModel, fastaLocationLocal, "Y", "./predictions/", 1, 50)


	
'''
chroName = "Y"
fastaLocs = ""
totalName = fastaLocs + "chr" + chroName + ".fasta"
chroSeq = SeqIO.read(totalName, "fasta")
sliceY = chroSeq.seq[10000:10000 + 202]
start_time = time.time()
print ("START TIME: ")
outVals = runPASTA(sliceY, loadedModel)
print ("ELAPSED TIME: ", time.time() - start_time)
print ("number outputs: ", len(outVals))
print ("sequence length: ", len(sliceY))
print ("ALL OUTPUTS: ")
print (outVals)
'''


	

'''
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
			

'''
