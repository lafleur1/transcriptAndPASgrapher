import re, sys, math, random, os, time
import numpy as np
import signal, subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Alphabet import IUPAC

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
stderr = sys.stderr
sys.stderr = open(os.devnull, 'w')
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




def runPASTA(sequence):
	#make fasta from sequence slice
	outputs = []
	current_output = 0
	loadedModel = compileModel()
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

chroName = "Y"
fastaLocs = ""
totalName = fastaLocs + "chr" + chroName + ".fasta"
chroSeq = SeqIO.read(totalName, "fasta")
sliceY = chroSeq.seq[500000:500000 + 1199]
print ("TOTAL SEQ: ", len(chroSeq.seq))
print (len(sliceY))
start_time = time.time()
print ("START TIME: ")
outVals = runPASTA(sliceY)
print ("ELAPSED TIME: ", time.time() - start_time)
print ("number outputs: ", len(outVals))
print ("sequence length: ", len(sliceY))
print ("ALL OUTPUTS: ")
print (outVals)





