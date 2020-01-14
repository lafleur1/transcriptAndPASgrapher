#extracting sequences near the true values in same regions for true negative datasets for APA model evaluation
#APARENT: 206 nt, overlapping
#DeepPASTA: 200 nts, polyA site in middle
#
#using same method as Leung et al (see  Inference of the human polyadenylation code supplementary material):
 #same padding requirement (i.e. that spacing between two PAS signals must allow four negative regions 
 #shift left and right by 50 nts but still within the region boundaries 
 
import pyensembl
import pandas as pd 
from Bio import SeqIO

#Note: Ensembl is 1-based, polyAsite is 0-based

#rules used for database: 
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


#open all PAS dataframes

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

#create new dataframe with negative examples


#save for later when verifying that new false negative falls in the region range

'''

def getLastTranscriptExon(transcript_strand, transcript_exons):
	#extracts last exon coordinates from PyEnsembl transcript object
	transcript_exons = sorted(transcript_exons)
	if transcript_strand == "+":
		return transcript_exons[-1]
	else: 
		return transcript_exons[0]

def checkOtherExons(transcript_strand, transcript_exons, pas_start, pas_end):
	winningRange = 0
	case = 0
	for trange in transcript_exons:
		if pas_start >= trange[0] and pas_end <= trange[1]: #exact match
			winningRange = trange
			case = "exact match"
			break
		elif pas_start >= trange[0] and pas_start <= trange[1]: #start in range and end is not
			winningRange = trange
			case = "start match"
		elif pas_end >= trange[0] and pas_end <= trange[1]: #end in range start is not
			winningRange = trange
			case = "end match"
		else:
			case = "no match"
	return winningRange, case

def doubleCheckIdentities(location, pas_type, pas_start, pas_end, pas_strand, pas_chromosome):
	matches = False
	newType = ""
	genes_in_range = ensembl.genes_at_locus(contig = pas_chromosome, position = pas_start, end = pas_end)
	for gene in genes_in_range:
		for transcript in gene.transcripts:
			exons = [[exon.start, exon.end] for exon in transcript.exons]
			last_exon = getLastTranscriptExon(transcript.strand, exons)
			if pas_start >= last_exon[0] and pas_end <= pas_end:
				newType = "TE"
			else:
				#check other exon
				if 
	if pas_type == "TE":
		#find the exon that it is in

		
	elif pas_type == "EX":
	elif pas_type == "IN":
	elif pas_type == "DS":
	elif pas_type == "AE":
	elif pas_type == "AI":
	elif pas_type == "AU":
	elif pas_type == "IG":
	else:
		print ("error!!!! Invalid PAS type")
		
'''
'''
print (auPAS.iloc[0])
row1 = auPAS.iloc[0]
print (row1)
print (row1['type'])

'''




#fails the pos signal if there are "N" in the true region, or in the region used to predict the true pas region [start - 205, end + 205]
def shiftLeftAndRightNegativesFilterPositives(pasRow, allPAS, spacingValue, shiftValue, fastaSeq, strictNFiltering = True, strictFilteringValue = 205):
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
		leftPassed = False
		
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
		rightPassed = False
		### LEFT
		if leftNCount == 0:
			leftCheck = pasRow['start'] - checkRange
			filterVals1 = (allPAS['strand'] == pasRow['strand']) & (allPAS['seqName'] == pasRow['seqName']) & (allPAS['end'] >= leftCheck) & (allPAS['start'] >= leftCheck) & (allPAS['start'] < pasRow['start'])
			filterVals2 = (allPAS['strand'] == pasRow['strand']) & (allPAS['seqName'] == pasRow['seqName']) & (allPAS['end'] >= leftCheck) & (allPAS['start'] <= leftCheck)
			otherPAS1Left = allPAS[filterVals1]
			otherPAS2Left = allPAS[filterVals2]
			sumInForbiddenRange = otherPAS1Left.shape[0] + otherPAS2Left.shape[0]
			if sumInForbiddenRange == 0: 
				leftPassed = True
		else:
			print ("Sequence: ", leftSlice)
			leftStart = -1
			leftEnd = -1
		#####
		##### RIGHT 
		if rightNCount ==0:
			rightCheck = pasRow['end'] + checkRange
			filterVals1 = (allPAS['strand'] == pasRow['strand']) & (allPAS['seqName'] == pasRow['seqName']) & (allPAS['start'] <= rightCheck) & (allPAS['end'] <= rightCheck) & (allPAS['start'] >= pasRow['end']) 
			filterVals2 = (allPAS['strand'] == pasRow['strand']) & (allPAS['seqName'] == pasRow['seqName']) & (allPAS['start'] <= rightCheck) & (allPAS['end'] >= rightCheck)
			otherPAS1Right = allPAS[filterVals1]
			otherPAS2Right = allPAS[filterVals2]
			sumInForbiddenRangeRight = otherPAS1Right.shape[0] + otherPAS2Right.shape[0]
			if sumInForbiddenRangeRight == 0:
				rightPassed = True
		else:
			print ("Sequence: ", rightSlice)
			rightStart = -1
			rightEnd = -1
		#####
		'''
		if not rightPassed and not leftPassed:
			print ("---------------------------------")
			print (pasRow)
			print ("RANGE :", leftCheck, " ", rightCheck)
			print ("LEFT: ")
			print (otherPAS1Left)
			print (otherPAS2Left)
			print ("RIGHT: ")
			print (otherPAS1Right)
			print (otherPAS2Right)
			print ("----------------------------------")
		'''
		return leftPassed, leftStart, leftEnd, rightPassed, rightStart, rightEnd
	else:
		print ("Sequence: ", filteringTrueSequence)
		return False, -2, -2, False, -2, -2
	

def createNegativeDataSet(trueVals, spacingValue, shiftValue, fileName, fastaSeq):
	#create dictionary for negative valeus left and right 
	#["seqName",  "start" , "end",  "clusterID",  "avgTPM",  "strand",   "percentSupporting",   "protocolsSupporting",  "avgTPM2",   "type",   "upstreamClusters"]
	blankDict = {"seqName":[], "start":[], "end":[], "clusterID":[], "strand":[], "type":[], "side":[]}
	copyTrueValues = trueVals.copy(deep = True) #make deep copy of dataframe
	numberFailedBoth = 0
	passedBoth = 0
	passedLeftOnly = 0
	passedRightOnly = 0
	total = 0
	for index, row in trueVals.iterrows():
		total += 1
		leftPassed, leftStart, leftEnd, rightPassed, rightStart, rightEnd = shiftLeftAndRightNegativesFilterPositives(row, trueVals, spacingValue, shiftValue, fastaSeq)
		if leftPassed:
			blankDict['seqName'].append(row['seqName'])
			blankDict['start'].append(leftStart)
			blankDict['end'].append(leftEnd)
			newID = row['clusterID'] + "_LeftNegative"
			blankDict['clusterID'].append(newID)
			blankDict['strand'].append(row['strand'])
			blankDict['type'].append(row['type'])
			blankDict['side'].append("Left")
		if rightPassed:
			blankDict['seqName'].append(row['seqName'])
			blankDict['start'].append(rightStart)
			blankDict['end'].append(rightEnd)
			newID = row['clusterID'] + "_RightNegative"
			blankDict['clusterID'].append(newID)
			blankDict['strand'].append(row['strand'])
			blankDict['type'].append(row['type'])
			blankDict['side'].append("Right")
		if leftStart == leftEnd == -1:
			print ("Failed due to N's in left sequence")
		if rightStart == rightEnd == -1:
			print ("Failed due to N's in right sequence")
		if leftStart == leftEnd == rightStart == rightEnd == -2:
			print ("Failed due to N's in true positive sequence")	
		if leftPassed and rightPassed:
			passedBoth += 1
		if leftPassed and not rightPassed:
			passedLeftOnly += 1
		if rightPassed and not leftPassed:
			passedRightOnly += 1
		if not leftPassed and not rightPassed:
			numberFailedBoth += 1
			#delete the row from the copy
			copyTrueValues.drop(index = index, inplace = True)
			print ("DROPPED: ", row['clusterID'])
	print ("total rows: ", total)
	print ("True positives remaining: ", total - numberFailedBoth)
	print ("Number with left and right negatives: ", passedBoth)
	print ("Number with only left negative: ", passedLeftOnly)
	print ("Number with only right negative: ", passedRightOnly)
	#print ("total negatives: ", total)
	asDataFrame = pd.DataFrame(blankDict)
	filteredTrueName = fileName + "BalancedPositives.csv"
	balancedNegatives = fileName + "BalancedNegatives.csv"
	copyTrueValues.to_csv(filteredTrueName)
	asDataFrame.to_csv(balancedNegatives)
	


def openForwardReverse(stem, name):
	totalNameFor = stem + name + ".npy"
	print (totalNameFor)
	forward = np.load(totalNameFor)
	reverse = np.load(stem + name + "RC.npy")
	return forward, reverse

def openBalancedNegatives(name):
	negativesName = name + "BalancedNegatives.csv"
	colNames = ["seqName", "start", "end", "clusterID", "strand", "type", "side"]
	return pd.read_csv( negativesName, names = colnames, dtype = {"seqName": str}) 

def openBalancedPositives(name):
	positivesName = name + "BalancedPositives.csv"
	colNames = ["seqName", "start", "end", "clusterID", "strand", "type", "side"]
	return pd.read_csv( positivesName, names = colnames, dtype = {"seqName": str}) 
	
def extractPredictionValues(name):
	dummy = []
	dummyBools = []
	#contigSeq = SeqIO.read("chrY.fasta", "fasta")
	forward, reverse = openForwardReverse("", "chrY")
	for index, row in negatives.iterrows():
		if row['strand'] == "+":
			#forward strand
			#negatives are 0-indexed still
			dummy = 0
		elif row['strand'] == "-":
			dummy = 0
		else:
			print ("ERROR!  Strand not listed for negative example")
		
#trying to make Negatives:
openedPAS = openPASClustersForChromosome("Y")
#print (openedPAS)
contigSeq = SeqIO.read("chrY.fasta", "fasta")

print ("1: ")
createNegativeDataSet(openedPAS, 1, 50, "chrYNegatives", contigSeq.seq)

#print ("1 w/25:")
#createNegativeDataSet(openedPAS, 1, 25, "chrYNegatives.csv")
#print ("1 w/20")
#createNegativeDataSet(openedPAS, 1, 20, "chrYNegatives.csv")
	




	
		
	
	
	
		
		
	
	


