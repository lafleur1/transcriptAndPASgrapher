#extracting sequences near the true values in same regions for true negative datasets for APA model evaluation
#APARENT: 206 nt, overlapping
#DeepPASTA: 200 nts, polyA site in middle
#
#using same method as Leung et al (see  Inference of the human polyadenylation code supplementary material):
 #same padding requirement (i.e. that spacing between two PAS signals must allow four negative regions 
 #shift left and right by 50 nts but still within the region boundaries 
 
import pyensembl
import pandas as pd 

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


def shiftLeftAndRightNegatives(pasRow, allPAS, spacing, shiftValue):
	width = pasRow['end'] - pasRow['start']
	checkRange = width + shiftValue
	leftStart = pas_start - shiftValue
	leftEnd = pas_end - shiftValue
	leftPassed = False
	rightStart = pas_start - shiftValue
	rightEnd = pas_end - shiftValue
	rightPassed = False
	#checking left
	filterVals = (allPAS['strand'] == pasRow['strand']) & (allPAS['seqName'] == pasRow['seqName']) & (allPAS['end'] >= pasRow['start'] - checkRange) & (allPAS['end'] <= pasRow['end'])
	otherPAS = allPAS[filterVals]
	if otherPAS.shape[0] == 1: 
		#left false is okay
		leftPassed = True
	filterVals = (allPAS['strand'] == pasRow['strand']) & (allPAS['seqName'] == pasRow['seqName']) & (allPAS['start'] <= pasRow['end'] + checkRange) & (allPAS['start'] >= pasRow['start'])
	otherPAS = allPAS[filterVals]
	if otherPAS.shape[0] == 1:
		#passed 
		rightPassed = True
	return leftPassed, leftStart, leftEnd, rightPassed, rightStart, rightEnd



	
		
	
	
	
		
		
	
	


