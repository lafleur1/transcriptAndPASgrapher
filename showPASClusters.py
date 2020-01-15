#AML
#1/6/19
#PyEnsembl tutorial - https://www.hammerlab.org/2015/02/04/exploring-the-genome-with-ensembl-and-python/
#already set up local copy of release 96

import pyensembl
import pandas as pd
from transcriptGraphing import SpliceVariantPASDiagram, MultiGeneVariantPASDiagram
from Bio import SeqIO
import numpy as np
from scipy.signal import find_peaks

#pulls current release version for the human genome
ensembl = pyensembl.EnsemblRelease(release = '96')

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
			
	
	
#creates gene graph w/out PAS markers for a strand 
def graphGeneTranscripts(geneId, thicken):
	gene = ensembl.gene_by_id(geneId)
	graphTitle = gene.name + " " + gene.contig + ":" + str(gene.start) + "-" + str(gene.end) + "(" + gene.strand + ")"
	transcriptsPos = []
	transcriptNames = []
	for transcript in gene.transcripts:
		transcriptNames.append(transcript.id + " " + transcript.biotype)
		temp = [[exon.start, exon.end] for exon in transcript.exons]
		transcriptsPos.append(temp[::-1])
	transcriptGraph = SpliceVariantPASDiagram(transcriptsPos, transcript_names = transcriptNames, diagramTitle = graphTitle, thickenExons = thicken)
	transcriptGraph.show()



#filters for PAS Clusters in a given gene range
def findAllPASClustersInGene(gene, chromosomePAS):
	maskPAS = (chromosomePAS['start'] > gene.start - 1) & (chromosomePAS['end'] < gene.end+ 1) & (chromosomePAS['strand'] == gene.strand)
	maskedPAS = chromosomePAS[maskPAS]
	return maskedPAS
	

def findallPASClustersInRangeBothStrands(start, stop, chromosomePAS):
	maskPAS = (chromosomePAS['start'] > start - 1) & (chromosomePAS['end'] < stop+ 1)
	maskedPAS = chromosomePAS[maskPAS]
	return maskedPAS


#ENSEMBL is 1-based indexing
#PolyASite2.0 is 0-based indexing
#change polyAsite to Ensembl notation here
def correctPolyASiteIndexing(pasRange):
	return [pasRange[0]+1, pasRange[1]+1]
	


#correct indexing here to Ensembl 1-based
def PASClusterToListFormatSingleStrand(clusters):
	positions = []
	types = []
	for index, row in clusters.iterrows():
		positions.append(correctPolyASiteIndexing([row.start, row.end]))
		types.append(row.type)
	return positions, types

#correct indexing here to Ensembl 1-based
def PASClusterToListFormatDoubleStrand(clusters):
	positions = [[],[]]
	types = [[],[]]
	for index, row in clusters.iterrows():
		if row.strand == "+":
			positions[0].append(correctPolyASiteIndexing([row.start, row.end]))
			types[0].append(row.type)
		else:
			positions[1].append(correctPolyASiteIndexing([row.start, row.end]))
			types[1].append(row.type)
	return positions, types


def graphGeneandPAS(gene, chromosomePAS, dropdowns):
	inGene = findAllPASClustersInRangeOnStrand(gene, chromosomePAS)
	pasPositions, pasTypes = PASClusterToListFormatSingleStrand(inGene)
	graphTitle = gene.name + " " + gene.contig + ":" + str(gene.start) + "-" + str(gene.end) + "(" + gene.strand + ")"
	transcriptsPos = []
	transcriptNames = []
	print ("number transcripts: ", len(gene.transcripts))
	for transcript in gene.transcripts:
		transcriptNames.append(transcript.id + " " + transcript.biotype)
		temp = [[exon.start, exon.end] for exon in transcript.exons]
		transcriptsPos.append(temp[::-1])
	transcriptGraph = SpliceVariantPASDiagram(transcriptsPos, transcript_names = transcriptNames, pas_pos= pasPositions, pas_types = pasTypes, diagramTitle = graphTitle, dropDownPASMarkers = dropdowns)
	transcriptGraph.show()


def extractTotalSpanFromGenePositions(gps):
	total_span_beginning = float('inf')
	total_span_end = float('-inf')
	for g in gps:
		if g[0] <= total_span_beginning:
			total_span_beginning = g[0]
		if g[1] >= total_span_end:
			total_span_end = g[1]
	return total_span_beginning, total_span_end

def extractTotalSpanFromGenes(gps):
	total_span_beginning = float('inf')
	total_span_end = float('-inf')
	for g in gps:
		if g[0] <= total_span_beginning:
			total_span_beginning = g[0]
		if g[1] >= total_span_end:
			total_span_end = g[1]
	return total_span_beginning, total_span_end


def transformPyEnsemblToPASDiagram(geneList):
	positions = []
	names = []
	geneNames = []
	genePositions = []
	geneStrands = []
	geneAddress = []
	for g in geneList:
		dummyPos = [] #store all transcripts for the gene
		dummyNames = []
		for t in g.transcripts:
			dummyPos.append([[exon.start, exon.end] for exon in t.exons])
			dummyNames.append(t.id + " " + t.biotype)
		positions.append(dummyPos)
		names.append(dummyNames)
		geneNames.append(g.id + ": " + g.name + " "+ g.biotype)
		genePositions.append([g.start,g.end])
		#print (g.start)
		#print (g.end)
		geneStrands.append(g.strand)
		geneAddress.append(g.name + " " + g.contig + ":" + str(g.start) + "-" + str(g.end) + "(" + g.strand + ")")
	return positions, names, geneNames, genePositions, geneStrands, geneAddress


def openForwardReverse(stem, name):
	totalNameFor = stem + name + ".npy"
	print (totalNameFor)
	forward = np.load(totalNameFor)
	reverse = np.load(stem + name + "RC.npy")
	return forward, reverse


#use to correct the indexes for the RC predictions so that they can be matched with the true cluster labels (since those are indexed according to the forward strand)
def flipSequenceIndex(index, lenSeq):
	b = lenSeq - 1
	return -1 * index + b


def find_peaks_ChromosomeVersion(avgPreds, peak_min_height, peak_min_distance, peak_prominence):
	#forwardPeaks = find_peaks_ChromosomeVersion(forward, minh, dist, (0.01, None)) 
	peaks, _ = find_peaks(avgPreds, height=peak_min_height, distance=peak_min_distance, prominence=peak_prominence) 
	return peaks

onChrY = openPASClustersForChromosome("Y")
contigSeq = SeqIO.read("chrY.fasta", "fasta")
print ("total length of chromosome Y: ", len(contigSeq.seq))
forwardY, reverseY = openForwardReverse("", "chrY")
print ("size predictions: ", forwardY.size)
reverseYFlipped = np.flip(reverseY)
#examining areas with overlapping genes on chromosome Y
notFoundTwo = True
posCurrent = 0

while notFoundTwo:
	genesDemo = ensembl.genes_at_locus(contig = "Y", position = posCurrent)
	if len(genesDemo) >= 1:
		#print (posCurrent)
		#print (genesDemo)		
		p,tn, gn, gp, gs, ga = transformPyEnsemblToPASDiagram(genesDemo)
		start, stop = extractTotalSpanFromGenePositions(gp)
		#print ("width of slice ", stop - start)
		pasTable = findallPASClustersInRangeBothStrands(start, stop, onChrY)
		genesRound2 = ensembl.genes_at_locus(contig = "Y", position = start, end = stop)
		p,tn, gn, gp, gs, ga = transformPyEnsemblToPASDiagram(genesRound2)
		start, stop = extractTotalSpanFromGenePositions(gp)
		#start and stop are 1-indexed, sequence and predictions are 0-indexed
		seqSlice = contigSeq.seq[start-1:stop]
		#print ("len seq: ", len(seqSlice))
		predForwardSlice = forwardY[start-1:stop] #numpt array
		#print ("forward: ", predForwardSlice.size)
		predReverseSlice =  reverseYFlipped[start-1 :stop]
		#print ("reverse: ", predReverseSlice.size)
		#print (genesRound2)
		#print (pasTable)
		#print ("Using min H of 0.5")
		preppedPASLocs, preppedPASTypes = PASClusterToListFormatDoubleStrand(pasTable)
		peaksForward = find_peaks_ChromosomeVersion(predForwardSlice, 0.5, 50, (0.01, None))
		peaksReverse  = find_peaks_ChromosomeVersion(predReverseSlice, 0.5, 50, (0.01, None))
		#print ('forward peaks: ', peaksForward)
		#print ('peaksReverse: ', peaksReverse)
		if len(peaksForward) > 0 or len(peaksReverse) > 0 :
			mgraph = MultiGeneVariantPASDiagram(p, tn, gn, gp,  gs, pas_pos= preppedPASLocs, pas_types = preppedPASTypes, sequence = seqSlice, forwardValues = predForwardSlice, reverseValues = predReverseSlice, forwardPeaks = peaksForward, reversePeaks = peaksReverse)
			mgraph.show()
		posCurrent += (stop - start)
	else:
		posCurrent += 1
	if posCurrent >= len(contigSeq.seq):
		notFoundTwo = False 




'''
#examining IG pas's on chromosome Y
igPAS = onChrY[onChrY['type'] == 'IG']
for i, row in igPAS.iterrows():
	#graph genes around first AU PAS
	print (row)
	start = int(row['start']) - 10000
	stop = int(row['end']) + 10000
	genesRound2 = ensembl.genes_at_locus(contig = "Y", position = start, end = stop)
	#print (genesRound2)
	for g in genesRound2:
	#	print (g.name)
		for t in g.transcripts:
			print (t.name)
			if t.contains_start_codon:
				print (t.start_codon_positions)
			else:
				print ("No start codon annotation")
			
	p,tn, gn, gp, gs, ga = transformPyEnsemblToPASDiagram(genesRound2)
	pasPos, pasType = PASClusterToListFormatDoubleStrand(igPAS[igPAS['clusterID'] == row['clusterID']])
	#print (pasPos, pasType)
	mgraph = MultiGeneVariantPASDiagram(p, tn, gn, gp,  gs, pas_pos= pasPos, pas_types = pasType, startOverride = start, stopOverride = stop)
	mgraph.show()
'''

'''
#examining AU pases on chromosome Y
auPAS = onChrY[onChrY['type'] == 'AU']
for i, row in auPAS.iterrows():
	#graph genes around first AU PAS
	print (row)
	start = int(row['start']) - 2000
	stop = int(row['end']) + 2000
	genesRound2 = ensembl.genes_at_locus(contig = "Y", position = start, end = stop)
	#print (genesRound2)
	for g in genesRound2:
	#	print (g.name)
		for t in g.transcripts:
			print (t.name)
			if t.contains_start_codon:
				print (t.start_codon_positions)
			else:
				print ("No start codon annotation")
			
	p,tn, gn, gp, gs, ga = transformPyEnsemblToPASDiagram(genesRound2)
	pasPos, pasType = PASClusterToListFormatDoubleStrand(auPAS[auPAS['clusterID'] == row['clusterID']])
	#print (pasPos, pasType)
	mgraph = MultiGeneVariantPASDiagram(p, tn, gn, gp,  gs, pas_pos= pasPos, pas_types = pasType, startOverride = start, stopOverride = stop)
	mgraph.show()

'''



'''
#examining areas with overlapping genes on chromosome Y
notFoundTwo = True
posCurrent = 6900000
while notFoundTwo:
	genesDemo = ensembl.genes_at_locus(contig = "Y", position = posCurrent)
	if len(genesDemo) > 1:
		#notFoundTwo = False
		print (posCurrent)
		print (genesDemo)		
		#print (genesDemo)
		#findAllPASClustersInRangeOnStrand(genesDemo[0], onChrY)
		p,tn, gn, gp, gs, ga = transformPyEnsemblToPASDiagram(genesDemo)
		start, stop = extractTotalSpanFromGenePositions(gp) 
		pasTable = findallPASClustersInRangeBothStrands(start, stop, onChrY)
		genesRound2 = ensembl.genes_at_locus(contig = "Y", position = start, end = stop)
		p,tn, gn, gp, gs, ga = transformPyEnsemblToPASDiagram(genesRound2)
		print (genesRound2)
		print (pasTable)
		preppedPASLocs, preppedPASTypes = PASClusterToListFormatDoubleStrand(pasTable)
		#print (preppedPASLocs)
		#print (preppedPASTypes)
		mgraph = MultiGeneVariantPASDiagram(p, tn, gn, gp,  gs, pas_pos= preppedPASLocs, pas_types = preppedPASTypes)
		mgraph.show()
		posCurrent += (stop - start)
	else:
		posCurrent += 1
	if posCurrent >= len(contigSeq.seq):
		notFoundTwo = False 

'''



'''

#errors in database???
id1 = "Y:6919894:+"
genes = ensembl.genes_at_locus(contig = "Y", position = 6919894)
print (genes)
id2 = "Y:6920421:+"
genes = ensembl.genes_at_locus(contig = "Y", position = 6919894)
print (genes)
id3 = "Y:6920474:+"
genes = ensembl.genes_at_locus(contig = "Y", position = 6919894)
print (genes)
id4 = "Y:7102025:+"
genes = ensembl.genes_at_locus(contig = "Y", position = 6919894)
print (genes)
#graphGeneandPAS(genesDemo[0], onChrY)
'''
#finding gene with most PAS signals contained
'''
gene_ids = ensembl.gene_ids()
genes = [ensembl.gene_by_id(gene_id) for gene_id in gene_ids]
genesOnYChr = [gene for gene in genes if gene.contig == 'Y']
maxPAS = float('-inf')
bestGene = ""
totalTranscripts = 0
totalTranscriptsWithStart = 0
for gene in genesOnYChr:
	#inGene = findAllPASClustersInRangeOnStrand(gene, onChrY)
	#print ("Numbe PAS clusters: ", inGene.shape[0])
	#if inGene.shape[0] > maxPAS:
	#	maxPAS = inGene.shape[0]
	#	bestGene = gene.id
	for t in gene.transcripts:
		totalTranscripts += 1
		if t.contains_start_codon:
			print (t)
			totalTranscriptsWithStart += 1
print (totalTranscripts)
print (totalTranscriptsWithStart)
#print ("highest number in gene: ", maxPAS)
#print ("Gene is: ", bestGene)
#graphGeneTranscripts(bestGene, True)
#graphGeneTranscripts(bestGene, False)
#gene = ensembl.gene_by_id(bestGene)
#graphGeneandPAS(gene, onChrY, True)
'''
