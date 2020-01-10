#AML
#1/6/19
#PyEnsembl tutorial - https://www.hammerlab.org/2015/02/04/exploring-the-genome-with-ensembl-and-python/
#already set up local copy of release 96

import pyensembl
import pandas as pd
from transcriptGraphing import SpliceVariantPASDiagram, MultiGeneVariantPASDiagram
from Bio import SeqIO

#pulls current release version for the human genome
ensembl = pyensembl.EnsemblRelease()

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

def PASClusterToListFormatSingleStrand(clusters):
	positions = []
	types = []
	for index, row in clusters.iterrows():
		positions.append([row.start, row.end])
		types.append(row.type)
	return positions, types

def PASClusterToListFormatDoubleStrand(clusters):
	positions = [[],[]]
	types = [[],[]]
	for index, row in clusters.iterrows():
		if row.strand == "+":
			positions[0].append([row.start, row.end])
			types[0].append(row.type)
		else:
			positions[1].append([row.start, row.end])
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
		geneStrands.append(g.strand)
		geneAddress.append(g.name + " " + g.contig + ":" + str(g.start) + "-" + str(g.end) + "(" + g.strand + ")")
	return positions, names, geneNames, genePositions, geneStrands, geneAddress

onChrY = openPASClustersForChromosome("Y")
contigSeq = SeqIO.read("chrY.fasta", "fasta")
print ("total length of chromosome Y: ", len(contigSeq.seq))

#position w/ 3 genes overlapping in the annotation 
#finding genes which overlap 





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
for gene in genesOnYChr:
	inGene = findAllPASClustersInRangeOnStrand(gene, onChrY)
	#print ("Numbe PAS clusters: ", inGene.shape[0])
	if inGene.shape[0] > maxPAS:
		maxPAS = inGene.shape[0]
		bestGene = gene.id
print ("highest number in gene: ", maxPAS)
print ("Gene is: ", bestGene)
#graphGeneTranscripts(bestGene, True)
#graphGeneTranscripts(bestGene, False)
#gene = ensembl.gene_by_id(bestGene)
#graphGeneandPAS(gene, onChrY, True)
'''
