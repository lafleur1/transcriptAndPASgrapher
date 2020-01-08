#AML
#1/6/19
#PyEnsembl tutorial - https://www.hammerlab.org/2015/02/04/exploring-the-genome-with-ensembl-and-python/
#already set up local copy of release 96

import pyensembl
import pandas as pd
from transcriptGraphing import SpliceVariantPASDiagram

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
def graphGeneTranscripts(geneId):
	gene = ensembl.gene_by_id(geneId)
	graphTitle = gene.name + " " + gene.contig + ":" + str(gene.start) + "-" + str(gene.end) + "(" + gene.strand + ")"
	transcriptsPos = []
	transcriptNames = []
	for transcript in gene.transcripts:
		transcriptNames.append(transcript.id + " " + transcript.biotype)
		temp = [[exon.start, exon.end] for exon in transcript.exons]
		transcriptsPos.append(temp[::-1])
	transcriptGraph = SpliceVariantPASDiagram(transcriptsPos, transcript_names = transcriptNames, diagramTitle = graphTitle)
	transcriptGraph.show()

#filters for PAS Clusters in a given gene range
def findAllPASClustersInRangeOnStrand(gene, chromosomePAS):
	maskPAS = (chromosomePAS['start'] > gene.start - 1) & (chromosomePAS['end'] < gene.end+ 1) & (chromosomePAS['strand'] == gene.strand)
	maskedPAS = chromosomePAS[maskPAS]
	return maskedPAS
	
	
def PASClusterToListFormat(clusters):
	positions = []
	types = []
	for index, row in clusters.iterrows():
		positions.append([row.start, row.end])
		types.append(row.type)
	return positions, types


def graphGeneandPAS(gene, chromosomePAS):
	inGene = findAllPASClustersInRangeOnStrand(gene, chromosomePAS)
	pasPositions, pasTypes = PASClusterToListFormat(inGene)
	graphTitle = gene.name + " " + gene.contig + ":" + str(gene.start) + "-" + str(gene.end) + "(" + gene.strand + ")"
	transcriptsPos = []
	transcriptNames = []
	print ("number transcripts: ", len(gene.transcripts))
	for transcript in gene.transcripts:
		transcriptNames.append(transcript.id + " " + transcript.biotype)
		temp = [[exon.start, exon.end] for exon in transcript.exons]
		transcriptsPos.append(temp[::-1])
	transcriptGraph = SpliceVariantPASDiagram(transcriptsPos, transcript_names = transcriptNames, pas_pos= pasPositions, pas_types = pasTypes, diagramTitle = graphTitle)
	transcriptGraph.show()
	



onChrY = openPASClustersForChromosome("Y")

notFoundTwo = True
posCurrent = 2786855
while notFoundTwo:
	genesDemo = ensembl.genes_at_locus(contig = "Y", position = posCurrent)
	if len(genesDemo) > 2:
		notFoundTwo = False
	else:
		posCurrent += 1

print (posCurrent)
print (genesDemo)		
#print (genesDemo)
#findAllPASClustersInRangeOnStrand(genesDemo[0], onChrY)
#graphGeneandPAS(genesDemo[0], onChrY)
