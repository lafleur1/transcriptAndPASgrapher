#AML
#opening gencode release 32 GTF file
#instructions from: https://biopython.org/wiki/GFF_Parsing

from transcriptGraphing import UCSCAnnotationDiagram
import pandas as pd

annotname = "gencode.v32.annotation.gtf"
fullAnnotName = "gencode.v32.chr_patch_hapl_scaff.annotation.gtf"

#ucsc data uses 0-based indexing

ucsc_wholeGene = "annotatedhg38" 
ucsc_upstream = "annotatedhg38Upstream"
ucsc_downstream = "annotatedhg38Downstream"
ucsc_allExons = "annotatedhg38Exons"
ucsc_allIntrons = "annotatedhg38Introns"
ucsc_3UTRExons = "annotatedhg38_3_UTRExons" #last exon
ucsc_5UTRExons= "annotatedhg38TerminalExons" #first exon
ucsc_CDSExons = "annotatedhg38CodingExons"


#filters for PAS Clusters in a given gene range
def findAllPASClustersInGene(start, stop, contig, strand, chromosomePAS):
	maskPAS = (chromosomePAS['start'] > start - 10000) & (chromosomePAS['end'] < stop+ 10000) & (chromosomePAS['strand'] == strand) & (chromosomePAS['chromosome'] == contig)
	maskedPAS = chromosomePAS[maskPAS]
	return maskedPAS

def findallPASClustersInRangeBothStrands(start, stop, chromosomePAS):
	maskPAS = (chromosomePAS['start'] > start - 1) & (chromosomePAS['end'] < stop+ 1)
	maskedPAS = chromosomePAS[maskPAS]
	return maskedPAS


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

def openAllForGeneName(gene_name):
	bedNames = ["chromosome", "start", "end", "name", "score", "strand"]
	bedNamesWholeGene = ["chromosome", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "itemRGB", "blockCount", "blockSizes", "blockStarts"]
	annotations =pd.read_csv(ucsc_wholeGene,delimiter='\t', names = bedNamesWholeGene ) #ignores comment lines at beginning of the GTF file'
	print ("WHOLE GENE")
	whole_gene = annotations[annotations["name"] == gene_name] 
	print (annotations[annotations["name"] == gene_name])
	print ( "---------------------------------")
	print ("UPSTREAM")
	annotations =pd.read_csv(ucsc_upstream,delimiter='\t', names = bedNames)
	print (annotations[annotations["name"].str.contains(gene_name)])
	upstream = annotations[annotations["name"].str.contains(gene_name)]
	print ( "---------------------------------")
	print ("DOWNSTREAM")
	annotations =pd.read_csv(ucsc_downstream,delimiter='\t', names = bedNames)
	print (annotations[annotations["name"].str.contains(gene_name)])
	downstream = annotations[annotations["name"].str.contains(gene_name)]
	#types = list(set(annotations.feature_type))
	#print (types)
	print ( "---------------------------------")
	print ("ALL EXONS")
	annotations =pd.read_csv(ucsc_allExons,delimiter='\t', names = bedNames)
	print (annotations[annotations["name"].str.contains(gene_name)])
	exons = annotations[annotations["name"].str.contains(gene_name)]
	print ( "---------------------------------")
	print ("ALL INTRONS:")
	annotations =pd.read_csv(ucsc_allIntrons,delimiter='\t', names = bedNames)
	print (annotations[annotations["name"].str.contains(gene_name)])
	introns = annotations[annotations["name"].str.contains(gene_name)]
	print ( "---------------------------------")
	print ("CODING EXONS:")
	annotations =pd.read_csv(ucsc_CDSExons,delimiter='\t', names = bedNames)
	print (annotations[annotations["name"].str.contains(gene_name)])
	coding_exons = annotations[annotations["name"].str.contains(gene_name)]
	print ( "---------------------------------")
	print ("3' UTR EXONS:")
	annotations =pd.read_csv(ucsc_3UTRExons,delimiter='\t', names = bedNames)
	print (annotations[annotations["name"].str.contains(gene_name)])
	three_prime_utr = annotations[annotations["name"].str.contains(gene_name)]
	print ( "---------------------------------")
	print ("5' UTR EXONS:")
	annotations =pd.read_csv(ucsc_5UTRExons,delimiter='\t', names = bedNames)
	print (annotations[annotations["name"].str.contains(gene_name)])
	five_prime_utr = annotations[annotations["name"].str.contains(gene_name)]
	return whole_gene, upstream, downstream, exons, introns, coding_exons, three_prime_utr, five_prime_utr


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

#colnames = ["chromosome_name", "annot_source", "feature_type", "start", "end", "score", "strand", "phase", "additional_info"]
#annotations =pd.read_csv(ucsc_Stuff,delimiter='\t', names = colnames, comment = "#") #ignores comment lines at beginning of the GTF file

bedNamesWholeGene = ["chromosome", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "itemRGB", "blockCount", "blockSizes", "blockStarts"]
annotations =pd.read_csv(ucsc_wholeGene,delimiter='\t', names = bedNamesWholeGene )
print (findAllPASClustersInGene(6900000, 7200000, "chrY", "+", annotations))

chrYPAS = openPASClustersForChromosome("Y")
ypas = findallPASClustersInRangeBothStrands(6900000, 7200000, chrYPAS)
ypasPos, ypasTypes = PASClusterToListFormatDoubleStrand(ypas)

#pas_pos=[], pas_types = []

#translation start site -> translation end site
whole_gene, upstream, downstream, exons, introns, coding_exons, three_prime_utr, five_prime_utr = openAllForGeneName("NM_134258")
#attempting to graph gene with messed up TE pas signals 
diagram = UCSCAnnotationDiagram( whole_gene, upstream, downstream, exons, introns, coding_exons, three_prime_utr, five_prime_utr, pas_pos = ypasPos, pas_types = ypasTypes, diagramTitle ="NM_134258")
diagram.show()

whole_gene, upstream, downstream, exons, introns, coding_exons, three_prime_utr, five_prime_utr = openAllForGeneName("NM_134259")
#attempting to graph gene with messed up TE pas signals 
diagram = UCSCAnnotationDiagram( whole_gene, upstream, downstream, exons, introns, coding_exons, three_prime_utr, five_prime_utr,pas_pos = ypasPos, pas_types = ypasTypes, diagramTitle ="NM_134259")
diagram.show()

whole_gene, upstream, downstream, exons, introns, coding_exons, three_prime_utr, five_prime_utr = openAllForGeneName("NM_033284")
#attempting to graph gene with messed up TE pas signals 
diagram = UCSCAnnotationDiagram( whole_gene, upstream, downstream, exons, introns, coding_exons, three_prime_utr, five_prime_utr,pas_pos = ypasPos, pas_types = ypasTypes, diagramTitle ="NM_033284")
diagram.show()




