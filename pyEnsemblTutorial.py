#AML
#1/6/19
#PyEnsembl tutorial - https://www.hammerlab.org/2015/02/04/exploring-the-genome-with-ensembl-and-python/
#already set up local copy of release 96

import pyensembl
from transcriptGraphing import SpliceVariantPASDiagram
from operator import itemgetter
#pulls current release version for the human genome
ensembl = pyensembl.EnsemblRelease()




#find a gene at a location 
print ("Examining 1st gene")
print ( " ")
genes = ensembl.genes_at_locus(contig = "1", position = 1000000)
#print ("first gene: ", genes[0])
hes4 = genes[0]
print (genes)
print ("gene: ", hes4.name)
print ("id: ", hes4.id)
print ("start of gene: ", hes4.start)
print ("end of gene: ", hes4.end)

'''
print (hes4.exons)
transcriptsPos = []
for exon in hes4.exons:
	transcriptsPos.append([exon.start, exon.end])
print (transcriptsPos)
print (transcriptsPos.sort(key = itemgetter(1)))
transcriptGraph = SpliceVariantPASDiagram([transcriptsPos])
transcriptGraph.show()
'''
'''
transcriptsPos = []
transcriptNames = []
for transcript in hes4.transcripts:
	transcriptNames.append(transcript.id + " " + transcript.biotype)
	temp = [[exon.start, exon.end] for exon in transcript.exons]
	transcriptsPos.append(temp[::-1])

transcriptGraph = SpliceVariantPASDiagram(transcriptsPos, transcript_names = transcriptNames, diagramTitle = hes4.name)
transcriptGraph.show()
	'''



print ("number genes on Y chro: ")
gene_ids = ensembl.gene_ids()
genes = [ensembl.gene_by_id(gene_id) for gene_id in gene_ids]
coding_genes = [gene for gene in genes if gene.contig == 'Y']
print (len(coding_genes))

'''
#looking at all genes
print ("Examining all genes")
print (" ")
gene_ids = ensembl.gene_ids()
print ("total number of genes: ", len(gene_ids))
genes = [ensembl.gene_by_id(gene_id) for gene_id in gene_ids]
coding_genes = [gene for gene in genes if gene.biotype == 'protein_coding']
biotypes = [gene.biotype for gene in genes]
biotypes = list(set(biotypes))


print ("Number protien coding genes: ", len(coding_genes))
print ("Biotypes: ")
for bt in biotypes:
	print (bt)
'''
