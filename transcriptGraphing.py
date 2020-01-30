import matplotlib.pyplot as plt
import matplotlib.colors as clr
import numpy as np
from matplotlib.patches import Patch, Circle
from operator import itemgetter
import random
import pandas as pd





class SpliceVariantPASDiagram(object):
	
	def __init__(self, transcripts, strands = [], transcript_names = [],  pas_pos=[], pas_types = [], marker_heights=[], marker_size=50, marker_weight=1, exon_colors= [], intron_color="gray", intron_weight=1, intron_style='-', bar_color='gray', bg_color="white", diagramTitle = "", thickenExons = True, dropDownPASMarkers = False, modelPredictions = []):
		###transcripts
		self.numberTranscripts = len(transcripts) #number of splice variants on Ensembl for the gene 
		self.transcripts = transcripts 
		self.strands = strands 
		self.transcript_names = transcript_names
		#if colors are provided, use those.  If none are specified, set all exons to black rectangles
		if exon_colors:
			self.exonColors = exon_colors
		else:
			self.exonColors = ["black"] * self.numberTranscripts
		self.numExons = [ len(x) for x in self.transcripts] #number of exons in each transcript
		spans = []
		#find smallest and largest positions 
		smallestPosition = float('inf')
		largestPosition = float('-inf')
		for exons in self.transcripts:
			spans.append(exons[-1][1] - exons[0][0])
			if exons[-1][1] >= largestPosition:
				largestPosition = exons[-1][1]
			if exons[0][0] <= smallestPosition:
				smallestPosition = exons[0][0]
		self.smallestIndex = smallestPosition
		self.largestIndex = largestPosition
		self.totalSpan =  max(spans)
		###PAS markers
		self.pasPositions = pas_pos #marker positions to display (list of [start,stop]) for the PAS cluster data
		self.markerHeights = marker_heights 
		self.pasTypes = pas_types #marker color depends on PAS cluster type
		self.markerSize = marker_size
		self.MarkerWeight = marker_weight
		self.colorKey = {'TE': "red", "EX": "salmon", "IN": "deepskyblue", "DS": "forestgreen", "AE": "darkred", "AI": "midnightblue", "AU": "limegreen", "IG": "gold"} #colors for the 8 PAS cluster types 
		#model predictions
		self.modelPredictions = modelPredictions #numpy array of model predictions on that slice
		
		###graph stuff
		self.dropDownPASMarkers = dropDownPASMarkers
		self.thickenExons = thickenExons #if True, adjusts exon ranges so that very small ones still appear in overall graph.  If false, does not adjust
		self.diagramTitle = diagramTitle
		self.intronColor = intron_color
		self.intronWeight = intron_weight
		self.intronStyle = intron_style
		self.barColor= bar_color
		self.bgColor = bg_color
		self.minExonLen = self.totalSpan*0.005 #set a mimimum exon length to show when graphing
		self.ylims = {'exon_max': 2, 'exon_min':1} #max and min for graph exon bars
		self.figure, self.canvas = plt.subplots(figsize=(15,1.5)) 
		self.canvas.set_facecolor(self.bgColor) 
		self._draw() #draw the graph 

	#setting the other y mins and maxes for the intron line points and the overhead bar (exon min/max heights are set above)
	def _set_limits(self):
		self.ylims['intron_max'] = self.ylims['exon_max']*0.9
		self.ylims['intron_min'] = (self.ylims['exon_max'] + self.ylims['exon_min'])/2.0
		self.ylims['bar_min'] = self.ylims['exon_max']+0.2 #setting the height of the overhead bar
		self.ylims['bar_max'] = self.ylims['bar_min']+(self.ylims['exon_max']-self.ylims['exon_min'])/5.0 #setting the height of the overhead bar

	#if exons exist with a lenght lower than the min exon length to draw, expand the covered interval to be visible
	def _transform_spans(self):
		for index, exons in enumerate(self.transcripts): #transform each set of exons seperately
			span_lens = [x[1]-x[0] for x in exons] #find all exon lengths in the given sequence
			max_len = float(max(span_lens)) #find the max length exon
			transformed_intervals = []
			if max_len < self.minExonLen: #if the smallest exon is lower than the minimum exon lenght expand the length it covers to be visible
				span_ratios = [x/max_len for x in span_lens]
				expansion_factor = self.totalSpan*1e-11
				for i in range(1,10):
					ef = (2**i)*expansion_factor
					if max_len+ef > self.minExonLen:
						expansion_factor = ef
						break
				for i,j in zip(exons, span_ratios):
					mid = (i[0] + i[1])/2
					f = (expansion_factor*j)/2
					if mid+f - mid-f > self.minExonLen:
						transformed_intervals.append([mid-f, mid+f])
					else:
						transformed_intervals.append([mid-(self.minExonLen/2), mid+(self.minExonLen/2)])
			else:
				for i in range(self.numExons[index]):
					if span_lens[i] < self.minExonLen: #expand if hte exon is lower than the mimimum 
						mid = (exons[i][0] + exons[i][0])/2 
						transformed_intervals.append([mid-(self.minExonLen/2), mid+(self.minExonLen/2)])
					else: #do nothing if it is larger than the minimal exon length
						transformed_intervals.append(exons[i])
			self.transcripts[index] = transformed_intervals #overwrite exon lengths for the graph

	def _draw_exon(self, span, offset, currentExonColor): #draw a rectangle between the two limits in the the span 
		self.canvas.fill_between(span, self.ylims['exon_min'] - offset, self.ylims['exon_max'] - offset, edgecolor=self.bgColor, facecolor=currentExonColor)
		return True  #if successful in drawing the exon, return true

	def _draw_intron(self, span, offset): #span is between last exon end point x and next exon's starting point 
		mid = (span[0]+span[1])/2.0 #find midpoint between the two exons 
		#make two line segments for the intron 
		self.canvas.plot([span[0], mid], [self.ylims['intron_min'] - offset, self.ylims['intron_max'] - offset], c=self.intronColor, lw=self.intronWeight, ls=self.intronStyle)
		self.canvas.plot([mid, span[1]], [self.ylims['intron_max'] - offset, self.ylims['intron_min'] - offset], c=self.intronColor, lw=self.intronWeight, ls=self.intronStyle)
		return True #if successful in drawing the intron, return true
    
	def _draw_markers(self):
		if self.pasPositions:
			if self.markerHeights == []:
				self.markerHeights = [self.ylims['exon_max']-self.ylims['exon_min'] for x in self.pasPositions] #make the markers the same height as the exon blocks          
			if self.dropDownPASMarkers:
				#drop PAS marker lines down through transcript graphs
				markerBottom =self.ylims['exon_min'] - ( len(self.transcripts) * 2)
			else:
				markerBottom = self.ylims['bar_max']
			for p,h,t in zip(self.pasPositions, self.markerHeights, self.pasTypes): #for each marker make a bar and then place a scatter plot dot on the bar
				markerColor = self.colorKey[t] #retrieve color for the pas type
				if type(p) == type(1):
					self.canvas.plot((p, p), (markerBottom, self.ylims['bar_max']+h), linestyle='-', color=markerColor, linewidth=self.MarkerWeight, alpha=0.7)
					self.canvas.scatter(p, self.ylims['bar_max']+h+0.25, s=self.markerSize, marker='o', c=markerColor, edgecolor=markerColor, alpha=1)
				elif type(p) == type([]) and len(p) == 2:
					midpoint = (p[0] + p[1] )/ 2.0
					self.canvas.fill_between(p, markerBottom, self.ylims['bar_max'] + h, edgecolor= markerColor, facecolor= markerColor)
					self.canvas.scatter(midpoint, self.ylims['bar_max']+h+0.25, s=self.markerSize, marker='o', c=markerColor, edgecolor=markerColor, alpha=1)
			#if there are PAS markers add a legend for the types
			types = [Patch(facecolor='red', edgecolor='red', label='TE'), 
					Patch(facecolor='salmon', edgecolor='salmon', label='EX'), 
					Patch(facecolor='deepskyblue', edgecolor='deepskyblue', label='IN'), 
					Patch(facecolor='forestgreen', edgecolor='forestgreen', label='DS'), 
					Patch(facecolor='darkred', edgecolor='darkred', label='AE'), 
					Patch(facecolor='midnightblue', edgecolor='midnightblue', label='AI'), 
					Patch(facecolor='limegreen', edgecolor='limegreen', label='AU'), 
					Patch(facecolor='gold', edgecolor='gold', label='IG')]
			l1 = plt.legend(handles=types, loc = 'center right', bbox_to_anchor = (0,0.5))
			canvas.add_artist(l1)
			#ax.legend(loc='center left', bbox_to_anchor=(1, 0.5)))

        
	#place markers between smallest/largest
	def _clean_axes(self):
		self.canvas.set_yticks([], []) #hides all y-axis ticks
		self.canvas.get_xaxis().tick_top() 
		self.canvas.tick_params(axis='x', direction='out') 
		self.canvas.set_xticks([]) #remove x ticks 
		for o in ["top", "bottom", "left", "right"]: #removing the box around the graph 
			self.canvas.spines[o].set_visible(False)
		min_pos = int(self.smallestIndex - self.totalSpan * 0.1) 
		if min_pos < 0:
			min_pos = 0
		max_pos = int(self.largestIndex + self.totalSpan * 0.1) 
		stepSize = int((max_pos-min_pos)/20) #force stepsize to be an int for placing range markers 
		minortick_pos = [x for x in range(min_pos, max_pos, stepSize)][1:] #divide up the total range by 20 and place tick marks 
		for i in minortick_pos:
			self.canvas.axvline(i, alpha=0.2, c='black', ls='--') #add vertical line across the axes 
		self.canvas.text(minortick_pos[0], self.ylims['exon_min']-0.5, minortick_pos[0], fontsize=8, ha='center') #place left value marker
		self.canvas.text(minortick_pos[-1], self.ylims['exon_min']-0.5, minortick_pos[-1], fontsize=8, ha='center') #place right value marker
		self.canvas.set_xlim(minortick_pos[0]-(minortick_pos[1]-minortick_pos[0]), minortick_pos[-1]+(minortick_pos[-1]-minortick_pos[-2])) #limit graph x-value range 
        
	def _draw(self):
		self._set_limits() #setting y values for the intron lines and the overhead bar in the graph
		if self.thickenExons:
			self._transform_spans() #expanding very small exons to show on the visual for the region
		for j in range(0,len(self.transcripts)):
			offset = (j ) * 2
			currentExonColor = self.exonColors[j]
			for i in range(self.numExons[j]):
				if i > 0: #if its not the first exon, draw the intron first
					self._draw_intron([self.transcripts[j][i-1][1], self.transcripts[j][i][0]], offset) #between last exon end point and current exon start
				self._draw_exon(self.transcripts[j][i], offset, currentExonColor) #draw first exon or following exons
				if i == self.numExons[j] - 1 and self.transcript_names: #to the right of the last exon place the transcript name if on exists
					tagPosition = self.ylims['exon_max'] - offset
					self.canvas.text( self.transcripts[j][i][1] + 2, tagPosition, self.transcript_names[j], fontsize = 8, va = "top", ha = "left")				
		self.canvas.fill_between([self.smallestIndex, self.largestIndex], self.ylims['bar_min'], self.ylims['bar_max'], edgecolor=self.bgColor, facecolor=self.barColor)
		if self.diagramTitle:
			plt.title(self.diagramTitle)
		self._draw_markers() #draw markers with the default color for the diagram instance
		self._clean_axes()

	def show(self):
		plt.show()



class MultiGeneVariantPASDiagram(object):
	#similar to above, but change so it graphs all the transcripts for a set of genes on +/- strands at once
	#and/or just graphs multiple sets of genes
	
	#extracts needed info from list of PyEnsembl genes 
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
			geneNames.append(g.id + ": " + g.name)
			genePositions.append([g.start,g.end])
			geneStrands.append(g.strand)
			geneAddress.append(g.name + " " + g.contig + ":" + str(g.start) + "-" + str(g.end) + "(" + g.strand + ")")
			if g.end >= total_span_end:
				total_span_end = g.end
			if g.start <= total_span_beginning:
				total_span_beginning = g.start
		return positions, names, geneNames, genePositions, geneStrands, geneAddress, total_span_beginning, total_span_end
			
	
	def extractTotalSpanFromGenePositions(self, gps):
		total_span_beginning = float('inf')
		total_span_end = float('-inf')
		for g in gps:
			#print (g)
			#print (g[0])
			#print (g[1])
			if g[0] <= total_span_beginning:
				total_span_beginning = g[0]
			if g[1] >= total_span_end:
				total_span_end = g[1]
		return total_span_beginning, total_span_end
	
	def _generateRandomColors(self, numberColors):
		#using matplotlib's built in color mapping dict, randomnly select named colors to use
		number_done = 0
		colors = clr.get_named_colors_mapping()
		colorsOut = []
		while len(colorsOut) != numberColors:
			rColor = random.choice(list(colors))
			if rColor not in colorsOut:
				colorsOut.append(rColor)
		return colorsOut
			
			
	def getHAVANAPolyANotationEnsemblInRange(self, chromosome, start, stop, showH):
		#open polya signals, returns any signals with a start between the overall start/stop being graphed
		if showH:
			names = ["chr", "group", "type", "start", "end", "misc", "strand", "misc2", "other"]
			chroName = "chr" + chromosome
			annotations =pd.read_csv("gencode.v33.polyAs.gtf",delimiter='\t', names = names, comment = '#' ) #ignores comment lines at beginning of the GTF file'
			return annotations[(annotations['chr'] == chroName) & (annotations['start'] >= start) & (annotations['start'] <= stop)]
		else:
			return []
			
	def __init__(self, chromosome, transcript_positions, transcript_names, gene_names, gene_positions,  gene_strands, gene_colors = [], pas_pos=[], pas_types = [], startOverride = 0, stopOverride = 0, marker_heights=[], marker_size=100, marker_weight=1.5, intron_color="gray", intron_weight=1, intron_style='-', bar_color='gray', bg_color="white", diagramTitle = "", dropDownPASMarkers = True, sequence = "", forwardValues = [], reverseValues = [], forwardPeaks = [], reversePeaks = [], showHAVANAAnnotations = False):
		#chromosome sequence 
		self.chromosome = chromosome
		self.sequence = sequence #used to show location of N's on axis line
		self.forwardCleavagePredictions = forwardValues 
		self.reverseCleavagePredictions = reverseValues
		self.forwardPeaks = forwardPeaks
		self.reversePeaks = reversePeaks
		if forwardValues != [] and reverseValues != []:
			self.yValStart = 3
		else:
			self.yValStart = 1	
		###transcripts
		self.numberGenes = len(gene_names) #number of splice variants on Ensembl for the gene 
		self.transcriptPositions = transcript_positions
		self.transcriptNames = transcript_names
		self.geneNames = gene_names
		self.genePositions = gene_positions
		self.geneStrands = gene_strands
		if gene_colors:
			self.geneColors = gene_colors
		else:
			self.geneColors = self._generateRandomColors(self.numberGenes)
			print (self.geneColors)
		if startOverride == stopOverride == 0:
			self.smallestIndex, self.largestIndex = self.extractTotalSpanFromGenePositions(self.genePositions)
		elif startOverride != 0 and stopOverride != 0:
			self.smallestIndex = startOverride
			self.largestIndex = stopOverride
		elif startOverride != 0 and stopOverride == 0:
			self.smallestIndex = startOverride
			_, self.largestIndex = self.extractTotalSpanFromGenePositions(self.genePositions)
		elif startOverride == 0 and stopOverride != 0:
			self.smallestIndex, _ = self.extractTotalSpanFromGenePositions(self.genePositions)
			self.largestIndex = stopOverride
		#print ("smallest: ", self.smallestIndex, " ", "largest: ", self.largestIndex)
		self.totalSpan =  self.largestIndex - self.smallestIndex
		
		###PAS markers
		self.pasPositions = pas_pos #marker positions to display [[+],[-]]
		self.markerHeights = marker_heights 
		self.pasTypes = pas_types #marker color depends on PAS cluster type
		self.markerSize = marker_size
		self.MarkerWeight = marker_weight
		self.colorKey = {'TE': "red", "EX": "salmon", "IN": "deepskyblue", "DS": "forestgreen", "AE": "darkred", "AI": "midnightblue", "AU": "limegreen", "IG": "gold"} #colors for the 8 PAS cluster types 
		#['pseudo_polyA', 'polyA_signal', 'polyA_site']
		self.colorKeyHAVANA = {"pseudo_polyA": "c", "polyA_signal": "m", "polyA_site": "y"}
		self.showHAVANA = showHAVANAAnnotations
		self.havanaPolyA = self.getHAVANAPolyANotationEnsemblInRange(self.chromosome, self.smallestIndex, self.largestIndex, self.showHAVANA) #get pandas rows of the dataframe with the polyA signals for both strands 
		###graph stuff
		self.dropDownPASMarkers = dropDownPASMarkers
		self.diagramTitle = diagramTitle
		self.intronColor = intron_color
		self.intronWeight = intron_weight
		self.intronStyle = intron_style
		self.barColor= bar_color
		self.bgColor = bg_color
		self.minExonLen = self.totalSpan*0.005 #set a mimimum exon length to show when graphing
		self.ylims = {'exon_width': 1} #max and min for graph exon bars
		self.figure, self.canvas = plt.subplots(figsize=(15,1.5)) 
		self.canvas.set_facecolor(self.bgColor) 
		self._draw() #draw the graph 

	
	
	
	#setting the other y mins and maxes for the intron line points and the overhead bar (exon min/max heights are set above)
	def _set_limits(self):
		self.ylims['intron_max'] = self.ylims['exon_width']*0.75
		self.ylims['intron_min'] = self.ylims['exon_width']*0.5
		
	#if exons exist with a lenght lower than the min exon length to draw, expand the covered interval to be visible
	def _transform_spans(self):
		for index2, geneTranscripts in enumerate(self.transcriptPositions):
			dummy = []
			for index, exons in enumerate(geneTranscripts): #transform each set of exons seperately
				span_lens = [x[1]-x[0] for x in exons] #find all exon lengths in the given sequence
				max_len = float(max(span_lens)) #find the max length exon
				transformed_intervals = []
				if max_len < self.minExonLen: #if the smallest exon is lower than the minimum exon lenght expand the length it covers to be visible
					span_ratios = [x/max_len for x in span_lens]
					expansion_factor = self.totalSpan*1e-11
					for i in range(1,10):
						ef = (2**i)*expansion_factor
						if max_len+ef > self.minExonLen:
							expansion_factor = ef
							break
					for i,j in zip(exons, span_ratios):
						mid = (i[0] + i[1])/2
						f = (expansion_factor*j)/2
						if mid+f - mid-f > self.minExonLen:
							transformed_intervals.append([mid-f, mid+f])
						else:
							transformed_intervals.append([mid-(self.minExonLen/2), mid+(self.minExonLen/2)])
				else:
					for i in range(self.numExons[index]):
						if span_lens[i] < self.minExonLen: #expand if hte exon is lower than the mimimum 
							mid = (exons[i][0] + exons[i][0])/2 
							transformed_intervals.append([mid-(self.minExonLen/2), mid+(self.minExonLen/2)])
						else: #do nothing if it is larger than the minimal exon length
							transformed_intervals.append(exons[i])
				dummy.append(transformed_intervals)
				#self.transcripts[index] = transformed_intervals #overwrite exon lengths for the graph
			self.transcriptPositions[index2] = dummy
	
	def _draw_exon(self, span, yValStart, geneColor, strand): #draw a rectangle between the two limits in the the span 
		minY = yValStart
		maxY =  self.ylims['exon_width'] +yValStart
		#print ("placing exon at: ", span)
		#print ("y vals: ", minY, " ", maxY) 
		if strand == "+":
			self.canvas.fill_between(span, -minY, -maxY, edgecolor=self.bgColor, facecolor=geneColor)
		else:
			self.canvas.fill_between(span, minY, maxY, edgecolor=self.bgColor, facecolor=geneColor)
			
		return True  #if successful in drawing the exon, return true

	def _draw_intron(self, span, yValStart, strand): #span is between last exon end point x and next exon's starting point 
		mid = (span[0]+span[1])/2.0 #find midpoint between the two exons 
		#make two line segments for the intron 
		minY = yValStart + self.ylims['intron_min']
		maxY = yValStart + self.ylims['intron_max']
		#print ("drawing intron between: ", span)
		#print ("y vals: ", minY, " ", maxY)
		if strand == "+":
			self.canvas.plot([span[0], mid], [-minY, -maxY], c=self.intronColor, lw=self.intronWeight, ls=self.intronStyle)
			self.canvas.plot([mid, span[1]], [-maxY,-minY], c=self.intronColor, lw=self.intronWeight, ls=self.intronStyle)
		else:
			self.canvas.plot([span[0], mid], [minY, maxY], c=self.intronColor, lw=self.intronWeight, ls=self.intronStyle)
			self.canvas.plot([mid, span[1]], [maxY, minY], c=self.intronColor, lw=self.intronWeight, ls=self.intronStyle)
		return True #if successful in drawing the intron, return true		
	
	
	
	def _draw_markers(self, yMax):
		if self.pasPositions: #if they are given
			#plus strand
			for i in range(0,2):
				for p,t in zip(self.pasPositions[i], self.pasTypes[i]): #for each marker make a bar and then place a scatter plot dot on the bar
					h = 0.5 * (-1) ** (i + 1)
					if self.dropDownPASMarkers:
						h = 1.1 * yMax * (-1) ** (i + 1)
					markerColor = self.colorKey[t] #retrieve color for the pas type
					if type(p) == type(1):
						self.canvas.plot((p, p), (0, h), linestyle='-', color=markerColor, linewidth=self.MarkerWeight)
					elif type(p) == type([]) and len(p) == 2:
						self.canvas.fill_between(p, 0, h, edgecolor= markerColor, facecolor= markerColor)
				
			#if there are PAS markers add a legend for the types
			types = [Patch(facecolor='red', edgecolor='red', label='TE'), 
					Patch(facecolor='salmon', edgecolor='salmon', label='EX'), 
					Patch(facecolor='deepskyblue', edgecolor='deepskyblue', label='IN'), 
					Patch(facecolor='forestgreen', edgecolor='forestgreen', label='DS'), 
					Patch(facecolor='darkred', edgecolor='darkred', label='AE'), 
					Patch(facecolor='midnightblue', edgecolor='midnightblue', label='AI'), 
					Patch(facecolor='limegreen', edgecolor='limegreen', label='AU'), 
					Patch(facecolor='gold', edgecolor='gold', label='IG')]
			l2 = plt.legend(handles=types, loc = 'center right', bbox_to_anchor = (0,0.5))
			self.canvas.add_artist(l2)
			#ax.legend(loc='center left', bbox_to_anchor=(1, 0.5)))
			
	def _draw_HAVANA_markers(self, yMax):
		if self.showHAVANA: #if they are given
			print (self.havanaPolyA)
			for index, row in self.havanaPolyA.iterrows():
				#for each pas in the range
				print (row)
				print (row['type'])
				markerColor = self.colorKeyHAVANA[row['type']]
				if row['strand'] == "+": #lower half of graph
					h = 0.5 
					if self.dropDownPASMarkers:
						h = 1.1 * yMax * (-1)
					if row['start'] == row['end']:
						#point 
						self.canvas.plot((p, p), (0, h), linestyle='-', color=markerColor, linewidth=self.MarkerWeight)
					else:
						p = [row['start'], row['end']]
						self.canvas.fill_between(p, 0, h, edgecolor= markerColor, facecolor= markerColor)
							
				else: #minus strand upper half of graph
					h = 0.5 
					if self.dropDownPASMarkers:
						h = 1.1 * yMax
					if row['start'] == row['end']:
						#point 
						self.canvas.plot((p, p), (0, h), linestyle='-', color=markerColor, linewidth=self.MarkerWeight)
					else:
						p = [row['start'], row['end']]
						self.canvas.fill_between(p, 0, h, edgecolor= markerColor, facecolor= markerColor)
			#self.colorKeyHAVANA = {'pseudo_polyA': "c", 'polyA_signal': "m", "polyA_site": "y"}
			types = [Patch(facecolor='c', edgecolor='c', label='HAVANA psuedo_polyA'), 
					Patch(facecolor='m', edgecolor='m', label='HAVANA polyA_signal'), 
					Patch(facecolor='y', edgecolor='y', label='HAVANA polyA_site')]
			l2 = plt.legend(handles=types, loc = 'center left', bbox_to_anchor = (1,0.5))
			self.canvas.add_artist(l2)
			
	
	def _draw_Peak_Markers(self, yMax):
		#print ("forward peaks: ", self.forwardPeaks)
		#print ("reverse peaks: ", self.reversePeaks)
		for p in self.forwardPeaks: #for each marker make a bar and then place a scatter plot dot on the bar
			#print (p)
			p = p + self.smallestIndex #go from 0 indexing to 1 indexing
			h = -0.25
			if self.dropDownPASMarkers:
				h = 0.55 * yMax * (-1)
			self.canvas.plot((p, p), (0, h), linestyle='-', color='black', linewidth=self.MarkerWeight)
		for p in self.reversePeaks: #for each marker make a bar and then place a scatter plot dot on the bar
			#print (p)
			p = p + self.smallestIndex #go from 0 indexing to 1 indexing 
			h = 0.25 
			if self.dropDownPASMarkers:
				h = 0.55 * yMax 
			self.canvas.plot((p, p), (0, h), linestyle='-', color='black', linewidth=self.MarkerWeight)
			
	
			
	#place markers between smallest/largest
	def _clean_axes(self):
		self.canvas.set_yticks([], []) #hides all y-axis ticks
		self.canvas.get_xaxis().tick_top() 
		self.canvas.tick_params(axis='x', direction='out') 
		self.canvas.set_xticks([]) #remove x ticks 
		for o in ["top", "bottom", "left", "right"]: #removing the box around the graph 
			self.canvas.spines[o].set_visible(False)
		min_pos = int(self.smallestIndex - self.totalSpan * 0.1) 
		if min_pos < 0:
			min_pos = 0
		max_pos = int(self.largestIndex + self.totalSpan * 0.1) 
		stepSize = int((max_pos-min_pos)/20) #force stepsize to be an int for placing range markers 
		minortick_pos = [x for x in range(min_pos, max_pos, stepSize)][1:] #divide up the total range by 20 and place tick marks 
		for i in minortick_pos:
			self.canvas.axvline(i, alpha=0.2, c='black', ls='--') #add vertical line across the axes 
		#self.canvas.text(minortick_pos[0], -0.5, minortick_pos[0], fontsize=8, ha='center') #place left value marker
		#self.canvas.text(minortick_pos[-1], -0.5, minortick_pos[-1], fontsize=8, ha='center') #place right value marker
		self.canvas.set_xlim(minortick_pos[0]-(minortick_pos[1]-minortick_pos[0]), minortick_pos[-1]+(minortick_pos[-1]-minortick_pos[-2])) #limit graph x-value range 
        
	def _drawGeneBoundaries(self, largestY):
		#draw lines extending up/down to show gene overall boundaries 
		lengthLine = largestY * 1.1
		patches = []
		for g in range(0, self.numberGenes):
			if self.geneStrands[g] == "+":
				self.canvas.plot([self.genePositions[g][0], self.genePositions[g][0]], [0, -1 * lengthLine], c= self.geneColors[g], lw=self.intronWeight, ls=":")
				self.canvas.plot([self.genePositions[g][1], self.genePositions[g][1]], [0, -1 * lengthLine], c= self.geneColors[g], lw=self.intronWeight, ls=":")
			else:
				self.canvas.plot([self.genePositions[g][0], self.genePositions[g][0]], [0, lengthLine], c= self.geneColors[g], lw=self.intronWeight, ls=":")
				self.canvas.plot([self.genePositions[g][1], self.genePositions[g][1]], [0, lengthLine], c= self.geneColors[g], lw=self.intronWeight, ls=":")
			#make legend of genes
			patches.append(Patch(facecolor = self.geneColors[g], edgecolor = self.geneColors[g], label = self.geneNames[g]))	
		self.canvas.legend(handles=patches, bbox_to_anchor=(0.5, -0.05), loc = 'upper center', ncol = self.numberGenes)


	def _dotNPositionsOnAxis(self):
		nPositions = [i for i, ltr in enumerate(self.sequence) if ltr == "N"]
		for position in nPositions:
			self.canvas.scatter(0, position + 1, marker = 'o', c = 'red')



	def _plotPredictions(self):
		xAxis = list(range(self.smallestIndex, self.largestIndex + 1))
		#print ("xaxis: ", len(xAxis))
		#print ("forward: ", self.forwardCleavagePredictions.size)
		#print ("reverse: ", self.reverseCleavagePredictions.size)
		if self.forwardCleavagePredictions != []:
			fixedForward = self.yValStart * self.forwardCleavagePredictions * -1 
			self.canvas.plot(xAxis, fixedForward)
		if self.reverseCleavagePredictions != []:
			fixedReverse = self.yValStart * self.reverseCleavagePredictions
			self.canvas.plot(xAxis, fixedReverse)
		
	def _draw(self):
		self._set_limits() #setting y values for the intron lines and the overhead bar in the graph
		#if self.thickenExons:
		#	self._transform_spans() #expanding very small exons to show on the visual for the region
		maxY = float('-inf')
		for k in range(0,self.numberGenes):
			transcriptsNow = self.transcriptPositions[k]
			yValStart = self.yValStart
			for j in range(0,len(transcriptsNow)):
				currentExonColor = self.geneColors[k]
				currentExons = transcriptsNow[j]
				#print (currentExons)
				currentExons = sorted(currentExons)
				#print (currentExons)
				for i in range(0,len(currentExons)):
					if i > 0: #if its not the first exon, draw the intron first
						self._draw_intron([currentExons[i-1][1], currentExons[i][0]], yValStart, self.geneStrands[k]) #between last exon end point and current exon start
					self._draw_exon(currentExons[i], yValStart, currentExonColor, self.geneStrands[k]) #draw first exon or following exons
					if i == len(currentExons) - 1 and self.transcriptNames: #to the right of the last exon place the transcript name if on exists
						tagPos = 0
						if self.geneStrands[k] == "+":
							tagPos = -1 * yValStart
						else:
							tagPos = yValStart
						self.canvas.text( self.largestIndex, tagPos, self.transcriptNames[k][j], fontsize = 8, va = "top", ha = "left")	
				yValStart += self.ylims['exon_width']
			if yValStart >= maxY:
				maxY = yValStart		
		
		self.canvas.plot([self.smallestIndex, self.largestIndex], [0,0], c='black', lw=self.intronWeight, ls=self.intronStyle)
		if self.sequence != "":
			self._dotNPositionsOnAxis()
		if self.yValStart != 1:
			self._plotPredictions()
		if self.diagramTitle:
			plt.title(self.diagramTitle)
		self._draw_markers(maxY) #draw markers with the default color for the diagram instance
		self._draw_HAVANA_markers(maxY)
		self._draw_Peak_Markers(maxY)
		self._drawGeneBoundaries(maxY)
		self._clean_axes()

	def show(self):
		plt.show()
		
		
		
		
class UCSCAnnotationDiagram(object):
	
	def _dataframeToRegions(self, df):
		regions = []
		for index, row in df.iterrows():
			regions.append([row.start, row.end])
		regions = sorted(regions)
		return regions
		
			
	def __init__(self, whole_gene, upstream, downstream, exons, introns, coding_exons, three_prime_utr, five_prime_utr, pas_pos=[], pas_types = [], marker_heights=[], marker_size=100, marker_weight=1.5, intron_color="gray", intron_weight=1, intron_style='-', bar_color='gray', bg_color="white", diagramTitle = "", dropDownPASMarkers = True):
		###transcripts
		self.wholeGene = self._dataframeToRegions(whole_gene)
		self.strand = whole_gene.iloc[0]['strand']
		self.upstream = self._dataframeToRegions(upstream)
		self.downstream = self._dataframeToRegions(downstream)
		self.exons = self._dataframeToRegions(exons)
		self.introns = self._dataframeToRegions(introns)
		self.codingExons = self._dataframeToRegions(coding_exons)
		self.threePrimeUTR = self._dataframeToRegions(three_prime_utr)
		self.fivePrimeUTR = self._dataframeToRegions(five_prime_utr)
		self.allRegions = [self.wholeGene, self.upstream, self.downstream, self.exons, self.introns, self.codingExons, self.threePrimeUTR, self.fivePrimeUTR]
		self.smallestIndex = self.upstream[0][0]
		self.largestIndex = self.downstream[0][1]
		self.totalSpan =  self.largestIndex - self.smallestIndex
		self.colors = ['grey', 'green', 'purple', 'red', 'blue', 'salmon', 'orange', 'yellow'] #whole gene, upstream, downstream, exons, introns, cds exons, 3', 5' 
		self.names = ['gene', 'upstream 1000 nts', 'downstream 1000 nts', 'exons', 'introns', 'CDS exons', "3'UTR", "5'UTR"] #names for graphs
		
		###PAS markers
		self.pasPositions = pas_pos #marker positions to display [[+],[-]]
		self.markerHeights = marker_heights 
		self.pasTypes = pas_types #marker color depends on PAS cluster type
		self.markerSize = marker_size
		self.MarkerWeight = marker_weight
		self.colorKey = {'TE': "red", "EX": "salmon", "IN": "deepskyblue", "DS": "forestgreen", "AE": "darkred", "AI": "midnightblue", "AU": "limegreen", "IG": "gold"} #colors for the 8 PAS cluster types 
		###graph stuff
		self.dropDownPASMarkers = dropDownPASMarkers
		self.diagramTitle = diagramTitle
		self.barColor= bar_color
		self.bgColor = bg_color
		self.minExonLen = self.totalSpan*0.005 #set a mimimum exon length to show when graphing
		self.ylims = {'exon_width': 1} #max and min for graph exon bars
		self.figure, self.canvas = plt.subplots(figsize=(15,1.5)) 
		self.canvas.set_facecolor(self.bgColor) 
		self._draw() #draw the graph 

	
	def _draw_exon(self, span, yValStart, geneColor, strand): #draw a rectangle between the two limits in the the span 
		minY = yValStart
		maxY =  self.ylims['exon_width'] +yValStart
		if strand == "+":
			self.canvas.fill_between(span, -minY, -maxY, edgecolor=self.bgColor, facecolor=geneColor)
		else:
			self.canvas.fill_between(span, minY, maxY, edgecolor=self.bgColor, facecolor=geneColor)
			
		return True  #if successful in drawing the exon, return true
		
	def _draw_markers(self, yMax):
		if self.pasPositions: #if they are given
			#plus strand
			for i in range(0,2):
				for p,t in zip(self.pasPositions[i], self.pasTypes[i]): #for each marker make a bar and then place a scatter plot dot on the bar
					h = 0.5 * (-1) ** (i + 1)
					if self.dropDownPASMarkers:
						h = 1.1 * yMax * (-1) ** (i + 1)
					markerColor = self.colorKey[t] #retrieve color for the pas type
					if type(p) == type(1):
						self.canvas.plot((p, p), (0, h), linestyle='-', color=markerColor, linewidth=self.MarkerWeight)
					elif type(p) == type([]) and len(p) == 2:
						self.canvas.fill_between(p, 0, h, edgecolor= markerColor, facecolor= markerColor)
				
			#if there are PAS markers add a legend for the types
			types = [Patch(facecolor='red', edgecolor='red', label='TE'), 
					Patch(facecolor='salmon', edgecolor='salmon', label='EX'), 
					Patch(facecolor='deepskyblue', edgecolor='deepskyblue', label='IN'), 
					Patch(facecolor='forestgreen', edgecolor='forestgreen', label='DS'), 
					Patch(facecolor='darkred', edgecolor='darkred', label='AE'), 
					Patch(facecolor='midnightblue', edgecolor='midnightblue', label='AI'), 
					Patch(facecolor='limegreen', edgecolor='limegreen', label='AU'), 
					Patch(facecolor='gold', edgecolor='gold', label='IG')]
			l2 = plt.legend(handles=types, loc = 'center right', bbox_to_anchor = (0,0.5))
			self.canvas.add_artist(l2)
			#ax.legend(loc='center left', bbox_to_anchor=(1, 0.5)))
			
	#place markers between smallest/largest
	def _clean_axes(self):
		self.canvas.set_yticks([], []) #hides all y-axis ticks
		self.canvas.get_xaxis().tick_top() 
		self.canvas.tick_params(axis='x', direction='out') 
		self.canvas.set_xticks([]) #remove x ticks 
		for o in ["top", "bottom", "left", "right"]: #removing the box around the graph 
			self.canvas.spines[o].set_visible(False)
		min_pos = int(self.smallestIndex - self.totalSpan * 0.1) 
		if min_pos < 0:
			min_pos = 0
		max_pos = int(self.largestIndex + self.totalSpan * 0.1) 
		stepSize = int((max_pos-min_pos)/20) #force stepsize to be an int for placing range markers 
		minortick_pos = [x for x in range(min_pos, max_pos, stepSize)][1:] #divide up the total range by 20 and place tick marks 
		for i in minortick_pos:
			self.canvas.axvline(i, alpha=0.2, c='black', ls='--') #add vertical line across the axes 
		self.canvas.set_xlim(minortick_pos[0]-(minortick_pos[1]-minortick_pos[0]), minortick_pos[-1]+(minortick_pos[-1]-minortick_pos[-2])) #limit graph x-value range 
		
	def _draw(self):
		#self._set_limits() #setting y values for the intron lines and the overhead bar in the graph
		#if self.thickenExons:
		#	self._transform_spans() #expanding very small exons to show on the visual for the region
		maxY = float('-inf')
		yValStart = 1
		for j in range(0,len(self.allRegions)):
			currentExonColor = self.colors[j]
			currentExons = self.allRegions[j]
			currentExons = sorted(currentExons)
			for i in range(0,len(currentExons)):
				self._draw_exon(currentExons[i], yValStart, currentExonColor, self.strand) #draw first exon or following exons
				if i == len(currentExons) - 1: #to the right of the last exon place the transcript name if on exists
					tagPos = 0
					if self.strand == "+":
						tagPos = -1 * yValStart
					else:
						tagPos = yValStart
					self.canvas.text( self.largestIndex, tagPos, self.names[j], fontsize = 8, va = "top", ha = "left")	
			yValStart += self.ylims['exon_width']
			if yValStart >= maxY:
				maxY = yValStart		
		self.canvas.plot([self.smallestIndex, self.largestIndex], [0,0], c='black', lw=1, ls='-')
		if self.diagramTitle:
			plt.title(self.diagramTitle)
		self._draw_markers(maxY) #draw markers with the default color for the diagram instance
		self._clean_axes()

	def show(self):
		plt.show()
		
if  __name__ == "__main__":
	#exons positions
	exon_pos = [97543299,97544702], [97547885,97548026], [97564044,97564188], [97658624,97658804], [97700407,97700550], [97770814,97770934], [97771732,97771853], [97839116,97839200], [97847948,97848017], [97915614,97915779], [97981281,97981497], [98015115,98015300], [98039315,98039526], [98058773,98058943], [98060614,98060722], [98144650,98144738], [98157272,98157354], [98164906,98165103], [98187065,98187227], [98205947,98206035], [98293669,98293752], [98348819,98348930], [98386439,98386615]
	exon_pos2 = [97543299,97544702], [97547885,97548026], [97564044,97564188], [97658624,97658804], [97700407,97700550], [97770814,97770934], [97771732,97771853], [97839116,97839200], [97847948,97848017], [97915614,97915779], [97981281,97981497], [98015115,98015300], [98039315,98039526], [98058773,98058943], [98060614,98060722], [98144650,98144738], [98157272,98157354], [98164906,98165103]
	exon_pos3 = [97543299,97544702], [97547885,97548026], [97564044,97564188], [97658624,97658804], [97700407,97700550], [97770814,97770934], [97771732,97771853], [97839116,97839200], [97847948,97848017], [97915614,97915779], [97981281,97981497]
	exon_pos4 = [97543299,97544702], [97547885,97548026], [97564044,97564188], [97658624,97658804], [97700407,97700550], [97770814,97770934], [97771732,97771853], [97839116,97839200], [97847948,97848017], [97915614,97915779], [97981281,97981497]
	exon_pos5 = [97543299,97544702], [97547885,97548026], [97564044,97564188], [97658624,97658804], [97700407,97700550], [97770814,97770934], [97771732,97771853], [97839116,97839200], [97847948,97848017], [97915614,97915779], [97981281,97981497]
	#marker positions
	marker_pos = [97947885, 98247485, [97548026,97548506 ]]
	pas_types = ['IN', 'EX', 'TE']
	exonColors = ["red", "blue", "green", "black", "grey"]
	transcript_names = ["top", "middle", "bottom", "b2", "b3"]
	multiGenes = SpliceVariantPASDiagram([exon_pos, exon_pos2, exon_pos3,exon_pos4,exon_pos5], transcript_names = transcript_names, pas_pos = marker_pos, pas_types = pas_types, exon_colors = exonColors, diagramTitle = "DEMO")
	multiGenes.show()

