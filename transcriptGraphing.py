import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Circle
from operator import itemgetter

'''
class ExonDiagram(object):
	def __init__(self, exons, diagramTitle):
		self.exons = exons
		self.diagramTitle = diagramTitle
		self._sortExons()
	
	def _sortExons(self):
		#sorts exons by earliest starting position 
		self.exons.sort(key = itemgetter(0))
		print (self.exons)
		
'''			
		


class SpliceVariantPASDiagram(object):
	
	def __init__(self, transcripts, strands = [], transcript_names = [],  pas_pos=[], pas_types = [], marker_heights=[], marker_size=100, marker_weight=1.5, exon_colors= [], intron_color="gray", intron_weight=1, intron_style='-', bar_color='gray', bg_color="white", diagramTitle = "", thickenExons = True, dropDownPASMarkers = False):
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
			self.canvas.legend(handles=types, loc = 'center right', bbox_to_anchor = (0,0.5))
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
	
	def __init__(self, transcripts, strands = [], transcript_names = [],  pas_pos=[], pas_types = [], marker_heights=[], marker_size=100, marker_weight=1.5, exon_colors= [], intron_color="gray", intron_weight=1, intron_style='-', bar_color='gray', bg_color="white", diagramTitle = ""):
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
		###graph stuff
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
			for p,h,t in zip(self.pasPositions, self.markerHeights, self.pasTypes): #for each marker make a bar and then place a scatter plot dot on the bar
				markerColor = self.colorKey[t] #retrieve color for the pas type
				if type(p) == type(1):
					self.canvas.plot((p, p), (self.ylims['bar_max'], self.ylims['bar_max']+h), linestyle='-', color=markerColor, linewidth=self.MarkerWeight, alpha=0.7)
					self.canvas.scatter(p, self.ylims['bar_max']+h+0.25, s=self.markerSize, marker='o', c=markerColor, edgecolor=markerColor, alpha=1)
				elif type(p) == type([]) and len(p) == 2:
					midpoint = (p[0] + p[1] )/ 2.0
					self.canvas.fill_between(p, self.ylims['bar_max'], self.ylims['bar_max'] + h, edgecolor= markerColor, facecolor= markerColor)
					self.canvas.scatter(midpoint, self.ylims['bar_max']+h+0.25, s=self.markerSize, marker='o', c=markerColor, edgecolor=markerColor, alpha=1)
			#if there are PAS markers add a legend for the types
			#'TE': "red", "EX": "salmon", "IN": "deepskyblue", "DS": "forestgreen", "AE": "darkred", "AI": "midnightblue", "AU": "limegreen", "IG": "gold"} #colors for the 8 PAS cluster types 
			types = [Patch(facecolor='red', edgecolor='red', label='TE'), 
					Patch(facecolor='salmon', edgecolor='salmon', label='EX'), 
					Patch(facecolor='deepskyblue', edgecolor='deepskyblue', label='IN'), 
					Patch(facecolor='forestgreen', edgecolor='forestgreen', label='DS'), 
					Patch(facecolor='darkred', edgecolor='darkred', label='AE'), 
					Patch(facecolor='midnightblue', edgecolor='midnightblue', label='AI'), 
					Patch(facecolor='limegreen', edgecolor='limegreen', label='AU'), 
					Patch(facecolor='gold', edgecolor='gold', label='IG')]
			self.canvas.legend(handles=types, loc = 'center right', bbox_to_anchor = (0,0.5))
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

'''
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
'''
