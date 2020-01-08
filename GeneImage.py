import matplotlib.pyplot as plt

"""
MIT License
Copyright (c) [2016] [Parashar Dhapola]
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

'''
__author__ = "Parashar Dhapola"
__email__ = "parashar.dhapola@gmail.com"
Original author information 
'''

#AML
#1/6/19 
#downloaded from https://gist.github.com/parashardhapola/e095e3ebd068c90e7089b53e46b8e0bc
#Updated to current matplotlib methods
'''
Editing to:
	a) display multiple transcripts for the gene at once
	b) allow 'PAS clusters' to be shown as modified markers -> set each type to different colors to see them easier????
	c) make option to remove axes and show simplified exon view (with reduced intron distances)
'''


class GeneImage(object):
    #exon_intervals is massive tuple w/ lists of [start,end] for each exon in the transcript/gene.  Must be ordered w/ lowest number exon to highest number exon 
    #marker_pos is a list of marker positions 
    #marker_heights
    #marker_colors
    #marker_size
    #marker_weights
    #exon_color
    #intron color
    #intron weight: line weight
    #intron style: line style
    #bar_color: color of overhead bar
    #bg_color: background color
    def __init__(self, exon_intervals, marker_pos=[], marker_heights=[], marker_colors=[],
                 marker_size=100, marker_weight=1.5, exon_color="black", intron_color="grey",
                 intron_weight=2, intron_style='-', bar_color='cornflowerblue', bg_color="white"):
        self.exonIntervals = exon_intervals #variable length tuple of lists of start and end positions
        self.markerPositions = marker_pos #marker positions to display
        self.markerHeights = marker_heights 
        self.markerColors = marker_colors #TODO: change so that each subtype of PAS cluster has its own color
        self.markerSize = marker_size
        self.MarkerWeight = marker_weight
        self.exonColor = exon_color
        self.intronColor = intron_color
        self.intronWeight = intron_weight
        self.intronStyle = intron_style
        self.barColor= bar_color
        self.bgColor = bg_color
        self.markerDefaultColor = 'grey'
        self.numExons = len(self.exonIntervals) #number of exons in the exon tuple
        self.totalSpan = self.exonIntervals[-1][1] - self.exonIntervals[0][0] 
        self.minExonLen = self.totalSpan*0.005 #set a mimimum exon length to show when graphing
        self.ylims = {'exon_max': 2, 'exon_min':1} #max and min for graph exon bars
        self.figure, self.canvas = plt.subplots(figsize=(15,1.5)) #set up transcript graph to make
        self.canvas.set_facecolor(self.bgColor) #AML: changed set_axis_bgcolor to set_facecolor, set background color
        self._draw() #draw the graph 

    #setting the other y mins and maxes for the intron line points and the overhead bar (exon min/max heights are set above)
    def _set_limits(self):
        self.ylims['intron_max'] = self.ylims['exon_max']*0.9
        self.ylims['intron_min'] = (self.ylims['exon_max'] + self.ylims['exon_min'])/2.0
        self.ylims['bar_min'] = self.ylims['exon_max']+0.2 #setting the height of the overhead bar
        self.ylims['bar_max'] = self.ylims['bar_min']+(self.ylims['exon_max']-self.ylims['exon_min'])/5.0 #setting the height of the overhead bar
        
    #if exons exist with a lenght lower than the min exon length to draw, expand the covered interval to be visible
    def _transform_spans(self):
        span_lens = [x[1]-x[0] for x in self.exonIntervals] #find all exon lengths in the given sequence
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
            for i,j in zip(self.exonIntervals, span_ratios):
                mid = (i[0] + i[1])/2
                f = (expansion_factor*j)/2
                if mid+f - mid-f > self.minExonLen:
                    transformed_intervals.append([mid-f, mid+f])
                else:
                    transformed_intervals.append([mid-(self.minExonLen/2), mid+(self.minExonLen/2)])
        else:
            for i in range(self.numExons):
                if span_lens[i] < self.minExonLen: #expand if hte exon is lower than the mimimum 
                    mid = (self.exonIntervals[i][0] + self.exonIntervals[i][0])/2 
                    transformed_intervals.append([mid-(self.minExonLen/2), mid+(self.minExonLen/2)])
                else: #do nothing if it is larger than the minimal exon length
                    transformed_intervals.append(self.exonIntervals[i])
        self.exonIntervals = transformed_intervals[:] #overwrite exon lengths for the graph
        
    def _draw_exon(self, span): #draw a rectangle between the two limits in the the span 
        self.canvas.fill_between(span, self.ylims['exon_min'], self.ylims['exon_max'],
                                 edgecolor=self.bgColor, facecolor=self.exonColor)
        return True  #if successful in drawing the exon, return true
        
    def _draw_intron(self, span): #span is between last exon end point x and next exon's starting point 
        mid = (span[0]+span[1])/2.0 #find midpoint between the two exons 
        #make two line segments for the intron 
        self.canvas.plot([span[0], mid], [self.ylims['intron_min'], self.ylims['intron_max']],
                         c=self.intronColor, lw=self.intronWeight, ls=self.intronStyle)
        self.canvas.plot([mid, span[1]], [self.ylims['intron_max'], self.ylims['intron_min']],
                         c=self.intronColor, lw=self.intronWeight, ls=self.intronStyle)
        return True #if successful in drawing the intron, return true
    
    def _draw_markers(self):
        if self.markerHeights == []:
            self.markerHeights = [self.ylims['exon_max']-self.ylims['exon_min'] for x in self.markerPositions] #make the markers the same height as the exon blocks
        if self.markerColors == []:
            self.markerColors = [self.markerDefaultColor for x in self.markerPositions]           
        for p,h,c in zip(self.markerPositions, self.markerHeights, self.markerColors): #for each marker make a bar and then place a scatter plot dot on the bar
            self.canvas.plot((p, p), (self.ylims['bar_max'], self.ylims['bar_max']+h),
                             linestyle='-', color='black', linewidth=self.MarkerWeight, alpha=0.7)
            self.canvas.scatter(p, self.ylims['bar_max']+h+0.25, s=self.markerSize, marker='o', c=c,
                                edgecolor=c, alpha=1)
        
    
    def _clean_axes(self):
        #self.canvas.set_ylim((self.ylims['exon_min'], self.ylims['bar_max']))
        self.canvas.set_yticks([], []) #hides all y-axis ticks
        self.canvas.get_xaxis().tick_top() #move ticks and tick labels to the top of the axis
        self.canvas.tick_params(axis='x', direction='out') #move ticks outside of axis 
        self.canvas.set_xticks([]) #remove x ticks 
        for o in ["top", "bottom", "left", "right"]: #removing the box around the graph 
            self.canvas.spines[o].set_visible(False)
        min_pos = int(self.exonIntervals[0][0] - self.totalSpan * 0.1) #find bound for the left
        if min_pos < 0:
            min_pos = 0
        max_pos = int(self.exonIntervals[-1][1] + self.totalSpan * 0.1) #find bound for the right 
        #print (min_pos, max_pos, (max_pos-min_pos)/20)
        stepSize = int((max_pos-min_pos)/20) #force stepsize to be an int for placing range markers 
        minortick_pos = [x for x in range(min_pos, max_pos, stepSize)][1:] #divide up the total range by 20 and place tick marks 
        for i in minortick_pos:
            self.canvas.axvline(i, alpha=0.2, c='black', ls='--') #add vertical line across the axes 
        self.canvas.text(minortick_pos[0], self.ylims['exon_min']-0.5,
                         minortick_pos[0], fontsize=8, ha='center') #place left value marker
        self.canvas.text(minortick_pos[-1], self.ylims['exon_min']-0.5,
                         minortick_pos[-1], fontsize=8, ha='center') #place right value marker
        self.canvas.set_xlim(minortick_pos[0]-(minortick_pos[1]-minortick_pos[0]),
                             minortick_pos[-1]+(minortick_pos[-1]-minortick_pos[-2])) #limit graph x-value range 
        
    def _draw(self):
        self._set_limits() #setting y values for the intron lines and the overhead bar in the graph
        self._transform_spans() #expanding very small exons to show on the visual for the region
        for i in range(self.numExons):
            if i > 0: #if its not the first exon, draw the intron first
                self._draw_intron([self.exonIntervals[i-1][1], self.exonIntervals[i][0]]) #between last exon end point and current exon start
            self._draw_exon(self.exonIntervals[i]) #draw first exon or following exons
        self.canvas.fill_between([self.exonIntervals[0][0], self.exonIntervals[-1][1]],
                                  self.ylims['bar_min'], self.ylims['bar_max'],
                                  edgecolor=self.bgColor, facecolor=self.barColor)
        self._draw_markers() #draw markers with the default color for the diagram instance
        self._clean_axes()
    
    def show(self):
        plt.show()
        
        
class SpliceVariantPASDiagram(object):
    colorKey = {'TE': "red", "EX": "salmon", "IN": "deepskyblue", "DS": "forestgreen", "AE": "darkred", "AI": "midnightblue", "AU": "limegreen", "IG": "gold"} #colors for the 8 PAS cluster types 
    
    def __init__(self, transcripts, strands = [], strand_names = [],  pas_pos=[], pas_types = [], marker_heights=[],
                 marker_size=100, marker_weight=1.5, exon_colors= [], intron_color="gray",
                 intron_weight=2, intron_style='-', bar_color='gray', bg_color="white"):
        ###transcripts
        self.numberTranscripts = len(transcripts) #number of splice variants on Ensembl for the gene 
        self.transcripts = transcripts 
        self.strands = strands 
        self.strand_names = strand_names
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
        ###graph stuff
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
        span_lens = [x[1]-x[0] for x in self.exonIntervals] #find all exon lengths in the given sequence
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
            for i,j in zip(self.exonIntervals, span_ratios):
                mid = (i[0] + i[1])/2
                f = (expansion_factor*j)/2
                if mid+f - mid-f > self.minExonLen:
                    transformed_intervals.append([mid-f, mid+f])
                else:
                    transformed_intervals.append([mid-(self.minExonLen/2), mid+(self.minExonLen/2)])
        else:
            for i in range(self.numExons):
                if span_lens[i] < self.minExonLen: #expand if hte exon is lower than the mimimum 
                    mid = (self.exonIntervals[i][0] + self.exonIntervals[i][0])/2 
                    transformed_intervals.append([mid-(self.minExonLen/2), mid+(self.minExonLen/2)])
                else: #do nothing if it is larger than the minimal exon length
                    transformed_intervals.append(self.exonIntervals[i])
        self.exonIntervals = transformed_intervals[:] #overwrite exon lengths for the graph
        
    def _draw_exon(self, span, offset): #draw a rectangle between the two limits in the the span 
        self.canvas.fill_between(span, self.ylims['exon_min'] - offset, self.ylims['exon_max'] - offset,
                                 edgecolor=self.bgColor, facecolor=self.exonColor)
        return True  #if successful in drawing the exon, return true
        
    def _draw_intron(self, span, offset): #span is between last exon end point x and next exon's starting point 
        mid = (span[0]+span[1])/2.0 #find midpoint between the two exons 
        #make two line segments for the intron 
        self.canvas.plot([span[0], mid], [self.ylims['intron_min'] - offset, self.ylims['intron_max'] - offset],
                         c=self.intronColor, lw=self.intronWeight, ls=self.intronStyle)
        self.canvas.plot([mid, span[1]], [self.ylims['intron_max'] - offset, self.ylims['intron_min'] - offset],
                         c=self.intronColor, lw=self.intronWeight, ls=self.intronStyle)
        return True #if successful in drawing the intron, return true
    
    def _draw_markers(self):
        if self.markerHeights == []:
            self.markerHeights = [self.ylims['exon_max']-self.ylims['exon_min'] for x in self.markerPositions] #make the markers the same height as the exon blocks
        if self.markerColors == []:
            self.markerColors = [self.markerDefaultColor for x in self.markerPositions]           
        for p,h,t in zip(self.markerPositions, self.markerHeights, self.pasTypes): #for each marker make a bar and then place a scatter plot dot on the bar
            c = colorKey[t] #retrieve color for the pas type
            self.canvas.plot((p, p), (self.ylims['bar_max'], self.ylims['bar_max']+h),
                             linestyle='-', color='black', linewidth=self.MarkerWeight, alpha=0.7)
            self.canvas.scatter(p, self.ylims['bar_max']+h+0.25, s=self.markerSize, marker='o', c=c,
                                edgecolor=c, alpha=1)
        
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
        self.canvas.text(minortick_pos[0], self.ylims['exon_min']-0.5,
                         minortick_pos[0], fontsize=8, ha='center') #place left value marker
        self.canvas.text(minortick_pos[-1], self.ylims['exon_min']-0.5,
                         minortick_pos[-1], fontsize=8, ha='center') #place right value marker
        self.canvas.set_xlim(minortick_pos[0]-(minortick_pos[1]-minortick_pos[0]),
                             minortick_pos[-1]+(minortick_pos[-1]-minortick_pos[-2])) #limit graph x-value range 
        
    def _draw(self):
        self._set_limits() #setting y values for the intron lines and the overhead bar in the graph
        self._transform_spans() #expanding very small exons to show on the visual for the region
        for j in range(0,len(self.transcripts)):
        	offset = (j + 1) * 3
		for i in range(self.numExons[j]):
		    if i > 0: #if its not the first exon, draw the intron first
		        self._draw_intron([self.transcripts[j][i-1][1], self.transcripts[j][i][0]], offset) #between last exon end point and current exon start
		    self._draw_exon(self.transcripts[j][i], offset) #draw first exon or following exons
		self.canvas.fill_between([self.smallestIndex, self.largestIndex],
		                          self.ylims['bar_min'], self.ylims['bar_max'],
		                          edgecolor=self.bgColor, facecolor=self.barColor)
        self._draw_markers() #draw markers with the default color for the diagram instance
        self._clean_axes()
    
    def show(self):
        plt.show()
        
        
if  __name__ == "__main__":
    #exons positions
    exon_pos = [97543299,97544702], [97547885,97548026], [97564044,97564188], [97658624,97658804], [97700407,97700550], [97770814,97770934], [97771732,97771853], [97839116,97839200], [97847948,97848017], [97915614,97915779], [97981281,97981497], [98015115,98015300], [98039315,98039526], [98058773,98058943], [98060614,98060722], [98144650,98144738], [98157272,98157354], [98164906,98165103], [98187065,98187227], [98205947,98206035], [98293669,98293752], [98348819,98348930], [98386439,98386615]
    exon_pos2 = [97543299,97544702], [97547885,97548026], [97564044,97564188], [97658624,97658804], [97700407,97700550], [97770814,97770934], [97771732,97771853], [97839116,97839200], [97847948,97848017], [97915614,97915779], [97981281,97981497], [98015115,98015300], [98039315,98039526], [98058773,98058943], [98060614,98060722], [98144650,98144738], [98157272,98157354], [98164906,98165103], [98187065,98187227], [98205947,98206035], [98293669,98293752], [98348819,98348930], [98386439,98386615]
    exon_pos3 = [97543299,97544702], [97547885,97548026], [97564044,97564188], [97658624,97658804], [97700407,97700550], [97770814,97770934], [97771732,97771853], [97839116,97839200], [97847948,97848017], [97915614,97915779], [97981281,97981497], [98015115,98015300], [98039315,98039526], [98058773,98058943], [98060614,98060722], [98144650,98144738], [98157272,98157354], [98164906,98165103], [98187065,98187227], [98205947,98206035], [98293669,98293752], [98348819,98348930], [98386439,98386615]
    #marker positions
    marker_pos = [97947885, 98247485]
    gene = GeneImage(exon_pos, marker_pos) #construct new class instance w/example info
    gene.show()
    multiGenes = SpliceVariantPASDiagram([exon_pos, exon_pos2, exon_pos3])
    multiGenes.show()
