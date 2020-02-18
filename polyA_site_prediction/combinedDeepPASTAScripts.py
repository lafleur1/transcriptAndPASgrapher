import re

#functions from DeepPASTA shape_asign_per_nucleotide.py 

def replaceSupportingCharacter(finalShape):
	for i in range(len(finalShape)):
		if finalShape[i] == '<':
			finalShape[i] = 'L'
		elif finalShape[i] == '>':
			finalShape[i] = 'R'
	return finalShape

def updateForIAndM(finalShape, pattern):
	firstIndex = -1
	for i in range(len(finalShape)):
		if finalShape[i] == '(' and finalShape[i+1] != '(':
			firstIndex = (i+1)
		elif (finalShape[i] == '.' or finalShape[i] == '>') and finalShape[i+1] == ')' and firstIndex != -1:
			lastIndex = i
			modified = False
			if len(re.findall(pattern, ''.join(finalShape)[firstIndex:lastIndex])) == 1:
				l = firstIndex
				while (l <= lastIndex):
					if finalShape[l] == '<':
						finalShape[l] = 'L'
					elif finalShape[l] == '>':
						finalShape[l] = 'R'
					elif finalShape[l] == '.':
						finalShape[l] = 'I'
					modified = True
					l = l + 1
			elif len(re.findall(pattern, ''.join(finalShape)[firstIndex:lastIndex])) > 1:
				l = firstIndex
				while (l <= lastIndex):
					if finalShape[l] == '<':
						finalShape[l] = 'L'
					elif finalShape[l] == '>':
						finalShape[l] = 'R'
					elif finalShape[l] == '.':
						finalShape[l] = 'M'
					modified = True
					l = l + 1
			if modified == True:
				l = firstIndex - 1
				m = lastIndex + 1
				while (finalShape[l] == '(' and finalShape[m] == ')'):
					finalShape[l] = '<'
					finalShape[m] = '>'
					l = l - 1
					m = m + 1
			firstIndex = -1
		elif finalShape[i] == '(' and finalShape[i+1] == '(':
			firstIndex = -1
	return finalShape

def assignUnpairedPosition(shape):
	for i in range(len(shape)):
		if shape[i] == '.':
			shape[i] = 'U'
	return shape

def shapeAtEachPosition(abstractShape):
	beforeOrAfterAnyParanthesis = False
	finalShape = ['*'] * len(abstractShape)
	lastParanthesisIndex = 0
	for i in range(len(abstractShape)):
		if abstractShape[i] == ')':
			lastParanthesisIndex = i
		finalShape[i] = abstractShape[i]
	lastbracket = ''
	firstIndex = 0
	lastIndex = 0
	for i in range(len(abstractShape)):
		if beforeOrAfterAnyParanthesis == False and abstractShape[i] == '(':
			beforeOrAfterAnyParanthesis = True
		if beforeOrAfterAnyParanthesis == True and i > lastParanthesisIndex:
			beforeOrAfterAnyParanthesis = False
			
		if abstractShape[i] == "." and beforeOrAfterAnyParanthesis == False:
			finalShape[i] = 'E'
		elif (abstractShape[i] == '(' or abstractShape[i] == ')') and abstractShape[i+1] == '.':
			lastbracket = abstractShape[i]
			firstIndex = (i+1)
		elif abstractShape[i] == '.' and (abstractShape[i+1] == ')' or abstractShape[i+1] == '('):
			lastIndex = i
			if lastbracket == '(' and abstractShape[i+1] == ')':
				l = firstIndex
				while (l <= lastIndex):
					finalShape[l] = 'H'
					l = l+1
				l = firstIndex - 1
				m = lastIndex + 1
				while (abstractShape[l] == '(' and abstractShape[m] == ')'):
					finalShape[l] = '<'
					finalShape[m] = '>'
					l = l - 1
					m = m + 1
	
			lastbracket = ''
	
	finalShape = updateForIAndM(finalShape, '<+\w+>+')
	count = 0
	while ('(' in ''.join(finalShape)) or (')' in ''.join(finalShape)):
		newFinalShape = updateForIAndM(finalShape, '<+\w*>+')
		if newFinalShape == finalShape:
			if count < 30:
				count = count + 1
			else:
				break
		finalShape = newFinalShape
	finalShape = replaceSupportingCharacter(finalShape)	
	if '.' in ''.join(finalShape):
		finalShape = assignUnpairedPosition(finalShape)

	return finalShape


#RNA shapes: https://academic.oup.com/bioinformatics/article/22/4/500/184565 

#FIRST PY: COMBINING_SUBSTRUCTURE ######################################################################

#combines first portion of 1-100 

def processRNAShapesOutput(inputFile):	
	firstSubStructure = []
	secondSubStructure = []
	outputList = []
	current = 'none'
	if inputFile != "":
		for line in open(inputFile):
			if ">" in line: #first line in file (fasta name)
				###
				if current == 'second': #enters if there is another fasta, not important for me
					for first in firstSubStructure:
						for second in secondSubStructure:
							outputList.append(first+second)
					current = 'none'
				###
				outputList.append(line[:-1])
				firstSubStructure = []
				secondSubStructure = []
			elif ('101' in line) and ('200' in line): 
				current = 'second'
			elif ('1' in line) and ('100' in line): #second line in file
				current = 'first' #set that it is in 1-100 area
			elif ("(" in line or ")" in line or "." in line):
				if current == 'first':
					firstSubStructure.append(line[:-1]) #removing the end line character between lines
				else:
					secondSubStructure.append(line)

	if current == 'second':
		for first in firstSubStructure:
			for second in secondSubStructure:
				outputList.append(first+second)
		current = 'none'


	#now have the outputs of combining_substructure 1 
	#########################################################################################################
	#SECOND PY: FILTERING_NUMBER_OF_SS #########################
	passedOn = outputList[0:4]
	#literally all this file does is keep the first 3 lines of the previous created file
	###########################################################################################################
	#THIRD PY: Shape assign per nucleotide 

	count = 0
	outputList3 = []
	lastShape = ''
	for line in passedOn:
		if ">" in line:
			outputList3.append(line)
			count = 0
		elif "(" in line or ")" in line or "." in line:
			lastShape = ''.join(shapeAtEachPosition(line))
			outputList3.append(lastShape[:-1])
			count = count + 1
	return outputList3


print (processRNAShapesOutput("output1.txt"))
