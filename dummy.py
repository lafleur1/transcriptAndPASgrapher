def extractPredictionValuesSeparateArrays(name, dataframe, pasType = "", usingStrand = ""):
	#returns values separated by strand or not for one dataframe
	dummyPlusStrand = np.array([])
	dummyMinusStrand = np.array([])
	predName = "chr" + name
	forward, reverse = openForwardReverse("../../aparentGenomeTesting/chromosomePredictions50/", predName)
	flippedReverse = np.flip(reverse)
	if pasType != "":
		mask = dataframe['type'] == pasType
		dataframe = dataframe[mask]
	if usingStrand != "":
		mask = dataframe['strand'] == usingStrand
		dataframe = dataframe[mask]
	for index, row in dataframe.iterrows():
		if row['strand'] == "+":
			#add all on positive strand  to forward
			sliceVals = forward[startInt:endInt+ 1]
			dummyPlusStrand = np.concatenate((dummyPlusStrand,sliceVals)) 
		elif row['strand'] == "-":
			#add all on negative strand to reverse 
			sliceVals = flippedReverse[startInt:endInt+ 1]
			dummyMinusStrand = np.concatenate((dummyMinusStrand,sliceVals))
	return dummyPlusStrand, dummyMinusStrand
