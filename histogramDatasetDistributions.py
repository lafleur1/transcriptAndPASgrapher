#AML
#2/21/20
#making histogram of total dataset vs. balanced dataset sizes

################# graphing type of clusters in human genome ###################3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import collections  as mc

#open the total human dataset 

def openAllPositives():
	colnames = ["seqName",  "start" , "end",  "clusterID",  "avgTPM",  "strand",   "percentSupporting",   "protocolsSupporting",  "avgTPM2",   "type",   "upstreamClusters"]
	pas_stuff =pd.read_csv('atlas.clusters.hg38.2-0.bed',delimiter='\t', names = colnames) 
	return pas_stuff


def openBalancedPositives(name):
	positivesName = name + "BalancedPositives.csv"
	return pd.read_csv( positivesName, dtype = {"seqName": str}) 


def openAllPositivesAsDict():
	pas_stuff = openAllPositives()
	type_dist = pas_stuff.type.value_counts()
	ind = type_dist.index
	arr = type_dist.array
	zipped = zip(ind, arr)
	blank_dict = {}
	for z in zipped:
		blank_dict[z[0]] = z[1]
	return blank_dict
	

def openBalancedPositiesAsDict():
	names =  ["1","2","3","4","5","6","7","8","9","10","11","12","13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]
	blank_dict = {'IN':0, 'TE':0, 'IG':0, 'AI':0, 'EX':0, 'DS':0, 'AE':0, 'AU':0}
	datasetsLocation = "./datasets/"
	for name in names:
		fileName = datasetsLocation + "chro" + name + "_NegSpaces" + str(1) + "_shifted" + str(50) + "Nts"
		balanced = openBalancedPositives(fileName)
		type_dist = balanced.type.value_counts()
		ind = type_dist.index
		arr = type_dist.array
		zipped = zip(ind, arr)
		for z in zipped:
			blank_dict[z[0]] += z[1]
	print (blank_dict)
	return blank_dict
		
	
	
	
def plotBothDistributions():
    allPos = openAllPositivesAsDict()
    balancedPos = openBalancedPositiesAsDict()
    #plt.bar(*zip(*allPos.items()), label = "PolyASite2.0 Dataset")
    #plt.bar(*zip(*balancedPos.items()), label = "Balanced Dataset",  )
    X = np.arange(len(allPos))
    ax = plt.subplot(111)
    ax.bar(X-0.2, allPos.values(), width=0.2, color='r', align='center', label = "PolyASite2.0 Dataset")
    ax.bar(X, balancedPos.values(), width=0.2, color='orange', align='center',label = "Balanced Dataset" )
    plt.xticks(X, allPos.keys())
    plt.legend()
    plt.title("PolyASite2.0 and Balanced Dataset Distributions")
    plt.xlabel("Cluster Type")
    plt.ylabel("Count")
    plt.show()



'''
type_dist.plot.bar()
plt.title("Human polyA Site Clusters by Type")
plt.ylabel("Cluster Count")
plt.xlabel("Type")
plt.show()
'''
allPos = openAllPositivesAsDict()
balancedPos = openBalancedPositiesAsDict()

plotBothDistributions()



