 #AML 2/18/20
 #opening CSVs and generating curves
 
 
import pyensembl
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 
import sklearn.metrics as metrics
from Bio import SeqIO
import random
import time
import math
from scipy.signal import find_peaks
import sys

#format of saved CSV:
'''
columns_dict = {'name': [], 'pasType': [], 'bufferVal':[], 'spacing':[], 'minh':[], 'distance':[], 'tolerance':[], 'TruePositives':[], 'FalsePositives':[], 'FalseNegatives':[], 'TrueNegatives':[]}
'''


def openCMDataOneChromosome(name):
	return pd.read_csv( negName,  dtype = {"seqName": str}) 
