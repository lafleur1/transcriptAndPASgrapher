#AML
#1/27/20
#making bash script to setup and run DeepPASTA for a given input sequence



#!/usr/bin/env python3
import os, re, signal, subprocess, sys
#import itertools to create student name combinations
import itertools
import queue
import statistics
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Alphabet import IUPAC
import matplotlib.pyplot as plt
import pandas as pd
import tensorflow as tf
tf.get_logger().setLevel('INFO')
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2' 
import time
#take sequence -> break into fasta files format
#get SS 
#process SS
#make predictions 
#extract prediction values
#delete intermediate files


'''
./RNAshapes -f <input_file_in_fa_format> -s -c 5 -t 1 -w 100 -W 100 -O 'D{%s\n}' > <output_file_in_txt_format>
python combining_substructure.py -i <output_from_the_previous_command> -o <output_file_in_txt_format>
python filtering_number_of_ss.py -n 3 -i <output_from_the_previous_command> -o <output_file_in_txt_format>
python shape_assign_per_nucleotide.py -c 3 -i <output_from_the_previous_command> -o <output_file_in_txt_format>

'''

#python DeepPASTA_polyA_site_prediction_testing.py -testSeq sample_sequence_input.hg19.fa -testSS sample_secondary_structure_input.txt  

TIMEOUT = 100000
WORKDIR = "temporary_script_folder"
#from Andrey's grader
# returns out, exits on error; throws TimeoutExpired
def bash(cmd, cwd=None, expectErrors=False):
    process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=cwd)
    out, err = process.communicate(timeout=TIMEOUT)
    if err:
        if expectErrors:
            return "ERROR"
        print(err.decode("utf-8"), "Error has occurred, exiting...")
        fail()
    return out.decode("utf-8")
 
'''
 from Bio import SeqIO
sequences = ...  # add code here
SeqIO.write(sequences, "example.fasta", "fasta")
'''


def runPASTAShortSeq(shortSeq):
	#for each 200 length kmer
	#generate temporary fasta file
	#generate SS files
	#grader_out = os.popen(command).read()
	#print (sequence[i:i+200])
	strSeq = str(shortSeq)
	print (strSeq)
	dname = "1"
	f = open("dummy.fa", "w")
	f.write(">" + dname + "\n" + strSeq + "\n")
	f.close()
	#
	comm1 = "./RNAshapes -f dummy.fa -s -c 5 -t 1 -w 100 -W 100 -O 'D{%s\n}' > output1.txt"
	p = subprocess.Popen(comm1, shell = True, stdout = subprocess.PIPE)
	print ("RUNNING RNASHAPES")
	p.wait()
	print (p.returncode)
	#
	comm2 = "python2 combining_substructure.py -i output1.txt -o output2.txt"
	p2 = subprocess.Popen(comm2, shell = True, stdout = subprocess.PIPE)
	print ("RUNNING COMBINE")
	p2.wait()
	print (p2.returncode)
	#
	comm3 = "python2 filtering_number_of_ss.py -n 3 -i  output2.txt -o output3.txt"
	print ("RUNNING FILTERING")
	p3 = subprocess.Popen(comm3, shell = True, stdout = subprocess.PIPE)
	p3.wait()
	print (p3.returncode)
	#
	comm4 = "python2 shape_assign_per_nucleotide.py -c 3 -i output3.txt -o ssFinalOutput.txt"
	p4 = subprocess.Popen(comm4, shell = True, stdout = subprocess.PIPE)
	print ("RUNNING SHAPE ASSIGN")
	p4.wait()
	print (p4.returncode)
	#run network
	runDeepPASTACommand = "python2 DeepPASTA_polyA_site_prediction_testing.py -testSeq dummy.fasta -testSS ssFinalOutput.txt -o fileOut.txt"
	print ("RUNNING DEEPPASTA")
	p5 = subprocess.Popen(runDeepPASTACommand, shell = True, stdout = subprocess.PIPE)
	p5.wait()
	print (p5.returncode)
	#add output value to sequence outputs 
	#print (output)
	print ("DONE")
	fout = open("fileOut.txt", "r")
	lines = [float(x.strip().split('\t')[1]) for x in fout.readlines()]
	return lines[0]
	
	
	
   
def runPASTA(sequence):
	#make fasta from sequence slice
	outputs = []
	current_output = 0
	for i in range(0,len(sequence) - (200-1)):
		#for each 200 length kmer
		#generate temporary fasta file
		#generate SS files
		#grader_out = os.popen(command).read()
		#print (sequence[i:i+200])
		strSeq = str(sequence[i:i+200])
		#print (strSeq)
		dname = "1"
		f = open("dummy.fa", "w")
		f.write(">" + dname + "\n" + strSeq + "\n")
		f.close()
		#
		comm1 = "./RNAshapes -f dummy.fa -s -c 5 -t 1 -w 100 -W 100 -O 'D{%s\n}' > output1.txt"
		p = subprocess.Popen(comm1, shell = True, stdout = subprocess.PIPE)
		p.wait()
		comm2 = "python2 combining_substructure.py -i output1.txt -o output2.txt"
		p2 = subprocess.Popen(comm2, shell = True, stdout = subprocess.PIPE)
		p2.wait()
		comm3 = "python2 filtering_number_of_ss.py -n 3 -i  output2.txt -o output3.txt"
		p3 = subprocess.Popen(comm3, shell = True, stdout = subprocess.PIPE)
		p3.wait()
		comm4 = "python2 shape_assign_per_nucleotide.py -c 3 -i output3.txt -o ssFinalOutput.txt"
		p4 = subprocess.Popen(comm4, shell = True, stdout = subprocess.PIPE)
		p4.wait()
		#added tensorflow warning message removal code...see if it works....
		runDeepPASTACommand = "python2 DeepPASTA_polyA_site_prediction_testing.py -testSeq dummy.fasta -testSS ssFinalOutput.txt -o fileOut.txt"
		#print ("RUNNING DEEPPASTA")
		p5 = subprocess.Popen(runDeepPASTACommand, shell = True, stdout = subprocess.PIPE)
		p5.wait()
		fout = open("fileOut.txt", "r")
		lines = [float(x.strip().split('\t')[1]) for x in fout.readlines()]
		outputs.append(lines[0])
		print ("FINISHED ONE ROUND")
		input()
	return outputs


def openBalancedPositives(name):
	positivesName = name + "BalancedPositives.csv"
	return pd.read_csv( positivesName, dtype = {"seqName": str}) 

#predicting values around the true values (same slcies used in piecewise prediction functions) 
#can only do this for chromosomes 15-Y (because those are the chromosomes not used in training DeepPASTA)
#need to save range +/- 500 nts around the true positive cluster (added another 100 buffer for running DeepPASTA (since it predicts only the middle of the sequence))
def scanPositivesForClustersNearEdge(balancedPos, chroSeq):
	for index, row in balancedPos.iterrows():
		if row['strand'] == "+":
			if row['start'] - 600 < 0:
				start = 0 
				print ("FOUND PAS TOO CLOSE TO START")
			else:
				start = row['start'] - 600
			if row['end'] + 600 > len(chroSeq) -1:
				end = len(chroSeq)-1
				print ("FOUND PAS TOO CLOSE TO END")
			else:
				end = row['end'] + 600
		else:
			#correct indexes
			rcIndexEnd = (len(chroSeq) -1) - row['start']
			if rcIndexEnd + 600 > len(chroSeq) - 1:
				end = len(chroSeq) -1
				print ("FOUND PAS TOO CLOSE TO START")
			else:
				end = rcIndexEnd + 600
			rcIndexStart = (len(chroSeq) -1) - row['end']
			if rcIndexStart - 600 < 0 :
				start = 0
				print ("FOUND PAS TOO CLOSE TO END")
			else:
				start = rcIndexStart - 600

def scanPositivesForClustersNearEdge(balancedPos, chroSeq):
	for index, row in balancedPos.iterrows():
		if row['strand'] == "+":
			if row['start'] - 600 < 0:
				start = 0 
				print ("FOUND PAS TOO CLOSE TO START")
			else:
				start = row['start'] - 600
			if row['end'] + 600 > len(chroSeq) -1:
				end = len(chroSeq)-1
				print ("FOUND PAS TOO CLOSE TO END")
			else:
				end = row['end'] + 600
			forward = predictionsForward[start:end]
			predictions = runPASTA(forward) #list of predictions for the 1200 -> 1000 slice
			#print (
		else:
			if row['start'] - 600 < 0:
				start = 0 
				print ("FOUND PAS TOO CLOSE TO START")
			else:
				start = row['start'] - 600
			if row['end'] + 600 > len(chroSeq) -1:
				end = len(chroSeq)-1
				print ("FOUND PAS TOO CLOSE TO END")
			else:
				end = row['end'] + 600
				
			#make reverse complement string (make temporary Seq object, use biopython reverse complement method)	
			forwardSliceString = predictionsRC[start:end]
			temporarySeqObject = Seq(forwardSliceString, generic_dna)
			reverseComplementString = temporarySeqObject.complement()
			predictions = runPASTA(reverseComplementString) 
			
			

#checking that there are no true positives to close to edge of chromosome to worry about 
def checkForClustersNearEdgeOfChromosome():
	#try on 1000 nts of chrY
	names = ["15","16","17","18","19","20","21","22","X","Y"]
	for chroName in names:
		print ("ON NAME: ", chroName)
		fastaLocs = "../../../aparentGenomeTesting/fastas/"
		totalName = fastaLocs + "chr" + chroName + ".fasta"
		chroSeq = SeqIO.read(totalName, "fasta")
		datasetsLocation = "../datasets/"
		fileName = datasetsLocation + "chro" + chroName + "_NegSpaces" + str(1) + "_shifted" + str(50) + "Nts"
		balancedPosDataset = openBalancedPositives(fileName)
		scanPositives(balancedPosDataset, chroSeq)







chroName = "Y"
fastaLocs = "../../../aparentGenomeTesting/fastas/"
totalName = fastaLocs + "chr" + chroName + ".fasta"
chroSeq = SeqIO.read(totalName, "fasta")
sliceY = chroSeq.seq[500000:500000 + 1000]
#sliceY2 = chroSeq.seq[500000:500000 + 300]
#sliceY3 = chroSeq.seq[5000000:5000000 + 300]
print ("TOTAL SEQ: ", len(chroSeq.seq))


print (sliceY)
print (len(sliceY))
#print (sliceY2)
#print (sliceY3)
start_time = time.time()
print ("START TIME: ")
outVals = runPASTA(sliceY)
#outVals2 = runPASTA(sliceY2)
#outVals3 = runPASTA(sliceY3)
print ("ELAPSED TIME: ", time.time() - start_time)
print ("numbe outputs: ", len(outVals))
print ("ALL OUTPUTS: ")
#plt.savefig("testrangey.png")
print (outVals)
print ("UNIQUE OUTPUTS: ")
print (set(outVals))

'''
print ("ALL OUTPUTS: ")
#plt.savefig("testrangey.png")
print (outVals2)
print ("UNIQUE OUTPUTS: ")
print (set(outVals2))

print ("ALL OUTPUTS: ")
#plt.savefig("testrangey.png")
print (outVals3)
print ("UNIQUE OUTPUTS: ")
print (set(outVals3))

plt.plot(list(range(0,len(outVals))), outVals, label = 'run1')
plt.show()

'''

