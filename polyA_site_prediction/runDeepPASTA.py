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
from Bio.Alphabet import IUPAC
import matplotlib.pyplot as plt

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
   
def makeFasta(sequence):
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
		print (strSeq)
		dname = "1"
		f = open("dummy.fa", "w")
		f.write(">" + dname + "\n" + strSeq + "\n")
		f.close()
		#
		comm1 = "./RNAshapes -f dummy.fasta -s -c 5 -t 1 -w 100 -W 100 -O 'D{%s\n}' > output1.txt"
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
		outputs.append(lines[0])
	print ("ALL OUTPUTS: ")
	plt.plot(list(range(0,len(outputs))), outputs)
	plt.savefig("testrangey.png")
	print (outputs)
	print ("UNIQUE OUTPUTS: ")
	print (set(outputs))
	
#try on 1000 nts of chrY
chroName = "Y"
totalName = "chr" + chroName + ".fasta"
chroSeq = SeqIO.read(totalName, "fasta")

sliceY = chroSeq.seq[110000:110000 + 1000]
makeFasta(sliceY)




