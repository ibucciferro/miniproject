#import the necessary tools
import os
from Bio import Entrez
from Bio import SeqIO
from Bio.Blast import NCBIWW
from Bio.Blast import NCBIXML

#start off by setting up the directory, copying files to it, and creating the log file
os.system('mkdir miniProject_Isabella_Bucciferro')
os.system('touch miniProject.log')

miniprojdir = '.../miniProject_Isabella_Bucciferro'

#1. retrieve the transcriptones and convert to paired-end fastq files
#start by opening the file and reading it into a list
file = open(miniprojdir + 'testdata.txt', 'r')
testdata = list(file.read().strip().split('\n'))



#2. extract the CDS features from GenBank format (write to log file)

#3. quantify the TPM of each CDS using kallisto (write to log file)
callkallisto = 'time kallisto quant -i HCMV_cds.idx -o"

#close the files
f.close()