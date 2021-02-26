#import the necessary tools
import os
from Bio import Entrez
from Bio import SeqIO
from Bio.Blast import NCBIWW
from Bio.Blast import NCBIXML

#start off by making the directory (mkdir) and copying the files to it
os.system('mkdir miniProject_Isabella_Bucciferro')
miniprojdir = '.../miniProject_Isabella_Bucciferro'
os.system('cp testdata.txt ' + miniprojdir)
os.system('cp Rsleuth.R ' + miniprojdir)
#then, change the directory to the miniproject directory (to move into it)
os.chdir(miniprojdir)
os.system('touch miniProject.log')



#Step 1. retrieve the transcriptones and convert to paired-end fastq files
#start by opening the file and reading it into a list (and then close the file)
file = open(miniprojdir + 'testdata.txt', 'r')
testdata = list(file.read().strip().split('\n'))
file.close()

#using wget, get the data from NCBI
for item in testdata:
    os.system('wget ' + str(item))

#then retrieve the SRR numbers from the testdata and turn them into paired end fastq files
SRRnum = []
for item in testdata:
    SRRnum.append(item[-10:])
    
#create a command to get the paired end fastq files
pairedend_command = 'fastq-dump -I --split-files'
for i in SRRnum:
    os.system(pairedend_command + ' ' + i)

#turn the files into sample data by only using the first #### lines of the files
#use head -n ##### to only use the first #### lines of the files
pairedend1 = '_1.fastq'
pairedend2 = '_2.fastq'
for j in SRRnum:
    os.system('head -n 50000 ' + j + pairedend1 + ' > " + j + '_sample_1.fastq')
    os.system('head -n 50000 ' + j + pairedend2 + ' > " + j + '_sample_2.fastq')



#Step 2. extract the CDS features from GenBank format (and build the index using kallisto)
#use efetch the search the NCBI databases (using the id = EF999921)
Entrez.email = 'ibucciferro@luc.edu'
handle = Entrez.efetch(db='nucleotide',id='EF999921',rettype='gbwithparts', retmode='text')
item = SeqIO.read(handle, 'gb')
outputfile = open("cdsHCMV.fasta', 'w')

#loop through item.features (the CDS features) and write them to the file
counter = 0
for j in item.features:
    if j.type == 'CDS':
        #if it is a CDS, add to the counter and determine the name and sequence (write to file)
        counter += counter
        name = j.qualifiers['protein_id']
        sequence = j.extract(item.seq)
        outputfile.write(">" + str(name) + '\n' + str(sequence) + '\n')
outputfile.close()

#write the HCMV CDS number to the log file
#store the intended output in a variable and call it when writing to the miniproject.log
numCDS = 'The HCMV genome (EF999921) has ' + str(counter) + ' CDS. \n'
os.system('echo ' + " '" + numCDS + "' >> miniProject.log")

#finally, build the index using kallisto
kallisto_command = 'kallisto index -i cdsHCMV.idx cdsHCMS.fasta'
os.system(kallisto_command)



#Step 3. quantify the TPM of each CDS using kallisto and run the results using the R package sleuth (while will be in an R file to be called at the end)
#start by creating a command to call kallisto and label the sample paired end files
samplepair1 = '_sample_1.fastq'
samplepair2 = '_sample_2.fastq'
call_kallisto = 'time kallisto quant -i HCMV_cds.idx -o"

#loop through the SRR files and call kallisto
for record in SRRnum:
    os.system(call_kallisto + miniprojdir + 'results/' + record + ' -b 30 -t 4 ' + record + samplepair1 + ' ' + record + samplepair2)
    
#then run the R script for sleuth! 
os.system('Rscript Rsleuth.R')



#Step 4. use bowtie2 to create an index for HCMV and save reads to the map

#Step 5. use bowtie2 output reads to assembly transcriptomes to produce 1 assembly via spades 



#Step 6. find the number of contigs with a length greater than 1000 

#Step 7. find the total number of bp in all of the contigs greater than 1000 bp in length

#Step 8. take longest contig from spades and do a blast analysis 

