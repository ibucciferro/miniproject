#import the necessary tools
import os
from Bio import Entrez
from Bio import SeqIO
#import logging to be used to add the output to the log file for steps 8
import logging





#start by creating output file folder
os.system("mkdir miniProject_Isabella_Bucciferro")
#Step 1. retrieve the transcriptones and convert to paired-end fastq files
#start by opening the file with the SRR links and reading it into a list (and then close the file)
file = open('testdata.txt', 'r')
testdata = list(file.read().strip().split('\n'))
file.close()
#using wget, get the data from NCBI
for item in testdata:
    os.system('wget ' + str(item))
#then retrieve the SRR numbers from the testdata and turn them into paired end fastq files
SRRnum = ["SRR5660030","SRR5660033","SRR5660044","SRR5660045"]
    
#create a command to get the paired end fastq files (fastq-dump --split-files)
pairedend_command = 'fastq-dump -I --split-files'
for i in SRRnum:
    os.system(pairedend_command + ' ' + i+'.1')

#turn the files into sample data by only using the first 50000 lines of the files
#can change how much of the data is run through the wrapper by removing the "head -n 50000"
for j in SRRnum:
    os.system('head -n 50000 ' + j + '.1_1.fastq' + ' > ' + j + '_s_1.fastq')
    os.system('head -n 50000 ' + j + '.1_2.fastq' + ' > ' + j + '_s_2.fastq')
    
    
    
    
#Step 2. extract the CDS features from GenBank format (and build the index using kallisto)
#use efetch the search the NCBI databases (using the id = EF999921)
Entrez.email = 'ibucciferro@luc.edu'
handle = Entrez.efetch(db='nucleotide',id='EF999921',rettype='gbwithparts', retmode='text')
item = SeqIO.read(handle, 'gb')
#loop through item.features (the CDS features) and write them to the file
outputfile = open('cdsHCMV.fasta', 'w')
counter = 0
for j in item.features:
    if j.type == 'CDS':
        #if it is a CDS, add to the counter and determine the name and sequence (write to file)
        counter = counter + 1
        name = j.qualifiers['protein_id']
        sequence = j.extract(item.seq)
        outputfile.write(">" + str(name) + '\n' + str(sequence) + '\n')
outputfile.close()

#write the HCMV CDS number to the log file
#store the intended output in a variable and call it when writing to the miniproject.log
with open('miniProject.log', 'a') as output:
    output.write('The HCMV genome (EF999921) has ' + str(counter) + ' CDS.' + '\n'+ '\n')
    output.close()

#finally, build the index using kallisto
kallisto_command = 'kallisto index -i cdsHCMV.idx cdsHCMV.fasta'
os.system(kallisto_command)





#Step 3. quantify the TPM of each CDS using kallisto and run the results using the R package sleuth (while will be in an R file to be called at the end)
#start by creating a command to call kallisto
kallisto_call = 'time kallisto quant -i cdsHCMV.idx -o'
os.system('mkdir kallisto_results')

#create directories for each of the SRR files in the kallisto_results folder
os.system('cd kallisto_results')
for record1 in SRRnum:
    os.system('mkdir ' +record1)
os.system('cd ..')

#loop through the SRR files and call kallisto before running the R script for sleuth
for record in SRRnum: 
    os.system('time kallisto quant -i cdsHCMV.idx -o kallisto_results/'+record+' -b 30 -t 2 '+ record + '_s_1.fastq ' + record + '_s_2.fastq')

#then change the directory, run the R script, and then change back to the old directory!
os.system('cd ..')
os.system('Rscript Rsleuth.R')

#using the logging system, write the output to the log file
logging.basicConfig(filename = 'miniProject.log', level = logging.INFO)
#open the sleuth results file and read the lines into the log file
FDRfile = open('sleuth_results.txt', 'r')
for newline in FDRfile:
    logging.info(newline)





#Step 4. use bowtie2 to create an index for HCMV and save reads to the map
#begin by building the bowtie index and running the bowtie command
os.system('bowtie2-build cdsHCMV.fasta CDS_HCMV')

#then loop through the SRR numbers and map the reads using a bowtie command line call
for item in SRRnum:
    os.system('bowtie2 --quiet -x CDS_HCMV -1 ' + item + '_s_1.fastq -2 ' + item + '_s_2.fastq -S ' + item + 'mapping.sam --al-conc-gz ' + item + '_mapped_%.fastq')

#finally, write the output to the log file
for item in SRRnum:
    os.system('wc -l < ' + item+ '_s_1.fastq >> beforereadfile.txt')
    os.system('wc -l < ' + item+ '_mapped_1.fastq >> afterreadfile.txt')

#then go through the before/after files and determine the lengths
beforebow = open('beforereadfile.txt').read().rstrip()
bbow = beforebow.split('\n')
beforelength = list(map(int, bbow))
beforelength = [int(length/4) for length in beforelength]

afterbow = open('afterreadfile.txt').read().rstrip()
abow = afterbow.split('\n')
afterlength = list(map(int, abow))
afterlength = [int(length/4) for length in afterlength]

#write the different before and after lengths to the log file
logfile1 = open('miniProject.log', 'a')
logfile1.write('Donor 1 (2dpi) had ' + str(beforelength[0])+ ' read pairs before Bowtie2 filtering and '+ str(afterlength[0])+ ' read pairs after.' + '\n')
logfile1.write('Donor 1 (6dpi) had ' + str(beforelength[1])+ ' read pairs before Bowtie2 filtering and '+ str(afterlength[1])+ ' read pairs after.' + '\n')
logfile1.write('Donor 3 (2dpi) had ' + str(beforelength[2])+ ' read pairs before Bowtie2 filtering and '+ str(afterlength[2])+ ' read pairs after.' + '\n')
logfile1.write('Donor 3 (6dpi) had ' + str(beforelength[3])+ ' read pairs before Bowtie2 filtering and '+ str(afterlength[3])+ ' read pairs after.' + '\n')
logfile1.close()






#Step 5. use bowtie2 output reads to assembly transcriptomes to produce 1 assembly via spades 

#run spades and add the spades command to the log file
spades_command = 'spades -k 55,77,99,127 -t 2 --only-assembler --pe1-1 SRR5660030_mapped_1.fastq --pe1-2 SRR5660030_mapped_2.fastq --pe2-1 SRR5660033_mapped_1.fastq --pe2-2 SRR5660033_mapped_2.fastq --pe3-1 SRR5660044_mapped_1.fastq --pe3-2 SRR5660044_mapped_2.fastq --pe4-1 SRR5660045_mapped_1.fastq --pe4-2 SRR5660045_mapped_2.fastq -o SpadesAssembly/'
os.system(spades_command)


#write the spades command to the log file
with open('miniProject.log', 'a') as output:
    output.write('\n' + spades_command + '\n' + '\n')
    output.close()
    
    
    
    
    
#Step 6. find the number of contigs with a length greater than 1000 bp
#open the file with the contigs in it and add each of the records to a list
contigslist = []
handle1 = open('SpadesAssembly/contigs.fasta')
for items in SeqIO.parse(handle1, 'fasta'):
    contigslist.append(items)
handle1.close()

#then create a counting loop for the number of contigs with length > 1000 bp
longcontiglist = []
for contig in contigslist:
    if len(str(contig.seq))>1000:
          longcontiglist.append(str(contig.seq))

    else:
          continue
 
#write the contigcounter to the log file
with open('miniProject.log', 'a') as output:
    output.write('There are ' + str(len(longcontiglist)) + ' contigs > 1000 in the assembly.' + '\n'+ '\n')
    output.close()
    
    
    
    
    
#Step 7. find the total number of bp in all of the contigs greater than 1000 bp in length
#loop through the contig list made in step 6 and add the lengths of each of the contigs that are greater than 1000 bp to the total length
totallength = 0
for contig in longcontiglist:
    totallength = totallength + len(contig)
    
#write the totallength variable to the log file
with open('miniProject.log', 'a') as output:
    output.write('There are ' + str(totallength) + ' bp in the assembly.' + '\n'+ '\n')
    output.close()
    
    
    
    
    
#Step 8. take longest contig from spades and do a blast analysis
#start by sorting the list by length and finding the longest contig (should be the last one in the list)
longcontiglist.sort(key = len)
longestcon = longcontiglist[-1]

#then loop through the contigs.fasta and determine the fasta id and write it to a file
longestid = ''
file12 = open('SpadesAssembly/contigs.fasta')
for items in SeqIO.parse(file12, 'fasta'):
    sequence = str(items.seq)
    if sequence == longestcon:
        longestid = str(items.id)
#write the ID and sequence of the longest contig to the longestcontigfile
filelongcon = open('longestcontigfile.txt','w')
filelongcon.write('>' + longestid+ '\n')
filelongcon.write(str(longestcon))
file12.close()
filelongcon.close()

#create the blast database command and call it 
makeblastdb = 'makeblastdb -in dbsequence.fasta -out betaherpesvirinaedb -title betaherpesvirinaedb -dbtype nucl'
os.system(makeblastdb)

#create the blast command and call it to run the blastn
blast_command = 'blastn -query longestcontigfile.txt -db betaherpesvirinaedb -out blastn_results.csv -outfmt "10 sacc pident length qstart qend sstart send bitscore evalue stitle"'
os.system(blast_command)

#then parse the blast results and grab the 10 hits
import csv
#create a function to parse the blast hits
def parse_blast(filename, headers):
    newlist = []
    blast_results = open(filename, 'r')
    rows = csv.DictReader(blast_results, headers, delimiter=',')
    for row in rows:
        newlist.append(row)
    blast_results.close()
    return newlist

#label the file/headers, call the function and grab the top 10 hits
headernames = ['sacc', 'pident', 'length', 'qstart', 'qend', 'sstart', 'send', 'bitscore', 'evalue', 'stitle']
parsedblast = parse_blast("blastn_results.csv", headernames)
topten = parsedblast[:10]

#now, write the output to the log file
logging.basicConfig(filename = 'miniProject.log', level = logging.INFO)
#tab delimit the headers
tabhead = '\t'.join(headernames)
logging.info(tabhead)
#then loop through the topten list and write the output to the files
for item in topten:
    iteminfo = list(item.values())
    iteminfo = [str(i) for i in iteminfo]
    finaloutput = '\t'.join(iteminfo)
    logging.info(finaloutput)


#finally, move all of the output files into the miniProject_NAME file (miniProject_Isabella_Bucciferro)
for num in SRRnum:
    os.system('mv '+num+'.1 ~/miniProject_Isabella_Bucciferro')
    os.system('mv '+num+'.1.1 ~/miniProject_Isabella_Bucciferro')
    os.system('mv '+num+'.1_1.fastq ~/miniProject_Isabella_Bucciferro')
    os.system('mv '+num+'.1_2.fastq ~/miniProject_Isabella_Bucciferro')
    os.system('mv '+num+'_mapped_1.fastq ~/miniProject_Isabella_Bucciferro')
    os.system('mv '+num+'_mapped_2.fastq ~/miniProject_Isabella_Bucciferro')
    os.system('mv '+num+'_s_1.fastq ~/miniProject_Isabella_Bucciferro')
    os.system('mv '+num+'_s_2.fastq ~/miniProject_Isabella_Bucciferro')
    os.system('mv '+num+'mapping.sam ~/miniProject_Isabella_Bucciferro') 
os.system('mv kallisto_results miniProject_Isabella_Bucciferro/kallisto_results')
os.system('mv SpadesAssembly miniProject_Isabella_Bucciferro/SpadesAssembly')
os.system('mv sleuth_results.txt ~/miniProject_Isabella_Bucciferro')
os.system('mv longestcontigfile.txt ~/miniProject_Isabella_Bucciferro')
os.system('mv cdsHCMV.idx ~/miniProject_Isabella_Bucciferro')
os.system('mv cdsHCMV.fasta ~/miniProject_Isabella_Bucciferro')
os.system('mv blastn_results.csv ~/miniProject_Isabella_Bucciferro')
os.system('mv afterreadfile.txt ~/miniProject_Isabella_Bucciferro')
os.system('mv beforereadfile.txt ~/miniProject_Isabella_Bucciferro')
os.system('mv betaherpesvirinaedb.nhr ~/miniProject_Isabella_Bucciferro')
os.system('mv betaherpesvirinaedb.nin ~/miniProject_Isabella_Bucciferro')
os.system('mv betaherpesvirinaedb.nsq ~/miniProject_Isabella_Bucciferro')

#end of the python wrapper