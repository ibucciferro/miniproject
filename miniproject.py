#import the necessary tools
import os
from Bio import Entrez
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML



#Step 1. retrieve the transcriptones and convert to paired-end fastq files
#start by opening the file with the SRR links and reading it into a list (and then close the file)
file = open('testdata.txt', 'r')
testdata = list(file.read().strip().split('\n'))
file.close()

#using wget, get the data from NCBI
for item in testdata:
    os.system('wget ' + str(item))

#then retrieve the SRR numbers from the testdata and turn them into paired end fastq files
SRRnum = []
for item in testdata:
    SRRnum.append(item[-10:])
    
#create a command to get the paired end fastq files (fastq-dump --split-files)
pairedend_command = 'fastq-dump -I --split-files'
for i in SRRnum:
    os.system(pairedend_command + ' ' + i)

#turn the files into sample data by only using the first #### lines of the files
#use head -n ##### to only use the first #### lines of the files
#should #### = 50000? Ask in class!
pairedend1 = '_1.fastq'
pairedend2 = '_2.fastq'
for j in SRRnum:
    os.system('head -n 50000 ' + j + pairedend1 + ' > ' + j + '_sample_1.fastq')
    os.system('head -n 50000 ' + j + pairedend2 + ' > ' + j + '_sample_2.fastq')



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
    output.write('The HCMV genome (EF999921) has ' + str(counter) + ' CDS.' + '\n')
    output.close()

#finally, build the index using kallisto
kallisto_command = 'kallisto index -i cdsHCMV.idx cdsHCMV.fasta'
os.system(kallisto_command)



#ADD STEP 3 HERE



#Step 4. use bowtie2 to create an index for HCMV and save reads to the map
#begin by building the bowtie index
os.system('bowtie2-build cdsHCMV.fasta CDS_HCMV')

#then loop through the SRR numbers and map the reads
for item in SRRnum:
    os.system('bowtie1 --quiet --no-unal --al-conc Bowtie_' +item+ '.fastq -x CDS_HCMV -1 ' +item+samplepair1+ ' -2 ' +item+samplepair2+ ' -S ' +item+ '.sam')

#loop through the SRR numbers and assign to the donor name (donor 1/3 2/6dpi)
donorname = ''
for j in SRRnum:
    if j == SRRnum[0]:
        donor = 'Donor 1 (2dpi)'
    elif j == SRRnum[1]:
        donor = 'Donor 1 (6dpi)'
    elif j == SRRnum[2]:
        donor = 'Donor 3 (2dpi)'
    elif j == SRRnum[3]:
        donor = 'Donor 3 (6dpi)'
    else:
          break
    
    #determine the counter before running bowtie (and open the files before running bowtie and count the reads to determine the before count)
    bcounter1 = 0
    bcounter2 = 0
    bSRR1 = open(str(j)+pairedend1)
    bSRR2 = open(str(j)+pairedend2)
    for i in bSRR1:
        bcounter1 = bcounter1 + 1
    for k in bSRR2:
        bcounter2 = bcounter2 + 1
    beforeCounter = (bcounter1 + bcounter2)/8
    
    #do the same steps to count the reads after running bowtie 
    acounter1 = 0
    acounter2 = 0
    aSRR1 = open('Bowtie_'+j+'.1.fastq')
    aSRR2 = open('Bowtie_'+j+'.2.fastq')
    for l in aSRR1:
        acounter1 = 0
    for m in aSRR2:
        acounter2 = 0
    afterCounter = (acounter1 + acounter2)/8
    
    with open('miniProject.log', 'a') as output:
        output.write(str(donor) + ' had ' + str(beforeCounter) + ' read pairs before Bowtie2 filtering and ' + str(afterCounter) + ' read pairs after.' + '\n')
        output.close()



#Step 5. use bowtie2 output reads to assembly transcriptomes to produce 1 assembly via spades 
#start by labeling all of the SRRs to be used in the spades command
SRR1 = SRRnum[0]
SRR2 = SRRnum[1]
SRR3 = SRRnum[2]
SRR4 = SRRnum[3]

#run spades and add the spades command to the log file
spades_command = 'spades -k 55, 77, 99, 127 --only-assembler -t 2 --pe1-1 Bowtie_'+SRR1+'.1.fastq --pe1-2 Bowtie_'+SRR1+'.2.fastq --pe2-1 Bowtie_'+SRR2+'.1.fastq --pe2-2 Bowtie_'+SRR2+'.2.fastq --pe3-1 Bowtie_'+SRR3+'.1.fastq --pe3-2 Bowtie_'+SRR3+'.2.fastq --pe4-1 Bowtie_'+SRR4+'.1.fastq --pe4-2 Bowtie_'+SRR4+'.2.fastq -o SpadesAssembly/'

#write the spades command to the log file
with open('miniProject.log', 'a') as output:
    output.write(spades_command + '\n')
    output.close()



#Step 6. find the number of contigs with a length greater than 1000 bp
#open the file with the contigs in it and add each of the records to a list
contigslist = []
handlefile = open('SpadesAssembly/contigs.fasta')
for items in SeqIO.parse(handlefile, 'fasta'):
    contigslist.append(items)
handlefile.close()

#then create a counter for the number of contigs with length > 1000 bp
contigcounter = 0
for contig in contigslist:
    if len(contig.seq)>1000:
          contigcounter = contigcounter + 1
    else:
          continue
 
#write the contigcounter to the log file
with open('miniProject.log', 'a') as output:
    output.write('There are ' + str(contigcounter) + ' contigs > 1000 in the assembly.' + '\n')
    output.close()

#Step 7. find the total number of bp in all of the contigs greater than 1000 bp in length
#loop through the contig list made in step 6 and add the lengths of each of the contigs that are greater than 1000 bp to the total length
totallength = 0
for contig in contigslist:
    if len(contig.seq)>1000:
          totallength = totallength + len(contig.seq)
    else:
          continue

#write the totallength variable to the log file
with open('miniProject.log', 'a') as output:
    output.write('There are ' + str(totallength) + ' bp in the assembly.' + '\n')
    output.close()



#Step 8. take longest contig from spades and do a blast analysis 



