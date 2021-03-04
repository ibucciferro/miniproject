# Comp Bio Mini Project


### Packages to have installed before running
Unix (kallisto, SPAdes, bowtie2, blast+)
- kallisto (https://pachterlab.github.io/kallisto/download)
- SPAdes (http://cab.spbu.ru/files/release3.14.0/manual.html)
- bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#obtaining-bowtie-2)

Python (Biopython)

R (sleuth)
- sleuth (https://pachterlab.github.io/sleuth/download)

### Files for the MiniProject Found in the Repo
- README_md: documentation for the miniproject repo
- testdata.txt: the text file that contains the links from NCBI (the input reads)
- miniproject.log: the output from running the miniproject.py file
- miniproject.py: the python wrapper script 
- sleuthdata.txt: the text file that contains the samples, conditions, and paths used in the sleuth package
- dbsequence.fasta: the fasta file of Betaherpesvirinae subfamily (taken from NCBI)
- Rsleuth.R: R script used to analyze the kallisto results

### How to Run the Program
- After cloning the repository using ```git clone https://github.com/ibucciferro/miniproject.git``` in Unix, move into the miniproject folder using ```cd miniproject```
- To run the script, call ```python3 miniproject.py``` in the Unix command line to run the python wrapper 
- If you want to run the full code, follow the comment in Step 1 (remove the head label). Otherwise, run the code as is.

