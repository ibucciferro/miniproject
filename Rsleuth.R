#start by loading the sleuth and dplyr libraries and specifying the directory
library(sleuth)
library(dplyr)

#label directory before creating a list of paths for the kallisto results
projdir <- '.../miniProject_Isabella_Bucciferro'
sample_id <- dir(file.path(projdir, "results")) 
kall_dirs <- file.path(projdir, "results", sample_id) 

#then load the table with the sleuth data.txt file
#the sleuth data should include the accession numbers, conditions, and the sample name (donor1 or donor3)
stab <- read.table(file.path('sleuthdata.txt'), header = TRUE, stringsAsFactors = FALSE)

#then, using dplyr, select the sample and use mutate to add the path variables so that kallisto can run properly
stab <- dplyr::select(stab, sample = accession_number, condition)
stab <- dplyr::mutate(stab, path = kall_dirs)

#create sleuth object, fit the full model, and then fit the reduced model before testing
so <- sleuth_prep(stab, extra_bootstrap_summary = TRUE)
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')

#using dplyr, filter for qvalues <= 0.05 and write all of those transcripts to the log file
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all=FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)

tableoutput <- sleuth_significant[1:4]
write.tabe(tableoutput, 'miniProject.log', append = TRUE, sep='\t', dec = '.', row.names = TRUE, col.names = TRUE)
