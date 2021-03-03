#start by loading the sleuth and dplyr libraries and specifying the directory
library(sleuth)
library(dplyr)

#label directory before creating a list of paths for the kallisto results
projdir <- '.../miniProject_Isabella_Bucciferro'
projdir <- paste(projdir, '/kallisto_results', sep='')
sample_id <- dir(file.path(projdir)) 
kall_dirs <- file.path(projdir, sample_id) 

#then, create the table by setting the conditions and making a data frame before adding the column names
condition <- c("2 dpi", "6 dpi", "2 dpi", "6 dpi")
table1 <-as.data.frame(cbind(sample_id, condition, kall_dirs))
colnames(table1) <-c("sample", "condition", "path")

#then make sure to convert all of the characters
table1$sample<-as.character(table1$sample) 
table1$condition<-as.character(table1$condition)
table1$path<-as.character(table1$path)

#create sleuth object, fit the full model, and then fit the reduced model before testing
so <- sleuth_prep(table1, extra_bootstrap_summary = TRUE)
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')

#using dplyr, filter for qvalues <= 0.05 and write all of those transcripts to the log file
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all=FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)

tableoutput <- sleuth_significant[1:4]
write.tabe(tableoutput, 'miniProject.log', append = TRUE, sep='\t',  row.names = FALSE, col.names = TRUE, quote = FALSE)
