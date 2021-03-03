#start by loading the sleuth and dplyr libraries and specifying the directory
library(sleuth)
library(dplyr)

#read the table txt file that includes the table values
stab <- read.table('~/miniproject/sleuthdata.txt', header = TRUE, stringsAsFactors = FALSE)

#create sleuth object, fit the full model, and then fit the reduced model before testing
so <- sleuth_prep(stab)
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')

#using dplyr, filter for qvalues <= 0.05 and write all of those transcripts to the log file
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all=FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05) %>% dplyr::arrange(pval)

#finally, write the desired values to the log file
write.tabe(sleuth_significant, target_id, test_stat, pval, qval, file = 'tableresults.txt', quote = FALSE, row.names = FALSE, sep = '\t')
