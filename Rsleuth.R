library(sleuth)
library(dplyr)
library(data.table)

#start by reading the data table txt file
stab <- read.table("sleuthdata.txt", header = TRUE, stringsAsFactors = FALSE)

#create sleuth object, fit the full model, and then fit the reduced model before testing

so <- sleuth_prep(stab)
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')

#using dplyr, filter for qvalues <= 0.05 and write all of those transcripts to the log file
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all=FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval<=0.05) %>% dplyr::arrange(pval)

#write to a file that can be used to write to the log file
parameters <- dplyr::select(sleuth_significant, target_id, test_stat, pval, qval)
write.table(parameters, file = 'sleuth_results.txt', quote = FALSE, row.names = FALSE, sep = '\t')
