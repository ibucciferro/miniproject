#start by loading the sleuth and dplyr libraries and specifying the directory
library(sleuth)
library(dplyr)
miniprojdir <- '.../miniProject_Isabella_Bucciferro'

#
sample_id <- dir(file.path(miniprojdir))
kall_dir <- file.path(miniprojdir, sample_id)

