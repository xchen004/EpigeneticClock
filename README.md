# Three Epigenetic Clocks

We develop three epigenetic clocks: Array based human epigenetic clock, sequencing based human epigenetic clock and one mouse epigenetic clock.

## Array based human epigenetic clock

The clock is used to predict DNAm age of array based methylation data. Please run code in 'Predict_Age_Human_Array.R'.

Example:

R_code=/code_path/Predict_Age_Human_Array.R
log_path=/path_to_log_file
log_file=/log_path/get_DNAm_age_log.txt
mkdir -p $log_path
R CMD BATCH --no-save --no-restore '--args /input_path /ref_path/' $R_code $log_file

Description

R command includes 2 arguments
--args $input $ref_path

1. input

input folder contains beta methylation calling .rdata files (beta.rdata). The each row of beta.rdata file represents CpG,
each column of beta.rdata file represents sample.

2. ref_path

ref_path folder contains one reference data set which will be used for DNAm age prediction.
The reference data set is the train data including 3665 samples, their chronological ages and 21723 CpGs to develop model.
The data set can be downloaded from website.

