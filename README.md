# Three Epigenetic Clocks

We develop three epigenetic clocks: Array based human epigenetic clock, sequencing based human epigenetic clock and one mouse epigenetic clock.

## Array based human epigenetic clock

The clock is used to predict DNAm age of array based methylation data. Please run code in 'Predict_Age_Human_Array.R'.

### Example:

R_code=/code_path/Predict_Age_Human_Array.R

log_path=/path_to_log_file

log_file=/log_path/get_DNAm_age_log.txt

mkdir -p $log_path

R CMD BATCH --no-save --no-restore '--args /input_path /ref_path/' $R_code $log_file

### Description:

R command includes 2 arguments
--args $input $ref_path

1. input

input folder contains beta methylation calling .rdata files (beta.rdata). The each row of beta.rdata file represents CpG,
each column of beta.rdata file represents sample.

2. ref_path

ref_path folder contains one reference data set which will be used for DNAm age prediction.
The reference data set is the train data including 3665 samples, their chronological ages and 21723 CpGs to develop model.
The data set can be downloaded from website.

## Sequencing based human epigenetic clock

The clock is used to predict DNAm age of sequencing based methylation data. Please run code in 'Predict_Age_Human_Seq.R'.

### Example:

R_code=/code_path/Predict_Age_Human_Seq.R

log_path=/path_to_log_file

log_file=/log_path/get_DNAm_age_log.txt

mkdir -p $log_path

R CMD BATCH --no-save --no-restore '--args /input_path /ref_path/ simplified_norm GRCh37' $R_code $log_file

### Description:

R command includes 4 arguments
--args $input $ref_path $norm_method $assembly


1. input

input folder contains all methylation calling files (.txt or .txt.gz file) of each sample and chronological age file (.csv).
The name of methylation calling files should be sample_id.txt or sample_id.txt.gz.
sample_id should be the first column of chronological age file.

The methylation call file has the same format as the input file for methylKit analysis. 
Below is the example:

chrBase   chr 	base  	strand  	coverage  	freqC 	freqT

chr8.99651	  8	  99651 	F 	1 	100 	0

chr8.99661	  8	  99661	  F	  1	  100	  0

chr8.99675	  8 	99675 	F	  2 	50  	50

chr8.99679	  8	  99679	  F	  2	  50	  50


The chronological age file contains "sample_id" and "age" column showing the chronological age of each sample.


2. ref_path

ref_path folder contains two reference data sets which will be used for DNAm age prediction.
The first reference data set is the CpG annotation file of illumina 450k array.
The second reference data set is the 450k train data including 2557 samples, their chronological ages and 469279 CpGs to develop model.
The two data sets can downloaded from website.

3. norm_method
normalization methods
There are two options: simplified normalization method and quantile normalization method.

4. assembly
human genome assembly version
There are two assembly version: GRCh37 and GRCh38.





