# Three Epigenetic Clocks

We develop three epigenetic clocks: Array based human epigenetic clock, sequencing based human epigenetic clock and one mouse epigenetic clock. Please go to google drive to download reference data sets.
Google drive link:
https://drive.google.com/drive/folders/1-Q3iTEQeLOjSLTqZrmGOhuxJjafRydy2?usp=sharing

Please do not change the name of reference data sets.

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
--args input ref_path

1. input

input folder contains beta methylation calling .rdata files (beta.rdata). The each row of beta.rdata file represents CpG,
each column of beta.rdata file represents sample.

2. ref_path

ref_path folder contains one reference data set which will be used for DNAm age prediction.
The reference data set is the train data including 3665 samples, their chronological ages and 21723 CpGs to develop model.
The data set can be downloaded from above google drive link.

## Sequencing based human epigenetic clock

The clock is used to predict DNAm age of sequencing based human methylation data. Please run code in 'Predict_Age_Human_Seq.R'.

### Example:

R_code=/code_path/Predict_Age_Human_Seq.R

log_path=/path_to_log_file

log_file=/log_path/get_DNAm_age_log.txt

mkdir -p $log_path

R CMD BATCH --no-save --no-restore '--args /input_path /ref_path/ simplified_norm GRCh37' $R_code $log_file

### Description:

R command includes 4 arguments
--args input ref_path norm_method assembly


1. input

input folder contains all methylation calling files (.txt or .txt.gz file) of each sample and chronological age file (.csv).
The name of methylation calling files should be sample_id.txt or sample_id.txt.gz.
sample_id should be the first column of chronological age file.

The methylation call file has the same format as the input file for methylKit analysis. 
Below is the example:

chrBase      chr		base  	strand  	coverage  	freqC 	freqT

chr8.99651	  8	  99651 	F 	1 	100 	0

chr8.99661	  8	  99661	  F	  1	  100	  0

chr8.99675	  8 	99675 	F	  2 	50  	50

chr8.99679	  8	  99679	  F	  2	  50	  50


The chronological age file contains "sample_id" and "age" column showing the chronological age of each sample.


2. ref_path

ref_path folder contains two reference data sets which will be used for DNAm age prediction.
The first reference data set is the CpG annotation file of illumina 450k array.
The second reference data set is the 450k train data including 2557 samples, their chronological ages and 469279 CpGs to develop model.
The two data sets can downloaded from above google drive link.

3. norm_method

There are two options: simplified normalization method and quantile normalization method.

4. assembly

There are two assembly version: GRCh37 and GRCh38.


## Mouse epigenetic clock

The clock is used to predict DNAm age of sequencing based mouse methylation data. Please run code in 'Predict_Age_Mouse.R'.

### Example:

R_code=/code_path/Predict_Age_Mouse.R

log_path=/path_to_log_file

log_file=/log_path/get_DNAm_age_log.txt

mkdir -p $log_path

R CMD BATCH --no-save --no-restore '--args /input_path /ref_path/ 10 quantile_norm GRCh37' $R_code $log_file

### Description:

R command includes 5 arguments
--args input ref_path cov_min norm_method assembly


1. input

input folder contains all methylation calling files (.txt or .txt.gz file) of each sample and chronological age file (.csv).

The methylation call file has the same format as the input file for methylKit.

Below is the example:

chrBase	chr	base	strand	coverage	freqC	freqT

chr5.3012454	5	3012454	F	1	100	0

chr5.3020804	5	3020804	R	1	100	0

chr5.3021139	5	3021139	F	48	95.83	4.17

chr5.3021183	5	3021183	F	47	68.09	31.91

The chronological age information contains "age" column showing the chronological age of each sample.


2. ref_path

ref_path folder contains one reference data set which will be used for DNAm age prediction.
The reference data set is the RRBS train data including 280 samples, their chronological ages and 1180864 CpGs to develop model.
The data set can downloaded from above google drive link..


3. cov_min

minimum read coverage of CpG to be included in age prediction. We used 10 as an example.

4. norm_method

There are two options: simplified normalization method and quantile normalization method.

5. assembly

There is only one assembly version for mouse: GRCm38.



