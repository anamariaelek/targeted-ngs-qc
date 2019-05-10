# Targeted NGS QC

Inspect capture efficiency and sequencing coverage of the targeted regions from mapped 
sequencing reads.  

## Inputs

* bam file(s) with reads from targeted NGS (e.g. gene panels, exomes) mapped to the 
reference genome  
* bed file with target regions  
* bedtools genome file that defines the expected chromosome order in the input files  

Make sure to have the consistent chromosome naming in all the files (e.g. chr1 or 1).  

## Prerequisites

* bedtools  
* samtools  
* GNU parallel  
* R with the following packages: docopt, data.table, stringr, ggplot2, plotly, foreach, 
doParallel  

Make sure to have the prerequisites added to `PATH`.  

## Running the analysis

1. `mapping_qc.sh path/to/input_directory path/to/target.bed path/to/genome`  

The first script `mapping_qc.sh` takes three arguments, 
 * path to directory containing sorted and indexed bam files 
 (the files for each sample should be placed in a separate subdirectory and the path to 
 **parent** directory given as an argument!)
 * path to bed file with target regions
 * path to bedtools genome file
 
2. `Rscript parse_qc.R`  

The script can be run without any arguments, in which case it will assume the current 
directory to be both the input and output directory, containing the files produced by 
running the previous step (again, the files for multiple samples will be placed in the 
respective subdirectories, and it's the parent directory only that should be specified, 
if specifying the directories at all).
The optional arguments can be specified with the following flags:  
`-i` &nbsp;&nbsp;&nbsp;&nbsp; directory containing input files  
`-o` &nbsp;&nbsp;&nbsp;&nbsp; directory to save the output files to  
`-n` &nbsp;&nbsp;&nbsp;&nbsp; number of cores to use for paralellization 
(default: number of cores on the host-1)

## Outputs

* `flagstat_summary_on_target_reads.csv`  
Percentage of reads mapping to target regions
* `target_coverage_barplot.pdf`  
Barplots with percent of target bases at non-zero coverage
* `target_coverage_lineplot.html`  
Lineplots with percent of target bases at certain coverage
