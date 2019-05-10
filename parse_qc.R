#!/usr/bin/Rscript

'Generate mapping QC reports

This script parses mapping QC output files produced using bedtools coverageBed 
and samtools flagstat, summarizes the data and creates diagnostic plots.
Usage:
parse_qc.R [-i <input_dir> -o <output_dir> -n <num_cores> -h]

Options:
-i         Directory containing input files
-o         Directory to save the output files to
-n         Number of cores to use for paralellization (default: numer of cores on the host-1)
-h         Show this screen

' -> doc

# Arguments =======================================================================================
library(docopt)
opts <- docopt(doc)
library(data.table)
library(stringr)
library(ggplot2)
library(plotly)
library(foreach)
library(doParallel)

if (opts$i) {
  input_dir <- opts$input_dir 
} else {
  input_dir <- getwd()
}

if (opts$o) {
  out_dir <- opts$output_dir
} else {
  out_dir <- input_dir
}

if (opts$n) {
  numCores <- as.integer(opts$num_cores)
} else {
  numCores <- detectCores()-1
}

message(sprintf("\nUsing %s as an input directory.", input_dir))
message(sprintf("Using %s as an output directory.", out_dir))

# Samtools flagstat ===============================================================================
message("\nParsing samtools flagstat output files:")
flagstat_files <- list.files(input_dir, pattern=".*stats$", recursive=TRUE, full.names=TRUE)

# Functions to parse flagstat output
extract_flagstat_value <- function(str) {
  values <- as.integer(stringr::str_extract_all(str, "\\d+")[[1]])[1:2]
  names(values) <- c("QC-passed","QC-failed")
  return(values[1:2])
}
extract_flagstat_description <- function(str) {
  values <- stringr::str_trim(
    stringr::str_extract_all(
      str, "(?<=\\d{1,10}[[:space:]]\\+[[:space:]]\\d{1,10}[[:space:]])[[[:alnum:]],[[:space:]]]+(\\(mapQ>=\\d+\\))*"
    )[[1]]
  )
  return(values)
}

# Parse flagstat output
flagstat_list <- lapply(flagstat_files, function(flagstat_file) {
  
  message(sprintf('Parsing input file %s', flagstat_file))
  flagstat <- readLines(flagstat_file)
  lines_to_read <- grep('^##', flagstat, invert=TRUE)
  out <- sapply(lines_to_read, function(i) extract_flagstat_value(flagstat[i]))
  out <- t(out)
  rownames(out) <- sapply(lines_to_read, function(i) extract_flagstat_description(flagstat[i]))
  data.table::data.table(out, keep.rownames='description')[,pass:=1:.N,by=description]
  
})

names(flagstat_list) <- basename(flagstat_files)
collate_flagstat_DT <- data.table::rbindlist(flagstat_list, idcol='file')

on_target_DT <- collate_flagstat_DT[description=='mapped'][
  ,.(`% reads on target`=.SD[pass==2,`QC-passed`]/.SD[pass==1,`QC-passed`]*100),by=file]

flagstat_csv <- file.path(out_dir, "flagstat_summary_on_target_reads.csv")
data.table::fwrite(
  x=on_target_DT,
  file=flagstat_csv
)
message(sprintf("\nSaved output to %s.", flagstat_csv))

# Bedtools coverageBed ============================================================================
message("\nParsing bedtools coverageBed output files:")
bed_files <- list.files(input_dir, pattern=".*bed$", recursive=TRUE, full.names=TRUE)

# Parse bedtools coverageBed output
hist_list <- list()
gghist_list <- list()
registerDoParallel(numCores)
hist_list <- foreach(bed_file=bed_files) %dopar% {
  sample <- basename(bed_file) %>% str_replace(.,'.all.coverage.report.bed','')
  message(sprintf('Parsing input file %s', bed_file))
  hist <- fread(bed_file)
  setDT(hist)
  setnames(
    hist, 
    c('chr','start','end','ref','overlapping_features','non_0_cov_bases','length','non_0_cov_bases_fraction'))
  hist
}

samples <- basename(bed_files) %>% str_replace(.,'.all.coverage.report.bed','')
names(hist_list) <- samples
hist_DT <- rbindlist(hist_list, idcol='sample')[,.(non_0_cov_bases_fraction,sample)]
setnames(hist_DT, c("x","sample"))

# Plotting parameters
theme_set(theme_light())
theme_update(
  panel.border=element_rect(fill=NULL, colour='black')
)

# Function to plot bedtools coverage
gghist <- function(hist,type='sample') {
  if (type=='sample') {
    gp <- ggplot(hist[,.(x)], aes(x)) + 
      geom_histogram(fill='grey',colour='white') +
      labs(subtitle=unique(hist$sample))
  } else if (type=='summary') {
    gp <- ggplot(unique(hist[,.(x,mean_count,se_count)]), aes(x=x,y=mean_count)) + 
      geom_bar(stat='identity', fill='grey',colour='white') +
      geom_errorbar(aes(ymin=mean_count-se_count, ymax=mean_count+se_count), width=.01) +
      labs(title='Coverage of target regions', subtitle='Mean for all samples ± SE')
  }
  gp <- gp + #facet_wrap(~chr) +
    scale_x_continuous(labels=scales::percent, expand=c(0, 0)) + 
    scale_y_continuous(trans='log10', breaks=10^(1:10), labels=scales::comma, expand=expand_scale(mult=c(0,0.1))) + 
    labs(x='Percent of bases in target region with non-zero coverage', y='Number of target regions')
  return(gp)
} 

# Plot coverage per sample
hist_plots <- lapply(samples, function(s) gghist(hist_DT[sample==s]))

# Sumarrise coverage for all samples
hist_layers <- lapply(hist_plots, layer_data)
names(hist_layers) <- samples
hist_layers_DT <- rbindlist(hist_layers, idcol='sample')
hist_layers_DT[,':='(
  mean_count=mean(count),
  se_count=(sd(count)/sqrt(.N))
  ), by=x]

# Plot summarized coverage for all samples
hist_plot_summary <- gghist(hist_layers_DT, type='summary')

# Save as pdf
cov_pdf <- file.path(out_dir,"target_coverage_barplot.pdf")
pdf(file=cov_pdf, height=5)
suppressMessages(print(hist_plot_summary))
suppressMessages(print(hist_plots))
closecon <- dev.off()
message(sprintf("\nSaved output to %s.", cov_pdf))

# Bedtools coverageBed -hist ======================================================================
message("\nParsing bedtools coverageBed -hist output files:")
report_all_files <- list.files(input_dir, pattern=".*all\\.coverage\\.report$", recursive=TRUE, full.names=TRUE)

cov_list <- lapply(report_all_files, function(file){
  message(sprintf('Parsing input file %s', file))
  df <- fread(file)
  setnames(df, c('region','cov','num_bases_at_cov','length','fraction_cov'))
  df[,cov_cumul:=1-cumsum(fraction_cov)]
  df[2:.N,.(cov,cov_cumul)]
})
samples <- str_replace(basename(report_all_files), '\\.all\\.coverage\\.report', '')
names(cov_list) <- samples
cov_DT <- rbindlist(cov_list, idcol='sample')

depth_threshold <- 300

cov_plot <- ggplot(cov_DT[cov<=depth_threshold], aes(cov, cov_cumul, colour=sample)) +
  geom_line() +
  scale_y_continuous(limits=c(0,1), breaks=c(0.1,0.2,0.5,0.8,0.9,1)) +
  scale_x_continuous(breaks = c(20,50,80,seq(100,depth_threshold,100))) +
  labs(title='Target region coverage', x='Depth', y='Percent of bases in target regions sequenced at \u2265 depth')

cov_plotly <- ggplotly(cov_plot)
cov_html <- file.path(out_dir,"target_coverage_lineplot.html")
htmlwidgets::saveWidget(widget=cov_plotly, file=cov_html, selfcontained=FALSE)
message(sprintf("\nSaved output to %s.", cov_html))

message("Done.")