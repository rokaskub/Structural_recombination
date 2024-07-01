library(data.table)
library(Rsamtools)
library(circlize)
library(ggplot2)
library(viridis)
library(ggraph)
library(igraph)
library(plyr)
library(cli)

# SSA Annotation file
annotation <- fread("../ReferenceMaterials/Annotations/SSA_window_annotation_updated.csv")
table(annotation$Note)
annotation[, repeat_id := paste(Note, Number, sep = "_")]
annotation[, fragment_size := (end - start) + 1] 
annotation <- annotation[!is.na(start)]
#setkey(annotation, seqname, R1_site_start, R1_site_end)
setkey(annotation, seqname, window.start, window.end)

# Load samples information
names <- fread("../ReferenceMaterials/Annotations/Samples.csv", header = T)

# Sample
names[SampleName == "WT_control"]

readCount <- function(bamPath, annotation){
  cli_alert_info("Extracting reads from {.file {basename(bamPath)}}")
  gr <- makeGRangesFromDataFrame(annotation[, c("seqname", "Note", "R1_site_start", "R1_site_end", "repeat_id")], 
                                 start.field = "R1_site_start", end.field = "R1_site_end", strand.field = "strand", 
                                 keep.extra.columns = T)
  seqlevels(gr) <- "NC_037304.1:258128-270303"
  bamFile <- BamFile(bamPath)
  params <- ScanBamParam(which = gr, what = c("qname", "pos", "mpos", "isize"))
  aln <- scanBam(bamFile, param = params)
  names(aln) <- annotation$repeat_id
  res <- rbindlist(aln, idcol = "repeat_id", use.names = T)
  res[, seqname := "NC_037304.1"]
  res <- res[isize >= 0]
  return(res)
}

# bamFiles
is1000 <- list.files("../1_Mapping_bwa_full/bwa_strict/insert_size_below_1000", "*.ChrM_258128_270303.MD.IS1000.filtered.bam$", full.names = T)
names(is1000) <- gsub(".ChrM_258128_270303.MD.IS1000.filtered.bam", "", basename(is1000))

# Parse bamFiles
is_inf1000.dt <- rbindlist(lapply(is1000, readCount, annotation = annotation), idcol = "SampleName")

fwrite(is_inf1000.dt, "../Results/SSA_PE_reads_IS_inf1000.raw_reads.csv.gz", sep = ",")

is_inf1000.dt[, `:=` ("R1.start" = pos, "R1.end" = pos, "R2.start" = mpos, "R2.end" = mpos)]
setkey(is_inf1000.dt, seqname, R1.start, R1.end)

repeats <- annotation$repeat_id
names(repeats) <- repeats

# Keep reads falling in annotation intervals R1start-R1end & R2start-R2end
by_repeat <- rbindlist(lapply(repeats, function(id){
  cli_alert_info(id)
  R <- annotation[repeat_id == id]
  is_inf1000.dt[pos >= R$R1_site_start & pos <= R$R1_site_end & mpos >= R$R2_site_start & mpos <= R$R2_site_end]
}), idcol = "repeat_id")

summary.dt <- by_repeat[, .N, by = c("repeat_id", "SampleName")]
unique(summary.dt$SampleName)
summary.dt <- merge(summary.dt, names[, c("SampleName", "Sample", "Genotype", "Construct", "Suffix", "Name")], by = "SampleName", all.x = T)

fwrite(summary.dt, "../Results/SSA_PE_reads_IS_inf1000.count.csv.gz")

# Load back if already computed
summary.dt <- fread("../Results/SSA_PE_reads_IS_inf1000.count.csv.gz")

# Fraction of data --------------------------------------------------------

# https://stackoverflow.com/questions/40851328/compute-area-under-density-estimation-curve-i-e-probability
# Create the function
getFraction <- function(v, size){
  d_fun <- ecdf(v)
  round(d_fun(size), 3)
}

repeats <- annotation$repeat_id
names(repeats) <- repeats

# Loop over repeats, use the fragment size calculated (size + 2x150bp)
by_repeat <- rbindlist(lapply(repeats, function(id){
  fragment_size <- annotation[repeat_id == id]$fragment_size
  Note <- annotation[repeat_id == id]$Note
  try(is_inf1000.dt[, list("Fraction" = getFraction(isize, size = fragment_size)), by = SampleName])
}), idcol = "repeat_id")

by_repeat <- merge(by_repeat, names, by = "SampleName", all.x = T)
by_repeat[, Norm := 1 - Fraction]

# We should get 3 columns:  Repeat_ID, Sample_Name, Fraction
fwrite(x = by_repeat, file = "../Results/Fraction_SSA_inf1000.csv", sep = ",")

# Load back
by_repeat <- fread("../Results/Fraction_SSA_inf1000.csv")

# Get wider format
fraction.dt <- setDT(tidyr::pivot_wider(data = by_repeat, id_cols = "Sample", names_from = repeat_id, values_from = Fraction))
fwrite(x = fraction.dt, file = "../Results/REPEATS_fragment_size.proportion.wide.csv", sep = ",")

library(pheatmap)
mat <- as.matrix(fraction.dt, rownames=T)
pheatmap(t(mat), cluster_rows = TRUE)

# Annotate paired-end reads
is_inf1000.dt[, `:=` ("R1.start" = pos, "R1.end" = pos, "R2.start" = mpos, "R2.end" = mpos)]
setkey(is_inf1000.dt, seqname, R1.start, R1.end)

setkey(annotation, seqname, R1_site_start, R1_site_end)

A <- foverlaps(is_inf1000.dt, annotation[, c("seqname", "R1_site_start", "R1_site_end", "repeat_id", "Note", "Number")], 
               by.x = c("seqname", "R1.start", "R1.end"),
               by.y = c("seqname", "R1_site_start", "R1_site_end"))
setnames(A, c("Note", "Number"), c("R1.Note", "R1.Number"))
setkey(annotation, "repeat_id", "R2_site_start", "R2_site_end")
setkey(A, repeat_id, R2.start, R2.end)

# Add annotation now for read2
B <- foverlaps(A, annotation[, c("repeat_id", "R2_site_start", "R2_site_end", "Note", "Number")],
               by.x = c("repeat_id", "R2.start", "R2.end"),
               by.y = c("repeat_id", "R2_site_start", "R2_site_end"))
setnames(B, c("repeat_id", "Note", "Number"), c("read2.repeat_id", "read2.Note", 'read2.Number'))

# Function to extract the paired reads position in a table:  
# readname, seqnames, read1_position, read2_position and insert_size
readBamRegion <- function(bamPath, seqnames, start, end, min_isize){
  bamFile <- BamFile(bamPath)
  cli_alert_info("Loading {.file {bamPath}}.")
  gr <- GRanges(seqnames, ranges = IRanges(start, end))
  #params <- ScanBamParam(which = gr, what = scanBamWhat()) # To get the entire bam info. Needs more compute power and ram!
  params <- ScanBamParam(which = gr, what = c("qname", "pos", "mpos", "isize"))
  aln <- scanBam(bamFile, param = params)
  region_id <- paste0(seqnames, ":", start, "-", end)
  dt <- data.table(rname = aln[[region_id]]$qname, read1 = aln[[region_id]]$pos, read2 = aln[[region_id]]$mpos, isize = aln[[region_id]]$isize, seqnames = seqnames)
  dt[, `:=` ("read1.start" = read1, 'read1.end' = read1)]
  dt[, `:=` ("read2.start" = read2, 'read2.end' = read2)]
  setkey(dt, seqnames, read1.start, read1.end)
  dt[isize >= min_isize]
}

