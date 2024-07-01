library(data.table)
library(Rsamtools)
library(circlize)
library(ggplot2)
library(viridis)
library(ggraph)
library(igraph)
library(plyr)
library(cli)

source("utils.R")

# SSA Annotation file
annotation <- fread("../ReferenceMaterials/Annotations/SSA_window_annotation_updated.csv")
table(annotation$Note)
annotation[, repeat_id := paste(Note, Number, sep = "_")]
annotation[, fragment_size := (end - start) + 1] 
annotation <- annotation[!is.na(start)]
setkey(annotation, seqname, window.start, window.end)

# Load samples information
names <- fread("../ReferenceMaterials/Annotations/Samples.csv", header = T)
names$Genotype_construct

# MD Coverage files for CPM normalization  ---------------------------------------------------
norm <- list.files("../byGenotype/1_Mapping_bwa_full/bwa_normal", "*.ChrM.MD.bam.coverage.txt$", full.names = T)
names(norm) <- gsub(".ChrM.MD.bam.coverage.txt", "", basename(norm))
norm.dt <- rbindlist(lapply(norm, fread), idcol = "Sample")
#norm.dt[, factor := round(numreads / 2932880, 3)] # Based on Col0

readCount <- function(bamPath, annotation){
  cli_alert_info("Extracting reads from {.file {basename(bamPath)}}")
  gr <- makeGRangesFromDataFrame(annotation[, c("seqname", "Note", "R1_site_start", "R1_site_end", "repeat_id")], 
                                 start.field = "R1_site_start", end.field = "R1_site_end", strand.field = "strand", 
                                 keep.extra.columns = T)
  bamFile <- BamFile(bamPath)
  params <- ScanBamParam(which = gr, what = c("qname", "pos", "mpos", "isize"))
  aln <- scanBam(bamFile, param = params)
  names(aln) <- annotation$repeat_id
  res <- rbindlist(aln, idcol = "repeat_id", use.names = T)
  res[, seqname := "NC_037304.1"]
  res <- res[isize >= 0]
  return(res)
}

# Insert size >=1000bp -----------------------------------------------------
bam_files <- list.files("../byGenotype/1_Mapping_bwa_full/bwa_strict/insert_size", "*.ChrM_258128_270303.MD.IS1000.filtered.bam$", full = T)
names(bam_files) <- gsub(".ChrM_258128_270303.MD.IS1000.filtered.bam", "", basename(bam_files))

# Go throught all the files
filtered_bam.dt <- rbindlist(lapply(bam_files, function(bam){
  count.dt <- readBamRegion(bamPath = bam, min_isize = 0)
}), idcol = "Sample")

# Change NC_037304.1:start-end in --> NC_037304.1 to be able to merge in the end
filtered_bam.dt[, seqname := tstrsplit(seqname, ":", keep = 1)]
fwrite(filtered_bam.dt, "../Results_genotype/SSA_PE_reads_IS_sup1000.raw_reads.csv.gz", sep = ",")

filtered_bam.dt <- fread("../Results_genotype/SSA_PE_reads_IS_sup1000.raw_reads.csv.gz")
setkey(filtered_bam.dt, seqname, read1.start, read1.end)

A <- foverlaps(filtered_bam.dt, annotation[, c("seqname", "window.start", "window.end", "repeat_id", "Note", "Number")], 
               by.x = c("seqname", "read1.start", "read1.end"),
               by.y = c("seqname", "window.start", "window.end"))
setnames(A, c("repeat_id", "Note", "Number"), c("read1.repeat_id", "read1.Note", "read1.Number"))

# Add annotation now for read2
B <- foverlaps(A, annotation[, c("seqname", "window.start", "window.end", "repeat_id", "Note", "Number")],
               by.x = c("seqname", "read2.start", "read2.end"), 
               by.y = c("seqname", "window.start", "window.end"))
setnames(B, c("repeat_id", "Note", "Number"), c("read2.repeat_id", "read2.Note", 'read2.Number'))

# Add pairwise name to filter them easily
B[, Pairwise := paste(read1.repeat_id, read2.repeat_id, sep = "-")]
B[! Pairwise %like% "NA" & read1.Note == read2.Note, .N, by = c("Sample", "Pairwise", "read1.Note")][order(N)]
#B <- B[!(read1.Note == read2.Note & read1.Number == read2.Number) & isize > 0]

B[read1.repeat_id ==  read2.repeat_id]

B[, Copy := "others"]

# Add Copy information. SSA only has 2 and 3
B[read1.Note == read2.Note & read1.Number == 2 & read2.Number == 3, Copy := "CP2_CP3"]
B[read1.Note == read2.Note & read1.Number == 3 & read2.Number == 2, Copy := "CP2_CP3"]

table(B$Pairwise)

# Save OTHERS
others <- B[Copy == "others"]
#fwrite(others, "~/Project/Rokas/Results/SSA_PE_reads_IS_sup1000.others_count.csv")

# Remove others for pairwise analysis
pairwise_B <- B[Copy != "others"]
# Check
length(B$qname)
length(others$qname) + length(pairwise_B$qname)

pairwise_B[, .N, by = c("Sample", "Pairwise")][N > 1]

# Total number
C <- pairwise_B[, .N, by = c("Sample",  "read1.repeat_id", "read2.repeat_id", "Pairwise", "Copy")]

D.tmp <- merge(C, annotation[, c("seqname", "start", "end", "repeat_id", "Note")], by.x = "read1.repeat_id", by.y = "repeat_id")
setnames(D.tmp, c("start", "end"), c("CP1.start", "CP1.end"))

E <- merge(D.tmp, annotation[, c("start", "end", "repeat_id")], by.x = "read2.repeat_id", by.y = "repeat_id")
setnames(E, c("start", "end"), c("CP2.start", "CP2.end"))

# Add CPM using the norm.dt
summary.dt <- merge.data.table(E, norm.dt[, c("Sample", "numreads")], by = "Sample", all.x = T)
summary.dt[, CPM := round((N / numreads) * 1e6, 3)]
summary.dt[, .N, by = Copy]

fwrite(summary.dt, "../Results_genotype/SSA_PE_reads_IS_sup1000.count.csv.gz", sep = ",")
summary.dt <- fread("../Results_genotype/SSA_PE_reads_IS_sup1000.count.csv.gz")

#G <- setDT(tidyr::pivot_wider(summary.dt, names_from = "Copy", values_from = "CPM"))

# Add repeats connection_type
#summary.dt <- merge(summary.dt, repeats_connection[, c("Note", "Copy", "Type")], by = c("Note", "Copy"), all.x = TRUE)

# Save single file per sample
for(sample in unique(E$Sample)){
  subset = E[Sample == sample]
  print(sample)
  fwrite(subset, paste0("../Results_genotype/", sample, ".count_summary.csv"), sep = ",")
}
