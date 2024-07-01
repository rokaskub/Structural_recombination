library(data.table)
library(Rsamtools)
library(circlize)
library(ggplot2)
library(viridis)
library(ggraph)
library(igraph)
library(plyr)
library(cli)

# Load repeats information
annotation <- fread("../ReferenceMaterials/Annotations/Repeats_links.csv")
setkey(annotation, seqname, window.start, window.end)

# Load samples information
names <- fread("../ReferenceMaterials/Annotations/Samples.csv", header = T)

# Load repeats connection type (directed/inverted)
repeats_connection <- fread("../ReferenceMaterials/Annotations/Repeats_links_type.csv")
repeats_connection <- tidyr::pivot_longer(repeats_connection, cols = starts_with("CP"), names_to = "Copy", values_to = "Type")

# Function to extract the paired reads position in a table:  
# readname, seqnames, read1_position, read2_position and insert_size
readBamRegion <- function(bamPath, seqnames, start, end, min_isize){
  bamFile <- BamFile(bamPath)
  cli_alert_info("Loading {.file {bamPath}}.")
  gr <- GRanges(seqnames, ranges = IRanges(start, end))
  #params <- ScanBamParam(which = gr, what = scanBamWhat()) # To get the entire bam info, but needs more compute power and ram!
  params <- ScanBamParam(which = gr, what = c("qname", "pos", "mpos", "isize")) #, flag = scanBamFlag(isFirstMateRead=TRUE))
  aln <- scanBam(bamFile, param = params)
  region_id <- paste0(seqnames, ":", start, "-", end)
  dt <- data.table(rname = aln[[region_id]]$qname, read1 = aln[[region_id]]$pos, read2 = aln[[region_id]]$mpos, isize = aln[[region_id]]$isize, seqnames = seqnames)
  dt[, `:=` ("read1.start" = read1, 'read1.end' = read1)]
  dt[, `:=` ("read2.start" = read2, 'read2.end' = read2)]
  setkey(dt, seqnames, read1.start, read1.end)
  dt[isize >= min_isize]
}

# Bam coverage files for normalization  ---------------------------------------------------
norm <- list.files("../byGenotype/1_Mapping_bwa_full/bwa_normal", "*.ChrM.MD.bam.coverage.txt$", full.names = T)
names(norm) <- gsub(".ChrM.MD.bam.coverage.txt", "", basename(norm))
norm.dt <- rbindlist(lapply(norm, fread), idcol = "Name")
fwrite(norm.dt, "../Results_genotype/Coverage.stats.csv", sep = ",")

# Microhomologies ---------------------------------------------------------
# Get counts
norm.hom <- list.files("../byGenotype/1_Mapping_bwa_full/bwa_normal/microhomologies", "*ChrM.MD.MICROHOMO.NM2_SCLEN20.bowtie2.bam$", full.names = T)
names(norm.hom) <- gsub(".ChrM.MD.MICROHOMO.NM2_SCLEN20.bowtie2.bam", "", basename(norm.hom))

norm.hom.dt <- rbindlist(lapply(norm.hom, function(bam){
  count.dt <- readBamRegion(bamPath = bam, seqnames = "NC_037304.1", start = 200, end = 367608, min_isize = 0)
}), idcol = "Name")

# Filtering out this region 194197-198389 ??
norm.hom.dt <- norm.hom.dt[!(read1 >= 194197 & read1 <= 198389)]
norm.hom.dt

# Get count summary
norm.hom.dt.count <- norm.hom.dt[, .N, by = Name]
norm.hom.dt.count <- merge.data.table(norm.hom.dt.count, norm.dt[, c("Name", "numreads")], by = "Name", all.x = T)
norm.hom.dt.count[, CPM := round((N/numreads)*1e6, 3)]

ggplot(norm.hom.dt.count, aes(Name, CPM)) + geom_col() + theme_bw() + coord_flip()

# Specific region 
# 264156 - 268528
cut_site <- norm.hom.dt[(read1 >= 264156 & read1 <= 268528) | (read1 >= 265000 & read1 <= 265600) , .N, by = Name]
cut_site <- norm.hom.dt[(read1 >= 266500 & read1 <= 267200) | (read1 >= 265000 & read1 <= 265600) , .N, by = Name]
cut_site <- cut_site[, .N, by = Name]
cut_site <- merge.data.table(cut_site, norm.dt[, c("Name", "numreads")], by = "Name", all.x = T)
cut_site[, CPM := round((N/numreads)*1e6, 3)]

ggplot(cut_site, aes(Name, CPM)) + geom_col() + theme_bw() + coord_flip() 

cut_site[Name %like% "msh"]
norm.hom.dt.count[Name %like% "msh"]

# Microhomologies
# Do not count reads in these windows
repeats_annotations <- import("../ReferenceMaterials/Annotations/Repeats_region.gff3")
repeats_annotations <- as.data.table(promoters(repeats_annotations, upstream = 100, downstream= 100)) # Add 100bp down/upstream
repeats_annotations[, Note := tstrsplit(Note, "-", keep = 1)]
others_repeats_annotations <- repeats_annotations[!Note %in% annotation$Note]

setkey(repeats_annotations, seqnames, start, end) 
setnames(norm.hom.dt, c("read1.start", "read1.end"), c("start", "end"))
setkey(norm.hom.dt, seqnames, start, end) 

test <- foverlaps(norm.hom.dt, repeats_annotations[, c("seqnames", "start", "end", "width", "strand", "Note")])
test[is.na(Note), .N, by = Name]

# MD insertsize filtered bam ----------------------------------------------
# Extract annotations for read1 and read2. Count by pairwise annotations. 

bam_files <- list.files("../byGenotype/1_Mapping_bwa_full/bwa_strict/insert_size", "*.ChrM.MD.IS1000.filtered.bam$", full = T)
names(bam_files) <- gsub(".ChrM.MD.IS1000.filtered.bam", "", basename(bam_files))

filtered_bam.dt <- rbindlist(lapply(bam_files, function(bam){
  count.dt <- readBamRegion(bamPath = bam, seqnames = "NC_037304.1", start = 200, end = 367608, min_isize = 0)
}), idcol = "Name")

fwrite(filtered_bam.dt, "../Results_genotype/Repeats_PE_reads_IS_sup1000.raw_reads.csv.gz", nThread = 6)
filtered_bam.dt <- fread("../Results_genotype/Repeats_PE_reads_IS_sup1000.raw_reads.csv.gz")

mean_isize <- filtered_bam.dt[, mean(isize), by = Name]
setkey(mean_isize, Name)
fwrite(mean_isize, "../Results_genotype/Repeats_PE_reads_IS_sup1000.mean_isize.csv", nThread = 6)

setkey(filtered_bam.dt, seqnames, read1.start, read1.end)

A <- foverlaps(filtered_bam.dt, annotation[, c("seqname", "window.start", "window.end", "repeat_id", "Note", "Number")], 
               by.x = c("seqnames", "read1.start", "read1.end"),
               by.y = c("seqname", "window.start", "window.end"))
setnames(A, c("repeat_id", "Note", "Number"), c("read1.repeat_id", "read1.Note", "read1.Number"))

# Add annotation now for read2
B <- foverlaps(A, annotation[, c("seqname", "window.start", "window.end", "repeat_id", "Note", "Number")],
               by.x = c("seqnames", "read2.start", "read2.end"), 
               by.y = c("seqname", "window.start", "window.end"))
setnames(B, c("repeat_id", "Note", "Number"), c("read2.repeat_id", "read2.Note", 'read2.Number'))

# Add pairwise name to filter them easily
B[, Pairwise := paste(read1.repeat_id, read2.repeat_id, sep = "-")]
B[! Pairwise %like% "NA" & read1.Note == read2.Note, .N, by = c("Name", "Pairwise", "read1.Note")]

B <- B[!(read1.Note == read2.Note & read1.Number == read2.Number) & isize > 0]

# Add Copy information (CP1-5)
B[, Copy := "others"]
B[(read1.Note == read2.Note & read1.Number == read2.Number), Copy := "others"]

B[read1.Note == read2.Note & read1.Number == 1 & read2.Number == 2, Copy := "CP1_CP2"]
B[read1.Note == read2.Note & read1.Number == 2 & read2.Number == 1, Copy := "CP1_CP2"]
B[read1.Note == read2.Note & read1.Number == 1 & read2.Number == 3, Copy := "CP1_CP3"]
B[read1.Note == read2.Note & read1.Number == 3 & read2.Number == 1, Copy := "CP1_CP3"]
B[read1.Note == read2.Note & read1.Number == 1 & read2.Number == 4, Copy := "CP1_CP4"]
B[read1.Note == read2.Note & read1.Number == 4 & read2.Number == 1, Copy := "CP1_CP4"]
B[read1.Note == read2.Note & read1.Number == 1 & read2.Number == 5, Copy := "CP1_CP5"]
B[read1.Note == read2.Note & read1.Number == 5 & read2.Number == 1, Copy := "CP1_CP5"]

B[read1.Note == read2.Note & read1.Number == 2 & read2.Number == 3, Copy := "CP2_CP3"]
B[read1.Note == read2.Note & read1.Number == 3 & read2.Number == 2, Copy := "CP2_CP3"]
B[read1.Note == read2.Note & read1.Number == 2 & read2.Number == 4, Copy := "CP2_CP4"]
B[read1.Note == read2.Note & read1.Number == 4 & read2.Number == 2, Copy := "CP2_CP4"]
B[read1.Note == read2.Note & read1.Number == 2 & read2.Number == 5, Copy := "CP2_CP5"]
B[read1.Note == read2.Note & read1.Number == 5 & read2.Number == 2, Copy := "CP2_CP5"]

B[read1.Note == read2.Note & read1.Number == 3 & read2.Number == 4, Copy := "CP3_CP4"]
B[read1.Note == read2.Note & read1.Number == 4 & read2.Number == 3, Copy := "CP3_CP4"]
B[read1.Note == read2.Note & read1.Number == 3 & read2.Number == 5, Copy := "CP3_CP5"]
B[read1.Note == read2.Note & read1.Number == 5 & read2.Number == 3, Copy := "CP3_CP5"]

B[read1.Note == read2.Note & read1.Number == 4 & read2.Number == 5, Copy := "CP4_CP5"]
B[read1.Note == read2.Note & read1.Number == 5 & read2.Number == 4, Copy := "CP4_CP5"]

B[, ID := paste(Name, Pairwise, sep = ".")]
B[, .N, by = Copy]

table(B[Copy == "others"]$Pairwise)
table(B$Copy)

B[, .N, by = Copy]

# Other links
others <- B[Copy == "others"]

fwrite(others, file = "../Results_genotype/Repeats_PE_reads_IS_sup1000.others_count.csv.gz", sep = ",")
#others <- fread("~/Projects/Rokas_final/Results/Repeats_PE_reads_IS_sup1000.others_count.csv.gz")

fwrite(B, "../Results_genotype/Repeats_PE_reads_IS_sup1000.raw_reads.csv.gz")

B[, .N, by = c("Name",  "read1.repeat_id", "read2.repeat_id", "Pairwise", "Copy")]

# Remove the others for pairwise analysis
pairwise_B <- B[Copy != "others"]
#pairwise_B <- B

# Check
length(B$rname)
length(others$rname) + length(pairwise_B$rname)

# Still few of them are overlapping
#pairwise_B[duplicated(pairwise_B$rname)]
#pairwise_B[rname == "VH00615:1:AACG7C7M5:1:2502:42166:46056"]
pairwise_B[, .N, by = c("Name", "Pairwise", "Copy")][N > 1]
pairwise_B[Name == "Col-0_TAL_6-2_S58"][, .N, by = "rname"][order(N)]

# Total number
C <- pairwise_B[, .N, by = c("Name",  "read1.repeat_id", "read2.repeat_id", "Pairwise", "Copy")]

D.tmp <- merge(C, annotation[, c("seqname", "start", "end", "repeat_id", "Note")], by.x = "read1.repeat_id", by.y = "repeat_id")
setnames(D.tmp, c("start", "end"), c("CP1.start", "CP1.end"))

E <- merge(D.tmp, annotation[, c("start", "end", "repeat_id")], by.x = "read2.repeat_id", by.y = "repeat_id")
setnames(E, c("start", "end"), c("CP2.start", "CP2.end"))

# Add CPM using the norm.dt
summary.dt <- merge.data.table(E, norm.dt[, c("Name", "numreads")], by = "Name", all.x = T)
summary.dt[, CPM := round((N / numreads) * 1e6, 3)]
summary.dt[, .N, by = Copy]
summary.dt[read2.repeat_id %like% "MM"]

fwrite(summary.dt, "../Results_genotype/Repeats_PE_reads_IS_sup1000.count.csv.gz", sep = ",")
summary.dt <- fread("../Results_genotype/Repeats_PE_reads_IS_sup1000.count.csv.gz")

# Save single file per sample
for(name in unique(E$Name)){
  subset = E[Name == name]
  print(name)
  fwrite(subset, paste0("../Results_genotype/", name, ".count_summary.csv"), sep = ",")
}
