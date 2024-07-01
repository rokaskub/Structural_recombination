library(data.table)
library(ggplot2)

source("utils.R")

# Load samples info
names <- fread("../ReferenceMaterials/Annotations/Samples.csv", header = T)
names <- names[Genotype_construct != "_"]

# Repeat annotation
repeats_annotations <- as.data.table(import("../ReferenceMaterials/Annotations/Repeats_region.gff3"))
repeats_annotations[, Note := tstrsplit(Note, "-", keep = 1)]

# Deeptools coverage files
covfiles <- list.files("../byGenotype/1_Mapping_bwa_full/bwa_normal/insert_size", "*.MD.IS1000.filtered.bam.deeptools.1bp_bins.cpm.txt", full.names = T)
names(covfiles) <- gsub(".MD.IS1000.filtered.bam.deeptools.1bp_bins.cpm.txt", "", basename(covfiles))
coverage.dt <- rbindlist(lapply(covfiles, fread, col.names = c("seqname", "start", "end", "coverage")), idcol = "Name")

# TODO add SSA: 258128 -> 0 
SSA.coverage.dt <- coverage.dt[Name %like% 'ChrM_258128_270303']
coverage.dt <- coverage.dt[!Name %like% 'ChrM_258128_270303']

merge(coverage.dt, names[, c("Genotype", "Construct", "Genotype_construct")])

coverage.dt[, max(coverage)]

setkey(coverage.dt, seqname, start, end)
setkey(repeats_annotations, seqnames, start, end)

coverage.dt <- foverlaps(coverage.dt, repeats_annotations)
unique(coverage.dt$Name)

ggplot(coverage.dt[Name == "msh1-1_TAL-2.ChrM"], aes(x = i.start, y = coverage, color = Note)) + geom_point(size = 0.8) +
  theme_bw() + xlab("Genomic position") + ylab("Coverage (in CPM)") + #+ scale_color_discrete(na.value="black")
  scale_x_continuous(labels = format_si())

ggplot(coverage.dt[Name =="why2-1_TAL-1.ChrM" ], aes(x = i.start, y = coverage, fill = Note, color = Note)) + geom_col() +
  theme_bw() + xlab("Genomic position") + ylab("Coverage (in CPM)") +
  scale_x_continuous(labels = format_si(), breaks = seq(from = 0, 367808, by = 20000))
           
unique(coverage.dt$Name)

# 
lapply(unique(coverage.dt$Name), function(sample){
  message(sample)
  A <- ggplot(coverage.dt[Name == sample], aes(x = i.start, y = coverage, fill = Note, color = Note)) + geom_col() +
  theme_bw() + xlab("Genomic position") + ylab("Coverage (in CPM)") + scale_x_continuous(labels = format_si(), breaks = seq(from = 0, 367808, by = 20000)) + 
  ggtitle(sample)
  ggsave(plot = A, filename = paste0("../Results_genotype/byGenotype/Coverage_plots_byGenotypeConstruct/", sample, ".pdf"), width = 14, height = 8)
})

# Plot links
readBamRegion <- function(bamPath, seqnames, start, end, min_isize){
  bamFile <- BamFile(bamPath)
  cli_alert_info("Loading {.file {bamPath}}.")
  gr <- GRanges(seqnames, ranges = IRanges(start, end))
  params <- ScanBamParam(which = gr, what = c("qname", "pos", "mpos", "isize"))
  aln <- scanBam(bamFile, param = params)
  region_id <- paste0(seqnames, ":", start, "-", end)
  dt <- data.table(rname = aln[[region_id]]$qname, read1 = aln[[region_id]]$pos, read2 = aln[[region_id]]$mpos, isize = aln[[region_id]]$isize, seqnames = seqnames)
  dt[, `:=` ("read1.start" = read1, 'read1.end' = read1)]
  dt[, `:=` ("read2.start" = read2, 'read2.end' = read2)]
  setkey(dt, seqnames, read1.start, read1.end)
  #dt[abs(isize) >= min_isize]
  dt[isize >= min_isize]
}

# Distribution analysis of reads with insert size <= 1000
is1000 <- list.files("../1_Mapping_bwa_full/bwa_normal", "*ChrM.MD.bam$", full.names = T)
names(is1000) <- gsub(".ChrM.MD.bam", "", basename(is1000))

is1000.dt <- rbindlist(lapply(is1000[1], function(bam){
  count.dt <- readBamRegion(bamPath = bam, seqnames = "NC_037304.1", start = 200, end = 367608, min_isize = 0)
}), idcol = "Sample")

ggraph(is1000.dt[read1 >= 600 & read2 <= 367000], layout = 'linear', circular = TRUE) + 
  geom_edge_arc(alpha = 0.1) + coord_fixed()
