library(data.table)
library(ggplot2)
library(ggraph)
library(rtracklayer)

source("utils.R")

# Load samples info
names <- fread("../ReferenceMaterials/Annotations/Samples.csv", header = T)
names <- names[Genotype_construct != "_"]

# Repeat annotation
repeats_annotations <- as.data.table(import("../ReferenceMaterials/Annotations/Repeats_region.gff3"))
repeats_annotations[, Note := tstrsplit(Note, "-", keep = 1)]

# Load bwa normal coverage per bin
#covfiles <- list.files("~/Projects/Rokas_final/byGenotype/1_Mapping_bwa_full/bwa_normal", "*.ChrM.MD.bam.deeptools.100bp_bins.cpm.txt", full.names = T)

# Solo bam
covfiles <- list.files("../1_Mapping_bwa_full/bwa_normal/insert_size", "*.MD.IS1000.filtered.bam.deeptools.100bp_bins.cpm.txt", full.names = T)
names(covfiles) <- gsub(".MD.IS1000.filtered.bam.deeptools.100bp_bins.cpm.txt", "", basename(covfiles))

# Bam merged ByGenotype 
# covfiles <- list.files("~/Projects/Rokas_final/byGenotype/1_Mapping_bwa_full/bwa_normal/insert_size", "*.MD.IS1000.filtered.bam.deeptools.1bp_bins.cpm.txt", full.names = T)
# names(covfiles) <- gsub(".MD.IS1000.filtered.bam.deeptools.100bp_bins.cpm.txt", "", basename(covfiles))

coverage.dt <- rbindlist(lapply(covfiles, fread, col.names = c("seqname", "start", "end", "coverage")), idcol = "Name")
coverage.dt <- coverage.dt[!Name %like% 'ChrM_258128_270303']

coverage.dt[, max(coverage)]

setkey(coverage.dt, seqname, start, end)
setkey(repeats_annotations, seqnames, start, end)

coverage.dt <- foverlaps(coverage.dt, repeats_annotations)
unique(coverage.dt$Name)

lapply(unique(coverage.dt$Name), function(sample){
  message(sample)
  A <- ggplot(coverage.dt[Name == sample], aes(x = i.start, y = coverage, color = Note)) + geom_point(size = 0.8) +
    theme_bw() + xlab("Genomic position") + ylab("Coverage (in CPM)") + scale_x_continuous(labels = format_si()) + ggtitle(sample)
  #ggsave(plot = A, filename = paste0("../Results/byGenotype/Coverage_plots/", sample, ".pdf"), width = 14, height = 8)
})

ggplot(coverage.dt[Name == "WT_control.ChrM"], aes(x = i.start, y = coverage, color = Note)) + geom_point(size = 0.8) +
  theme_bw() + xlab("Genomic position") + ylab("Coverage (in CPM)") + #+ scale_color_discrete(na.value="black")
  scale_x_continuous(labels = format_si())
