library(rtracklayer)
library(data.table)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(pander)
library(knitr)
library(cli)

source("utils.R")

# Samples sheet
names <- fread("../ReferenceMaterials/Annotations/Samples.csv", header = T)

# Annotation file of repeats of interest 
annotation <- fread("../ReferenceMaterials/Annotations/Repeats_links.csv")
setkey(annotation, seqname, window.start, window.end)

# GFF of repeats 
repeats_annotations <- as.data.table(import("../ReferenceMaterials/Annotations/Repeats_region.gff3"))
repeats_annotations[, Note := tstrsplit(Note, "-", keep = 1)]
others_repeats_annotations <- repeats_annotations[!Note %in% annotation$Note]

# Annotation to draw in plot
# CutSite --> 6-2 & UP
cutsites <- data.table(seqname = c("NC_037304.1", "NC_037304.1"), 
                       start = c(265238, 266820), end = c( 265292, 266880), 
                       Repeat = c("6-2", "UP"))

setnames(annotation, "Note", "Repeat", skip_absent = TRUE)
setnames(others_repeats_annotations, c("Note", "seqnames"), c("Repeat", "seqname"), skip_absent = TRUE)

# DATA TABLE USED IN GGPLOT 
annotations_to_show <- rbindlist(list(cutsites,
                                      annotation[Repeat %in% c("V", "Q", "L", "G"), c("seqname", "start", "end", "Repeat")],
                                      others_repeats_annotations[Repeat %like% "Large", c("seqname", "start", "end", "Repeat")]))

# Coverage files
full_filtered_paths <- list.files("../1_Mapping_bwa_full/bwa_normal/microhomologies", "*.full_filtered_reads.coverage.txt", full.names = T)
names(full_filtered_paths) <- basename(gsub(".ChrM.SCLEN20.RLEN20.full_filtered_reads.coverage.txt", "", full_filtered_paths))

softclips_paths <- list.files("~/Projects/Rokas_final/1_Mapping_bwa_full/bwa_normal/microhomologies", "*.softclips.no_small_repeats.coverage.txt", full.names = T)
names(softclips_paths) <- basename(gsub(".ChrM.SCLEN20.RLEN20.softclips.no_small_repeats.coverage.txt", "", softclips_paths))

full_filtered.dt <- rbindlist(lapply(full_filtered_paths, fread, col.names = c("Chr", "Start", "Coverage")), id = "sample")
softclips_paths.dt <- rbindlist(lapply(softclips_paths, fread, col.names = c("Chr", "Start", "Coverage")), id = "sample")

# Total mapping stats
bam_coverage <- fread("../Results/Coverage.stats.csv")

# Make a named vector as dict
coverage_dict <- bam_coverage$numreads
names(coverage_dict) <- bam_coverage$SampleName

# Compute CPM with the dict Total read mapped per sample
full_filtered.dt[, full_filtered_CPM := (Coverage / coverage_dict[sample]) * 1e6, by = sample]
softclips_paths.dt[, softclip_CPM := (Coverage / coverage_dict[sample]) * 1e6, by = sample]


# Merge FULL and SOFT
coverage.dt <- merge(full_filtered.dt, softclips_paths.dt, by = c("sample", "Chr", "Start"))
coverage.dt[, sum := full_filtered_CPM + softclip_CPM]

# Remove 150bp at start & end region (HIGH coverage zone)
coverage.filt.dt <- coverage.dt[Start >= 150 & Start < 367658 ][!(Start >= 194197 & Start <= 194347)]
coverage.filt.dt[sample == "Col-0_TAL_6-2_S58"]

fwrite(coverage.filt.dt, file = "../Results/coverage.microhomology.txt.gz", nThread = 4)

# If already computed
coverage.filt.dt <- fread("../Results/coverage.microhomology.txt.gz")

ggplot(coverage.filt.dt[sample == "Col-0_TAL_6-2_S58"], aes(x = Start, y = softclip_CPM)) + geom_point(size = 0.1, alpha = 0.8) + theme_bw()

#coverage.filt.dt[sample == "NGS257-1E12_S4" & softclip_CPM >= 5]$softclip_CPM

for(Sample in unique(coverage.filt.dt$sample)) {
  print(Sample)
  ggplot(coverage.filt.dt[Name == Sample], aes(x = Start, y = softclip_CPM)) + geom_point(size = 0.1) + theme_bw()
  ggsave(filename = paste0(Sample, "_softclip_CPM.png"), path = "../Results/Plots/Softclip_plots", width = 12)
}

# Peak calling with IDPmisc -----------------------------------------------
# IDPmisc is the best option right now 
# The options of the peaks() function are:
# minPH:	Mimimum height of peak to be reported. We are using 3*(mean)
# minPW: Minimum width of peak at half maximum to be reported.
# thr: Threshold below which the signal is not processed.
library(IDPmisc)

# Get peaks from a distribution
getPeaks <- function(Sample, dt, minPW = 20, thr = 4, ...) {
  cli_alert_info("Analysing {Sample}")
  subdt <- dt[sample == Sample]
  sumdt <- summary(subdt[softclip_CPM > 0]$softclip_CPM)
  pander(sumdt)
  tryCatch({
    pts <- IDPmisc::peaks(subdt$Start, subdt$softclip_CPM, minPH = sumdt["Mean"]*3, minPW, thr, ...)
    cli_alert_success("{nrow(pts)} peaks detected!")
    return(pts)
  }, error = function(e) {
    cli_alert_danger("No peak detected!")
    return(NA)
  })
}

getPeaks(Sample = "Col-0_TAL_6-2_S1", dt = coverage.filt.dt)
getPeaks(Sample = "Col-0_TAL_6-2_S58", dt = coverage.filt.dt)

unique(coverage.dt$sample)
unique(coverage.filt.dt$sample)

# LOOP through all samples
SAMPLES <- unique(coverage.filt.dt$sample)
names(SAMPLES) <- unique(coverage.filt.dt$sample)

res <- lapply(SAMPLES, getPeaks, dt = coverage.filt.dt)
res[is.na(res)]
res.dt <- rbindlist(res[!is.na(res)], idcol = "Sample")
setnames(res.dt, c("x", "y", "w"), c("Position", "Height", "Width"))
fwrite(res.dt, "../Results/Peaks_detected.csv")

# LOAD back once computed
PEAKS <- fread("../Results/Peaks_detected.csv")
unique(PEAKS$Sample)

# Merge annotation to show with the coverage dt
setkey(annotations_to_show, seqname, start, end)
coverage.filt.dt[, End := Start]
setkey(coverage.filt.dt, Chr, Start, End)

coverage.filt.dt <- foverlaps(coverage.filt.dt, annotations_to_show, 
          by.x = c("Chr", "Start", "End"), 
          by.y = c("seqname", "start", "end"))


drawPeaks <- function(dt, peaks = PEAKS, zoom = F, ggarr = F, ...) {
  sample <- unique(dt$sample) 
  cli_alert_info("Analysing {sample}")
  sumdt <- summary(dt[softclip_CPM > 0]$softclip_CPM)
  pander(sumdt)
  pts <- peaks[Sample %in% sample]
  print(kable(pts, caption = "Peaks detected"))
  
  a <- ggplot(dt, aes(x = Start, y = softclip_CPM, color = Repeat)) + 
    geom_point(size = 0.2) + 
    theme_bw() + # geom_line() +
      scale_x_continuous(labels = format_si(), breaks = seq(0, 370000, 20000), 
                         minor_breaks = seq(0, 370000, 10000)) +

      scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50), 
                         labels = c("0", "10", "20", "30", "40", "50"), 
                                    limits = c(-5, 50)) +     
      #geom_hline(aes(yintercept = sumdt["Mean"], linetype = "Mean coverage"), col="steelblue") + 
      ylab("Softclip coverage (in CPM)") + xlab("NC_037304.1 genomic position") +
      geom_rect(data = annotations_to_show, aes(xmin = start, xmax = end, ymin = -1, ymax = -4, fill = Repeat), inherit.aes = F) 
      #scale_linetype_manual(name = "", values = c(2, 2), guide = guide_legend(override.aes = list(color = "steelblue")))
   
  b <- ggplot(dt, aes(x = Start, y = softclip_CPM, color = Repeat)) + 
    geom_point(size = 0.2) + 
    theme_bw() + #geom_col() +
    scale_x_continuous(labels = format_si(), breaks = seq(0, 370000, 20000), 
                       minor_breaks = seq(0, 370000, 10000)) +
    scale_y_continuous(breaks = c(50, 100, 150), 
                       labels = c("50", "100", "150"), 
                       limits = c(50, 150)) +
    ylab("Softclip coverage (in CPM)") + xlab("NC_037304.1 genomic position") +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

  if(sample %in% pts$Sample){
      a <- a + geom_point(data = pts, aes(x = Position, y = Height ), color = "red")
  }
  
  if(zoom){
      a <- a + scale_x_continuous(limits = c(250000, 300000), breaks = seq(250000, 300000, 5000), labels = format_si())
      b <- b + scale_x_continuous(limits = c(250000, 300000), breaks = seq(250000, 300000, 5000), labels = format_si())
  }
  c = b / a + plot_layout(guides = 'collect', heights = c(0.5, 0.8)) + labs(x = "NC_037304.1 genomic position", y = "Softclip coverage (in CPM)")
  if(ggarr) {
    c <- ggarrange(b + rremove("ylab"), a + rremove("ylab"), ncol = 1, common.legend = T, legend = "right", align = "v", heights = c(1, 2)) 
    d <- annotate_figure(c, left = text_grob("Softclip coverage", rot = 90, vjust = 1))
    return(d)
  }
  return(c)
}


coverage.filt.dt[sample == "why2-1_TAL_6-2_S1" & softclip_CPM > 0]

drawPeaks(coverage.filt.dt[sample == "why2-1_TAL_6-2_S1" & softclip_CPM > 0], ggarr = T)
drawPeaks(coverage.filt.dt[sample == "NGS354_Why2-1_BM1_NA_D1_S5" & softclip_CPM > 0], zoom = F, ggarr = F)
drawPeaks(coverage.filt.dt[sample == "why2-1_TAL_UP_S5" & softclip_CPM > 0], ggarr = T)

for(Sample in SAMPLES) {
  message(Sample)
  tryCatch({
    A <- drawPeaks(coverage.filt.dt[softclip_CPM >= 1 & sample == Sample], ggarr = T)
    print(A)
    ggsave(plot = A, filename = paste0(Sample, ".softclip_CPM_FULL.pdf"), 
           path = "../Results/Plots/Softclip_coverage_with_peaks", width = 12)
  }, error = function(e) print(e))
  #B <- drawPeaks(coverage.filt.dt[softclip_CPM >= 1 & sample == Sample], ggplot = T, zoom = T)
  #print(B)
  #ggsave(plot = B, filename = paste0(Sample, ".softclip_CPM_ZOOM.pdf"), path = "../Results/Plots/Softclip_coverage_with_peaks", width = 12)
}

for(Sample in SAMPLES) {
  message(Sample)
  tryCatch({
    A <- drawPeaks(coverage.filt.dt[softclip_CPM >= 1 & sample == Sample], ggarr = T, zoom = T)
    ggsave(plot = A, filename = paste0(Sample, ".softclip_CPM_ZOOM.pdf"), 
           path = "../Results/Plots/Softclip_coverage_with_peaks", width = 12)
  }, error = function(e) print(e))
 }






