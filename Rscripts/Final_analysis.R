library(rtracklayer)
library(data.table)
library(circlize)
library(ggplot2)
library(scales)
library(tidyr)
library(cli)

source("utils.R")

# Retrieve sample names
names <- fread("../ReferenceMaterials/Annotations/Samples.csv", header = T)
names <- names[Genotype != ""]

names[SampleName == "WT_control"]

# Paired-end data
res <- fread("../Results/Results_fraction_COPY.csv") 

unique(res$Note)
res[Note == "Y"]

unique(res$SampleName) # 126
res[SampleName == "WT_control"]

# Mapping stats
mapping_stats <- fread("../Results/Coverage.stats.csv")
mapping_stats <- merge(mapping_stats, names[, c("SampleName", "Name")], by = "SampleName")

# Annotation
annotation <- fread("../ReferenceMaterials/Annotations/Repeats_links.csv")

# Cut Sites
cut_sites <- data.table(seqname = c("NC_037304.1", "NC_037304.1"), start = c(265238, 266820), end = c(265292, 266880), Note = c("6-2", "UP"))

# Others repeats
repeats_annotations <- as.data.table(import("../ReferenceMaterials/Annotations/Repeats_region.gff3"))
repeats_annotations[, Note := tstrsplit(Note, "-", keep = 1)]
others_repeats_annotations <- repeats_annotations[!Note %in% annotation$Note]

# Large repeats
large_repeats <- others_repeats_annotations[Note %like% "Large"]

# Load bwa normal coverage per bin
covfiles <- list.files("../1_Mapping_bwa_full/bwa_normal", "*.ChrM.MD.bam.deeptools.100bp_bins.cpm.txt", full.names = T)
names(covfiles) <- gsub(".ChrM.MD.bam.deeptools.100bp_bins.cpm.txt", "", basename(covfiles))
coverage.dt <- rbindlist(lapply(covfiles, fread, col.names = c("seqname", "start", "end", "coverage")), idcol = "SampleName")
coverage.dt <- merge(coverage.dt, names[, c("SampleName", "Name")], by = "SampleName")
coverage.dt[, max(coverage)]

setkey(coverage.dt, seqname, start, end)
setkey(repeats_annotations, seqnames, start, end)

coverage.dt <- foverlaps(coverage.dt, repeats_annotations)
unique(coverage.dt$Name)

lapply(unique(coverage.dt$Name), function(sample){
  message(sample)
  A <- ggplot(coverage.dt[Name == sample], aes(x = i.start, y = coverage, color = Note)) + geom_point(size = 0.8) +
    theme_bw() + xlab("Genomic position") + ylab("Coverage (CPM)") + scale_x_continuous(labels = format_si()) + ggtitle(sample)
  #ggsave(plot = A, filename = paste0("../Results/Coverage_plots/", sample, ".pdf"), width = 14, height = 8)
  print(A)
})

ggplot(coverage.dt[Name == "Col-0_TAL-1_S29"], aes(x = i.start, y = coverage, color = Note)) + geom_point(size = 0.8) +
  theme_bw() + xlab("Genomic position") + ylab("Coverage (in CPM)") + #+ scale_color_discrete(na.value="black")
  scale_x_continuous(labels = format_si())

# Repeat connection file
repeats_connection <- fread("../ReferenceMaterials/Annotations/Repeats_links_type.csv")
repeats_connection <- tidyr::pivot_longer(repeats_connection, cols = starts_with("CP"), names_to = "Copy", values_to = "Type")
repeats_connection <- setDT(repeats_connection)

repeats_connection <- repeats_connection[!is.na(Type)]
repeats_connection[, tmp := gsub("CP", "", Copy)]
repeats_connection[, c("CP1", "CP2") := tstrsplit(tmp, "_")]

# Pairwise
repeats_connection[, Pairwise := paste0(Note, "_", CP1, "-", Note, "_", CP2)]

all_pairwise <- setDT(expand.grid(Pairwise = repeats_connection$Pairwise, Name = unique(res$Name)))
all_pairwise[Name == "Col-0_no-TAL1_mix"]
all_pairwise <- merge(all_pairwise, names, by = "Name")
all_pairwise <- merge(all_pairwise, repeats_connection[, c("Pairwise", "Repeat_size", "Type", "Note")], by = "Pairwise")
unique(all_pairwise$Pairwise)

res2 <- setDT(merge(all_pairwise, res[, c("Name", "Pairwise", "N_sup1000", "read2.normal", "read1.normal", "available_norm", "fraction")], 
                    by = c("Pairwise", "Name"), all.x = TRUE))

# Add fraction
res2[is.na(fraction), fraction := 0]

# Active?
res2[fraction >= 0.5, Active := "yes"]
res2[fraction < 0.5, Active := "no"]
res2[, .N, by = Active]

unique(res2$SampleName)

fwrite(res2, "../Results/Fraction_table.csv")
res2 <- fread("../Results/Fraction_table.csv")

ggplot(res2[Construct %in% c("no-TAL", "TAL-1", 'TAL-2') & Note %in% c("Y", "V", "Q", "M", "L", "G") & !Genotype %in% c("recA3_3_WT")][, .N, by = c("Genotype", "Active", "Construct", "Note")], 
       aes(x = Genotype, y = N, fill = Active)) + geom_col(position = "fill") + facet_grid(Construct ~ Note) + theme_bw() 

ggplot(res2[Construct %in% c("TAL-1", 'TAL-2') & Note %in% c("Y", "V", "Q", "M", "L", "G") & !Genotype %in% c("recA3_3_WT")][, .N, by = c("Genotype", "Active", "Construct", "Note")], 
       aes(x = Genotype, y = N, fill = Active)) + geom_col(position = "fill") + facet_grid(Construct ~ Note) + theme_bw() + 
  scale_fill_brewer(palette = "Dark2") + scale_y_continuous(labels = scales::percent)

ggplot(res2[Construct %in% c("TAL-1", 'TAL-2') & Note %in% c("Y", "V", "Q", "M", "L", "G") & !Genotype %in% c("recA3_3_WT")][, .N, by = c("Genotype_construct", "Active", "Note")], 
       aes(x = Note, y = N, fill = Active)) + geom_col(position = "fill") + facet_grid(~ Genotype_construct) + theme_bw() + 
  scale_fill_brewer(palette = "Dark2") + scale_y_continuous(labels = scales::percent)

# Loop through all the genotype
ggplot(res2[Genotype == "odb1-1" & Construct %in% c("TAL-1", "TAL-2", "no-TAL1", "no-TAL") & 
              (Note %in% c("Q", "M") | Pairwise %in% c("Y_2-Y_3", "V_1-V_3", "L_1-L_2", "G_1-G_2")) & !Genotype %in% c("recA3_3_WT") &
              !(Construct == "TAL-2" & Note %in% c("V", "Y"))], 
       aes(x = Note, y = fraction, fill = Note)) + geom_boxplot() + facet_grid(~ Construct) + theme_bw() + ggtitle("Col-0") +
  scale_fill_brewer(palette = "Dark2") + 
  scale_y_continuous(labels = scales::percent) + 
  scale_x_discrete(limits = c("Y", "V", "Q", "M", "L", "G"))
unique(res2$Genotype)
  
lapply(unique(res2$Genotype), function(genotype){
  message(genotype)
  A <- ggplot(res2[Genotype == genotype & Construct %in% c("TAL-1", "TAL-2", "no-TAL1", "no-TAL") & 
                     (Note %in% c("Q", "M") | Pairwise %in% c("Y_2-Y_3", "V_1-V_3", "L_1-L_2", "G_1-G_2")) & !Genotype %in% c("recA3_3_WT") &
                     !(Construct == "TAL-2" & Note %in% c("V", "Y"))], 
              aes(x = Note, y = fraction, fill = Note)) + geom_boxplot() + facet_grid(~ Construct) + theme_bw() + ggtitle(genotype) +
    scale_fill_brewer(palette = "Dark2") + 
    scale_y_continuous(labels = scales::percent) + 
    scale_x_discrete(limits = c("Y", "V", "Q", "M", "L", "G"))
  try(ggsave(plot = A, filename = paste0("../Results/Boxplot_repeats_byGenotype/", genotype, ".Boxplot_repeat_genotype.pdf"), width = 10, height = 6))
})
  
res2[Genotype == "Col-0" & Construct %in% c("TAL-1", 'TAL-2') & Note %in% c("Q") & !Genotype %in% c("recA3_3_WT")]

#tidyr::pivot_wider(res[, c("Sample", "Pairwise", "N_sup1000")], names_from = c("Sample"), values_from = c("N_sup1000"), values_fill = 0)
tmp <- melt(dcast(res[, c("Name", "Construct", "Genotype", "Pairwise", "N_sup1000")], 
                  Name + Construct + Genotype ~ Pairwise, value.var = 'N_sup1000', fill = 0), 
            id.vars = c("Name", "Construct", "Genotype"), variable.name = "Pairwise", value.name = "N_sup1000")
tmp1 <- merge(tmp, res[, c("Name", "Pairwise", "fraction")], by = c("Name", "Pairwise"), all.x = T)

tmp1[N_sup1000 == 0, fraction := 0]
# Active?
tmp1[fraction >= 0.5, Active := "yes"]
tmp1[fraction < 0.5, Active := "no"]

tmp1[Name == "Col-0_TAL-2_S1"]
merge(tmp1, res[, c("Name", "Note", "Type", "Pairwise")], by = c("Name", "Pairwise"), all.x = TRUE)

# Add Construct, Note, Genotype, Type, Size of repeat to tmp1
# Update the plot code
ggplot(res2[Construct %in% c("no-TAL", "TAL-1", "TAL-2") & Note %in% c("Y", "V", "Q", "M", "L", "G") & !Genotype %in% c("recA3_3_WT")][, .N, by = c("Genotype", "Active", "Construct", "Note")], 
       aes(x = Genotype, y = N, fill = Active)) + geom_col(position = "fill") + facet_grid(Construct~Note) + theme_bw() + scale_fill_brewer(palette = "Dark2")

# Circos plot -------------------------------------------------------------

# Circos fraction ---------------------------------------------------------
drawFractionCircos <- function(dt, size = 367808, min_perc = 0, min_pairs = 0, others = FALSE) {
  if(others) {
    dt <- dt[Copy == "others"]
  } else {
    dt <- dt[Copy != "others"]
  }
  sample_name <- unique(dt$Name)
  cli_alert_info(sample_name)
  cli_alert_info("Total row: {nrow(dt)}, ➚ max. reads: {max(dt$N_sup1000)} ➘ min. reads: {min(dt$N_sup1000)}")
  #cli_alert_info("Removing everything below max*{min_perc} -> {max(dt$N_sup1000)*min_perc}")
  
  # Removes everything below n% of max value & by min_pairs
  dt <- dt[N_sup1000 > max(N_sup1000) * min_perc]
  dt <- dt[N_sup1000 >= min_pairs]
  
  # Init circle track
  #circos.par("track.height" = 0.1) # 10% of space for each track
  circos.initialize(sectors = unique(dt$seqname), xlim = c(0, size)) # init the sectors
  
  # Title
  # Get mean depth for the specific sample
  meandepth <- round(mapping_stats[Name == sample_name]$meandepth)
  title(paste0(unique(dt$Name), "   (Mean depth = ", meandepth, ")"))
  
  # Annotation
  subset_annotation <- annotation[Note %in% c(dt$Note, "Y"), c("seqname", "start", "end", "Note")] 
  # TODO Remove the notes that was not present in the original 
  # repeats_annotations <- import("/home/dpflieger/Project/Rokas/ReferenceMaterials/all_repeat_region.gff3")
  # repeats_annotations <- as.data.table(promoters(repeats_annotations, upstream = 100, downstream= 100)) # Add 100bp down/upstream
  # repeats_annotations[, Note := tstrsplit(Note, "-", keep = 1)]
  # others_repeats_annotations <- repeats_annotations[!Note %in% annotation$Note]
  
  subset_annotation[, `:=` (color = "black", line_col = "black", cex = 0.7)]
  cut_sites[, `:=` (color = "red", line_col = "red", cex = 1)]
  
  # ADD the large repeats 
  subset_annotation <- rbind(subset_annotation, cut_sites)
  
  my_palette <- colorRampPalette(c("grey85", "orange", "red1"))(100)
  dt$colors <- my_palette[cut(dt$fraction, 100)]
  print(dt)
  
  # Draw label
  #circos.track(unique(dt$seqname), ylim = c(0, 1), track.height = 2)
  circos.labels(subset_annotation$seqname, x = subset_annotation$start,
                labels = subset_annotation$Note, col = subset_annotation$color,
                line_col = subset_annotation$line_col, side = "outside",
                cex = subset_annotation$cex)
  
  
  circos.track(unique(dt$seqname), ylim = c(0, 1), track.height = 0.05)
  mapply(circos.rect, subset_annotation$start, 0, subset_annotation$end, 1, col = "grey30", border = "grey30")
  mapply(circos.rect, large_repeats[Note == "Large1"]$start, 0, large_repeats[Note == "Large1"]$end, 1, col = "aquamarine1", border = "black")
  mapply(circos.rect, large_repeats[Note == "Large2"]$start, 0, large_repeats[Note == "Large2"]$end, 1, col = "aquamarine4", border = "black")
  mapply(circos.rect, cut_sites$start, 0, cut_sites$end, 1, col = "red", border = "red")
  
  subset_annotation <- annotation[Note %in% c(dt$Note, "Y")]
  
  #circos.track(ylim = c(0, 1), track.height = 0.1)
  
  # Coverage as point
  circos.track(unique(dt$seqname), ylim = c(0, 1), track.height = 0.2)
  # mapply(circos.points,
  #        coverage.dt[Sample == sample_name]$start,
  #        rescale(coverage.dt[Sample == sample_name]$coverage),
  #        col = "steelblue", pch = 20, cex = 0.05)
  
  # Coverage barplot
  mapply(circos.barplot,
         rescale(coverage.dt[Name == sample_name]$coverage,
                 from = c(0, max(coverage.dt[Name == sample_name]$coverage))),
         coverage.dt[Name == sample_name]$i.start,
         col = "steelblue", border = "steelblue")
  
  
  # Draw a line example
  #mean_rescaled = rescale(meandepth, from = c(0, max(coverage.dt[Name == sample_name]$coverage)))
  #circos.lines(x = c(0, size), y = c(mean_rescaled, mean_rescaled), col = "red")
  
  circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
    circos.genomicAxis(
      h = "bottom",
      direction = "inside",
      labels.cex=0.5,            
      labels.font=1 # labels font (1, 2, 3, 4)
    )
  })
  
  #circos.track(unique(dt$seqname), ylim = c(0, 1), track.height = 0.06)
  setorder(dt, N_sup1000)
  for(i in 1:nrow(dt)) {
    circos.link(dt[i]$seqname, c(dt[i]$CP1.start, dt[i]$CP1.end),
                dt[i]$seqname, c(dt[i]$CP2.start, dt[i]$CP2.end),
                col = dt[i]$color, lwd = 2)
  }
}

# For the NGS246* the min_pairs == 3, because they are low covered
drawFractionCircos(res[Name == "Col-0_TAL-1_S58"], min_perc = 0, min_pairs = 5)

# Plot fraction by size
test <- merge.data.table(res[, c("Note", "SampleName", "fraction", "Type")], annotation[, c("Note", "Size")], by = "Note", all.x = T, allow.cartesian=TRUE)
ggplot(test[SampleName == "Col-0_TAL_6-2_S58"], aes(as.factor(Size), fraction)) + geom_boxplot() 

othersCircos <- function(dt, size = 367808) {
  sample_name <- unique(dt$Sample)
  cli_alert_info(sample_name)
  circos.par("track.height" = 0.05) # 10% of space for each track
  circos.initialize(sectors = unique(dt$seqnames), xlim = c(0, size))
  
  # Title
  meandepth <- round(mapping_stats[Sample == sample_name]$meandepth)
  title(paste0(unique(dt$Name), "   (Mean depth = ", meandepth, ")"))
  # Repeat annotations
  subset_annotation <- annotation[Note %in% c(unique(dt$read2.Note), unique(dt$read1.Note)), c("seqname", "start", "end", "Note")] 
  subset_annotation[, `:=` (color = "black", line_col = "black", cex = 0.7)]
  cut_sites[, `:=` (color = "red", line_col = "red", cex = 1)]
  subset_annotation <- rbind(subset_annotation, cut_sites)
  
  # circos.track(factors = unique(dt$seqname), ylim = c(0, 1), panel.fun = function(x, y) {
  #   circos.axis(
  #     h="top",                   # x axis on the inner or outer part of the track?
  #     labels=TRUE,               # show the labels of the axis?
  #     major.tick=TRUE,           # show ticks?
  #     labels.cex=0.5,            # labels size (higher=bigger)
  #     labels.font=1,             # labels font (1, 2, 3 , 4)
  #     direction="outside",       # ticks point to the outside or inside of the circle ?
  #     minor.ticks=4,             # Number of minor (=small) ticks
  #     major.tick.length=0.1, # The size of the ticks in percentage of the track height
  #     lwd=1                      # thickness of ticks and x axis.
  #   )
  # })
  
  # Draw repeat annotations
  # Draw label
  circos.labels(subset_annotation$seqname, x = subset_annotation$start, 
                labels = subset_annotation$Note, col = subset_annotation$color, 
                line_col = subset_annotation$line_col, side = "outside", 
                cex = subset_annotation$cex)
  
  circos.track(unique(dt$seqname), ylim = c(0, 1))
  mapply(circos.rect, subset_annotation$start, 0, subset_annotation$end, 1, col = "grey30", border = "grey30")
  #mapply(circos.rect, large_repeats$start, 0, large_repeats$end, 1, col = "steelblue", border = "black")
  mapply(circos.rect, cut_sites$start, 0, cut_sites$end, 1, col = "red", border = "red")
  
  dt$color <- "#20809F20"
  #mapply(circos.link, dt$seqnames, dt$read1, dt$seqnames, dt$read2, col = dt$color, lwd = 2)
  for(i in 1:nrow(dt)){
    circos.link(dt[i]$seqnames, dt[i]$read1, dt[i]$seqnames, dt[i]$read2, col = dt[i]$color, lwd = 2)
  }
}

othersCircos(dt = others[Sample == "Col-0_TAL_6-2_S58"])

drawCircos <- function(dt){
  # Function to add a legend on a circlize plot
  library(ComplexHeatmap)
  library(gridBase)
  library(grid)
  
  # Based on https://jokergoo.github.io/circlize_book/book/legends.html
  plot.new()
  circle_size = unit(1, "snpc") # snpc unit gives you a square region
  pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                        just = c("left", "center")))
  par(omi = gridOMI(), new = TRUE)
  drawFractionCircos(dt = dt, min_perc = 0, min_pairs = 5)
  upViewport()
  # Create the legend colorpalette
  col_fun = colorRamp2(c(0, 0.5, 1), c("grey85", "orange", "red1"))
  
  # Legend Fraction
  lgd_links = Legend(col_fun = col_fun, title_position = "topleft", title = "Fraction")
  
  # Legend Large Repeat
  lgd_large1 = Legend(at = c("Large1 repeat"), type = "line",background = "aquamarine1", border = "black", 
                      title_position = "topleft" 
                      #title = "Repeats"
  )
  
  lgd_large2 = Legend(at = c("Large2 repeat"), type = "line", background = "aquamarine4", border = "black", title_position = "topleft")
  
  lgd_list_vertical = packLegend(lgd_links,lgd_large1, lgd_large2)
  draw(lgd_list_vertical, x = circle_size, just = "left")
}

drawCircos(dt = res[Name == "Col-0_TAL-1_S58"])

# DRAW CIRCOS

# Draw single circos
plyr::l_ply(unique(res$Name), function(i){
  pdf(file = paste0("../Results/Plots/Plots_with_coverage/", i, ".plots.pdf"), width = 10, height = 8)
  try(drawCircos(dt = res[Name == i]))
  dev.off()
})

plyr::l_ply(unique(res$Name), function(i){
  pdf(file = paste0("../Results/Plots/Circos_plots/", i, ".plots.pdf"), width = 10, height = 8)
  try(drawCircos(dt = res[Name == i])) # min_pairs == 5
  dev.off()
})

# DRAW OTHERS
# All in one pdf
pdf(file = "../Results/Plots/All_samples_OTHERS_circos.plots.pdf", width = 10, height = 8)
plyr::l_ply(unique(others$Sample), function(i){
  try(othersCircos(dt = others[Sample == i]))
})
dev.off()

# One pdf per sample
plyr::l_ply(unique(res$Sample), function(i){
  pdf(file = paste0("../Results/Plots/Circos_others/", i, ".others.plots.pdf"), width = 10, height = 8)
  try(othersCircos(dt = others[Sample == i]))
  dev.off()
})
