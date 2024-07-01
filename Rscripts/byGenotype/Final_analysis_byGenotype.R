library(rtracklayer)
library(data.table)
library(circlize)
library(ggplot2)
library(scales)
library(tidyr)
library(cli)

# Retrieve sample names
names <- fread("../ReferenceMaterials/Annotations/Samples.csv", header = T)
names <- names[Genotype_construct != "_"]

unique(names$Genotype_construct)

# Paired-end data
res <- fread("../Results_genotype/Results_fraction_COPY.csv") # Fraction 

# Mapping stats
mapping_stats <- fread("../Results_genotype/Coverage.stats.csv")

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
covfiles <- list.files("../byGenotype/1_Mapping_bwa_full/bwa_normal", "*.ChrM.MD.bam.deeptools.100bp_bins.cpm.txt", full.names = T)
names(covfiles) <- gsub(".ChrM.MD.bam.deeptools.100bp_bins.cpm.txt", "", basename(covfiles))
coverage.dt <- rbindlist(lapply(covfiles, fread, col.names = c("seqname", "start", "end", "coverage")), idcol = "Name")
coverage.dt[, max(coverage)]

# Repeat connection file
repeats_connection <- fread("../ReferenceMaterials/Annotations/Repeats_links_type.csv")
repeats_connection <- tidyr::pivot_longer(repeats_connection, cols = starts_with("CP"), names_to = "Copy", values_to = "Type")
repeats_connection <- setDT(repeats_connection)

repeats_connection <- repeats_connection[!is.na(Type)]
repeats_connection[, tmp := gsub("CP", "", Copy)]
repeats_connection[, c("CP1", "CP2") := tstrsplit(tmp, "_")]

# Pairwise
repeats_connection[, Pairwise := paste0(Note, "_", CP1, "-", Note, "_", CP2)]

all_pairwise <- expand.grid(Pairwise = repeats_connection$Pairwise, Name = unique(res$Name))

all_pairwise <- merge(all_pairwise, repeats_connection[, c("Pairwise", "Repeat_size", "Type", "Note")], by = "Pairwise")
res2 <- setDT(merge(all_pairwise, res[, c("Name", "Pairwise", "N_sup1000", "read2.normal", "read1.normal", "available_norm", "fraction")], 
                    by = c("Pairwise", "Name"), all.x = TRUE))

res2[is.na(fraction), fraction := 0]

# Active?
res2[fraction >= 0.5, Active := "yes"]
res2[fraction < 0.5, Active := "no"]
res2[, .N, by = Active]

res2[Note == "Q"]
View(res2[Note == "H"])

fwrite(res2, "../Results_genotype/Fraction_table.csv")
res2 <- fread("../Results_genotype/Fraction_table.csv")

res2[, c("Name", "Construct") := tstrsplit(Name, "_")]

TAL1 <- res2[Construct == "TAL-1"]

ggplot(res2[Note %in% c("Y", "V", "Q", "H", "L", "G")], 
       aes(x = Name, y = fraction, fill = Name)) + 
  geom_boxplot() + facet_wrap(~Note) + theme_bw()

ggplot(TAL1[Note %in% c("Y", "V", "Q", "H", "L", "G")], 
       aes(x = Name, y = fraction, fill = Name)) + 
  geom_boxplot() + facet_wrap(~Note) + theme_bw()


tmp <- melt(dcast(res[, c("Name", "Pairwise", "N_sup1000")], Name ~ Pairwise, value.var = 'N_sup1000', fill = 0), 
            id.vars = c("Name"), variable.name = "Pairwise", value.name = "N_sup1000")
tmp1 <- merge(tmp, res[, c("Name", "Pairwise", "fraction")], by = c("Name", "Pairwise"), all.x = T)
#full_res <- merge.data.table(tmp1, unique(res[, c("Sample", "Construct", "Genotype")]), by = "Sample")

tmp1[N_sup1000 == 0, fraction := 0]
# Active?
tmp1[fraction >= 0.5, Active := "yes"]
tmp1[fraction < 0.5, Active := "no"]

table(tmp1$Active)

merge(tmp1, res[, c("Name", "Note", "Type", "Pairwise")], by = c("Name", "Pairwise"), all.x = TRUE)

ggplot(res2[Note %in% c("Y", "V", "Q", "H", "L", "G")][, .N, by = c("Name", "Active", "Note")], 
       aes(x = Name, y = N, fill = Active)) + geom_col(position = "fill") + facet_wrap(~Note) + theme_bw()

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
         coverage.dt[Name == sample_name]$start,
         col = "steelblue", border = "steelblue")
  
  mean_rescaled = rescale(meandepth, from = c(0, max(coverage.dt[Name == sample_name]$coverage)))
  
  # Draw a line example
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

# Draw the circos, maybe change the min pair since we merge by construct and genotype
drawFractionCircos(res[Name == "Col-0_TAL-1"], min_perc = 0.1, min_pairs = 5)
drawFractionCircos(res[Name == "Col-0_TAL-2"], min_perc = 0.1, min_pairs = 5)
drawFractionCircos(res[Name == "WT_control"], min_perc = 0.1, min_pairs = 5)

# Plot fraction by size
# Draw by sample the boxplot of the 
ggplot(test, aes(Note, fraction)) + geom_boxplot() + facet_wrap(~ Name)

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

drawCircos(dt = res[Name == "Col-0_TAL-1"])
drawCircos(dt = res[Name == "WT_control"])

# DRAW CIRCOS

# Draw single circos
plyr::l_ply(unique(res$Name), function(i){
  pdf(file = paste0("../Results_genotype/Circos_with_coverage/", i, ".plots.pdf"), width = 10, height = 8)
  try(drawCircos(dt = res[Name == i]))
  dev.off()
})

plyr::l_ply(unique(res$Name), function(i){
  pdf(file = paste0("../Results_genotype/Plots_with_coverage", i, ".plots.pdf"), width = 10, height = 8)
  try(drawCircos(dt = res[Name == i]))
  dev.off()
})
