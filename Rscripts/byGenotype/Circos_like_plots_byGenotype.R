library(rtracklayer)
library(data.table)
library(circlize)
library(ggplot2)
library(viridis)
library(ggraph)
library(igraph)
library(plyr)
library(cli)

# Load annotations files --------------------------------------------------

# Annotation file of repeats of interest 
annotation <- fread("../ReferenceMaterials/Annotations/Repeats_links.csv")

# SSA Annotations
SSA_annotation <- fread("../ReferenceMaterials/Annotations/SSA_window_annotation_updated.csv")
table(SSA_annotation$Note)
SSA_annotation[, repeat_id := paste(Note, Number, sep = "_")]
SSA_annotation[, fragment_size := (end - start) + 1] 
SSA_annotation <- SSA_annotation[!is.na(start)]
setkey(SSA_annotation, seqname, R1_site_start, R1_site_end)

# GFF of repeats to draw specific stuff
repeats_annotations <- as.data.table(import("../ReferenceMaterials/Annotations/Repeats_region.gff3"))
repeats_annotations[, Note := tstrsplit(Note, "-", keep = 1)]
others_repeats_annotations <- repeats_annotations[!Note %in% annotation$Note]

# Retrieve the repeats connection type (directed/inverted)
repeats_connection <- fread("../ReferenceMaterials/Annotations/Repeats_links_type.csv")
repeats_connection <- tidyr::pivot_longer(repeats_connection, cols = starts_with("CP"), names_to = "Copy", values_to = "Type")

# "Large" repeats for circos plot 
large_repeats_annotations <- others_repeats_annotations[Note %like% "Large"]

# To draw the cut sites
cut_sites <- data.table(seqname = c("NC_037304.1", "NC_037304.1"), start = c(265238, 266820), end = c(265292, 266880), Note = c("6-2", "UP"))
 
## Merge SSA results with the others results 
SSA_inf1000 <- fread("../Results_genotype/SSA_PE_reads_IS_inf1000.count.csv.gz", col.names = c("repeat_id", "Name", "N"))
SSA_sup1000 <- fread("../Results_genotype/SSA_PE_reads_IS_sup1000.count.csv.gz") %>% setnames("Sample", "Name")
Repeats_inf1000 <- fread("../Results_genotype/Repeats_PE_reads_IS_inf1000.count.csv.gz")
Repeats_sup1000 <- fread("../Results_genotype/Repeats_PE_reads_IS_sup1000.count.csv.gz")

inf1000 <- rbindlist(list(SSA_inf1000, Repeats_inf1000))
sup1000 <- rbindlist(list(SSA_sup1000, Repeats_sup1000))

names(SSA_sup1000)
names(Repeats_sup1000)

setnames(inf1000, "N", "N_inf1000")
setnames(sup1000, "N", "N_sup1000")

names(sup1000)
names(inf1000)

tmp <- merge(sup1000, inf1000, 
             by.x = c("Name", "read2.repeat_id"), 
             by.y = c("Name", "repeat_id"), 
             all.x = T)
setnames(tmp, "N_inf1000", "read2.normal")

tmp <- merge(tmp, inf1000, 
             by.x = c("Name","read1.repeat_id"), 
             by.y = c("Name", "repeat_id"), 
             all.x = T)
setnames(tmp, "N_inf1000", "read1.normal")

tmp[is.na(read2.normal), read2.normal := 0]
tmp[is.na(read1.normal), read1.normal := 0]
tmp[, available_norm := ifelse(read2.normal > read1.normal, read1.normal, read2.normal)]
tmp[, fraction := round(N_sup1000  / (N_sup1000 + available_norm), 3)]

# Add repeats connection information
res <- merge(tmp, repeats_connection[, c("Note", "Copy", "Type")], by = c("Note", "Copy"), all.x = T)
fwrite(res, "../Results_genotype/Results_fraction_COPY.csv")

ggplot(res[Note == "L"], aes(x = Name, y = fraction)) + geom_point(position = "dodge") #+ coord_flip()
ggplot(res[Note == "L"], aes(x = Sample, y = fraction, fill = Name)) + geom_boxplot() + facet_wrap(~Construct)

# Circlize plots ----------------------------------------------------------
drawLinkCircos <- function(dt, size = 367808) {
  circos.par("track.height" = 0.1) # 10% of space for each track
  circos.initialize(sectors = unique(dt$seqname), xlim = c(0, size))
  subset_annotation <- annotation[Note %in% dt$read2.Note]
  print(subset_annotation)
  circos.labels(subset_annotation$seqname, x = subset_annotation$start, labels = subset_annotation$Note, side = "outside")
  circos.track(ylim = c(0, 1))
  
  mapply(circos.rect, subset_annotation$start, 0, subset_annotation$end, 1)
  dt$color <- "#20809F10"
  #mapply(circos.link, dt$seqname, dt$read1, dt$seqname, dt$read2, col = dt$color, lwd = 2)
  for(i in 1:nrow(dt)) {
    circos.link(dt[i]$seqname, dt[i]$read1, dt[i]$seqname, dt[i]$read2, col = dt[i]$color, lwd = 1)
  }
}

drawLinkCountCircos <- function(dt, size = 367808, min_perc = 0.1) {
  cli_alert_info("Total row: {nrow(dt)}, ➚ max. reads: {max(dt$N_sup1000)} ➘ min. reads: {min(dt$N_sup1000)}")
  cli_alert_info("Removing everything below max*{min_perc} -> {max(dt$N)*min_perc}")
  # Removes everything below 10% of max value  
  subset_dt <- dt[N_sup1000 > max(N_sup1000) * min_perc]
  cli_alert_info("Drawing {nrow(subset_dt)} row{?s}")
  my_palette <- colorRampPalette(c("lightgrey", "orange", "red"))(max(dt$N_sup1000))
  subset_dt$colors <- my_palette[cut(subset_dt$N_sup1000, max(dt$N_sup1000))]
  print(subset_dt)
  
  circos.par("track.height" = 0.1) # 10% of space for each track
  circos.initialize(sectors = unique(dt$seqname), xlim = c(0, size))
  subset_annotation <- annotation[Note %in% subset_dt$Note] # TODO: Filter name 
  circos.labels(subset_annotation$seqname, x = subset_annotation$start, labels = subset_annotation$Note, side = "outside")
  circos.track(ylim = c(0, 1))
  title(unique(subset_dt$Name))
  
  mapply(circos.rect, subset_annotation$start, 0, subset_annotation$end, 1)
  #mapply(circos.link, dt$seqname, dt$read1, dt$seqname, dt$read2, col = dt$color, lwd = 2)
  for(i in 1:nrow(subset_dt)) {
    # circos.link(sector.index1, c(0, 1), sector.index2, 0, col, lwd, lty, border)
    #circos.link(dt[i]$seqname, dt[i]$read1, dt[i]$seqname, dt[i]$read2, col = dt[i]$color, lwd = 2)
    circos.link(subset_dt[i]$seqname, c(subset_dt[i]$CP1.start, subset_dt[i]$CP1.end),
                subset_dt[i]$seqname, c(subset_dt[i]$CP2.start, subset_dt[i]$CP2.end),
                col = subset_dt[i]$color, lwd = 2)
  }
}

drawLinkCPMCircos <- function(dt, size = 367808, min_perc = 0.1) {
  cli_alert_info("Total row: {nrow(dt)}, ➚ max. reads: {max(dt$CPM)} ➘ min. reads: {min(dt$CPM)}")
  cli_alert_info("Removing everything below max*{min_perc} -> {max(dt$CPM)*min_perc}")
  # Removes everything below 10% of max value  
  subset_dt <- dt[CPM > (max(CPM) * min_perc)]
  cli_alert_info("Drawing {nrow(subset_dt)} row{?s}")
  my_palette <- colorRampPalette(c("lightgrey", "orange", "red"))(max(dt$CPM))
  subset_dt$colors <- my_palette[cut(subset_dt$CPM, max(dt$CPM))]
  print(subset_dt)
  
  circos.par("track.height" = 0.1) # 10% of space for each track
  circos.initialize(sectors = unique(dt$seqname), xlim = c(0, size))
  subset_annotation <- annotation[Note %in% subset_dt$Note] # TODO: Filter name 
  circos.labels(subset_annotation$seqname, x = subset_annotation$start, labels = subset_annotation$Note, side = "outside")
  circos.track(ylim = c(0, 1))
  title(unique(subset_dt$Name))
  
  mapply(circos.rect, subset_annotation$start, 0, subset_annotation$end, 1)
  #mapply(circos.link, dt$seqname, dt$read1, dt$seqname, dt$read2, col = dt$color, lwd = 2)
  for(i in 1:nrow(subset_dt)) {
    # circos.link(sector.index1, c(0, 1), sector.index2, 0, col, lwd, lty, border)
    #circos.link(dt[i]$seqname, dt[i]$read1, dt[i]$seqname, dt[i]$read2, col = dt[i]$color, lwd = 2)
    circos.link(subset_dt[i]$seqname, c(subset_dt[i]$CP1.start, subset_dt[i]$CP1.end),
                subset_dt[i]$seqname, c(subset_dt[i]$CP2.start, subset_dt[i]$CP2.end),
                col = subset_dt[i]$color, lwd = 2)
  }
}

drawLinkFractionCircos <- function(dt, size = 367808, min_perc = 0.1) {
  cli_alert_info("Total row: {nrow(dt)}, ➚ max. reads: {max(dt$N_sup1000)} ➘ min. reads: {min(dt$N_sup1000)}")
  cli_alert_info("Removing everything below max*{min_perc} -> {max(dt$N_sup1000)*min_perc}")
  # Removes everything below 10% of max value  
  subset_dt <- dt[N_sup1000 > max(N_sup1000) * min_perc]
  
  cli_alert_info("Drawing {nrow(subset_dt)} row{?s}")
  
  my_palette <- colorRampPalette(c("lightgrey", "orange", "red"))(100)
  subset_dt$colors <- my_palette[cut(subset_dt$fraction, 100)]
  
  circos.par("track.height" = 0.1) # 10% of space for each track
  
  circos.initialize(sectors = unique(dt$seqname), xlim = c(0, size))
  title(paste0(unique(subset_dt$Name), "   (N = ", nrow(subset_dt), ")"))

  circos.track(factors = unique(dt$seqname), ylim = c(0, 1), panel.fun = function(x, y) {
    circos.axis(
      h="top",                   # x axis on the inner or outer part of the track?
      labels=TRUE,               # show the labels of the axis?
      major.tick=TRUE,           # show ticks?
      labels.cex=0.5,            # labels size (higher=bigger)
      labels.font=1,             # labels font (1, 2, 3 , 4)
      direction="outside",       # ticks point to the outside or inside of the circle ?
      minor.ticks=4,             # Number of minor (=small) ticks
      major.tick.percentage=0.1, # The size of the ticks in percentage of the track height
      lwd=1                      # thickness of ticks and x axis.
    )
  })
  
  subset_annotation <- annotation[Note %in% c(subset_dt$Note, "Y")]  

  #todraw <- rbind(cut_sites, subset_annotation[, c("seqnames", "start", "end", "Note")])
  circos.labels(subset_annotation$seqname, x = subset_annotation$start, labels = subset_annotation$Note, side = "outside")
  #circos.labels(cut_sites$seqname, x = cut_sites$start, labels = cut_sites$Note, col = "blue", side = "outside")
  
  setorder(subset_dt, N_sup1000)
  print(subset_dt)
  mapply(circos.rect, subset_annotation$start, 0, subset_annotation$end, 1)
  mapply(circos.rect, cut_sites$start, 0, cut_sites$end, 1, col = "red")
  #mapply(circos.link, dt$seqname, dt$read1, dt$seqname, dt$read2, col = dt$color, lwd = 2)
  circos.track(ylim = c(0, 1))
  for(i in 1:nrow(subset_dt)) {
    # circos.link(sector.index1, c(0, 1), sector.index2, 0, col, lwd, lty, border)
    #circos.link(dt[i]$seqname, dt[i]$read1, dt[i]$seqname, dt[i]$read2, col = dt[i]$color, lwd = 2)
    circos.link(subset_dt[i]$seqname, c(subset_dt[i]$CP1.start, subset_dt[i]$CP1.end),
                subset_dt[i]$seqname, c(subset_dt[i]$CP2.start, subset_dt[i]$CP2.end),
                col = subset_dt[i]$color, lwd = 2)
  }
}

testCircos <- function(dt, size = 367808, min_perc = 0.1) {
  circos.par("track.height" = 0.1) # 10% of space for each track
  circos.initialize(sectors = unique(dt$seqname), xlim = c(0, size))
  title(paste0(unique(dt$Name), "   (N = ", nrow(dt), ")"))
  
  circos.track(factors = unique(dt$seqname), ylim = c(0, 1), panel.fun = function(x, y) {
    subset_annotation <- annotation[Note %in% c(dt$Note, "Y")]  
    #circos.labels(subset_annotation$seqname, x = subset_annotation$start, labels = subset_annotation$Note, side = "outside")
    circos.axis(
      h="top",                   # x axis on the inner or outer part of the track?
      labels=TRUE,               # show the labels of the axis?
      major.tick=TRUE,           # show ticks?
      labels.cex=0.5,            # labels size (higher=bigger)
      labels.font=1,             # labels font (1, 2, 3 , 4)
      direction="outside",       # ticks point to the outside or inside of the circle ?
      minor.ticks=4,             # Number of minor (=small) ticks
      major.tick.length=0.1, # The size of the ticks in percentage of the track height
      lwd=1                      # thickness of ticks and x axis.
    )
    mapply(circos.rect, subset_annotation$start, 0, subset_annotation$end, 1)
    mapply(circos.rect, cut_sites$start, 0, cut_sites$end, 1, col = "red", border = NA)
  })
  
  #setorder(dt, N_sup1000)
  #mapply(circos.link, dt$seqname, dt$read1, dt$seqname, dt$read2, col = dt$color, lwd = 2)
  
  #circos.track(ylim = c(0, 1))
  
}

drawLinkFractionCircos(res[Name == "Col-0_TAL-1"])
drawLinkCPMCircos(res[Name == "Col-0_TAL-1"], min_perc = 0.1)

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
  drawLinkCountCircos(dt = dt)
  upViewport()
  # Create the legend colorpalette
  col_fun = colorRamp2(c(min(dt$N), max(dt$N)/2, max(dt$N)), c("lightgrey", "yellow", "red"))
  # Legend
  lgd_links = Legend(col_fun = col_fun, title_position = "topleft", title = "# links")
  lgd_list_vertical = packLegend(lgd_links)
  draw(lgd_list_vertical, x = circle_size, just = "left")
}

drawCircosCPM <- function(dt){
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
  drawLinkCPMCircos(dt = dt)
  upViewport()
  # Create the legend colorpalette
  col_fun = colorRamp2(c(min(dt$CPM), max(dt$CPM)/2, max(dt$CPM)), c("lightgrey", "yellow", "red"))
  # Legend
  lgd_links = Legend(col_fun = col_fun, title_position = "topleft", title = "Read pairs (in CPM)")
  lgd_list_vertical = packLegend(lgd_links)
  draw(lgd_list_vertical, x = circle_size, just = "left")
}

pdf(file = "../Results_genotype/Plots/Circos_with_legend.pdf", width = 8, height = 6)
l_ply(unique(summary.dt$Sample), function(i){
  cli_alert_info(i)
  #pdf(file = paste0("/home/dpflieger/Project/Rokas/Plots/Circos_with_legend/", i, ".circos.pdf"))
  try(drawCircosCPM(dt = summary.dt[Sample == i]))
})
dev.off()
