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

# Load repeats information
annotation <- fread("../ReferenceMaterials/Annotations/Repeats_links.csv")
setkey(annotation, seqname, R1_site_start, R1_site_end)

# Load samples information
names <- fread("../ReferenceMaterials/Annotations/Samples.csv", header = T)
names <- names[Genotype != ""]

unique(names$SampleName)

# Load repeats connection type (directed/inverted)
repeats_connection <- fread("../ReferenceMaterials/Annotations/Repeats_links_type.csv")
repeats_connection <- tidyr::pivot_longer(repeats_connection, cols = starts_with("CP"), names_to = "Copy", values_to = "Type")

# Create GRanges object to extract reads from specific regions  
gr1 <- makeGRangesFromDataFrame(annotation, start.field = "R1_site_start", end.field = "R1_site_end", strand.field = "strand")

# Load bam files with insert size <= 1000bp 
is1000 <- list.files("../1_Mapping_bwa_full/bwa_strict/insert_size_below_1000", "*ChrM.MD.IS1000.filtered.bam$", full.names = T)
names(is1000) <- gsub(".ChrM.MD.IS1000.filtered.bam", "", basename(is1000))

# ~ 20-30 min of computation
is1000.dt <- rbindlist(lapply(is1000, readCount, annotation = annotation), idcol = "SampleName")

# Save results
fwrite(is1000.dt, file = "../Results/Repeats_PE_reads_IS_inf1000.raw_reads.csv.gz", nThread = 8, sep = ",")

# Load back if already computed
is1000.dt <- fread("../Results/Repeats_PE_reads_IS_inf1000.raw_reads.csv.gz", nThread = 8)

# Add samples information
is1000.dt <- merge(is1000.dt, names[, c("Genotype", "Construct", "Suffix", "Name", "SampleName")], by = "SampleName", all.x = T)

# Compute mean
mean_isize <- is1000.dt[, .(mean_isize = mean(isize)), by = Name]
setkey(mean_isize, Name)

fwrite(mean_isize, "../Results/Repeats_PE_reads_IS_inf1000.mean_isize.csv", sep = ";")
mean_isize <- fread("../Results/Repeats_PE_reads_IS_inf1000.mean_isize.csv")

#ggplot(is1000.dt, aes(x = isize, fill = Name)) + geom_histogram(binwidth = 100) + theme_bw()

ggplot(is1000.dt, aes(x = isize, color = Name)) + geom_density() + theme_bw() + facet_wrap(~ Genotype) + theme(legend.position = "none") + xlab("Insert size (bp)")
ggsave("../Results/Density.IS_inf1000.pdf", height = 12, width = 16)

ggplot(is1000.dt, aes(x = isize, fill = Name)) + geom_density(alpha = 0.2) + theme_bw() + facet_wrap(~ Genotype) + theme(legend.position = "none") + xlab("Insert size (bp)")
ggsave("../Results/Density_fill.IS_inf1000.pdf", height = 12, width = 16)

is1000.dt[, `:=` ("R1.start" = pos, "R1.end" = pos, "R2.start" = mpos, "R2.end" = mpos)]
setkey(is1000.dt, seqname, R1.start, R1.end)

is1000.dt[repeat_id == "D_1"]
# Create list of repeats to loop on
repeats <- annotation$repeat_id
names(repeats) <- repeats

# Keep reads falling in annotation intervals
by_repeat <- rbindlist(lapply(repeats, function(id){
  R <- annotation[repeat_id == id]
  is1000.dt[pos >= R$R1_site_start & pos <= R$R1_site_end & mpos >= R$R2_site_start & mpos <= R$R2_site_end]
}), idcol = "repeat_id")

summary.dt <- by_repeat[, .N, by = c("repeat_id", "SampleName")]
summary.dt <- merge(summary.dt, names[, c("Genotype", "Construct", "Suffix", "Name", "SampleName")], by = "SampleName", all.x = T)

unique(summary.dt$SampleName)

names[SampleName == "Col-0_TAL_6-2_S1"]

fwrite(summary.dt, "../Results/Repeats_PE_reads_IS_inf1000.count.csv.gz")
summary.dt <- fread("../Results/Repeats_PE_reads_IS_inf1000.count.csv.gz")

# Microhomologies ---------------------------------------------------------
# Get counts
norm.hom <- list.files("../1_Mapping_bwa_full/bwa_normal/microhomologies", "*ChrM.MD.MICROHOMO.NM2_SCLEN20.bowtie2.bam$", full.names = T)
names(norm.hom) <- gsub(".ChrM.MD.MICROHOMO.NM2_SCLEN20.bowtie2.bam", "", basename(norm.hom))

norm.hom.dt <- rbindlist(lapply(norm.hom, function(bam){
  count.dt <- readBamRegion(bamPath = bam, seqnames = "NC_037304.1",  min_isize = 0)
}), idcol = "Name")

norm.hom.dt <- norm.hom.dt[pos >= 200 & pos <= 367608]
# Filtering out this region 194197-198389
norm.hom.dt <- norm.hom.dt[!(read1 >= 194197 & read1 <= 198389)]
norm.hom.dt

# Get count summary
norm.hom.dt.count <- norm.hom.dt[, .N, by = Name]

# MD Coverage files for normalization  ---------------------------------------------------
norm <- list.files("../1_Mapping_bwa_full/bwa_normal", "*.ChrM.MD.bam.coverage.txt$", full.names = T)
names(norm) <- gsub(".ChrM.MD.bam.coverage.txt", "", basename(norm))
norm.dt <- rbindlist(lapply(norm, fread), idcol = "Name")
norm.dt[, factor := round(numreads / 2932880, 3)]

res <- merge(norm.dt, norm.hom.dt.count, by = "Name", all.x = T)
res <- res[!is.na(N)]
res[, CPM := round((N / numreads) * 1e6, 2)]

fwrite(res, "../Results/Microhomologies_norm_table_all_samples.csv", sep = ",")

# MD insertsize filtered bam ----------------------------------------------
# Extract annotations for read1 and read2. Count by pairwise annotations. 
bam_files <- list.files("../1_Mapping_bwa_full/bwa_strict/insert_size_below_1000", "*.ChrM.MD.IS1000.filtered.bam$", full = T)
names(bam_files) <- gsub(".ChrM.MD.IS1000.filtered.bam", "", basename(bam_files))

filtered_bam.dt <- rbindlist(lapply(bam_files, function(bam){
  count.dt <- readBamRegion(bamPath = bam, seqnames = "NC_037304.1", min_isize = 0)
}), idcol = "Name")


fwrite(filtered_bam.dt, "~../Results/Repeats_PE_reads_IS_inf1000.raw_reads.csv.gz")
filtered_bam.dt <- fread("../Results/Repeats_PE_reads_IS_inf1000.raw_reads.csv.gz", nThread = 8)

# Set keys for foverlaps
setkey(filtered_bam.dt, seqname, read1.start, read1.end)
setkey(annotation, seqname, window.start, window.end)

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
unique(B$Pairwise)

#B <- B[!(read1.Note == read2.Note & read1.Number == read2.Number) & isize > 0]
B[read1.Note == read2.Note & read1.Number == 1 & read2.Number == 2]
# Add Copy information (CP1-5)
B[, Copy := "others"]
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

table(B$Copy)

others <- B[Copy == "others"]

# Remove others for pairwise analysis
pairwise_B <- B[Copy != "others"]

# Check
length(B$rname)
length(others$rname) + length(pairwise_B$rname)

pairwise_B <- merge(pairwise_B, names, by.x = "Sample", by.y = "SampleName", all.x = TRUE)
pairwise_B[is.na(Name), Name := Sample]

pairwise_B[, .N, by = c("Sample", "Pairwise")][N > 1]
pairwise_B[Sample == "Col-0-6-2-1E12_S58"][, .N, by = "rname"][order(N)]

View(B[Sample == "Col-0-6-2-1E12_S58"][, .N, by = c("Sample", "Pairwise")])

#B[rname == "VH00615:1:AACG7C7M5:1:2606:57920:13306"]

# Unique B = 717000
pairwise_B[Sample == "Sloan_MSH1_rep1" & isize > 0][, .N, by = c("Sample", "Pairwise")][order(Pairwise)]

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

# Add new names 
summary.dt <- merge.data.table(summary.dt, names, by = "Sample", all.x = T)
# Some are missing so just keep the Sample name
summary.dt[is.na(Name), Name := Sample]

G <- setDT(tidyr::pivot_wider(summary.dt, names_from = "Copy", values_from = "CPM"))

fwrite(G, "../Results/PAIRED_END_READS_IN_REPEATS.inRPM.csv", sep = ",")

# Save single file per sample
for(sample in unique(E$Sample)){
  subset = E[Sample == sample]
  print(sample)
  fwrite(subset, paste0("../Results/", sample, ".count_summary.csv"), sep = ",")
}

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
  cli_alert_info("Total row: {nrow(dt)}, ➚ max. reads: {max(dt$N)} ➘ min. reads: {min(dt$N)}")
  cli_alert_info("Removing everything below max*{min_perc} -> {max(dt$N)*min_perc}")
  # Removes everything below 10% of max value  
  subset_dt <- dt[N > max(N) * min_perc]
  cli_alert_info("Drawing {nrow(subset_dt)} row{?s}")
  my_palette <- colorRampPalette(c("#F8F8F850", "orange", "red"))(max(dt$N))
  subset_dt$colors <- my_palette[cut(subset_dt$N, max(dt$N))]
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
  my_palette <- colorRampPalette(c("#F8F8F850", "orange", "red"))(max(dt$CPM))
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

drawLinkCountCircos(E[Sample == "msh1-6-2-4D10_S6"], min_perc = 0.1)
drawLinkCountCircos(E[Sample == "NGS257-1A7_S27"], min_perc = 0.1)

library(devtools)
install_github("jokergoo/ComplexHeatmap")

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
  col_fun = colorRamp2(c(min(dt$N), max(dt$N)/2, max(dt$N)), c("#F8F8F850", "yellow", "red"))
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

# Draw circos with the legend
drawCircos(dt = summary.dt[Sample == "NGS257-2B11_S49"])
drawCircosCPM(dt = summary.dt[Sample == "NGS257-2B11_S49"])

# To keep legend in the right position use : width = 654 and height = 571 when saving the plot
drawCircos(dt = E[Sample == "NGS257-2B11_S49"])

pdf(file = "../Results/Plots/Circos_with_legend.pdf", width = 8, height = 6)
l_ply(unique(summary.dt$Sample), function(i){
  cli_alert_info(i)
  #pdf(file = paste0("../Results/Plots/Circos_with_legend/", i, ".circos.pdf"))
  try(drawCircosCPM(dt = summary.dt[Sample == i]))
})
dev.off()

