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
setkey(annotation, seqname, R1_site_start, R1_site_end)

# Load samples information
names <- fread("../ReferenceMaterials/Annotations/Samples.csv", header = T)

# Load repeats connection type (directed/inverted)
repeats_connection <- fread("../ReferenceMaterials/Annotations/Repeats_links_type.csv")
repeats_connection <- tidyr::pivot_longer(repeats_connection, cols = starts_with("CP"), names_to = "Copy", values_to = "Type")

# Create GRanges object to extract reads from specific regions  
gr1 <- makeGRangesFromDataFrame(annotation, start.field = "R1_site_start", end.field = "R1_site_end", strand.field = "strand")

readCount <- function(bamPath, annotation){
  cli_alert_info("Extracting reads from {.file {basename(bamPath)}}")
  gr <- makeGRangesFromDataFrame(annotation[, c("seqname", "Note", "R1_site_start", "R1_site_end", "repeat_id")], 
                                 start.field = "R1_site_start", end.field = "R1_site_end", strand.field = "strand", 
                                 keep.extra.columns = T)
  bamFile <- BamFile(bamPath)
  params <- ScanBamParam(which = gr, what = c("qname", "pos", "mpos", "isize")) #, flag = scanBamFlag(isFirstMateRead=TRUE))
  aln <- scanBam(bamFile, param = params)
  names(aln) <- annotation$repeat_id
  res <- rbindlist(aln, idcol = "repeat_id", use.names = T)
  res <- res[isize >= 0] # insert size >= 0
  res[, seqname := "NC_037304.1"]
  return(res)
}

# Load bam files with insert <= 1000bp is1000 <- list.files("~/Projects/Rokas_final/byGenotype/1_Mapping_bwa_full/bwa_strict/insert_size_below_1000", "*ChrM.MD.IS1000.filtered.bam$", full.names = T)
names(is1000) <- gsub(".ChrM.MD.IS1000.filtered.bam", "", basename(is1000))

# ~ 20-30 min of computation
is1000.dt <- rbindlist(lapply(is1000, readCount, annotation = annotation), idcol = "Name")

fwrite(is1000.dt, file = "../Results_genotype/Repeats_PE_reads_IS_inf1000.raw_reads.csv.gz", nThread = 8, sep = ",")
is1000.dt <- fread("../Results_genotype/Repeats_PE_reads_IS_inf1000.raw_reads.csv.gz")

groups <- seq(0, 1000, by=100)
is1000.dt[, isize_group := cut(isize, breaks = groups), by = isize]

ggplot(is1000.dt, aes(x = Name, y = isize)) + geom_boxplot() + theme_bw()
ggsave("../Results_genotype/Plots/Insert_size.inf1000.boxplot.pdf", height = 8, width = 12)
ggplot(is1000.dt, aes(x = isize, fill = Name)) + geom_density(alpha = 0.2) + theme_bw()
ggsave("../Results_genotype/Plots/Insert_size.inf1000.density.pdf", height = 8, width = 12)

mean_isize <- is1000.dt[, .(mean_isize = mean(isize)), by = Name]
setkey(mean_isize, Name)

fwrite(mean_isize, "../Results_genotype/Repeats_PE_reads_IS_inf1000.mean_isize.csv", sep = ";")

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

summary.dt <- by_repeat[, .N, by = c("repeat_id", "Name")]

fwrite(summary.dt, "../Results_genotype/Repeats_PE_reads_IS_inf1000.count.csv.gz")

# A <- foverlaps(is1000.dt, annotation[, c("seqname", "R1_site_start", "R1_site_end", "repeat_id", "Note", "Number")], 
#                by.x = c("seqname", "R1.start", "R1.end"),
#                by.y = c("seqname", "R1_site_start", "R1_site_end"))
# 
# setnames(A, c("Note", "Number"), c("R1.Note", "R1.Number"))
# 
# setkey(annotation, "repeat_id", "R2_site_start", "R2_site_end")
# setkey(A, repeat_id, R2.start, R2.end)
# 
# # Add annotation now for read2
# B <- foverlaps(A, annotation[, c("repeat_id", "R2_site_start", "R2_site_end", "Note", "Number")],
#                by.x = c("repeat_id", "R2.start", "R2.end"),
#                by.y = c("repeat_id", "R2_site_start", "R2_site_end"))
# setnames(B, c("repeat_id", "Note", "Number"), c("read2.repeat_id", "read2.Note", 'read2.Number'))
# 
# 
# 
# # Function to extract the paired reads position in a table:  
# # readname, seqnames, read1_position, read2_position and insert_size
readBamRegion <- function(bamPath, seqnames, start, end, min_isize){
  bamFile <- BamFile(bamPath)
  cli_alert_info("Loading {.file {bamPath}}.")
  gr <- GRanges(seqnames, ranges = IRanges(start, end))
  params <- ScanBamParam(which = gr, what = c("qname", "pos", "mpos", "isize")) #, flag = scanBamFlag(isFirstMateRead=TRUE))
  aln <- scanBam(bamFile, param = params)
  region_id <- paste0(seqnames, ":", start, "-", end)
  dt <- data.table(rname = aln[[region_id]]$qname, 
                   read1 = aln[[region_id]]$pos, 
                   read2 = aln[[region_id]]$mpos, 
                   isize = aln[[region_id]]$isize, 
                   seqnames = seqnames)
  dt[, `:=` ("read1.start" = read1, 'read1.end' = read1)]
  dt[, `:=` ("read2.start" = read2, 'read2.end' = read2)]
  setkey(dt, seqnames, read1.start, read1.end)
  #dt[abs(isize) >= min_isize]
  dt[isize >= min_isize]
}

# Microhomologies ---------------------------------------------------------
# Get counts
norm.hom <- list.files("../byGenotype/1_Mapping_bwa_full/bwa_normal/microhomologies", "*ChrM.MD.MICROHOMO.NM2_SCLEN20.bowtie2.bam$", full.names = T)
names(norm.hom) <- gsub(".ChrM.MD.MICROHOMO.NM2_SCLEN20.bowtie2.bam", "", basename(norm.hom))

norm.hom.dt <- rbindlist(lapply(norm.hom, function(bam){
  count.dt <- readBamRegion(bamPath = bam, seqnames = "NC_037304.1", start = 200, end = 367608, min_isize = 0)
}), idcol = "Name")

# Filtering out this region 194197-198389
norm.hom.dt <- norm.hom.dt[!(read1 >= 194197 & read1 <= 198389)]
norm.hom.dt

# Get count summary
norm.hom.dt.count <- norm.hom.dt[, .N, by = Name]

# MD Coverage files for normalization  ---------------------------------------------------
norm <- list.files("../byGenotype/1_Mapping_bwa_full/bwa_normal", "*.ChrM.MD.bam.coverage.txt$", full.names = T)
names(norm) <- gsub(".ChrM.MD.bam.coverage.txt", "", basename(norm))
norm.dt <- rbindlist(lapply(norm, fread), idcol = "Name")

res <- merge(norm.dt, norm.hom.dt.count, by = "Name", all.x = T)
res <- res[!is.na(N)]
res[, CPM := round((N / numreads) * 1e6, 2)]

fwrite(res, "../Results_genotype/Microhomologies_norm_table_all_samples.csv", sep = ",")

# MD insertsize filtered bam ----------------------------------------------
# Extract annotations for read1 and read2. Count by pairwise annotations. 
bam_files <- list.files("../byGenotype/1_Mapping_bwa_full/bwa_strict/insert_size_below_1000", "*.ChrM.MD.IS1000.filtered.bam$", full = T)
names(bam_files) <- gsub(".ChrM.MD.IS1000.filtered.bam", "", basename(bam_files))

filtered_bam.dt <- rbindlist(lapply(bam_files, function(bam){
  count.dt <- readBamRegion(bamPath = bam, seqnames = "NC_037304.1", start = 200, end = 367608, min_isize = 0)
}), idcol = "Name")

setkey(filtered_bam.dt, seqnames, read1.start, read1.end)
setkey(annotation, seqname, window.start, window.end)

fwrite(filtered_bam.dt, "../Results_genotype/Repeats_PE_reads_IS_inf1000.raw_reads.csv.gz")

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

# Other
# TODO: draw the circos for the others (link only)
others <- B[Copy == "others"]

# Remove others for pairwise analysis
pairwise_B <- B[Copy != "others"]

# Check
length(B$rname)
length(others$rname) + length(pairwise_B$rname)

pairwise_B <- merge(pairwise_B, names, by = "Sample", all.x = TRUE)
pairwise_B[is.na(Name), Name := Sample]

# Still few of them are overlapping
#pairwise_B[duplicated(pairwise_B$rname)]
#pairwise_B[rname == "VH00615:1:AACG7C7M5:1:2502:42166:46056"]

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

# TODO merge with the all_repeats_connection_type
# 1. Mean of the fraction + try boxplot of inverted vs direct 
# 2. CPM normalized * (1 - norm_factor)


fwrite(G, "../Results_genotype/PAIRED_END_READS_IN_REPEATS.inRPM.csv", sep = ",")

# Repeat1 id # repeat2 id 2
summary.dt[Sample == "msh1-6-2-4D10_S6"]

# Save single file per sample
for(sample in unique(E$Sample)){
  subset = E[Sample == sample]
  print(sample)
  fwrite(subset, paste0("../Results_genotype/", sample, ".count_summary.csv"), sep = ",")
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


#library(devtools)
#install_github("jokergoo/ComplexHeatmap")
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

pdf(file = "../Results_genotype/Plots/Circos_with_legend.pdf", width = 8, height = 6)
l_ply(unique(summary.dt$Sample), function(i){
  cli_alert_info(i)
  #pdf(file = paste0("/home/dpflieger/Project/Rokas/Plots/Circos_with_legend/", i, ".circos.pdf"))
  try(drawCircosCPM(dt = summary.dt[Sample == i]))
})
dev.off()

# TEST ZONE ---------------------------------------------------------------
getCircos <- function(dt, size = 367808) {
  dt <- dt[!is.na(read1)]
  circos.par("track.height" = 0.1) # 10% of space for each track
  circos.initialize(unique(dt$seqnames), xlim = c(0, size))
  
  #circos.labels(annotation$seqname, x = annotation$start, labels = annotation$Note, de = "outside")
  # circos.track(dt$seqnames, y = dt$read2,
  #              panel.fun = function(x, y) {
  #                circos.text(CELL_META$xcenter,
  #                            CELL_META$cell.ylim[2] + mm_y(5),
  #                            CELL_META$sector.index)
  #                circos.axis(labels.cex = 0.6)
  #              })
  circos.labels(annotation$seqname, x = annotation$start, labels = annotation$Note, side = "outside")
  circos.track(ylim = c(0, 1))
  mapply(circos.rect, annotation$start, 0, annotation$end, 1)
  #circos.trackPoints(dt$sectors, dt$read1, dt$read2, pch = 16, cex = 0.5)
  for(i in 1:nrow(dt)) {
    circos.link(dt[i]$seqname, dt[i]$read1, dt[i]$seqname, dt[i]$read2, col = dt[i]$color, lwd = 5)
  }
}

bamFiles <- list.files("../byGenotype/1_Mapping_bwa_full/bwa_normal/insert_size", "*.ChrM.MD.IS1000.filtered.bam$", full.names = T)
names(bamFiles) <- gsub(".ChrM.MD.IS1000.filtered.bam", "", basename(bamFiles))

full_res <- lapply(bamFiles, readBamRegion, seqnames = "NC_037304.1", start = 1, end = 367808, min_isize = 1000)

# Test circos-like plot ------------------------------------------------------
l_ply(names(full_res), function(i){
  cli_alert_info(i)
  full_res[[i]][, color := "#20809F05"]
  pdf(file = paste0("../Results_genotype/Plots/", i, ".withSmallLabels.pdf"))
  #A = ggraph(full_res[[i]][read1 >= 600 & read2 <= 367000], layout = 'linear') + geom_edge_arc(alpha = 0.05) + ggtitle(i)
  #print(A)
  getCircos(full_res[[i]])
  dev.off()
})

