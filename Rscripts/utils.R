library(data.table)
library(Rsamtools)
library(cli)

# # Function to extract the paired reads position in a table:  
# # readname, seqnames, read1_position, read2_position and insert_size
# readBamRegion <- function(bamPath, seqnames, start, end, min_isize){
#   bamFile <- BamFile(bamPath)
#   cli_alert_info("Loading {.file {bamPath}}.")
#   gr <- GRanges(seqnames, ranges = IRanges(start, end))
#   #params <- ScanBamParam(which = gr, what = scanBamWhat()) # To get the entire bam info, but needs more compute power and ram!
#   params <- ScanBamParam(which = gr, what = c("qname", "pos", "mpos", "isize")) #, flag = scanBamFlag(isFirstMateRead=TRUE))
#   aln <- scanBam(bamFile, param = params)
#   region_id <- paste0(seqnames, ":", start, "-", end)
#   dt <- data.table(rname = aln[[region_id]]$qname, read1 = aln[[region_id]]$pos, read2 = aln[[region_id]]$mpos, isize = aln[[region_id]]$isize, seqnames = seqnames)
#   dt[, `:=` ("read1.start" = read1, 'read1.end' = read1)]
#   dt[, `:=` ("read2.start" = read2, 'read2.end' = read2)]
#   setkey(dt, seqnames, read1.start, read1.end)
#   dt[isize >= min_isize]
# }

# Function to extract the paired reads position in a table:  
# seqname, readname, read1_position, read2_position and insert_size
readBamRegion <- function(bamPath, seqnames, min_isize){
  bamFile <- BamFile(bamPath)
  cli_alert_info("Loading {.file {bamPath}}.")
  params <- ScanBamParam(what = c("rname", "qname", "pos", "mpos", "isize"))
  aln <- scanBam(bamFile, param = params)
  dt <- rbindlist(aln)
  dt[, `:=` ("read1.start" = pos, 'read1.end' = pos)]
  dt[, `:=` ("read2.start" = mpos, 'read2.end' = mpos)]
  setnames(dt, "rname", "seqname")
  setkey(dt, seqname, read1.start, read1.end)
  dt[isize >= min_isize]
}

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


format_si <- function(...) {
  # Format a vector of numeric values according
  # to the International System of Units.
  # http://en.wikipedia.org/wiki/SI_prefix
  #
  # Based on code by Ben Tupper
  # https://stat.ethz.ch/pipermail/r-help/2012-January/299804.html
  # Args:
  #   ...: Args passed to format()
  #
  # Returns:
  #   A function to format a vector of strings using
  #   SI prefix notation
  #
  
  function(x) {
    limits <- c(1e-24, 1e-21, 1e-18, 1e-15, 1e-12,
                1e-9,  1e-6,  1e-3,  1e0,   1e3,
                1e6,   1e9,   1e12,  1e15,  1e18,
                1e21,  1e24)
    prefix <- c("y",   "z",   "a",   "f",   "p",
                "n",   "µ",   "m",   " ",   "Kb",
                "Mb",   "Gb",   "Tb",   "Pb",   "E",
                "Z",   "Y")
    
    # Vector with array indices according to position in intervals
    i <- findInterval(abs(x), limits)
    
    # Set prefix to " " for very small values < 1e-24
    i <- ifelse(i==0, which(limits == 1e0), i)
    
    paste(format(round(x/limits[i], 1),
                 trim=TRUE, scientific=FALSE, ...),
          prefix[i])
  }
}

