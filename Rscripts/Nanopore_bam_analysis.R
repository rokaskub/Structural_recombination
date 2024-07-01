library(Rsamtools)
library(data.table)

source("utils.R")

bed <- fread("../Nanopore/rokas.nanopore.minimap2.ChrM.sorted.sclen200.bam.bed", col.names = c("Chr", "Start", "End", "qname", "X", "strand", "cigar"))
bam <- BamFile("../Nanopore/rokas.nanopore.minimap2.ChrM.sorted.sclen200.bam")

params <- ScanBamParam(what = c("rname", "qname", "qwidth", "pos", "cigar"))
aln <- scanBam(bam, param = params)
res.dt <- data.table("rname" = aln[[1]]$rname, 
                     "qname" = aln[[1]]$qname, 
                     "qwidth" = aln[[1]]$qwidth, 
                     "pos" = aln[[1]]$pos, 
                     "cigar" = aln[[1]]$cigar)

# Function to extract softclip sizes
extract_softclip_sizes <- function(cigar_string) {

  # Match softclips at the beginning
  first_softclip <- regmatches(cigar_string, regexec("(\\d+)S", cigar_string))

  # If first softclip is found -> first_softclip = list(c("15S" "15"))
  # If no softclip found -> list(charactor(0))
  if (length(first_softclip[[1]]) > 1) {
    first_softclip_size <- as.integer(first_softclip[[1]][2])
  } else {
    first_softclip_size <- 0
  }
  
  # Match softclips at the end
  end_softclip <- regmatches(cigar_string, regexec(".*?(\\d+)S$", cigar_string))
  
  # If end softclip is found
  if (length(end_softclip[[1]]) > 1) {
    end_softclip_size <- as.integer(end_softclip[[1]][2])
  } else {
    end_softclip_size <- 0
  }
  
  return(list(first_softclip = first_softclip_size, end_softclip = end_softclip_size))
}

# Function to extract softclip sizes
extract_softclip_sizes <- function(cigar_string) {
  # Pattern to match softclips
  softclip_pattern <- "^(\\d+)S"
  
  # Match softclips at the beginning
  first_softclip <- regmatches(cigar_string, regexec(softclip_pattern, cigar_string))
  
  # Initialize sizes
  first_softclip_size <- 0
  end_softclip_size <- 0
  
  # If first softclip is found
  if (length(first_softclip[[1]]) > 1) {
    first_softclip_size <- as.integer(first_softclip[[1]][2])
  }
  
  # Match softclips at the end
  end_softclip <- regmatches(cigar_string, regexec(paste0("(\\d+)S$"), cigar_string))
  
  # If end softclip is found
  if (length(end_softclip[[1]]) > 1) {
    end_softclip_size <- as.integer(end_softclip[[1]][2])
  }
  
  return(list(first_softclip = first_softclip_size, end_softclip = end_softclip_size))
}


# Function to calculate mapped read length from CIGAR string
calculate_mapped_read_length <- function(cigar_string) {
  # Pattern to match alignment segments (M) and deletions (D)
  alignment_pattern <- "(\\d+)(M|D)"
  
  # Extract alignment segments and deletions
  alignment_segments <- regmatches(cigar_string, gregexpr(alignment_pattern, cigar_string))[[1]]
  
  # Initialize mapped read length
  mapped_read_length <- 0
  
  # Iterate through alignment segments and deletions
  for (segment in alignment_segments) {
    # Extract length and operation
    segment_info <- regmatches(segment, regexec(alignment_pattern, segment))[[1]]
    length <- as.integer(segment_info[2])
    operation <- segment_info[3]
    
    # Update mapped read length based on operation
    if (operation == "M" || operation == "D") {
      mapped_read_length <- mapped_read_length + length
    }
  }
  
  return(mapped_read_length)
}

# Example CIGAR string (511 length mapped)
cigar_string <- "164S11M5I92M2D120M2D143M2D3M1D22M1D2M1D44M2I17M1D14M1D32M9S"
cigar_string <- "180S18M1D45M"
cigar_string <- "18M1D100M15S"

# Calculate mapped read length
calculate_mapped_read_length(cigar_string)
extract_softclip_sizes(cigar_string)

bed[, c("first_softclip", "end_softclip") := extract_softclip_sizes(cigar), by = cigar]
bed[, "mapped_part" := calculate_mapped_read_length(cigar), by = cigar]

bed[qname == "2e7b9c99-aec4-471e-bf1e-9a6edf0fed42"]
res.dt[qname == "2e7b9c99-aec4-471e-bf1e-9a6edf0fed42"]

# Binwidth 500bp
ggplot() + 
  geom_histogram(data = bed[first_softclip >= 200,], aes(x = Start, fill = "5'"), binwidth = 500, alpha = 0.5)+
  geom_histogram(data = bed[end_softclip >= 200,], aes(x = End, fill = "3'"), binwidth = 500, alpha = 0.5) +
  scale_x_continuous(labels = format_si()) + theme_bw() + ylab("Nanopore read counts") + xlab("Genomic positions") +
  ggtitle(label = "Nanopore softclip positions (3' and 5') in reads after mapping") +
  scale_fill_manual(values=c("5'" = "blue", "3'" = "red"), 'Position')

# Bp resolution using geom line 
tmp <- merge(bed[first_softclip >= 200, .N, by = c("Start")], bed[end_softclip >= 200, .N, by = c("End")], by.x = "Start", by.y = "End", all = T)
tmp[, coverage := sum(N.x, N.y, na.rm = T), by = c("N.x", "N.y")]
ggplot(tmp, aes(x = Start, y = coverage)) + geom_line() + theme_bw()

