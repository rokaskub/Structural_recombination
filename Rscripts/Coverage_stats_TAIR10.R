library(data.table)
library(ggplot2)

samples <- fread("../ReferenceMaterials/Annotations/Samples.csv")

myfiles <- list.files("../1_Mapping_bwa_full/bwa_normal", glob2rx("*.TAIR10.bam.coverage.txt"), full = T)
names(myfiles) <- gsub(".ChrM.TAIR10.bam.coverage.txt", "", basename(myfiles))

# Load all the files 
cov.dt <- rbindlist(lapply(myfiles, fread), idcol = "Sample") 
cov.dt <- merge(cov.dt, samples[, c("SampleName", "Name")], by.x = "Sample", by.y = "SampleName", all.x = TRUE)
setnames(cov.dt, "#rname", "Chromosome")

# Save file
fwrite(cov.dt[, c("Sample", "Name", "Chromosome", "numreads", "meandepth", "meanbaseq")], "../Results/Coverage.stats.csv")

# Stats by chromo
res.dt <- dcast(cov.dt[, c("Sample", "Name", "Chromosome", "numreads")], Sample + Name ~ Chromosome, value.var = "numreads")
res.dt[, Total_reads := Chr1 + Chr2 + Chr3 + Chr4 + Chr5 + ChrC + ChrM]
fwrite(res.dt, "../Results/Coverage.stats_chromo.csv")




