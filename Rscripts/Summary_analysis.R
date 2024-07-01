library(data.table)
library(ggplot2)

# Name / Size
A <- fread("../Results/Fraction_table.csv")

unique(A$Note)

AA <- fread("../Results/byGenotype/Fraction_table.csv")

names <- fread("../ReferenceMaterials/Annotations/Samples.csv", header = T)
names <- names[SampleName != ""]

names[SampleName == "WT_control"]
A[SampleName == "WT_control"]

max(A$Repeat_size)

A[Construct == "no-TAL2"]

# Group by size
groups <- seq(0, 600, by=100)
A[, Repeat_group := cut(Repeat_size, breaks = groups), by = Repeat_size]

# by Name and Group 
A[, sum(N_sup1000, na.rm = T), by = .(Name, Genotype, Construct, Suffix, Type, Repeat_group)]

A[Name %in% c("Col-0_no-TAL1_mix")]
A[SampleName %in% c("WT_control")]

# by Name and Type
total <- A[, .(total_sum = sum(N_sup1000, na.rm = T)), by = .(Name, Genotype, Construct, Suffix)]
by_type <- A[, .(type_sum = sum(N_sup1000, na.rm = T)), by = .(Name, Genotype, Construct, Suffix, Type)]
by_repeat_group <- A[, .(repeat_sum = sum(N_sup1000, na.rm = T)), by = .(Name, Genotype, Construct, Suffix, Repeat_group)]

## Merge both table for the ratio
# by_type
B <- merge(total, by_type, by = c("Name", "Genotype", "Construct", "Suffix"))
B[, ratio := round(type_sum / total_sum, 2)]

# by_repeat_group
C <- merge(total, by_repeat_group, by = c("Name", "Genotype", "Construct", "Suffix"))
C[, ratio := round(repeat_sum / total_sum, 2)]

# Boxplot direct / inverted 
ggplot(B[Construct %in% c("no-TAL", "no-TAL1", "no-TAL-2", "TAL-1-WT", "TAL-1", "TAL-2")], aes(Genotype, ratio, fill = Type)) +
  geom_boxplot() + facet_wrap(~ Type, scales = "free_x", ncol = 1) + theme_bw() + 
  scale_y_continuous(labels = scales::percent, limits = c(0, 1))

ggsave("../Results/ByGenotype.pdf", width = 12, height = 10)

ggplot(B[Construct %in% c("no-TAL", "no-TAL1", "no-TAL-2", "TAL-1-WT", "TAL-1", "TAL-2") & !Genotype %in% c("recA3-3-WT")], aes(Genotype, ratio, fill = Construct)) +
  geom_boxplot() + facet_wrap(~ Type, scales = "free_x", ncol = 1) + theme_bw() + 
  scale_y_continuous(labels = scales::percent, limits = c(0, 1))
ggsave("../Results/ByConstruct.pdf", width = 12, height = 10)

ggplot(C[Construct %in% c("no-TAL", "no-TAL1", "no-TAL-2", "TAL-1-WT", "TAL-1", "TAL-2") & Genotype == "Col-0"], aes(Repeat_group, ratio)) +
  geom_boxplot() + theme_bw() + facet_wrap(~ Construct, nrow = 1) + 
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) + ggtitle("Col-0")
ggsave("../Results/ByRepeats_Col0.pdf", width = 14, height = 6)

C[Genotype == "Col-0"]

byNote <- A[Note %in% c("Y", "V", "Q", "L", "G")][, .N, by = c("Note", "Name", "Genotype", "Construct")]
total <- A[Note %in% c("Y", "V", "Q", "L", "G")][, .(total = .N), by = c("Name", "Genotype", "Construct")]
D <- merge(byNote, total, by = c("Name", "Genotype", "Construct"))
D[, perc := round(N / total * 100, 2)]

ggplot(D, aes(x = Note, y = N, fill = Construct)) + geom_boxplot() + facet_wrap(~ Genotype)
