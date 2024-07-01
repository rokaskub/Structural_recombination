library(data.table)
library(ggplot2)

# Name / Size
A <- fread("../Results_genotype/byGenotype/Fraction_table.csv")

A[Name == "WT_control"]

length(unique(A$Sample))
max(A$Repeat_size)

# Group by size
groups <- seq(0, 600, by=100)
A[, Repeat_group := cut(Repeat_size, breaks = groups), by = Repeat_size]

# by Name and Group 
A[, sum(N_sup1000, na.rm = T), by = .(Name, Type, Repeat_group)]

# by Name and Type
total <- A[, .(total_sum = sum(N_sup1000, na.rm = T)), by = .(Name)]
by_type <- A[, .(type_sum = sum(N_sup1000, na.rm = T)), by = .(Name, Type)]
by_repeat_group <- A[, .(repeat_sum = sum(N_sup1000, na.rm = T)), by = .(Name, Repeat_group)]

# Merge both table for the percentage.
A[, perc := round(type_sum / total_sum * 100, 2)]

# 
C <- merge(total, by_repeat_group, by = c("Name", "Genotype", "Construct", "Suffix"))
C[, perc := round(repeat_sum / total_sum * 100, 2)]

table(B$Construct)

ggplot(B[Construct %in% c("no-TAL", "TAL-1-WT", "TAL-1", "TAL-2")], aes(Genotype, perc, fill = Type)) + 
  geom_boxplot() + facet_wrap(~ Type, scales = "free_x", ncol = 1) + theme_bw()
ggsave("../Results_genotype/Plots/ByGenotype.pdf", width = 12, height = 10)

# OLD PART
# Boxplot x-axis = (Genotype + Type) + Construct, y-axis "N" 
# By Genotype
ggplot(B[Construct %in% c("Col-0", "TAL_6-2", "TAL_UP")], aes(Genotype, perc, fill = Type)) + 
  geom_boxplot() + facet_wrap(~ Type, scales = "free_x", ncol = 1) + theme_bw()
ggsave("../Results_genotype/Plots/ByGenotype.pdf", width = 12, height = 10)


ggplot(B[Construct %in% c("Col-0", "TAL_6-2", "TAL_UP")], aes(Genotype, perc, fill = Construct)) + 
  geom_boxplot() + facet_wrap(~ Type, scales = "free_x", ncol = 1) + theme_bw()
ggsave("../Results_genotype/Plots/ByType.pdf", width = 12, height = 10)

ggplot(B[Construct %in% c("Col-0", "TAL_6-2", "TAL_UP")], aes(Genotype, perc, fill = Construct)) + 
  geom_boxplot() + facet_wrap(~ Type, scales = "free_x", ncol = 1) + theme_bw()


# TODO loop for Name
lapply(unique(C$Name), function(name){
  message(name)
  A <- ggplot(C[Construct %in% c("Col-0", "TAL_6-2", "TAL_UP") & Name == name], aes(Repeat_group, perc)) + 
    geom_col()+ theme_bw() + ggtitle(name)
  try(ggsave(plot = A, filename = paste0("../Results_genotype/Plots/", name, ".Repeat_group.byName.pdf"), width = 10, height = 5))
})

# Loop for Genotype
lapply(unique(C$Genotype), function(genotype){
  message(genotype)
  A <- ggplot(C[Construct %in% c("Col-0", "TAL_6-2", "TAL_UP") & Genotype == genotype], aes(Repeat_group, perc)) + 
    geom_boxplot()+ theme_bw() + facet_grid(Genotype ~ Construct) 
  try(ggsave(plot = A, filename = paste0("../Results_genotype/Plots/", genotype, ".byGenomeConstruct.pdf"), width = 10, height = 5))
})

byNote <- A[Note %in% c("Y", "V", "Q", "L", "G")][, .N, by = c("Note", "Name", "Genotype", "Construct")]
total <- A[Note %in% c("Y", "V", "Q", "L", "G")][, .(total = .N), by = c("Name", "Genotype", "Construct")]
D <- merge(byNote, total, by = c("Name", "Genotype", "Construct"))
D[, perc := round(N / total * 100, 2)]

#ggplot(D, aes(x = Note, y = N, fill = Construct)) + geom_boxplot() + facet_wrap(~ Genotype)

