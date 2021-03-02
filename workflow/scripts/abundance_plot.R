#the genera to look for:
gen <- c("Roseomonas", "Nocardia", "Anabaena", "Paracoccus", "Rhodococcus")

###################
##FOR SHORT READS##
###################

#sourmash

sour <- read.csv(snakemake@input[["sourmash"]][1],
                 header = T, stringsAsFactors = F)

sel.sour <- sour[sour$genus %in% gen,]
sel.sour <- sel.sour[!duplicated(sel.sour[,"genus"]),]

sour.df <- sel.sour[,c("count", "genus")]

if (any(!gen %in% sour.df$genus)) {
  add <- data.frame("count" = 0,
                    genus = gen[!gen %in% sour.df$genus])
  sour.fin <- rbind(sour.df, add)
  sour.fin$method <- "sourmash"
} else {
  sour.fin <- sour.df
  sour.fin$method <- "sourmash"
}

#kraken2

kra <- read.table(snakemake@input[["kraken2"]][1],
                  header = F, sep = "\t")

kra.g <- kra[kra$V4 == "G",]
kra.g$V6 <- gsub( " ", "", kra.g$V6) 
sel.kra <- kra.g[kra.g$V6 %in% gen,] #v2 is the count
kra.fin <- sel.kra[,c("V2", "V6")]
colnames(kra.fin) <- c("count", "genus")
kra.fin$method <- "kraken2"

#kallisto - (waiting)

#final table
fin <- rbind(sour.fin, kra.fin)
#fin$count <- log2(fin$count + 1)

#plot
library(ggplot2)
p <- ggplot(data = fin, aes(x = genus, y = count, group = count, fill = genus)) + ylim(0, 15000)
p <- p + geom_bar(stat = "identity", width = 0.5, position = "dodge")
p <- p + facet_grid(. ~ method)
p <- p + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 90)) 

# ggsave(filename = paste0("results/final_abundance/","Sample",as.character(snakemake@params[["n_samples"]]),
#                          "_fraction",snakemake@params[["fraction"]],"_total_abundance.pdf"))
ggsave(filename = paste0("results/final_abundance/second_try/sr/","Sample",as.character(snakemake@params[["n_samples"]]),
                         "_fraction",snakemake@params[["fraction"]],"_total_abundance.pdf"))

#ggplot(fin, aes(x=method, y=count, color=genus)) + geom_point() + theme_bw()

###################
##FOR LONG READS##
###################

#sourmash - no match at genus level

#kraken2

# kra <- read.table(snakemake@input[["kraken2"]][1],
#                   header = F, sep = "\t")
# kra <- read.table("/home/uzuner/Documents/ReadSimulation/Tumor-microbiome-calling/workflow/results/kraken2/lr/bacterial-db/evol1_Sample1_fraction10000_bracken_species",header = F, sep = "\t")
# 
# kra.g <- kra[kra$V4 == "G",]
# kra.g$V6 <- gsub( " ", "", kra.g$V6) 
# sel.kra <- kra.g[kra.g$V6 %in% gen,] #v2 is the count
# kra.fin <- sel.kra[,c("V2", "V6")]
# colnames(kra.fin) <- c("count", "genus")
# kra.fin$method <- "kraken2"
# 
# #kallisto
# 
# #final table
# fin <- rbind(sour.fin, kra.fin)
# #fin$count <- log2(fin$count + 1)
# 
# #plot
# library(ggplot2)
# p <- ggplot(data = fin, aes(x = genus, y = count, group = count, fill = genus)) + ylim(0, 15000)
# p <- p + geom_bar(stat = "identity", width = 0.5, position = "dodge")
# p <- p + facet_grid(. ~ method)
# p <- p + theme_bw()
# p <- p + theme(axis.text.x = element_text(angle = 90)) 
# 
# # ggsave(filename = paste0("results/final_abundance/","Sample",as.character(snakemake@params[["n_samples"]]),
# #                          "_fraction",snakemake@params[["fraction"]],"_total_abundance.pdf"))
# ggsave(filename = paste0("results/final_abundance/second_try/sr/","Sample",as.character(snakemake@params[["n_samples"]]),
#                          "_fraction",snakemake@params[["fraction"]],"_total_abundance.pdf"))
