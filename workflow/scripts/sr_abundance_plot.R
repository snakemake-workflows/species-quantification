library(ggplot2)

#this script creates an abundance plot for all samples and compares existing methods.
#the species to look for:
#removed "Candidatus Carsonella ruddii"

sp <- c("Yersinia pestis",
        "Myxococcus xanthus",
        "Deinococcus radiodurans",
        "Salmonella enterica",
        "Pseudomonas aeruginosa",
        "Roseomonas mucosa")

#sourmash#

fin.s <- data.frame()
for (i in snakemake@input[["sourmash"]]){
  #read and match the genera
  sour <- read.csv(i)
  sel.sour <- sour[sour$species %in% sp,]
  sel.sour <- sel.sour[!duplicated(sel.sour[,"species"]),]
  sour.df <- sel.sour[,c("count", "species")]
  
  #if the genus doesn't have a hit, get rid of NAs
  if (any(!sp %in% sour.df$species)) {
    add <- data.frame("count" = 0,
                      species = sp[!sp %in% sour.df$species])
    sour.fin <- rbind(sour.df, add)
    sour.fin$method <- "sourmash_k51"
  } else {
    sour.fin <- sour.df
    sour.fin$method <- "sourmash_k51"
  }
  
  #calculate observed and real fractions
  
  total_n <- setNames(c(98195609, 98223609, 98253929, 98283929),
                      c("results/sourmash/sr/lca-class/Scaled_2000_mixed_sample1_1000_R1.csv",
                        "results/sourmash/sr/lca-class/Scaled_2000_mixed_sample1_5000_R1.csv",
                        "results/sourmash/sr/lca-class/Scaled_2000_mixed_sample1_10000_R1.csv",
                        "results/sourmash/sr/lca-class/Scaled_2000_mixed_sample1_15000_R1.csv"))
  #calculate observed fraction
  sour.fin$o_fraction <- sour.fin$count/total_n[i]
  
  #calculate real fraction
  r_fraction <- unlist(strsplit(i, split = "_"))[5] 
  
  #divide the total n reads by the fraction
  sour.fin$r_fraction <- as.numeric(r_fraction)/total_n[i]
  
  #add sample names
  sour.fin$sample <- unlist(strsplit(i, split = "_"))[5]

  #have the final table for sourmash
  fin.s <- rbind(fin.s, sour.fin)
}

#kraken2#

fin.k <- data.frame()
for (j in snakemake@input[["kraken2"]]){
  print(j)
  
  #read and match the genera
  kra <- read.table(j, header = F, sep = "\t", strip.white = T)
  kra.s <- kra[kra$V4 == "S",]
  sel.kra <- kra.s[kra.s$V6 %in% sp,] #v2 is the count
  kra.fin <- sel.kra[,c("V2", "V6")]
  colnames(kra.fin) <- c("count", "species")
  kra.fin$method <- "kraken2"
  
  #all inserted genera have hits for kraken2, no need to take care of NAs.
  
  #BU KISMI SOR!
  total_n <- setNames(c(98195609, 98223609, 98253929, 98283929),
                      c("results/kraken2/sr/sb/evol1_Sample1_fraction1000",
                        "results/kraken2/sr/sb/evol1_Sample1_fraction5000",
                        "results/kraken2/sr/sb/evol1_Sample1_fraction10000",
                        "results/kraken2/sr/sb/evol1_Sample1_fraction15000"))
  
  #calculate the observed fraction
  kra.fin$o_fraction <- kra.fin$count/total_n[j]
  
  #calculate the real fraction
  tmp <- unlist(strsplit(j, split = "_"))[3] 
  r_fraction <- unlist(strsplit(tmp, split = "fraction"))[2]
  
  #divide the total n reads by the fraction
  kra.fin$r_fraction <- as.numeric(r_fraction)/total_n[j] #divide the total n reads by the fraction
  
  #add sample names
  kra.fin$sample <- unlist(strsplit(j, split = "_"))[3]
  
  #have the final table for kraken2
  fin.k <- rbind(fin.k, kra.fin)
}  

#bracken#

fin.b <- data.frame()
for (k in snakemake@input[["bracken"]]){
  #read and match the genera
  bra <- read.table(k, header = T, sep = "\t", strip.white = T)
  sel.bra <- bra[bra[,"name"] %in% sp,] 
  bra.fin <- sel.bra[,c("name", "new_est_reads")]
  colnames(bra.fin) <- c("species", "count")
  bra.fin$method <- "bracken"
  
  #all inserted genera have hits for bracken, no need to take care of NAs.
  
  #BU KISMI SOR!
  total_n <- setNames(c(98195609, 98223609, 98253929, 98283929),
                      c("results/bracken/sr/sb/evol1_Sample1_fraction1000.bracken",
                        "results/bracken/sr/sb/evol1_Sample1_fraction5000.bracken",
                        "results/bracken/sr/sb/evol1_Sample1_fraction10000.bracken",
                        "results/bracken/sr/sb/evol1_Sample1_fraction15000.bracken"))
  
  #calculate the observed fraction
  bra.fin$o_fraction <- bra.fin$count/total_n[k]
  
  #calculate the real fraction
  tmp <- unlist(strsplit(k, split = "_"))[3] 
  r_fraction <- unlist(strsplit(tmp, split = "fraction"))[2]
  r_fraction <- gsub( ".bracken", "", r_fraction)
  
  #divide the total n reads by the fraction
  bra.fin$r_fraction <- as.numeric(r_fraction)/total_n[k] #divide the total n reads by the fraction
  
  #add sample names
  bra.fin$sample <- unlist(strsplit(k, split = "_"))[3]
  
  #have the final table for kraken2
  fin.b <- rbind(fin.b, bra.fin)
  fin.b$sample <- gsub( ".bracken", "", fin.b$sample)
}  
print(fin.b)

#final table
fin <- rbind(fin.s, fin.k)
fin <- rbind(fin, fin.b)

#modify sample names
fin$sample[1:24] <- paste0("fraction", fin$sample)
fin$sample_f <- factor(fin$sample, levels=c("fraction1000", "fraction5000", "fraction10000", "fraction15000"))

print(fin)
#scatter plot 
p <- ggplot(fin, aes(x=r_fraction, y=o_fraction, shape=species, color=species))+
  geom_point()+
  facet_grid(method~sample_f) + coord_fixed() + geom_abline(linetype="dashed")+
  scale_shape_manual(values = seq(0,6))+
  theme_classic() + annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)

#save the pdf file containing the scatter plot
ggsave(p, filename = paste0("results/final_abundance/scatter_plot/sr/sr_final_abundance_all_samples_coord_fixed.pdf"),
       width = 11,
       height = 8.5,)

#output the table
write.csv(fin, "results/final_abundance/scatter_plot/sr/sr_final_abundance_all_samples.csv", row.names = F)
