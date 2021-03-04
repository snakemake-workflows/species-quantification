library(ggplot2)

#this script creates an abundance plot for all samples and compares existing methods.

#the genera to look for:
#discarded paracoccus because kraken2 estimates absurdly high and sourmash doesnt identfy it at all.

gen <- c("Roseomonas", "Nocardia", "Anabaena", "Rhodococcus")

#sourmash#

fin.s <- data.frame()
for (i in snakemake@input[["sourmash"]]){
  
  #read and match the genera
  sour <- read.csv(i)
  sel.sour <- sour[sour$genus %in% gen,]
  sel.sour <- sel.sour[!duplicated(sel.sour[,"genus"]),]
  sour.df <- sel.sour[,c("count", "genus")]
  
  #if the genus doesn't have a hit, get rid of NAs
  if (any(!gen %in% sour.df$genus)) {
    add <- data.frame("count" = 0,
                      genus = gen[!gen %in% sour.df$genus])
    sour.fin <- rbind(sour.df, add)
    sour.fin$method <- "sourmash"
  } else {
    sour.fin <- sour.df
    sour.fin$method <- "sourmash"
  }
  
  #calculate observed and real fractions
  
  #number of reads, mixed_10000_R1: 98236994
  #number of reads, mixed_1000_R1: 98191994
  
  #calculate observed fraction
  sour.fin$o_fraction <- sour.fin$count/98000000
  
  #calculate real fraction
  r_fraction <- unlist(strsplit(i, split = "_"))[6] 
 
  #divide the total n reads by the fraction
  sour.fin$r_fraction <- as.numeric(r_fraction)/98000000 
  
  #add sample names
  sour.fin$sample <- unlist(strsplit(i, split = "_"))[6]
  
  #have the final table for sourmash
  fin.s <- rbind(fin.s, sour.fin)
}

#kraken2#

fin.k <- data.frame()
for (j in snakemake@input[["kraken2"]]){

  #read and match the genera
  kra <- read.table(j, header = F, sep = "\t")
  kra.g <- kra[kra$V4 == "G",]
  kra.g$V6 <- gsub( " ", "", kra.g$V6) 
  sel.kra <- kra.g[kra.g$V6 %in% gen,] #v2 is the count
  kra.fin <- sel.kra[,c("V2", "V6")]
  colnames(kra.fin) <- c("count", "genus")
  kra.fin$method <- "kraken2"
  
  #all inserted genera have hits for kraken2, no need to take care of NAs.
  
  #calculate the observed fraction
  kra.fin$o_fraction <- kra.fin$count/98000000
  
  #calculate the real fraction
  tmp <- unlist(strsplit(j, split = "_"))[3] 
  r_fraction <- unlist(strsplit(tmp, split = "fraction"))[2]
  
  #divide the total n reads by the fraction
  kra.fin$r_fraction <- as.numeric(r_fraction)/98000000 #divide the total n reads by the fraction

  #add sample names
  kra.fin$sample <- unlist(strsplit(j, split = "_"))[3]
  
  #have the final table for kraken2
  fin.k <- rbind(fin.k, kra.fin)
}  

#final table
fin <- rbind(fin.s, fin.k)

#modify sample names
fin$sample[1:16] <- paste0("fraction", fin$sample)
fin$sample_f <- factor(fin$sample, levels=c("fraction1000", "fraction5000", "fraction10000", "fraction15000"))

#scatter plot 
p <- ggplot(fin, aes(x=r_fraction, y=o_fraction, shape=genus, color=genus))+
  geom_point()+
  theme_minimal() + facet_grid(method~sample_f) + coord_fixed() + geom_abline(linetype="dashed")

#save the pdf file containing the scatter plot
ggsave(p, filename = paste0("results/final_abundance/scatter_plot/sr/final_abundance_all_samples_coord_fixed.pdf"),
       width = 11,
       height = 8.5,)


