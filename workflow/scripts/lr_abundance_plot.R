library(ggplot2)

#this script creates an abundance plot for all samples and compares 3 quantification methods; sourmash, kraken2, bracken.
#the species to look for:

# sp <- c("Yersinia pestis",
#         "Myxococcus xanthus",
#         "Deinococcus radiodurans",
#         "Salmonella enterica",
#         "Pseudomonas aeruginosa",
#         "Roseomonas mucosa")

sp <- snakemake@params[["species"]]
number_of_reads <- snakemake@params[["total_reads"]]

#sourmash#

fin.s <- data.frame()
for (i in snakemake@input[["sourmash"]]){
  
  #read and match the genera
  sour <- read.csv(i)
  substr(sour$species,start = 4, stop = 1000000L)
  sel.sour <- sour[sour$species %in% sp,]
  sel.sour <- sel.sour[!duplicated(sel.sour[,"species"]),]
  sour.df <- sel.sour[,c("count", "species")]
  
  #if the genus doesn't have a hit, get rid of NAs
  if (any(!sp %in% sour.df$species)) {
    add <- data.frame("count" = 0,
                      species = sp[!sp %in% sour.df$species])
    sour.fin <- rbind(sour.df, add)
    sour.fin$method <- "sourmash_k21"
  } else {
    sour.fin <- sour.df
    sour.fin$method <- "sourmash_k21"
  }
  
  #match the total number of reads with samples
  total_n <- setNames(number_of_reads, snakemake@input[["sourmash"]]) #suspicous, if matched or not?
  
  #calculate observed fraction
  sour.fin$o_fraction <- sour.fin$count/total_n[i]

    #calculate real fraction
  r_fraction <- unlist(strsplit(i, split = "_"))[5] 
  
  #divide the total n reads by the fraction
  sour.fin$r_fraction <- as.numeric(r_fraction)/total_n[i]
  
  #add sample names
  sour.fin$sample <- unlist(strsplit(i, split = "_"))[5]
  sour.fin$sample <- paste0("fraction", sour.fin$sample)
  
  #have the final table for sourmash
  fin.s <- rbind(fin.s, sour.fin)
}

#kraken2#

fin.k <- data.frame()
for (j in snakemake@input[["kraken2"]]){
 
   #read and match the genera
  kra <- read.table(j, header = F, sep = "\t", strip.white = T)
  kra.s <- kra[kra$V4 == "S",]
  
  sel.kra <- kra.s[kra.s$V6 %in% sp,]#v2 is the count
  kra.df <- sel.kra[,c("V2", "V6")]
  colnames(kra.df) <- c("count", "species")

  #if the genus doesn't have a hit, get rid of NAs
  if (any(!sp %in% kra.df$species)) {
    add <- data.frame("count" = 0,
                        species = sp[!sp %in% kra.df$species])
    kra.fin <- rbind(kra.df, add)
    kra.fin$method <- "kraken2"
  } else {
    kra.fin <- kra.df
    kra.fin$method <- "kraken2"
  }
  
  #match the total number of reads with samples
  total_n <- setNames(number_of_reads, snakemake@input[["kraken2"]])
  
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
  bra.df <- sel.bra[,c("name", "new_est_reads")]
  colnames(bra.df) <- c("species", "count")
  
  #if the genus doesn't have a hit, get rid of NAs
  if (any(!sp %in% bra.df$species)) {
    add <- data.frame("count" = 0,
                      species = sp[!sp %in% bra.df$species])
    bra.fin <- rbind(bra.df, add)
    bra.fin$method <- "bracken"
  } else {
    bra.fin <- bra.df
    bra.fin$method <- "bracken"
  }
  
  #match the total number of reads with samples
  total_n <- setNames(number_of_reads, snakemake@input[["bracken"]])
  
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

#final table
fin <- rbind(fin.s, fin.k)
fin <- rbind(fin, fin.b)

#modify sample names
fin$sample_f <- factor(fin$sample, levels=c("fraction1000", "fraction5000", "fraction10000", "fraction15000"))

print("checkpoint")
#scatter plot 
p <- ggplot(fin, aes(x=r_fraction, y=o_fraction, shape=species, color=species))+
  geom_point()+
  facet_grid(method~sample_f) + coord_fixed() + geom_abline(linetype="dashed")+
  scale_shape_manual(values = seq(0,6))+
  theme_classic() + annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)

#save the pdf file containing the scatter plot
ggsave(p, filename = paste0("results/final_abundance/scatter_plot/lr/lr_final_abundance_all_samples_coord_fixed.pdf"),
       width = 11,
       height = 8.5,)

#output the table
write.csv(fin, "results/final_abundance/scatter_plot/lr/lr_final_abundance_all_samples.csv", row.names = F)
