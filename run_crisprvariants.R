library("CrispRVariants")
library("GenomicRanges")
library(rtracklayer)
library(reshape2)
library(stringr)
## Analyze vector-free BAM files with crisprvariants R package and add previously quantified vector insertions

# Focusing range
# Select coordinates to match the figures by Ibon
blat.strand <- "-"
blat.T.name <- "NC_000068.7"
blat.T.start <- 134548201
blat.T.end <- 134548300

# Focusing range of 100 bp
gd = GRanges(blat.T.name, blat.strand, ranges = IRanges(blat.T.start, blat.T.end))

# Run this command in linux enviroment with samtools to get the sequence
command=sprintf("samtools faidx blastdb/GCF_000001635.26_GRCm38.p6_genomic.fna %s:%s-%s",
                seqnames(gd)[1], start(gd)[1], end(gd)[1])

# Focusing range sequences
# 100bp
focusing.sequence.char <- "TGCTGACTCTCTGTCCTAAAACAGAAGTTGACAGATCGATATCAGCAACGTTGCGAAGCATCCGTGGATAGAGCTTCCATctggaatttaaaataatatt"
focusing.sequence.char=Biostrings::reverseComplement(Biostrings::DNAString(focusing.sequence.char))

# Cut point in the focusing range (5'->3')
target.location <- 22

# Flag to determine how chimeric reads are treated. One of "ignore", "exclude", "merge" or "count".
treat.chimeras <- "count"

# Open sample tables
md <- read.table("samples.csv", sep = "\t", header = F)
md.mef <- md[md$V4=="MEFs",]
md.sacas9 <- md[md$V4=="Sacas9",]
md.tbg <- md[md$V4=="TBG",]
md.aat <- md[md$V4=="AAT",]

# Sample names of this analysis
dataname.char.mef <- md.mef$V3
dataname.char.sacas9 <- md.sacas9$V3
dataname.char.tbg <- md.tbg$V3
dataname.char.aat <- md.aat$V3

# Place of bam files
bam.file.path.mef <- file.path(md.mef$V6, md.mef$V5)
bam.file.path.sacas9 <- file.path(md.sacas9$V13, md.sacas9$V10)
bam.file.path.tbg <- file.path(md.tbg$V13, md.tbg$V10)
bam.file.path.aat <- file.path(md.aat$V13, md.aat$V10)

# Open viral read count table (from multiqc on vector-only fastq files)
vcounts <- read.table("multiqc_fastqc.txt", sep = "\t", header=T)
vcounts$Sample <- gsub("_virus_mapped_1", "", vcounts$Sample)

# Merge tables and fill missing values of samples without viral mapping
md <- merge(md, vcounts, by.x = "V1", by.y = "Sample", all.x = T)
rownames(md) <- md$V3
# Fill empty cells with zeros (no virus found)
md[["Total.Sequences"]][is.na(md[["Total.Sequences"]])] <- 0

cols <- c("Deletion" = "#fc8d62", "Insertion" = "#8da0cb", "Vector" = "#b3b3b3")

### FUNCTIONS ##################################################################

get_deletions <- function(var_counts, byType){
  del <- cbind.data.frame(var_counts, byType)
  del <- del[del$byType == "deletion",]
  del[c('Loc', 'Size')] <- str_split_fixed(rownames(del), ':', 2)
  del$Size <- mapply(gsub, "D", "", del$Size)
  del$byType <- NULL
  del$Loc <- NULL
  del <- aggregate(.~Size,data=del,FUN=sum)
  del <- melt(del, id = c("Size"))
  del$Size = as.numeric(as.character(del$Size))
  del <- del[rep(row.names(del), del$value), 1:ncol(del)]
  del$Type <- "Deletion"
  return(del)
}

get_insertions <- function(var_counts, byType){
  ins <- cbind.data.frame(var_counts, byType)
  ins <- ins[ins$byType == "insertion",]
  ins[c('Loc', 'Size')] <- str_split_fixed(rownames(ins), ':', 2)
  ins$Size <- mapply(gsub, "I", "", ins$Size)
  ins$byType <- NULL
  ins$Loc <- NULL
  ins <- aggregate(.~Size,data=ins,FUN=sum)
  ins <- melt(ins, id = c("Size"))
  ins$Size = as.numeric(as.character(ins$Size))
  ins <- ins[rep(row.names(ins), ins$value), 1:ncol(ins)]
  ins$Type <- "Insertion"
  return(ins)
}

get_vector <- function(sampletable) {
  dir <- "C:/Users/Usuario/Desktop/CIMA/PH1_Torella/crisprvariants/"
  columns = c("Size","variable","value","Type") 
  vdf = data.frame(matrix(nrow = 0, ncol = length(columns))) 
  colnames(vdf) = columns
  for (i in 1:length(rownames(sampletable)) ){
    vfiledir <- paste(dir,"v_pipeline_insertions/",sampletable$V1[i], "_virus.tsv", sep = "")
    vfile <- read.table(vfiledir, sep = "\t", header=T)
    if (length(rownames(vfile)) != 0) {
      vfile$Size <- vfile$Length
      vfile$variable <- sampletable$V3[i]
      vfile$value <- 0
      vfile$Type <- "Vector"
      vfile <- vfile[ , columns]
      vdf <- rbind(vdf, vfile)
    }
  }
  return(vdf)
}


### Analyse ####################################################################
################################################################################
## MEFs
crispr.set.mef <- readsToTarget(bam.file.path.mef, gd,
                            reference = focusing.sequence.char,
                            target.loc = target.location,
                            names = dataname.char.mef,
                            chimeras = treat.chimeras,
                            collapse.pairs = TRUE,
                            chimera.to.target = 200)

eff.mef <- mutationEfficiency(crispr.set.mef)
var_counts.mef <- variantCounts(crispr.set.mef)
#sqs <- consensusSeqs(crispr.set)
write.table(var_counts.mef, "in_vitro_MEFs.results.csv", sep = "\t", quote = FALSE, row.names = TRUE)

pdf("in_vitro_MEFs.pdf",
    width = 12, height = 8, bg = "white")
plotVariants(crispr.set.mef, col.wdth.ratio = c(4,1),
                  plotFreqHeatmap.args = list(top.n = 40, type="proportions"),
                  plotAlignments.args = list(pam.start = c(14,90),
                                             pam.end = c(19,95),
                                             target.loc = c(22, 86),
                                             max.insertion.size = 100,
                                             top.n = 40,
                                             legend.cols = 4,
                                             plot.text.size = 1,
                                             guide.loc = IRanges::IRanges(c(20, 68),c(41, 89))))
dev.off() 

p <- plotVariants(crispr.set2,
                  plotAlignments.args = list(pam.start = c(64,140),
                                             pam.end = c(69,145),
                                             target.loc = c(72, 136),
                                             max.insertion.size = 100,
                                             guide.loc = IRanges::IRanges(c(70, 118),c(91, 139))))

# Variant type plot
byType.mef <- crispr.set.mef$classifyVariantsByType()
unique(byType.mef)
byType.mef <- unlist(byType.mef)

var_counts_plot.mef <- apply(var_counts.mef,2,function(x){x/sum(x)*100})
var_counts_plot.mef <- cbind.data.frame(var_counts_plot.mef, byType.mef)
var_counts_plot.mef$byType.mef[var_counts_plot.mef$byType.mef == "SNV"] <- "no variant"
var_counts_plot.mef <- aggregate(.~byType.mef,data=var_counts_plot.mef,FUN=sum)
var_counts_plot.mef <- melt(var_counts_plot.mef, id = c("byType.mef"))
var_counts_plot.mef$byType.mef <- factor(var_counts_plot.mef$byType.mef, levels = unique(byType.mef)) 

pdf("in_vitro_MEFs.types.pdf",
    width = 6, height = 2.5, bg = "white")
ggplot(var_counts_plot.mef, aes(fill=byType.mef, x = variable, y = value)) + 
  geom_bar(position='stack', stat='identity', alpha=0.75) + coord_flip() +
  theme_linedraw() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
  scale_fill_brewer(palette="Set2") +
  labs(x = "Sample", y = "% mapping reads", fill = "Type")
dev.off()

var_means.mef <-  merge(md[c("V3","V8")], var_counts_plot.mef, by.x = "V3", by.y = "variable")
var_means.mef$V3 <- NULL
var_means.mef$value <- as.numeric(as.character(var_means.mef$value))
var_means.mef <- aggregate(.~V8+byType.mef,data=var_means.mef,FUN=mean)

pdf("in_vitro_MEFs.means.pdf",
    width = 6, height = 2.5, bg = "white")
ggplot(var_means.mef, aes(fill=byType.mef, x = V8, y = value)) + 
  geom_bar(position='stack', stat='identity', alpha=0.75) + coord_flip() +
  theme_linedraw() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
  scale_fill_brewer(palette="Set2") +
  scale_x_discrete(limits = rev(md.mef$V8)) +
  labs(x = "Group", y = "% mapping reads", fill = "Type")
dev.off()

del.mef <- get_deletions(var_counts.mef, byType.mef)
ins.mef <- get_insertions(var_counts.mef, byType.mef)
indels.mef <- rbind(del.mef, ins.mef)

pdf("in_vitro_MEFs.deletions.pdf",
    width = 6, height = 3.5, bg = "white")
ggplot(indels.mef, aes(color=NULL, fill=Type , x = Size)) + 
  geom_histogram(binwidth=5, alpha=0.7, position = "identity") +
  scale_fill_manual(values=cols) + theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  labs(x = "Indel Size", y = "Indel count") +
  facet_wrap(~variable,  ncol=2, strip.position = "top") 
dev.off()

indels.mef.ut <- indels.mef[indels.mef$variable=="UT-MEF",]
indels.mef.g1 <- indels.mef[indels.mef$variable=="SaCas9g1-MEF",]
indels.mef.g2 <- indels.mef[indels.mef$variable=="SaCas9g2-MEF",]
indels.mef.g1g2 <- indels.mef[indels.mef$variable=="SaCas9g1+g2-MEF",]

pdf("SaCas9g1-MEF.deletions.pdf",
    width = 4, height = 2.5, bg = "white")
ggplot(indels.mef.g1, aes(color=NULL, fill=Type , x = Size)) + xlim(-5, 100) +
  geom_histogram(aes(y = after_stat(count / sum(count))), binwidth=5, alpha=0.7, position = "identity") +
  scale_fill_manual(values=cols) + theme_linedraw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        #legend.position = c(0.7, 0.7),
        plot.title = element_text(vjust = 0, hjust=0.5, size = 10)) +
  labs(x = "Indel Size", y = "Indel fraction (%)") +
  scale_y_continuous(labels = scales::percent) +
  ggtitle("SaCas9g1-MEF")
dev.off()

pdf("SaCas9g2-MEF.deletions.pdf",
    width = 4, height = 2.5, bg = "white")
ggplot(indels.mef.g2, aes(color=NULL, fill=Type , x = Size)) + xlim(-5, 100) +
  geom_histogram(aes(y = after_stat(count / sum(count))), binwidth=5, alpha=0.7, position = "identity") +
  scale_fill_manual(values=cols) + theme_linedraw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        #legend.position = c(0.7, 0.7),
        plot.title = element_text(vjust = 0, hjust=0.5, size = 10)) +
  labs(x = "Indel Size", y = "Indel fraction (%)") +
  scale_y_continuous(labels = scales::percent) +
  ggtitle("SaCas9g2-MEF")
dev.off()

pdf("SaCas9g1g2-MEF.deletions.pdf",
    width = 4, height = 2.5, bg = "white")
ggplot(indels.mef.g1g2, aes(color=NULL, fill=Type , x = Size)) + xlim(-5, 100) +
  geom_histogram(aes(y = after_stat(count / sum(count))), binwidth=5, alpha=0.7, position = "identity") +
  scale_fill_manual(values=cols) + theme_linedraw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        #legend.position = c(0.7, 0.7),
        plot.title = element_text(vjust = 0, hjust=0.5, size = 10)) +
  labs(x = "Indel Size", y = "Indel fraction (%)") +
  scale_y_continuous(labels = scales::percent) +
  ggtitle("SaCas9g1+g2-MEF")
dev.off()


################################################################################
## Sacas9 samples
crispr.set.sacas9 <- readsToTarget(bam.file.path.sacas9, gd,
                                  reference = focusing.sequence.char,
                                  target.loc = target.location,
                                  names = dataname.char.sacas9,
                                  chimeras = treat.chimeras,
                                  collapse.pairs = TRUE,
                                  chimera.to.target = 200)

eff.sacas9 <- mutationEfficiency(crispr.set.sacas9)
var_counts.sacas9 <- variantCounts(crispr.set.sacas9)
# Add vector integrations
vcounts.sacas9 <- t(md[colnames(var_counts.sacas9),]["Total.Sequences"])
rownames(vcounts.sacas9) <- c("Vector")
var_counts.sacas9 <- rbind(var_counts.sacas9, vcounts.sacas9)
# Save variant table
write.table(var_counts.sacas9, "Sacas9.results.vector.csv", sep = "\t", quote = FALSE, row.names = TRUE)

pdf("Sacas9.pdf",
    width = 15, height = 6, bg = "white")
plotVariants(crispr.set.sacas9, col.wdth.ratio = c(4,1),
                  plotFreqHeatmap.args = list(top.n = 25, type="proportions"),
                  plotAlignments.args = list(pam.start = c(14,90),
                                             pam.end = c(19,95),
                                             target.loc = c(22, 86),
                                             max.insertion.size = 100,
                                             top.n = 25,
                                             legend.cols = 4,
                                             plot.text.size = 1,
                                             guide.loc = IRanges::IRanges(c(20, 68),c(41, 89))))
dev.off() 

# Plot top 3 variants: perfect deletion plus next two
pdf("Sacas9_top3_width8.pdf",
    width = 8, height = 1.5, bg = "white")
plotAlignments(crispr.set.sacas9, pam.start = c(14,90),
                                  pam.end = c(19,95),
                                  target.loc = c(22, 86),
                                  max.insertion.size = 100,
                                  top.n = 5,
                                  legend.cols = 4,
                                  plot.text.size = 1.5,
                                  add.other=FALSE,
                                  guide.loc = IRanges::IRanges(c(20, 68),c(41, 89)))
dev.off() 

# Variant type plot
byType.sacas9 <- crispr.set.sacas9$classifyVariantsByType()
unique(byType.sacas9)
byType.sacas9 <- unlist(byType.sacas9)
byType.sacas9["Vector"] <- "vector"
byType.sacas9["Other"] <- "other"
var_counts.sacas9["Vector",]

var_counts_plot.sacas9 <- apply(var_counts.sacas9,2,function(x){x/sum(x)*100})
var_counts_plot.sacas9 <- cbind.data.frame(var_counts_plot.sacas9, byType.sacas9)
var_counts_plot.sacas9$byType.sacas9[var_counts_plot.sacas9$byType.sacas9 == "SNV"] <- "no variant"
var_counts_plot.sacas9 <- aggregate(.~byType.sacas9,data=var_counts_plot.sacas9,FUN=sum)
var_counts_plot.sacas9 <- melt(var_counts_plot.sacas9, id = c("byType.sacas9"))
var_counts_plot.sacas9$byType.sacas9 <- factor(var_counts_plot.sacas9$byType.sacas9, levels = unique(byType.sacas9)) 

pdf("Sacas9.types.virus.pdf",
    width = 6, height = 3, bg = "white")
ggplot(var_counts_plot.sacas9, aes(fill=byType.sacas9, x = variable, y = value)) + 
  geom_bar(position='stack', stat='identity', alpha=0.75) + coord_flip() +
  theme_linedraw() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
  scale_fill_brewer(palette="Set2") +
  labs(x = "Sample", y = "% mapping reads", fill = "Type")
dev.off()

var_means.sacas9 <-  merge(md[c("V3","V8")], var_counts_plot.sacas9, by.x = "V3", by.y = "variable")
var_means.sacas9$V3 <- NULL
var_means.sacas9$value <- as.numeric(as.character(var_means.sacas9$value))
var_means.sacas9 <- aggregate(.~V8+byType.sacas9,data=var_means.sacas9,FUN=mean)

del.sacas9 <- get_deletions(var_counts.sacas9, byType.sacas9)
ins.sacas9 <- get_insertions(var_counts.sacas9, byType.sacas9)
vec.sacas9 <- get_vector(md.sacas9)
  
indels.sacas9 <- rbind(del.sacas9, ins.sacas9, vec.sacas9)

pdf("Sacas9.deletions.pdf",
    width = 4, height = 2.5, bg = "white")
ggplot(indels.sacas9, aes(color=NULL, fill=Type , x = Size)) + xlim(-5, 100) +
  geom_histogram(aes(y = after_stat(count / sum(count))), binwidth=5, alpha=0.7, position = "identity") +
  scale_fill_manual(values=cols) + theme_linedraw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        #legend.position = c(0.7, 0.7),
        plot.title = element_text(vjust = 0, hjust=0.5, size = 10)) +
  labs(x = "Indel Size", y = "Indel fraction (%)") +
  scale_y_continuous(labels = scales::percent) +
  ggtitle("TBG−SaCas9g1+g2")
dev.off()

################################################################################
## TBG + PBS samples
crispr.set.tbg <- readsToTarget(bam.file.path.tbg, gd,
                                   reference = focusing.sequence.char,
                                   target.loc = target.location,
                                   names = dataname.char.tbg,
                                   chimeras = treat.chimeras,
                                   collapse.pairs = TRUE,
                                   chimera.to.target = 200)

crispr.set.tbgonly <- readsToTarget(bam.file.path.tbg[9:13], gd,
                                reference = focusing.sequence.char,
                                target.loc = target.location,
                                names = dataname.char.tbg[9:13],
                                chimeras = treat.chimeras,
                                collapse.pairs = TRUE,
                                chimera.to.target = 200)

eff.tbg <- mutationEfficiency(crispr.set.tbg)
var_counts.tbg <- variantCounts(crispr.set.tbg)
var_counts.tbgonly <- variantCounts(crispr.set.tbgonly)


# Add vector integrations
vcounts.tbg <- t(md[colnames(var_counts.tbg),]["Total.Sequences"])
rownames(vcounts.tbg) <- c("Vector")
var_counts.tbg <- rbind(var_counts.tbg, vcounts.tbg)
vcounts.tbgonly <- t(md[colnames(var_counts.tbgonly),]["Total.Sequences"])
rownames(vcounts.tbgonly) <- c("Vector")
var_counts.tbgonly <- rbind(var_counts.tbgonly, vcounts.tbgonly)
# Save variant table
write.table(var_counts.tbg, "TBG.results.vector.csv", sep = "\t", quote = FALSE, row.names = TRUE)

png("TBG.png",
    width = 14, height = 8, units = "in", bg = "white", res=300)
plotVariants(crispr.set.tbg, col.wdth.ratio = c(1.75,1),
                  plotFreqHeatmap.args = list(top.n = 40, type="proportions"),
                  plotAlignments.args = list(pam.start = c(14,90),
                                             pam.end = c(19,95),
                                             target.loc = c(22, 86),
                                             max.insertion.size = 100,
                                             top.n = 40,
                                             legend.cols = 4,
                                             plot.text.size = 2,
                                             guide.loc = IRanges::IRanges(c(20, 68),c(41, 89))))
dev.off() 

# Plot top 20 variants without SNVs and only treatment samples
pdf("TBG_top20_width8.pdf",
    width = 8, height = 3, bg = "white")
plotAlignments(crispr.set.tbgonly,pam.start = c(14,90),
               pam.end = c(19,95),
               target.loc = c(22, 86),
               max.insertion.size = 100,
               min.freq = 0.17,  # Limit to keep SNVs out of the plot
               top.n = 50,
               legend.cols = 3,
               plot.text.size = 1.5,
               add.other=FALSE,
               style="all",
               guide.loc = IRanges::IRanges(c(20, 68),c(41, 89)))
dev.off()

# Plot top 20 variants without SNVs and only treatment samples WITH heatmap
pdf("TBG_top20_with_heatmap.pdf",
    width = 12, height = 6, bg = "white")
plotVariants(crispr.set.tbgonly, col.wdth.ratio = c(3,1),
             plotFreqHeatmap.args = list(top.n = 40, type="proportions"),
             plotAlignments.args = list(pam.start = c(14,90),
                                        pam.end = c(19,95),
                                        target.loc = c(22, 86),
                                        max.insertion.size = 100,
                                        min.freq = 0.17,  # Limit to keep SNVs out of the plot
                                        top.n = 50,
                                        legend.cols = 4,
                                        plot.text.size = 2,
                                        guide.loc = IRanges::IRanges(c(20, 68),c(41, 89))))
dev.off()

# Variant type plot
byType.tbg <- crispr.set.tbg$classifyVariantsByType()
unique(byType.tbg)
byType.tbg <- unlist(byType.tbg)
byType.tbg["Vector"] <- "vector"
byType.tbg["Other"] <- "other"

byType.tbgonly <- crispr.set.tbgonly$classifyVariantsByType()
unique(byType.tbgonly)
byType.tbgonly <- unlist(byType.tbgonly)
byType.tbgonly["Vector"] <- "vector"
byType.tbgonly["Other"] <- "other"


var_counts_plot.tbg <- apply(var_counts.tbg,2,function(x){x/sum(x)*100})
var_counts_plot.tbg <- cbind.data.frame(var_counts_plot.tbg, byType.tbg)
var_counts_plot.tbg$byType.tbg[var_counts_plot.tbg$byType.tbg == "SNV"] <- "no variant"
var_counts_plot.tbg <- aggregate(.~byType.tbg,data=var_counts_plot.tbg,FUN=sum)
var_counts_plot.tbg <- melt(var_counts_plot.tbg, id = c("byType.tbg"))
var_counts_plot.tbg$byType.tbg <- factor(var_counts_plot.tbg$byType.tbg, levels = unique(byType.tbg)) 


pdf("TBG.types.vector.pdf",
    width = 6, height = 3, bg = "white")
ggplot(var_counts_plot.tbg, aes(fill=byType.tbg, x = variable, y = value)) + 
  geom_bar(position='stack', stat='identity', alpha=0.75) + coord_flip() +
  theme_linedraw() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
  scale_fill_brewer(palette="Set2") +
  labs(x = "Sample", y = "% mapping reads", fill = "Type")
dev.off()

var_means.tbg <-  merge(md[c("V3","V8")], var_counts_plot.tbg, by.x = "V3", by.y = "variable")
var_means.tbg$V3 <- NULL
var_means.tbg$value <- as.numeric(as.character(var_means.tbg$value))
var_means.tbg <- aggregate(.~V8+byType.tbg,data=var_means.tbg,FUN=mean)

pdf("TBG.means.vector.pdf",
    width = 6, height = 3, bg = "white")
ggplot(var_means.tbg, aes(fill=byType.tbg, x = V8, y = value)) + 
  geom_bar(position='stack', stat='identity', alpha=0.75) + coord_flip() +
  theme_linedraw() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
  scale_fill_brewer(palette="Set2") +
  scale_x_discrete(limits = rev(unique(md.tbg$V8))) +
  labs(x = "Group", y = "% mapping reads", fill = "Type")
dev.off()

del.tbg <- get_deletions(var_counts.tbgonly, byType.tbgonly)
ins.tbg <- get_insertions(var_counts.tbgonly, byType.tbgonly)
vec.tbg <- get_vector(md.tbg[9:13,])
indels.tbg <- rbind(del.tbg, ins.tbg, vec.tbg)

pdf("TBG.deletions.pdf",
    width = 4, height = 2.5, bg = "white")
ggplot(indels.tbg, aes(color=NULL, fill=Type , x = Size)) + xlim(-5, 100) +
  geom_histogram(aes(y = after_stat(count / sum(count))), binwidth=5, alpha=0.7, position = "identity") +
  scale_fill_manual(values=cols) + theme_linedraw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        #legend.position = c(0.7, 0.7),
        plot.title = element_text(vjust = 0, hjust=0.5, size = 10)) +
  labs(x = "Indel Size", y = "Indel fraction (%)") +
  scale_y_continuous(labels = scales::percent) +
  ggtitle("TBG-D10Ag1+g2")
dev.off()

################################################################################
## Join Sacas9 + TBG + Saline
colnames(var_counts_plot.sacas9) <- colnames(var_counts_plot.tbg)
var_counts_plot.tbg_sacas9 <- rbind(var_counts_plot.tbg, var_counts_plot.sacas9)

pdf("TBG_Sacas9.types.vector.pdf",
    width = 6, height = 3, bg = "white")
ggplot(var_counts_plot.tbg_sacas9, aes(fill=byType.tbg, x = variable, y = value)) + 
  geom_bar(position='stack', stat='identity', alpha=0.75) + coord_flip() +
  theme_linedraw() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
  scale_fill_brewer(palette="Set2") +
  labs(x = "Sample", y = "% mapping reads", fill = "Type")
dev.off()

var_means.tbg_sacas9 <-  merge(md[c("V3","V8")], var_counts_plot.tbg_sacas9, by.x = "V3", by.y = "variable")
var_means.tbg_sacas9$V3 <- NULL
var_means.tbg_sacas9$value <- as.numeric(as.character(var_means.tbg_sacas9$value))
var_means.tbg_sacas9 <- aggregate(.~V8+byType.tbg,data=var_means.tbg_sacas9,FUN=mean)

pdf("TBG_Sacas9.means.vector.pdf",
    width = 5.5, height = 2.5, bg = "white")
ggplot(var_means.tbg_sacas9, aes(fill=byType.tbg, x = V8, y = value)) + 
  geom_bar(position='stack', stat='identity', alpha=0.75) + coord_flip() +
  theme_linedraw() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
  scale_fill_brewer(palette="Set2") +
  scale_x_discrete(limits = c("TBG-SaCas9g1+g2",rev(unique(md.tbg$V8)))) +
  labs(x = "Group", y = "% mapping reads", fill = "Type")
dev.off()

# Add results from single guides sacas9 (zabaleta 2018)
var_means.zab <- read.table("Zabaletag1g2.means.vector.csv", sep = "\t")
colnames(var_means.zab) <- colnames(var_means.tbg_sacas9)
var_means.tbg_sacas9_zab <- rbind(var_means.tbg_sacas9, var_means.zab)
var_means.tbg_sacas9_zab <- var_means.tbg_sacas9_zab[var_means.tbg_sacas9_zab$V8 != "Control_g1",]
var_means.tbg_sacas9_zab <- var_means.tbg_sacas9_zab[var_means.tbg_sacas9_zab$V8 != "Control_g2",]

pdf("TBG_Sacas9_zab2018.means.vector3.pdf",
    width = 4, height = 4, bg = "white")
ggplot(var_means.tbg_sacas9_zab, aes(fill=byType.tbg, x = V8, y = value)) + 
  geom_bar(position='stack', stat='identity', alpha=0.75) + coord_flip() +
  theme_linedraw() + 
  theme(legend.position = "bottom",
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.key.size = unit(0.5, "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  guides(fill=guide_legend(nrow=4,byrow=TRUE)) +
  scale_fill_brewer(palette="Set2") +
  scale_x_discrete(limits = c("TBG-SaCas9g2","TBG-SaCas9g1","TBG-SaCas9g1+g2",rev(unique(md.tbg$V8)))) +
  labs(x = "", y = "% mapping reads", fill = "Type")
dev.off()

################################################################################
## AAT samples
crispr.set.aat <- readsToTarget(bam.file.path.aat, gd,
                                reference = focusing.sequence.char,
                                target.loc = target.location,
                                names = dataname.char.aat,
                                chimeras = treat.chimeras,
                                collapse.pairs = TRUE,
                                chimera.to.target = 200)

crispr.set.aatonly <- readsToTarget(bam.file.path.aat[5:20], gd,
                                reference = focusing.sequence.char,
                                target.loc = target.location,
                                names = dataname.char.aat[5:20],
                                chimeras = treat.chimeras,
                                collapse.pairs = TRUE,
                                chimera.to.target = 200)

eff.aat <- mutationEfficiency(crispr.set.aat)
var_counts.aat <- variantCounts(crispr.set.aat)

# Add vector integrations
vcounts.aat <- t(md[colnames(var_counts.aat),]["Total.Sequences"])
rownames(vcounts.aat) <- c("Vector")
var_counts.aat <- rbind(var_counts.aat, vcounts.aat)
# Save variant table
write.table(var_counts.aat, "AAT.results.vector.csv", sep = "\t", quote = FALSE, row.names = TRUE)

pdf("AAT.pdf",
    width = 16, height = 8, bg = "white")
plotVariants(crispr.set.aatonly, col.wdth.ratio = c(1.5,1),
                  plotFreqHeatmap.args = list(top.n = 40, plot.text.size = 2, type="proportions"),
                  plotAlignments.args = list(pam.start = c(14,90),
                                             pam.end = c(19,95),
                                             target.loc = c(22, 86),
                                             max.insertion.size = 100,
                                             top.n = 40,
                                             legend.cols = 4,
                                             plot.text.size = 2,
                                             guide.loc = IRanges::IRanges(c(20, 68),c(41, 89))))
dev.off() 

# Plot top 20 variants without SNVs and only treatment samples
pdf("AAT_top20_width6.pdf",
    width = 8, height = 5, bg = "white")
plotAlignments(crispr.set.aatonly,pam.start = c(14,90),
               pam.end = c(19,95),
               target.loc = c(22, 86),
               max.insertion.size = 100,
               min.freq = 0.15,  # Limit to keep SNVs out of the plot
               top.n = 40,
               legend.cols = 3,
               plot.text.size = 1.5,
               add.other=FALSE,
               style="all",
               guide.loc = IRanges::IRanges(c(20, 68),c(41, 89)))
dev.off()

# Variant type plot
byType.aat <- crispr.set.aat$classifyVariantsByType()
unique(byType.aat)
byType.aat <- unlist(byType.aat)
byType.aat["Vector"] <- "vector"
byType.aat["Other"] <- "other"

var_counts_plot.aat <- apply(var_counts.aat,2,function(x){x/sum(x)*100})
var_counts_plot.aat <- cbind.data.frame(var_counts_plot.aat, byType.aat)
var_counts_plot.aat$byType.aat[var_counts_plot.aat$byType.aat == "SNV"] <- "no variant"
var_counts_plot.aat <- aggregate(.~byType.aat,data=var_counts_plot.aat,FUN=sum)
var_counts_plot.aat <- melt(var_counts_plot.aat, id = c("byType.aat"))
var_counts_plot.aat$byType.aat <- factor(var_counts_plot.aat$byType.aat, levels = unique(byType.aat)) 

pdf("AAT.types.vector.pdf",
    width = 6, height = 5, bg = "white")
ggplot(var_counts_plot.aat, aes(fill=byType.aat, x = variable, y = value)) + 
  geom_bar(position='stack', stat='identity', alpha=0.75) + coord_flip() +
  theme_linedraw() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
  scale_fill_brewer(palette="Set2") +
  labs(x = "Sample", y = "% mapping reads", fill = "Type")
dev.off()

var_means.aat <-  merge(md[c("V3","V8")], var_counts_plot.aat, by.x = "V3", by.y = "variable")
var_means.aat$V3 <- NULL
var_means.aat$value <- as.numeric(as.character(var_means.aat$value))
var_means.aat <- aggregate(.~V8+byType.aat,data=var_means.aat,FUN=mean)

pdf("AAT.means.vector.pdf",
    width = 4.85, height = 2.5, bg = "white")
ggplot(var_means.aat, aes(fill=byType.aat, x = V8, y = value)) + 
  geom_bar(position='stack', stat='identity', alpha=0.75) + coord_flip() +
  theme_linedraw() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
  scale_fill_brewer(palette="Set2") +
  scale_x_discrete(limits = rev(unique(md.aat$V8))) +
  labs(x = "Group", y = "% mapping reads", fill = "Type")
dev.off()
