library("CrispRVariants")
library("GenomicRanges")
library(rtracklayer)
library(reshape2)
library(stringr)
## Analyze vector-free BAM files with crisprvariants R package and add previously quantified vector insertions

# Focusing range
blat.strand <- "-"
blat.T.name <- "NC_000068.7"

blat.T.start1 <- 134548261 - 20
blat.T.end1 <- 134548281 + 20

blat.T.start2 <- 134548212 - 20
blat.T.end2 <- 134548232 + 20

amp.start1 <- 134548178
amp.end1 <- 134548458
amp.start2 <- 134548069
amp.end2 <- 134548283

# Focusing range of 61 bp around each guide
gd1 = GRanges(blat.T.name, blat.strand, ranges = IRanges(blat.T.start1, blat.T.end1))
gd2 = GRanges(blat.T.name, blat.strand, ranges = IRanges(blat.T.start2, blat.T.end2))
amplicon1 = GRanges(blat.T.name, blat.strand, ranges = IRanges(amp.start1, amp.end1))
amplicon2 = GRanges(blat.T.name, blat.strand, ranges = IRanges(amp.start2, amp.end2))

# Run this command in linux enviroment with samtools to get the sequence
command1=sprintf("samtools faidx blastdb/GCF_000001635.26_GRCm38.p6_genomic.fna %s:%s-%s",
                seqnames(gd1)[1], start(gd1)[1], end(gd1)[1])
command2=sprintf("samtools faidx blastdb/GCF_000001635.26_GRCm38.p6_genomic.fna %s:%s-%s",
                seqnames(gd2)[1], start(gd2)[1], end(gd2)[1])
command_amp1=sprintf("samtools faidx blastdb/GCF_000001635.26_GRCm38.p6_genomic.fna %s:%s-%s",
                 seqnames(amplicon1)[1], start(amplicon1)[1], end(amplicon1)[1])
command_amp2=sprintf("samtools faidx blastdb/GCF_000001635.26_GRCm38.p6_genomic.fna %s:%s-%s",
                 seqnames(amplicon2)[1], start(amplicon2)[1], end(amplicon2)[1])

# Focusing range sequences
# Guide1
focusing.sequence.char1 <- "ATCAGCAACGTTGCGAAGCATCCGTGGATAGAGCTTCCATCTGGAATTTAAAATAATATTA"
focusing.sequence.char1=Biostrings::reverseComplement(Biostrings::DNAString(focusing.sequence.char1))
focusing.sequence.char2 <- "ATATTGGCATGCTGACTCTCTGTCCTAAAACAGAAGTTGACAGATCGATATCAGCAACGTT"
focusing.sequence.char2=Biostrings::reverseComplement(Biostrings::DNAString(focusing.sequence.char2))

# Cut point in the focusing range (5'->3')
target.location1 <- 23
target.location2 <- 38

# Flag to determine how chimeric reads are treated. One of "ignore", "exclude", "merge" or "count".
treat.chimeras <- "count"

# Open sample tables
md <- read.table("samples.csv", sep = "\t", header = F)
md.g1 <- md[md$V8=="TBG-SaCas9g1",]
md.g2 <- md[md$V8=="TBG-SaCas9g2",]
md.c1 <- md[md$V8=="Control_g1",]
md.c2 <- md[md$V8=="Control_g2",]
md.g1 <- rbind(md.g1,md.c1)
md.g2 <- rbind(md.g2,md.c2)
md.g1only <- md[md$V8=="TBG-SaCas9g1",]
md.g2only <- md[md$V8=="TBG-SaCas9g2",]

# Sample names of this analysis
dataname.char.g1 <- md.g1$V3
dataname.char.g2 <- md.g2$V3
dataname.char.g1only <- md.g1only$V3
dataname.char.g2only <- md.g2only$V3

# Place of bam files
bam.file.path.g1 <- file.path(md.g1$V13, md.g1$V10)
bam.file.path.g2 <- file.path(md.g2$V13, md.g2$V10)
bam.file.path.g1only <- file.path(md.g1only$V13, md.g1only$V10)
bam.file.path.g2only <- file.path(md.g2only$V13, md.g2only$V10)

# Open viral read count table (from multiqc)
vcounts <- read.table("multiqc_fastqc_zabaleta.txt", sep = "\t", header=T)
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


################################################################################
## Guide1 samples
crispr.set.g1 <- readsToTarget(bam.file.path.g1, gd1,
                              reference = focusing.sequence.char1,
                              target.loc = target.location1,
                              names = dataname.char.g1,
                              chimeras = treat.chimeras,
                              collapse.pairs = TRUE,
                              chimera.to.target = 200)
crispr.set.g1only <- readsToTarget(bam.file.path.g1only, gd1,
                               reference = focusing.sequence.char1,
                               target.loc = target.location1,
                               names = dataname.char.g1only,
                               chimeras = treat.chimeras,
                               collapse.pairs = TRUE,
                               chimera.to.target = 200)

eff.g1 <- mutationEfficiency(crispr.set.g1)
var_counts.g1 <- variantCounts(crispr.set.g1)
# Add vector integrations
vcounts.g1 <- t(md[colnames(var_counts.g1),]["Total.Sequences"])
rownames(vcounts.g1) <- c("Vector")
var_counts.g1 <- rbind(var_counts.g1, vcounts.g1)

eff.g1only <- mutationEfficiency(crispr.set.g1only)
var_counts.g1only <- variantCounts(crispr.set.g1only)
# Add vector integrations
vcounts.g1only <- t(md[colnames(var_counts.g1only),]["Total.Sequences"])
rownames(vcounts.g1only) <- c("Vector")
var_counts.g1only <- rbind(var_counts.g1only, vcounts.g1only)
# Save variant table
write.table(var_counts.g1, "zabaleta_g1.results.vector.csv", sep = "\t", quote = FALSE, row.names = TRUE)

pdf("zabaleta_g1.pdf",
    width = 15, height = 7, bg = "white")
plotVariants(crispr.set.g1, col.wdth.ratio = c(2,1),
                  plotFreqHeatmap.args = list(top.n = 30, type="proportions"),
                  plotAlignments.args = list(pam.start = c(15),
                                             pam.end = c(20),
                                             target.loc = c(23),
                                             max.insertion.size = 100,
                                             top.n = 30,
                                             legend.cols = 5,
                                             plot.text.size = 2,
                                             guide.loc = IRanges::IRanges(c(21),c(42))))
dev.off()

# Plot top 3 variants: perfect deletion plus next two
pdf("zabaleta_g1_top3.pdf",
    width = 4, height = 1.5, bg = "white")
plotAlignments(crispr.set.g1, pam.start = c(15),
               pam.end = c(20),
               target.loc = c(23),
               max.insertion.size = 100,
               top.n = 4,
               legend.cols = 4,
               plot.text.size = 1,
               add.other = FALSE,
               guide.loc = IRanges::IRanges(c(21),c(42)))
dev.off() 

# Variant type plot
byType.g1 <- crispr.set.g1$classifyVariantsByType()
unique(byType.g1)
byType.g1 <- unlist(byType.g1)
byType.g1["Vector"] <- "vector"
byType.g1["Other"] <- "other"
var_counts.g1["Vector",]

byType.g1only <- crispr.set.g1only$classifyVariantsByType()
unique(byType.g1only)
byType.g1only <- unlist(byType.g1only)
byType.g1only["Vector"] <- "vector"
byType.g1only["Other"] <- "other"
var_counts.g1only["Vector",]

var_counts_plot.g1 <- apply(var_counts.g1,2,function(x){x/sum(x)*100})
var_counts_plot.g1 <- cbind.data.frame(var_counts_plot.g1, byType.g1)
var_counts_plot.g1$byType.g1[var_counts_plot.g1$byType.g1 == "SNV"] <- "no variant"
var_counts_plot.g1 <- aggregate(.~byType.g1,data=var_counts_plot.g1,FUN=sum)
var_counts_plot.g1 <- melt(var_counts_plot.g1, id = c("byType.g1"))
var_counts_plot.g1$byType.g1 <- factor(var_counts_plot.g1$byType.g1, levels = unique(byType.g1)) 

pdf("Zabaleta_g1.types.virus.pdf",
    width = 6, height = 3, bg = "white")
ggplot(var_counts_plot.g1, aes(fill=byType.g1, x = variable, y = value)) + 
  geom_bar(position='stack', stat='identity', alpha=0.75) + coord_flip() +
  theme_linedraw() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
  scale_fill_brewer(palette="Set2") +
  labs(x = "Sample", y = "% mapping reads", fill = "Type")
dev.off()

del.g1 <- get_deletions(var_counts.g1only, byType.g1only)
ins.g1 <- get_insertions(var_counts.g1only, byType.g1only)
vec.g1 <- get_vector(md.g1only)
indels.g1 <- rbind(del.g1, ins.g1, vec.g1)


pdf("g1.deletions_samples.pdf",
    width = 5, height = 6, bg = "white")
ggplot(indels.g1, aes(color=NULL, fill=Type , x = Size)) + 
  geom_histogram(binwidth=5, alpha=0.5, position = "identity") +
  scale_fill_brewer(palette="Set1") + theme_linedraw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  labs(x = "Indel Size", y = "Indel count") +
  facet_wrap(~variable,  ncol=2, strip.position = "left") 
dev.off()

pdf("g1.deletions.pdf",
    width = 4, height = 2.5, bg = "white")
ggplot(indels.g1, aes(color=NULL, fill=Type , x = Size)) + xlim(-5,100) +
  geom_histogram(aes(y = after_stat(count / sum(count))), binwidth=5, alpha=0.7, position = "identity") +
  scale_fill_manual(values=cols) + theme_linedraw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        #legend.position = c(0.7, 0.7),
        plot.title = element_text(vjust = 0, hjust=0.5, size = 10)) +
  labs(x = "Indel Size", y = "Indel fraction (%)") +
  scale_y_continuous(labels = scales::percent) +
  ggtitle("TBG−SaCas9g1")
dev.off()


################################################################################
## Guide2 samples
crispr.set.g2 <- readsToTarget(bam.file.path.g2, gd2,
                               reference = focusing.sequence.char2,
                               target.loc = target.location2,
                               names = dataname.char.g2,
                               chimeras = treat.chimeras,
                               collapse.pairs = TRUE,
                               chimera.to.target = 200)
crispr.set.g2only <- readsToTarget(bam.file.path.g2only, gd2,
                                   reference = focusing.sequence.char2,
                                   target.loc = target.location2,
                                   names = dataname.char.g2only,
                                   chimeras = treat.chimeras,
                                   collapse.pairs = TRUE,
                                   chimera.to.target = 200)

eff.g2 <- mutationEfficiency(crispr.set.g2)
var_counts.g2 <- variantCounts(crispr.set.g2)
# Add vector integrations
vcounts.g2 <- t(md[colnames(var_counts.g2),]["Total.Sequences"])
rownames(vcounts.g2) <- c("Vector")
var_counts.g2 <- rbind(var_counts.g2, vcounts.g2)
# Save variant table
write.table(var_counts.g2, "zabaleta_g2.results.vector.csv", sep = "\t", quote = FALSE, row.names = TRUE)

eff.g2only <- mutationEfficiency(crispr.set.g2only)
var_counts.g2only <- variantCounts(crispr.set.g2only)
# Add vector integrations
vcounts.g2only <- t(md[colnames(var_counts.g2only),]["Total.Sequences"])
rownames(vcounts.g2only) <- c("Vector")
var_counts.g2only <- rbind(var_counts.g2only, vcounts.g2only)

pdf("zabaleta_g2.pdf",
    width = 15, height = 7, bg = "white")
plotVariants(crispr.set.g2, col.wdth.ratio = c(1.5,1),
             plotFreqHeatmap.args = list(top.n = 30, type="proportions"),
             plotAlignments.args = list(pam.start = c(42),
                                        pam.end = c(47),
                                        target.loc = c(38), # 86 lehen
                                        max.insertion.size = 100,
                                        top.n = 30,
                                        legend.cols = 5,
                                        plot.text.size = 2,
                                        guide.loc = IRanges::IRanges(c(20),c(41))))
dev.off()

# Plot top 3 variants: perfect deletion plus next two
pdf("zabaleta_g2_top3.pdf",
    width = 4, height = 0.9, bg = "white")
plotAlignments(crispr.set.g2,
               pam.start = c(42),
               pam.end = c(47),
               target.loc = c(38), # 86 lehen
               max.insertion.size = 100,
               top.n = 4,
               legend.cols = 4,
               plot.text.size = 1,
               add.other = FALSE,
               guide.loc = IRanges::IRanges(c(20),c(41)))
dev.off()

# Variant type plot
byType.g2 <- crispr.set.g2$classifyVariantsByType()
unique(byType.g2)
byType.g2 <- unlist(byType.g2)
byType.g2["Vector"] <- "vector"
byType.g2["Other"] <- "other"
var_counts.g2["Vector",]

byType.g2only <- crispr.set.g2only$classifyVariantsByType()
unique(byType.g2only)
byType.g2only <- unlist(byType.g2only)
byType.g2only["Vector"] <- "vector"
byType.g2only["Other"] <- "other"
var_counts.g2only["Vector",]

var_counts_plot.g2 <- apply(var_counts.g2,2,function(x){x/sum(x)*100})
var_counts_plot.g2 <- cbind.data.frame(var_counts_plot.g2, byType.g2)
var_counts_plot.g2$byType.g2[var_counts_plot.g2$byType.g2 == "SNV"] <- "no variant"
var_counts_plot.g2 <- aggregate(.~byType.g2,data=var_counts_plot.g2,FUN=sum)
var_counts_plot.g2 <- melt(var_counts_plot.g2, id = c("byType.g2"))
var_counts_plot.g2$byType.g2 <- factor(var_counts_plot.g2$byType.g2, levels = unique(byType.g2)) 

pdf("Zabaleta_g2.types.virus.pdf",
    width = 6, height = 3, bg = "white")
ggplot(var_counts_plot.g2, aes(fill=byType.g2, x = variable, y = value)) + 
  geom_bar(position='stack', stat='identity', alpha=0.75) + coord_flip() +
  theme_linedraw() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
  scale_fill_brewer(palette="Set2") +
  labs(x = "Sample", y = "% mapping reads", fill = "Type")
dev.off()

del.g2 <- get_deletions(var_counts.g2only, byType.g2only)
ins.g2 <- get_insertions(var_counts.g2only, byType.g2only)
vec.g2 <- get_vector(md.g2only)
indels.g2 <- rbind(del.g2, ins.g2, vec.g2)

pdf("g2.deletions_samples.pdf",
    width = 5, height = 6, bg = "white")
ggplot(indels.g2, aes(color=NULL, fill=Type , x = Size)) + 
  geom_histogram(binwidth=5, alpha=0.5, position = "identity") +
  scale_fill_brewer(palette="Set1") + theme_linedraw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  labs(x = "Indel Size", y = "Indel count") +
  facet_wrap(~variable,  ncol=2, strip.position = "left") 
dev.off()

pdf("g2.deletions.pdf",
    width = 4, height = 2.5, bg = "white")
ggplot(indels.g2, aes(color=NULL, fill=Type , x = Size)) + xlim(-5,100) +
  geom_histogram(aes(y = after_stat(count / sum(count))), binwidth=5, alpha=0.7, position = "identity") +
  scale_fill_manual(values=cols) + theme_linedraw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        #legend.position = c(0.7, 0.7),
        plot.title = element_text(vjust = 0, hjust=0.5, size = 10)) +
  labs(x = "Indel Size", y = "Indel fraction (%)") +
  scale_y_continuous(labels = scales::percent) +
  ggtitle("TBG−SaCas9g2")
dev.off()

################################################################################
## Join g1 and g2
var_counts_plot.g1_copy <- var_counts_plot.g1
colnames(var_counts_plot.g1_copy) <- colnames(var_counts_plot.g2)
var_counts_plot.g1g2 <- rbind(var_counts_plot.g1_copy, var_counts_plot.g2)
write.table(var_counts_plot.g1g2, "Zabaletag1g2.types.vector.csv", sep = "\t", quote = FALSE, row.names = TRUE)


pdf("Zabaletag1g2.types.vector.pdf",
    width = 6, height = 4, bg = "white")
ggplot(var_counts_plot.g1g2, aes(fill=byType.g2, x = variable, y = value)) + 
  geom_bar(position='stack', stat='identity', alpha=0.75) + coord_flip() +
  theme_linedraw() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
  scale_fill_brewer(palette="Set2") +
  labs(x = "Sample", y = "% mapping reads", fill = "Type")
dev.off()

var_means.g1g2 <-  merge(md[c("V3","V8")], var_counts_plot.g1g2, by.x = "V3", by.y = "variable")
var_means.g1g2$V3 <- NULL
var_means.g1g2$value <- as.numeric(as.character(var_means.g1g2$value))
var_means.g1g2 <- aggregate(.~V8+byType.g2,data=var_means.g1g2,FUN=mean)
write.table(var_means.g1g2, "Zabaletag1g2.means.vector.csv", sep = "\t", quote = FALSE, row.names = TRUE)

pdf("Zabaletag1g2.means.vector.pdf",
    width = 5, height = 2.5, bg = "white")
ggplot(var_means.g1g2, aes(fill=byType.g2, x = V8, y = value)) + 
  geom_bar(position='stack', stat='identity', alpha=0.75) + coord_flip() +
  theme_linedraw() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
  scale_fill_brewer(palette="Set2") +
  labs(x = "Group", y = "% mapping reads", fill = "Type")
dev.off()
