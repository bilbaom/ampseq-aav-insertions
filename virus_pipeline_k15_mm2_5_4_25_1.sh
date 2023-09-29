#!/bin/bash

# Preprocessing pipeline that is used before indel characterisation.
# 1. It first keeps the reads that map to the amplified region
# 2. Quantify the reads that harbour AAV vector insertions and remove them
# 3. Vector-free reads are mapped with Minimap2 with alignment parameters -A5 -B4 -O25 -E1
# Necesary packages in conda environment: BWA, samtools, flash, minimap2, bbtools.

conda activate alignments
# Index genome for alignment (once)
bwa index -p /home/bilbaom/ph1/ref/GCF_000001635.26_GRCm38.p6_genomic  /home/bilbaom/ph1/ref/GCF_000001635.26_GRCm38.p6_genomic.fna
minimap2 -d /home/bilbaom/ph1/ref/GCF_000001635.26_GRCm38.p6_genomic.mmi /home/bilbaom/ph1/ref/GCF_000001635.26_GRCm38.p6_genomic.fna -k21 -w11

# index virus genome for alignment (once)
bwa index -p /home/bilbaom/ph1/ref/aat-d10ag1_N /home/bilbaom/ph1/ref/aat-d10ag1_N.fasta
bwa index -p /home/bilbaom/ph1/ref/aat-d10ag1g2_N /home/bilbaom/ph1/ref/aat-d10ag1g2_N.fasta
bwa index -p /home/bilbaom/ph1/ref/tbg-d10ag1_N /home/bilbaom/ph1/ref/tbg-d10ag1_N.fasta
bwa index -p /home/bilbaom/ph1/ref/tbg-sacas9g1_N /home/bilbaom/ph1/ref/tbg-sacas9g1_N.fasta

mkdir -p bam
mkdir -p bam_fractions
mkdir -p seq_fractions
mkdir -p bamv
mkdir -p bamv_final
mkdir -p bam_final
mkdir -p readlists
mkdir -p flash


function analizefiles {
	echo "Starting sample ${sample}"

	# Map to genome (to ensure that both reads come from targeted amplicon)
	bwa mem -t 4 /home/bilbaom/ph1/ref/GCF_000001635.26_GRCm38.p6_genomic \
		${fastq_r1} \
		${fastq_r2} | samtools sort -o bam/${sample}.fastq.bam

	# Get mapped to amplicon (R1 mapped and R2 mapped) 
	samtools view -b -h -F 12 bam/${sample}.fastq.bam | samtools sort -o bam_fractions/${sample}_map_map.bam -
	
	# Filter by location (ON-TARGET: NC_000068.7:134548000-134548500)
	samtools index bam_fractions/${sample}_map_map.bam
	samtools view -b -h bam_fractions/${sample}_map_map.bam NC_000068.7:134548000-134548500 > bam_fractions/${sample}_map_map.ont.bam
	
	# Sort by read name for fastq
	samtools sort -n bam_fractions/${sample}_map_map.ont.bam -o bam_fractions/${sample}_map_map.ont.sort.bam
	
	# Get mapped FASTQ
	samtools fastq -@ 4 \
		-1 seq_fractions/${sample}_mapped_1.fastq \
		-2 seq_fractions/${sample}_mapped_2.fastq \
		-0 /dev/null -s /dev/null -n bam_fractions/${sample}_map_map.ont.sort.bam

	gzip -f seq_fractions/${sample}_mapped*
	
	# Merge paired end reads
	flash --min-overlap 10 --max-overlap 250 -t 4 -o ${sample} -d flash/ \
		seq_fractions/${sample}_mapped_1.fastq.gz \
		seq_fractions/${sample}_mapped_2.fastq.gz > flash/${sample}.flash.log 2>&1
	
	# Map MERGED to genome (for indel detection with all reads)
	minimap2 -A5 -B4 -O25 -E1 \
		-a /home/bilbaom/ph1/ref/GCF_000001635.26_GRCm38.p6_genomic.mmi \
		flash/${sample}.extendedFrags.fastq | samtools sort -o bam_final/${sample}_mapped.bam

	# Map  MERGED to virus genome
	bwa mem -t 4 -k 15 /home/bilbaom/ph1/ref/${virus_seq} \
		flash/${sample}.extendedFrags.fastq | samtools sort -o bamv/${sample}_virus.bam
	
	# Get MERGED mapped to virus
	samtools view -h -F 4 bamv/${sample}_virus.bam > bamv/${sample}_virus_mapped.bam
	samtools sort -n bamv/${sample}_virus_mapped.bam -o bamv_final/${sample}_virus_mapped.sort.bam
	
	# Map UNMERGED to virus genome
	bwa mem -t 4 -k 15 /home/bilbaom/ph1/ref/${virus_seq} \
		flash/${sample}.notCombined_1.fastq \
		flash/${sample}.notCombined_2.fastq | samtools sort -o bamv/${sample}_virus_notCombined.bam
		
	# Get UNMERGED mapped to virus (R1 or R2)
	samtools view -b -F 4 -f 8 bamv/${sample}_virus_notCombined.bam > tmps1.bam  # first mate mapped
	samtools view -b -F 8 -f 4 bamv/${sample}_virus_notCombined.bam > tmps2.bam  # second mate mapped
	samtools view -b -F 12 bamv/${sample}_virus_notCombined.bam > tmps3.bam  # both mates mapped
	samtools merge -u - tmps[123].bam | samtools sort -n -o bamv_final/${sample}_virus_notCombined_mapped.sort.bam
	
	# Get virus mapped FASTQ
	samtools fastq -@ 4 \
		bamv_final/${sample}_virus_mapped.sort.bam > seq_fractions/${sample}_virus_mapped.fastq
	
	gzip -f seq_fractions/${sample}_virus*
	
	# Map merged reads with virus to amplicon
	minimap2 -A5 -B4 -O25 -E1 \
		-a /home/bilbaom/ph1/ref/GCF_000001635.26_GRCm38.p6_genomic.mmi \
		seq_fractions/${sample}_virus_mapped.fastq.gz | samtools sort -o bam_final/${sample}_v.fastq.bam

	# Filter amplicon mapped to get virus-free BAM (for indel detection without vector reads)
	# Get read names of reads that map to the vector
	samtools view bam_final/${sample}_v.fastq.bam | cut -f1 | sort | uniq > readlists/${sample}_vir_reads.txt
	samtools view bamv_final/${sample}_virus_notCombined_mapped.sort.bam | cut -f1 | sort | uniq > readlists/${sample}_vir_reads_notCombined.txt
	# Filter out those reads from all the mapped reads
	filterbyname.sh \
		in=flash/${sample}.extendedFrags.fastq \
		out=seq_fractions/${sample}_map_map.vfree.fastq.gz \
		names=readlists/${sample}_vir_reads.txt
	
	minimap2 -A5 -B4 -O25 -E1 \
		-a /home/bilbaom/ph1/ref/GCF_000001635.26_GRCm38.p6_genomic.mmi \
		seq_fractions/${sample}_map_map.vfree.fastq.gz | samtools sort -o bam_final/${sample}_map_map.vfree.sort.bam
}


# Nickase and Sacas9 samples
grep "home" /home/bilbaom/ph1/all_samples.csv | while read sample fastq_r1 fastq_r2 virus_seq guide_seq guide_name; do
  analizefiles
done

# Zabaleta2018 samples
grep "home" /home/bilbaom/ph1/zabaleta_samples_k11.csv | while read sample fastq_r1 fastq_r2 virus_seq guide_seq guide_name; do
  analizefiles
done

conda deactivate