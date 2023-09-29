# ampseq-aav-insertions

This repository contains the scripts used to analyse gene editing from amplicon sequencing data in the manuscript: [under revision]

## Scripts
* [virus_pipeline_k15_mm2_5_4_25_1.sh](/virus_pipeline_k15_mm2_5_4_25_1.sh): **Preprocessing pipeline that is used before indel characterisation.** It first keeps the reads that map to the amplified region, and then quantifies the reads that harbour AAV vector insertions and removes them. Finally vector-free reads are mapped with Minimap2 with alignment parameters -A5 -B4 -O25 -E1 (based on parameters used in Amplican and CRISPECTOR tools).
  * Necesary packages in conda environment: BWA, samtools, flash, minimap2, bbtools.

* [run_crisprvariants.R](/run_crisprvariants.R): **Indel characterisation.** Analyze vector-free BAM files with crisprvariants R package and add previously quantified vector insertions.
  
* [run_crisprvariants_zabaleta2018.R](/run_crisprvariants_zabaleta2018.R): **Indel characterisation of Zabaleta et al. (2018) data.** Analyze vector-free BAM files with crisprvariants R package and add previously quantified vector insertions.
