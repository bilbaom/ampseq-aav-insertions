# ampseq-aav-insertions

This repository contains the scripts used to analyse gene editing and AAV vector insertion from amplicon sequencing data in the manuscript:
"Efficient and safe therapeutic use of paired Cas9-nickases for primary hyperoxaluria type 1" [https://doi.org/10.1038/s44321-023-00008-8](https://doi.org/10.1038/s44321-023-00008-8)

## Scripts
* [virus_pipeline_k15_mm2_5_4_25_1.sh](/virus_pipeline_k15_mm2_5_4_25_1.sh): **Preprocessing pipeline that is used before indel characterisation.** It first keeps the reads that map to the amplified region, and then quantifies the reads that harbour AAV vector insertions and removes them. Finally vector-free reads are mapped with Minimap2 with alignment parameters -A5 -B4 -O25 -E1 (based on parameters used in Amplican and CRISPECTOR tools).
  * Necesary packages in conda environment: BWA, samtools, flash, minimap2, bbtools.

* [run_crisprvariants.R](/run_crisprvariants.R): **Indel characterisation.** Analyze vector-free BAM files with crisprvariants R package and add previously quantified vector insertions. Reproduce manuscript plots.
  
* [run_crisprvariants_zabaleta2018.R](/run_crisprvariants_zabaleta2018.R): **Indel characterisation of Zabaleta et al. (2018) data.** Analyze vector-free BAM files with crisprvariants R package and add previously quantified vector insertions. Reproduce manuscript plots.
