# Dominican Republic Coral Population Genomics 2bRAD
### Repository Overview
This is an accompanying repository to the paper O'Donnell et al. "Reproductive Strategy of Caribbean Corals Influences Population Structure on a Microgeographic Scale". Here, we use reduced representation sequencing (2b-RAD) to analyze the population structure of regionally abundant brooding (Agaricia agaricites, Porites astreoides) and broadcasting (Orbicella faveolata) species along the southeastern coast of the Dominican Republic. This repository contains the bioinformatic and statistical scripts, R-scripts, necessary to re-create analyses, as well as the final datasets used to generate figures. Raw sequence read files can be obtained from NCBI's SRA under BioProject PRJNA1093451   

### Files in this Repository

#### DR_GitHub_Meta.csv

Metadata for sample files contained in this repository

#### DR_processing

Bioinformatic workflow for processing 2bRAD sample reads. Custom scripts used in this analysis are listed below and in the repository. This pipeline was written for a cluster which uses a SLURM scheduler. A denovo reference genome created for Agaricia agaricites is also listed in this repository and necessary for mapping Agaricia samples. 

- Perl scripts and .fasta files
  - trim2bRAD_2barcodes_dedup_N2.pl
  - adaptors.fasta
- Shell Scripts
  - batch.sh
  - RefComp.sh
- Agaricia reference genome (created following Wang Z et al (2012) The genome of flax (Linum usitatissimum) assembled de novo from short shotgun sequence reads. Plant J 72:461–473)  
  - AagarBcgIdenovo_symABCDlongref.fasta


#### DR_ANGSD

ANGSD workflow and commands used to anaylze population structure by species. Notes within the code indicate differences in workflow by species. This script is run in tandem with R-scripts to visualize population strucuture and plot final figures. Final bam sample lists, .ibsMat files, and Fst values used for plotting final paper figures for each species are also listed.

- Accompanying R-scripts
  - Fst_95Intervals.R
  - plotQC.R
  - DR_ibs_pcoa_admix.R
- Final datasets and bam files
  - bams_Ofav_final
  - bams_Past_final
  - bams_Agar_final
  - Ofav_final_LD.ibsMat
  - Past_final.ibsMat
  - Agar_final_.ibsMat
  - Ofav_HardCall_FstFinal_LD.txt
  - Past_HardCall_FstFinal_LD.txt
  - Agar_HardCall_FstFinal.txt
  - Agar_FstGlobal.txt
  - Past_FstGlobal.txt
  - Ofav_FstGlobal.txt
  - Agar_boot_fst_results.txt
  - Past_boot_fst_results.txt
  - Ofav_boot_fst_results.txt
 
  







