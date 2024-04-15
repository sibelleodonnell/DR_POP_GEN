#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4GB
#SBATCH --time=00:45:00
#SBATCH --error=%x_%j.err
#SBATCH --output=%x_%j.out

export SCRIPTS='/project/ckenkel_26/scripts'
export REFS='/project/ckenkel_26/RefSeqs'


module load gcc/8.3.0
module load gcc/9.2.0
module load intel/18.0.4
module load intel/19.0.4
module load gcc/11.3.0
module load openblas/0.3.20
module load zlib
module load bzip2
module load htslib
module load jdk
module load perl
module load python
module load py-cython
module load bowtie2
module load curl
module load fastqc
module load git
module load gsl/2.5
module load ncbi-rmblastn/2.11.0
module load parallel
module load picard
module load plink
module load r/4.2.3
module load samtools
module load screen
module load slurm/20.02.4
module load sra-toolkit/2.9.6
module load structure
module load tar
module load vcftools
module load bcftools
module load vim

export PATH=/project/ckenkel_26/scripts:$PATH
export PATH=/project/ckenkel_26/software:$PATH
export PATH=/project/ckenkel_26/software/binAIMS:$PATH
export PATH=/project/ckenkel_26/software/bin:$PATH
export PATH=/project/ckenkel_26/software/bbmap:$PATH
export PATH=/project/ckenkel_26/software/SHRiMP_2_2_3:$PATH
export PATH=/project/ckenkel_26/software/SHRiMP_2_2_3/bin:$PATH
export PATH=/project/ckenkel_26/software/cdhit:$PATH
export PATH=/project/ckenkel_26/software/ngsPopGen:$PATH
export PATH=/project/ckenkel_26/software/dist/admixture_linux-1.3.0:$PATH
export PATH=/project/ckenkel_26/software/PopGenTools:$PATH
export PATH=/project/ckenkel_26/software/PGDSpider_2.0.7.1:$PATH
export PATH=/project/ckenkel_26/software/standard-RAxML-master:$PATH
export PATH=/project/ckenkel_26/software/ncbi-blast-2.10.1+/bin:$PATH
export PATH=/project/ckenkel_26/software/angsd/bin:$PATH

export PYTHONPATH=$PYTHONPATH:/project/ckenkel_26/software/moments

#job parameters
cat $NAME".tr0" | fastq_quality_filter -q 20 -p 90 > $NAME".trim2" #relaxed percent of reads mapping to increase retention of reads
bbduk.sh in=$NAME".trim2" ref=/project/ckenkel_26/RefSeqs/adaptors.fasta k=12 stats=stats2.txt out=$NAME".clean2"
GENOME_FASTA=/project/ckenkel_26/RefSeqs/PastGenome/V2lessRepeats/Past_symABCD_genomev2.fasta
bowtie2 --threads 16 --no-unal --score-min L,16,1 --local -L 16 -x $GENOME_FASTA -U $NAME".clean2" -S $NAME"bt2.sam2"


