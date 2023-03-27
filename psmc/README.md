## PSMC
```bash
# hisat2: Genome fasta file index
# (base) kang1234@celia-PowerEdge-T640 Sat Mar 11 23:36:11 ~/genome/Gene_annotation
hisat2-build -p 22 -f Cleaner_wrasse_softmasked_ChaHeader_final.fasta Ldim
```
### Pre PSMC
```bash
# Pre_PSMC.sh
hisat2 -p 22 --dta -q -x Ldim -1 ../DTG-OmniC-105_R1_001.fastq.gz -2 ../DTG-OmniC-105_R2_001.fastq.gz -S Ldim_genome.sam --summary Ldim_genome.txt
samtools view -bS Ldim_genome.sam >Ldim_genome.bam
samtools sort Ldim_genome.bam -o Ldim_genome_sorted.bam
samtools index Ldim_genome_sorted.bam
samtools mpileup -Q 30 -q 30 -uf Cleaner_wrasse_softmasked_ChaHeader_final.fasta Ldim_genome_sorted.bam|bcftools call -c - |vcfutils.pl vcf2fq -d 5 -D 100 -Q 30 > Ldim_genome.fq
```
```bash
# (base) kang1234@celia-PowerEdge-T640 Sat Mar 11 23:59:27 ~/genome/Gene_annotation
nohup bash -x Pre_PSMC.sh > Pre_PSMC.process 2>&1 &
# [1] 31138
```
