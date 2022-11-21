## DEGs enrichment
### social behavior
```bash
# FB_enrichment.txt; HB_enrichment.txt; MB_enrichment.txt
# kangjingliang@kangjingliangdeMacBook-Pro 三 11 16 22:15:48 ~/Documents/2021/Cleaner_wrasse/gene_expression/enrichment
cp ~/Documents/2022/Ldim_genome_Restart/PSGs/Gene_annotation.final.txt ./
cp ~/Documents/2022/Ldim_genome_Restart/PSGs/*_GOs.txt ./
extract_gene_functions -i FB_enrichment.txt -a Gene_annotation.final.txt --gene_column 1 --func_column 3 --functions social_GOs.txt --output FB_DEGs_social_GOs
less FB_DEGs_social_GOs.txt|perl -alne 'my @a=split /\t/;print "$a[2]\t$a[-2]\t$a[-1]"'|sort -u >FB_DEGs_social_GOs_uniq.txt
extract_gene_functions -i HB_enrichment.txt -a Gene_annotation.final.txt --gene_column 1 --func_column 3 --functions social_GOs.txt --output HB_DEGs_social_GOs
less HB_DEGs_social_GOs.txt|perl -alne 'my @a=split /\t/;print "$a[2]\t$a[-2]\t$a[-1]"'|sort -u >HB_DEGs_social_GOs_uniq.txt
extract_gene_functions -i MB_enrichment.txt -a Gene_annotation.final.txt --gene_column 1 --func_column 3 --functions social_GOs.txt --output MB_DEGs_social_GOs
less MB_DEGs_social_GOs.txt|perl -alne 'my @a=split /\t/;print "$a[2]\t$a[-2]\t$a[-1]"'|sort -u >MB_DEGs_social_GOs_uniq.txt
perl Combine_DEGs.pl --FB FB_DEGs_social_GOs_uniq.txt --HB HB_DEGs_social_GOs_uniq.txt --MB MB_DEGs_social_GOs_uniq.txt > Total_DEGs_social_GOs_uniq.txt
```
### immune
```bash
# kangjingliang@kangjingliangdeMacBook-Pro 三 11 16 22:15:48 ~/Documents/2021/Cleaner_wrasse/gene_expression/enrichment
extract_gene_functions -i FB_enrichment.txt -a Gene_annotation.final.txt --gene_column 1 --func_column 3 --functions immune_GOs.txt --output FB_DEGs_immune_GOs
less FB_DEGs_immune_GOs.txt|perl -alne 'my @a=split /\t/;print "$a[2]\t$a[-2]\t$a[-1]"'|sort -u >FB_DEGs_immune_GOs_uniq.txt
extract_gene_functions -i HB_enrichment.txt -a Gene_annotation.final.txt --gene_column 1 --func_column 3 --functions immune_GOs.txt --output HB_DEGs_immune_GOs
less HB_DEGs_immune_GOs.txt|perl -alne 'my @a=split /\t/;print "$a[2]\t$a[-2]\t$a[-1]"'|sort -u >HB_DEGs_immune_GOs_uniq.txt
extract_gene_functions -i MB_enrichment.txt -a Gene_annotation.final.txt --gene_column 1 --func_column 3 --functions immune_GOs.txt --output MB_DEGs_immune_GOs
less MB_DEGs_immune_GOs.txt|perl -alne 'my @a=split /\t/;print "$a[2]\t$a[-2]\t$a[-1]"'|sort -u >MB_DEGs_immune_GOs_uniq.txt
perl Combine_DEGs.pl --FB FB_DEGs_immune_GOs_uniq.txt --HB HB_DEGs_immune_GOs_uniq.txt --MB MB_DEGs_immune_GOs_uniq.txt > Total_DEGs_immune_GOs_uniq.txt
```
### sensory
```bash
# kangjingliang@kangjingliangdeMacBook-Pro 三 11 16 22:15:48 ~/Documents/2021/Cleaner_wrasse/gene_expression/enrichment
extract_gene_functions -i FB_enrichment.txt -a Gene_annotation.final.txt --gene_column 1 --func_column 3 --functions sensory_GOs.txt --output FB_DEGs_sensory_GOs
less FB_DEGs_sensory_GOs.txt|perl -alne 'my @a=split /\t/;print "$a[2]\t$a[-2]\t$a[-1]"'|sort -u >FB_DEGs_sensory_GOs_uniq.txt
extract_gene_functions -i HB_enrichment.txt -a Gene_annotation.final.txt --gene_column 1 --func_column 3 --functions sensory_GOs.txt --output HB_DEGs_sensory_GOs
less HB_DEGs_sensory_GOs.txt|perl -alne 'my @a=split /\t/;print "$a[2]\t$a[-2]\t$a[-1]"'|sort -u >HB_DEGs_sensory_GOs_uniq.txt
extract_gene_functions -i MB_enrichment.txt -a Gene_annotation.final.txt --gene_column 1 --func_column 3 --functions sensory_GOs.txt --output MB_DEGs_sensory_GOs
less MB_DEGs_sensory_GOs.txt|perl -alne 'my @a=split /\t/;print "$a[2]\t$a[-2]\t$a[-1]"'|sort -u >MB_DEGs_sensory_GOs_uniq.txt
perl Combine_DEGs.pl --FB FB_DEGs_sensory_GOs_uniq.txt --HB HB_DEGs_sensory_GOs_uniq.txt --MB MB_DEGs_sensory_GOs_uniq.txt > Total_DEGs_sensory_GOs_uniq.txt
```
## Build the gff for the selected gene families
```bash
### ORs
# Kang@fishlab3 Fri Nov 18 16:01:48 /media/HDD/cleaner_fish/genome/OR_detection/Cleaner_wrasse/group
perl Create_ORsgtf_genewise.pl > ORs_phy.gtf
scp ORs_phy.gtf kang1234@147.8.76.177:~/genome/Gene_annotation/RNA-seq/RNA-align

### Opsins
# Kang@fishlab3 Mon Nov 21 12:08:51 /media/HDD/cleaner_fish/genome/Opsin_new/Labroides_dimidiatus/filtering
perl Create_gtf_genewise.pl >Opsins_phy.gtf

### Glutamates
# Kang@fishlab3 Mon Nov 21 13:53:41 /media/HDD/cleaner_fish/genome/Glutamates/Labroides_dimidiatus/filtering
cp /media/HDD/cleaner_fish/genome/Opsin_new/Labroides_dimidiatus/filtering/Create_gtf_genewise.pl ./
perl Create_gtf_genewise.pl >Glutamates_phy.gtf

### Dopamines
# Kang@fishlab3 Mon Nov 21 13:57:51 /media/HDD/cleaner_fish/genome/Dopamines/Labroides_dimidiatus/filtering
cp /media/HDD/cleaner_fish/genome/Opsin_new/Labroides_dimidiatus/filtering/Create_gtf_genewise.pl ./
perl Create_gtf_genewise.pl >Dopamines_phy.gtf

### Crystallins
# Kang@fishlab3 Mon Nov 21 14:08:05 /media/HDD/cleaner_fish/genome/Crystallins/Labroides_dimidiatus/filtering
cp /media/HDD/cleaner_fish/genome/Opsin_new/Labroides_dimidiatus/filtering/Create_gtf_genewise.pl ./
perl Create_gtf_genewise.pl >Crystallins_phy.gtf

### Hox
# Kang@fishlab3 Mon Nov 21 14:11:56 /media/HDD/cleaner_fish/genome/Hox/Labroides_dimidiatus/filtering
cp /media/HDD/cleaner_fish/genome/Opsin_new/Labroides_dimidiatus/filtering/Create_gtf_genewise.pl ./
perl Create_gtf_genewise.pl >Hox_phy.gtf

### BMPs
# Kang@fishlab3 Mon Nov 21 14:16:36 /media/HDD/cleaner_fish/genome/BMPs/Labroides_dimidiatus/filtering
cp /media/HDD/cleaner_fish/genome/Opsin_new/Labroides_dimidiatus/filtering/Create_gtf_genewise.pl ./
perl Create_gtf_genewise.pl >BMPs_phy.gtf

### Cytochrome P450 enzymes
# Kang@fishlab3 Mon Nov 21 14:20:51 /media/HDD/cleaner_fish/genome/Cytochrome/Labroides_dimidiatus/filtering
cp /media/HDD/cleaner_fish/genome/Opsin_new/Labroides_dimidiatus/filtering/Create_gtf_genewise.pl ./
perl Create_gtf_genewise.pl >CP450_phy.gtf

### Toll-like receptors
# Kang@fishlab3 Mon Nov 21 14:23:02 /media/HDD/cleaner_fish/genome/TLRs/Labroides_dimidiatus/filtering
cp /media/HDD/cleaner_fish/genome/Opsin_new/Labroides_dimidiatus/filtering/Create_gtf_genewise.pl ./
perl Create_gtf_genewise.pl >TLRs_phy.gtf

### RIG-like receptors
# Kang@fishlab3 Mon Nov 21 14:26:16 /media/HDD/cleaner_fish/genome/RLRs/Labroides_dimidiatus/filtering
cp /media/HDD/cleaner_fish/genome/Opsin_new/Labroides_dimidiatus/filtering/Create_gtf_genewise.pl ./
perl Create_gtf_genewise.pl >RLRs_phy.gtf

### NLRs
# Kang@fishlab3 Mon Nov 21 14:27:55 /media/HDD/cleaner_fish/genome/NLRs/Labroides_dimidiatus/filtering
cp /media/HDD/cleaner_fish/genome/Opsin_new/Labroides_dimidiatus/filtering/Create_gtf_genewise.pl ./
perl Create_gtf_genewise.pl >NLRs_phy.gtf

### Taste Receptors
# Kang@fishlab3 Mon Nov 21 14:30:23 /media/HDD/cleaner_fish/genome/TRs/Labroides_dimidiatus/filtering
cp /media/HDD/cleaner_fish/genome/Opsin_new/Labroides_dimidiatus/filtering/Create_gtf_genewise.pl ./
perl Create_gtf_genewise.pl >TRs_phy.gtf

### Immunoglobulin
# Kang@fishlab3 Mon Nov 21 14:31:35 /media/HDD/cleaner_fish/genome/ImmuHeavy/Labroides_dimidiatus/filtering
perl Create_gtf_genewise.pl >ImmuHeavy_phy.gtf
less ImmuHeavy_phy.gtf

### Immediate early response genes
# Kang@fishlab3 Mon Nov 21 14:32:58 /media/HDD/cleaner_fish/genome/IERs/Labroides_dimidiatus/filtering
cp /media/HDD/cleaner_fish/genome/Opsin_new/Labroides_dimidiatus/filtering/Create_gtf_genewise.pl ./
perl Create_gtf_genewise.pl >IERs_phy.gtf
```
## Combine gtf of all predicted genes
```bash
# Kang@fishlab3 Mon Nov 21 14:56:29 /media/HDD/cleaner_fish/genome
cat */Labroides_dimidiatus/filtering/*_phy.gtf >predicted_genes_phy.gtf
cat predicted_genes_phy.gtf /media/HDD/cleaner_fish/genome/OR_detection/Cleaner_wrasse/group/ORs_phy.gtf >predicted_genes_phy_1.gtf
mv predicted_genes_phy_1.gtf predicted_genes_phy.gtf

scp predicted_genes_phy.gtf kang1234@147.8.76.177:~/genome/Gene_annotation/RNA-seq/RNA-align/
# (base) kang1234@celia-PowerEdge-T640 Mon Nov 21 14:54:33 ~/genome/Gene_annotation/RNA-seq/RNA-align/
featureCounts -a predicted_genes_phy.gtf -o predicted_genes_phy_reads_nb.txt -T 24 LD5FB.sorted.bam LD6FB.sorted.bam LD15FB.sorted.bam LD16FB.sorted.bam LD25FB.sorted.bam LD26FB.sorted.bam LD3FB.sorted.bam LD4FB.sorted.bam LD13FB.sorted.bam LD14FB.sorted.bam LD23FB.sorted.bam LD24FB.sorted.bam LD5HB.sorted.bam LD6HB.sorted.bam LD15HB.sorted.bam LD16HB.sorted.bam LD25HB.sorted.bam LD26HB.sorted.bam LD3HB.sorted.bam LD4HB.sorted.bam LD13HB.sorted.bam LD14HB.sorted.bam LD23HB.sorted.bam LD24HB.sorted.bam LD5MB.sorted.bam LD6MB.sorted.bam LD15MB.sorted.bam LD16MB.sorted.bam LD25MB.sorted.bam LD26MB.sorted.bam LD3MB.sorted.bam LD4MB.sorted.bam LD13MB.sorted.bam LD14MB.sorted.bam LD23MB.sorted.bam LD24MB.sorted.bam
cp predicted_genes_phy_reads_nb.txt read_matrix/

# (base) kang1234@celia-PowerEdge-T640 Mon Nov 21 15:56:11 ~/genome/Gene_annotation/RNA-seq/RNA-align/read_matrix
perl Creat_matrix_for_coefficient.pl >all_inds_read_nb_for_coefficient.txt
```
**TPM-transform**   
```R
# we just used subset genes, thus we used the colSums(rpk) from the data of all individuals
# (base) kang1234@celia-PowerEdge-T640 Fri Nov 18 16:35:51 ~/genome/Gene_annotation/RNA-seq/RNA-align/read_matrix
countdata<-read.table("all_inds_read_nb_for_coefficient.txt",skip = 1,sep="\t",header = T,row.names = 1)
metadata <- countdata[,1:5]
countdata <- countdata[,6:ncol(countdata)]
kb <- metadata$Length / 1000
rpk1 <- countdata / kb

countdata<-read.table("predicted_genes_phy_reads_nb.txt",skip = 1,sep="\t",header = T,row.names = 1)
metadata <- countdata[,1:5]
countdata <- countdata[,6:ncol(countdata)]
prefix<-"predicted_genes_phy_reads_nb"
kb <- metadata$Length / 1000
rpk <- countdata / kb
tpm <- t(t(rpk)/colSums(rpk1) * 1000000)
write.csv(tpm,paste0(prefix,"_tpm.csv"))
```
**TPM value: round number to 2 digits after decimal point**   
```bash
# (base) kang1234@celia-PowerEdge-T640 Mon Nov 21 17:04:00 ~/genome/Gene_annotation/RNA-seq/RNA-align/read_matrix
less predicted_genes_phy_reads_nb_tpm.csv|perl -alne 's/\.sorted\.bam//g;$i++;@F=split /\,/;print if $i==1;if ($i>1){my $info=$F[0].",";for(my $i = 1; $i < @F; $i++){my $a=$F[$i]; $a=sprintf("%.2f",$a);$info.=$a.","};$info=~s/\,$//;print "$info"}' >predicted_genes_phy_reads_nb_tpm_1.csv
mv predicted_genes_phy_reads_nb_tpm_1.csv predicted_genes_phy_reads_nb_tpm.csv
```
