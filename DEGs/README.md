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

# (base) kang1234@celia-PowerEdge-T640 Wed Nov 23 10:16:10 ~/genome/Gene_annotation/RNA-seq/RNA-align/read_matrix
less all_inds_read_nb_for_coefficient_tpm.txt|perl -alne 's/\.sorted\.bam//g;$i++;print if $i==1;if ($i>1){my $info=$F[0]."\t";for(my $i = 1; $i < @F; $i++){my $a=$F[$i]; $a=sprintf("%.2f",$a);$info.=$a."\t"};$info=~s/\s+$//;print "$info"}' >all_inds_read_nb_for_coefficient_tpm_1.txt
mv all_inds_read_nb_for_coefficient_tpm_1.txt all_inds_read_nb_for_coefficient_tpm.txt

# social
# kangjingliang@kangjingliangdeMacBook-Pro 三 11 23 10:08:43 ~/Documents/2021/Cleaner_wrasse/gene_expression/enrichment
less Total_DEGs_social_GOs_uniq.txt|cut -f 1 >DEGs_social_GOs.txt
scp kang1234@147.8.76.177:~/genome/Gene_annotation/RNA-seq/RNA-align/read_matrix/all_inds_read_nb_for_coefficient_tpm.txt ./
extract_reads_nb --matrix all_inds_read_nb_for_coefficient_tpm.txt --genes DEGs_social_GOs.txt --samples coldata.txt >DEGs_social_reads_nb_tpm.txt

# sensory
less Total_DEGs_sensory_GOs_uniq.txt|cut -f 1 >DEGs_sensory_GOs.txt
extract_reads_nb --matrix all_inds_read_nb_for_coefficient_tpm.txt --genes DEGs_sensory_GOs.txt --samples coldata.txt >DEGs_sensory_reads_nb_tpm.txt

# All DEGs
# kangjingliang@kangjingliangdeMacBook-Pro 四 11 24 10:58:55 ~/Documents/2021/Cleaner_wrasse/gene_expression
perl temp1.pl gtf_FB_Interaction_Solo.DEGs.txt >gtf_FB_Interaction_Solo.DEGs.ano.txt
perl temp1.pl gtf_HB_Interaction_Solo.DEGs.txt >gtf_HB_Interaction_Solo.DEGs.ano.txt
perl temp1.pl gtf_MB_Interaction_Solo.DEGs.txt >gtf_MB_Interaction_Solo.DEGs.ano.txt
cp gtf_*_Interaction_Solo.DEGs.ano.txt enrichment/
# kangjingliang@kangjingliangdeMacBook-Pro 四 11 24 11:02:33 ~/Documents/2021/Cleaner_wrasse/gene_expression/enrichment
perl Combine_DEGs.pl --FB gtf_FB_Interaction_Solo.DEGs.ano.txt --HB gtf_HB_Interaction_Solo.DEGs.ano.txt --MB gtf_MB_Interaction_Solo.DEGs.ano.txt > Total_DEGs_uniq.txt
# Check the regulation tendency of all DEGs across FB, HB, MB
# kangjingliang@kangjingliangdeMacBook-Pro 五 11 25 09:35:02 ~/Documents/2021/Cleaner_wrasse/gene_expression/enrichment
less Total_DEGs_uniq.txt|perl -alne 'my @F=split /\t/;my $a;for (my $i = 3; $i < @F; $i+=2){my @k=split /\;/, $F[$i];$a.=$k[-1]."\t"};$a=~s/\s+$//;$hash{$a}++;END{foreach my $key (sort keys %hash){print "$key\t$hash{$key}"}}'
# FB-DEG	HB-DEG	MB-DEG	65
# FB-DEG	HB-DEG	MB-NONDEG	407
# FB-DEG	HB-NONDEG	MB-DEG	253
# FB-DEG	HB-NONDEG	MB-NONDEG	2010
# FB-NONDEG	HB-DEG	MB-DEG	35
# FB-NONDEG	HB-DEG	MB-NONDEG	1075
# FB-NONDEG	HB-NONDEG	MB-DEG	159

# kangjingliang@kangjingliangdeMacBook-Pro 五 11 25 09:48:49 ~/Documents/2021/Cleaner_wrasse/gene_expression/enrichment
less Total_DEGs_uniq.txt|perl -alne 'my @F=split /\t/;my $a;for (my $i = 4; $i < @F; $i+=2){my @k=split /\;/, $F[$i];$a.=$k[-1]."\t"};$a=~s/\s+$//;$hash{$a}++;END{foreach my $key (sort keys %hash){print "$key\t$hash{$key}"}}'
# DOWN	DOWN	DOWN	1008
# DOWN	DOWN	UP	64
# DOWN	UP	DOWN	256
# DOWN	UP	UP	138
# UP	DOWN	DOWN	51
# UP	DOWN	UP	169
# UP	UP	DOWN	40
# UP	UP	UP	2278

# Functional enrichment for DEGs in at least one tissue
# kangjingliang@kangjingliangdeMacBook-Pro 五 11 25 09:57:45 ~/Documents/2021/Cleaner_wrasse/gene_expression
cat gtf_FB_Interaction_Solo.DEGs.txt gtf_HB_Interaction_Solo.DEGs.txt gtf_MB_Interaction_Solo.DEGs.txt|sort -u >Total_DEGs.txt
# Functional enrichenment result: Total_DEGs_enrichment.txt
# kangjingliang@kangjingliangdeMacBook-Pro 五 11 25 15:10:24 ~/Documents/2021/Cleaner_wrasse/gene_expression/enrichment
perl temp3.pl > Total_DEGs_enrichment_regulation_ratio.txt
# Plot: Total_DEGs_enrichment_regulation_ratio_plot.csv (GO_Category="BIOLOGICAL_PROCESS"; Total_Nb >= 20)

# Functional enrichment reduced
perl temp4.pl > Total_DEGs_enrichment_reduced_regulation_ratio.txt
# Plot: Func_regulation.R
# Extract the genes of GOs in the plot: Plot_GOs.txt
# kangjingliang@kangjingliangdeMacBook-Pro 五 11 25 17:08:44 ~/Documents/2021/Cleaner_wrasse/gene_expression/enrichment
mv Total_DEGs_enrichment_reduced.txt Total_DEGs_reduced_enrichment.txt
extract_gene_functions -i Total_DEGs_reduced_enrichment.txt -a Gene_annotation.final.txt --gene_column 1 --func_column 3 --functions Plot_GOs.txt --output Total_DEGs_reduced_enrichment_GOs
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

### Protocadherins
# Kang@fishlab3 Tue Nov 22 23:06:18 /media/HDD/cleaner_fish/genome/Protocadherin_new/Labroides_dimidiatus/filtering
cp /media/HDD/cleaner_fish/genome/Opsin_new/Labroides_dimidiatus/filtering/Create_gtf_genewise.pl ./
perl Create_gtf_genewise.pl >Protocadherins_phy.gtf
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

# Kang@fishlab3 Tue Nov 22 23:09:36 /media/HDD/cleaner_fish/genome/Protocadherin_new/Labroides_dimidiatus/filtering
scp Protocadherins_phy.gtf kang1234@147.8.76.177:~/genome/Gene_annotation/RNA-seq/RNA-align/
featureCounts -a Protocadherins_phy.gtf -o predicted_protocadherins_phy_reads_nb.txt -T 24 LD5FB.sorted.bam LD6FB.sorted.bam LD15FB.sorted.bam LD16FB.sorted.bam LD25FB.sorted.bam LD26FB.sorted.bam LD3FB.sorted.bam LD4FB.sorted.bam LD13FB.sorted.bam LD14FB.sorted.bam LD23FB.sorted.bam LD24FB.sorted.bam LD5HB.sorted.bam LD6HB.sorted.bam LD15HB.sorted.bam LD16HB.sorted.bam LD25HB.sorted.bam LD26HB.sorted.bam LD3HB.sorted.bam LD4HB.sorted.bam LD13HB.sorted.bam LD14HB.sorted.bam LD23HB.sorted.bam LD24HB.sorted.bam LD5MB.sorted.bam LD6MB.sorted.bam LD15MB.sorted.bam LD16MB.sorted.bam LD25MB.sorted.bam LD26MB.sorted.bam LD3MB.sorted.bam LD4MB.sorted.bam LD13MB.sorted.bam LD14MB.sorted.bam LD23MB.sorted.bam LD24MB.sorted.bam
mv predicted_protocadherins_phy_reads_nb.txt read_matrix/
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
write.table(tpm,paste0(prefix,"_tpm.txt"),quote = F, sep = "\t", col.names = NA, row.names = TRUE)
```
**TPM value: round number to 2 digits after decimal point**   
```bash
# (base) kang1234@celia-PowerEdge-T640 Mon Nov 21 17:04:00 ~/genome/Gene_annotation/RNA-seq/RNA-align/read_matrix
less predicted_genes_phy_reads_nb_tpm.txt|perl -alne 's/\.sorted\.bam//g;$i++;print if $i==1;if ($i>1){my $info=$F[0]."\t";for(my $i = 1; $i < @F; $i++){my $a=$F[$i]; $a=sprintf("%.2f",$a);$info.=$a."\t"};$info=~s/\s+$//;print "$info"}' >predicted_genes_phy_reads_nb_tpm_1.txt
mv predicted_genes_phy_reads_nb_tpm_1.txt predicted_genes_phy_reads_nb_tpm.txt
```
## Plot the heatmap of gene family expression one by one
**extract the matrix**   
### ORs
```bash
# kangjingliang@kangjingliangdeMacBook-Pro 二 11 22 11:17:04 ~/Documents/2021/Cleaner_wrasse/gene_expression/Gene_families
cat ../coldata_*.txt|perl -alne '$i++;s/\r//g;if (/^\s+/ && $i==1){print}elsif(/^LD/ && $i>1){print}' >coldata.txt
less predicted_genes_phy_reads_nb_tpm.txt|perl -alne 'print $F[0] if /L\.dimidiatus/'|sort >ORs.txt
extract_reads_nb --matrix predicted_genes_phy_reads_nb_tpm.txt --genes ORs.txt --samples coldata.txt >ORs_reads_nb_tpm.txt

# Count the number of genes with "expression == 0" and the maximun expression per gene
# kangjingliang@kangjingliangdeMacBook-Pro 二 11 29 10:30:39 ~/Documents/2021/Cleaner_wrasse/gene_expression/Gene_families
less ORs_reads_nb_tpm.txt|perl -alne '$k++;print "Gene\tNbeq0\tNbeq0_ratio\tMaximun_Nb\tMean_Nb" if $k==1;next if /^\s+/;my ($j, $max, $mean);for (my $i = 1; $i < @F; $i++){$j++ if $F[$i]==0;$mean+=$F[$i];($F[$i]<$max)?($max=$max):($max=$F[$i])};($j>0)?($j=$j):($j=0);my $rat=$j/36;$rat=sprintf("%.3f",$rat);$mean=$mean/36;$mean=sprintf("%.3f",$mean);print "$F[0]\t$j\t$rat\t$max\t$mean"' >ORs_Expression_info.txt

# Plot the expression of ORs if its expression equals to 0 in less than 80% individuals (< 28 individuals: 36 individuals in total)
cp ~/Documents/2021/Cleaner_wrasse/gene_expression/enrichment/Create_reads_nb_gene.pl ./
# kangjingliang@kangjingliangdeMacBook-Pro 二 11 29 11:18:51 ~/Documents/2021/Cleaner_wrasse/gene_expression/Gene_families
vi ORs_expression_less_80_0.txt
perl Create_reads_nb_gene.pl --matr ORs_reads_nb_tpm.txt --gene L.dimidiatus.Kappa-1,L.dimidiatus.Delta-7,L.dimidiatus.Delta-13,L.dimidiatus.Delta-9,L.dimidiatus.Delta-8,L.dimidiatus.Delta-10,L.dimidiatus.Delta-6,L.dimidiatus.Delta-32,L.dimidiatus.Delta-11,L.dimidiatus.Delta-14,L.dimidiatus.Beta-1,L.dimidiatus.Delta-19,L.dimidiatus.Delta-29,L.dimidiatus.Delta-15,L.dimidiatus.Eta-1,L.dimidiatus.Delta-31,L.dimidiatus.Delta-30,L.dimidiatus.Delta-12 > ORs_expression_less_80_0_tpm_plot.txt
```
**Plot**   
```R
setwd("~/Documents/2021/Cleaner_wrasse/gene_expression/Gene_families")
library(pheatmap)
data<-read.table("ORs_reads_nb_tpm.txt", header = TRUE, row.names = 1)
coldata <- read.table(file="coldata.txt") 
df<-coldata
ann_colors = list(
  Group = c(Solo = "#8491B4", Interaction = "#3C5488")
)
pheatmap(data, cluster_rows=FALSE, cluster_cols=FALSE, 
         show_rownames=TRUE, show_colnames = FALSE, 
         gaps_col = c(12, 24), gaps_row = c(3, 37, 40, 42, 43),
         color = colorRampPalette(c("#fcf9cb", "#ea8544", "#6a1928"))(50),
         annotation_col=df, annotation_colors = ann_colors,border_color=NA,
         fontfamily= "Times New Roman",fontsize = 8)
```
### Opsins
```bash
# Kang@fishlab3 Tue Nov 22 14:02:59 /media/HDD/cleaner_fish/genome/Opsin_new
less Predict_Opsins.phy|perl -alne 'next if /^\d+/;if ($F[0]=~/(.*?)_(.*?)_.*/){$hash{$2}->{$1}++};END{foreach my $ge (sort keys %hash){my $info; foreach $sp (qw(Ldim Tbif Smel Tads Spul Lber Ncel Cund)){my $nb; $hash{$ge}->{$sp}?($nb=$hash{$ge}->{$sp}):($nb=0);$info.=$nb.":";};$info=~s/\:$//;print "$ge\t[$info]"}}' >Opsins_nb.txt
cp /media/HDD/cleaner_fish/genome/BMPs/temp2.pl ./
perl temp2.pl Opsins_nb.txt >Opsins_nb_2.txt
# kangjingliang@kangjingliangdeMacBook-Pro 二 11 22 14:26:21 ~/Documents/2021/Cleaner_wrasse/gene_expression/Gene_families
extract_reads_nb --matrix predicted_genes_phy_reads_nb_tpm.txt --genes Opsins.txt --samples coldata.txt >Opsins_reads_nb_tpm.txt

# Barplot
# kangjingliang@kangjingliangdeMacBook-Pro 二 11 29 12:01:50 ~/Documents/2021/Cleaner_wrasse/gene_expression/Gene_families
perl Create_reads_nb_gene.pl --matr Opsins_reads_nb_tpm.txt --gene OPSR_1C,OPSV_1C,OPSB_1T,OPSD_1C,OPSD_2C,OPSG_3T,OPSG_2C,OPSG4_1C > Opsins_reads_nb_tpm_plot.txt
```

### Glutamates
```bash
# Kang@fishlab3 Tue Nov 22 15:49:12 /media/HDD/cleaner_fish/genome/Glutamates
less Predict_Glutamates.phy|perl -alne 'next if /^\d+/;if ($F[0]=~/(.*?)_(.*?)_.*/){$hash{$2}->{$1}++};END{foreach my $ge (sort keys %hash){my $info; foreach $sp (qw(Ldim Tbif Smel Tads Spul Lber Ncel Cund)){my $nb; $hash{$ge}->{$sp}?($nb=$hash{$ge}->{$sp}):($nb=0);$info.=$nb.":";};$info=~s/\:$//;print "$ge\t[$info]"}}' >Glutamates_nb.txt
cp /media/HDD/cleaner_fish/genome/BMPs/temp2.pl ./
perl temp2.pl Glutamates_nb.txt >Glutamates_nb_2.txt
# kangjingliang@kangjingliangdeMacBook-Pro 二 11 22 14:54:36 ~/Documents/2021/Cleaner_wrasse/gene_expression/Gene_families
extract_reads_nb --matrix predicted_genes_phy_reads_nb_tpm.txt --genes Glutamates.txt --samples coldata.txt >Glutamates_reads_nb_tpm.txt
```

### Dopamine receptors
```bash
# kangjingliang@kangjingliangdeMacBook-Pro 二 11 22 15:51:24 ~/Documents/2021/Cleaner_wrasse/gene_expression/Gene_families
less predicted_genes_phy_reads_nb_tpm.txt|grep -i 'DRD'|cut -f 1|sort -u >Dopamines.txt
extract_reads_nb --matrix predicted_genes_phy_reads_nb_tpm.txt --genes Dopamines.txt --samples coldata.txt >Dopamines_reads_nb_tpm.txt
```

### Toll-like receptors; RIG-like receptors; NLRs
```bash
# Kang@fishlab3 Tue Nov 22 17:53:57 /media/HDD/cleaner_fish/genome/TLRs
less Predict_TLRs.phy|perl -alne 'next if /^\d+/;if ($F[0]=~/(.*?)_(.*?)_.*/){$hash{$2}->{$1}++};END{foreach my $ge (sort keys %hash){my $info; foreach $sp (qw(Ldim Tbif Smel Tads Spul Lber Ncel Cund)){my $nb; $hash{$ge}->{$sp}?($nb=$hash{$ge}->{$sp}):($nb=0);$info.=$nb.":";};$info=~s/\:$//;print "$ge\t[$info]"}}' >TLRs_nb.txt
cp /media/HDD/cleaner_fish/genome/BMPs/temp2.pl ./
perl temp2.pl TLRs_nb.txt >TLRs_nb_2.txt

# Kang@fishlab3 Tue Nov 22 18:00:34 /media/HDD/cleaner_fish/genome/RLRs
less Predict_RLRs.phy|perl -alne 'next if /^\d+/;if ($F[0]=~/(.*?)_(.*?)_.*/){$hash{$2}->{$1}++};END{foreach my $ge (sort keys %hash){my $info; foreach $sp (qw(Ldim Tbif Smel Tads Spul Lber Ncel Cund)){my $nb; $hash{$ge}->{$sp}?($nb=$hash{$ge}->{$sp}):($nb=0);$info.=$nb.":";};$info=~s/\:$//;print "$ge\t[$info]"}}' >RLRs_nb.txt
cp /media/HDD/cleaner_fish/genome/BMPs/temp2.pl ./
perl temp2.pl RLRs_nb.txt >RLRs_nb_2.txt
# kangjingliang@kangjingliangdeMacBook-Pro 二 11 22 18:15:30 ~/Documents/2021/Cleaner_wrasse/gene_expression/Gene_families
perl temp1.pl|sort -u
extract_reads_nb --matrix predicted_genes_phy_reads_nb_tpm.txt --genes Immune.txt --samples coldata.txt >Immune_reads_nb_tpm.txt

less Immune_reads_nb_tpm.txt|perl -alne '$k++;print "Gene\tNbeq1\tNbeq1_ratio\tMaximun_Nb\tMean_Nb" if $k==1;next if /^\s+/;my ($j, $max, $mean);for (my $i = 1; $i < @F; $i++){$j++ if $F[$i]<1;$mean+=$F[$i];($F[$i]<$max)?($max=$max):($max=$F[$i])};($j>0)?($j=$j):($j=0);my $rat=$j/36;$rat=sprintf("%.3f",$rat);$mean=$mean/36;$mean=sprintf("%.3f",$mean);print "$F[0]\t$j\t$rat\t$max\t$mean"' > Immune_Expression_info.txt
vi Immune_expression_less_80_1.txt
# Expression TPM >=1 in 80% indviduals
perl Create_reads_nb_gene.pl --matr Immune_reads_nb_tpm.txt --gene DHX58_1C,IFIH1_1C,NAL12_4C,NAL12_5C,NLRX1_1C,NOD2_1C,TLR13_3C,TLR1_1C,TLR21_1C,TLR3_1C > Immune_reads_nb_tpm_plot.txt
```

### Protocadherins: alpha; gamma
```bash
# kangjingliang@kangjingliangdeMacBook-Pro 二 11 22 23:31:16 ~/Documents/2021/Cleaner_wrasse/gene_expression/Gene_families
perl temp2.pl|sort -u >protocadherins.txt
extract_reads_nb --matrix predicted_protocadherins_phy_reads_nb_tpm.txt --genes protocadherins.txt --samples coldata.txt >protocadherins_reads_nb_tpm.txt

# kangjingliang@kangjingliangdeMacBook-Pro 三 11 30 10:49:51 ~/Documents/2021/Cleaner_wrasse/gene_expression/Gene_families
less protocadherins_reads_nb_tpm.txt|perl -alne '$k++;print "Gene\tNbeq1\tNbeq1_ratio\tMaximun_Nb\tMean_Nb" if $k==1;next if /^\s+/;my ($j, $max, $mean);for (my $i = 1; $i < @F; $i++){$j++ if $F[$i]<1;$mean+=$F[$i];($F[$i]<$max)?($max=$max):($max=$F[$i])};($j>0)?($j=$j):($j=0);my $rat=$j/36;$rat=sprintf("%.3f",$rat);$mean=$mean/36;$mean=sprintf("%.3f",$mean);print "$F[0]\t$j\t$rat\t$max\t$mean"' >Protocadherins_Expression_info.txt
vi Protocadherins_expression_less_80_1.txt

# Expression TPM >=1 in 80% indviduals
perl Create_reads_nb_gene.pl --matr protocadherins_reads_nb_tpm.txt --gene PCDA2_1T,PCDA3_6C,PCDA6_3C,PCDA8_4C,PCDC2_1C,PCDC2_2C,PCDG2_1C,PCDG2_5T,PCDG3_2C,PCDG8_2C,PCDGM_2C,PCDGM_3C,PCDGM_5C,PCDGM_6C > protocadherins_reads_nb_tpm_plot.txt
```

## Plot specific gene reads nb
**TPM_plot.R**    
```bash
# Neurexins
# kangjingliang@kangjingliangdeMacBook-Pro 三 11 23 15:46:37 ~/Documents/2021/Cleaner_wrasse/gene_expression/enrichment
less Total_DEGs_social_GOs_uniq.txt|grep -i 'Neurexin'
perl Create_reads_nb_gene.pl --matr DEGs_social_reads_nb_tpm.txt --gene Ldim_g17130,Ldim_g17131,Ldim_g22189,Ldim_g25203,Ldim_g25204,Ldim_g557,Ldim_g7336 > Neurexin_tpm_plot.txt

# use the gene name
perl temp1.pl Total_DEGs_social_GOs_uniq.txt Neurexin_tpm_plot.txt > Neurexin_tpm_plot_update.txt

# according the print to mark "*"
perl temp2.pl Total_DEGs_social_GOs_uniq.txt Ldim_g17130 Ldim_g17131 Ldim_g22189 Ldim_g25203 Ldim_g25204 Ldim_g557 Ldim_g7336

# Glutamates involved in social behavior
# kangjingliang@kangjingliangdeMacBook-Pro 三 11 23 23:25:08 ~/Documents/2021/Cleaner_wrasse/gene_expression/enrichment
less Total_DEGs_social_GOs_uniq.txt|grep -i 'glutamate'|perl -alne '$in.=$F[0].",";END{$in=~s/\,$//;print $in}'
perl Create_reads_nb_gene.pl --matr DEGs_social_reads_nb_tpm.txt --gene Ldim_g10133,Ldim_g20661,Ldim_g21576,Ldim_g3641,Ldim_g7869 > Glutamate_social_tpm_plot.txt
perl temp1.pl Total_DEGs_social_GOs_uniq.txt Glutamate_social_tpm_plot.txt > Glutamate_social_tpm_plot_update.txt

less Total_DEGs_social_GOs_uniq.txt|grep -i 'glutamate'|perl -alne '$in.=$F[0]." ";END{$in=~s/\s+$//;print $in}'
perl temp2.pl Total_DEGs_social_GOs_uniq.txt Ldim_g10133 Ldim_g20661 Ldim_g21576 Ldim_g3641 Ldim_g7869

# Plot Neurexins and Glutamates involved in social behavior together
cat Neurexin_tpm_plot_update.txt Glutamate_social_tpm_plot_update.txt > Neurexin_Glutamate_social_tpm_plot_update.txt

# Glutamate receptor involved in sensory
less Total_DEGs_sensory_GOs_uniq.txt |grep -i 'glutamate receptor'|perl -alne '$in.=$F[0].",";END{$in=~s/\,$//;print $in}'
perl Create_reads_nb_gene.pl --matr DEGs_sensory_reads_nb_tpm.txt --gene Ldim_g10930,Ldim_g12359,Ldim_g12891,Ldim_g15257,Ldim_g16063,Ldim_g17696,Ldim_g17780,Ldim_g18806,Ldim_g19424,Ldim_g19799,Ldim_g5621,Ldim_g690 > Glutamate_sensory_tpm_plot.txt
perl temp1.pl Total_DEGs_sensory_GOs_uniq.txt Glutamate_sensory_tpm_plot.txt > Glutamate_sensory_tpm_plot_update.txt
less Total_DEGs_sensory_GOs_uniq.txt|grep -i 'glutamate receptor'|perl -alne '$in.=$F[0]." ";END{$in=~s/\s+$//;print $in}'
perl temp2.pl Total_DEGs_sensory_GOs_uniq.txt Ldim_g10930 Ldim_g12359 Ldim_g12891 Ldim_g15257 Ldim_g16063 Ldim_g17696 Ldim_g17780 Ldim_g18806 Ldim_g19424 Ldim_g19799 Ldim_g5621 Ldim_g690

# For all Glutamate receptors
# kangjingliang@kangjingliangdeMacBook-Pro 一 11 28 14:07:20 ~/Documents/2021/Cleaner_wrasse/gene_expression/enrichment
less Total_DEGs_uniq.txt|grep -i 'glutamate receptor'|perl -alne '$in.=$F[0].",";END{$in=~s/\,$//;print $in}'|less
perl temp2.pl Total_DEGs_uniq.txt Ldim_g10831 Ldim_g10930 Ldim_g12359 Ldim_g12891 Ldim_g13699 Ldim_g13736 Ldim_g15257 Ldim_g16055 Ldim_g16062 Ldim_g16063 Ldim_g17696 Ldim_g17780 Ldim_g18806 Ldim_g19110 Ldim_g19424 Ldim_g19799 Ldim_g20661 Ldim_g21258 Ldim_g21576 Ldim_g25970 Ldim_g3040 Ldim_g418 Ldim_g5526 Ldim_g5612 Ldim_g5621 Ldim_g6091 Ldim_g686 Ldim_g690 Ldim_g7121
perl Create_reads_nb_gene.pl --matr all_inds_read_nb_for_coefficient_tpm.txt --gene Ldim_g10831,Ldim_g10930,Ldim_g12359,Ldim_g12891,Ldim_g13699,Ldim_g13736,Ldim_g15257,Ldim_g16055,Ldim_g16062,Ldim_g16063,Ldim_g17696,Ldim_g17780,Ldim_g18806,Ldim_g19110,Ldim_g19424,Ldim_g19799,Ldim_g20661,Ldim_g21258,Ldim_g21576,Ldim_g25970,Ldim_g3040,Ldim_g418,Ldim_g5526,Ldim_g5612,Ldim_g5621,Ldim_g6091,Ldim_g686,Ldim_g690,Ldim_g7121  > All_glutamate_DEGs_tpm_plot.txt
perl temp1.pl Total_DEGs_uniq.txt All_glutamate_DEGs_tpm_plot.txt > All_glutamate_DEGs_tpm_plot_update.txt
less Total_DEGs_uniq.txt|grep -i 'glutamate receptor'|perl -alne '$in.=$F[0]." ";END{$in=~s/\,$//;print $in}'|less
perl temp2.pl Total_DEGs_uniq.txt Ldim_g10831 Ldim_g10930 Ldim_g12359 Ldim_g12891 Ldim_g13699 Ldim_g13736 Ldim_g15257 Ldim_g16055 Ldim_g16062 Ldim_g16063 Ldim_g17696 Ldim_g17780 Ldim_g18806 Ldim_g19110 Ldim_g19424 Ldim_g19799 Ldim_g20661 Ldim_g21258 Ldim_g21576 Ldim_g25970 Ldim_g3040 Ldim_g418 Ldim_g5526 Ldim_g5612 Ldim_g5621 Ldim_g6091 Ldim_g686 Ldim_g690 Ldim_g7121|less
```

```R
library(ggpubr)
setwd("~/Documents/2021/Cleaner_wrasse/gene_expression/enrichment")
data<-read.table("Neurexin_tpm_plot_update.txt", header = TRUE, sep = "\t")
names(data)
ggbarplot(data, x = "Tissue", y = "TPM", facet.by = "Gene", scales = "free",
          add = c("mean_se", "jitter"),font.family="Times",
          color = "Type", palette = c("#4DBBD5", "#E64B35"),
          position = position_dodge(0.8)) +
  theme(axis.text.x=element_text(colour="black",family="Times",size=16), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Times",size=16,face="plain",colour = "black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 20,face="plain"), #设置y轴标题的字体属性
        panel.border = element_blank(),axis.line = element_line(colour = "black"), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="plain", family="Times", colour="black", size=18), #设置图例的子标题的字体属性
        legend.title=element_blank(), #设置图例的总标题的字体属性
        axis.title.x=element_blank() )+
  theme(strip.text = element_text(face="plain", family="Times", colour="black", size=15),
        strip.background = element_rect(color = "white", fill = "white"))
```

## Combine the predicted gtf with the original gtf together for DEGs detection
```bash
# (base) kang1234@celia-PowerEdge-T640 Mon Feb 06 15:02:08 ~/genome/Gene_annotation/RNA-seq/RNA-align
cat ../../combined/Ldim_original_name.gtf predicted_genes_phy.gtf >Original_genewise_combine.gtf
featureCounts -a Original_genewise_combine.gtf -o Original_genewise_combine_reads_nb.txt -T 24 LD5FB.sorted.bam LD6FB.sorted.bam LD15FB.sorted.bam LD16FB.sorted.bam LD25FB.sorted.bam LD26FB.sorted.bam LD3FB.sorted.bam LD4FB.sorted.bam LD13FB.sorted.bam LD14FB.sorted.bam LD23FB.sorted.bam LD24FB.sorted.bam LD5HB.sorted.bam LD6HB.sorted.bam LD15HB.sorted.bam LD16HB.sorted.bam LD25HB.sorted.bam LD26HB.sorted.bam LD3HB.sorted.bam LD4HB.sorted.bam LD13HB.sorted.bam LD14HB.sorted.bam LD23HB.sorted.bam LD24HB.sorted.bam LD5MB.sorted.bam LD6MB.sorted.bam LD15MB.sorted.bam LD16MB.sorted.bam LD25MB.sorted.bam LD26MB.sorted.bam LD3MB.sorted.bam LD4MB.sorted.bam LD13MB.sorted.bam LD14MB.sorted.bam LD23MB.sorted.bam LD24MB.sorted.bam
# (base) kang1234@celia-PowerEdge-T640 Mon Feb 06 15:19:16 ~/genome/Gene_annotation/RNA-seq/RNA-align
cp Original_genewise_combine_reads_nb.txt read_matrix/
# (base) kang1234@celia-PowerEdge-T640 Mon Feb 06 15:20:56 ~/genome/Gene_annotation/RNA-seq/RNA-align/read_matrix
perl Extract_ind_reads_nb.pl Original_genewise_combine_reads_nb.txt
# kangjingliang@kangjingliangdeMacBook-Pro 一  2 06 15:23:36 ~/Documents/2023/genome/DEGs
scp kang1234@147.8.76.177:~/genome/Gene_annotation/RNA-seq/RNA-align/read_matrix/Original_genewise_combine_reads_nb_*.txt ./
# kangjingliang@kangjingliangdeMacBook-Pro 一  2 06 15:26:14 ~/Documents/2023/genome/DEGs
cp ~/Documents/2021/Cleaner_wrasse/gene_expression/coldata_*.txt ./

# FB: 2747 DEGs
DESeq --matrix Original_genewise_combine_reads_nb_FB.txt --samples coldata_FB.txt --column Group --prefix FB

# MB: 512 DEGs
DESeq --matrix Original_genewise_combine_reads_nb_MB.txt --samples coldata_MB.txt --column Group --prefix MB

# HB: 1559 DEGs
DESeq --matrix Original_genewise_combine_reads_nb_HB.txt --samples coldata_HB.txt --column Group --prefix HB
```
