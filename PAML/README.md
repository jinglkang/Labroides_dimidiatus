# PAML: evolution analysis
## swiss-prot blastp results
```bash
# SNORLAX
~/genome/Gene_annotation/swis-blastp.result # six reference fish spceis
~/genome/Gene_annotation/combined/Ldim_blastp.result.best.txt # Labroides_dimidiatus
~/genome/Gene_annotation/Cheilinus_undulatus/swis-blastp.result # Cheilinus_undulatus
~/genome/Gene_annotation/Labrus_bergylta/swis-blastp.result # Labrus_bergylta
~/genome/Gene_annotation/Notolabrus_celidotus/swis-blastp.result # Notolabrus_celidotus
~/genome/Gene_annotation/Semicossyphus_pulcher/swis-blastp.result # Semicossyphus_pulcher
~/genome/Gene_annotation/Symphodus_melops/swis-blastp.result # Symphodus_melops
~/genome/Gene_annotation/Tautogolabrus_adspersus/swis-blastp.result # Tautogolabrus_adspersus
~/genome/Gene_annotation/Thalassoma_bifasciatum/swis-blastp.result # Thalassoma_bifasciatum
```

## cat all blastp results together
```bash
# (base) kang1234@celia-PowerEdge-T640 Tue Nov 01 16:24:58 ~/genome/Gene_annotation
cat ~/genome/Gene_annotation/swis-blastp.result \
~/genome/Gene_annotation/combined/Ldim_blastp.result.best.txt \
~/genome/Gene_annotation/Cheilinus_undulatus/swis-blastp.result \
~/genome/Gene_annotation/Labrus_bergylta/swis-blastp.result \
~/genome/Gene_annotation/Notolabrus_celidotus/swis-blastp.result \
~/genome/Gene_annotation/Semicossyphus_pulcher/swis-blastp.result \
~/genome/Gene_annotation/Symphodus_melops/swis-blastp.result \
~/genome/Gene_annotation/Tautogolabrus_adspersus/swis-blastp.result \
~/genome/Gene_annotation/Thalassoma_bifasciatum/swis-blastp.result > allspe_blastp.result

# Kang@fishlab3 Tue Nov 01 16:26:04 /media/HDD/cleaner_fish/genome/gene_family_3/OrthoFinder/Results_May09/sub_orth
scp kang1234@147.8.76.177:~/genome/Gene_annotation/allspe_blastp.result ./
# Extract gene list for PAML
# 1. keep the genes in the group have the same gene name
# 2. keep the one has the highest blast score if the species has more than two genes in this group
# Kang@fishlab3 Wed Nov 02 09:52:18 /media/HDD/cleaner_fish/genome/gene_family_3/OrthoFinder/Results_May09/sub_orth
perl Generate_genelist_paml.pl > Genelist_paml.txt
# Extract the annotation for each gene list
scp kang1234@147.8.76.177:~/genome/Gene_annotation/combined/braker2+3_combined_renamed.aa.long.anno.final.txt Ldim_ano.txt
perl Extract_ano.pl > Genelist_ano.txt

mkdir paml_input; cd paml_input
# Kang@fishlab3 Wed Nov 02 09:58:28 /media/HDD/cleaner_fish/genome/gene_family_3/OrthoFinder/Results_May09/sub_orth/paml_input
cp /media/HDD/cleaner_fish/genome/gene_family_2/paml_input/prepare_input_paml.pl ./
cp /media/HDD/cleaner_fish/genome/gene_family_2/paml_input/prepare_input_paml_parallel.pl ./
cp /media/HDD/cleaner_fish/genome/gene_family_3/OrthoFinder/Results_May09/paml_input/correlation.txt ./
cp /media/HDD/cleaner_fish/genome/gene_family_3/OrthoFinder/Results_May09/paml_input/*.fasta ./
nohup perl prepare_input_paml_parallel.pl Genelist_paml.txt >prepare_input_paml.process 2>&1 &
# [1] 23011
```
## Run PAML
### 1. free-ratio
```bash
# Kang@fishlab3 Thu Nov 03 23:29:23 /media/HDD/cleaner_fish/genome/gene_family_3/OrthoFinder/Results_May09/sub_orth/paml_input
cp /media/HDD/cleaner_fish/genome/gene_family_2/paml_input/codeml.pl ./
cp /media/HDD/cleaner_fish/genome/gene_family_2/paml_input/codeml_parallel.pl ./
vi spe.tre
# ((((Fugu,(Stickleback,(Spul,((Cund,((Smel,Tads),Lber)),(Ncel,(Ldim,Tbif)))))),(Platyfish,Medaka)),Zebrafish),Spottedgar);
# perl codeml.pl --input final_orth_input_paml.txt --model free-ratio --dir . --tree spe.tre --icode 0 --omega 1.2
# revise the command in codeml_parallel.pl
# "my $cmd="perl codeml.pl --input temp/$temp --model free-ratio --dir . --tree spe.tre --icode 0 --omega 1.2";"
nohup perl codeml_parallel.pl final_orth_input_paml.txt >free_ratio.process 2>&1 &
# [1] 31944
```
### 2. Positive selection
#### 2.1 Ldim
```bash
# (base) kang1234@celia-PowerEdge-T640 Thu Nov 03 23:50:11 ~/genome
mkdir paml_new; cd paml_new
# Kang@fishlab3 Thu Nov 03 23:47:32 /media/HDD/cleaner_fish/genome/gene_family_3/OrthoFinder/Results_May09/sub_orth
nohup scp -r paml_input/ kang1234@147.8.76.177:~/genome/paml_new/
# ctrl+z
bg
# [2]+ nohup scp -r paml_input/ kang1234@147.8.76.177:~/genome/paml_new/ &
# (base) kang1234@celia-PowerEdge-T640 Fri Nov 04 00:08:04 ~/genome/paml_new/paml_input
vi spe_Ldim.tre
# ((((Fugu,(Stickleback,(Spul,((Cund,((Smel,Tads),Lber)),(Ncel,(Ldim #1,Tbif)))))),(Platyfish,Medaka)),Zebrafish),Spottedgar);
# perl codeml.pl --input temp/$temp --model branch-site --dir . --output_suf Ldim --tree spe_Ldim.tre --icode 0 --omega 1.2
nohup perl codeml_parallel.pl final_orth_input_paml.txt >codeml.process 2>&1 &
# [1] 3530
# (base) kang1234@celia-PowerEdge-T640 Mon Nov 07 11:02:09 ~/genome/paml_new/paml_input
perl Extract_PSGs.pl Genelist_ano.txt > Ldim_PSGs.txt # Ldim has 162 PSGs
```
### 2.2 All Labridae fish species
```bash
# (base) kang1234@celia-PowerEdge-T640 Mon Nov 07 11:02:44 ~/genome/paml_new/paml_input
vi Ancestor_wrasses.tre
# ((((Fugu,(Stickleback,(Spul,((Cund,((Smel,Tads),Lber)),(Ncel,(Ldim,Tbif)))) #1)),(Platyfish,Medaka)),Zebrafish),Spottedgar);

# revise codeml_parallel.pl
# "my $cmd="perl codeml.pl --input temp/$temp --model branch-site --dir . --output_suf Ldim --tree spe_Ldim.tre --icode 0 --omega 1.2";" =>
# "my $cmd="perl codeml.pl --input temp/$temp --model branch-site --dir . --output_suf Ancestor_wrasses --tree Ancestor_wrasses.tre --icode 0 --omega 1.2";"
nohup perl codeml_parallel.pl final_orth_input_paml.txt >codeml.process 2>&1 &
# [1] 30809
```

## 3. Plot the pep sequences (PSGs in Ldim)
### GRIA3
```bash
# OG0000065_OG8   GRIA3   Glutamate receptor 3
# (base) kang1234@celia-PowerEdge-T640 Mon Nov 07 11:52:07 ~/genome/paml_new/paml_input/OG0000065_OG8
cds2pep.pl final_alignment.fa > final_alignment_pep.fa
# kangjingliang@kangjingliangdeMacBook-Pro 一 11 07 11:53:11 ~/Desktop
scp kang1234@147.8.76.177:~/genome/paml_new/paml_input/OG0000065_OG8/final_alignment_pep.fa GRIA3_alignment_pep.fa
```
### BMPR2
```bash
# OG0013224       BMPR2   Bone morphogenetic protein receptor type-2
# (base) kang1234@celia-PowerEdge-T640 Mon Nov 07 12:40:19 ~/genome/paml_new/paml_input/OG0013224
cds2pep.pl final_alignment.fa > final_alignment_pep.fa
# kangjingliang@kangjingliangdeMacBook-Pro 一 11 07 12:42:20 ~/Desktop
scp kang1234@147.8.76.177:~/genome/paml_new/paml_input/OG0013224/final_alignment_pep.fa BMPR2_alignment_pep.fa
```
### TLR5
```bash
# OG0001476_OG1	Tlr5	Toll-like receptor 5
# (base) kang1234@celia-PowerEdge-T640 Wed Nov 09 16:29:23 ~/genome/paml_new/paml_input/OG0001476_OG1
cds2pep.pl final_alignment.fa > final_alignment_pep.fa
# kangjingliang@kangjingliangdeMacBook-Pro 三 11 09 16:39:58 ~/Documents/2022/Ldim_genome_Restart/PSGs
scp kang1234@147.8.76.177:~/genome/paml_new/paml_input/OG0001476_OG1/final_alignment_pep.fa ./TLR5_alignment_pep.fa
```
### cyp26b1
```bash
# OG0008849_OG0	cyp26b1	Cytochrome P450 26B1
# (base) kang1234@celia-PowerEdge-T640 Wed Nov 09 16:47:15 ~/genome/paml_new/paml_input/OG0008849_OG0
cds2pep.pl final_alignment.fa > final_alignment_pep.fa
# kangjingliang@kangjingliangdeMacBook-Pro 三 11 09 16:40:59 ~/Documents/2022/Ldim_genome_Restart/PSGs
scp kang1234@147.8.76.177:~/genome/paml_new/paml_input/OG0008849_OG0/final_alignment_pep.fa ./cyp26b1_alignment_pep.fa
```
## 4. Reconstruct the sequences on each ancestral node
```bash
# (base) kang1234@celia-PowerEdge-T640 Thu Nov 10 10:09:11 ~/genome/paml_new/paml_input
cp OG0000065_OG8/1.ctr Ancestral.crl
# (base) kang1234@celia-PowerEdge-T640 Thu Nov 10 11:23:59 ~/genome/paml_new/paml_input
nohup perl Rebuild_seq_ancetral_parallel.pl final_orth_input_paml.txt Ancestral.crl >Rebuild_seq_ancetral.process 2>&1 &
# [1] 26271
```
