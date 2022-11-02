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
