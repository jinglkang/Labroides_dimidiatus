# The other gene families: Glutamate; Dopamine; Crystallins; Hox; BMPs; Cytochrome P450 enzymes; Toll-like receptors; RIG-like receptors; NLRs; Taste Receptors; Immunoglobulin; Immediate early response genes
## Glutamate receptor prediction
```bash
# (base) kang1234@celia-PowerEdge-T640 Fri Oct 21 15:07:57 ~/genome/Gene_annotation
perl Create_query_seqs.pl Glutamate_names.txt >Query_glutamate.fasta
# Kang@fishlab3 Fri Oct 21 15:28:34 /media/HDD/cleaner_fish/genome
mkdir Glutamates; cd Glutamates
# Kang@fishlab3 Fri Oct 21 15:30:29 /media/HDD/cleaner_fish/genome/Glutamates
scp kang1234@147.8.76.177:~/genome/Gene_annotation/Query_glutamate.fasta ./
scp kang1234@147.8.76.177:~/genome/Gene_annotation/Glutamate_names.txt ./
cp ../Opsin_new/run_Fmdetect.pl ./
cp ../Opsin_new/Fmdetect.pl ./
# cover ratio set to 50%, remove other limitations, just confirm the name is same and cover ratio >= 0.5
cp ../TLRs/Fmdetect.pl ./
nohup perl run_Fmdetect.pl /media/HDD/cleaner_fish/genome/Glutamates/Query_glutamate.fasta 50 /media/HDD/cleaner_fish/genome/Glutamates/Glutamate_names.txt >Detect_glutamate.process 2>&1 &
#[1] 15741
# Build the phylogenetic tree for glutamate receptors
# Kang@fishlab3 Sat Oct 22 12:00:22 /media/HDD/cleaner_fish/genome/Glutamates
cp /media/HDD/cleaner_fish/genome/Protocadherin_new/temp6.pl ./
perl temp1.pl
mv Predict_Opsins.fa Predict_Glutamates.fa
mv Predict_Opsins.phy Predict_Glutamates.phy
# (base) kang1234@celia-PowerEdge-T640 Sat Oct 22 12:15:55 ~/genome/gene_family
mkdir Glutamates; cd Glutamates
# Kang@fishlab3 Sat Oct 22 12:15:14 /media/HDD/cleaner_fish/genome/Glutamates
scp Predict_Glutamates.phy kang1234@147.8.76.177:~/genome/gene_family/Glutamates
# (base) kang1234@celia-PowerEdge-T640 Sat Oct 22 23:44:56 ~/genome/gene_family/Glutamates
nohup raxmlHPC -T 24 -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s Predict_Glutamates.phy -n Glutamates >Glutamates.process 2>&1 &
# [1] 20600
# Kang@fishlab3 Mon Oct 24 11:08:13 /media/HDD/cleaner_fish/genome/Glutamates
cp /media/HDD/cleaner_fish/genome/Opsin_new/label_color.pl ./
cp /media/HDD/cleaner_fish/genome/Opsin_new/node_symbol.pl ./
# Kang@fishlab3 Mon Oct 24 11:09:27 /media/HDD/cleaner_fish/genome/Glutamates
perl label_color.pl Predict_Glutamates.fa|less
perl node_symbol.pl Predict_Glutamates.fa|less
less Predict_Glutamates.phy|perl -alne 'next if /^\d+/;if ($F[0]=~/(.*?)_(.*?)_.*/){$hash{$1}++ if ($2 eq "GRM1" || $2 eq "GRM5")};END{foreach $sp (qw(Ldim Tbif Smel Tads Spul Ncel Lber Cund)){print "$sp\t$hash{$sp}"}}'
less Predict_Glutamates.phy|perl -alne 'next if /^\d+/;unless ($F[0]=~/(.*?)_(.*?)_.*/ && $2 eq "PCDGM"){$hash{$1}++};END{foreach $sp (qw(Ldim Tbif Smel Tads Spul Ncel Lber Cund)){print "$sp\t$hash{$sp}"}}'
```

## Dopamine receptors
```bash
# (base) kang1234@celia-PowerEdge-T640 Sat Oct 22 16:59:17 ~/genome/Gene_annotation
less ref.fasta.ano.final|grep -i 'Dopamine'|perl -alne '($nm)=$_=~/Name=\"(.*)\s+\[.*\]\"/;print $nm'|sort -u
perl Create_query_seqs.pl Dopamines_names.txt >Query_Dopamines.fasta
# Kang@fishlab3 Sat Oct 22 17:09:17 /media/HDD/cleaner_fish/genome
mkdir Dopamines; cd Dopamines
# Kang@fishlab3 Sat Oct 22 17:10:16 /media/HDD/cleaner_fish/genome/Dopamines
scp kang1234@147.8.76.177:~/genome/Gene_annotation/Query_Dopamines.fasta ./
scp kang1234@147.8.76.177:~/genome/Gene_annotation/Dopamines_names.txt ./
cp ../Glutamates/Fmdetect.pl ./
cp ../Glutamates/run_Fmdetect.pl ./
# cover ratio set to 50%, remove other limitations, just confirm the name is same and cover ratio >= 0.5
cp ../TLRs/Fmdetect.pl ./
nohup perl run_Fmdetect.pl /media/HDD/cleaner_fish/genome/Dopamines/Query_Dopamines.fasta 50 /media/HDD/cleaner_fish/genome/Dopamines/Dopamines_names.txt >Detect_Dopamines.process 2>&1 &
# [1] 29258
# Kang@fishlab3 Sat Oct 22 23:41:02 /media/HDD/cleaner_fish/genome/Dopamines
cp /media/HDD/cleaner_fish/genome/Protocadherin_new/temp6.pl ./
cp /media/HDD/cleaner_fish/genome/Glutamates/temp1.pl ./
perl temp1.pl
mv Predict_Opsins.fa Predict_Dopamines.fa
mv Predict_Opsins.phy Predict_Dopamines.phy
scp Predict_Dopamines.phy kang1234@147.8.76.177:~/genome/gene_family/Glutamates
# (base) kang1234@celia-PowerEdge-T640 Sat Oct 22 23:44:56 ~/genome/gene_family/Glutamates
nohup raxmlHPC -T 24 -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s Predict_Dopamines.phy -n Dopamines >Dopamines.process 2>&1 &
# [1] 25251
cp ../NLRs/temp2.pl ./
less Predict_Dopamines.phy|perl -alne 'next if /^\d+/;if ($F[0]=~/(.*?)_(.*?)_.*/){$hash{$2}->{$1}++};END{foreach my $ge (sort keys %hash){my $info; foreach $sp (qw(Ldim Tbif Smel Tads Spul Ncel Lber Cund)){my $nb; $hash{$ge}->{$sp}?($nb=$hash{$ge}->{$sp}):($nb=0);$info.=$nb.":";};$info=~s/\:$//;print "$ge\t[$info]"}}' >Dopamines_nb.txt
perl temp2.pl Dopamines_nb.txt > Dopamines_nb_2.txt
```

## Crystallins
```bash
# (base) kang1234@celia-PowerEdge-T640 Sat Oct 22 16:59:17 ~/genome/Gene_annotation
less ref.fasta.ano.final|grep -i 'crystallin'|perl -alne '($nm)=$_=~/Name=\"(.*)\s+\[.*\]\"/;print $nm'|sort -u
perl Create_query_seqs.pl Crystallins_names.txt > Query_Crystallins.fasta
# Kang@fishlab3 Sat Oct 22 17:02:09 /media/HDD/cleaner_fish/genome
mkdir Crystallins;cd Crystallins
# Kang@fishlab3 Mon Oct 24 00:02:01 /media/HDD/cleaner_fish/genome/Crystallins
scp kang1234@147.8.76.177:~/genome/Gene_annotation/Query_Crystallins.fasta ./
scp kang1234@147.8.76.177:~/genome/Gene_annotation/Crystallins_names.txt ./
cp ../Glutamates/Fmdetect.pl ./
cp ../Glutamates/run_Fmdetect.pl ./
# cover ratio set to 50%, remove other limitations, just confirm the name is same and cover ratio >= 0.5
cp ../TLRs/Fmdetect.pl ./
nohup perl run_Fmdetect.pl /media/HDD/cleaner_fish/genome/Crystallins/Query_Crystallins.fasta 50 /media/HDD/cleaner_fish/genome/Crystallins/Crystallins_names.txt >Detect_Crystallins.process 2>&1 &
# [1] 20378
# Build the phylogenetic tree
cp /media/HDD/cleaner_fish/genome/Protocadherin_new/temp6.pl ./
cp /media/HDD/cleaner_fish/genome/Glutamates/temp1.pl ./
perl temp1.pl
mv Predict_Opsins.fa Predict_Crystallins.fa
mv Predict_Opsins.phy Predict_Crystallins.phy
scp Predict_Crystallins.phy kang1234@147.8.76.177:~/genome/gene_family/Glutamates
# (base) kang1234@celia-PowerEdge-T640 Mon Oct 24 10:06:06 ~/genome/gene_family/Glutamates
nohup raxmlHPC -T 24 -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s Predict_Crystallins.phy -n Crystallins >Crystallins.process 2>&1 &
# [1] 32061
# Kang@fishlab3 Mon Oct 24 14:56:22 /media/HDD/cleaner_fish/genome/Crystallins
cp /media/HDD/cleaner_fish/genome/Opsin_new/label_color.pl ./
cp /media/HDD/cleaner_fish/genome/Opsin_new/node_symbol.pl ./
less Alpha_allspe.phy|perl -alne 'next if /^\d+/;if ($F[0]=~/(.*?)_(.*?)_.*/ && $2 eq "PCDC2"){$hash{$1}++};END{foreach $sp (qw(Ldim Tbif Smel Tads Spul Ncel Lber Cund)){print "$sp\t$hash{$sp}"}}'
```


## Hox genes
```bash
# (base) kang1234@celia-PowerEdge-T640 Sun Oct 23 23:48:29 ~/genome/Gene_annotation
less ref.fasta.ano.final|grep -i 'hox'|perl -alne '($nm)=$_=~/Name=\"(.*)\s+\[.*\]\"/;print $nm'|sort -u
perl Create_query_seqs.pl Hox_names.txt > Query_Hox.fasta
# Kang@fishlab3 Sun Oct 23 23:55:59 /media/HDD/cleaner_fish/genome
mkdir Hox; cd Hox
# Kang@fishlab3 Sun Oct 23 23:57:50 /media/HDD/cleaner_fish/genome/Hox
scp kang1234@147.8.76.177:~/genome/Gene_annotation/Query_Hox.fasta ./
scp kang1234@147.8.76.177:~/genome/Gene_annotation/Hox_names.txt ./
cp ../Glutamates/Fmdetect.pl ./
cp ../Glutamates/run_Fmdetect.pl ./
# cover ratio set to 50%, remove other limitations, just confirm the name is same and cover ratio >= 0.5
cp ../TLRs/Fmdetect.pl ./
nohup perl run_Fmdetect.pl /media/HDD/cleaner_fish/genome/Hox/Query_Hox.fasta 50 /media/HDD/cleaner_fish/genome/Hox/Hox_names.txt >Detect_Hox.process 2>&1 &
# [1] 17828
cp /media/HDD/cleaner_fish/genome/Protocadherin_new/temp6.pl ./
cp /media/HDD/cleaner_fish/genome/Glutamates/temp1.pl ./
perl temp1.pl
mv Predict_Opsins.fa Predict_Hox.fa
mv Predict_Opsins.phy Predict_Hox.phy
scp Predict_Hox.phy kang1234@147.8.76.177:~/genome/gene_family/Glutamates
# (base) kang1234@celia-PowerEdge-T640 Mon Oct 24 10:06:06 ~/genome/gene_family/Glutamates
nohup raxmlHPC -T 22 -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s Predict_Hox.phy -n Hox >Hox.process 2>&1 &
# [1] 16875
# Check the gene number
# Kang@fishlab3 Mon Oct 24 18:13:42 /media/HDD/cleaner_fish/genome/Hox
less Predict_Hox.phy|perl -alne 'next if /^\d+/;if ($F[0]=~/(.*?)_(.*?)_.*/){$hash{$2}->{$1}++};END{foreach my $ge (sort keys %hash){foreach $sp (qw(Ldim Tbif Smel Tads Spul Ncel Lber Cund)){print "$ge\t$sp\t$hash{$ge}->{$sp}"}print}}'|less
# Kang@fishlab3 Wed Oct 26 10:29:07 /media/HDD/cleaner_fish/genome/Hox
less Predict_Hox.phy|perl -alne 'next if /^\d+/;if ($F[0]=~/(.*?)_(.*?)_.*/){$hash{$2}->{$1}++};END{foreach my $ge (sort keys %hash){my $info; foreach $sp (qw(Ldim Tbif Smel Tads Spul Ncel Lber Cund)){my $nb; $hash{$ge}->{$sp}?($nb=$hash{$ge}->{$sp}):($nb=0);$info.=$nb.":";};$info=~s/\:$//;print "$ge\t[$info]"}}' >Hox_nb.txt
# Kang@fishlab3 Wed Oct 26 11:20:16 /media/HDD/cleaner_fish/genome/Hox
perl temp2.pl > Hox_nb_2.txt
```


## BMPs (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2857416/)
```bash
# Embryos infected with RCAS-Bmp-2 or RCAS-Bmp-4 at this early time point exhibited severe facial malformations at day 9. The upper and lower components of the jaw were short and clefts between the maxillary and frontonasal process were evident
# (base) kang1234@celia-PowerEdge-T640 Mon Oct 24 10:51:38 ~/genome/Gene_annotation
less ref.fasta.ano.final|grep -i 'bone'|perl -alne '($nm)=$_=~/Name=\"(.*)\s+\[.*\]\"/;print $nm'|sort -u >BMPs_names.txt
perl Create_query_seqs.pl BMPs_names.txt > Query_BMPs.fasta
# Kang@fishlab3 Mon Oct 24 10:58:40 /media/HDD/cleaner_fish/genome
mkdir BMPs; cd BMPs
# Kang@fishlab3 Mon Oct 24 11:00:26 /media/HDD/cleaner_fish/genome/BMPs
scp kang1234@147.8.76.177:~/genome/Gene_annotation/Query_BMPs.fasta ./
scp kang1234@147.8.76.177:~/genome/Gene_annotation/BMPs_names.txt ./
cp ../Glutamates/Fmdetect.pl ./
cp ../Glutamates/run_Fmdetect.pl ./
# cover ratio set to 50%, remove other limitations, just confirm the name is same and cover ratio >= 0.5
cp ../TLRs/Fmdetect.pl ./
nohup perl run_Fmdetect.pl /media/HDD/cleaner_fish/genome/BMPs/Query_BMPs.fasta 50 /media/HDD/cleaner_fish/genome/BMPs/BMPs_names.txt >Detect_BMPs.process 2>&1 &
# [1] 14086
cp /media/HDD/cleaner_fish/genome/Protocadherin_new/temp6.pl ./
cp /media/HDD/cleaner_fish/genome/Glutamates/temp1.pl ./
perl temp1.pl
mv Predict_Opsins.fa Predict_BMPs.fa
mv Predict_Opsins.phy Predict_BMPs.phy
scp Predict_BMPs.phy kang1234@147.8.76.177:~/genome/gene_family/Glutamates
# (base) kang1234@celia-PowerEdge-T640 Mon Oct 24 10:06:06 ~/genome/gene_family/Glutamates
nohup raxmlHPC -T 22 -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s Predict_BMPs.phy -n BMPs >BMPs.process 2>&1 &
# iTol plot
# Kang@fishlab3 Fri Oct 28 10:29:15 /media/HDD/cleaner_fish/genome/BMPs
perl ../Opsin_new/label_color.pl Predict_BMPs.fa|less
cp ../NLRs/temp2.pl ./
less Predict_BMPs.phy|perl -alne 'next if /^\d+/;if ($F[0]=~/(.*?)_(.*?)_.*/){$hash{$2}->{$1}++};END{foreach my $ge (sort keys %hash){my $info; foreach $sp (qw(Ldim Tbif Smel Tads Spul Lber Ncel Cund)){my $nb; $hash{$ge}->{$sp}?($nb=$hash{$ge}->{$sp}):($nb=0);$info.=$nb.":";};$info=~s/\:$//;print "$ge\t[$info]"}}' >BMPs_nb.txt
perl temp2.pl BMPs_nb.txt > BMPs_nb_2.txt
```

## Cytochrome P450 enzymes
```bash
# (base) kang1234@celia-PowerEdge-T640 Mon Oct 24 21:52:21 ~/genome/Gene_annotation
less ref.fasta.ano.final|grep -i 'cytochrome p450'|perl -alne '($nm)=$_=~/Name=\"(.*)\s+\[.*\]\"/;print $nm'|sort -u|less
perl Create_query_seqs.pl Cytochrome_names.txt > Query_Cytochrome.fasta
# Kang@fishlab3 Mon Oct 24 21:49:10 /media/HDD/cleaner_fish/genome
mkdir Cytochrome; cd Cytochrome
# Kang@fishlab3 Mon Oct 24 21:59:07 /media/HDD/cleaner_fish/genome/Cytochrome
scp kang1234@147.8.76.177:~/genome/Gene_annotation/Query_Cytochrome.fasta ./
scp kang1234@147.8.76.177:~/genome/Gene_annotation/Cytochrome_names.txt ./
cp ../Glutamates/Fmdetect.pl ./
cp ../Glutamates/run_Fmdetect.pl ./
# cover ratio set to 50%, remove other limitations, just confirm the name is same and cover ratio >= 0.5
cp ../TLRs/Fmdetect.pl ./
nohup perl run_Fmdetect.pl /media/HDD/cleaner_fish/genome/Cytochrome/Query_Cytochrome.fasta 50 /media/HDD/cleaner_fish/genome/Cytochrome/Cytochrome_names.txt >Detect_Cytochrome.process 2>&1 &
# [1] 13901
cp /media/HDD/cleaner_fish/genome/Protocadherin_new/temp6.pl ./
cp /media/HDD/cleaner_fish/genome/Glutamates/temp1.pl ./
perl temp1.pl
mv Predict_Opsins.fa Predict_CP450.fa
mv Predict_Opsins.phy Predict_CP450.phy
scp Predict_CP450.phy kang1234@147.8.76.177:~/genome/gene_family/Glutamates
# (base) kang1234@celia-PowerEdge-T640 Mon Oct 24 10:06:06 ~/genome/gene_family/Glutamates
nohup raxmlHPC -T 22 -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s Predict_CP450.phy -n CP450 >CP450.process 2>&1 &
# iTol plot
# Kang@fishlab3 Fri Oct 28 17:44:32 /media/HDD/cleaner_fish/genome/Cytochrome
perl ../Opsin_new/label_color.pl Predict_CP450.fa|less
cp ../Dopamines/temp2.pl ./
less Predict_CP450.phy|perl -alne 'next if /^\d+/;if ($F[0]=~/(.*?)_(.*?)_.*/){$hash{$2}->{$1}++};END{foreach my $ge (sort keys %hash){my $info; foreach $sp (qw(Ldim Tbif Smel Tads Spul Ncel Lber Cund)){my $nb; $hash{$ge}->{$sp}?($nb=$hash{$ge}->{$sp}):($nb=0);$info.=$nb.":";};$info=~s/\:$//;print "$ge\t[$info]"}}' >CP450_nb.txt
perl temp2.pl CP450_nb.txt > CP450_nb_2.txt
```

## Toll-like receptors
```bash
# (base) kang1234@celia-PowerEdge-T640 Mon Oct 24 22:01:32 ~/genome/Gene_annotation
less ref.fasta.ano.final|grep -i 'toll-like'|perl -alne '($nm)=$_=~/Name=\"(.*)\s+\[.*\]\"/;print $nm'|sort -u >TLRs_names.txt
perl Create_query_seqs.pl TLRs_names.txt > Query_TLRs.fasta
# Kang@fishlab3 Mon Oct 24 23:02:19 /media/HDD/cleaner_fish/genome
mkdir TLRs; cd TLRs
# Kang@fishlab3 Mon Oct 24 23:02:37 /media/HDD/cleaner_fish/genome/TLRs
scp kang1234@147.8.76.177:~/genome/Gene_annotation/Query_TLRs.fasta ./
scp kang1234@147.8.76.177:~/genome/Gene_annotation/TLRs_names.txt ./
cp ../Glutamates/Fmdetect.pl ./
cp ../Glutamates/run_Fmdetect.pl ./
# cover ratio set to 50%, remove other limitations, just confirm the name is same and cover ratio >= 0.5
nohup perl run_Fmdetect.pl /media/HDD/cleaner_fish/genome/TLRs/Query_TLRs.fasta 50 /media/HDD/cleaner_fish/genome/TLRs/TLRs_names.txt >Detect_TLRs.process 2>&1 &
# [1] 11820
cp /media/HDD/cleaner_fish/genome/Protocadherin_new/temp6.pl ./
cp /media/HDD/cleaner_fish/genome/Glutamates/temp1.pl ./
perl temp1.pl
mv Predict_Opsins.fa Predict_TLRs.fa
mv Predict_Opsins.phy Predict_TLRs.phy
scp Predict_TLRs.phy kang1234@147.8.76.177:~/genome/gene_family/Glutamates
# (base) kang1234@celia-PowerEdge-T640 Mon Oct 24 10:06:06 ~/genome/gene_family/Glutamates
nohup raxmlHPC -T 22 -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s Predict_TLRs.phy -n TLRs >TLRs.process 2>&1 &
```

## RIG-like receptors
```bash
# RLRs: retinoic acid‐inducible gene I (RIG‐I) ‐like receptors (RLRs)
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5382327/
# Run in SNORLAX
# Kang@fishlab3 Tue Oct 25 15:33:46 ~/software/bin
scp *.pl kang1234@147.8.76.177:~/software/bin
scp solar kang1234@147.8.76.177:~/software/bin
# Kang@fishlab3 Tue Oct 25 15:36:22 /media/HDD/cleaner_fish/genome/TLRs
scp /media/HDD/cleaner_fish/genome/All_genomes/*.fasta kang1234@147.8.76.177:~/genome/gene_family/gene_family_prediction
# (base) kang1234@celia-PowerEdge-T640 Tue Oct 25 15:35:03 ~/genome/gene_family
mkdir gene_family_prediction; cd gene_family_prediction
# (base) kang1234@celia-PowerEdge-T640 Tue Oct 25 15:58:38 ~/genome/gene_family/gene_family_prediction
mkdir All_genomes; mv *.fasta All_genomes/
# uniprot: /home/kang1234/swiss-prot/uniprot-filtered-reviewed_yes.fasta
# genome: /home/kang1234/genome/gene_family/gene_family_prediction/All_genomes
# Kang@fishlab3 Thu Oct 27 10:57:50 /media/HDD/cleaner_fish/genome
mkdir RLRs; cd RLRs
# kangjingliang@kangjingliangdeMacBook-Pro 四 10 27 11:00:25 ~/Documents/2022/Ldim_genome_Restart/RLRs
scp Query_RLRs.fasta Kang@147.8.76.231:/media/HDD/cleaner_fish/genome/RLRs
blastp -outfmt 6 -query Query_RLRs.fasta -out Query_RLRs.bla -db ~/Desktop/Annotation_database/swiss-prot/uniprot-filtered-reviewed_yes.fasta -num_threads 30
cat Query_RLRs.bla |perl -lane 'print if $F[0] ne $flag; $flag=$F[0]'|perl -aF'\t' -lne '$h=join "\t", @F;print $h' |sort -k2 >Query_RLRs.bla.best
# Kang@fishlab3 Thu Oct 27 11:12:19 /media/HDD/cleaner_fish/genome/RLRs
less Query_RLRs.bla.best|perl -alne '(my $ne)=$F[1]=~/sp\|.*\|(.*)_.*/;print $ne'|sort -u
# DDX58
# DHX58
# IFIH1
# Search the name in swiss-prot and add the gene description in RLRs_names.txt
cp ../NLRs/Fmdetect.pl ./
cp ../NLRs/run_Fmdetect.pl ./
nohup perl run_Fmdetect.pl /media/HDD/cleaner_fish/genome/RLRs/Query_RLRs.fasta 50 /media/HDD/cleaner_fish/genome/RLRs/RLRs_names.txt >Detect_RLRs.process 2>&1 &
# [1] 6166
cp /media/HDD/cleaner_fish/genome/Protocadherin_new/temp6.pl ./
cp /media/HDD/cleaner_fish/genome/Glutamates/temp1.pl ./
perl temp1.pl
mv Predict_Opsins.fa Predict_RLRs.fa
mv Predict_Opsins.phy Predict_RLRs.phy
# Kang@fishlab3 Fri Oct 28 13:33:41 /media/HDD/cleaner_fish/genome/RLRs
scp Predict_RLRs.phy kang1234@147.8.76.177:~/genome/gene_family/Glutamates
nohup raxmlHPC -T 22 -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s Predict_RLRs.phy -n RLRs >RLRs.process 2>&1 &
```

## NLRs
```bash
# Run in SNORLAX
# (base) kang1234@celia-PowerEdge-T640 Mon Oct 24 22:21:44 ~/genome/Gene_annotation
less ref.fasta.ano.final|grep -i 'NLR family'|perl -alne '($nm)=$_=~/Name=\"(.*)\s+\[.*\]\"/;print $nm'|sort -u
less ref.fasta.ano.final|grep -i 'NACHT'|perl -alne '($nm)=$_=~/Name=\"(.*)\s+\[.*\]\"/;print $nm'|sort -u
less ref.fasta.ano.final|grep -i 'MHC class'|perl -alne '($nm)=$_=~/Name=\"(.*)\s+\[.*\]\"/;print $nm'|sort -u
less ref.fasta.ano.final|grep -i 'Baculoviral IAP'|perl -alne '($nm)=$_=~/Name=\"(.*)\s+\[.*\]\"/;print $nm'|sort -u
less ref.fasta.ano.final|grep -i 'Nucleotide-binding oligomerization'|perl -alne '($nm)=$_=~/Name=\"(.*)\s+\[.*\]\"/;print $nm'|sort -u
less ref.fasta.ano.final|grep -i 'NLRC'|perl -alne '($nm)=$_=~/Name=\"(.*)\s+\[.*\]\"/;print $nm'|sort -u
perl Create_query_seqs.pl NLRs_names.txt > Query_NLRs.fasta
# (base) kang1234@celia-PowerEdge-T640 Tue Oct 25 15:59:00 ~/genome/gene_family/gene_family_prediction
mkdir NLRs; cd NLRs
# (base) kang1234@celia-PowerEdge-T640 Tue Oct 25 16:01:16 ~/genome/gene_family/gene_family_prediction/NLRs
mv ~/genome/Gene_annotation/Query_NLRs.fasta ./
mv ~/genome/Gene_annotation/NLRs_names.txt ./
# Kang@fishlab3 Tue Oct 25 15:58:35 /media/HDD/cleaner_fish/genome/TLRs
scp /media/HDD/cleaner_fish/genome/Glutamates/Fmdetect.pl kang1234@147.8.76.177:~/genome/gene_family/gene_family_prediction/NLRs
scp /media/HDD/cleaner_fish/genome/Glutamates/run_Fmdetect.pl kang1234@147.8.76.177:~/genome/gene_family/gene_family_prediction/NLRs
# (base) kang1234@celia-PowerEdge-T640 Tue Oct 25 16:05:51 ~/genome/gene_family/gene_family_prediction/NLRs
nohup perl run_Fmdetect.pl /home/kang1234/genome/gene_family/gene_family_prediction/NLRs/Query_NLRs.fasta 50 /home/kang1234/genome/gene_family/gene_family_prediction/NLRs/NLRs_names.txt >Detect_NLRs.process 2>&1 &
# [1] 28400
# Run in my own workstation
# Kang@fishlab3 Wed Oct 26 00:11:03 /media/HDD/cleaner_fish/genome
mkdir NLRs; cd NLRs
# Kang@fishlab3 Wed Oct 26 00:19:41 /media/HDD/cleaner_fish/genome/NLRs
scp kang1234@147.8.76.177:~/genome/gene_family/gene_family_prediction/NLRs/NLRs_names.txt ./
scp kang1234@147.8.76.177:~/genome/gene_family/gene_family_prediction/NLRs/Query_NLRs.fasta ./
nohup perl run_Fmdetect.pl /media/HDD/cleaner_fish/genome/NLRs/Query_NLRs.fasta 50 /media/HDD/cleaner_fish/genome/NLRs/NLRs_names.txt >Detect_NLRs.process 2>&1 &
# [1] 30205
cp /media/HDD/cleaner_fish/genome/Protocadherin_new/temp6.pl ./
cp /media/HDD/cleaner_fish/genome/Glutamates/temp1.pl ./
perl temp1.pl
mv Predict_Opsins.fa Predict_NLRs.fa
mv Predict_Opsins.phy Predict_NLRs.phy
scp Predict_NLRs.phy kang1234@147.8.76.177:~/genome/gene_family/Glutamates
# (base) kang1234@celia-PowerEdge-T640 Thu Oct 27 10:45:59 ~/genome/gene_family/Glutamates
nohup raxmlHPC -T 22 -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s Predict_NLRs.phy -n NLRs >NLRs.process 2>&1 &
# [1] 5374
# Kang@fishlab3 Thu Oct 27 17:33:39 /media/HDD/cleaner_fish/genome/NLRs
less Predict_NLRs.phy|perl -alne 'next if /^\d+/;if ($F[0]=~/(.*?)_(.*?)_.*/){$hash{$2}->{$1}++};END{foreach my $ge (sort keys %hash){my $info; foreach $sp (qw(Ldim Tbif Smel Tads Spul Ncel Lber Cund)){my $nb; $hash{$ge}->{$sp}?($nb=$hash{$ge}->{$sp}):($nb=0);$info.=$nb.":";};$info=~s/\:$//;print "$ge\t[$info]"}}' >NLRs_nb.txt
# Kang@fishlab3 Thu Oct 27 17:36:07 /media/HDD/cleaner_fish/genome/NLRs
cp ../Hox/temp2.pl ./
# Kang@fishlab3 Thu Oct 27 17:56:57 /media/HDD/cleaner_fish/genome/NLRs
perl temp2.pl > NLRs_nb_2.txt
# Extract NLRC3 for the phylogenetic tree construction, revise temp1.pl
perl temp1.pl
mv Predict_Opsins.fa Predict_NLRC3.fa
mv Predict_Opsins.phy Predict_NLRC3.phy
scp Predict_NLRC3.phy kang1234@147.8.76.177:~/genome/gene_family/Glutamates
# (base) kang1234@celia-PowerEdge-T640 Fri Oct 28 01:54:53 ~/genome/gene_family/Glutamates
nohup raxmlHPC -T 22 -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s Predict_NLRC3.phy -n NLRC3 >NLRC3.process 2>&1 &
# [1] 11273
```


## Taste Receptors
```bash
# kangjingliang@kangjingliangdeMacBook-Pro 四 10 27 15:18:36 ~/Documents/2022/Ldim_genome_Restart/TRs
cat *.fasta|perl -alne 'print unless /^\s*$/' >Query_TRs.fasta
scp Query_TRs.fasta Kang@147.8.76.231:/media/HDD/cleaner_fish/genome/TRs
# Kang@fishlab3 Thu Oct 27 15:19:49 /media/HDD/cleaner_fish/genome
mkdir TRs; cd TRs
# Kang@fishlab3 Thu Oct 27 15:32:42 /media/HDD/cleaner_fish/genome/TRs
blastp -outfmt 6 -query Query_TRs.fasta -out Query_TRs.bla -db ~/Desktop/Annotation_database/swiss-prot/uniprot-filtered-reviewed_yes.fasta -num_threads 30
cat Query_TRs.bla |perl -lane 'print if $F[0] ne $flag; $flag=$F[0]'|perl -aF'\t' -lne '$h=join "\t", @F;print $h' |sort -k2 >Query_TRs.bla.best
# Kang@fishlab3 Thu Oct 27 11:12:19 /media/HDD/cleaner_fish/genome/RLRs
less Query_TRs.bla.best|perl -alne '(my $ne)=$F[1]=~/sp\|.*\|(.*)_.*/;print $ne'|sort -u
# SCNNA
# SCNNB
# SCNNG
# T2R40
# TS1R1
# TS1R3
# Search the name in swiss-prot and add the gene description in TRs_names.txt
# (base) kang1234@celia-PowerEdge-T640 Thu Oct 27 15:58:14 ~/genome/Gene_annotation
perl Create_query_seqs.pl TRs_names.txt > Query_TRs.fasta
# Kang@fishlab3 Thu Oct 27 16:00:03 /media/HDD/cleaner_fish/genome/TRs
scp kang1234@147.8.76.177:~/genome/Gene_annotation/Query_TRs.fasta ./
cp ../NLRs/Fmdetect.pl ./
cp ../NLRs/run_Fmdetect.pl ./
nohup perl run_Fmdetect.pl /media/HDD/cleaner_fish/genome/TRs/Query_TRs.fasta 50 /media/HDD/cleaner_fish/genome/TRs/TRs_names.txt >Detect_TRs.process 2>&1 &
# [1] 10174
# Kang@fishlab3 Fri Oct 28 13:39:47 /media/HDD/cleaner_fish/genome/TRs
cp /media/HDD/cleaner_fish/genome/Protocadherin_new/temp6.pl ./
cp /media/HDD/cleaner_fish/genome/Glutamates/temp1.pl ./
perl temp1.pl
mv Predict_Opsins.fa Predict_TRs.fa
mv Predict_Opsins.phy Predict_TRs.phy
# Kang@fishlab3 Fri Oct 28 13:40:12 /media/HDD/cleaner_fish/genome/TRs
scp Predict_TRs.phy kang1234@147.8.76.177:~/genome/gene_family/Glutamates
nohup raxmlHPC -T 22 -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s Predict_TRs.phy -n TRs >TRs.process 2>&1 &
# Kang@fishlab3 Mon Oct 31 10:04:32 /media/HDD/cleaner_fish/genome/TRs
perl ../Opsin_new/label_color.pl Predict_TRs.fa|less
cp ../NLRs/temp2.pl ./
less Predict_TRs.phy|perl -alne 'next if /^\d+/;if ($F[0]=~/(.*?)_(.*?)_.*/){$hash{$2}->{$1}++};END{foreach my $ge (sort keys %hash){my $info; foreach $sp (qw(Ldim Tbif Smel Tads Spul Ncel Lber Cund)){my $nb; $hash{$ge}->{$sp}?($nb=$hash{$ge}->{$sp}):($nb=0);$info.=$nb.":";};$info=~s/\:$//;print "$ge\t[$info]"}}' >TRs_nb.txt
perl temp2.pl TRs_nb.txt > TRs_nb_2.txt
```

## Immunoglobulin
```bash
# (base) kang1234@celia-PowerEdge-T640 Fri Oct 28 00:25:17 ~/genome/Gene_annotation
less ref.fasta.ano.final|grep -i 'Immunoglobulin heavy variable'|perl -alne '($nm)=$_=~/Name=\"(.*)\s+\[.*\]\"/;print $nm'|sort -u >ImmuHeavy_names.txt
perl Create_query_seqs.pl ImmuHeavy_names.txt > Query_ImmuHeavy.fasta
# Kang@fishlab3 Fri Oct 28 00:28:36 /media/HDD/cleaner_fish/genome
mkdir ImmuHeavy; cd ImmuHeavy
# Kang@fishlab3 Fri Oct 28 00:29:24 /media/HDD/cleaner_fish/genome/ImmuHeavy
scp kang1234@147.8.76.177:~/genome/Gene_annotation/Query_ImmuHeavy.fasta ./
scp kang1234@147.8.76.177:~/genome/Gene_annotation/ImmuHeavy_names.txt ./
cp ../NLRs/Fmdetect.pl ./
cp ../NLRs/run_Fmdetect.pl ./
nohup perl run_Fmdetect.pl /media/HDD/cleaner_fish/genome/ImmuHeavy/Query_ImmuHeavy.fasta 50 /media/HDD/cleaner_fish/genome/ImmuHeavy/ImmuHeavy_names.txt >Detect_ImmuHeavy.process 2>&1 &
# [1] 822
cp /media/HDD/cleaner_fish/genome/Protocadherin_new/temp6.pl ./
cp /media/HDD/cleaner_fish/genome/Glutamates/temp1.pl ./
perl temp1.pl
mv Predict_Opsins.fa Predict_ImmuHeavy.fa
mv Predict_Opsins.phy Predict_ImmuHeavy.phy
scp Predict_ImmuHeavy.phy kang1234@147.8.76.177:~/genome/gene_family/Glutamates
# (base) kang1234@celia-PowerEdge-T640 Thu Oct 27 10:45:59 ~/genome/gene_family/Glutamates
nohup raxmlHPC -T 22 -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s Predict_ImmuHeavy.phy -n ImmuHeavy > ImmuHeavy.process 2>&1 &
# [1] 14207
# Kang@fishlab3 Mon Oct 31 10:38:56 /media/HDD/cleaner_fish/genome/ImmuHeavy
cp ../Opsin_new/node_symbol.pl ./
cp ../TRs/temp2.pl ./
less Predict_ImmuHeavy.phy|perl -alne 'next if /^\d+/;if ($F[0]=~/(.*?)_(.*?)_.*/){$hash{$2}->{$1}++};END{foreach my $ge (sort keys %hash){my $info; foreach $sp (qw(Ldim Tbif Smel Tads Spul Ncel Lber Cund)){my $nb; $hash{$ge}->{$sp}?($nb=$hash{$ge}->{$sp}):($nb=0);$info.=$nb.":";};$info=~s/\:$//;print "$ge\t[$info]"}}' >ImmuHeavy_nb.txt
perl temp2.pl ImmuHeavy_nb.txt > ImmuHeavy_nb_2.txt
```

## Immediate early response genes
```bash
# (base) kang1234@celia-PowerEdge-T640 Fri Oct 28 10:13:54 ~/genome/Gene_annotation
less ref.fasta.ano.final|grep -i 'Immediate early response'|perl -alne '($nm)=$_=~/Name=\"(.*)\s+\[.*\]\"/;print $nm'|sort -u >IERs_names.txt
perl Create_query_seqs.pl IERs_names.txt > Query_IERs.fasta
# Kang@fishlab3 Fri Oct 28 10:17:36 /media/HDD/cleaner_fish/genome
mkdir IERs; cd IERs
# Kang@fishlab3 Fri Oct 28 10:17:53 /media/HDD/cleaner_fish/genome/IERs
scp kang1234@147.8.76.177:~/genome/Gene_annotation/Query_IERs.fasta ./
scp kang1234@147.8.76.177:~/genome/Gene_annotation/IERs_names.txt ./
cp ../NLRs/Fmdetect.pl ./
cp ../NLRs/run_Fmdetect.pl ./
nohup perl run_Fmdetect.pl /media/HDD/cleaner_fish/genome/IERs/Query_IERs.fasta 50 /media/HDD/cleaner_fish/genome/IERs/IERs_names.txt >Detect_IERs.process 2>&1 &
# [1] 7561
cp /media/HDD/cleaner_fish/genome/Protocadherin_new/temp6.pl ./
cp /media/HDD/cleaner_fish/genome/TRs/temp1.pl ./
perl temp1.pl
mv Predict_Opsins.fa Predict_IERs.fa
mv Predict_Opsins.phy Predict_IERs.phy
scp Predict_IERs.phy kang1234@147.8.76.177:~/genome/gene_family/Glutamates
# (base) kang1234@celia-PowerEdge-T640 Thu Oct 27 10:45:59 ~/genome/gene_family/Glutamates
nohup raxmlHPC -T 22 -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s Predict_IERs.phy -n IERs > IERs.process 2>&1 &
# [1] 14207
# Kang@fishlab3 Mon Oct 31 11:26:38 /media/HDD/cleaner_fish/genome/IERs
perl ../Opsin_new/label_color.pl Predict_IERs.fa|less
cp ../TRs/temp2.pl ./
less Predict_IERs.phy|perl -alne 'next if /^\d+/;if ($F[0]=~/(.*?)_(.*?)_.*/){$hash{$2}->{$1}++};END{foreach my $ge (sort keys %hash){my $info; foreach $sp (qw(Ldim Tbif Smel Tads Spul Ncel Lber Cund)){my $nb; $hash{$ge}->{$sp}?($nb=$hash{$ge}->{$sp}):($nb=0);$info.=$nb.":";};$info=~s/\:$//;print "$ge\t[$info]"}}' >IERs_nb.txt
perl temp2.pl IERs_nb.txt > IERs_nb_2.txt
```
## Change the color of node symbol: Lber to "#999999"
### 1. OR
```bash
# kangjingliang@kangjingliangdeMacBook-Pro 五 11 04 11:07:44 ~/Desktop
scp kang1234@147.8.76.177:~/genome/ORs_phy/RAxML_bipartitions.*-C ./

# Delta
# Kang@fishlab3 Fri Nov 04 11:05:12 /media/HDD/cleaner_fish/genome/OR_detection/ORs_class
perl label_color.pl Delta-C.fasta|less
perl node_symbol.pl Delta-C.fasta|less
```
### 2. BMPs
```bash
# Kang@fishlab3 Fri Nov 04 15:47:12 /media/HDD/cleaner_fish/genome/BMPs
less Predict_BMPs.phy|perl -alne 'next if /^\d+/;if ($F[0]=~/(.*?)_(.*?)_.*/){$hash{$2}->{$1}++};END{foreach my $ge (sort keys %hash){my $info; foreach $sp (qw(Ldim Tbif Smel Tads Spul Lber Ncel Cund)){my $nb; $hash{$ge}->{$sp}?($nb=$hash{$ge}->{$sp}):($nb=0);$info.=$nb.":";};$info=~s/\:$//;print "$ge\t[$info]"}}' >BMPs_nb.txt
perl temp2.pl BMPs_nb.txt > BMPs_nb_2.txt
```
### 3. Crystallins
```bash
# Kang@fishlab3 Fri Nov 04 15:59:25 /media/HDD/cleaner_fish/genome/Crystallins
cp ../BMPs/temp2.pl ./
less Predict_Crystallins.phy|perl -alne 'next if /^\d+/;if ($F[0]=~/(.*?)_(.*?)_.*/){$hash{$2}->{$1}++};END{foreach my $ge (sort keys %hash){my $info; foreach $sp (qw(Ldim Tbif Smel Tads Spul Lber Ncel Cund)){my $nb; $hash{$ge}->{$sp}?($nb=$hash{$ge}->{$sp}):($nb=0);$info.=$nb.":";};$info=~s/\:$//;print "$ge\t[$info]"}}' >Crystallins_nb.txt
perl temp2.pl Crystallins_nb.txt > Crystallins_nb_2.txt
```

### 3. Opsin
```bash
# Kang@fishlab3 Sun Nov 06 15:49:14 /media/HDD/cleaner_fish/genome/Opsin_new
cp ../BMPs/temp2.pl ./ temp4.pl
less Predict_Opsins.phy|perl -alne 'next if /^\d+/;if ($F[0]=~/(.*?)_(.*?)_.*/){$hash{$2}->{$1}++};END{foreach my $ge (sort keys %hash){my $info; foreach $sp (qw(Ldim Tbif Smel Tads Spul Lber Ncel Cund)){my $nb; $hash{$ge}->{$sp}?($nb=$hash{$ge}->{$sp}):($nb=0);$info.=$nb.":";};$info=~s/\:$//;print "$ge\t[$info]"}}' >Opsins_nb.txt
perl temp4.pl Opsins_nb.txt > Opsins_nb_2.txt
```
