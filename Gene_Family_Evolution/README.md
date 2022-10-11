# Detection of the gene family size variation
## make sure the gene can be mapped to swiss-prot and refseq, and keep the gene if the subject gene with coverage (mapped region) more than 70%
```bash
# (base) kang1234@celia-PowerEdge-T640 Sun Oct 02 23:28:23 ~/genome/Gene_annotation/Semicossyphus_pulcher
less Semicossyphus_pulcher.pep.all.1.conca.long.fasta.ano.final|perl -alne 'if (/scov=\"(.*?)\";/){print if $1>=70}'|wc -l # 19,075

# (base) kang1234@celia-PowerEdge-T640 Sun Oct 02 23:34:10 ~/genome/Gene_annotation/Symphodus_melops
less Symphodus_melops.pep.all.1.conca.long.fasta.ano.final|perl -alne 'if (/scov=\"(.*?)\";/){print if $1>=70}'|wc -l # 17,619

# (base) kang1234@celia-PowerEdge-T640 Sun Oct 02 23:35:23 ~/genome/Gene_annotation/Tautogolabrus_adspersus
less Tautogolabrus_adspersus.pep.all.1.conca.long.fasta.ano.final|perl -alne 'if (/scov=\"(.*?)\";/){print if $1>=70}'|wc -l # 19,113

# (base) kang1234@celia-PowerEdge-T640 Sun Oct 02 23:36:42 ~/genome/Gene_annotation/Thalassoma_bifasciatum
less Thalassoma_bifasciatum.pep.all.1.conca.long.fasta.ano.final|perl -alne 'if (/scov=\"(.*?)\";/){print if $1>=70}'|wc -l # 21,109

# (base) kang1234@celia-PowerEdge-T640 Sun Oct 02 23:38:38 ~/genome/Gene_annotation/Notolabrus_celidotus
less Notolabrus_celidotus.pep.all.1.conca.long.fasta.ano.final|perl -alne 'if (/scov=\"(.*?)\";/){print if $1>=70}'|wc -l # 18,010

# (base) kang1234@celia-PowerEdge-T640 Sun Oct 02 23:39:36 ~/genome/Gene_annotation/Labrus_bergylta
less Labrus_bergylta.pep.all.1.conca.long.fasta.ano.final|perl -alne 'if (/scov=\"(.*?)\";/){print if $1>=70}'|wc -l # 19,313

# Ldim
# (base) kang1234@celia-PowerEdge-T640 Sun Oct 02 23:40:54 ~/genome/Gene_annotation/combined
less braker2+3_combined_renamed.aa.long.anno.final.txt|perl -alne 'if (/scov=\"(.*?)\";/){print if $1>=70}'|wc -l # 17,437

# (base) kang1234@celia-PowerEdge-T640 Sun Oct 02 23:42:26 ~/genome/Gene_annotation/Cheilinus_undulatus
less Cheilinus_undulatus.pep.all.1.conca.long.fasta.ano.final|perl -alne 'if (/scov=\"(.*?)\";/){print if $1>=70}'|wc -l # 18,806

# Annotate the pep seqs of other 6 ref species
# (base) kang1234@celia-PowerEdge-T640 Sun Oct 02 23:55:41 ~/genome/Gene_annotation
nohup annotate --fasta ref.fasta > ano_ref.process 2>&1 &
# [1] 20138
# revise few for the name
nohup perl annotate.pl --fasta ref.fasta > ano_ref.process 2>&1 &
# [1] 19568
# (base) kang1234@celia-PowerEdge-T640 Mon Oct 03 10:06:26 ~/genome/Gene_annotation
less ref.fasta.ano.final|perl -alne 'if (/scov=\"(.*?)\";/){print if $1>=70}'|perl -alne '($spe)=$F[0]=~/(.*?)\_/;$hash{$spe}++;END{foreach my $key (sort keys %hash){print "$key\t$hash{$key}"}}'
# Fugu	16395
# Medaka	17217
# Platyfish	17842
# Spottedgar	14683
# Stickleback	15678
# Zebrafish	22897
```
### Extract the qualified pep sequences per species
```bash
# Kang@fishlab3 Mon Oct 03 10:17:10 /media/HDD/cleaner_fish
mkdir Genome_analysis_restart; cd Genome_analysis_restart
# Kang@fishlab3 Mon Oct 03 10:20:20 /media/HDD/cleaner_fish/Genome_analysis_restart
mkdir ano_files; cd ano_files
scp kang1234@147.8.76.177:~/genome/Gene_annotation/Semicossyphus_pulcher/Semicossyphus_pulcher.pep.all.1.conca.long.fasta.ano.final ./
scp kang1234@147.8.76.177:~/genome/Gene_annotation/Symphodus_melops/Symphodus_melops.pep.all.1.conca.long.fasta.ano.final ./
scp kang1234@147.8.76.177:~/genome/Gene_annotation/Tautogolabrus_adspersus/Tautogolabrus_adspersus.pep.all.1.conca.long.fasta.ano.final ./
scp kang1234@147.8.76.177:~/genome/Gene_annotation/Thalassoma_bifasciatum/Thalassoma_bifasciatum.pep.all.1.conca.long.fasta.ano.final ./
scp kang1234@147.8.76.177:~/genome/Gene_annotation/Notolabrus_celidotus/Notolabrus_celidotus.pep.all.1.conca.long.fasta.ano.final ./
scp kang1234@147.8.76.177:~/genome/Gene_annotation/Labrus_bergylta/Labrus_bergylta.pep.all.1.conca.long.fasta.ano.final ./
scp kang1234@147.8.76.177:~/genome/Gene_annotation/combined/braker2+3_combined_renamed.aa.long.anno.final.txt ./
scp kang1234@147.8.76.177:~/genome/Gene_annotation/Cheilinus_undulatus/Cheilinus_undulatus.pep.all.1.conca.long.fasta.ano.final ./
scp kang1234@147.8.76.177:~/genome/Gene_annotation/ref.fasta.ano.final ./

# Kang@fishlab3 Mon Oct 03 10:38:04 /media/HDD/cleaner_fish/Genome_analysis_restart
mkdir gene_family_seq_input; cd gene_family_seq_input

# the seq length should be more than 50 && can be mapped to swiss-prot
# Kang@fishlab3 Mon Oct 03 11:02:00 /media/HDD/cleaner_fish/Genome_analysis_restart/gene_family_seq_input
perl extract_seq.pl
# Kang@fishlab3 Mon Oct 03 11:00:39 /media/HDD/cleaner_fish/Genome_analysis_restart/gene_family_seq_input
for fa in *.fasta;do echo -n "$fa   ";grep '>' ${fa}|wc -l;done
# Cheilinus_undulatus.fasta   16797
# Fugu.fasta   15909
# Labroides_dimidiatus.fasta   16621
# Labrus_bergylta.fasta   17223
# Medaka.fasta   16600
# Notolabrus_celidotus.fasta   16118
# Patyfish.fasta   17175
# Semicossyphus_pulcher.fasta   17344
# Spottedgar.fasta   14365
# Stickleback.fasta   15368
# Symphodus_melops.fasta   16313
# Tautogolabrus_adspersus.fasta   17150
# Thalassoma_bifasciatum.fasta   19294
# Zebrafish.fasta   21570
```

### orthofinder detect the orthogroups
```bash
# Kang@fishlab3 Mon Oct 03 11:08:36 /media/HDD/cleaner_fish/Genome_analysis_restart
nohup orthofinder -f gene_family_seq_input -a 32 >Orthofinder.process 2>&1 &
# [1] 6631

# the gene family should exist in at least three ref species and then annotate the gene family with species that have the most gene number in this gene family
# Kang@fishlab3 Mon Oct 03 14:07:05 /media/HDD/cleaner_fish/Genome_analysis_restart/ano_files
# cat *.final|perl -alne 'if (/scov=\"(.*?)\";/){print if $1>=70}'>all_qualified.anno.txt
# Kang@fishlab3 Mon Oct 03 14:25:00 /media/HDD/cleaner_fish/Genome_analysis_restart/gene_family_seq_input/OrthoFinder/Results_Oct03/Orthogroups
# perl filter_ano_gene_family.pl # output: fm_nb_ano.txt; fm_nb.txt; 13733 gene families
# Kang@fishlab3 Mon Oct 03 16:05:14 /media/HDD/cleaner_fish/Genome_analysis_restart/gene_family_seq_input/OrthoFinder/Results_Oct03/Orthogroups
# python clade_and_size_filter.py -i fm_nb.txt -o fm_nb_filtered.txt -s # fm_nb_filtered.txt: 13733 gene families
# scp fm_nb_filtered.txt kang1234@147.8.76.177:~/genome/gene_family

# the orthogroups were divided into gene family by gene name to make sure only one gene name in each gene family
# Kang@fishlab3 Tue Oct 11 10:40:57 /media/HDD/cleaner_fish/Genome_analysis_restart/ano_files
cat *.final|perl -alne 'if (/swiss-prot/ && /scov=\"(.*?)\";.*qlenth=\"(.*?)\"/){print if $1>=70}'>all_qualified.anno.txt
# Kang@fishlab3 Tue Oct 11 11:40:52 /media/HDD/cleaner_fish/Genome_analysis_restart/gene_family_seq_input/OrthoFinder/Results_Oct11/Orthogroups
mv ../../../OrthoFinder_old/Results_Oct03/Orthogroups/temp1.pl filter_ano_gene_family.pl
perl filter_ano_gene_family.pl
# in each gene family, at least 7 species with gene number > 0
less fm_nb.txt|perl -alne 'if (/^Desc/){print}else{my $k;for (my $i = 2; $i < @F; $i++){$k++ if $F[$i]>0};print if $k>=7}' >fm_nb_conserve.txt
python clade_and_size_filter.py -i fm_nb_conserve.txt -o fm_nb_conserve_filtered.txt -s # 14390 gene families
less fm_nb_conserve_filtered.txt |perl -alne 's/\_//g;print' >fm_nb_conserve_filtered.txt.1
mv fm_nb_conserve_filtered.txt.1 fm_nb_conserve_filtered.txt
scp fm_nb_conserve_filtered.txt kang1234@147.8.76.177:~/genome/gene_family
scp fm_nb_ano.txt kang1234@147.8.76.177:~/genome/gene_family/Restart
```
### Run cafe
```Restart.cafe
#! cafe
load -i fm_nb_conserve_filtered.txt -t 20 -l Restart/log_run1.txt -p 0.01 -r 10000
tree ((((((((((Symphodusmelops:12.207213,Tautogolabrusadspersus:12.207213):4.220335,Labrusbergylta:16.427548):16.178243,Cheilinusundulatus:32.605791):3.790756,((Labroidesdimidiatus:14.975052,Thalassomabifasciatum:14.975052):11.654002,Notolabruscelidotus:26.629053):9.767494):3.822500,Semicossyphuspulcher:40.219047):26.542146,Stickleback:66.761193):7.420018,Fugu:74.181211):5.877738,(Platyfish:68.523066,Medaka:68.523066):11.535883):76.471051,Zebrafish:156.530000):94.757281,Spottedgar:251.287281)
lambda -s -t ((((((((((1,1)1,1)1,1)1,((1,1)1,1)1)1,1)1,1)1,1)1,(1,1)1)1,1)1,1)
report Restart/report_run1
```

```bash
mkdir Restart/
cafe ./Restart.cafe
python2 report_analysis.py -i Restart/report_run1.cafe -o Restart/summary_run1
# Plot the result: Rapid, Expansions, Contractions
python2 draw_tree.py -i Restart/summary_run1_node.txt -t '((((((((((Symphodusmelops:12.2072,Tautogolabrusadspersus:12.2072):4.22034,Labrusbergylta:16.4275):16.1782,Cheilinusundulatus:32.6058):3.79076,((Labroidesdimidiatus:14.9751,Thalassomabifasciatum:14.9751):11.654,Notolabruscelidotus:26.6291):9.76749):3.8225,Semicossyphuspulcher:40.219):26.5421,Stickleback:66.7612):7.42002,Fugu:74.1812):5.87774,(Platyfish:68.5231,Medaka:68.5231):11.5359):76.4711,Zebrafish:156.53):94.7573,Spottedgar:251.287)' -d '((((((((((Symphodusmelops<0>,Tautogolabrusadspersus<2>)<1>,Labrusbergylta<4>)<3>,Cheilinusundulatus<6>)<5>,((Labroidesdimidiatus<8>,Thalassomabifasciatum<10>)<9>,Notolabruscelidotus<12>)<11>)<7>,Semicossyphuspulcher<14>)<13>,Stickleback<16>)<15>,Fugu<18>)<17>,(Platyfish<20>,Medaka<22>)<21>)<19>,Zebrafish<24>)<23>,Spottedgar<26>)<25>' -o Restart/summary_run1_tree_rapid.png -y Rapid
python2 draw_tree.py -i Restart/summary_run1_node.txt -t '((((((((((Symphodusmelops:12.2072,Tautogolabrusadspersus:12.2072):4.22034,Labrusbergylta:16.4275):16.1782,Cheilinusundulatus:32.6058):3.79076,((Labroidesdimidiatus:14.9751,Thalassomabifasciatum:14.9751):11.654,Notolabruscelidotus:26.6291):9.76749):3.8225,Semicossyphuspulcher:40.219):26.5421,Stickleback:66.7612):7.42002,Fugu:74.1812):5.87774,(Platyfish:68.5231,Medaka:68.5231):11.5359):76.4711,Zebrafish:156.53):94.7573,Spottedgar:251.287)' -d '((((((((((Symphodusmelops<0>,Tautogolabrusadspersus<2>)<1>,Labrusbergylta<4>)<3>,Cheilinusundulatus<6>)<5>,((Labroidesdimidiatus<8>,Thalassomabifasciatum<10>)<9>,Notolabruscelidotus<12>)<11>)<7>,Semicossyphuspulcher<14>)<13>,Stickleback<16>)<15>,Fugu<18>)<17>,(Platyfish<20>,Medaka<22>)<21>)<19>,Zebrafish<24>)<23>,Spottedgar<26>)<25>' -o Restart/summary_run1_tree_Expansions.png -y Expansions
python2 draw_tree.py -i Restart/summary_run1_node.txt -t '((((((((((Symphodusmelops:12.2072,Tautogolabrusadspersus:12.2072):4.22034,Labrusbergylta:16.4275):16.1782,Cheilinusundulatus:32.6058):3.79076,((Labroidesdimidiatus:14.9751,Thalassomabifasciatum:14.9751):11.654,Notolabruscelidotus:26.6291):9.76749):3.8225,Semicossyphuspulcher:40.219):26.5421,Stickleback:66.7612):7.42002,Fugu:74.1812):5.87774,(Platyfish:68.5231,Medaka:68.5231):11.5359):76.4711,Zebrafish:156.53):94.7573,Spottedgar:251.287)' -d '((((((((((Symphodusmelops<0>,Tautogolabrusadspersus<2>)<1>,Labrusbergylta<4>)<3>,Cheilinusundulatus<6>)<5>,((Labroidesdimidiatus<8>,Thalassomabifasciatum<10>)<9>,Notolabruscelidotus<12>)<11>)<7>,Semicossyphuspulcher<14>)<13>,Stickleback<16>)<15>,Fugu<18>)<17>,(Platyfish<20>,Medaka<22>)<21>)<19>,Zebrafish<24>)<23>,Spottedgar<26>)<25>' -o Restart/summary_run1_tree_Contractions.png -y Contractions
```
# Extract the sig. expansion and contraction gene family
```bash
# (base) kang1234@celia-PowerEdge-T640 Tue Oct 11 16:40:42 ~/genome/gene_family/Restart
perl extract_sig_fm.pl
```
