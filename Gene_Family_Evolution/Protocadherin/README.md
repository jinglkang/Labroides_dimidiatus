# Protocadherin gene family detection
## Download the query protein sequences
### Select the protocadherin sequences of human in Uniprot, and also select the genes in uniprot if the gene name is same with the human protocadherin gene name   
```bash
# Kang@fishlab3 Fri Oct 21 00:02:12 /media/HDD/cleaner_fish/genome/Protocadherin_new
cat ../Protocadherin_alpha/query_protein.fasta ../Protocadherin_beta/query_protein.fasta ../Protocadherin_gamma/query_protein.fasta ../Protocadherin_non/query_protein.fasta >query_protocadherin_protein.fasta
# Select the target predicted genes: e-value<=10e-20 and gene decription inculding keyword "protocadherin"
nohup perl run_Fmdetect.pl /media/HDD/cleaner_fish/genome/Protocadherin_new/query_protocadherin_protein.fasta 50 >Detect_protocadherin.process 2>&1 &
# [1] 9367
```
### Domain detection (remove "!" if the predicted gene is a pseudogene)     
```temp1.pl
#!/usr/bin/perl
use strict;
use warnings;

my $query="query_protocadherin_protein.fasta";
my %que;
my $quenm;
open QUERY, $query or die "can not open $query\n$!\n";
while (<QUERY>) {
        chomp;
        if (/^\>/) {
                s/\>//;
                my @a=split;
                $quenm=$a[0];
        } else {
                $que{$quenm}.=$_;
        }
}

foreach my $nm (sort keys %que) {
        if ($nm=~/\_HUMAN/) {
                #               print ">$nm\n$que{$nm}\n";
        }
}

my @species=qw(Semicossyphus_pulcher Cheilinus_undulatus Labrus_bergylta Tautogolabrus_adspersus Symphodus_melops Notolabrus_celidotus Thalassoma_bifasciatum Labroides_dimidiatus);
foreach my $spe (@species) {
        my $fa="$spe/filtering/3_phy.bla.fa";
        my $fa1="$spe/filtering/3_phy.bla.fa1";
        # Remove "!" in the pep sequences
        open FA, $fa or die "can not open $fa\n$!\n";
        open FA1,">$fa1" or die "can not open $fa1\n$!\n";
        while (<FA>) {
                chomp;
                if (/^\>/) {
                        print FA1 "$_\n";
                } else {
                        s/\!//g;
                        print FA1 "$_\n";
                }
        }
        # pfam detect domain
        close FA; close FA1;

        my %pepseq; my $pepnm;
        open FA1, "$fa1" or die "can not open $fa1\n$!\n";
        while (<FA1>) {
                chomp;
                if (/^\>/) {
                        s/\>//;
                        my @a=split;
                        $pepnm=$a[0];
                } else {
                        $pepseq{$pepnm}.=$_;
                }
        }

        my $pfamR="$spe/filtering/filter_pfam.txt";
        my %pfam;
        system("~/software/PfamScan/pfam_scan.pl -fasta $fa1 -dir /media/HDD/pfam/ >$pfamR");
        # confirm that the predicted peps have at least six extracellular Cadherin domains
        open PFAM, $pfamR or die "can not open $pfamR\n$!\n";
        while (<PFAM>) {
        chomp;
        if (/^#/ || /^\s*$/) {
                next;
        } else {
                my @F=split;
                $pfam{$F[0]}.=$F[6]."\t";
        }
        }

        foreach my $pep (sort keys %pfam) {
                if ($pfam{$pep}=~/Cadherin/i) {
                        (my $nb)=$pfam{$pep}=~s/Cadherin/Cadherin/g;
                        #               print ">$pep\n$pepseq{$pep}\n" if $nb>=6;
                }
        }
}
```

```bash
# Kang@fishlab3 Fri Oct 21 00:23:03 /media/HDD/cleaner_fish/genome/Protocadherin_new
nohup perl temp1.pl &
```
### Protocadherin gamma genes
```temp7.pl
#!/usr/bin/perl
use strict;
use warnings;

my @species=qw(Semicossyphus_pulcher Cheilinus_undulatus Labrus_bergylta Tautogolabrus_adspersus Symphodus_melops Notolabrus_celidotus Thalassoma_bifasciatum Labroides_dimidiatus);
foreach my $spe (@species) {
        my %blast;
        my $bla="$spe/filtering/3_phy.bla";
        open BLA, $bla or die "can not open $bla\n$!\n";
        while (<BLA>) {
                chomp;
                if (/Protocadherin gamma/) {
                        my @a=split /\t/;
                        (my $nm)=$a[1]=~/sp\|.*\|(.*)\_/;
                        $blast{$nm}++;
                }
        }

        my $fa1="$spe/filtering/3_phy.bla.fa1";
        my %pepseq; my $pepnm;
        open FA1, "$fa1" or die "can not open $fa1\n$!\n";
        while (<FA1>) {
                chomp;
                if (/^\>/) {
                        s/\>//;
                        my @a=split;
                        $pepnm=$a[0];
                } else {
                        s/x//ig;
                        $pepseq{$pepnm}.=$_;
                }
        }

        my $pfamR="$spe/filtering/filter_pfam.txt";
        my %pfam;
        # confirm that the predicted peps have at least six extracellular Cadherin domains
        open PFAM, $pfamR or die "can not open $pfamR\n$!\n";
        while (<PFAM>) {
        chomp;
        if (/^#/ || /^\s*$/) {
                next;
        } else {
                my @F=split;
                $pfam{$F[0]}.=$F[6]."\t";
        }
        }

        foreach my $pep (sort keys %pfam) {
                (my $nm)=$pep=~/.*\_(.*)\_.*/;
                if ($blast{$nm} && $pfam{$pep}=~/Cadherin/i) {
                        (my $nb)=$pfam{$pep}=~s/Cadherin/Cadherin/g;
                        print ">$pep\n$pepseq{$pep}\n" if $nb>=6 && $pep !~ /F$/;
                }
        }
}
```

```bash
# revise the keyword "/Protocadherin gamma/" && not contain the predicted pseudogene genes
perl temp7.pl >1.fa
muscle -in 1.fa -out 1_align.fa
mv 1.fa Predict_Gamma_allspe.fa
trimal -in 1_align.fa -out 1_align_trim.fa -gt 0.8 -st 0.001 -cons 60
perl temp6.pl 1_align_trim.fa >1_align_trim_conc.fa
fasta2phy.pl 1_align_trim_conc.fa >1_align_trim_conc.phy
mv 1_align_trim_conc.phy Gamma_allspe.phy
scp Gamma_allspe.phy kang1234@147.8.76.177:~/genome/gene_family/Cadherin/
# (base) kang1234@celia-PowerEdge-T640 Fri Oct 21 00:39:16 ~/genome/gene_family/Cadherin
nohup raxmlHPC -T 24 -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s Gamma_allspe.phy -n Gamma_allspe >Gamma_allspe.process 2>&1 &
# [1] 5480
```
### Protocadherin alpha genes
```bash
# revise the keyword "/Protocadherin alpha/" && not contain the predicted pseudogene genes
perl temp7.pl >1.fa
mv 1.fa Predict_Alpha_allspe.fa
muscle -in 1.fa -out 1_align.fa
trimal -in 1_align.fa -out 1_align_trim.fa -gt 0.8 -st 0.001 -cons 60
perl temp6.pl 1_align_trim.fa >1_align_trim_conc.fa
fasta2phy.pl 1_align_trim_conc.fa >1_align_trim_conc.phy
mv 1_align_trim_conc.phy Alpha_allspe.phy
scp Alpha_allspe.phy kang1234@147.8.76.177:~/genome/gene_family/Cadherin/
nohup raxmlHPC -T 24 -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s Alpha_allspe.phy -n Alpha_allspe >Alpha_allspe.process 2>&1 &
# [1] 4779
```
### iTol set the node shape and label colour
```bash
# Kang@fishlab3 Fri Oct 21 10:51:43 /media/HDD/cleaner_fish/genome/Protocadherin_new
perl temp4.pl 1.fa|less
# Kang@fishlab3 Fri Oct 21 10:51:53 /media/HDD/cleaner_fish/genome/Protocadherin_new
perl temp5.pl 1.fa|less
```
### Summary the genes number per species
```bash
# Gamma
# (base) kang1234@celia-PowerEdge-T640 Fri Oct 21 10:36:28 ~/genome/gene_family/Cadherin
less Gamma_allspe.phy|perl -alne 'next if /^\d+/;if ($F[0]=~/(.*?)_(.*?)_.*/ && $2 eq "PCDGM"){$hash{$1}++};END{foreach $sp (qw(Ldim Tbif Smel Tads Spul Ncel Lber Cund)){print "$sp\t$hash{$sp}"}}'
less Gamma_allspe.phy|perl -alne 'next if /^\d+/;unless ($F[0]=~/(.*?)_(.*?)_.*/ && $2 eq "PCDGM"){$hash{$1}++};END{foreach $sp (qw(Ldim Tbif Smel Tads Spul Ncel Lber Cund)){print "$sp\t$hash{$sp}"}}'
# Alpha
# (base) kang1234@celia-PowerEdge-T640 Fri Oct 21 10:46:58 ~/genome/gene_family/Cadherin
less Alpha_allspe.phy|perl -alne 'next if /^\d+/;if ($F[0]=~/(.*?)_(.*?)_.*/ && $2 eq "PCDC2"){$hash{$1}++};END{foreach $sp (qw(Ldim Tbif Smel Tads Spul Ncel Lber Cund)){print "$sp\t$hash{$sp}"}}'
less Alpha_allspe.phy|perl -alne 'next if /^\d+/;unless ($F[0]=~/(.*?)_(.*?)_.*/ && $2 eq "PCDC2"){$hash{$1}++};END{foreach $sp (qw(Ldim Tbif Smel Tads Spul Ncel Lber Cund)){print "$sp\t$hash{$sp}"}}'

less Alpha_allspe.phy|perl -alne 'next if /^\d+/;if ($F[0]=~/(.*?)_(.*?)_.*/){$hash{$2}->{$1}++};END{foreach my $ge (sort keys %hash){my $info; foreach $sp (qw(Ldim Tbif Smel Tads Spul Ncel Lber Cund)){my $nb; $hash{$ge}->{$sp}?($nb=$hash{$ge}->{$sp}):($nb=0);$info.=$nb.":";};$info=~s/\:$//;print "$ge\t[$info]"}}' >Alpha_nb.txt
perl temp8.pl Alpha_nb.txt > Alpha_nb_2.txt
less Alpha_nb_2.txt
# Kang@fishlab3 Tue Nov 01 09:33:29 /media/HDD/cleaner_fish/genome/Protocadherin_new
less Gamma_allspe.phy|perl -alne 'next if /^\d+/;if ($F[0]=~/(.*?)_(.*?)_.*/){$hash{$2}->{$1}++};END{foreach my $ge (sort keys %hash){my $info; foreach $sp (qw(Ldim Tbif Smel Tads Spul Ncel Lber Cund)){my $nb; $hash{$ge}->{$sp}?($nb=$hash{$ge}->{$sp}):($nb=0);$info.=$nb.":";};$info=~s/\:$//;print "$ge\t[$info]"}}' >Gamma_nb.txt
perl temp8.pl Gamma_nb.txt > Gamma_nb_2.txt

# Extract PCDA2 (Protocadherin alpha-2: PCDHA2) for phylogeny
# Kang@fishlab3 Tue Nov 01 10:51:45 /media/HDD/cleaner_fish/genome/Protocadherin_new
# mv temp9.pl extract_geneseqs.pl
perl extract_geneseqs.pl PCDA2
mv Predict_Opsins.fa Predict_PCDA2.fa
mv Predict_Opsins.phy Predict_PCDA2.phy
scp Predict_PCDA2.phy kang1234@147.8.76.177:~/genome/gene_family/Glutamates
# (base) kang1234@celia-PowerEdge-T640 Tue Nov 01 11:19:54 ~/genome/gene_family/Glutamates
nohup raxmlHPC -T 24 -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s Predict_PCDA2.phy -n PCDA2 >PCDA2.process 2>&1 &
# [1] 11596

# Extract PCDGB (Protocadherin gamma-A11: PCDHGA11) for phylogeny
# Kang@fishlab3 Tue Nov 01 11:22:39 /media/HDD/cleaner_fish/genome/Protocadherin_new
perl extract_geneseqs.pl PCDGB
mv Predict_Opsins.fa Predict_PCDGB.fa
mv Predict_Opsins.phy Predict_PCDGB.phy
scp Predict_PCDGB.phy kang1234@147.8.76.177:~/genome/gene_family/Glutamates
# (base) kang1234@celia-PowerEdge-T640 Tue Nov 01 11:19:54 ~/genome/gene_family/Glutamates
nohup raxmlHPC -T 24 -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s Predict_PCDGB.phy -n PCDGB >PCDGB.process 2>&1 &
# [1] 11878
```
