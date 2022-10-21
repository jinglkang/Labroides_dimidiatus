# Opsins prediction
## Create query sequences
```bash
# 1. extract the gene names containing "opsin" from six reference fish species and search the names one by one to confirm the gene name is "opsin"
# save the target name in "Opsin_names.txt"
# (base) kang1234@celia-PowerEdge-T640 Thu Oct 20 14:00:54 ~/genome/Gene_annotation
less ref.fasta.ano.final|grep -i 'opsin'|perl -alne '($nm)=$_=~/Name=\"(.*)\s+\[.*\]\"/;print $nm'|sort -u
```
Extract the sequences from six reference fish species    
```Create_query_seqs.pl
#!/usr/bin/perl
use strict;
use warnings;

my $genenm=$ARGV[0];
my %genm;
open GENENM, $genenm or die "can not open $genenm\n";
while (<GENENM>) {
        chomp;
        $genm{$_}++;
}

my %tapep;
my $refano="ref.fasta.ano.final";
open REFANO, $refano or die "can not open $refano\n";
while (<REFANO>) {
        chomp;
        my $nm;
        if (/Name=\"(.*)\s+\[.*\]\"/) {
                $nm=$1;
                my @a=split /\t/;
                $tapep{$a[0]}++ if $genm{$nm};
        }
}

my %pepseq;
my $refseq="ref.fasta";
my $pepnm;
open REFSEQ, $refseq or die "can not open $refseq\n";
while (<REFSEQ>) {
        chomp;
        if (/^>/) {
                s/\>//;
                $pepnm=$_;
        } else {
                $pepseq{$pepnm}.=$_;
        }
}

foreach my $pep (sort keys %tapep) {
        my $seq=$pepseq{$pep};
        print ">$pep\n$seq\n";
}
```

```bash
# (base) kang1234@celia-PowerEdge-T640 Thu Oct 20 14:08:56 ~/genome/Gene_annotation
perl Create_query_seqs.pl Opsin_names.txt >Query_opsins.fasta
```
## Start prediction
```bash
# Kang@fishlab3 Thu Oct 20 14:16:12 /media/HDD/cleaner_fish/genome
mkdir Opsin_new; cd Opsin_new
# Kang@fishlab3 Thu Oct 20 14:17:28 /media/HDD/cleaner_fish/genome/Opsin_new
scp kang1234@147.8.76.177:~/genome/Gene_annotation/Query_opsins.fasta ./
scp kang1234@147.8.76.177:~/genome/Gene_annotation/Opsin_names.txt ./
# revise the filtering, the mapped name should be included by "Opsin_names.txt"
nohup perl run_Fmdetect.pl /media/HDD/cleaner_fish/genome/Opsin_new/Query_opsins.fasta 50 /media/HDD/cleaner_fish/genome/Opsin_new/Opsin_names.txt >Detect_opsin.process 2>&1 &
# [1] 13971
```
## Phylogenetic tree construction
```temp1.pl
#!/usr/bin/perl
use strict;
use warnings;

my $allseq="Predict_Opsins.fasta";
system("cat */filtering/3_phy.bla.fa > $allseq");
my %pepseq;
my $pepnm;
open ALLSEQ, $allseq or die "can not open $allseq\n";
while (<ALLSEQ>) {
        chomp;
        if (/^>/) {
                s/>//;
                my @a=split;
                $pepnm=$a[0];
        } else {
                $pepseq{$pepnm}.=$_;
        }
}

my $physeq="Predict_Opsins.fa";
open PHYSEQ, ">$physeq" or die "can not open $physeq\n";
foreach my $pep (sort keys %pepseq) {
        my $seq=$pepseq{$pep};
        unless ($pep=~/F$/) {
                print PHYSEQ ">$pep\n$seq\n";
        }
}

my $align="Phy_Opsins_align.fa";
my $trim ="Phy_Opsins_align_trim.fa";
my $conc ="Phy_Opsins_align_trim_conc.fa";
my $phy  ="Predict_Opsins.phy";

system("muscle -in $physeq -out $align");
system("trimal -in $align -out $trim -gt 0.8 -st 0.001 -cons 60");
system("perl temp6.pl $trim > $conc");
system("fasta2phy.pl $conc > $phy");
system("rm $allseq $align $trim $conc");
```

```bash
# Select the predicted Opsins except the predicted pseudogene gene, and then align, trim, concatenate, phylip file for RaxMl
# Kang@fishlab3 Fri Oct 21 11:47:57 /media/HDD/cleaner_fish/genome/Opsin_new
cp /media/HDD/cleaner_fish/genome/Protocadherin_new/temp6.pl ./
# Kang@fishlab3 Fri Oct 21 12:05:26 /media/HDD/cleaner_fish/genome/Opsin_new
perl temp1.pl
# (base) kang1234@celia-PowerEdge-T640 Fri Oct 21 12:08:33 ~/genome/gene_family
mkdir Opsins; cd Opsins/
# Kang@fishlab3 Fri Oct 21 12:08:22 /media/HDD/cleaner_fish/genome/Opsin_new
scp Predict_Opsins.phy kang1234@147.8.76.177:~/genome/gene_family/Opsins
# (base) kang1234@celia-PowerEdge-T640 Fri Oct 21 12:10:07 ~/genome/gene_family/Opsins
nohup raxmlHPC -T 24 -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s Predict_Opsins.phy -n Opsins >Opsins.process 2>&1 &
# [1] 5718
```
## Extract Pinopsins sequences as outgroup
```temp3.pl
#!/usr/bin/perl
use strict;
use warnings;

my $allseq="Predict_Opsins.fasta";
my $outgro="Pinopsins.fa";
my $physeq="Predict_Opsins.fa";

system("cat $outgro $physeq > $allseq");

my $align="Phy_Opsins_align.fa";
my $trim ="Phy_Opsins_align_trim.fa";
my $conc ="Phy_Opsins_align_trim_conc.fa";
my $phy  ="Predict_Opsins.phy";

system("muscle -in $allseq -out $align");
system("trimal -in $align -out $trim -gt 0.8 -st 0.001 -cons 60");
system("perl temp6.pl $trim > $conc");
system("fasta2phy.pl $conc > $phy");
system("rm $allseq $align $trim $conc");
```

```bash
# "1.txt": Pinopsin
# (base) kang1234@celia-PowerEdge-T640 Fri Oct 21 12:57:26 ~/genome/Gene_annotation
perl Create_query_seqs.pl 1.txt >Pinopsins.fa
# Kang@fishlab3 Fri Oct 21 13:00:56 /media/HDD/cleaner_fish/genome/Opsin_new
scp kang1234@147.8.76.177:~/genome/Gene_annotation/Pinopsins.fa ./
# Kang@fishlab3 Fri Oct 21 13:10:57 /media/HDD/cleaner_fish/genome/Opsin_new
perl temp3.pl
scp Predict_Opsins.phy kang1234@147.8.76.177:~/genome/gene_family/Opsins
nohup raxmlHPC -T 24 -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s Predict_Opsins.phy -o Fugu_ENSTRUG00000004747,Fugu_ENSTRUG00000013247,Spottedgar_ENSLOCG00000001991,Spottedgar_ENSLOCG00000004577 -n Opsins >Opsins.process 2>&1 &
# [1] 7813
```
