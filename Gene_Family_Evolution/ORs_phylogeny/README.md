# Build the phylogenetic tree for each OR group
## If the predicted pep sequences include "!", remove it; detect the domian; keep the seq if contain "7tm"; and include the Zebrafish ORs and Non-ORs to build a phylogenetic tree to confirm the classification
```temp1.pl
#!/usr/bin/perl
use strict;
use warnings;

my %pep;
my $pepseq="filter.out.1.info.fasta1";
my $pepnm;
open PEP, $pepseq or die "can not open $pepseq\n$!\n";
while (<PEP>) {
        chomp;
        if (/>/) {
                s/>//;
                my @a=split;
                $pepnm=$a[0];
        } else {
                $pep{$pepnm}.=$_;
        }
}

# domain should contain "7tm"
my %pfam;
my @genes;
my $pfamR="filter_pfam.txt";
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

# include the zebrafish pep seq for classify
my $Rule="ruler.fa";
my $clas_nm;
my @clas_nms;
my %Class;
open RULE, $Rule or die "can not open $Rule\n$!\n";
while (<RULE>) {
        chomp;
        if (/>/) {
                s/>//;
                $clas_nm=$_;
                push @clas_nms, $clas_nm;
        } else {
                $Class{$clas_nm}.=$_;
        }
}

foreach my $gene (@clas_nms) {
        if ($gene=~/NonOR/ || $gene=~/D\.rerio.*/) {
                print ">$gene\n$Class{$gene}\n";
        }
}

foreach my $name (sort keys %pep) {
        if ($pfam{$name}) {
                print ">$name\n$pep{$name}\n";
        }
}
```

```1.sh
less filter.out.1.info.fasta|perl -alne 's/\!//g;print'>filter.out.1.info.fasta1 # remove "!" in frameshift pep sequences
~/software/PfamScan/pfam_scan.pl -fasta filter.out.1.info.fasta1 -dir /media/HDD/pfam/ >filter_pfam.txt
perl temp1.pl >OR_phy.fa
muscle -in OR_phy.fa -out OR_phy_align.fas
less OR_phy_align.fas|perl -alne 's/\:/|/g;print' >OR_phy_align.fas1
FastTreeMP OR_phy_align.fas1 >OR_phy_align.phy
```

```bash
# Kang@fishlab3 Wed Oct 12 02:22:51 /media/HDD/cleaner_fish/genome/OR_detection/Labrus_bergylta/group
cp /media/HDD/cleaner_fish/genome/OR_detection/Cheilinus_undulatus/group/temp1.pl ./
cp /media/HDD/cleaner_fish/genome/OR_detection/Cheilinus_undulatus/group/1.sh ./
sh 1.sh
```
