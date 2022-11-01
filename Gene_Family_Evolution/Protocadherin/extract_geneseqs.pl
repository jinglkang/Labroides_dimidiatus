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
        s/x//g;
        $pepseq{$pepnm}.=$_;
    }
}

my $physeq="Predict_Opsins.fa";
open PHYSEQ, ">$physeq" or die "can not open $physeq\n";
foreach my $pep (sort keys %pepseq) {
    my $seq=$pepseq{$pep};
    my ($nm, $str);
    if ($pep=~/.*\_(.*)\_\d+(\D+)/) {
        $nm=$1; $str=$2;
    }
    unless ( $str eq "F") {
        my $nm1=uc($nm);
        my $nm2=uc($ARGV[0]);
        print PHYSEQ ">$pep\n$seq\n" if $nm1 eq $nm2;
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
