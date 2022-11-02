#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;

my $bla="allspe_blastp.result";
my %blast;
open BLA, $bla or die "can not open $bla\n";
while (<BLA>) {
        chomp;
        my @a=split /\t/;
        my $gene=$a[0];
        $gene=~s/\.?t\d+$//;
        (my $unid)=$a[1]=~/sp\|(.*)\|.*\_.*/;
        (my $genm)=$a[1]=~/sp\|.*\|(.*)\_.*/;
        my $eval=$a[10];
        my $scor=$a[11];
        $blast{$gene}={
                'UNID' => $unid, # use this id to extract the gene description
                'NAME' => $genm, # gene name
                'EVAL' => $eval, # E-value
                'SCOR' => $scor  # blast score
        };
}

my $orthnm="sub_orth_genecount.txt";
my %nm;
open ORTHNM, $orthnm or die "can not open $orthnm\n";
while (<ORTHNM>) {
        chomp;
        my @a=split /\t/;
        next if /^Orthogroup/;
#       print "$_\n" if /^Orthogroup/;
        my @b;
        for (my $i = 1; $i < @a; $i++) {
                push @b, $a[$i] if $a[$i]>=1;
        }
        $nm{$a[0]}++ if @b==14;
}

my @speall=qw(Spottedgar Zebrafish Medaka Platyfish Fugu Stickleback Spul Cund Lber Tads Smel Ncel Tbif Ldim);
my $orth="sub_orth_id.txt";
open ORTH, $orth or die "can not open $orth\n";
while (<ORTH>) {
        chomp;
        s/\,//g;
        my @a=split;
        if ($nm{$a[0]}) { # all species in this orthogroup
                my %hash; my @spes;
                for (my $i = 1; $i < @a; $i++) {
                        if ($blast{$a[$i]}->{'NAME'}) {
                                my $name=$blast{$a[$i]}->{'NAME'};
#                               print "$name\n";
                                $hash{$name}=[] unless exists $hash{$name};
                                push @{$hash{$name}}, $a[$i];
                        }
                }

                foreach my $name (keys %hash) {
                        my @genes=@{$hash{$name}};
                        my (%hash1, %hash2); my @spes;
                        foreach my $gene (@genes) {
                                (my $spe)=$gene=~/(.*)\_.*/;
                                $hash1{$spe}++;
                                push @spes, $spe if $hash1{$spe}==1;
                                $hash2{$spe}=[] unless exists $hash2{$spe};
                                push @{$hash2{$spe}}, $gene;
                        }
                        my $number=@spes;
#                       print "$name\t$number\n";
                        if (@spes==14) {
                                my $info=$a[0]."\t";
                                foreach my $spe (@speall) {
                                        my @genes=@{$hash2{$spe}};
                                        if (@genes==1) {
                                                $info.=$genes[0]."\t";
                                        } else {
                                                my %hash3;
                                                foreach my $gene (@genes) {
                                                        my $score=$blast{$gene}->{'SCOR'};
                                                        if ($hash3{$spe}) {
                                                                my $oldscore=$hash3{$spe}->{'SCORE'};
                                                                if ($oldscore < $score) {
                                                                        $hash3{$spe}={
                                                                                'SCORE' => $score,
                                                                                'GENE'  => $gene
                                                                        };
                                                                }
                                                        } else {
                                                                $hash3{$spe}={
                                                                        'SCORE' => $score,
                                                                        'GENE'  => $gene
                                                                };
                                                        }
                                                }
                                                $info.=$hash3{$spe}->{'GENE'}."\t";
                                        }
                                }
                                $info=~s/\s+$//;
                                print "$info\n";
                        }
                }
        }
}
