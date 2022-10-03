#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;

my @anos="/media/HDD/cleaner_fish/Genome_analysis_restart/ano_files/*.final";
my $ano ="all_ano.txt";
system("cat @anos > $ano");
my %anno;
open ANO, $ano or die "can not open $ano\n";
while (<ANO>) {
        chomp;
        my @a=split /\t/;
        if (/scov=\"(.*?)\";/) {
                $anno{$a[0]}=$_ if $1>=70;
        }
}

my @fas=</media/HDD/cleaner_fish/genome/gene_family_2/longest_pep/*.fasta>;
foreach my $fas (@fas) {
        my ($fasta)=basename($fas);
        open FASTA, $fas or die "can not open $fas\n";
        my %seq; my $gene; my @genes;
        while (<FASTA>) {
                chomp;
                if (/>/) {
                        s/>//;
                        $gene=$_;
                        push @genes, $gene;
                } else {
                        $seq{$gene}.=$_;
                }
        }
        close FASTA;

        open FILE, ">$fasta" or die "can not create $fasta\n";
        foreach my $gene (@genes) {
                my $len=length($seq{$gene});
                if ($len>=50 && $seq{$gene} !~ /\*\S/ && $anno{$gene}) {
                        print FILE ">$gene\n$seq{$gene}\n";
                }
        }
        close FILE;
}
