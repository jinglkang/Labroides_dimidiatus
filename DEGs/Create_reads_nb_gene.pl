#!/usr/bin/perl
use strict;
use warnings;
use File::Basename; 
use Cwd qw(cwd);
use Getopt::Long 'HelpMessage';

my $usage=<<_EOH_;;
----------------------------------------------------------------------------------------
This script is used to create the dataframe for reads nb plot of a specific gene
Usage:
perl Create_reads_nb_gene.pl --matr DEGs_social_reads_nb_tpm.txt --gene Ldim_g17130

                                        Kang 2022-11-17
----------------------------------------------------------------------------------------
_EOH_
;

GetOptions(
    'matr:s', \my $matr,
    'gene:s', \my $gene,
    );

unless ($matr && $gene) {
    die $usage;
}

my %inds=(
    'LD5' =>'solo',
    'LD6' =>'solo',
    'LD15'=>'solo',
    'LD16'=>'solo',
    'LD25'=>'solo',
    'LD26'=>'solo',
    'LD3' =>'inte',
    'LD4' =>'inte',
    'LD13'=>'inte',
    'LD14'=>'inte',
    'LD23'=>'inte',
    'LD24'=>'inte',
    );

my %hash; my @heads;
open MATR, $matr or die "can not open $matr\n";
while (<MATR>) {
    chomp;
    my @a=split /\t/;
    if (/^\s+/) {
        @heads=@a;
    } else {
        for (my $i = 1; $i < @a; $i++) {
            $hash{$heads[$i]}->{$a[0]}=$a[$i];
        }
    }
}

print "Gene\tInd\tTissue\tType\tTPM\n";
foreach my $ind (@heads) {
    next unless $ind=~/LD/;
    my ($nm, $tiss)=$ind=~/(LD\d+)(\D{2})/;
    my $type;
    ($inds{$nm} eq "solo")?($type="1No-interaction"):($type="2Interaction");
    my $tpm=$hash{$ind}->{$gene};
    print "$gene\t$ind\t$tiss\t$type\t$tpm\n";
}
