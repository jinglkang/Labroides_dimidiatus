#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;

my $ano ="Ldim_ano.txt";
my %ann;
open ANO, $ano or die "can not open $ano\n";
while (<ANO>) {
        chomp;
        my @a=split /\t/;
        my $gene=$a[0];
        (my $symb)=$a[2]=~/symbol=\"(.*?)\"/;
        (my $decr)=$a[3]=~/Name=\"(.*)\s+\[.*\]\"/;
        $ann{$gene}={
                'SYMB' => $symb,
                'DESR' => $decr
        };
}

my $list="Genelist_paml.txt";
open LIST, $list or die "can not open $list\n";
while (<LIST>) {
        chomp;
        my @a=split;
        my $gene=$a[-1];
        my $orth=$a[0];
        my $symb=$ann{$gene}->{'SYMB'};
        my $decr=$ann{$gene}->{'DESR'};
        print "$orth\t$symb\t$decr\n";
}
