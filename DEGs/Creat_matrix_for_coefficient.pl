#!/usr/bin/perl -w
use strict;
use warnings;

my %bams;
my @bams=qw(LD5FB.sorted.bam LD6FB.sorted.bam LD15FB.sorted.bam LD16FB.sorted.bam LD25FB.sorted.bam LD26FB.sorted.bam LD3FB.sorted.bam LD4FB.sorted.bam LD13FB.sorted.bam LD14FB.sorted.bam LD23FB.sorted.bam LD24FB.sorted.bam LD5HB.sorted.bam LD6HB.sorted.bam LD15HB.sorted.bam LD16HB.sorted.bam LD25HB.sorted.bam LD26HB.sorted.bam LD3HB.sorted.bam LD4HB.sorted.bam LD13HB.sorted.bam LD14HB.sorted.bam LD23HB.sorted.bam LD24HB.sorted.bam LD5MB.sorted.bam LD6MB.sorted.bam LD15MB.sorted.bam LD16MB.sorted.bam LD25MB.sorted.bam LD26MB.sorted.bam LD3MB.sorted.bam LD4MB.sorted.bam LD13MB.sorted.bam LD14MB.sorted.bam LD23MB.sorted.bam LD24MB.sorted.bam);
my $head="Geneid\tChr\tStart\tEnd\tStrand\tLength\t";
foreach my $bam (@bams) {
	$head.=$bam."\t";
	$bams{$bam}++;
}
$head=~s/\s+$//;

my %pos;
my $read="all_inds_read_nb.txt";
open READ, $read or die "can not open $read\n";
while (<READ>) {
	chomp;
	my @a=split;
	if (/^#/) {
		print "$_\n";
	} elsif (/^Geneid/) {
		print "$head\n";
		my @heads=split;
		for (my $i = 0; $i < @heads; $i++) {
			if ($bams{$heads[$i]}) {
				$pos{$i}++;
			}
		}
	} else {
		my @a=split;
		my $info;
		for (my $i = 0; $i < @a; $i++) {
			if ($i<=5) {
				$info.=$a[$i]."\t";
			} elsif ($pos{$i}) {
				$info.=$a[$i]."\t";
			}
		}
		$info=~s/\s+$//;
		print "$info\n";
	}
}
