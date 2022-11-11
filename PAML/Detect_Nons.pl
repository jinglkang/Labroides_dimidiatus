#!/usr/bin/perl
use strict;
use warnings;

# Detect_Nons.pl
my %seq; my $name;
my $orth=$ARGV[0];
my $fas="$orth/final_alignment_pep_nodes.fa";
open FAS, $fas or die "can not open $fas\n";
while (<FAS>) {
	chomp;
	if (/^>/) {
		s/\>//;
		$name=$_;
	} else {
		$seq{$name}.=$_;
	}
}

my %cleaner=(
	'Spul'=> 1,
	'Smel'=> 1,
	'Tads'=> 1,
	'Lber'=> 1,
	'Ldim'=> 1,
	'Tbif'=> 1,
	);

# compare the nonsynonymous position pep sequences one by one
my (%hash1, %hash2);
my @grp1=qw(Spul); # compare to "Node20"
my @grp2=qw(Smel Tads Lber); # compare to "Node22"
my @grp3=qw(Ldim Tbif); # compare to "Node25"

my @nocl1=qw(Cund); # compare to "Node22"
my @nocl2=qw(Ncel); # compare to "Node25"

&Detect_Nons(\@grp1, "Node20");
&Detect_Nons(\@grp2, "Node22");
&Detect_Nons(\@grp3, "Node25");
&Detect_Nons(\@nocl1, "Node22");
&Detect_Nons(\@nocl2, "Node25");

my @nodes=qw(Node20 Spul Node22 Smel Tads Lber Cund Node25 Ldim Tbif Ncel);
my @facul=qw(Smel Tads Lber Tbif Spul);

foreach my $pos (sort {$a <=> $b} keys %hash2) {
	if ($hash2{$pos} == 6) {
		&Print_Nons($pos);
	}
}

sub Detect_Nons {
	my ($grp, $ref)=@_;
	my @grp=@{$grp};
	my $refseq=$seq{$ref};
	foreach my $spe (@grp) {
		my $seq=$seq{$spe};
		my $len=length($seq);
		for (my $i = 0; $i < $len; $i++) {
			my $refpos=substr($refseq,$i,1);
			$hash1{$ref}->{$i}=$refpos;
			my $spepos=substr($seq,$i,1);
			$hash1{$spe}->{$i}=$spepos;
			if ($cleaner{$spe} && ($spepos ne $refpos) ) {
				$hash2{$i}++;
			}
		}
	}
}

sub Print_Nons {
	my ($pos)=@_;
	my $info=$pos."[$hash2{$pos}]: ";
	foreach my $spe (@nodes) {
		my $spepos=$hash1{$spe}->{$pos};
		$info.=$spe."($spepos);";
	}
	$info=~s/\;$//;
	my %hash3;
	foreach my $spe (@facul) {
		my $spepos=$hash1{$spe}->{$pos};
		$hash3{$spepos}++;
	}
	my $nb=keys %hash3;
	if ( $nb==1 && ($hash1{'Tbif'}->{$pos} ne $hash1{'Ncel'}->{$pos}) && ($hash1{'Smel'}->{$pos} ne $hash1{'Cund'}->{$pos}) ) {
		print "$orth\t$info\n";
	}
}
