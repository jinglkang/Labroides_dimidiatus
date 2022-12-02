#!/usr/bin/perl
use strict;
use warnings;
use Array::Utils qw(:all);

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
my %hash1;
my @poss;
my @cleas=qw(Smel Tads Lber Tbif Spul Ldim);
my @nocls=qw(Spottedgar Zebrafish Medaka Platyfish Fugu Stickleback Cund Ncel);
my @aspes=qw(Spottedgar Zebrafish Medaka Platyfish Fugu Stickleback Cund Ncel Smel Tads Lber Tbif Spul Ldim);

&Build_pos_hash(\@aspes);

sub Build_pos_hash {
        my ($grp)=@_;
        my @grp=@{$grp};
        my $len;
        foreach my $spe (@grp) {
                my $seq=$seq{$spe};
                $len=length($seq);
                for (my $i = 0; $i < $len; $i++) {
                        my $spepos=substr($seq,$i,1);
                        $hash1{$spe}->{$i}=$spepos;
                }
        }

        for (my $i = 0; $i < $len; $i++) {
                my %hash2;
                my $pos=$i;
                my $newp=$pos+1;
                my $info=$newp.":";
                foreach my $spe (@aspes) {
                        my $spepos=$hash1{$spe}->{$pos};
                        $info.=$spe."($spepos);";
                }

                my (@cleas_pos, @nocls_pos);
                foreach my $spe (@cleas) {
                        my $spepos=$hash1{$spe}->{$pos};
                        $hash2{$spepos}++;
                        push @cleas_pos, $spepos;
                }
                foreach my $spe (@nocls) {
                        my $spepos=$hash1{$spe}->{$pos};
                        push @nocls_pos, $spepos;
                }
                my @isect = intersect(@cleas_pos, @nocls_pos);
                my $numb=keys %hash2;
                unless (@isect) {
                        print "$orth\t$numb\t$info\n";
                }
        }
}
