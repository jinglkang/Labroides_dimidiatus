#!/usr/bin/perl
use strict;
use warnings;

#node 100379 will also have a filled green rectangle in the middle of the branch, same size as the circle defined above (size is 10)
#100379,1,10,#00ff00,1,0.5

my %clean_clor=(
        'Ldim'  => '#000000',
        'Tbif' => '#999999',
        'Ncel'   => '#e6e6e6',
        'Smel'      => '#999999',
        'Tads'   => '#999999',
        'Lber'    => '#e6e6e6',
        'Cund'   => '#e6e6e6',
        'Spul'     => '#999999'
        );

my %spe_color=(
        'Ldim'  => '#f44336',
        'Tbif' => '#744700',
        'Ncel'   => '#ce7e00',
        'Smel'      => '#8fce00',
        'Tads'   => '#2986cc',
        'Lber'    => '#16537e',
        'Cund'   => '#6a329f',
        'Spul'     => '#c90076'
        );

open FAS, $ARGV[0] or die "can not open $ARGV[0]\n";
while (<FAS>) {
        chomp;
        if (/\>/) {
#                s/\>//;
                next if /NonOR/;
                (my $spe)=$_=~/>(.*?)\_/;
                if ($clean_clor{$spe} && $spe_color{$spe}) {
                        my $clean=$clean_clor{$spe};
                        my $speco=$spe_color{$spe};
                        s/\>//;
                        my $info ="$_,1,10,$clean,1,1";
                        print "$info\n";
                } else {
                        s/\>//;
                        print "$_,1,10,#e6e6e6,1,1\n";
                }
        }
}
