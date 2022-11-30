#!/usr/bin/perl
use strict;
use warnings;

#node 100379 will also have a filled green rectangle in the middle of the branch, same size as the circle defined above (size is 10)
#100379,1,10,#00ff00,1,0.5

my %clean_clor=(
        'L.dimidiatus'  => '#000000',
        'T.bifasciatum' => '#999999',
        'N.celidotus'   => '#e6e6e6',
        'S.melops'      => '#999999',
        'T.adspersus'   => '#999999',
        'L.bergylta'    => '#999999',
        'C.undulatus'   => '#e6e6e6',
        'S.pulcher'     => '#999999'
        );

my %spe_color=(
        'L.dimidiatus'  => '#f44336',
        'T.bifasciatum' => '#744700',
        'N.celidotus'   => '#ce7e00',
        'S.melops'      => '#8fce00',
        'T.adspersus'   => '#2986cc',
        'L.bergylta'    => '#16537e',
        'C.undulatus'   => '#6a329f',
        'S.pulcher'     => '#c90076'
        );

open FAS, $ARGV[0] or die "can not open $ARGV[0]\n";
while (<FAS>) {
        chomp;
        if (/\>/) {
                s/\>//;
                if (/NonOR/) {
                        s/\>//;
                        print "$_,1,10,#e6e6e6,1,1\n";
                } else {
                        (my $spe)=$_=~/(\D\..*?)\./;
                        my $clean=$clean_clor{$spe};
                        my $speco=$spe_color{$spe};
                        my $info ="$_,1,10,$clean,1,1";
                        print "$info\n";
                }
        }
}
