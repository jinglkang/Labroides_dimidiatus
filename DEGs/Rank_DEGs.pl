#!/usr/bin/perl
use strict;
use warnings;
use File::Basename; 
use Cwd qw(cwd);
use Getopt::Long 'HelpMessage';

my $usage=<<_EOH_;;
----------------------------------------------------------------------------------------
This script is used to combine the DEGs from FB, HB, MB
Usage:
perl Rank_DEGs.pl --DEG gtf_FB_Interaction_Solo.DEGs.txt --TYP FB

                                        Kang 2022-11-17
----------------------------------------------------------------------------------------
_EOH_
;

GetOptions(
    'DEG:s', \my $DEG,
    'TYP:s', \my $typ
    );

unless ($DEG && $typ) {
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

my %nbs;
my $fbnbm="gtf_read_nb_tpm_FB.txt";
my $hbnbm="gtf_read_nb_tpm_HB.txt";
my $mbnbm="gtf_read_nb_tpm_MB.txt";
&phase_matrix($fbnbm, "FB");
&phase_matrix($hbnbm, "HB");
&phase_matrix($mbnbm, "MB");

my $ano="Gene_annotation.final.txt";
my %ANO;
open ANO, $ano or die "can not open $ano\n";
while (<ANO>) {
    chomp;
    my @a=split /\t/;
    $ANO{$a[0]}=$a[2]."\t".$a[3];
}

my (%solo, %inte);
open DEG, $DEG or die "can not open $DEG\n";
while (<DEG>) {
    chomp;
    my $gene=$_;
    my $solo_nb=$nbs{$gene}->{$typ}->{'solo'};
    my $inte_nb=$nbs{$gene}->{$typ}->{'inte'};
    $solo{$gene}=$solo_nb;
    $inte{$gene}=$inte_nb;
}

my ($sn, $tn); my (%ranksolo, %rankinte);
my (@topsolo, @topinte);
foreach my $gene (sort {$solo{$b} <=> $solo{$a}} keys %solo) {
    $sn++;
    my $numb=$solo{$gene};
    $ranksolo{$gene}=$sn.":".$numb;
    push @topsolo, $gene; #if $sn<=100;
}

foreach my $gene (sort {$inte{$b} <=> $inte{$a}} keys %inte) {
    $tn++;
    my $numb=$inte{$gene};
    $rankinte{$gene}=$tn.":".$numb;
    push @topinte, $gene; # if $tn<=100;
}

foreach my $gene (@topsolo) {
    my $anoo=$ANO{$gene};
    my $solo_nb=$solo{$gene};
    my $inte_nb=$inte{$gene};
    if ($solo_nb > $inte_nb) {
        print "$gene\tsolo:$ranksolo{$gene}\tinte:$rankinte{$gene}\tUP\t$anoo\n";
    } else {
        print "$gene\tsolo:$ranksolo{$gene}\tinte:$rankinte{$gene}\tDOWN\t$anoo\n";
    }
}


sub phase_matrix {
    my ($reads, $type)=@_;
    my %newinds;
    foreach my $ind (keys %inds) {
        my $newind=$ind.$type;
        $newinds{$newind}=$inds{$ind};
#        print "$newind\t$newinds{$newind}\n";
    }
    my @headers;
    open READS, $reads or die "can not open $reads\n";
    while (<READS>) {
        chomp;
        my @a=split;
        if (/^Geneid/) {
            @headers=@a;
        } else {
            my $gene=$a[0];
            my ($st, $it, $sm, $im);
            for (my $i = 1; $i < @a; $i++) {
                if ($newinds{$headers[$i]} eq "solo") {
                    $st+=$a[$i];
                } else {
                    $it+=$a[$i];
                }
            }
            $sm=$st/6; $im=$it/6;
            $sm=sprintf("%.2f",$sm);
            $im=sprintf("%.2f",$im);
            $nbs{$gene}->{$type}={
                'solo'=> $sm,
                'inte'=> $im
            };
        }
    }
    return %nbs;
}
