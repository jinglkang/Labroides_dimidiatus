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
perl Combine_DEGs.pl --FB FB_DEGs_social_GOs_uniq.txt
         --HB HB_DEGs_social_GOs_uniq.txt
         --MB MB_DEGs_social_GOs_uniq.txt

                                        Kang 2022-11-17
----------------------------------------------------------------------------------------
_EOH_
;

GetOptions(
    'FB:s', \my $fb,
    'HB:s', \my $hb,
    'MB:s', \my $mb,
    );

unless ($fb && $hb && $mb) {
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

my %DEGs;
&phase_DEGs($fb, "FB");
&phase_DEGs($hb, "HB");
&phase_DEGs($mb, "MB");
foreach my $DEG (sort keys %DEGs) {
    my $ano;
    foreach my $tiss (qw(FB HB MB)) {
        if ($DEGs{$DEG}->{$tiss}) {
            $ano=$DEGs{$DEG}->{$tiss};
            last;
        }
    }

    my $info;
    foreach my $tiss (qw(FB HB MB)) {
        my $info1;
        my $solo_nb=$nbs{$DEG}->{$tiss}->{'solo'};
        my $inte_nb=$nbs{$DEG}->{$tiss}->{'inte'};
        my $code;
        ($solo_nb > $inte_nb)?($code="UP"):($code="DOWN");
        $info1="solo:".$solo_nb.";inte:".$inte_nb;
        if ($DEGs{$DEG}->{$tiss}) {
            $info.=$tiss."-DEG"."\t$info1;$code\t";
        } else {
            $info.=$tiss."-NONDEG"."\t$info1\t";
        }
    }
    $info=~s/\s+$//;
    print "$DEG\t$ano\t$info\n";
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

sub phase_DEGs {
    my ($DEG, $type)=@_;
    open DEG, $DEG or die "can not open $DEG\n";
    while (<DEG>) {
        chomp;
        my @a=split /\t/;
        $DEGs{$a[0]}->{$type}=$a[1]."\t".$a[2];
    }
}
