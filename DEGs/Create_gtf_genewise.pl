#!/usr/bin/perl -w
use strict;
use warnings;

# Kang@fishlab3 Sat Nov 19 18:14:58 /media/HDD/cleaner_fish/genome/Opsin_new/Labroides_dimidiatus/filtering

# change "chr" to "Scx22uW_"
my %chro;
my $chro="/media/HDD/cleaner_fish/genome/Ldim_genome.info.change.txt";
open CHRO, $chro or die "can not open $chro\n";
while (<CHRO>) {
    chomp;
    my @a=split;
    $chro{$a[1]}=$a[0];
}

# the pep sequences of opsins for phylogeny
my %ORre;
my $fast="3_phy.bla.fa";
open FAST, $fast or die "can not open $fast\n";
while (<FAST>) {
    chomp;
    if (/^>/) {
        s/>//;
        my @a=split;
        (my $geid)=$a[0]=~/Ldim_(.*)/;
        my @b=split /\:/, $a[1];
        my $gnpo="$b[0]:$b[1]:$b[2]";
        my $orie=$b[3];
        my $stru=$b[4];
        $ORre{$gnpo}={
        'geid' => $geid, # gene id
#        'desp' => $anno, # gene description
#        'posi' => $gnpo_new, # position
        'orie' => $orie, # orientation
        'stru' => $stru  # structure
    };
    }
}

my $blas="3_phy.bla";
open BLAS, $blas or die "can not open $blas\n";
while (<BLAS>) {
    chomp;
    my @a=split /\t/;
    my @b=split /\:/, $a[0];
    my $gnpo="$b[0]:$b[1]:$b[2]";
    $ORre{$gnpo}->{'desp'}=$a[-1];
}

$/ = "\/\/";
# parse the gff file: ../genewise/query.fa.bla.solar.besthit.lt250.wise.best.1.gff
my $GFF="../genewise/2_wise.best.1.gff";
my $gene;
open GFF, $GFF or die "can not open $GFF\n";
while (<GFF>) {
    chomp;
    s/^\s+//; s/\s+$//;
    my @a=split /\n/;
    my ($strand, $gnpo);
    if (@a>0) {
        for (my $i = 0; $i < @a; $i++) {
            my @b=split /\t/, $a[$i];
            my $chr=$chro{$b[0]};
            my $info=$b[5]."\t".$b[6]."\t".$b[7];
            $strand=$b[6];
            my @c=split /\:/, $b[-1];
            $gnpo="$c[0]:$c[1]:$c[2]";
            next unless $ORre{$gnpo}->{'geid'};
            next if $ORre{$gnpo}->{'stru'} eq "F";
            my $gene_id=$ORre{$gnpo}->{'geid'};
            my $tran_id=$gene_id.".t1";
            my $ano=$ORre{$gnpo}->{'desp'};
            my $gene_info="gene_id \"$gene_id\"\; gene_description \"$ano\"\;";
            my $tran_info="transcript_id \"$tran_id\"\; gene_id \"$gene_id\"\;";
            
            if ($i==0) {
                if ($strand eq "+") {
                    print "$chr\tGeneWise\tgene\t$b[3]\t$b[4]\t$info\t$gene_info\n";
                    print "$chr\tGeneWise\ttranscript\t$b[3]\t$b[4]\t$info\t$tran_info\n";
                    my $start=$b[3]; my $stop=$b[3]+2;
                    my $info1=".\t+\t0"; #my $info2="transcript_id \"$tran\"\; gene_id \"$gene\"\;";
                    print "$chr\tGeneWise\tstart_codon\t$start\t$stop\t$info1\t$tran_info\n";
                } else {
                    print "$chr\tGeneWise\tgene\t$b[4]\t$b[3]\t$info\t$gene_info\n";
                    print "$chr\tGeneWise\ttranscript\t$b[4]\t$b[3]\t$info\t$tran_info\n";
                    my $start=$b[4]; my $stop=$b[4]+2;
                    my $info1=".\t-\t0"; #my $info2="transcript_id \"$tran\"\; gene_id \"$gene\"\;";
                    print "$chr\tGeneWise\tstop_codon\t$start\t$stop\t$info1\t$tran_info\n";
                }
            } elsif ($b[2] eq "cds") {
                if ($strand eq "+") {
                    print "$chr\tGeneWise\tCDS\t$b[3]\t$b[4]\t$info\t$tran_info\n";
                    print "$chr\tGeneWise\texon\t$b[3]\t$b[4]\t$info\t$tran_info\n";
                } else {
                    print "$chr\tGeneWise\tCDS\t$b[4]\t$b[3]\t$info\t$tran_info\n";
                    print "$chr\tGeneWise\texon\t$b[4]\t$b[3]\t$info\t$tran_info\n";
                }
            } elsif ($b[2] eq "intron") {
                if ($strand eq "+") {
                    print "$chr\tGeneWise\tintron\t$b[3]\t$b[4]\t$info\t$tran_info\n";
                } else {
                    print "$chr\tGeneWise\tintron\t$b[3]\t$b[4]\t$info\t$tran_info\n";
                }
            }
            if ($i==@a-1) {
                if ($strand eq "+") {
                    my $start=$b[3]; my $stop=$b[3]+2;
                    my $info1=".\t+\t0";
                    my $end_d=$b[4]; my $end_s=$b[4]-2;
                    print "$chr\tGeneWise\tstop_codon\t$end_s\t$end_d\t$info1\t$tran_info\n";
                } else {
                    my $start=$b[4]; my $stop=$b[4]+2;
                    my $info1=".\t-\t0";
                    my $end_d=$b[3]; my $end_s=$b[3]-2;
                    print "$chr\tGeneWise\tstart_codon\t$end_s\t$end_d\t$info1\t$tran_info\n";
                }
            }
        }
    }
}
