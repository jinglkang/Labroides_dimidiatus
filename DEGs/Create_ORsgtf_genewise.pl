#!/usr/bin/perl -w
use strict;
use warnings;

# Kang@fishlab3 Fri Nov 18 14:29:18 /media/HDD/cleaner_fish/genome/OR_detection/Cleaner_wrasse/group
my $filt="filter.out.1.info";
my %ANOT;
open FILT, $filt or die "can not open $filt\n";
while (<FILT>) {
    chomp;
    my @a=split /\t/;
    my ($gnpo, $stru)=$a[0]=~/(.*)\-(\D)/;
    $ANOT{$gnpo}={
        'genm' => $a[2],
        'geid' => $a[1],
    };
}

# The rename txt of all ORs for all species: ~/Documents/2022/Ldim_genome_Restart/ORs/all_ORs_spes_rename.txt
my $ORre="../../all_ORs_spes_rename.txt";
my %ORre;
open ORRE, $ORre or die "can not open $ORre\n";
while (<ORRE>) {
    chomp;
    my @a=split;
    next unless /L\.dimidiatus/;
    my $genm=$a[2]; # the new gene name

    # $gnpo for the gene location in genome
    # $stru for the "intact; partial" information of predicted OR genes
    my ($gnpo, $stru)=$a[0]=~/(.*)\-(\D)/;
    my $anno=$ANOT{$gnpo}->{'genm'};
    my $geid=$ANOT{$gnpo}->{'geid'};
    my @b=split /\:/, $a[0];
    # the orientation
    my ($orie, $gnpo_new);
    if ($b[1]<$b[1]) {
        $orie="+";
        $gnpo_new=$b[0]."\t".$b[1]."\t".$b[2];
    } else {
        $orie="-";
        $gnpo_new=$b[0]."\t".$b[2]."\t".$b[1];
    }

    $ORre{$gnpo}={
        'genm' => $genm, # gene name
        'geid' => $geid, # gene id
        'desp' => $anno, # gene description
        'posi' => $gnpo_new, # position
        'orie' => $orie, # orientation
        'stru' => $stru  # structure
    };
}

$/ = "\/\/";
# parse the gff file: ../genewise/query.fa.bla.solar.besthit.lt250.wise.best.1.gff
my $GFF="../genewise/query.fa.bla.solar.besthit.lt250.wise.best.1.gff";
my $gene;
open GFF, $GFF or die "can not open $GFF\n";
while (<GFF>) {
    chomp;
    s/^\s+//; s/\s+$//;
    my @a=split /\n/;
    my ($strand, $id);
    if (@a>0) {
        for (my $i = 0; $i < @a; $i++) {
            my @b=split /\t/, $a[$i];
            my $info=$b[5]."\t".$b[6]."\t".$b[7];
            $strand=$b[6];
            my @c=split /\-/, $b[-1];
            $id=$c[0];
            next unless $ORre{$id}->{'genm'} && ($ORre{$id}->{'stru'} eq "C");
            my $gene_id=$ORre{$id}->{'genm'};
            my $tran_id=$gene_id.".t1";
            my $gene_info="gene_id \"$gene_id\"";
            my $tran_info="transcript_id \"$tran_id\"\; gene_id \"$gene_id\"\;";
            
            if ($i==0) {
                if ($strand eq "+") {
                    print "$b[0]\tGeneWise\tgene\t$b[3]\t$b[4]\t$info\t$gene_info\n";
                    print "$b[0]\tGeneWise\ttranscript\t$b[3]\t$b[4]\t$info\t$tran_info\n";
                    my $start=$b[3]; my $stop=$b[3]+2;
                    my $info1=".\t+\t0"; #my $info2="transcript_id \"$tran\"\; gene_id \"$gene\"\;";
                    print "$b[0]\tGeneWise\tstart_codon\t$start\t$stop\t$info1\t$tran_info\n";
#                    print "$b[0]\tGeneWise\tCDS\t$b[3]\t$b[4]\t$info\t$tran_info\n";
#                    print "$b[0]\tGeneWise\texon\t$b[3]\t$b[4]\t$info\t$tran_info\n";
#                    my $end_d=$b[4]; my $end_s=$b[4]-2;
#                    print "$b[0]\tGeneWise\tstop_codon\t$end_s\t$end_d\t$info1\t$tran_info\n";
                } else {
                    print "$b[0]\tGeneWise\tgene\t$b[4]\t$b[3]\t$info\t$gene_info\n";
                    print "$b[0]\tGeneWise\ttranscript\t$b[4]\t$b[3]\t$info\t$tran_info\n";
                    my $start=$b[4]; my $stop=$b[4]+2;
                    my $info1=".\t-\t0"; #my $info2="transcript_id \"$tran\"\; gene_id \"$gene\"\;";
                    print "$b[0]\tGeneWise\tstop_codon\t$start\t$stop\t$info1\t$tran_info\n";
#                    print "$b[0]\tGeneWise\tCDS\t$b[4]\t$b[3]\t$info\t$tran_info\n";
#                    print "$b[0]\tGeneWise\texon\t$b[4]\t$b[3]\t$info\t$tran_info\n";
#                    my $end_d=$b[3]; my $end_s=$b[3]-2;
#                    print "$b[0]\tGeneWise\tstart_codon\t$end_s\t$end_d\t$info1\t$tran_info\n";
                }
            } elsif ($b[2] eq "cds") {
                if ($strand eq "+") {
                    print "$b[0]\tGeneWise\tCDS\t$b[3]\t$b[4]\t$info\t$tran_info\n";
                    print "$b[0]\tGeneWise\texon\t$b[3]\t$b[4]\t$info\t$tran_info\n";
                } else {
                    print "$b[0]\tGeneWise\tCDS\t$b[4]\t$b[3]\t$info\t$tran_info\n";
                    print "$b[0]\tGeneWise\texon\t$b[4]\t$b[3]\t$info\t$tran_info\n";
                }
            } elsif ($b[2] eq "intron") {
                if ($strand eq "+") {
                    print "$b[0]\tGeneWise\tintron\t$b[3]\t$b[4]\t$info\t$tran_info\n";
                } else {
                    print "$b[0]\tGeneWise\tintron\t$b[3]\t$b[4]\t$info\t$tran_info\n";
                }
            }
            if ($i==@a-1) {
                if ($strand eq "+") {
                    my $start=$b[3]; my $stop=$b[3]+2;
                    my $info1=".\t+\t0";
                    my $end_d=$b[4]; my $end_s=$b[4]-2;
                    print "$b[0]\tGeneWise\tstop_codon\t$end_s\t$end_d\t$info1\t$tran_info\n";
                } else {
                    my $start=$b[4]; my $stop=$b[4]+2;
                    my $info1=".\t-\t0";
                    my $end_d=$b[3]; my $end_s=$b[3]-2;
                    print "$b[0]\tGeneWise\tstart_codon\t$end_s\t$end_d\t$info1\t$tran_info\n";
                }
            }
        }
    }
}
