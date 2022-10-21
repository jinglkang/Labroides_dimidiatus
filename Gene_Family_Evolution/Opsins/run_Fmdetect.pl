#!/usr/bin/perl
use strict;
use warnings;

my %species=(
        'Cheilinus_undulatus'  =>  'Cund',
        'Fugu'  =>  'Fugu',
        'Labroides_dimidiatus' =>  'Ldim',
        'Labrus_bergylta'  =>  'Lber',
        'Medaka'  =>  'Medaka',
        'Notolabrus_celidotus' =>  'Ncel',
        'Platyfish'  =>  'Platyfish',
        'Semicossyphus_pulcher' =>  'Spul',
        'Spottedgar'  =>  'Spottedgar',
        'Stickleback' =>  'Stickleback',
        'Symphodus_melops' =>  'Smel',
        'Tautogolabrus_adspersus'=> 'Tads',
        'Thalassoma_bifasciatum' => 'Tbif',
        'Zebrafish'  =>  'Zebrafish'
        );

my @species=qw(Cheilinus_undulatus Fugu Labroides_dimidiatus Labrus_bergylta Medaka Notolabrus_celidotus Platyfish Semicossyphus_pulcher Spottedgar Stickleback Symphodus_melops Tautogolabrus_adspersus Thalassoma_bifasciatum Zebrafish);
my $swiss="~/Desktop/Annotation_database/swiss-prot/uniprot-filtered-reviewed_yes.fasta";
my $query=$ARGV[0]; # /media/HDD/cleaner_fish/genome/Protocadherin_gamma/query_protein.fasta
my $len=$ARGV[1]; # 579
my $name=$ARGV[2];
foreach my $spe (sort keys %species) {
    my $dir  ="/media/HDD/cleaner_fish/genome/All_genomes";
    my $short=$species{$spe};
    my $fasta="$spe.fasta";
    my $genom="$dir/$fasta";
    my $pep  ="$dir/longest_pep/$fasta";
    #    my $name ="/media/HDD/cleaner_fish/genome/Opsin_new/Opsin_names.txt";
    system("perl Fmdetect.pl --genome $genom --species $spe --query $query --len $len --uniprot $swiss --shortnm $short --quenm $name\n");
}
