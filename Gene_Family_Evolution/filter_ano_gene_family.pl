#!/usr/bin/perl -w
use strict;
use warnings;



my $header="Family_id\tOrthogroup_id\t";
my $header1;
my @spes=qw(Spottedgar Zebrafish Medaka Platyfish Fugu Stickleback Semicossyphus_pulcher Cheilinus_undulatus Labrus_bergylta Tautogolabrus_adspersus Symphodus_melops Notolabrus_celidotus Thalassoma_bifasciatum Labroides_dimidiatus);
foreach my $spe (@spes) {
	$header1.=$spe."\t";
}
$header.=$header1;
$header=~s/\s+$//;

my %hash_ref=(
	'Spottedgar' =>'Spottedgar',
	'Zebrafish'=>'Zebrafish',
	'Medaka'=>'Medaka',
	'Platyfish'=>'Platyfish',
	'Fugu'=>'Fugu',
	'Stickleback'=>'Stickleback',
	'Semicossyphus_pulcher'=>'Spul',
	'Cheilinus_undulatus'=>'Cund',
	'Labrus_bergylta'=>'Lber',
	'Symphodus_melops'=>'Smel',
	'Tautogolabrus_adspersus'=>'Tads',
	'Notolabrus_celidotus'=>'Ncel',
	'Thalassoma_bifasciatum'=>'Tbif',
	'Labroides_dimidiatus'=>'Ldim'
	);

my @refs=qw(Spottedgar Zebrafish Medaka Platyfish Fugu Stickleback);

my $orthog="Orthogroups.GeneCount.tsv";
my (@heads, @orths);
my (%nb, %nb_least);
open ORTH, $orthog or die "can not open $orthog\n$!\n";
while (<ORTH>) {
	chomp;
	my @a=split;
	if (/^Orthogroup/) {
		@heads=@a;
	} else {
		push @orths, $a[0];
		for (my $i = 1; $i < @a-1; $i++) {
			$nb{$heads[$i]}->{$a[0]}=$a[$i];
			if ($nb_least{$a[0]}) {
				if ($nb_least{$a[0]}->{'NB'} < $a[$i]) {
					$nb_least{$a[0]}={
						'NB'  => $a[$i],
						'spe' => $heads[$i]
					};
				}
			} else {
				$nb_least{$a[0]}={
					'NB'  => $a[$i],
					'spe' => $heads[$i]
				};
			}
		}
	}
}

my $zeb="/media/HDD/cleaner_fish/Genome_analysis_restart/ano_files/all_qualified.anno.txt";
my %ano;
open ZEB, $zeb or die "can not open $zeb\n";
while (<ZEB>) {
	chomp;
#	s/\s+\[.*\]$//;
	my @a=split /\t/;
	my $name=$a[0];
	my ($gene_nm, $des);
	if (/symbol=\"(.*?)\"\s+Name=\"(.*?)\s+\[(.*?)\]\"$/) {
		($gene_nm, $des)=($1, $2);
		$ano{$name}=$gene_nm."(".$des.")";
	}
}

my $list="Orthogroups.tsv";
my %ano1;
open LIST, $list or die "can not open $list";
while (<LIST>) {
	next if /^Orthogroup/;
	chomp;
	s/\,//g;
    my @a=split;
    my $fm_id=$a[0];
    my (%nb1, %hash);
    if ($nb_least{$fm_id}) {
    	my $ano_spe=$nb_least{$fm_id}->{'spe'};
    	my $ano_spe1  =$hash_ref{$ano_spe};
    	for (my $i = 1; $i < @a; $i++) {
    		(my $spe)=$a[$i]=~/(.*)\_/;
	        $nb1{$spe}++;
    	    if ($spe eq $ano_spe1 && $ano{$a[$i]}) {
        		my $fm_name=$ano{$a[$i]};
            	$hash{$fm_name}++;
            	$ano1{$fm_id}.=$fm_name.";" if $hash{$fm_name}==1;
	        }
    	}
    }
    $ano1{$fm_id}=~s/\;$// if $ano1{$fm_id};
}

my $fmnb_ano="fm_nb_ano.txt";
my $fmnb    ="fm_nb.txt";
open FMNB_ANO, ">$fmnb_ano" or die "can not create $fmnb_ano\n$!\n";
open FMNB,     ">$fmnb"     or die "can not create $fmnb\n$!\n";

print FMNB_ANO "$header\tGene_name(Description)\n";
print FMNB     "Desc\tFamily ID\t$header1\n";

foreach my $orth (@orths) {
	(my $id)=$orth=~/OG(\d+)/;
	$id="FM".$id;
	my $nb_info;
	foreach my $spe (@spes) {
		$nb_info.=$nb{$spe}->{$orth}."\t";
	}
	$nb_info=~s/\s+$//;
	my @nbs;
	foreach my $ref (@refs) {
		my $nb=$nb{$ref}->{$orth};
		if ($nb>=1) {
			push @nbs, $nb;
		}
	}
	my $annot;
	$ano1{$orth}?($annot=$ano1{$orth}):($annot="-");
	print FMNB_ANO "$id\t$orth\t$nb_info\t$annot\n" if @nbs>=4; # \t$nb_least{$orth}->{'spe'}\t$nb_least{$orth}->{'NB'}
	print FMNB     "(null)\t$id\t$nb_info\n" if @nbs>=4; # \t$nb_least{$orth}->{'spe'}\t$nb_least{$orth}->{'NB'}
}
