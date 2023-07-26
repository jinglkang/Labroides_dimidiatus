# Lidm: revision
## Try to run hyphy for a gene
```bash
# hyphy busted --alignment final_alignment.fa --tree spe.tre --multiple-hit 8
# download "BUSTED-MH.bf" from https://github.com/veg/hyphy-analyses/tree/master/BUSTED-MH
# Kang@fishlab3 Tue Jun 06 14:44:16 /media/HDD/white_island/Compevo/orth16_new/paml_files/OG0026786
hyphy BUSTED-MH.bf --alignment final_alignment.fa --tree spe.tre > 1.txt

# Select specific branch
# Select "Common" for positive selection analysis
# (((Ocomp:0.03585542443421290554,((Stickleback:0.06974035184994364922,Fugu:0.09586551543927242236):0.00721328629797507014,(((Daru:0.01122384596377027331,((Acura:0.00363603864681824899,Apoly:0.00990330563475714035):0.00131655138479240194,(Pmol:0.00064302488464996338,Padel:0.00126163438902896332):0.00430604231776857919):0.00585625365821386065):0.00902505599906567452,(Blenny:0.04482948203984413876,((Yaldwyn:0.00234697253924803784,Blueeyed:0.00508064272075202123):0.00799614071893736157,Common{Foreground}:0.00832227325936589837):0.03669599295005335216):0.01248602538660074229):0.00216716981416960200,(Platyfish:0.07805852278148092682,Medaka:0.09036308277882683371):0.01658547148651793798):0.00395400639645793802):0.00500320730692958848):0.08192758751906517589,Zebrafish:0.16783918724349181084):0.09545358185444878518,Spottedgar:0.09545358185444878518);
nohup hyphy BUSTED-MH.bf --alignment final_alignment.fa --tree spe_len.tre --branches Foreground  > 1.txt 2>&1  &

# Check "1.txt" for P value of Likelihood ratio test
# run hyphy in: 1) my workstation; 2) snorlax; 3) hpc
# The orth id for anaysis: the first column of "Genelist_ano.txt"
# The 2,915 single copy genes used for phylogenetic tree construction
# Kang@fishlab3 Tue Jun 06 21:58:45 /media/HDD/cleaner_fish/genome/gene_family_3/OrthoFinder/Results_May09/Orthogroups/fm_no_care_single_copy

# The phylogenetic tree "RAxML_bestTree.fm_no_care_single_copy_concatenated"
# (base) kang1234@celia-PowerEdge-T640 Tue Jun 06 22:05:17 ~/genome/gene_family

# Kang@fishlab3 Tue Jun 06 19:03:54 /media/HDD/cleaner_fish/genome/gene_family_3/OrthoFinder/Results_May09/sub_orth
nohup scp -r paml_input/ jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/Ldim_revision
# ctrl+z
bg
scp Genelist_ano.txt jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/Ldim_revision
# jlkang@hpc2021 Tue Jun 06 19:06:41 /lustre1/g/sbs_schunter/Kang
mkdir Ldim_revision
# jlkang@hpc2021 Tue Jun 06 22:08:23 /lustre1/g/sbs_schunter/Kang/Ldim_revision
vi spe_hyphy.tre
# ((((Fugu:0.11214230077313369627,(Stickleback:0.08183290882938504263,(Spul:0.02616206046091430065,((Cund:0.03376526434538235089,((Smel:0.01856753434553917032,Tads:0.00824419255946880376):0.00447539782457445374,Lber:0.01808109499928191310):0.01658219238652142158):0.00350401472400454616,(Ncel:0.02867723539260350757,(Ldim{Foreground}:0.01545412627275501334,Tbif:0.01913100688699834531):0.01367612238276721785):0.01088509588524225018):0.00345154831846404640):0.01107897686089813309):0.00376303496703858523):0.00518919920431696984,(Platyfish:0.08177360295703335613,Medaka:0.09610027472192766984):0.01683889120830171435):0.09226451000309993100,Zebrafish:0.16679757628902702749):0.10764984973932918699,Spottedgar:0.10764984973932918699);

# jlkang@hpc2021 Tue Jun 06 22:12:47 /lustre1/g/sbs_schunter/Kang/Ldim_revision
module load hyphy
module load perl/5.34.0 perl-lib/5.34.0

###########################################
# Try
# jlkang@hpc2021 Tue Jun 06 22:16:02 /lustre1/g/sbs_schunter/Kang/Ldim_revision/paml_input/OG0000027_OG8
hyphy /lustre1/g/sbs_schunter/Kang/Ldim_revision/BUSTED-MH.bf --alignment final_alignment.fa --tree /lustre1/g/sbs_schunter/Kang/Ldim_revision/spe_hyphy.tre --branches Foreground
###########################################

# jlkang@hpc2021 Tue Jun 06 22:22:59 /lustre1/g/sbs_schunter/Kang/Ldim_revision
split -l 1155 Genelist_ano.txt Genelist_ano_split_
```

## Run hyphy for positive selection of Ldim
```run_hyphy.pl
#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Parallel::ForkManager;

# hyphy /lustre1/g/sbs_schunter/Kang/Ldim_revision/BUSTED-MH.bf --alignment final_alignment.fa --tree /lustre1/g/sbs_schunter/Kang/Ldim_revision/spe_hyphy.tre --branches Foreground
my $list=$ARGV[0]; # The list
my $tree=$ARGV[1]; # The tree
my $spe =$ARGV[2]; # The species for positive selection analysis
my $outd="Hyphy";
unless (-d $outd) {
	mkdir $outd;
}

my @cmds;
open LIST, $list or die "can not open $list\n";
while (<LIST>) {
	chomp;
	my @a=split;
	my $orth=$a[0];
	my $alig="paml_input/$orth/final_alignment.fa"; # the alignment
	my $outl="$outd/$orth"."_".$spe.".txt";
	my $cmd ="hyphy busted --alignment $alig --tree $tree --multiple-hits Double+Triple --starting-points 5 --branches Foreground > $outl";
#	print "$cmd\n";
	push @cmds, $cmd;
}

my $manager = new Parallel::ForkManager(5);
foreach my $cmd (@cmds) {
	$manager->start and next;
    system($cmd);
    $manager->finish;
}
$manager -> wait_all_children;
```

```script_aa.cmd
#!/bin/bash
#SBATCH --job-name=psg_aa        # 1. Job name
#SBATCH --mail-type=BEGIN,END,FAIL    # 2. Send email upon events (Options: NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jlkang@hku.hk     #    Email address to receive notification
#SBATCH --partition=amd               # 3. Request a partition
#SBATCH --qos=normal                  # 4. Request a QoS
#SBATCH --nodes=1                     #    Request number of node(s)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=1G
#SBATCH --time=7-00:00:00             # 7. Job execution duration limit day-hour:min:sec
#SBATCH --output=%x_%j.out            # 8. Standard output log as $job_name_$job_id.out
#SBATCH --error=%x_%j.err             #    Standard error log as $job_name_$job_id.err

# print the start time
date
perl run_hyphy.pl Genelist_ano_split_aa spe_hyphy.tre Ldim
# print the end time
date
```

```bash
# check the p value
# jlkang@hpc2021 Wed Jun 07 09:45:39 /lustre1/g/sbs_schunter/Kang/Ldim_revision/Hyphy
for i in *.txt;do tail -n 1 ${i};done

# run too slowly
# run the left orthologous genes
# jlkang@hpc2021 Wed Jun 07 10:18:46 /lustre1/g/sbs_schunter/Kang/Ldim_revision/Hyphy
ll -h|perl -alne 'if (/\-/){$F[4]=~s/K//;$F[-1]=~s/_Ldim\.txt//;print "$F[-1]" if $F[4]>=6}' > ../finish_orth.txt
perl temp1.pl > left_orth.txt # divide into 12 files
split -l 541 left_orth.txt left_orth_split_

# just run hyphy for the positive selected genes detected by paml
# Kang@fishlab3 Wed Jun 07 18:46:29 /media/HDD/cleaner_fish/genome/gene_family_3/OrthoFinder/Results_May09/sub_orth
vi PSGs_paml.txt
nohup perl run_hyphy.pl PSGs_paml.txt spe_hyphy.tre Ldim > hyphy_process.txt 2>&1 &
# [1] 29676

# cancel all and use the new command; and better to update the version to the latest HYPHY 2.5.51
# revise "run_hyphy.pl"
# hyphy busted --alignment paml_input/OG0000065_OG8/final_alignment.fa --tree spe_hyphy.tre --multiple-hits Double+Triple --starting-points 5 --branches Foreground

# Kang@fishlab3 Thu Jun 08 22:27:32 /media/HDD/cleaner_fish/genome/gene_family_3/OrthoFinder/Results_May09/sub_orth
nohup perl run_hyphy.pl PSGs_paml.txt spe_hyphy.tre Ldim > hyphy_process.txt 2>&1 &
# [1] 14990

# rerun in HPC
split -l 594 Genelist_ano.txt Genelist_ano_split_

sbatch script_aa.cmd
sbatch script_ab.cmd
sbatch script_ac.cmd
sbatch script_ad.cmd
sbatch script_ae.cmd
sbatch script_af.cmd
sbatch script_ag.cmd
sbatch script_ah.cmd
sbatch script_ai.cmd
sbatch script_aj.cmd
sbatch script_ak.cmd
sbatch script_al.cmd

# jlkang@hpc2021 Mon Jun 12 09:01:10 /lustre1/g/sbs_schunter/Kang/Ldim_revision/Hyphy
perl temp1.pl > hyphy_PSGs.txt # Hyphy result
```

```R
# fdr correct
p_apoly<-read.table(file="hyphy_PSGs.txt")
p_apoly$fdr<- p.adjust(p_apoly$V2,method="fdr",length(p_apoly$V2))
write.table(p_apoly, file="hyphy_PSGs_fdr.txt",row.names=F,col.names=F,quote=F,sep="\t")
# 
```

```bash
# jlkang@hpc2021 Mon Jun 12 09:14:22 /lustre1/g/sbs_schunter/Kang/Ldim_revision/Hyphy
less hyphy_PSGs_fdr.txt |perl -alne 'print if $F[-1] <= 0.05'|wc -l # 45 PSGs detected by hyphy
vi paml_PSGs_fdr.txt
# get the overlap between hyphy and paml
perl temp2.pl > PSG_hyphy_paml.txt # 41 common
```

## Convergence
```bash
# Kang@fishlab3 Fri Jun 09 23:06:59 ~/software/paml4.9j
evolver 7 MCaa.dat
```

```MCaa.dat
0        * 0: paml format (mc.paml); 1:paup format (mc.nex)
13147       * random number seed (odd number)

5 10000 5   * <# seqs>  <# sites>  <# replicates>

-1         * <tree length, use -1 if tree below has absolute branch lengths>

(((Human:0.06135, Chimpanzee:0.07636):0.03287, Gorilla:0.08197):0.11219, Orangutan:0.28339, Gibbon:0.42389);

.5 8        * <alpha; see notes below>  <#categories for discrete gamma>
2 dat/mtmam.dat * <model> [aa substitution rate file, need only if model=2 or 3]

0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05
0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05

 A R N D C Q E G H I
 L K M F P S T W Y V

// end of file

=============================================================================
Notes for using the option in evolver to simulate amino acid sequences.
Change values of parameters, but do not delete them.  It is o.k. to add
empty lines, but do not break down the same line into two or more lines.

  model = 0 (poisson), 1 (proportional), 2 (empirical), 3 (empirical_F)
  Use 0 for alpha to have the same rate for all sites.
  Use 0 for <#categories for discrete gamma> to use the continuous gamma
  <aa substitution rate file> can be dayhoff.dat, jones.dat, and so on.
  <aa frequencies> have to be in the right order, as indicated.
=================!! Check screen output carefully!! =====================
```

## Sequences simulation
### codeml: branch lengths, amino acid frequencies and the best shape parameter for variable rates among sites
```bash
# based on the amino acid sequences
# concatenate all AA sequences of orthologous genes
# Translate nucleotide acid to amino acid
# Kang@fishlab3 Sun Jun 11 16:01:27 /media/HDD/cleaner_fish/genome/gene_family_3/OrthoFinder/Results_May09/sub_orth
mkdir simulation; cd simulation
# Kang@fishlab3 Sun Jun 11 16:18:49 /media/HDD/cleaner_fish/genome/gene_family_3/OrthoFinder/Results_May09/sub_orth/simulation
perl extract_all_pep.pl > all_pep.fa # 2765503 anino acid
```

```codeml_est.ctr
seqfile = all_pep.fa        * sequence data file name
outfile = mlc_est               * main result file name
treefile = spe.tre    * tree structure file name
noisy = 9                   * 0,1,2,3,9: how much rubbish on the screen
verbose = 0                 * 1: detailed output, 0: concise output
runmode = 0                 * 0: user tree;  1: semi-automatic;  2: automatic
                            * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise
seqtype = 2                 * 1:codons; 2:AAs; 3:codons-->AAs
CodonFreq = 2               * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
* ndata = 10
clock = 0                   * 0:no clock, 1:clock; 2:local clock
aaDist = 0                  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a * 7:AAClasses
aaRatefile = /home/Kang/software/paml4.9j/dat/jones.dat        * only used for aa seqs with model=empirical(_F)
                            * dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own
model = 2
NSsites = 0                 * models for codons:
                            * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
                            * models for AAs or codon-translated AAs:
                            * 0:poisson, 1:proportional,2:Empirical,3:Empirical+F * 6:FromCodon, 8:REVaa_0, 9:REVaa(nr=189)
                            * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;
                            * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
                            * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
                            * 13:3normal>0
icode = 0                   * 0:universal code; 1:mammalian mt; 2-11:see below
Mgene = 0                   * 0:rates, 1:separate;
fix_kappa = 0               * 1: kappa fixed, 0: kappa to be estimated
kappa = 2                   * initial or fixed kappa
fix_omega = 0               * 1: omega or omega_1 fixed, 0: estimate
omega = .4                  * initial or fixed omega, for codons or codon-based AAs
fix_alpha = 0               * 0: estimate gamma shape parameter; 1: fix it at alpha
alpha =                     * initial or fixed alpha, 0:infinity (constant rate)
Malpha = 0                  * different alphas for genes
ncatG = 3                   * # of categories in dG of NSsites models
fix_rho = 1                 * 0: estimate rho; 1: fix it at rho
rho = 0.                    * initial or fixed rho,   0:no correlation
getSE = 0                   * 0: don't want them, 1: want S.E.s of estimates RateAncestor = 0 * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
Small_Diff = .5e-6
*   cleandata = 0           * remove sites with ambiguity data (1:yes, 0:no)?
* fix_blength = 0           * 0: ignore, -1: random, 1: initial, 2: fixed
method = 0                  * 0: simultaneous; 1: one branch at a time
```


```bash
# Kang@fishlab3 Sun Jun 11 16:25:39 /media/HDD/cleaner_fish/genome/gene_family_3/OrthoFinder/Results_May09/sub_orth/simulation
nohup codeml codeml_est.ctr > codeml_est.process 2>&1 &
# low speed: estimate genes one by one and then for the average

```

## Forget about the simulation, use other guys method
```bash
# (base) kang1234@celia-PowerEdge-T640 Sun Jun 11 22:53:48 ~/genome/paml_new/paml_input
vi Detect_Nons_revision.pl
vi Detect_Nons_all_revision.pl # make sure all the non-cleaner shared the same character
perl Detect_Nons_all_revision.pl > convergent_evo_genes_revision.txt # 38 convergent sites
# plot
# http://vision.hyphy.org/BUSTED to vision the positively selected site
# kangjingliang@kangjingliangdeMacBook-Pro 一  6 12 14:21:28 ~/Desktop
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/Ldim_revision/paml_input/OG0002026_OG1/final_alignment.fa.BUSTED.json ./OG0002026_OG1_ADCY1.BUSTED.json
# kangjingliang@kangjingliangdeMacBook-Pro 一  6 12 14:49:51 ~/Desktop
scp kang1234@147.8.76.177:~/genome/paml_new/paml_input/OG0002026_OG1/final_alignment_pep.fa ~/Documents/2022/Ldim_genome_Restart/PSGs/ADCY1_alignment_pep.fa
```

## Estimate the genome directly by BUSCO
```bash
# Kang@fishlab3 Tue Jun 13 01:15:31 /media/HDD/cleaner_fish/genome/All_genomes
nohup python ~/software/Busco/scripts/run_BUSCO.py -i Labroides_dimidiatus.fasta -o Ldim_genome -l ~/software/Busco/lineage/actinopterygii_odb9 -m geno -c 32 > Busco_Ldim_genome.process 2>&1 &
# [1] 14709
```

## Opsin
```bash
# Kang@fishlab3 Wed Jun 14 23:38:54 /media/HDD/cleaner_fish/genome/Opsin_new
perl temp5.pl > Predict_Opsins_cleaner.fa
cp temp3.pl temp7.pl
perl temp7.pl # Predict_Opsins_cleaner.phy
scp Predict_Opsins_cleaner.phy kang1234@147.8.76.177:~/genome/gene_family/Opsins
# (base) kang1234@celia-PowerEdge-T640 Wed Jun 14 23:47:22 ~/genome/gene_family/Opsins
nohup raxmlHPC -T 24 -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s Predict_Opsins_cleaner.phy -o Fugu_ENSTRUG00000004747,Fugu_ENSTRUG00000013247,Spottedgar_ENSLOCG00000004577 -n Opsins_cleaner >Opsins.process 2>&1 &
# [1] 3874
```


## Glutamate
```bash
## Glutamate receptor
```bash
# do the filtering once again
# 1. the "Glutamate_names.txt" didn't include the name of all glutamate receptors in Uniprot
# 2. the predicted genes are too short
perl run_Fmdetect.pl /media/HDD/cleaner_fish/genome/Glutamates/Query_glutamate.fasta 50 /media/HDD/cleaner_fish/genome/Glutamates/Glutamate_names.txt
# domain detection
Kang@fishlab3 Thu Jun 15 10:58:23 /media/HDD/cleaner_fish/genome/Glutamates/Labroides_dimidiatus/filtering
less 3_phy.bla.fa|perl -alne 's/\!//g;print' >3_phy.bla.fa.2
~/software/PfamScan/pfam_scan.pl -fasta 3_phy.bla.fa.2 -dir /media/HDD/pfam/ > pfam.txt
perl temp2.pl|perl -alne 'my $nb=@F;print if $nb>2'|less # more than 2 domain
perl temp3.pl

# Check query domain
# Kang@fishlab3 Thu Jun 15 10:56:01 /media/HDD/cleaner_fish/genome/Glutamates
cp /media/HDD/cleaner_fish/genome/Glutamates/Labroides_dimidiatus/filtering/temp2.pl temp3.pl
# Check domain for each query
~/software/PfamScan/pfam_scan.pl -fasta Query_glutamate.fasta -dir /media/HDD/pfam/ > query_pfam.txt
perl temp3.pl|grep -i 'zebrafish'
# Check the length
less Query_glutamate.fasta|perl -alne 'if (/>/){$info=$_}else{$len=length($_);$info.="\t".$len;print $info;my $info}'

# Error: Error reading input stream at line 24: Invalid character (!) in sequence
# Kang@fishlab3 Thu Jun 15 12:41:26 /media/HDD/cleaner_fish/genome/Glutamates/Labroides_dimidiatus/genewise
diamond blastp -q 2_wise.best.1.aa -e 1e-5 --sensitive -k 1 -d ~/Desktop/Annotation_database/swiss-prot/uniprot-filtered-reviewed_yes.fasta --out 1.txt
blastp -outfmt 6 -query 2_wise.best.1.aa -out 1.txt -db ~/Desktop/Annotation_database/swiss-prot/uniprot-filtered-reviewed_yes.fasta -num_threads 30
```

```bash
# 1.sh
# based on the previous result
cp /media/HDD/cleaner_fish/genome/Glutamates/Labroides_dimidiatus/filtering/temp3.pl ./
less 3_phy.bla.fa|perl -alne 's/\!//g;print' >3_phy.bla.fa.2
~/software/PfamScan/pfam_scan.pl -fasta 3_phy.bla.fa.2 -dir /media/HDD/pfam/ > pfam.txt
perl temp3.pl > glutamate_revision.txt
```

```bash
# Kang@fishlab3 Thu Jun 15 12:53:25 /media/HDD/cleaner_fish/genome/Glutamates/Cheilinus_undulatus/filtering
cp /media/HDD/cleaner_fish/genome/Glutamates/Symphodus_melops/filtering/1.sh ./
sh 1.sh
```

```bash
# build a phylogeny according these final glutamate receptors
# and remove if many "x"
# Kang@fishlab3 Thu Jun 15 13:18:18 /media/HDD/cleaner_fish/genome/Glutamates/Thalassoma_bifasciatum/filtering
perl temp11.pl > glutamate_revision.fa
```

```bash
# 2.bash
cp /media/HDD/cleaner_fish/genome/Glutamates/Thalassoma_bifasciatum/filtering/temp11.pl ./
perl temp11.pl > glutamate_revision.fa

# Kang@fishlab3 Thu Jun 15 13:39:47 /media/HDD/cleaner_fish/genome/Glutamates
perl temp4.pl
mv Predict_Opsins.fa Predict_Glutamates_revision.fa
mv Predict_Opsins_revison.phy Predict_Glutamates_revision.phy
scp Predict_Glutamates_revision.phy kang1234@147.8.76.177:~/genome/gene_family/Glutamates
# (base) kang1234@celia-PowerEdge-T640 Thu Jun 15 13:46:28 ~/genome/gene_family/Glutamates
nohup raxmlHPC -T 24 -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s Predict_Glutamates_revision.phy -n Glutamates_revison >Glutamates.process 2>&1 &
# [1] 8310
```

```bash
# the number of glutamate receptor
# Kang@fishlab3 Thu Jun 15 14:25:42 /media/HDD/cleaner_fish/genome/Glutamates
perl temp7.pl
###############################################################
# mGluRs
###############################################################
# mGluR group1
# GRM1	Ldim:Tbif:Smel:Tads:Spul:Ncel:Lber:Cund	3:5:4:3:3:3:2:4
# GRM5	Ldim:Tbif:Smel:Tads:Spul:Ncel:Lber:Cund	3:3:2:2:3:3:1:5

# mGluR group2
# GRM2	Ldim:Tbif:Smel:Tads:Spul:Ncel:Lber:Cund	2:2:2:1:1:2:2:2
# GRM3	Ldim:Tbif:Smel:Tads:Spul:Ncel:Lber:Cund	0:2:0:3:1:1:2:1

# mGluR group3
# GRM4	Ldim:Tbif:Smel:Tads:Spul:Ncel:Lber:Cund	5:3:4:5:3:4:4:4
# GRM8	Ldim:Tbif:Smel:Tads:Spul:Ncel:Lber:Cund	5:9:3:4:6:5:4:4
# GRM7	Ldim:Tbif:Smel:Tads:Spul:Ncel:Lber:Cund	3:3:1:3:2:2:1:3


###############################################################
# iGluRs
###############################################################
# GRID
# GRID1	Ldim:Tbif:Smel:Tads:Spul:Ncel:Lber:Cund	2:3:1:2:2:3:2:2
# GRID2	Ldim:Tbif:Smel:Tads:Spul:Ncel:Lber:Cund	1:1:1:1:1:1:2:1

# AMPA
# GRIA4	Ldim:Tbif:Smel:Tads:Spul:Ncel:Lber:Cund	5:5:4:4:4:5:4:5
# GRIA1	Ldim:Tbif:Smel:Tads:Spul:Ncel:Lber:Cund	2:3:3:3:3:1:3:3
# GRIA2	Ldim:Tbif:Smel:Tads:Spul:Ncel:Lber:Cund	1:1:2:3:3:1:2:3
# GRIA3	Ldim:Tbif:Smel:Tads:Spul:Ncel:Lber:Cund	2:3:2:2:3:4:2:3

# Kainate
# GRIK2	Ldim:Tbif:Smel:Tads:Spul:Ncel:Lber:Cund	3:2:4:3:3:3:3:6
# GRIK3	Ldim:Tbif:Smel:Tads:Spul:Ncel:Lber:Cund	1:1:1:1:1:1:1:2
# GRIK1	Ldim:Tbif:Smel:Tads:Spul:Ncel:Lber:Cund	1:2:2:2:1:2:2:2
# GRIK5	Ldim:Tbif:Smel:Tads:Spul:Ncel:Lber:Cund	4:7:1:3:1:1:3:4
# GRIK4	Ldim:Tbif:Smel:Tads:Spul:Ncel:Lber:Cund	1:1:1:1:1:0:1:1

# NMDA
# NMDE1	Ldim:Tbif:Smel:Tads:Spul:Ncel:Lber:Cund	3:3:3:3:3:3:2:5
# NMDE3	Ldim:Tbif:Smel:Tads:Spul:Ncel:Lber:Cund	2:2:3:3:2:4:1:5
# NMDE4	Ldim:Tbif:Smel:Tads:Spul:Ncel:Lber:Cund	2:4:1:2:4:1:3:2
# NMD3A	Ldim:Tbif:Smel:Tads:Spul:Ncel:Lber:Cund	3:3:1:2:3:2:4:3
# NMD3B	Ldim:Tbif:Smel:Tads:Spul:Ncel:Lber:Cund	0:1:1:0:0:0:0:0
# NMDE2	Ldim:Tbif:Smel:Tads:Spul:Ncel:Lber:Cund	4:4:2:4:4:4:2:6
# NMDZ1	Ldim:Tbif:Smel:Tads:Spul:Ncel:Lber:Cund	3:3:1:2:2:0:3:4

# construct the phylogeny according to the classify
perl temp8.pl # Predict_iGluRs_revision.fa; Predict_mGluRs_revision.fa

# Kang@fishlab3 Thu Jun 15 15:43:36 /media/HDD/cleaner_fish/genome/Glutamates
perl temp9.pl Predict_mGluRs_revision.fa Predict_mGluRs # Predict_mGluRs_revison.phy
perl temp9.pl Predict_iGluRs_revision.fa Predict_iGluRs # Predict_iGluRs_revison.phy
scp Predict_iGluRs_revison.phy kang1234@147.8.76.177:~/genome/gene_family/Glutamates
# (base) kang1234@celia-PowerEdge-T640 Thu Jun 15 13:46:28 ~/genome/gene_family/Glutamates
nohup raxmlHPC -T 24 -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s Predict_iGluRs_revison.phy -n iGluRs_revison > iGluRs.process 2>&1 &
# [1] 8880

# Kang@fishlab3 Thu Jun 15 16:07:56 /media/HDD/cleaner_fish/genome/Glutamates
scp Predict_mGluRs_revison.phy jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/Ldim_revision
# jlkang@hpc2021 Thu Jun 15 16:14:04 /lustre1/g/sbs_schunter/Kang/Ldim_revision
sbatch script_aa.cmd
# Submitted batch job 1242921
```

```bash
# 
# Kang@fishlab3 Thu Jun 15 18:39:16 /media/HDD/cleaner_fish/genome/Glutamates/Labroides_dimidiatus/filtering
cp /media/HDD/cleaner_fish/genome/Glutamates/temp6.pl ./
cp /media/HDD/cleaner_fish/genome/Glutamates/temp9.pl ./
perl temp9.pl Predict_iGluRs_revision.fa Predict_iGluRs

# iGluRs contain 'SYTANLAAF' motif
# Kang@fishlab3 Thu Jun 15 19:15:07 /media/HDD/cleaner_fish/genome/Glutamates/Cheilinus_undulatus/filtering
less 3_phy.bla.fa.2 |grep -i 'SYTANLAAF'|wc -l # 26 iGluRs

# mGluRs contain seven transmembrane region
less glutamate_revision.txt|grep -i '7tm_3'|wc -l # 15 mGluRs

less 3_phy.bla.fa.2 |grep -i 'SYTANLAAF'|wc -l
less glutamate_revision.txt|grep -i '7tm_3'|wc -l

# Extract the predicted iGluRs and mGluRs for each species
# Kang@fishlab3 Fri Jun 16 10:54:07 /media/HDD/cleaner_fish/genome/Glutamates
perl temp10.pl # Predict_iGluRs_revision.fa and Predict_mGluRs_revision.fa will be the dir "filtering" for each species
perl temp11.pl Predict_iGluRs_revision.fa iGluRs # output: iGluRs_predicted_all.fasta; iGluRs.phy
perl temp11.pl Predict_mGluRs_revision.fa mGluRs # output: mGluRs_predicted_all.fasta; mGluRs.phy

# summany the number for each genes
# Kang@fishlab3 Fri Jun 16 11:27:02 /media/HDD/cleaner_fish/genome/Glutamates
perl temp7.pl iGluRs.phy
# NMDA
NMDE1	2	2	2	2	2	2	2	2
NMDE2	2	2	1	2	2	1	2	2
NMDE3	2	2	2	2	2	1	2	2
NMDE4	2	2	1	0	2	1	1	1
NMDZ1	2	1	1	2	1	1	0	2

# GRID
GRID2	1	1	1	1	1	2	0	0

# AMPA
GRIA1	1	2	2	1	2	2	1	2
GRIA2	1	1	2	3	2	2	1	2
GRIA3	2	2	2	2	2	2	2	3
GRIA4	2	2	2	2	2	2	2	1

# Kainate
GRIK1	1	2	2	2	1	2	2	2
GRIK2	2	2	2	2	2	2	2	3
GRIK3	1	1	1	1	1	1	1	0
GRIK4	1	1	1	1	1	1	0	1
GRIK5	2	3	1	2	1	1	1	2
# Kang@fishlab3 Fri Jun 16 11:29:32 /media/HDD/cleaner_fish/genome/Glutamates
perl temp7.pl mGluRs.phy
# mGluR group1
GRM1	2	4	3	2	2	1	2	2
GRM5	2	2	1	1	2	0	2	2

# mGluR group2
GRM2	2	2	2	1	1	2	2	2
GRM3	0	1	0	2	1	0	1	1

# mGluR group3
GRM4	2	2	2	2	1	2	2	2
GRM7	2	2	1	3	2	1	1	2
GRM8	4	6	3	3	4	3	4	4


##############################################################################
# Differentially expressed GluRs
# iGluRs
# 1. AMPA
Ldim_g3040	GRIA1	Glutamate receptor 1
Ldim_g19110	GRIA2	Glutamate receptor 2
Ldim_g418	GRIA3	Glutamate receptor 3
Ldim_g5612	GRIA4	Glutamate receptor 4

# 2. GRID
Ldim_g21576	GRID1	Glutamate receptor ionotropic, delta-1

# 3. Kainate
Ldim_g6091	GRIK1	Glutamate receptor ionotropic, kainate 1
Ldim_g686	GRIK1	Glutamate receptor ionotropic, kainate 1
Ldim_g7121	GRIK2	Glutamate receptor ionotropic, kainate 2
Ldim_g13736	GRIK3	Glutamate receptor ionotropic, kainate 3
Ldim_g5526	GRIK4	Glutamate receptor ionotropic, kainate 4
Ldim_g10831	GRIK5	Glutamate receptor ionotropic, kainate 5
Ldim_g13699	GRIK5	Glutamate receptor ionotropic, kainate 5

# 4. NMDA
Ldim_g20661	GRIN1	Glutamate receptor ionotropic, NMDA 1
Ldim_g12359	GRIN2A	Glutamate receptor ionotropic, NMDA 2A
Ldim_g19799	GRIN2A	Glutamate receptor ionotropic, NMDA 2A
Ldim_g12891	GRIN2B	Glutamate receptor ionotropic, NMDA 2B
Ldim_g19424	GRIN2B	Glutamate receptor ionotropic, NMDA 2B
Ldim_g10930	GRIN2D	Glutamate receptor ionotropic, NMDA 2D
Ldim_g16055	GRIN3A	Glutamate receptor ionotropic, NMDA 3A
Ldim_g16062	GRIN3A	Glutamate receptor ionotropic, NMDA 3A
Ldim_g21258	GRIN3A	Glutamate receptor ionotropic, NMDA 3A
Ldim_g16063	GRIN3B	Glutamate receptor ionotropic, NMDA 3B

# mGluRs
# 1. Group1
Ldim_g5621	GRM5	Metabotropic glutamate receptor 5
Ldim_g690	GRM5	Metabotropic glutamate receptor 5
# 2. Group3
Ldim_g15257	GRM4	Metabotropic glutamate receptor 4
Ldim_g17696	GRM4	Metabotropic glutamate receptor 4
Ldim_g18806	GRM7	Metabotropic glutamate receptor 7

# GluR interacting
Ldim_g25970	GRIP1	Glutamate receptor-interacting protein 1
Ldim_g17780	GRIP2	Glutamate receptor-interacting protein 2

# make a heatmap
# kangjingliang@kangjingliangdeMacBook-Pro 五  6 16 14:53:48 ~/Documents/2021/Cleaner_wrasse/gene_expression/enrichment
less Glutamate_DEGs.txt|perl -alne 'if (/^\D+/){print $F[0] unless /^#/}' >Glutamate_DEGs_heatmap_plot.txt
extract_reads_nb --matrix gtf_read_nb_tpm_FB.txt --genes Glutamate_DEGs_heatmap_plot.txt --samples ../coldata_FB.txt >Glutamate_DEGs_FB_reads_nb_tpm.txt
extract_reads_nb --matrix gtf_read_nb_tpm_HB.txt --genes Glutamate_DEGs_heatmap_plot.txt --samples ../coldata_HB.txt >Glutamate_DEGs_HB_reads_nb_tpm.txt
extract_reads_nb --matrix gtf_read_nb_tpm_MB.txt --genes Glutamate_DEGs_heatmap_plot.txt --samples ../coldata_MB.txt >Glutamate_DEGs_MB_reads_nb_tpm.txt

mv Glutamate_DEGs_FB_reads_nb_tpm.txt Glutamate_DEGs_HB_reads_nb_tpm.txt Glutamate_DEGs_MB_reads_nb_tpm.txt ~/Documents/2021/Cleaner_wrasse/gene_expression/Gene_families/
##############################################################################
```
