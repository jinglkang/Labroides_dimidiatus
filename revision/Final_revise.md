## Squences simmulaion for convergence estimation
```bash
# estimate the average codons length of paml input
# Kang@fishlab3 Mon Jul 24 12:03:09 /media/HDD/cleaner_fish/genome/gene_family_3/OrthoFinder/Results_May09/sub_orth/paml_input
wc -l final_orth_input_paml.txt # 6353 orthologous paml input
perl temp2.pl
# print "$nm\t$total\t$aver\n";
# 6353	7720731	1215 # 1215/3=405; average length of amino acid sequences of all paml input orthologous genes
```

```MCaa.dat
 0        * 0: paml format (mc.paml); 1:paup format (mc.nex)
13147       * random number seed (odd number)

14 405 6353   * <# seqs>  <# sites>  <# replicates>

-1         * <tree length, use -1 if tree below has absolute branch lengths>

((((Fugu: 0.127813, (Stickleback: 0.096377, (Spul: 0.030647, ((Cund: 0.040846, ((Smel: 0.022863, Tads: 0.010059): 0.005857, Lber: 0.021824): 0.022094): 0.005489, (Ncel: 0.035071, (Ldim: 0.018691, Tbif: 0.024406): 0.017995): 0.014671): 0.005980): 0.019867): 0.007304): 0.009708, (Platyfish: 0.094498, Medaka: 0.109604): 0.019513): 0.121733, Zebrafish: 0.176199): 0.000005, Spottedgar: 0.229246);

.37114 3        * <alpha; see notes below>  <#categories for discrete gamma>
2 /home/Kang/software/paml4.9j/dat/jones.dat * <model> [aa substitution rate file, need only if model=2 or 3]

0.06499 0.05651 0.03858 0.05233 0.02262 0.04549 0.06706 0.06060 0.02632 0.04729
0.10194 0.05840 0.02517 0.04017 0.05174 0.07847 0.05283 0.01278 0.03099 0.06572

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

```bash
# Kang@fishlab3 Mon Jul 24 14:03:38 /media/HDD/cleaner_fish/genome/gene_family_3/OrthoFinder/Results_May09/sub_orth/simulation/Evolver_2step
evolver 7
# EVOLVER in paml version 4.9j, February 2020
# Reading options from data file MCaa.dat
# Kang@fishlab3 Mon Jul 24 14:10:21 /media/HDD/cleaner_fish/genome/gene_family_3/OrthoFinder/Results_May09/sub_orth/simulation/Evolver_2step
less ancestral.txt # Ancestral sequences generated during simulation (check against mc.txt)
less mc.txt # Species sequences generated during simulation

# make a directory for each simulate pep sequences
perl temp1.pl
```

```temp1.pl
#!/usr/bin/perl
use strict;
use warnings;

my $seqs="mc.txt";
my ($nb, $nb2, $dir, $file);
my $phy="pep.phy";
open SEQS, $seqs or die "can not open $seqs\n";
while (<SEQS>) {
  chomp;
  if (/^\D+/) {
    $nb++;
    $nb2++ if $nb==1;
    $dir="Replicate".$nb2;
    mkdir $dir;
    $file=$dir."/".$phy;
    open FILE, ">>$file";
    my @a=split;
    my $info;
    for (my $i = 1; $i < @a; $i++) {
      $info.=$a[$i];
    }
    if ($nb==1) {
      print FILE "14 405\n";
      print FILE "$a[0]   $info\n";
    } elsif ($nb<14) {
      print FILE "$a[0]   $info\n";
    } elsif ($nb==14) {
      print FILE "$a[0]   $info\n";
      $nb=0;
    }
  }
}
```

```bash
# Ancestral State Inference
# Kang@fishlab3 Tue Jul 25 22:38:40 /media/HDD/cleaner_fish/genome/gene_family_3/OrthoFinder/Results_May09/sub_orth/simulation/Evolver_2step/Replicate1435
cp ../spe.tre ./
codeml Ancestral.ctr
# Make sure to separate the sequence from its name by 2 or more spaces
```

```Ancestral.ctr
seqfile = pep.phy
CodonFreq = 0
clock = 0
aaDist = 0
noisy = 0
Mgene = 0
RateAncestor = 1
fix_alpha = 1
alpha = .37114
verbose = 1
seqtype = 2
aaRatefile = /home/Kang/software/paml4.9j/dat/jones.dat 
fix_kappa = 0
kappa = 2
Malpha = 0
ncatG = 3
fix_rho = 1
rho = 0.
getSE = 0
Small_Diff = .5e-6
cleandata = 0
fix_blength = 0
method = 0
runmode = 0
model = 2
NSsites = 5
fix_omega = 1
omega = 1.2
outfile = 1.txt
treefile = spe.tre
icode = 0
estFreq = 0
```

```bash
# Ancestral State Inference
# Kang@fishlab3 Wed Jul 26 10:10:51 /media/HDD/cleaner_fish/genome/gene_family_3/OrthoFinder/Results_May09/sub_orth/simulation/Evolver_2step
nohup perl temp2.pl > estimate_ancestral.process 2>&1 &
# [1] 17603
# combine the ancestral sequences and species sequences to the file "pep_nodes.phy"
perl temp3.pl
```

```temp2.pl
#!/usr/bin/perl
use strict;
use warnings;
use Parallel::ForkManager;

my @dirs=<Replicate*>;
my @cmds;
foreach my $dir (@dirs) {
    my $cmd1="cp Ancestral.ctr spe.tre $dir";
    system($cmd1);
    my $cmd2="chdir $dir;codeml Ancestral.ctr;chdir ../";
    push @cmds, $cmd2;
}

my $manager = new Parallel::ForkManager(18);
foreach my $cmd (@cmds) {
    $manager->start and next;
    system($cmd);
    $manager->finish;
}
$manager -> wait_all_children;
```

```temp3.pl
#!/usr/bin/perl
use strict;
use warnings;
use Parallel::ForkManager;

my @dirs=<Replicate*>;
foreach my $dir (@dirs) {
        chdir $dir;
        system("rm pep_nodes.phy");
        system("cp pep.phy pep_nodes.phy");
        open PEP, ">>pep_nodes.phy" or die "can not open pep_nodes.phy\n";
        open RST, "rst" or die "There is no rst in $dir/\n";
        while (<RST>) {
                chomp;
                if (/^(Node\s+#\d+)\s+(.*)/i) {
                        my ($nm, $seq)=($1, $2);
                        $nm =~s/\s+//g;
                        $nm =~s/\#//g;
                        $seq=~s/\s+//g;
                        print PEP "$nm   $seq\n";
                }
        }
        chdir("../");
}
```


```bash
# the convergent sites detected in my study
# Kang@fishlab3 Wed Jul 26 15:46:35 /media/HDD/cleaner_fish/genome/gene_family_3/OrthoFinder/Results_May09/sub_orth/simulation/Evolver_2step
ll|perl -alne 'print $F[-1] if $F[-1]=~/Replicate/;'|perl -alne 's/\/$//;print' > final_orth_input_paml.txt
scp kang1234@147.8.76.177:~/genome/paml_new/paml_input/Detect_Nons_revision.pl ./
scp kang1234@147.8.76.177:~/genome/paml_new/paml_input/Detect_Nons_all_revision.pl ./
perl Detect_Nons_all_revision.pl > convergent_evo_genes_revision.txt

# detect the convergent sites among six cleaner fish species
# the sites in cleaners (C1 - C6), at least three are identical and different with non-cleaners
# Kang@fishlab3 Thu Jul 27 10:32:18 /media/HDD/cleaner_fish/genome/gene_family_3/OrthoFinder/Results_May09/sub_orth/simulation/Evolver_2step
vi Detect_convergences_cleaners.pl
vi Detect_convergences_cleaners_all.pl
perl Detect_convergences_cleaners_all.pl > convergences_cleaners.txt # 12785 convergent sites in total
# random convergences in 12785 convergent sites
vi Detect_random_convergences_cleaners.pl
vi Detect_random_convergences_cleaners_all.pl
perl Detect_random_convergences_cleaners_all.pl > random_convergences_cleaners.txt # 10277 random convergent sites in total (80.3%)
# false convergences in 12785 convergent sites: if the amino acid of node15 is not consistent with spotted gar
vi Detect_false_convergences_cleaners.pl
vi Detect_false_convergences_cleaners_all.pl
perl Detect_false_convergences_cleaners_all.pl > false_convergences_cleaners.txt # 1040 false convergent sites in total (8.1%)
# Incorrect
vi Detect_incorresct.pl
vi Detect_incorresct_all.pl
perl Detect_incorresct_all.pl > incorresct.txt # 113731 incorrest of 7720731 total amino acids (98.5%)

# based on conservative sites among non-cleaners, and cleaners occurs different amino aicds at least one
vi detect_conservative_noncleaners.pl
vi detect_conservative_noncleaners_all.pl
perl detect_conservative_noncleaners_all.pl > conservative_noncleaners.txt # 856 conservative sites in total
# random convergences in 856 convergent sites
vi detect_conservative_noncleaners_random_convergence.pl
vi detect_conservative_noncleaners_random_convergence_all.pl
perl detect_conservative_noncleaners_random_convergence_all.pl > conservative_noncleaners_random_convergence.txt # 687 random (80.2%)
# false convergence in 856 convergent sites
vi detect_conservative_noncleaners_false_convergence.pl
vi detect_conservative_noncleaners_false_convergence_all.pl
perl detect_conservative_noncleaners_false_convergence_all.pl > conservative_noncleaners_flase_convergence.txt # 0 false convergence
# incorrect
vi detect_conservative_noncleaners_incorresct.pl
vi detect_conservative_noncleaners_incorresct_all.pl
perl detect_conservative_noncleaners_incorresct_all.pl > incorresct_conservative_noncleaners.txt # 0 incorrect of 1473948 total amino acid
perl detect_conservative_noncleaners_incorresct_all.pl > conservative_noncleaners_all.txt # keep non-cleaners are the same: 1473948
```
