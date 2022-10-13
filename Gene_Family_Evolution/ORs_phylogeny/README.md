# Build the phylogenetic tree for each OR group
## If the predicted pep sequences include "!", remove it; detect the domian; keep the seq if contain "7tm"; and include the Zebrafish ORs and Non-ORs to build a phylogenetic tree to confirm the classification
```temp1.pl
#!/usr/bin/perl
use strict;
use warnings;

my %pep;
my $pepseq="filter.out.1.info.fasta1";
my $pepnm;
open PEP, $pepseq or die "can not open $pepseq\n$!\n";
while (<PEP>) {
        chomp;
        if (/>/) {
                s/>//;
                my @a=split;
                $pepnm=$a[0];
        } else {
                $pep{$pepnm}.=$_;
        }
}

# domain should contain "7tm"
my %pfam;
my @genes;
my $pfamR="filter_pfam.txt";
open PFAM, $pfamR or die "can not open $pfamR\n$!\n";
while (<PFAM>) {
        chomp;
        if (/^#/ || /^\s*$/) {
                next;
        } else {
                my @F=split;
                $pfam{$F[0]}.=$F[6]."\t";
        }
}

# include the zebrafish pep seq for classify
my $Rule="ruler.fa";
my $clas_nm;
my @clas_nms;
my %Class;
open RULE, $Rule or die "can not open $Rule\n$!\n";
while (<RULE>) {
        chomp;
        if (/>/) {
                s/>//;
                $clas_nm=$_;
                push @clas_nms, $clas_nm;
        } else {
                $Class{$clas_nm}.=$_;
        }
}

foreach my $gene (@clas_nms) {
        if ($gene=~/NonOR/ || $gene=~/D\.rerio.*/) {
                print ">$gene\n$Class{$gene}\n";
        }
}

foreach my $name (sort keys %pep) {
        if ($pfam{$name}) {
                print ">$name\n$pep{$name}\n";
        }
}
```

```1.sh
less filter.out.1.info.fasta|perl -alne 's/\!//g;print'>filter.out.1.info.fasta1 # remove "!" in frameshift pep sequences
~/software/PfamScan/pfam_scan.pl -fasta filter.out.1.info.fasta1 -dir /media/HDD/pfam/ >filter_pfam.txt
perl temp1.pl >OR_phy.fa
muscle -in OR_phy.fa -out OR_phy_align.fas
less OR_phy_align.fas|perl -alne 's/\:/|/g;print' >OR_phy_align.fas1
FastTreeMP OR_phy_align.fas1 >OR_phy_align.phy
```

```bash
# Kang@fishlab3 Wed Oct 12 02:22:51 /media/HDD/cleaner_fish/genome/OR_detection/Labrus_bergylta/group
cp /media/HDD/cleaner_fish/genome/OR_detection/Cheilinus_undulatus/group/temp1.pl ./
cp /media/HDD/cleaner_fish/genome/OR_detection/Cheilinus_undulatus/group/1.sh ./
sh 1.sh
```

## Based on the classification, build the phylogenetic tree for each OR
```temp1.pl
#!/usr/bin/perl
use strict;
use warnings;

my %species=(
    'Spottedgar' => 'Spottedgar',
    'Zebrafish' => 'Zebrafish',
    'Medaka'=> 'Medaka',
    'Platyfish'=>'Platyfish',
    'Fugu'=>'Fugu',
    'Stickleback'=>'Stickleback',
    'Spul'=>'S.pulcher',
    'Cund'=>'C.undulatus',
    'Lber'=>'L.bergylta',
    'Tads'=>'T.adspersus',
    'Smel'=>'S.melops',
    'Ncel'=>'N.celidotus',
    'Tbif'=>'T.bifasciatum',
    'Ldim'=>'L.dimidiatus',
    );

my @groups=<*_ORs_group.txt>;
foreach my $group (@groups) {
    (my $spe)=$group=~/(.*?)_ORs_group\.txt/;
#    my $new=$spe."_ORs_group_name.txt";
    my $new_spe=$species{$spe};
    my %hash;
    open GROUP, $group or die "can not open $group\n$!\n";
#    open NEW,  ">$new" or die "can not open $new\n$!\n";
    while (<GROUP>) {
        chomp;
        if (/^#/) {
            next;
        } else {
            my @a=split /\t/;
            my ($str, $grp);
            if ($a[0]=~/-(\D)$/) {
                $str=$1;
            } elsif ($a[0]=~/\|(\D)$/) {
                $str=$1;
            }
            ($grp)=$a[1];
            $hash{$grp}->{$str}++;
            if ($str eq "C") {
                my $name=$new_spe.".$grp\-$hash{$grp}->{$str}";
                $_=~s/\|/\:/g;
                print "$_\t$name\n";
            } else {
                my $name=$new_spe.".$grp\-$str\-$hash{$grp}->{$str}";
                $_=~s/\|/\:/g;
                print "$_\t$name\n";
            }
        }
    }
}
```

```temp8.pl
#!/usr/bin/perl
use strict;
use warnings;

my %pep;
my $ORnm;
my $ORpep="all_ORs_original_predict.fasta";
open ORPEP, $ORpep or die "can not open $!\n$ORpep\n";
while (<ORPEP>) {
        chomp;
        if (/>/) {
                s/^>//;
                my @a=split;
                $ORnm=$a[0];
        } else {
                $pep{$ORnm}.=$_;
        }
}

# include the zebrafish pep seq for classify
my $Rule="Cheilinus_undulatus/group/ruler.fa";
my $clas_nm;
my @clas_nms;
my %Class;
open RULE, $Rule or die "can not open $Rule\n$!\n";
while (<RULE>) {
        chomp;
        if (/>/) {
                s/>//;
                $clas_nm=$_;
                push @clas_nms, $clas_nm;
        } else {
                $Class{$clas_nm}.=$_;
        }
}

my $NonOR;
foreach my $gene (@clas_nms) {
        if ($gene=~/NonOR/) {
                $NonOR.=">$gene\n$Class{$gene}\n";
        }
}

my $dir="ORs_class";
my $fnm="all_ORs_rename.txt";
open FNM, $fnm or die "can not open $!\n$fnm\n";
while (<FNM>) {
        chomp;
        my @a=split /\t/;
        my ($str, $grp);
        $grp=$a[1];
        if ($a[0]=~/\-(\D)$/) {
                $str=$1;
        } elsif ($a[0]=~/\:(\D)$/) {
                $str=$1;
        } else {
                die "no gene structure detected\n";
        }
        my $file=$grp."-".$str.".fasta";
        if (-e "$dir/$file") {
                open FILE, ">>$dir/$file" or die "can not open $!\n$file\n";
                print FILE ">$a[2]\n$pep{$a[0]}\n";
        } else {
                open FILE, ">$dir/$file" or die "can not open $!\n$file\n";
                print FILE "$NonOR";
                print FILE ">$a[2]\n$pep{$a[0]}\n";
        }
}
```

```bash
# prepare the sequences
# kangjingliang@kangjingliangdeMacBook-Pro 三 10 12 17:52:39 ~/Desktop
perl temp1.pl >all_ORs_rename.txt
scp all_ORs_rename.txt Kang@147.8.76.231:/media/HDD/cleaner_fish/genome/OR_detection
# Kang@fishlab3 Wed Oct 12 17:55:25 /media/HDD/cleaner_fish/genome/OR_detection
cat */group/filter.out.1.info.fasta >all_ORs_original_predict.fasta
# kangjingliang@kangjingliangdeMacBook-Pro 三 10 12 17:57:59 ~/Desktop
scp all_ORs_rename.txt Kang@147.8.76.231:/media/HDD/cleaner_fish/genome/OR_detection/
# Kang@fishlab3 Wed Oct 12 18:21:50 /media/HDD/cleaner_fish/genome/OR_detection
mkdir ORs_class
# Kang@fishlab3 Wed Oct 12 22:43:14 /media/HDD/cleaner_fish/genome/OR_detection
perl temp8.pl
```

```/media/HDD/cleaner_fish/genome/OR_detection/ORs_class/temp2.pl
#!/usr/bin/perl
use strict;
use warnings;

# muscle -in Beta-C.fasta -out Beta-C_align.fasta
# trimal -in Beta-C_align.fasta -out Beta-C_align_trim.fasta -gt 0.8 -st 0.001 -cons 60
# perl temp1.pl Beta-C_align_trim.fasta >Beta-C_align_trim_conc.fasta
# fasta2phy.pl Beta-C_align_trim_conc.fasta >Beta-C_align_trim_conc.phy
# scp Beta-C_align_trim_conc.phy kang1234@147.8.76.177:~/genome/ORs_phy

my $fas=$ARGV[0];
(my $name)=$fas=~/(.*)\.fasta/;
my $align=$name."_align.fasta";
my $trim =$name."_align_trim.fasta";
my $conc =$name."_align_trim_conc.fasta";
my $phy  =$name."_align_trim_conc.phy";
system("muscle -in $fas -out $align");
system("trimal -in $align -out $trim -gt 0.8 -st 0.001 -cons 60");
system("perl temp1.pl $trim > $conc");
system("fasta2phy.pl $conc > $phy");
system("rm $align $trim $conc");
system("scp $phy kang1234\@147.8.76.177:~/genome/ORs_phy");
```

```bash
# Delta
# Kang@fishlab3 Wed Oct 12 23:42:01 /media/HDD/cleaner_fish/genome/OR_detection/ORs_class
perl temp2.pl Delta-C.fasta
# (base) kang1234@celia-PowerEdge-T640 Wed Oct 12 23:26:11 ~/genome/ORs_phy
nohup raxmlHPC -T 24 -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s Delta-C_align_trim_conc.phy -o NP_001471.2.NonOR_8,NP_001043.2.NonOR_9,NP_000667.1.NonOR_3,NP_000675.1.NonOR_2,NP_000862.1.NonOR_7,NP_071640.1.NonOR_4,NP_000671.2.NonOR_1,NP_000854.1.NonOR_5,NP_000857.1.NonOR_6 -n Delta-C >Delta-C.process 2>&1 &
# [1] 6894
```

## prepare the data for itol plot
```/media/HDD/cleaner_fish/genome/OR_detection/ORs_class/temp3.pl
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
        'L.bergylta'    => '#e6e6e6',
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
                next if /NonOR/;
                (my $spe)=$_=~/(\D\..*?)\./;
                my $clean=$clean_clor{$spe};
                my $speco=$spe_color{$spe};
                my $info ="$_,1,10,$clean,1,1";
                print "$info\n";
        }
}
```

```/media/HDD/cleaner_fish/genome/OR_detection/ORs_class/temp4.pl
#!/usr/bin/perl
use strict;
use warnings;

#node 100379 will also have a filled green rectangle in the middle of the branch, same size as the circle defined above (size is 10)
#100379,1,10,#00ff00,1,0.5

#leaf label for node 9606 will be displayed in green, bold and twice the regular font size
#9606 label #00ff00 bold 2

my %clean_clor=(
        'L.dimidiatus'  => '#000000',
        'T.bifasciatum' => '#999999',
        'N.celidotus'   => '#e6e6e6',
        'S.melops'      => '#999999',
        'T.adspersus'   => '#999999',
        'L.bergylta'    => '#e6e6e6',
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
                next if /NonOR/;
                (my $spe)=$_=~/(\D\..*?)\./;
                my $clean=$clean_clor{$spe};
                my $speco=$spe_color{$spe};
                my $info ="$_ label $speco bold 2";
                print "$info\n";
        }
}
```

```bash 
# Example
# Kang@fishlab3 Thu Oct 13 15:54:57 /media/HDD/cleaner_fish/genome/OR_detection/ORs_class
perl temp3.pl Beta-C.fasta # For node shape: dataset_symbols_template.txt
perl temp4.pl Beta-C.fasta # For label colour: colors_styles_template.txt
```
