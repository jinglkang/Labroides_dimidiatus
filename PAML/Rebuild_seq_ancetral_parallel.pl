use strict;
use warnings;
use Parallel::ForkManager;

# Run Rebuild_seq_ancetral.pl parallel
# Rebuild_seq_ancetral_parallel.pl
my $fm=$ARGV[0];
my $conf=$ARGV[1]; # the control file that should be copied to the target directory
die "There is no $conf\n" unless (-e $conf);

my @cmds;
open FM, $fm or die "can not open $fm\n";
while (<FM>) {
    chomp; my @a=split;
    my $dir=$a[0];
    system("cp $conf $dir/");
    my $cmd="perl Rebuild_seq_ancetral.pl $dir";
    push @cmds, $cmd;
}

my $manager = new Parallel::ForkManager(18);
foreach my $cmd (@cmds) {
    $manager->start and next;
    system($cmd);
    $manager->finish;
}
$manager -> wait_all_children;
