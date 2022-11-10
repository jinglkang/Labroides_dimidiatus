#!/usr/bin/perl
use strict;
use warnings;

# Rebuild_seq_ancetral.pl
# RUN it in the dir including all paml-ortho dir
my $dir=$ARGV[0];
chdir $dir;
if (-e "final_alignment.fa") {
	system("cds2pep.pl final_alignment.fa > final_alignment_pep.fa");
	system("codeml Ancestral.crl");
	if (-e "rst") {
		system("cp final_alignment_pep.fa final_alignment_pep_nodes.fa");
		open PEP, ">>final_alignment_pep_nodes.fa" or die "can not open final_alignment_pep_nodes.fa\n";
		open RST, "rst" or die "There is no rst in $dir/\n";
		while (<RST>) {
			chomp;
			if (/^(Node\s+#\d+)\s+(.*)/) {
				my ($nm, $seq)=($1, $2);
				$nm =~s/\s+//g;
				$nm =~s/\#//g;
				$seq=~s/\s+//g;
				print PEP ">$nm\n$seq\n";
			}
		}
	} else {
		print "No rst in $dir/\n";
	}
}
