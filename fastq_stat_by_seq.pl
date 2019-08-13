#!/usr/bin/perl
use strict;
use 5.010;

my $file = $ARGV[0];
my $fa_fh;
if ($file && $file ne '-') {
    die "file: $file does not exists\n" unless -e $file;
    open $fa_fh, "<", $file or die "cannot open file $file:$!\n";
}else {
    $fa_fh = *STDIN;
}

my %seqs;
while(<$fa_fh>) {
	chomp;
	if ($. % 4 == 2) {
		if ( exists $seqs{$_}) {
			$seqs{$_}++;
		}else{
			$seqs{$_} = 1;
		}
	}
}
close $fa_fh;

print "SEQ\tCount\n";
foreach my $seq (sort {$seqs{$b} <=> $seqs{$a}} keys %seqs) {
	print "$seq\t$seqs{$seq}\n";
}
