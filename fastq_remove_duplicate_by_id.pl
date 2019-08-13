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

my %ids;
my $flag_print;
while(<$fa_fh>) {
	chomp;
	if ($. % 4 == 1) {
		my $header = $1 if /^@(.*)\b/;
		if ( exists $ids{$header}) {
			$flag_print = 0;
		}else{
			$flag_print = 1;
			$ids{$header} = 1;
		}
	}
	print $_,"\n" if $flag_print;
}
close $fa_fh;

