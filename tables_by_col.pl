#!/usr/bin/perl -w

use strict;
use 5.010;
use Getopt::Long;
use Data::Printer;

#usage information
my $usage = <<USAGE;
tables_by_col V1.0, written by corephi

----------------------------------------------------------
Usage: grep_gtf_by_id.pl -i in.txt -c 2 
 Options:
  -i    gtf file name, the first line must be id, default STDIN
  -o    reverse match, output the different line
  -h 	gtf file has a header, default 0;
  -c    columns, 1:3,4:4,5:8, default first column, 1 based.
  

USAGE
my $in_file     = '-';
my $out_file    =  '-';
my $reverse_flag = 0;
my $header_flag  = 0;
my $col_str = '1';
die $usage
  unless GetOptions(
    "i:s" => \$in_file,
    "o=s" => \$out_file,
    "h"   => \$header_flag,
    'r'   => \$reverse_flag,
    'c'   => \$col_str,
  );

  
#groups
my @columns = ();
if ($col_str){
	my @group_strs = split /,/, $col_str;
	foreach my $group (@group_strs) {
		my ($start, $end ) = split /:/, $group ;
		$end = $end ? $end : $start;	
		push @columns, $start .. $end;
	}
}else{
	die "No column speficified\n";
}

@columns = map {$_ -1 } @columns;  

my $in_fh;
if ($in_file && $in_file ne '-') {
    die "in file does not exists\n" unless -e $in_file;
    open $in_fh, "<", $in_file or die "cannot open file $in_file:$!\n";
}
else {
	if (@ARGV) {
		if ($ARGV[0] eq '-') {
			$in_fh = *STDIN;
		}elsif (-e $ARGV[0]) {
			open $in_fh, "<", $ARGV[0] or die "cannot open file $ARGV[0]:$!\n";		
		}else{
			die "$ARGV[0] does not exists\n";
		}
	}else{
		$in_fh = *STDIN;				
	}
}

my $out_fh;
if ($out_file && $out_file ne '-') {
    open $out_fh, ">", $out_file or die "cannot open file $out_file:$!\n";
}
else {
    $out_fh = *STDOUT;
}

	
if ($header_flag) {
    my $header = readline $in_fh;
}

my %data;
while(<$in_fh>) {
	s/\r?\n//;
	my @cols = split /\t/;
	foreach my $col (@columns) {
		if (exists $data{$cols[$col]}) {
			$data{$cols[$col]}++;
		}else{
			$data{$cols[$col]} = 1;	
		}
	}
}

foreach my $key (sort keys %data) {
	my $count = $data{$key};
	print $out_fh "$key\t$count\n";
}
