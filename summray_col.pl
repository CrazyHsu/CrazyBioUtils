#!/usr/bin/perl -w

use strict;
use 5.010;
use Getopt::Long;

#usage information
my $usage = <<USAGE;
grep_file_by_id V1.2, written by corephi
This program is equal to linux bash "grep -f", the differece
is: this program use hash to inprove speed.
----------------------------------------------------------
Usage: grep_file_by_id.pl -a afile -b bfile -o outfile
 Options:
  -a    afile name, the first line must be id, default STDIN
  -r    reverse match, output the different line
  -h 	afile has a header, default 0;
  -b    idlist file name, it must be one line per id
		such as:
		AT1G01010
		AT1G01020
		AT1G01030
Note: This scrpit can also be used in pipeline.such as
cat a.txt | grep_file_by_id.pl -b id.list > greped.txt
if a.txt has a header so, use
cat a.txt | grep_file_by_id.pl -b id.list -h > greped.txt
istead.

USAGE
my $afile        = '';
my $bfile        = '';
my $outfile      = '-';
my $reverse_flag = 0;
my $header_flag  = 0;
die $usage
  unless GetOptions(
    "a:s" => \$afile,
    "b=s" => \$bfile,
    "o:s" => \$outfile,
    "h"   => \$header_flag,
    'r'   => \$reverse_flag,
  );

die "a and b can not get from STDIN at the same time\n" 
	if $afile && $bfile && $afile eq '-' && $bfile eq '-';
die "a and b can not get from STDIN at the same time\n" 
	unless $afile || $bfile;

  
my $bfile_fh;
if ($bfile && $bfile ne '-') {
    die "b file does not exists\n" unless -e $bfile;
    open $bfile_fh, "<", $bfile or die "cannot open file $bfile:$!\n";
}
else {
    $bfile_fh = *STDIN;
}

my $outfile_fh;
if ($outfile && $outfile ne '-') {
	if (-e $outfile) {
		warn "$outfile exists, appended...\n";
		open $outfile_fh, ">>", $outfile or die "cannot open file $outfile:$!\n";
	}else{
	    open $outfile_fh, ">", $outfile or die "cannot open file $outfile:$!\n";

	}
}else {
    $outfile_fh = *STDOUT;
}

#store bfile in to hash
my %lists;
while (<$bfile_fh>) {
    next if /^#/;
    s/\r\n/\n/;
    chomp;
    my @ids = split /\t/;
    foreach my $id (@ids) {
        $lists{$id} = 1;
    }
}

#append content to a file
my $afile_fh;
if ($afile && $afile ne '-') {
    die "a file does not exists\n" unless -e $afile;
    open $afile_fh, "<", $afile or die "cannot open file $afile:$!\n";
}
else {
    $afile_fh = *STDIN;
}
if ($header_flag) {
    my $header = readline $afile_fh;
    print $outfile_fh $header if $header_flag;
}

while (<$afile_fh>) {
    my $line = $_;
    next if /^#/;
    s/\r\n/\n/;
    chomp;
    my @lines = split /\t/;
    my $id    = shift @lines;
    if ($reverse_flag) {
        print $outfile_fh $line unless exists $lists{$id};
    }
    else {
        print $outfile_fh $line if exists $lists{$id};
    }

}

