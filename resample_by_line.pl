#!/usr/bin/perl
use 5.010001;
use strict;
use warnings;
use File::Basename;
use Algorithm::Numerical::Sample qw(sample);
use Getopt::Long;

my $usage = <<USAGE;
resample_by_line V1.0, written by corephi
----------------------------------------------------------
More scripts? Join "bioinformatics*CN" QQ group: 276151571
If you have any probel or suggestions about this program, 
please mail to: hpyu\@genetics.ac.cn
----------------------------------------------------------
Usage: resample_by_line -s 10 -i infile > outfile
 Options:
  -i    infile name, default STDIN
  -o    output file name
  -h 	infile has a header, default 0;
  -s    sample size, 0 for all, default 0
Note: This scrpit can also be used in pipeline.such as
cat infile | resample_by_line -s 10 > outfile
resample_by_line -s 10 -i infile > outfile
resample_by_line -s 10 -i infile -o outfile
USAGE
my $in_file = "-";
my $out_file = "";
my $sample_size = 0;
my $header = 0;
die $usage
  unless GetOptions(
    "s:i"        => \$sample_size,
    "i:s"        => \$in_file,
    "o:s"        => \$out_file,
    "h"        => \$header,
  );
  
my $in_fh;
if ($in_file && $in_file ne '-') {
    die "in file does not exists\n" unless -e $in_file;
    open $in_fh, "<", $in_file || die "cannot open file $in_file:$!\n";
}else {
    $in_fh = *STDIN;
}

my $out_fh;
if ($out_file) {
    open $out_fh, ">", $out_file || die "cannot create file $out_file:$!\n";
}else {
    $out_fh = *STDOUT;
}

#get the total number of fastq
my $total_num = 0;
while (<$in_fh>) {
	$total_num++;
}
$total_num-- if $header;
close $in_fh;

#sample;
my %sample;
if ( $sample_size == 0 || $sample_size >= $total_num ) {
	$sample{$_} = 1 foreach 1 .. $total_num;
}
else {
	$sample{$_} = 1 foreach sample(
		-set         => [ 1 .. $total_num],
		-sample_size => $sample_size
	);
}

open $in_fh, "<", $in_file || die "cannot open file $in_file:$!\n";
if ($header) {
	my $head = readline $in_fh;
	print $out_fh $head;
}
while (<$in_fh>) {
	my $line_num = $header ? $. - 1 : $. - 1;
	print $out_fh $_ if exists $sample{$line_num};
}
close $in_fh;
close $out_fh;


