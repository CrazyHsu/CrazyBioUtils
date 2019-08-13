#!/usr/bin/perl -w

use strict;
use 5.010;
use Getopt::Long;
# use Data::Printer;

#usage information
my $usage = <<USAGE;
transpose.pl V1.0, written by corephi
This program is similar to linux bash "cut -f", the differece
is this can not only output by ordered col numbers, but also 
can output by ordered col names
----------------------------------------------------------
Usage: grep_file_by_col.pl -a in.txt -c col -o outfile
 Options:
  -i    input file name, default STDIN
  -o    output file name, default STDOUT

Note: This scrpit can also be used in pipeline.
cat in.txt | transpose.pl > out.txt
     

USAGE
my $in_file        = '';
my $column_str      = '';
my $outfile      = '-';
my $mode         = 'number';
die $usage
  unless GetOptions(
    "i:s" => \$in_file,
    "o:s" => \$outfile,
  );

die "a and b can not get from STDIN at the same time\n" 
	if $in_file && $in_file eq '-';

 
my $out_fh;
if ($outfile && $outfile ne '-') {

	open $out_fh, ">", $outfile or die "cannot open file $outfile:$!\n";

}else {
    $out_fh = *STDOUT;
}


#append content to a file
my $in_fh;
if ($in_file && $in_file ne '-') {
    die "a file does not exists\n" unless -e $in_file;
    open $in_fh, "<", $in_file or die "cannot open file $in_file:$!\n";
}
else {
    $in_fh = *STDIN;
}

my @matrix = ();
my ($row_num, $col_num) = (0, 0);
while (<$in_fh>) {
    s/\r?\n//;
	next unless $_;
    my @lines = split /\t/;
	push @matrix, \@lines;
	$col_num = @lines if @lines > $col_num;
	$row_num ++;
}
foreach my $col (0 .. $col_num - 1) {
	my @lines = ();
	foreach my $row (0 .. $row_num - 1) {
		push @lines, $matrix[$row][$col];
	}
	print $out_fh join "\t", @lines;
	print $out_fh "\n";
}




