#!/usr/bin/perl -w

use strict;
use 5.010;
use Getopt::Long;
# use Data::Printer;
use File::Basename;

#usage information
my $usage = <<USAGE;
merge_file_by_col V1.0.1, written by corephi
This program is used to merge files by specific lines. 
----------------------------------------------------------
More scripts? Join "bioinformatics*CN" QQ group: 276151571
If you have any problems or suggestions about this program, 
please mail to: hpyu\@genetics.ac.cn
----------------------------------------------------------
Usage: merge_file_by_col.pl -a -c 2 -s .txt *.txt > merged.txt

Options:
  -m    missing value, default "NA"
  -a    files containing header, default off
  -c    columns specified to merge, 0 based
  -s    suffix used to calculate basename
  -o    output filename, default STDOUT
  
Note: 
  make sure your first column is the identifer, and all the 
  files has same identifer such as "ID", "GID" etc, so that
  this program can handle the header as well.

USAGE

my $missing_value   = "NA";
my $outfile         = '-';
my $header_flag         = 0;
my $colmn = 0;
my $suffix = '';
die $usage
  unless GetOptions(
    "s:s" => \$suffix,
    "c:i" => \$colmn,
    "m:s" => \$missing_value,
    "a" => \$header_flag,
    "o:s" => \$outfile,
  );

my @files = grep {-s $_} @ARGV;
die "No files speificied\n" unless @files;
my @suffixlist = split /,/, $suffix;
die "colmn number cannot be setted as 0\n" unless $colmn;

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

my %data;
my @basenames = ();
foreach my $file (@files) {
	warn "Reading $file...\n";
	my $basename = basename($file, @suffixlist);
	push @basenames, $basename;
	my %rpkm = read_file($file);
	foreach my $gid (keys %rpkm) {
		$data{$gid}{$basename}{$_} = $rpkm{$gid}{$_}
			foreach keys %{$rpkm{$gid}}
	}
}
warn "Outputt...\n";
print $outfile_fh join "\t", ("GID", @basenames);
print $outfile_fh "\n";

foreach my $gid (keys %data) {
	my @values = ();
	foreach my $basename (@basenames) {
		push @values, exists $data{$gid}{$basename}{$colmn} ? $data{$gid}{$basename}{$colmn}: "NA"
	}
	print $outfile_fh join "\t", ($gid, @values);
	print $outfile_fh "\n";
}


sub read_file {
	my ($file, $col) = @_;
	open my $fh, "<" ,$file;
	readline $fh if $header_flag;
	my %rpkm;
	while (<$fh>) {
		my @lines = split /\t/;
		my $gid = $lines[0];
		$rpkm{$gid}{$_} = $lines[$_] foreach 1..$#lines;
	}
	return %rpkm;
}
