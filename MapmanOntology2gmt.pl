#!/usr/bin/perl -w

use strict;
use 5.010;
use YAML;
use Getopt::Long;
use Data::Printer;
use File::Basename;

my $usage = <<USAGE;
SYSNOPSIS
Mapman2gmt.pl mapman.annoation.txt > mapman.gmt
 Options:
   -i     isoform mode, default off;
   
USAGE

my $isoform_mode = '0';
die $usage
  unless GetOptions(
    "i"      => \$isoform_mode,
  );
  
my $in_file = shift;

#read file to information
open my $in_fh, "<", $in_file or die "cannot open file:$!\n";
readline $in_fh;
my %data;
while (<$in_fh>) {
	chomp;
	s/\'//ig;
	my ($bin_code, $name, $identifier, $desc) = split /\t/;
	$data{$bin_code}{name} = $name;
	if ($identifier) {
		my ($gid , $isoform_num) = split /\./, $identifier;
		if ($isoform_mode) {
			$data{$bin_code}{ids}{$identifier} = $desc;		
		}else{
			$data{$bin_code}{ids}{$gid} = $desc;
		}		
	}

}
close $in_fh;

foreach my $bin_code (keys %data) {
	my $name = $data{$bin_code}{name};
	my $col1 = "$bin_code $name";
	my @cols2 = ();
	foreach my $bin (keys %data) {
		if (length $bin >= length $bin_code && $bin =~ /^$bin_code/) {
			my @ids = keys %{$data{$bin}{ids}};
			push @cols2, @ids;
		}
	}
	my %tmp = map {$_ => 1} @cols2;
	@cols2 = map {lc} grep {$_} sort keys %tmp;
	print "$col1\t";
	print join "\t", @cols2;
	print "\n";
}
