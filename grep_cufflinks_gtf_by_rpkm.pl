#!/usr/bin/perl -w

use strict;
use 5.010;
use Getopt::Long;
use File::Basename;

my $usage = <<USAGE;
SYSNOPSIS
grep_cufflinks_gtf_by_rpkm.pl -i cufflinks.gtf > filterd.gtf
written by corephi, group:276151571
this program is fliter the transcripts by rpkm

 Options:
  -i   in.gtf, default STDIN
  -f   fpkm  threashold to kepp, default 0
  -m   compare mode: 'gt', 'ge', 'ls', 'le'
  -o   out.gtf, default STDOUT
  
Examples:
grep_cufflinks_gtf_by_rpkm.pl cufflinks.gtf > filterd.gtf
grep_cufflinks_gtf_by_rpkm.pl < cufflinks.gtf > filterd.gtf
grep_cufflinks_gtf_by_rpkm.pl -i cufflinks.gtf -o filterd.gtf
cat in.gtf | grep_cufflinks_gtf_by_rpkm.pl >  filterd.gtf

Note:
gtf file must be fix by fix_cufflinks_gtf.pl


USAGE

my $in_gtf = '-';
my $out_gtf = '-';
my $threashold = 0;
my $mode = 'gt';
die $usage
  unless GetOptions(
    "i:s" => \$in_gtf,
	"f:f" => \$threashold,
	"m:s" => \$mode,
    "o:s" => \$out_gtf,
  );
#check mode
if ($mode ne 'gt' && $mode ne "ge" && $mode ne "ls" && $mode ne 'le') {
	die "$mode is not supported\n". $usage;
}

my @in_fhs;
if ($in_gtf && $in_gtf ne '-') {
    open my $in_fh, "<", $in_gtf or die "cannot open file $in_gtf:$!\n";
	push @in_fhs, $in_fh;
}elsif (@ARGV) {
	foreach my $in_file (@ARGV) {
		if (-s $in_file) {
			open my $in_fh, "<", $in_file or die "cannot open file $in_file:$!\n";
			push @in_fhs, $in_fh;		
		}
	}
}else{
	my $in_fh = *STDIN;
	push @in_fhs, $in_fh;
}
 
my $out_fh;
if ($out_gtf ne '-') {
    open $out_fh, ">", $out_gtf or die "cannot open file $out_gtf:$!\n";
}
else {
    $out_fh = *STDOUT;
}


#read the gtf info
foreach my $in_fh (@in_fhs) {
	while (<$in_fh>) {
		my $line = $_;
		if (/FPKM "(.*?)";/) {
			my $fpkm = $1;
			my $print_flag = 0;
			if ($mode eq 'gt') {
				$print_flag = 1 if $fpkm > $threashold;
			}elsif ($mode eq 'ge'){
				$print_flag = 1 if $fpkm >= $threashold;	
			}elsif ($mode eq 'ls'){
				$print_flag = 1 if $fpkm < $threashold;	
			}elsif ($mode eq 'le'){
				$print_flag = 1 if $fpkm <= $threashold;	
			}
			print $out_fh $line if $print_flag;			
		}
	}	
}

close $out_fh;

