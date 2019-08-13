#!/usr/bin/perl -w

use strict;
use 5.010;
use Getopt::Long;
use File::Basename;

my $usage = <<USAGE;
SYSNOPSIS
grep_gid_by_rpkm.pl -i total.count > out.count
written by corephi, group:276151571
----------------------------------------------------------
More scripts? Join "bioinformatics*CN" QQ group: 276151571
If you have any probel or suggestions about this program, 
please mail to: hpyu\@genetics.ac.cn
----------------------------------------------------------

 Options:
  -i   in.count, default STDIN
  -t   count threashold set as expressed,default 1
  -o   exprees.count, default STDOUT

  
Examples:
grep_gid_by_rpkm.pl total.rpkm > out.count
grep_gid_by_rpkm.pl < total.rpkm > out.count
grep_gid_by_rpkm.pl -i total.rpkm -o out.count
cat total.rpkm | grep_gid_by_rpkm.pl >  out.count

USAGE

my $in_count_file = '-';
my $out_count_file = '-';
my $count_threshold = 10;

die $usage
  unless GetOptions(
    "i:s" => \$in_count_file,
    "t:f" => \$count_threshold,
    "o:s" => \$out_count_file,
  );
  
   
my $out_fh;
if ($out_count_file ne '-') {
    open $out_fh, ">", $out_count_file or die "cannot open file $out_count_file:$!\n";
}
else {
    $out_fh = *STDOUT;
}
my $in_fh;
if ($in_count_file && $in_count_file ne '-') {
    open $in_fh, "<", $in_count_file or die "cannot open file $in_count_file:$!\n";
}else{
	$in_fh = *STDIN;
}
 
#Check the column number and ,print the header;
my $header = readline $in_fh;
print $out_fh $header;

while(<$in_fh>) {
	my $line = $_;
	chomp;
	my @tmp = split /\t/;
	shift @tmp;
	my @kept = grep {$_ >= $count_threshold} @tmp;
	print $line if @kept;		
}
close $in_fh;
close $out_fh;

