#!/usr/bin/perl -w

use strict;
use 5.010;
use Getopt::Long;
use File::Basename;

my $usage = <<USAGE;
SYSNOPSIS
grep_gid_by_rpkm.pl -i total.rpkm > out.fpkm
written by corephi, group:276151571
----------------------------------------------------------
More scripts? Join "bioinformatics*CN" QQ group: 276151571
If you have any probel or suggestions about this program, 
please mail to: hpyu\@genetics.ac.cn
----------------------------------------------------------

 Options:
  -i   in.rpkm, default STDIN
  -f   fpkm threashold set as expressed,default 1
  -o   exprees.fpkm, default STDOUT
  -c   columns, 1:3,4:4,5:8, 1 start with the first data
       matrix, default total. ',' - group, ':' - range
  -n   average fpkm group number,default 1
  
Examples:
grep_gid_by_rpkm.pl total.rpkm > out.fpkm
grep_gid_by_rpkm.pl < total.rpkm > out.fpkm
grep_gid_by_rpkm.pl -i total.rpkm -o out.fpkm
cat total.rpkm | grep_gid_by_rpkm.pl >  out.fpkm

USAGE

my $in_fpkm_file = '-';
my $out_fpkm_file = '-';
my $fpkm_threshold = 1;
my $group_str = '';
my $goup_number_threshold = 1;

die $usage
  unless GetOptions(
    "i:s" => \$in_fpkm_file,
    "f:f" => \$fpkm_threshold,
    "c:s" => \$group_str,
    "n:i" => \$goup_number_threshold,
    "o:s" => \$out_fpkm_file,
  );
  
   
my $out_fh;
if ($out_fpkm_file ne '-') {
    open $out_fh, ">", $out_fpkm_file or die "cannot open file $out_fpkm_file:$!\n";
}
else {
    $out_fh = *STDOUT;
}

my @in_fhs;
if ($in_fpkm_file && $in_fpkm_file ne '-') {
    open my $in_fh, "<", $in_fpkm_file or die "cannot open file $in_fpkm_file:$!\n";
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
 
#Check the column number and ,print the header;
my $header = '';
my $col_num = 0;
foreach my $in_fh (@in_fhs) {
	$header = readline $in_fh;
	my @tmp = split /\t/, $header;
	if ($col_num) {
		die "Col number is not the same!" if $col_num != @tmp
	}else{
		$col_num = @tmp;
	}
}
print $out_fh $header;

#groups
my %groups = ();
if ($group_str){
	my @group_strs = split /,/, $group_str;
	foreach my $group (@group_strs) {
		my ($start, $end ) = split /:/, $group ;
		$groups{$group}{start} = $start;
		$groups{$group}{end} = $end ? $end : $start;	
	}
}else{
	my $end = $col_num - 1;
	$groups{"1:$end"}{end} = $end;	
	$groups{"1:$end"}{start} = 1;	
}




foreach my $in_fh (@in_fhs) {
	while(<$in_fh>) {
		my $line = $_;
		chomp;
		my @tmp = split /\t/;
		my $fpkm_group = 0;
		foreach my $group (keys %groups) {
			my $start = $groups{$group}{start};
			my $end = $groups{$group}{end};
			my @numbers = @tmp[$start..$end];
			$fpkm_group ++ if avg(@numbers) >= $fpkm_threshold;
		}
		
		print $out_fh $line if $fpkm_group >= $goup_number_threshold;		
	}
	close $in_fh;
}

close $out_fh;

sub avg {
	my @numbers = @_;
	my $sum = 0;
	$sum += $_ foreach @numbers;
	return  @numbers ?  $sum / @numbers : 0;
}
