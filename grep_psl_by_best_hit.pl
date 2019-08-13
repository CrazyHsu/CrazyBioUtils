#!/usr/bin/perl -w

use strict;
use 5.010;
use Getopt::Long;
use File::Basename;

my $usage = <<USAGE;
SYSNOPSIS
grep_psl_by_best_hit.pl -i in.psl > best_hits.psl
written by corephi, group:276151571
----------------------------------------------------------
More scripts? Join "bioinformatics*CN" QQ group: 276151571
If you have any probel or suggestions about this program, 
please mail to: hpyu\@genetics.ac.cn
----------------------------------------------------------

 Options:
  -i   in.count, default STDIN
  -t   identity_threshold ,default 0
  -o   exprees.count, default STDOUT

  
Examples:
grep_psl_by_best_hit.pl in.psl > best_hits.psl
grep_psl_by_best_hit.pl < in.psl > best_hits.psl
grep_psl_by_best_hit.pl -i in.psl -o best_hits.psl
cat in.psl | grep_psl_by_best_hit.pl >  best_hits.psl

USAGE

my $in_count_file = '-';
my $out_count_file = '-';
my $identity_threshold = 0;

die $usage
  unless GetOptions(
    "i:s" => \$in_count_file,
    "t:f" => \$identity_threshold,
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
 
my %alignments = ();
while(<$in_fh>) {
	my $line = $_;
	chomp;
	
	next unless /^\d+/;
	my ($match, $mismatch, $rep_match, $ambiguious_num, $q_gap_count, $q_gap_bases, $t_gap_count, $t_gap_bases, $strand, $q_name, $q_size, $q_start, $q_end, $t_name, $t_size, $t_start, $t_end, $block_count, $block_size, $q_starts, $t_starts) = split /\t/;	
	
	my $t_locus = "$t_name:$t_start-$t_end/$strand";		

	my $homology = $match / $q_size;
	
	
	if ($homology >= $identity_threshold ) {
		$alignments{$q_name}{$t_locus}{line} = $line;
		$alignments{$q_name}{$t_locus}{homology} = $homology;
	}else{
		next;
	}			
	
}
close $in_fh;

my %best_hits;
foreach my $q_name (keys %alignments) {
	my @t_loci = sort {$alignments{$q_name}{$a}{homology} <=> $alignments{$q_name}{$b}{homology}} 
				keys %{$alignments{$q_name}};
	my $best_hits = $t_loci[-1];
	print $out_fh $alignments{$q_name}{$best_hits}{line};
}

close $out_fh;







