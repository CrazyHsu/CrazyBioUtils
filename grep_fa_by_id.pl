#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use 5.010;

my $usage = <<USAGE;
grep_fa_by_id.pl -i in.fa -l id.txt > out.fa
written by corephi, group:276151571
----------------------------------------------------------
More scripts? Join "bioinformatics*CN" QQ group: 276151571
If you have any probel or suggestions about this program, 
please mail to: hpyu\@genetics.ac.cn
----------------------------------------------------------
Usage: grep_gtf_by_id.pl -i transcripts.gtf -l id.txt
 Options:
  -i    fa file name, the first line must be id, default STDIN
  -o    out.fa, default STDOUT
  -r    reverse match, output the different line
  -l    idlist file name, it must be one line per id
  

USAGE

my $in_fasta  = '';
my $id_list_file = '';
my $out_fasta = '-';
my $reverse_flag = 0;
die $usage
  unless GetOptions(
    "i:s" => \$in_fasta,
    "l=s" => \$id_list_file,
    "o=s" => \$out_fasta,
    'r'   => \$reverse_flag,
  );
  
die "-f and -l can not get from STDIN at the same time\n" 
	if $in_fasta && $id_list_file && $in_fasta eq '-' && $id_list_file eq '-';
die "-f and -l can not get from STDIN at the same time\n" 
	unless $in_fasta || $id_list_file;


my $id_list_file_fh;
if ($id_list_file && $id_list_file ne '-') {
    open $id_list_file_fh, "<", $id_list_file or die "cannot open file $id_list_file:$!\n";
}else{

		$id_list_file_fh = *STDIN;				
}


#store id_list_file in to hash
my %lists;
while (<$id_list_file_fh>) {
    next if /^#/;
    s/\r?\n$//;
    my @ids = split /\t/;
    foreach my $id (@ids) {
        $lists{$id} = 1;
    }
}

#append content to a file
my $in_fasta_fh;
if ($in_fasta && $in_fasta ne '-') {
    die "a file does not exists\n" unless -e $in_fasta;
    open $in_fasta_fh, "<", $in_fasta or die "cannot open file $in_fasta:$!\n";
}
else {
	if (@ARGV) {
		if ($ARGV[0] eq '-') {
			$in_fasta_fh = *STDIN;
		}elsif (-e $ARGV[0]) {
			open $in_fasta_fh, "<", $ARGV[0] or die "cannot open file $ARGV[0]:$!\n";		
		}else{
			die "$ARGV[0] does not exists\n";
		}
	}else{
		$in_fasta_fh = *STDIN;				
	}
}

my $out_fh;
if ($out_fasta ne '-') {
    open $out_fh, ">", $out_fasta or die "cannot open file $out_fasta:$!\n";
}
else {
    $out_fh = *STDOUT;
}

my $out_flag = 0;
while (<$in_fasta_fh>) {
    s/\r?\n//;
	if (/^>/) {
		my ($tmp,$id) = split />|\s+/;
		if ($reverse_flag) {
			$out_flag = exists $lists{$id} ? 0 : 1;
		}else{
			$out_flag = exists $lists{$id} ? 1 : 0;			
		}
	}else {
	}
	print $out_fh "$_\n" if $out_flag;
}
