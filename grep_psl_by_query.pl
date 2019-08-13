#!/usr/bin/perl -w

use strict;
use 5.010;
use Getopt::Long;

#usage information
my $usage = <<USAGE;
grep_psl_by_query V1.0, written by corephi

----------------------------------------------------------
Usage: grep_psl_by_query.pl -i in.psl -l list.txt -o out.psl
 Options:
  -i    in.psl default STDIN
  -l    idlist file name, it must be one line per id
  -r    reverse match, output the different line
  -o    out.psl, default STDOUT


Note: This scrpit can also be used in pipeline

cat in.psl | grep_psl_by_query.pl -i - -l id.list > out.psl
cat id.list | grep_psl_by_query.pl -i in.psl -l - > out.psl


USAGE
my $afile        = '';
my $bfile        = '';
my $outfile      = '';
my $reverse_flag = 0;
my $header_flag  = 0;
die $usage
  unless GetOptions(
    "i:s" => \$afile,
    "l=s" => \$bfile,
    "o:s" => \$outfile,
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
	    open $outfile_fh, ">", $outfile or die "cannot open file $outfile:$!\n";
}else {
    $outfile_fh = *STDOUT;
}

#store bfile in to hash
my %lists;
while (<$bfile_fh>) {
    next if /^#/;
    s/\r?\n//;
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
	next unless /^\d+/;
    my $line = $_;
    s/\r?\n//;
	my ($match, $mismatch, $rep_match, $ambiguious_num, $q_gap_count, $q_gap_bases, $t_gap_count, $t_gap_bases, $strand, $q_name, $q_size, $q_start, $q_end, $t_name, $t_size, $t_start, $t_end, $block_count, $block_size, $q_starts, $t_starts) = split /\t/;	
	
    if ($reverse_flag) {
        print $outfile_fh $line unless exists $lists{$q_name};
    }
    else {
        print $outfile_fh $line if exists $lists{$q_name};
    }

}

