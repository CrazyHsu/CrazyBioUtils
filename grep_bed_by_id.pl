#!/usr/bin/perl -w

use strict;
use 5.010;
use Getopt::Long;
use Data::Printer;

#usage information
my $usage = <<USAGE;
grep_bed_by_id V1.0, written by corephi

----------------------------------------------------------
Usage: grep_bed_by_id.pl -i transcripts.bed -l id.txt
 Options:
  -i    bed file name, the first line must be id, default STDIN
  -r    reverse match, output the different line
  -a 	bed file has a header, default 0;
  -l    idlist file name, it must be one line per id
  -m    match by 'gid' or 'tid', default 'tid'
  -o    out.bed, default STDOUT
  

USAGE
my $bed_file     = '';
my $id_list_file      = '';
my $reverse_flag = 0;
my $out_file      = '';
my $header_flag  = 0;
my $mode         = '';
die $usage
  unless GetOptions(
    "i:s" => \$bed_file,
    "l=s" => \$id_list_file,
    "a"   => \$header_flag,
    'r'   => \$reverse_flag,
	"o:s" => \$out_file,
    'm=s' => \$mode,
  );
$mode = 'tid' unless $mode;

my $id_list_file_fh;
if ($id_list_file && $id_list_file ne '-') {
    open $id_list_file_fh, "<", $id_list_file or die "cannot open file $id_list_file:$!\n";
}else{
	$id_list_file_fh = *STDIN;				
}


#store id_file in to hash
my %lists;
while (<$id_list_file_fh>) {
    next if /^#/;
    s/\r?\n$//;
    my @ids = split /\t/;
    foreach my $id (@ids) {
        $lists{$id} = 1;
    }
}

my $bed_file_fh;
if ($bed_file && $bed_file ne '-') {
    die "bed file does not exists\n" unless -e $bed_file;
    open $bed_file_fh, "<", $bed_file or die "cannot open file $bed_file:$!\n";
}
else {
	if (@ARGV) {
		if ($ARGV[0] eq '-') {
			$bed_file_fh = *STDIN;
		}elsif (-e $ARGV[0]) {
			open $bed_file_fh, "<", $ARGV[0] or die "cannot open file $ARGV[0]:$!\n";		
		}else{
			die "$ARGV[0] does not exists\n";
		}
	}else{
		$bed_file_fh = *STDIN;	
	}
}


my $out_fh;
if ($out_file && $out_file ne '-') {
    open $out_fh, ">", $out_file or die "cannot open file $out_file:$!\n";
}
else {
    $out_fh = *STDOUT;
}


if ($header_flag) {
    my $header = readline $bed_file_fh;
    print $out_fh $header if $header_flag;
}

while (<$bed_file_fh>) {
    my $line = $_;
    next if /^#/;
    s/\r?\n//;
	my @tmp = split /\t/;
	my $tid = $tmp[3];
	my $gid = get_gid($tid);
	
	my $id = $tid;
	if ($mode eq 'gid') {
		$id = $gid;
	}else{
		$id = $tid;
	}
    if ($reverse_flag) {
        print $out_fh $line unless exists $lists{$id};
    }
    else {
        print $out_fh $line if exists $lists{$id};
    }

}

sub get_gid {
	my $tid = shift;
	my $gid = $tid;
	if ($tid =~ /(.*)\.\d+/) {
		$gid = $1 
	}else{
		$gid = $tid;
	}
}

