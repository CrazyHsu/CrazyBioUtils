#!/usr/bin/perl -w

use strict;
use 5.010;

use Bio::DB::Sam;
use Getopt::Long;

my $usage = <<USAGE;
Usage: grep_bam_by_read_name.pl [options] -i bamfile -l readnames.lst 

SYSNOPSIS
grep_bam_by_read_name.pl [options] -i bamfile -l readnames.lst  

 Options:
   -i         input bam
   -l         readnames fils, one name per line
   -r         reverse match, output the different line
   -o         output prefix of the reslut 

USAGE

#pare the arguments
my $bam_file = '';
my $id_list_file = '';
my $prefix = 'out';
my $min_intron = 60;
my $reverse_flag = 0;

die $usage unless 
GetOptions ("i=s" => \$bam_file,
			"l=s" => \$id_list_file,
			'r'   => \$reverse_flag,
			"o|out-prefix:s" => \$prefix,
			);
			
my $result_bam = $prefix.".filtered.bam";
#open file to proceed;
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

#open the in bam file
my $inbam = '';
my $bam_file;
if ($bam_file) {
    die "bam file does not exists\n" unless -e $bam_file;
	$inbam = Bio::DB::Bam->open($bam_file,"r");
}
else {
	if (@ARGV) {
		if ($ARGV[0] eq '-') {
			die "STDIN is not supported\n";
		}elsif (-e $ARGV[0]) {
			$inbam = Bio::DB::Bam->open($bam_file,"r");		
		}else{
			die "$ARGV[0] does not exists\n";
		}
	}else{
		die "No input bam file\n";		
	}
}


my $inheader = $inbam->header;
my $target_names = $inheader->target_name;

my $outbam = Bio::DB::Bam->open($result_bam,"w");
$outbam->header_write($inheader);	

while (my $align = $inbam->read1) {
	my $read_name = $align->qname;
	$outbam->write1($align) if exists $reads{$read_name};
}

