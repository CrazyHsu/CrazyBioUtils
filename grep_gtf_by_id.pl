#!/usr/bin/perl -w

use strict;
use 5.010;
use Getopt::Long;
# use Data::Printer;

#usage information
my $usage = <<USAGE;
grep_gtf_by_id.pl -i in.gtf -l tid.txt > out.gtf
written by corephi, group:276151571
----------------------------------------------------------
More scripts? Join "bioinformatics*CN" QQ group: 276151571
If you have any probel or suggestions about this program, 
please mail to: hpyu\@genetics.ac.cn
----------------------------------------------------------
Usage: grep_gtf_by_id.pl -i transcripts.gtf -l id.txt
 Options:
  -i    gtf file name, the first line must be id, default STDIN
  -o    out.gtf, default STDOUT
  -r    reverse match, output the different line
  -h 	gtf file has a header, default 0;
  -l    idlist file name, it must be one line per id
  -m    match by 'gid' , 'acc' or 'tid', default 'tid'
  

USAGE
my $gtf_file     = '';
my $id_list_file = '';
my $reverse_flag = 0;
my $out_gtf = '';

my $header_flag  = 0;
my $mode         = 'tid';
die $usage
  unless GetOptions(
    "i:s" => \$gtf_file,
    "l=s" => \$id_list_file,
    "h"   => \$header_flag,
    'r'   => \$reverse_flag,
    'm=s' => \$mode,
    "o:s" => \$out_gtf,

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

my $out_fh;
if ($$out_gtf && $out_gtf ne '-') {
    open $out_fh, ">", $out_gtf or die "cannot open file $out_gtf:$!\n";
}
else {
    $out_fh = *STDOUT;
}

#append content to a file
my $gtf_file_fh;
if ($gtf_file && $gtf_file ne '-') {
    die "gtf file does not exists\n" unless -e $gtf_file;
    open $gtf_file_fh, "<", $gtf_file or die "cannot open file $gtf_file:$!\n";
}
else {
	if (@ARGV) {
		if ($ARGV[0] eq '-') {
			$gtf_file_fh = *STDIN;
		}elsif (-e $ARGV[0]) {
			open $gtf_file_fh, "<", $ARGV[0] or die "cannot open file $ARGV[0]:$!\n";		
		}else{
			die "$ARGV[0] does not exists\n";
		}
	}else{
		$gtf_file_fh = *STDIN;				
	}
}
	
if ($header_flag) {
    my $header = readline $gtf_file_fh;
    print $out_fh $header if $header_flag;
}

while (<$gtf_file_fh>) {
    my $line = $_;
    next if /^#/;
    s/\r?\n//;
    my @lines     = split /\t/;
    my $attribute = pop @lines;
    my $type      = lc $lines[2];

    # say $attribute;
    my @temp = split ';', $attribute;
    my %attribttes;
    foreach my $attr (@temp) {
        $attr =~ s/^\s+//;
        next unless $attr;
        my ( $key, $value );
        if ( $attr =~ /(\S+) "(.*)"/ ) {
            ( $key, $value ) = ( $1, $2 );

            # say "$key\t\t\t$value";
            $attribttes{$key} = $value;
        }
        else {
            warn "$attr:is not well formatted\n";
        }
    }
    my $gid = $attribttes{gene_id};
    my $tid = $attribttes{transcript_id};
    my $acc = $attribttes{accession_numer};
	my $id;
	if ($mode eq 'gid') {
		$id = $gid;
	}elsif ($mode eq 'acc'){
		$id = $acc;
	}else{
		$id = $tid;
	}
	next unless $tid;
    if ($reverse_flag) {
        print $out_fh $line unless exists $lists{$id};
    }
    else {
        print $out_fh $line if exists $lists{$id};
    }

}

