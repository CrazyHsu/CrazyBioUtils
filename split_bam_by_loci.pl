#!/usr/bin/perl -w

use strict;
use 5.010;

use Bio::DB::Sam;
use Getopt::Long;

my $usage = <<USAGE;
Usage: split_bam_by_loci.pl [options...] bamfile
written by corephi, group:276151571
this program is used to split bamfile by multi hits, splice junction etc.

 Options:
  -l|--loci-file        the loci u want to keep
  -o|--out-prefix       prefix of output files. optional, default "output".

USAGE

my $loci_file  = '';
my $out_prefix = 'out';

die $usage
  unless GetOptions(
    "l|loci-file:s"  => \$loci_file,
    "o|out-prefix:s" => \$out_prefix,
  );

my $inbam_file = shift or die $usage;
die "loci-file:$loci_file dos not exists\n" unless -e $loci_file;
die "bam file:$inbam_file dos not exists\n" unless -e $inbam_file;

open my $loci_fh, "<", $loci_file
  or die "cannot open file:$loci_file\n";
my $insam = Bio::DB::Sam->new(
    -bam       => $inbam_file,
    -autoindex => 1,
);
my %kept_reads;
while (<$loci_fh>) {
    chomp;
    my ( $chr, $start, $end, $attribute ) = split /\t/;
    my @alignments = $insam->get_features_by_location(
        -seq_id => $chr,
        -start  => $start,
        -end    => $end
    );
    for my $align (@alignments) {
        my $read_name = $align->query->name;
        my $aln_num   = $align->aux_get("NH");
        $kept_reads{$read_name} = $aln_num;
    }
}
$insam = 0;

#open the in bam file
my $inbam = Bio::DB::Bam->open( $inbam_file, "r" );
die $usage unless $inbam;
my $inheader = $inbam->header;
my $out_left = Bio::DB::Bam->open( $out_prefix . "_left.bam", "w" );
$out_left->header_write($inheader);
my $out_discard = Bio::DB::Bam->open( $out_prefix . "_discard.bam", "w" );
$out_discard->header_write($inheader);

my $total_record_num   = 0;
my $kept_record_num    = 0;
my $discard_record_num = 0;

while ( my $align = $inbam->read1 ) {
    $total_record_num++;
    my $read_name = $align->qname;
    if ( exists $kept_reads{$read_name} ) {
        $kept_record_num++;
        $out_left->write1($align);
    }
    else {
        $discard_record_num++;
        $out_discard->write1($align);
    }
}

print "Total records: $total_record_num\n";
print "Left records: $kept_record_num\n";
print "Discard records: $discard_record_num\n";
