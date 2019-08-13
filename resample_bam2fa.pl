#!/usr/bin/perl -w

use strict;
use 5.010;
use Bio::DB::Sam;
use Algorithm::Numerical::Sample qw(sample);
use Getopt::Long;

my $usage = <<USAGE;
resample_bam2fa V1.1, written by corephi
This program is used to random pickup specified number of reads
from a bam, and store to fasta file
----------------------------------------------------------
More scripts? Join "bioinformatics*CN" QQ group: 276151571
If you have any probel or suggestions about this program, 
please mail to: hpyu\@genetics.ac.cn
----------------------------------------------------------
Usage: resample_bam2fa [options...] bamfile


 Options:
  -n|--number             number of reads to ramdom sample
                          optional, default "20"
  -o|--out-prefix		  prefix of output files.
                          optional, default "output"

USAGE
my $prefix = "output";
my $num    = 20;
die $usage
  unless GetOptions(
    "n|number:s"     => \$num,
    "o|out-prefix:s" => \$prefix
  );
my $bamfile = shift or die $usage;
my $fafile = "${prefix}_$num.fa";
SampleBam( $bamfile, $fafile, $num );

sub SampleBam {
    my ( $bamfile, $fafile, $num ) = @_;

    #get the number of alignments;
    my $sam = Bio::DB::Sam->new( -bam => $bamfile );
    my $total_num;
    my $iterator = $sam->features(
        -type     => "match",
        -iterator => 1
    );
    $total_num++ while $iterator->next_seq;

    #sample;
    my @sample_array;
    if ( $num >= $total_num ) {
        @sample_array = 1 .. $total_num;
    }
    else {
        @sample_array = sample(
            -set         => [ 1 .. $total_num ],
            -sample_size => $num
        );
    }
    my %sample = map { $_ => 1 } @sample_array;

    #open file to read;
    open my $out_fh, ">", $fafile or die "Cannot Create file $fafile";

    #open bam to output sample
    my $bam    = Bio::DB::Bam->open($bamfile);
    my $header = $bam->header;
    my $line_count;
    while ( my $aln = $bam->read1() ) {
        $line_count++;
        my $dna = $aln->qseq;
        my $id  = $aln->qname;
        print $out_fh ">$id\n$dna\n" if exists $sample{$line_count};
    }

    close $out_fh;
}

