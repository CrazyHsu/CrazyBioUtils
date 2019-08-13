#!/usr/bin/perl -w

use strict;
use 5.010;
use Bio::DB::Sam;
use autodie;
use YAML;
use Getopt::Long;

#usage information
my $usage = <<USAGE;
Usage: chrom_mapped_reads.pl [options...] bamfile 
This program is used to calculate the number of reads mapped to each chromosome

 Options:
  -o --out-prefix    Prefix of output pdf file, optional
  Note:
  The default parameters are optimized for arabidopsis.
USAGE
my $bam_file = shift or die $usage;
my $prefix = 'out';
die $usage
  unless GetOptions( "o|out-prefix:s" => \$prefix, );

open my $data_fh, ">", "$prefix.dat";
print $data_fh "chrom\tchrom_len\tchrom_reads\n";

my $bam    = Bio::DB::Bam->open($bam_file);
my $header = $bam->header;

my $target_names_rf = $header->target_name;
my %chrs_reads = map { $_ => 0 } @$target_names_rf;

my $chr_length_rf = $header->target_len;
my %chr_length;
warn "======================================\n";
warn "Calculating chromosome length...\n";
for my $i ( 0 .. @$target_names_rf - 1 ) {
    $chr_length{ $target_names_rf->[$i] } = $chr_length_rf->[$i];
}
warn "Calculating chromosome reads...\n";
while ( my $align = $bam->read1 ) {
    my $seqid   = $target_names_rf->[ $align->tid ];
    my $hit_num = $align->aux_get("NH");
    $chrs_reads{$seqid}++ if $hit_num == 1;
}
warn "Done!\n";
warn "======================================\n";
foreach my $chr (@$target_names_rf) {
    print $data_fh "$chr\t$chr_length{$chr}\t$chrs_reads{$chr}\n";
    print "$chr\t$chr_length{$chr}\t$chrs_reads{$chr}\n";
}
