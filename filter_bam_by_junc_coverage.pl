#!/usr/bin/perl -w

use strict;
use 5.010;

use Bio::DB::Sam;
use Getopt::Long;
use YAML;

my $usage = <<USAGE;
Usage: filter_bam_by_junc_coverage.pl [options...] -i aln.bam -p 3

SYSNOPSIS
filter_bam_by_junc_coverage.pl [options] -i aln.bam -p 3

 Options:
   -i --in-bam           bam file wich spliced reads
   -p --pos-count        count of aligment start postion, 
                         default on, value 3.
   -r --read-count       count of reads span one junction,
                         default off, value 0.
   -o --out-prefix       output prefix of the reslut
   -m --min-intron   

USAGE

#pare the arguments
my $bam_file   = '';
my $read_count = 0;
my $pos_count  = 3;
my $prefix     = 'out';
my $min_intron = 60;
die $usage
  unless GetOptions(
    "i|in-bam=s"     => \$bam_file,
    "p|pos-count:i"  => \$pos_count,
    "r|read-count:i" => \$read_count,
    "o|out-prefix:s" => \$prefix,
    "m|min-intorn:i" => \$min_intron,
  );

die "only one of pos-count or read-count can be activitated.\n" . $usage
  unless ( $read_count xor $pos_count );
my $result_bam = $prefix . "_filterd_junction.bam";

my %introns;

#open the in bam file
my $count_bam          = Bio::DB::Bam->open( $bam_file, "r" );
my $count_inheader     = $count_bam->header;
my $count_target_names = $count_inheader->target_name;

while ( my $align = $count_bam->read1 ) {

    my $chr   = $count_target_names->[ $align->tid ];
    my $start = $align->pos + 1;                        #1-based
    my $end   = $align->calend;
    die if $end < $start;
    my $loci     = "$chr:$start-$end";
    my $cigar_rf = $align->cigar_array;

    #calculate the referece region
    my @starts;
    my @n_loci;
    my $start_t = $start;
    push @starts, $start_t;
    for my $i ( 0 .. $#{$cigar_rf} ) {
        if ( $cigar_rf->[$i][0] =~ /M|N|D/i ) {
            $start_t += $cigar_rf->[$i][1];
            push @starts, $start_t;
            push @n_loci, $i
              if ( $cigar_rf->[$i][0] eq 'N' )
              && ( $cigar_rf->[$i][1] >= $min_intron );

        }
        elsif ( $cigar_rf->[$i][0] =~ /I/i ) {
            push @starts, $start_t;
        }
        else {
            warn "ecounted:$cigar_rf->[$i][0]\n";
        }
    }

    foreach my $i (@n_loci) {
        my $start_t = $starts[$i];               #1-based;
        my $end_t   = $starts[ $i + 1 ] - 1;
        my $loci_t  = "$chr:$start_t-$end_t";    #1 - based;

        if ( exists $introns{$loci_t} ) {
            $introns{$loci_t}{read_count}++;
            if ( exists $introns{$loci_t}{details}{$start} ) {
                $introns{$loci_t}{details}{$start}++;
            }
            else {
                $introns{$loci_t}{details}{$start} = 1;
                $introns{$loci_t}{pos_count}++;
            }

        }
        else {
            $introns{$loci_t}{read_count}      = 1;
            $introns{$loci_t}{pos_count}       = 1;
            $introns{$loci_t}{details}{$start} = 1;
        }
    }
}

#check the results
foreach my $loci_t ( keys %introns ) {
    my $pos_count = keys %{ $introns{$loci_t}{details} };
    warn "error in constructing intron database: pos_count is not right\n"
      if $introns{$loci_t}{pos_count} != $pos_count;
    my $sum = 0;
    foreach my $pos ( keys %{ $introns{$loci_t}{details} } ) {
        $sum += $introns{$loci_t}{details}{$pos};
    }
    warn "error in constructing intron database:read_count is not right\n"
      if $introns{$loci_t}{read_count} != $sum;
}

#open the in bam file
my $inbam        = Bio::DB::Bam->open( $bam_file, "r" );
my $inheader     = $inbam->header;
my $target_names = $inheader->target_name;

my $outbam = Bio::DB::Bam->open( $result_bam, "w" );
$outbam->header_write($inheader);

my $t_count = 0;
my $d_count = 0;
my $l_count = 0;
my %discard;
while ( my $align = $inbam->read1 ) {

    my $chr   = $target_names->[ $align->tid ];
    my $start = $align->pos + 1;                  # 1-based
    my $end   = $align->calend;
    die if $end < $start;
    my $loci     = "$chr:$start-$end";
    my $cigar_rf = $align->cigar_array;

    my $debug = 0;

    #calculate the referece region
    my @starts;
    my @n_loci;
    my $start_t = $start;
    push @starts, $start_t;
    for my $i ( 0 .. $#{$cigar_rf} ) {
        if ( $cigar_rf->[$i][0] =~ /M|N|D/i ) {
            $start_t += $cigar_rf->[$i][1];
            push @starts, $start_t;
            push @n_loci, $i
              if ( $cigar_rf->[$i][0] eq 'N' )
              && ( $cigar_rf->[$i][1] >= $min_intron );

            # $debug = 1 if  $cigar_rf->[$i][0] eq 'N'
            # && $cigar_rf->[$i][1] < $min_intron;
        }
        elsif ( $cigar_rf->[$i][0] =~ /I/i ) {
            push @starts, $start_t;

            # $debug = 1;
        }
        else {
            warn "ecounted:$cigar_rf->[$i][0]\n";
        }
    }

    # $debug = 1 if @n_loci > 1;
    if ($debug) {
        say $loci;
        say $align->cigar_str;
        say Dump $cigar_rf;
        say Dump \@n_loci;
        say Dump \@starts;
    }

    #fetch the intron start and end
    my $hit_flag = 0;

    my @junctions = ();
    foreach my $i (@n_loci) {
        my $start_t = $starts[$i];               #1-based;
        my $end_t   = $starts[ $i + 1 ] - 1;     #1-based;
        my $loci_t  = "$chr:$start_t-$end_t";    #1 - based;
        push @junctions, $loci_t;
        say $loci_t if $debug;
        if ( exists $introns{$loci_t} ) {
            if ($pos_count) {
                $hit_flag = 1 if $introns{$loci_t}{pos_count} >= $pos_count;
            }
            elsif ($read_count) {
                $hit_flag = 1 if $introns{$loci_t}{read_count} >= $read_count;
            }
        }
        else {

        }

    }
    if ($hit_flag) {
        $outbam->write1($align);
        $l_count++;
    }
    else {
        $d_count++;
        foreach my $junction (@junctions) {
            if ( exists $discard{$junction} ) {
                $discard{$junction}++;
            }
            else {
                $discard{$junction} = 1;
            }
        }
    }
    $t_count++;
}
my $discard_intron_num = keys %discard;
my $total_intron_num   = keys %introns;
print "Total records: $t_count\n";
print "Discard records: $d_count\n";
print "Total introns: $total_intron_num\n";
print "Discard introns: $discard_intron_num\n";
print "Left records: $l_count\n";
