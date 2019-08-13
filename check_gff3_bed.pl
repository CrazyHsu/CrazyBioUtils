#!/usr/bin/perl -w

use strict;
use 5.010;
use Getopt::Long;
use YAML;

my $usage = <<USAGE;
Usage: check_gff3_bed.pl [options...] -b annotian.bed -g annotaion.gff

This program is used to compare gff annotatioin and bed format annotation

 Options:
  -b --bed        bed format annotation file
  -g --gff        gff3 format annotaion file
  -o --out-prefix		  prefix of output files.
                      optional, default "output"

USAGE

#pare the arguments
my $prefix = "output";
my ( $bed_file, $gff_file );
die $usage
  unless GetOptions(
    "b|bed=s"        => \$bed_file,
    "o|out-prefix:s" => \$prefix,
    "g|gff=s"        => \$gff_file,
  );
my $debug = 1;

my %cds;
my %exons;
my %gff;
open my $gff_fh, "<", $gff_file or die "cannot open file $gff_file: $!";
while (<$gff_fh>) {
    chomp;
    my (
        $chr,   $source, $type,  $start, $end,
        $score, $strand, $frame, $attribute
    ) = split /\t/;

    #store attribute to hash.
    my @temp = split /;/, $attribute;
    my %attributes;
    foreach my $tmp (@temp) {

        # print "atribute:$tmp\t" if $debug;
        my ( $key, $value ) = split /=/, $tmp;

        # print "key:$key\tvalue:$value\n" if $debug;
        my @values = split /,/, $value;

        # print @values,"\n" if $debug;
        $attributes{$key} = \@values;
    }

    #
    if ( $type =~ /CDS/ ) {
        my $transcript = $attributes{Parent}->[0];

        # say $transcript if $debug;
        if ( exists $cds{$transcript} ) {
            push @{ $cds{$transcript} }, [ $start, $end ];
        }
        else {
            $cds{$transcript} = [ [ $start, $end ], ];
        }
    }

    if ( $type =~ /exon/ ) {
        my $transcript = $attributes{Parent}->[0];

        # say $transcript if $debug;
        $exons{$transcript}{"$chr:$start-$end"} = 1;
    }

}

#find the cds
foreach my $transcript ( keys %cds ) {
    my @CDS = @{ $cds{$transcript} };

    # print Dump \@CDS;
    my ( $start, $end ) = @{ $CDS[0] };
    for my $i ( 1 .. $#CDS ) {
        $start = $CDS[$i][0] < $start ? $CDS[$i][0] : $start;
        $end   = $CDS[$i][1] > $end   ? $CDS[$i][1] : $end;
    }
    $gff{$transcript}{CDS} = "$start-$end";

    # say "$transcript\t$start\t$end" if $debug;
}

#fetch the exon information
foreach my $transcript ( keys %exons ) {
    $gff{$transcript}{exon} = $exons{$transcript};
}

# print Dump \%gff;

open my $bed_fh, "<", $bed_file or die "cannot open file $bed_file: $!";
my %bed;
while (<$bed_fh>) {
    chomp;
    next if /^track/;
    my (
        $chrom,   $start,      $end,        $name,
        $score,   $strand,     $thickStart, $thickEnd,
        $itemRgb, $blockCount, $blockSizes, $blockStarts
    ) = split /\t/;

    #cds
    my $cds_start = $thickStart + 1;
    my $cds_end   = $thickEnd;
    $bed{$name}{CDS} = "$cds_start-$cds_end" unless $thickStart == $thickEnd;

    #exon
    my @blockSize  = split ',', $blockSizes;
    my @blockStart = split ',', $blockStarts;

    for my $i ( 0 .. $#blockSize ) {
        my $exon_start = $start + $blockStart[$i] + 1;
        my $exon_end   = $start + $blockStart[$i] + $blockSize[$i];
        $bed{$name}{exon}{"$chrom:$exon_start-$exon_end"} = 1;
    }
}

#check the transcript
my ( @gff_transcript_only, @bed_transcipt_only );
my ( @gff_exon_only,       @bed_exon_only );
my (%cds_diff);

#check gff only
foreach my $transcript ( keys %gff ) {

    #check the transcipt
    if ( exists $bed{$transcript} ) {

        #check the exon
        foreach my $exon ( keys %{ $gff{$transcript}{exon} } ) {
            if ( exists $bed{$transcript}{exon}{$exon} ) {

            }
            else {
                push @bed_exon_only, "$transcript\t$exon";
            }
        }

        #check the cds
        if ( exists $gff{$transcript}{CDS} && $bed{$transcript}{CDS} ) {
            if ( $gff{$transcript}{CDS} eq $bed{$transcript}{CDS} ) {

            }
            else {
                $cds_diff{$transcript}{gff} = $gff{$transcript}{CDS};
                $cds_diff{$transcript}{bed} = $bed{$transcript}{CDS};
            }
        }
        else {
            if ( exists $gff{$transcript}{CDS} ) {
                $cds_diff{$transcript}{gff} = $gff{$transcript}{CDS};
                $cds_diff{$transcript}{bed} = 'NA';
            }
            elsif ( exists $bed{$transcript}{CDS} ) {
                $cds_diff{$transcript}{bed} = $bed{$transcript}{CDS};
                $cds_diff{$transcript}{gff} = 'NA';
            }
            else {

            }
        }
    }
    else {
        push @gff_transcript_only, $transcript;
    }
}

#check bed only
foreach my $transcript ( keys %bed ) {

    #check the transcipt
    if ( exists $gff{$transcript} ) {

        #check the exon
        foreach my $exon ( keys %{ $bed{$transcript}{exon} } ) {
            if ( exists $gff{$transcript}{exon}{$exon} ) {

            }
            else {
                push @bed_exon_only, "$transcript\t$exon";
            }
        }

        #check the cds
        if ( exists $gff{$transcript}{CDS} && $bed{$transcript}{CDS} ) {
            if ( $gff{$transcript}{CDS} eq $bed{$transcript}{CDS} ) {

            }
            else {
                $cds_diff{$transcript}{gff} = $gff{$transcript}{CDS};
                $cds_diff{$transcript}{bed} = $bed{$transcript}{CDS};
            }
        }
        else {
            if ( exists $gff{$transcript}{CDS} ) {
                $cds_diff{$transcript}{gff} = $gff{$transcript}{CDS};
                $cds_diff{$transcript}{bed} = 'NA';
            }
            elsif ( exists $bed{$transcript}{CDS} ) {
                $cds_diff{$transcript}{bed} = $bed{$transcript}{CDS};
                $cds_diff{$transcript}{gff} = 'NA';
            }
            else {

            }
        }

    }
    else {
        push @bed_transcipt_only, $transcript;
    }
}

#bed output
print "bed transcript only:\n";
foreach my $transcript (@bed_transcipt_only) {
    print "$transcript\n";
}

print "bed exon only:\n";
foreach my $exon (@bed_exon_only) {
    print "$exon\n";
}

#gff output
print "gff transcript only:\n";
foreach my $transcript (@gff_transcript_only) {
    print "$transcript\n";
}

print "gff exon only:\n";
foreach my $exon (@gff_exon_only) {
    print "$exon\n";
}

#cds
print "CDS:\n";
print "transcript_id\tgff\tbed\n";
foreach my $transcript ( keys %cds_diff ) {
    print
      "$transcript\t$cds_diff{$transcript}{gff}\t$cds_diff{$transcript}{bed}\n";
}
