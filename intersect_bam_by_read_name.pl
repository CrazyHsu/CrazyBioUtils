#!/usr/bin/perl -w

use strict;
use 5.010;

use Bio::DB::Sam;
use Getopt::Long;
# use Data::Printer;

my $usage = <<USAGE;
Usage: intersect_bam_by_read_name.pl [options...] file1.bam file2.bam
written by corephi, group:276151571

 Options:
  -o|--out-file       prefix of output files. optional, default "out.bam".

USAGE

#pare the arguments
my $out_file = 'output.bam';
die $usage
  unless GetOptions(
    "o|out-prefix:s"    => \$out_file,
  );
#check files
die "You must specify at least 2 files" if @ARGV < 2;
my @bams = @ARGV; 
foreach my $bam (@bams) {
	die "$bam does not exists\n" unless -e $bam
}

#proceed the bam file
my $first_bam_file = shift @bams;
warn "Reading $first_bam_file\n";
my $shared_reads_rf = get_reads_name($first_bam_file);
# p $shared_reads_rf;
my $temp_in_bam = $bams[0];
my @temp_out_bams = ();
foreach my $bam (@bams) {
	my $temp_out_bam = "$bam.tmp";
	write_intersect_bam($temp_in_bam, $shared_reads_rf, $temp_out_bam);
	$shared_reads_rf = get_reads_name($temp_out_bam);
	$temp_in_bam = $temp_out_bam;
	push @temp_out_bams, $temp_out_bam;
}

my $last_temp = pop @temp_out_bams;
rename $last_temp, $out_file;
unlink $_ foreach @temp_out_bams;

sub get_reads_name {
	my $inbam_file = shift;
	my %reads = ();
	my $inbam = Bio::DB::Bam->open( $inbam_file, "r" );
	die $usage unless $inbam;
	my $inheader     = $inbam->header;
	my $target_names = $inheader->target_name;
	while ( my $align = $inbam->read1 ) {
		my $seqname = $align->query->name;
		$reads{$seqname} = 0;
	}
	return \%reads;
}

sub write_intersect_bam {
	my ($inbam_file, $reads_rf, $outbam_file) = @_;
	my $inbam = Bio::DB::Bam->open( $inbam_file, "r" );
	warn "Reading $inbam_file\n";
	my $outbam = Bio::DB::Bam->open( $outbam_file, "w" );
	die $usage unless $inbam;
	my $inheader     = $inbam->header;
	$outbam->header_write($inheader);
	my $target_names = $inheader->target_name;
	while ( my $align = $inbam->read1 ) {
		my $seqname = $align->query->name;
		$outbam->write1($align) if exists $reads_rf->{$seqname};
	}
}
