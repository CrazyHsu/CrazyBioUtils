#!/usr/bin/perl -w
use strict;
use 5.010;
use Getopt::Long;
use Data::Printer;
use File::Basename;
use IPC::Cmd qw[can_run run];
use File::Copy;

my $usage = <<USAGE;
SYSNOPSIS
remove_low_quality_reads.pl file|glob

 Options:
   -m|--mode          'pe' or 'se' end reads.
   -a|--adaptor-lib   default "N", see NGSQCToolkit
   -l|--cut-length    cutoff length, default 70 
   -s|--cut-score     cutoff score, default 20
   -z|--zip           Output zipped results, default off
   -o|--output        output folder for removed reads 
   -p|--progress      progress, default 8
USAGE
my $mode       = 'pe';
my $out_folder = dirname './';
my $progress   = 8;
my $cutoff_len = 70;
my $cutoff_score = 20;
my $zipped = 0;
my $primer_lib = "N";
die $usage
  unless GetOptions(
    "m|mode:s"       => \$mode,
    "o|output=s"   => \$out_folder,
	"a|adapter-lib:s"       =>  \$primer_lib,
	"l|cut-length:i"       =>  \$cutoff_len,
	"s|cut-score:i"       =>  \$cutoff_score,
	"z|zip"       =>  \$zipped,
    "p|progress=i" => \$progress,
  );
  
my @fqs  = ();
push @fqs, glob $_ foreach @ARGV;
die $usage unless @fqs;

can_run('IlluQC_PRLL.pl') or die 'NGSQCToolkit is not installed!';
$out_folder =~ s/[\/|\|]+$//;
$zipped = $zipped  ? 'g' : 't';

my %fastqs = ();
my @suffx  = qw (_1.fq.gz _2.fq.gz _1.fastq.gz _2.fastq.gz _1.fq _2.fq _1.fastq _2.fastq
  .1.fq.gz .2.fq.gz .1.fastq.gz .2.fastq.gz .1.fq .2.fq .1.fastq .2.fastq
  .fq.gz .fq .fastq.gz .fastq);
if ( $mode eq 'pe' ) {
    foreach my $fq ( sort @fqs ) {
        if ( $fq =~ /(.*)[\.|_]([1|2])\.(fq\.gz|fq|fastq\.gz|fastq)/ ) {
            my $basename = basename( $fq, @suffx );
            my $mate = $2;
            $fastqs{$basename}{$mate} = $fq;
        }
        else {
            warn "Not supported format:$fq\n";
        }
    }
}else {
    foreach my $fq ( sort @fqs ) {
        if ( $fq =~ /(.*)\.(fq\.gz|fq|fastq\.gz|fastq)/ ) {
            my $basename = basename( $fq, @suffx );
            $fastqs{$basename} = $fq;
        }
        else {
            warn "Not supported format:$fq\n";
        }
    }
}
p %fastqs;

mkdir $out_folder unless -e $out_folder;
my $log_foder = $out_folder."/logs";
mkdir $log_foder unless -e $log_foder;


my $fq_in_param          = "";
print "Sample\tTotal\tLQ\tContaminated\tLeft\tPct\n";
foreach my $basename ( sort keys %fastqs ) {
    warn "Removing $basename:\n";
	my $remove_foder = $log_foder."/$basename";
	mkdir $remove_foder unless -e $remove_foder;
	my $stat_file = '';
    if ( $mode eq "pe" ) {
        $fq_in_param  = "-pe $fastqs{$basename}{1} $fastqs{$basename}{2}";
		$stat_file = basename($fastqs{$basename}{1})."_".basename($fastqs{$basename}{2})."_stat";
    }
    else {
        $fq_in_param  = "-se $fastqs{$basename}";
		$stat_file = basename($fastqs{$basename})."_stat";
    }
	$fq_in_param .= " $primer_lib A";
	my $command = "IlluQC_PRLL.pl $fq_in_param -o $remove_foder -l $cutoff_len -s $cutoff_score -c $progress -z g -t 2 1> $log_foder/$basename.stdout.log 2> $log_foder/$basename.stderr.log";
    warn "\t$command\n";
    my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
      run( command => $command, verbose => 0 );
    if ($success) {
        warn "\tDone!\n";
    }
    else {
        warn "Something went wrong:\n$stderr_buf";
    }
	
	#counts
    my %counts = summarize( $remove_foder."/$stat_file", $mode );
    print
"$basename\t$counts{total}\t$counts{LQ}\t$counts{Contaminated}\t$counts{HQ}\t$counts{Pct}\n";

	#move files
    if ( $mode eq "pe" ) {
        my $sourcefile1  = "$remove_foder/".basename($fastqs{$basename}{1})."_filtered";
		$sourcefile1 .= ".gz" if $zipped;
		my $destifile1 = "$out_folder/${basename}_1.fq";
		$destifile1 .= ".gz" if $zipped;
        my $sourcefile2  = "$remove_foder/".basename($fastqs{$basename}{2})."_filtered";
		$sourcefile2 .= ".gz" if $zipped;
		my $destifile2 = "$out_folder/${basename}_2.fq";
		$destifile2 .= ".gz" if $zipped;
		move($sourcefile1, $destifile1);
		move($sourcefile2, $destifile2);
    }
    else {
        my $sourcefile  = "$remove_foder/".basename($fastqs{$basename})."_filtered";
		$sourcefile .= ".gz" if $zipped;
		move($sourcefile, $out_folder);
    }	

}

sub summarize {
    my $bowite_file = shift;
    my $mode        = shift;
    $mode = 'se' if $mode ne 'pe';
    open my $fh, "<", $bowite_file or die "cannot open $bowite_file\n$!";
    my %counts;
    while (<$fh>) {
        chomp;
        if ( $mode eq 'pe' ) {

        }
        else {

        }
		my @lines = split /\s+/;
		given ( $. ) {
			when (13) {
				$counts{total} = $lines[-1];
			}
			when (14){
				$counts{LQ} = $counts{total} -$lines[-1];
			}
			when (20){
				$counts{Contaminated} = $lines[-1];
			}			
			when (21){
				$counts{HQ} = $lines[-1];
			}
			when (22){
				$counts{Pct} = $lines[-1];
			}
		}
		
    }
    return %counts;
}

