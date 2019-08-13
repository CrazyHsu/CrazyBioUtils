#!/usr/bin/perl -w
use strict;
use 5.010;
# use Data::Printer;
use File::Basename;
use IPC::Cmd qw[can_run run];
use Getopt::Long;

my $usage = <<USAGE;

remove_reads_by_bowtie2db V1.1, written by corephi
This program is used to random pickup specified number of reads
from a bam, and store to fasta file
----------------------------------------------------------
More scripts? Join "bioinformatics*CN" QQ group: 276151571
If you have any probel or suggestions about this program, 
please mail to: hpyu\@genetics.ac.cn
----------------------------------------------------------
SYSNOPSIS:

remove_reads_by_bowtie2db.pl [options] bw2db file|glob

 Options:
   -m|--mode          'pe' or 'se' end reads.
   -o|--output        output folder for removed reads 
   -p|--progress      progress, default 8

Note: 

Suppored formats:
_1.fq.gz _2.fq.gz _1.fastq.gz _2.fastq.gz _1.fq _2.fq _1.fastq _2.fastq
.1.fq.gz .2.fq.gz .1.fastq.gz .2.fastq.gz .1.fq .2.fq .1.fastq .2.fastq
.fq.gz .fq .fastq.gz .fastq
   
USAGE
my $mode       = 'pe';
my $out_folder = dirname './';
my $progress   = 8;
die $usage
  unless GetOptions(
    "m|mode:s"       => \$mode,
    "o|output:s"   => \$out_folder,
    "p|progress=i" => \$progress,
  );
my $xRNA_bowtiedb = shift;
if ($xRNA_bowtiedb) {
	-e "$xRNA_bowtiedb.1.bt2"
	  || -e "$xRNA_bowtiedb.1.bt2l"
	  || die "Cannot detect bowtie2 db:$xRNA_bowtiedb";	
}else{
	die "No bowtie2 db was specified\n$usage\n";
}


my @fqs;
foreach my $glob_str (@ARGV) {
	push @fqs, grep {-s $_} glob $glob_str;
}
die $usage unless @fqs;
can_run('bowtie2') or die 'bowtie2 is not installed!';
$out_folder =~ s/[\/|\|]+$//;
mkdir $out_folder unless -e $out_folder;

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
}
else {
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

my $bowtie2_dexRNA_param = "-k 1 -t -p $progress";
my $fq_in_param          = "";
my $fq_out_param         = "";
print "Sample\tTotal\tMapped\tLeft\n";
foreach my $basename ( sort keys %fastqs ) {
    warn "Removing $basename:\n";
    if ( $mode eq "pe" ) {
        $fq_in_param  = "-1 $fastqs{$basename}{1} -2 $fastqs{$basename}{2}";
        $fq_out_param = "--un-conc-gz $out_folder/${basename}_%.fq.gz";
    }
    else {
        $fq_in_param  = "-U $fastqs{$basename}";
        $fq_out_param = "--un-gz $out_folder/${basename}.fq.gz";
    }
    my $log_filename = "$out_folder/Remove_xRNA_$basename.log";
    my $command =
"bowtie2 $bowtie2_dexRNA_param -x $xRNA_bowtiedb $fq_in_param $fq_out_param >/dev/null 2> $log_filename";
    warn "\t$command\n";
    my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
      run( command => $command, verbose => 0 );
    if ($success) {
        warn "\tDone!\n";
    }
    else {
        warn "Something went wrong:\n$stderr_buf";
    }
    my %counts = summarize_bowtie2( $log_filename, $mode );
    print
"$basename\t$counts{concordant}{total}\t$counts{concordant}{unique}\t$counts{concordant}{unmapped}\n";
}

sub summarize_bowtie2 {
    my $bowite_file = shift;
    my $mode        = shift;
    $mode = 'se' if $mode ne 'pe';
    open my $fh, "<", $bowite_file or die "cannot open $bowite_file\n$!";
    my %counts;
    while (<$fh>) {
        chomp;
        if ( $mode eq 'pe' ) {
            $counts{concordant}{total} = $1 if /(\d+).*were paired; of these:$/;
            $counts{concordant}{unmapped} = $1
              if /(\d+).*aligned concordantly 0 times$/;
            $counts{concordant}{unique} = $1
              if /(\d+).*aligned concordantly exactly 1 time$/;
            $counts{concordant}{multiple} = $1
              if /(\d+).*aligned concordantly >1 times$/;
            $counts{discordantly}{mapped} = $1
              if /(\d+).*aligned discordantly 1 time$/;
        }
        else {
            $counts{concordant}{total} = $1 if /(\d+).*reads; of these:$/;
            $counts{concordant}{unmapped} = $1
              if /(\d+).*aligned 0 times$/;
            $counts{concordant}{unique} = $1
              if /(\d+).*aligned exactly 1 time$/;
            $counts{concordant}{multiple} = $1
              if /(\d+).*aligned >1 times$/;
            $counts{discordantly}{mapped} = 0;
        }
    }
    return %counts;
}

