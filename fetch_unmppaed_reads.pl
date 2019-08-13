#!/usr/bin/perl -w
use strict;
use 5.010;
# use Data::Printer;
use File::Basename;
use IPC::Cmd qw[can_run run];
use Getopt::Long;

my $usage = <<USAGE;
SYSNOPSIS
tophat2_batch.pl [options] bw2db file|glob

 Options:
   -m|--mode          'pe' or 'se' end reads.
   -s|--summarize     summarize only, default off
   -a|--args          argments used for tophat
   -g|--gtf           gtf file
   -o|--output        output folder for removed reads 
   -p|--progress      progress, default 8
USAGE
my $mode       = 'pe';
my $gtf_file = '';
my $out_folder = dirname './';
my $summarize = 0;
my $progress   = 8;
my $args_file = '';
die $usage
  unless GetOptions(
    "m|mode:s"       => \$mode,
    "a|args:s"       => \$args_file,
    "g|gtf:s"       => \$gtf_file,
    "s|summarize"       => \$summarize,
    "o|output:s"   => \$out_folder,
    "p|progress=i" => \$progress,
  );
my $bowtie2db = shift;
die "Does not specifed bowtie2 DB\n" unless $bowtie2db;
-e "$bowtie2db.1.bt2"
  || -e "$bowtie2db.1.bt2l"
  || die "Cannot detect bowtie2 db:$bowtie2db\n";
my $glob = shift;
my @fqs  = glob "$glob";
die $usage unless @fqs;
can_run('bowtie2') or die 'bowtie2 is not installed!';
can_run('tophat') or die 'tophat is not installed!';
if ($gtf_file) {
	die "$gtf_file does not exits\n" unless -e $gtf_file;
}
$out_folder = dirname $out_folder;
mkdir $out_folder unless -d $out_folder;

my %fastqs = ();
my @suffx  = qw (_1.fq.gz _2.fq.gz _1.fastq.gz _2.fastq.gz _1.fq _2.fq _1.fastq _2.fastq
  .1.fq.gz .2.fq.gz .1.fastq.gz .2.fastq.gz .1.fq .2.fq .1.fastq .2.fastq
  .fq.gz .fq .fastq.gz .fastq);
if ( $mode eq 'pe' ) {
    foreach my $fq ( sort @fqs ) {
        if ( $fq =~ /(.*)[\.|_]([1|2])\.(fq\.gz|fq)/ ) {
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
        if ( $fq =~ /(.*)\.(fq\.gz|fq)/ ) {
            my $basename = basename( $fq, @suffx );
            $fastqs{$basename} = $fq;
        }
        else {
            warn "Not supported format:$fq\n";
        }
    }
}

my $tophat_param = "-p $progress ";
$tophat_param .= "-G $gtf_file " if $gtf_file;
if ($args_file) {
	open my $args_fh, "<", $args_file or die "cannot open $args_file\n";
	while (<$args_fh>) {
		chomp;
		$tophat_param .= "$_ ";		
	}
}
my $fq_in_param          = "";
print "Sample\tLeft\tLeft_Multi\tLeft_Uniq\tRight\tRight_Multi\tRight_Uniq\tPairs_Multi\tPairs_Uniq\n";
foreach my $basename ( sort keys %fastqs ) {
    if ( $mode eq "pe" ) {
        $fq_in_param  = " $fastqs{$basename}{1}  $fastqs{$basename}{2}";
    }
    else {
        $fq_in_param  = " $fastqs{$basename}";
    }
    my $stdout = "$basename.stdout.log";
    my $stderr = "$basename.stderr.log";
    my $command =
"tophat -o $out_folder/$basename $tophat_param $bowtie2db $fq_in_param >$stdout 2> $stderr";
	if ($summarize) {
		
	}else{
		warn "\t$command\n";
		my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
		  run( command => $command, verbose => 0 );
		if ($success) {
			warn "\tDone!\n";
		}
		else {
			warn "Something went wrong:\n$stderr_buf";
		}		
	}

    my %counts = summarize_tophat2( "$out_folder/$basename/align_summary.txt", $mode );
	# p %counts;
    print
"$basename\t$counts{left}{input}\t$counts{left}{multiple}\t$counts{left}{unique}\t$counts{right}{input}\t$counts{right}{multiple}\t$counts{right}{unique}\t$counts{right}{multiple}\t$counts{right}{unique}\n";
}

sub summarize_tophat2 {
    my $tophat2_file = shift;
    my $mode        = shift;
    $mode = 'se' if $mode ne 'pe';
    open my $fh, "<", $tophat2_file or die "cannot open $tophat2_file\n$!";
    my %counts;
    while (<$fh>) {
        chomp;
        if ( $mode eq 'pe' ) {
			#left
			$counts{left}{input} = $1 if $. == 2  && /Input\s+:\s+(\d+)/;
			$counts{left}{mapped} = $1 if $. == 3 && /Mapped\s+:\s+(\d+)\s+.*of input\)/;
			$counts{left}{multiple} = $1 if $. == 4 && /of these:\s+(\d+)\s+.*have multiple alignments/;

			#right
			$counts{right}{input} = $1 if $. == 6  && /Input\s+:\s+(\d+)/;
			$counts{right}{mapped} = $1 if $. == 7 && /Mapped\s+:\s+(\d+)\s+.*of input\)/;
			$counts{right}{multiple} = $1 if $. == 8 && /of these:\s+(\d+)\s+.*have multiple alignments/;
			
			#pairs
			$counts{pair}{mapped} = $1 if $. == 11  && /Aligned pairs:\s+(\d+)/;
			$counts{pair}{multiple} = $1 if $. == 12 && /of these:\s+(\d+)\s+.*have multiple alignments/;
        }
        else {

        }
    }
	$counts{left}{unique} = $counts{left}{mapped} - $counts{left}{multiple};
	$counts{right}{unique} = $counts{right}{mapped} - $counts{right}{multiple};
	$counts{pair}{unique} = $counts{pair}{mapped} - $counts{pair}{multiple};

    return %counts;
}

