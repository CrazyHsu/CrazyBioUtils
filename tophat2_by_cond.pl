#!/usr/bin/perl -w
use strict;
use 5.010;
# use Data::Printer;
use File::Basename;
use IPC::Cmd qw[can_run run];
use Getopt::Long;

my $usage = <<USAGE;
SYSNOPSIS
tophat2.pl [options] bw2db file|glob

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
my @fqs;
foreach my $glob_str (@ARGV) {
	push @fqs, grep {-s $_} glob $glob_str;
}
die $usage unless @fqs;
can_run('bowtie2') or die 'bowtie2 is not installed!';
can_run('tophat') or die 'tophat is not installed!';
if ($gtf_file) {
	die "$gtf_file does not exits\n" unless -e $gtf_file;
}
$out_folder =~ s/[\/|\|]+$//;
mkdir $out_folder unless -d $out_folder;
my $log_folder = $out_folder."/logs";
mkdir $log_folder unless -d $log_folder;


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
my $header = $mode eq 'pe' ? 
	"Sample\tLeftTotal\tLeftMultiple\tLeftUnique\tRightTotal\tRightMultiple\tRightUnique\tPairsMultiple\tPairsUnique\tPairsConcordant\tPairsDiscordant\tParisOveralRate\tConcordantRate\n" :
	"Sample\tInput\tUnique\tMultiple\tOveralRate\n";
print $header ;
# p %fastqs;

my %groups;
foreach my $basename ( sort keys %fastqs ) {
	my ($group, $replicate);
	if ( $basename =~ /^(.+)_R(\d+)/ ) {
		$group = $1;
		$replicate = $2;
	}
	else {
		warn "Not supported group format:$basename\nTeat as samples without replicates\n";
		$group = $basename;
		$replicate = 1;
	}
	
    if ( $mode eq "pe" ) {
		$groups{$group}{$replicate}{1} = $fastqs{$basename}{1};
		$groups{$group}{$replicate}{2} = $fastqs{$basename}{2};
    }else {
		$groups{$group}{$replicate} = fastqs{$basename};
    }	
}







warn "Staring Tophat2...\n";
foreach my $group ( sort keys %groups ) {
	my $fq_in_param;
    if ( $mode eq "pe" ) {
		my @mate1 = ();
		my @mate2 = ();
		foreach my $replicate (sort keys %{$groups{$group}}) {
			push @mate1, $groups{$group}{$replicate}{1};
			push @mate2, $groups{$group}{$replicate}{2};
		}
		my $mate1_files = join ",", @mate1;
		my $mate2_files = join ",", @mate2;
        $fq_in_param  = " $mate1_files $mate2_files ";
    }
    else {
		my @fastqs;
		foreach my $replicate (sort keys %{$groups{$group}}) {
			push @fastqs, $groups{$group}{$replicate};
		}
        $fq_in_param  = " ". join ",", @fastqs;
    }
    my $stdout = "$group.stdout.log";
    my $stderr = "$group.stderr.log";
    my $command =
"tophat -o $out_folder/$group $tophat_param $bowtie2db $fq_in_param >$log_folder/$stdout 2> $log_folder/$stderr";
	if ($summarize) {
		warn "\t$command\n";		
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

    my %counts = summarize_tophat2( "$out_folder/$group/align_summary.txt", $mode );
	# p %counts;
	my $data = $mode eq 'pe' ? 
		"$group\t$counts{left}{input}\t$counts{left}{multiple}\t$counts{left}{unique}\t$counts{right}{input}\t$counts{right}{multiple}\t$counts{right}{unique}\t$counts{pair}{multiple}\t$counts{pair}{unique}\t$counts{pair}{concordant}\t$counts{pair}{discordant}\t$counts{rate}{overal}\t$counts{rate}{concordant}\n" :
		"$group\t$counts{pair}{input}\t$counts{pair}{multiple}\t$counts{pair}{unique}\t$counts{rate}{concordant}\n";
    print $data;
;
}

sub summarize_tophat2 {
    my $tophat2_file = shift;
    my $mode        = shift;
    $mode = 'se' if $mode ne 'pe';
    open my $fh, "<", $tophat2_file or die "cannot open $tophat2_file\n$!";
    my %counts;
    if ( $mode eq 'pe' ) {
		while (<$fh>) {
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
			$counts{pair}{discordant} = $1 if $. == 13 && /\s+(\d+)\s+.*are discordant alignments/;
		

			#mapping rate
			$counts{rate}{overal} = $1 if $. == 9 &&  /(.*)% overall read mapping rate/;
			$counts{rate}{concordant} = $1 if $. == 14 &&  /(.*)% concordant pair alignment rate/;
		}
		$counts{pair}{unique} = $counts{pair}{mapped} - $counts{pair}{multiple};
		$counts{pair}{concordant} = $counts{pair}{mapped}  - $counts{pair}{discordant};
		$counts{left}{unique} = $counts{left}{mapped} - $counts{left}{multiple};
		$counts{right}{unique} = $counts{right}{mapped} - $counts{right}{multiple};	

	}else{
		while (<$fh>) {
			#pairs
			$counts{pair}{input} = $1 if $. == 2  && /Input:\s+(\d+)/;
			$counts{pair}{mapped} = $1 if $. == 3  && /Mapped:\s+(\d+)/;
			$counts{pair}{multiple} = $1 if $. == 4 && /of these:\s+(\d+)\s+.*have multiple alignments/;	
			$counts{rate}{overal} = $1 if $. == 5 &&  /$(.*)% overall read mapping rate/;
		}	
		$counts{pair}{unique} = $counts{pair}{mapped} - $counts{pair}{multiple};

	}
	close $fh;
	
    return %counts;
}

