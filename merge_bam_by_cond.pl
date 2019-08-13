#!/usr/bin/perl -w
use strict;
use 5.010;
use Data::Printer;
use File::Basename;
use IPC::Cmd qw[can_run run];
use Cwd 'abs_path';
use Getopt::Long;
use threads;
use Thread::Queue;

my $usage = <<USAGE;
SYSNOPSIS
merge_reads_by_cond.pl V1.0, written by corephi

This program is used to run merge the bams files in sample
to level to condtion level


----------------------------------------------------------
More scripts? Join "bioinformatics*CN" QQ group: 276151571
If you have any probel or suggestions about this program, 
please mail to: hpyu\@genetics.ac.cn
----------------------------------------------------------
Usage: merge_reads_by_cond.pl [options]  file|glob

 Options:
   -b|--bam-flag           default NA, can be setted as'accepted_hits'
   -o|--output             output folder for removed reads 
   -p|--progress           progress, default 8
   
Notes:
This program will automatic guess the sample name for each bam files with
the format of "Material_Treatment_ReplicatNumber", such as "Leaf_Day3_R1",
"Leaf_Day5_R2", "mutation_drought_R1".
To compatible with tophat, this grogram can also guess the sample name
from the directory. For example file "../Leaf_Day3_R1/accepted_hits.bam", 
will be named as "Leaf_Day3_R1"

Samples:

merge_reads_by_cond.pl -o counts -r genes.gtf ../uniquelyAligned/*.bam
merge_reads_by_cond.pl -o counts -b accepted_hits -r genes.gtf ../*.bam


USAGE
my $out_folder = dirname './';
my $progress   = 16;
my $length_flag = 0;
my $bam_flag = '';
my $help = 0;
die $usage
  unless GetOptions(
    "b|bam-flag:s"   => \$bam_flag,
    "o|output:s"   => \$out_folder,
    "p|progress=i" => \$progress,
    "h|help" => \$help,
  );
  
die $usage if $help;
#############################################################################
#Parameters
############################################################################# 

#check running enviroment
can_run('samtools') or die 'samtools is not installed!';

#create the output folder
$out_folder =~ s/[\/|\|]+$//;
mkdir $out_folder unless -d $out_folder;

my $log_folder = $out_folder."/logs";
mkdir $log_folder unless -d $log_folder;

#checking the bam files
my @bams = ();
foreach my $glob_str(@ARGV) {
	#pick up directory
	my @bam_tmp = ();
	if ($bam_flag) {
		@bam_tmp = map {abs_path("$_/$bam_flag.bam")} 
				grep {-d $_ and -s "$_/$bam_flag.bam"} 
				glob $glob_str;			
	}else{
		@bam_tmp = map {abs_path($_)} 
			grep {-s $_} 
				glob $glob_str;	
	}
	
	#pick up bams	
	push @bams, @bam_tmp if @bam_tmp;
}

die "No bams were found, please the the glob string\n" 
	unless @bams;

#############################################################################
#Getting the experimental design
#This is sutiable for differentical test
#Material_Treatment_ReplicatNumber
############################################################################# 
my %groups;
foreach my $bam ( sort @bams ) {
	my @paths = split /\//, $bam;
	
	my $basename = $bam;
	if ($bam_flag && $bam =~ /$bam_flag.bam/) {
		$basename = $paths[-2];
	}else{
		$basename = basename($bam, ".uniq.bam", "_uniq.bam", ".bam");
		warn "$bam does not match bam flag, sample name set as:$basename \n";
	}
	
	if ( $basename =~ /^(.+)_R(\d+)/ ) {
		my $group = $1;
		my $replicate = $2;
		$groups{$group}{$replicate} = $bam;
	}
	else {
		warn "Not supported group format:$basename\nTeat as samples without replicates\n";
		my $group = $basename;
		my $replicate = 1;
		$groups{$group}{$replicate} = $bam;
	}
}


#############################################################################
#run samtools multiple-progress
#############################################################################
my @commonds = (); 
foreach my $group (sort keys %groups) {
	my @group_bams = ();
	foreach my $replicate(sort {$a <=>$b } keys %{$groups{$group}}) {
		push @group_bams, $groups{$group}{$replicate};
	}
	my $bams =  join " ", @group_bams;
	my $stdout = "$group.stdout.log";
    my $stderr = "$group.stderr.log";
    my $command =
	"samtools merge $out_folder/$group.bam  $bams >$log_folder/$stdout 2> $log_folder/$stderr";
	push @commonds, $command;
}

warn "Staring HTseq-count ...\n";
run_parallel($progress, @commonds);


#===  FUNCTION  ================================================================
#         NAME: run_parallel
#      PURPOSE: given commands, run them in multiple threads
#   PARAMETERS: $process_num: 
#				$missions: commonds
#      RETURNS: NA
#===============================================================================
sub run_parallel{
	my ($process_num, @missions) = @_;
	my $stream = Thread::Queue->new(@missions,undef);
	my $mission_num = scalar @missions;

	#assgn the task
	my @running = ();
	my @Threads;
	while (@Threads < @missions) {
	    @running = threads->list(threads::running);

	    if ( @running < $process_num) {
			my $command = $stream->dequeue();
	        my $thread = threads->new(\&sub_run,$command);
       		push (@Threads, $thread);
	        my $tid = $thread->tid;
    	}
    	@running = threads->list(threads::running);
    	foreach my $thr (@Threads) {
        	if ($thr->is_running()) {
            		my $tid = $thr->tid;
        	}
        	elsif ($thr->is_joinable()) {
            		my $tid = $thr->tid;
            		$thr->join;
        	}
    	}
 
    	@running = threads->list(threads::running);
	}

	#join the threads
	while (@running) {
    		foreach my $thr (@Threads) {
        		$thr->join if ($thr->is_joinable());
    		}
    		@running = threads->list(threads::running);
    		sleep(3);
	}
	return 0;
}


sub sub_run {
	my $command = shift;
	warn "\t$command\n";
	my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
	  run( command => $command, verbose => 0 );
	if ($success) {
		warn "\tDone: $command!\n";
	}
	else {
		my @stderrs = @$stderr_buf;
		warn "Something went wrong:\n@stderrs";
	}		
}

