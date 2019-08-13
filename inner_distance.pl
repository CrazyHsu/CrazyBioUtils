#!/usr/bin/perl -w
use strict;
use 5.010;
# use Data::Printer;
use File::Basename;
use IPC::Cmd qw[can_run run];
use Getopt::Long;

my $usage = <<USAGE;
SYSNOPSIS
inner_distance_batch.pl [options] tophat_batch_folder

 Options:
   -g|--bed           bed file
   -l|--lower         lower bound size, default -250
   -u|--upper         upper bound size, default 250
   -s|--step          step size
   -o|--output        output folder for removed reads 
USAGE
my $bed_file = '';
my $lower = -250;
my $upper = 250;
my $step = 10;
my $out_folder = dirname './';
die $usage
  unless GetOptions(
    "g|gtf:s"       => \$bed_file,
    "o|output:s"   => \$out_folder,
    "l|lower:i"   => \$lower,
    "s|step:i"   => \$step,
    "u|upper:i"   => \$upper,
  );
can_run('inner_distance.py') or die 'RseQC is not installed!';
my $tophat_folder = shift;
$tophat_folder =~ s/[\/|\|]+$//;
die "$tophat_folder does not exits\n;" 
	unless -e $tophat_folder && -d $tophat_folder;
if ($bed_file) {
	die "$bed_file does not exits\n" unless -e $bed_file;
}
$out_folder =~ s/[\/|\|]+$//;
mkdir $out_folder unless -d $out_folder;
my $log_folder = "$out_folder/logs";
mkdir $log_folder unless -d $log_folder;


opendir(my $dh, $tophat_folder) || die "can't opendir $tophat_folder: $!";
# my @dots = grep { /^\./ && -f "$tophat_folder/$_" } readdir($dh);
my @samples = grep { !/^\.+/ && -d "$tophat_folder/$_" } readdir($dh);
closedir $dh;

foreach my $sample (@samples) {
	my $bam_file;
	if ( $sample =~ /(.*)_R(\d+)/ ) {
		my $cond = $1;
		my $rep = $2;
		$bam_file = "$tophat_folder/$sample/accepted_hits.bam";
	}
	else {
		warn "Directory $sample seem not like sample data, skipping...\n";
	}
	warn "Staring calculating inner_distance of $sample...\n";


	my $command =
	"inner_distance.py -o $out_folder/$sample -l $lower -u $upper -s $step -r $bed_file -i $bam_file >$log_folder/$sample.inner_distance.stdout.log 2> $log_folder/$sample.inner_distance.stderr.log";

	warn "\t$command\n";
	my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
	run( command => $command, verbose => 0 );
	if ($success) {
	    warn "\tDone!\n";
	}
	else {
		warn "Something went wrong:\n@{$stderr_buf}";
	}	
}

