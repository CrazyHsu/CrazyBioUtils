#!/usr/bin/perl -w
use strict;
use 5.010;
use Data::Printer;
use File::Basename;
use IPC::Cmd qw[can_run run];
use Cwd 'abs_path';
use Getopt::Long;

my $usage = <<USAGE;
SYSNOPSIS
cuffdiff.pl [options] gtffile file|glob

 Options:
   -b|--bias-correct       genome sequence fasta file, optional 
   -m|--mask-file          genome sequence fasta file, optional
   -u|--multi-read-correct multiple aligned read count correct, default on
   -n|--no-diff            no differentail test, default off
   -s|--bam-flag           default NA, can be setted as'accepted_hits'
   -c|–-min-align-count    minimum alignment count, default 10
   -q|--FDR                FDR, default 0.05
   -f|--min-isfm-fraction  multiple aligned read count correct, default 0.1
   -a|--args               argments used for cufflinks
   -o|--output             output folder for removed reads 
   -p|--progress           progress, default 8
USAGE
my $out_folder = dirname './';
my $progress   = 8;
my $no_diff = 0;
my $bias_correct = '';
my $bam_flag = '';
my $multi_read_correct = 0;
my $min_isoform_fraction = 0.05;
my $fdr = 0.05;
my $mask_file = '';
my $args_file = '';
die $usage
  unless GetOptions(
    "a|args:s"       => \$args_file,
    "n|no-diff"   => \$no_diff,
    "s|bam-flag:s"   => \$bam_flag,
    "b|bias-correct:s"   => \$bias_correct,
    "M|mask-file:s"   => \$mask_file,
    "u|multi-read-correct"   => \$multi_read_correct,
    "f|min-isoform-fraction:f"   => \$min_isoform_fraction,
    "q|FDR:f"   => \$fdr,
    "o|output:s"   => \$out_folder,
    "p|progress=i" => \$progress,
  );

my $gtf_file = shift @ARGV;
# file_check($gtf_file) if $gtf_file;
# file_check($bias_correct) if $bias_correct;
# file_check($mask_file) if $mask_file;

can_run('cuffdiff') or die 'cuffdiff is not installed!';
$out_folder =~ s/[\/|\|]+$//;
mkdir $out_folder unless -d $out_folder;

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
#constructing group and replicates
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

my $cuffdiff_param = "-p $progress --no-update-check -F $min_isoform_fraction ";
$cuffdiff_param .= "-u " if $multi_read_correct;
$cuffdiff_param .= "-M $mask_file " if $mask_file;
$cuffdiff_param .= "-b $bias_correct " if $bias_correct;

my ($label_param, $bams_param);
my @labels;
my @bam_by_group;
foreach my $group (sort keys %groups) {
	push @labels, $group;
	my @group_bams = ();
	foreach my $replicate(sort {$a <=>$b } keys %{$groups{$group}}) {
		push @group_bams, $groups{$group}{$replicate};
	}
	push @bam_by_group, join ",", @group_bams;
}
$label_param = join ",", @labels;
$bams_param = join " ", @bam_by_group;


if ($args_file) {
	open my $args_fh, "<", $args_file or die "cannot open $args_file\n";
	while (<$args_fh>) {
		chomp;
		$cuffdiff_param .= "$_ ";		
	}
}


warn "Staring cuffdiff...\n";
my $stdout = "cuffdiff.stdout.log";
my $stderr = "cuffdiff.stderr.log";
my $command =
"cuffdiff -o $out_folder/ $cuffdiff_param -L $label_param $gtf_file $bams_param >$out_folder/$stdout 2> $out_folder/$stderr";

warn "\t$command\n";
my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
  run( command => $command, verbose => 0 );
if ($success) {
	warn "\tDone!\n";
}
else {
	my @stderrs = @$stderr_buf;
	warn "Something went wrong:\n@stderrs";
}		




sub file_check {
	my $file = shift;
die "File $file does not exists\n" 
	unless -e  $file;
die "File:$file is empty\n" 
	unless -s $file;	
}