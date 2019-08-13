#!/usr/bin/perl -w
use strict;
use 5.010;
# use Data::Printer;
use File::Basename;
use IPC::Cmd qw[can_run run];
use Cwd 'abs_path';
use Getopt::Long;

my $usage = <<USAGE;
SYSNOPSIS
cuffmerge.pl [options] bw2db file|glob

 Options:
   -f|--gtf-flag      default NA, can be setted as'accepted_hits'
   -g|--ref-gtf       An optional “reference” annotation GTF
   -s|–-ref-seq       genomic DNA sequences for the reference
   -o|--output        output folder for removed reads 
   -p|--progress      progress, default 8
USAGE
my $mode       = 'pe';
my $ref_gtf_file = '';
my $genome_seq_file = '';
my $out_folder = dirname './';
my $progress   = 8;
my $gtf_flag = '';
die $usage
  unless GetOptions(
    "f|gtf-flag:s"       => \$gtf_flag,
    "g|ref-gtf:s"       => \$ref_gtf_file,
    "s|ref-seq:s"       => \$genome_seq_file,
    "o|output:s"   => \$out_folder,
    "p|progress=i" => \$progress,
  );

#############################################################################
#Parameters
#############################################################################  
file_check($genome_seq_file) if $genome_seq_file;
file_check($ref_gtf_file) if $ref_gtf_file;

can_run('cuffmerge') or die 'cuffmerge is not installed!';
$out_folder =~ s/[\/|\|]+$//;
mkdir $out_folder unless -d $out_folder;
my $log_folder = $out_folder."/logs";
mkdir $log_folder unless -d $log_folder;



my @gtfs = ();
foreach my $glob_str(@ARGV) {
	
	#pick up directory
	my @gtf_tmp = ();
	if ($gtf_flag) {
		@gtf_tmp = map {abs_path("$_/$gtf_flag.gtf")} 
				grep {-d $_ and -s "$_/$gtf_flag.gtf"} 
				glob $glob_str;			
	}else{
		@gtf_tmp = map {abs_path($_)} 
			grep {-s $_} 
				glob $glob_str;	
	}
	
	#pick up gtfs	
	push @gtfs, @gtf_tmp if @gtf_tmp;
}

die "No gtfs were found, please check the glob string\n" 
	unless @gtfs;
  
my $assembly_list_file = "$out_folder/assembly_list.txt";
open my $assembly_list_fh, ">", $assembly_list_file
	or die "Cannot create assembly list:$assembly_list_file$!\n";
foreach my $gtf (@gtfs) {
	print $assembly_list_fh "$gtf\n";
}
close $assembly_list_fh;

my $cuffmerge_parameter = "-p $progress";
$cuffmerge_parameter .= " -o $out_folder";
$cuffmerge_parameter .= " -g $ref_gtf_file" if $ref_gtf_file;
$cuffmerge_parameter .= " -s $genome_seq_file" if $genome_seq_file;


my $stdout = "cuffmerge.stdout.log";
my $stderr = "cuffmerge.stderr.log";
my $command =
"cuffmerge  $cuffmerge_parameter $assembly_list_file >$log_folder/$stdout 2> $log_folder/$stderr";

warn "\t$command\n";
my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
  run( command => $command, verbose => 0 );
if ($success) {
	warn "\tDone!\n";
}
else {
	warn "Something went wrong:\n$stderr_buf";
}
		



sub file_check {
	my $file = shift;
die "File $file does not exists\n" 
	unless -e  $file;
die "File:$file is empty\n" 
	unless -s $file;	
}

