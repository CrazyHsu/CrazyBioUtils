#!/usr/bin/perl -w
use strict;
use 5.010;
# use Data::Printer;
use Getopt::Long;
use File::Copy;
use File::Basename;


my $usage = <<USAGE;
SYSNOPSIS
brename.pl V1.0.0, written by corephi

This program is bath rename files by given replacement rules.
replacement rules must be perl regular exprssion. This can be
used for rename the fastq files.


----------------------------------------------------------
More scripts? Join "bioinformatics*CN" QQ group: 276151571
If you have any probel or suggestions about this program, 
please mail to: hpyu\@genetics.ac.cn
----------------------------------------------------------
Usage: brename.pl [options] file|glob

 Options:
   -f|--file          file store batch rename pattern,batch 
                      rename pattern must be stored as two 
					  column: source_pattern \t replacement
   -s|--pattern       source string, regular expression such 'ab'.
   -r|--replacement   replacement string.
   -p|print           print only
   -h|--help               

Samples:

brename.pl -s 'SRR0001' -r 'Heart' *.gz
brename.pl -f string *.gz


USAGE
my $pattern_file = '';
my $source_reg = '';
my $replacement = '';
my $print_flag = 0;
my $help_flag = 0;

die $usage
  unless GetOptions(
    "f|file:s"   => \$pattern_file,  
    "s|pattern:s"   => \$source_reg,
    "r|replacement:s"   => \$replacement,
    "p|print"   => \$print_flag,
    "h|help" => \$help_flag,
  );
  
die $usage if $help_flag;
#############################################################################
#Parameters
############################################################################# 
#check file existance
my @patterns;
if ($pattern_file) {
	file_check($pattern_file);	
	@patterns = read_pattern($pattern_file)
}else{
	die "source pattern must be specified" unless $source_reg;
	@patterns = ([$source_reg, $replacement])
}

#checking the bam files
my @files = ();
foreach my $glob_str(@ARGV) {
	#pick up directory
	my @file_tmp = glob $glob_str;	
				
	#pick up bams	
	push @files, @file_tmp if @file_tmp;
}

die "No files were found, please the the glob string\n" 
	unless @files;
	
#Preprocess the patterns
my %replacement_files = ();	
my %file_mapping = ();
foreach my $file (@files) {
	my $replacement_file = '';
	my ($name,$path,$suffix) = fileparse($file,());
	my $source = $path.$name;
	foreach my $i (0 .. $#patterns) {
		my ($pattern, $replacement) = @{$patterns[$i]};
		if ($name =~ /$pattern/) {
			my $replacement_name = $name;
			# p $replacement_name;
			# p $replacement;
			$replacement_name =~ s/$pattern/$replacement/;
			# p $replacement_name;
			$replacement_file = $path . $replacement_name;

			if (exists $replacement_files{$replacement_file} ) {
				die "'$pattern' to '$replacement' may conflict with $replacement_files{$replacement_file} ";
			}else{
				$replacement_files{$replacement_file} = "'$pattern' to '$replacement'";
			}
			last;
		}
	}

	if ($replacement_file) {
		$file_mapping{$source}{source} = $source ;		
		$file_mapping{$source}{dest} = $replacement_file ;		
	}
	
}
# p %file_mapping;
warn "No files math to the pattern\n" unless %file_mapping;
foreach my $file (sort keys %file_mapping) {
	my $replacement = $file_mapping{$file}{dest};
	my $source = $file_mapping{$file}{source};
	warn "$source" . " -> ".  $replacement. "\n";
	move($source, $replacement) unless $print_flag;
}



#===  FUNCTION  ================================================================
#         NAME: Read pattern
#      PURPOSE: 
#   PARAMETERS: $file: string
#      RETURNS: @patterns: two-dimintional array
#			         $patterns[$i][0] = $pattern
#			         $patterns[$i][1] = $replacement
#===============================================================================
sub read_pattern {
	my $file = shift;
	open my $in_fh, "<", $file
		or die "Cannot open file:$file\n";
	my @patterns;
	while (<$in_fh>) {
		chomp;
		my ($pattern , $replacement) = split /\t/;
		push @patterns, [$pattern , $replacement];
	}
	return @patterns;
}



#===  FUNCTION  ================================================================
#         NAME: file_check
#      PURPOSE: check wheather the file exists
#   PARAMETERS: $file: string
#      RETURNS: NA
#===============================================================================
sub file_check {
	my $file = shift;
die "File $file does not exists\n" 
	unless -e  $file;
die "File:$file is empty\n" 
	unless -s $file;	
}

