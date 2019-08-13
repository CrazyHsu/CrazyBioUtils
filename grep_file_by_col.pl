#!/usr/bin/perl -w

use strict;
use 5.010;
use Getopt::Long;

#usage information
my $usage = <<USAGE;
grep_file_by_col.pl V1.0, written by corephi
This program is similar to linux bash "cut -f", the differece
is this can not only output by ordered col numbers, but also 
can output by ordered col names
----------------------------------------------------------
Usage: grep_file_by_col.pl -a in.txt -c col -o outfile
 Options:
  -i    input file name, default STDIN
  -o    output file name, default STDOUT
  -c    col number or name, number such as '1:3,4:4,5:8' is 
        accepted,  ',' - separater, ':' - range  
  -m    'name' or 'number', default 'number'
Note: This scrpit can also be used in pipeline.
cat in.txt | grep_file_by_col.pl -c 3,2,1 
     

USAGE
my $in_file        = '';
my $column_str      = '';
my $outfile      = '-';
my $mode         = 'number';
die $usage
  unless GetOptions(
    "i:s" => \$in_file,
    "c=s" => \$column_str,
    "m:s" => \$mode,
    "o:s" => \$outfile,
  );

die "a and b can not get from STDIN at the same time\n" 
	if $in_file && $in_file eq '-';

 
my $outfile_fh;
if ($outfile && $outfile ne '-') {
	open $outfile_fh, ">", $outfile or die "cannot open file $outfile:$!\n";
}else {
    $outfile_fh = *STDOUT;
}


#append content to a file
my $in_fh;
if ($in_file && $in_file ne '-') {
    die "a file does not exists\n" unless -e $in_file;
    open $in_fh, "<", $in_file or die "cannot open file $in_file:$!\n";
}
else {
    $in_fh = *STDIN;
}



#groups

my @columns = ();   #0-based
my %columns_tmp = ();
if ($column_str){
	if ($mode eq 'number') {
		my @column_strs = split /,/, $column_str;
		foreach my $column (@column_strs) {
			my ($start, $end ) = split /:/, $column ;
			$end = $start unless $end;
			my @tmp = ();
			my $revere_flag = 0;
			if ($start > $end) {
				($start, $end) = ($end, $start);
				$revere_flag = 1;
			}
			foreach my $col ($start .. $end) {
				if (exists $columns_tmp{$col}) {
					next;
				}else{
					$columns_tmp{$col} = 1;
					push @tmp, $col - 1;
				}
			}
			@tmp = reverse @tmp if $revere_flag;
			push @columns, @tmp;			

		}
	}else{
		my $header = readline $in_fh;
		$header =~ s/\r?\n//;
		my @file_col_names = split /\t/, $header;
		my %file_cols = ();
		foreach my $col (0 .. $#file_col_names) {
			$file_cols{$file_col_names[$col]} = $col;
		}
		
		my @seed_col_names = split /,/, $column_str;
		foreach my $col_name (@seed_col_names) {
			if (exists $file_cols{$col_name}) {
				push @columns, $file_cols{$col_name} ;			
			}else{
				die "column name:$col_name does not exists in input file\n";
			}		
			if (exists $columns_tmp{$col_name}) {
				next;
			}else{
				$columns_tmp{$col_name} = 1;

			}			
		}
		
		print $outfile_fh join "\t", @file_col_names[@columns];
		print $outfile_fh "\n";

	}

}else{
	die "no column is specified\n";
}


while (<$in_fh>) {
    s/\r?\n//;
    my @lines = split /\t/;
	print $outfile_fh join "\t", @lines[@columns];
	print $outfile_fh "\n";
}


