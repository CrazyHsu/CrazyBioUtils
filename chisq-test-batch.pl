#!/usr/bin/perl -w

use strict;
use 5.010;
use Getopt::Long;
# use Data::Printer;

my $usage = <<USAGE;
Usage: chisq-test.pl [options...] [options]files
written by corephi, group:276151571
this program is used to do chisq-test one by one.

 Options:
  -a|--header           input files containd header, default 0.
  -m|--mode             mode used for comparsion, 'seq' means to 
                        compare one by one, for example 1-2, 2-3,
                        3-4. 'combination' means to compare
                        by all possible pairs, for example 1-2, 1
                        -3, 1-4, 2-3, 2-4, 3-4
						
Note: the input file must containd 3 colums
condition  ture   false
con1       21     19
con2       77     214
con3       156    1071
USAGE

#pare the arguments
my $header_flag          = 0;
my $mode = 'seq';
die $usage
  unless GetOptions(
    "a|header"            => \$header_flag,
    "m|mode:s"            => \$mode,
  );
# p @ARGV;
die "no files inputed!\n$usage\n" unless @ARGV;
foreach my $file (@ARGV) {
	die "$file does not exits\n" unless -e $file;
}

foreach my $file (@ARGV) {
	open my $in_fh, "<", $file;
	readline $in_fh if $header_flag;
	my %data;
	my @conds;
	while (<$in_fh>) {
		chomp;
		next if /^\s/;
		my ($cond, $true,$false) = split /\t/;
		push @conds, $cond;
		$data{$cond}{true} = $true;
		$data{$cond}{false} = $false;
	}
	close $in_fh;
	# p %data;
	
	open my $out_fh, ">", $file.".test.txt";
	print $out_fh "Comparasioin\tPvalue\tCorretion\n";
	if ($mode eq 'seq') {
		foreach my $i(0 .. $#conds - 1) {
			my $cond1 = $conds[$i];
			my $cond2 = $conds[$i+1];
			my $cond1_t = $data{$cond1}{true};
			my $cond1_f = $data{$cond1}{false};
			my $cond2_t = $data{$cond2}{true};
			my $cond2_f = $data{$cond2}{false};
			my $correction = 1;
			$correction = 0 if $cond1_t > 5 && $cond1_f > 5 && $cond2_t > 5 && $cond2_f > 5; 
			my $pvalue = chisq_test($cond1_t, $cond1_f, $cond2_t, $cond2_f , $correction);
			print $out_fh "$cond1 vs $cond2\t$pvalue\t$correction\n";
		}
	}else{
		foreach my $i(0 .. $#conds - 1)	{
			foreach my $j ($i + 1 .. $#conds) {
				my $cond1 = $conds[$i];
				my $cond2 = $conds[$j];
				my $cond1_t = $data{$cond1}{true};
				my $cond1_f = $data{$cond1}{false};
				my $cond2_t = $data{$cond2}{true};
				my $cond2_f = $data{$cond2}{false};
				my $correction = 1;
				$correction = 0 if $cond1_t > 5 && $cond1_f > 5 && $cond2_t > 5 && $cond2_f > 5; 
				my $pvalue = chisq_test($cond1_t, $cond1_f, $cond2_t, $cond2_f , $correction);
				print $out_fh "$cond1 vs $cond2\t$pvalue\t$correction\n";
			}
		}		
	}
	close $out_fh;
}


sub chisq_test {
	use Statistics::R;
	my ($con1_t, $con1_f , $con2_t, $con2_f, $correct) = @_;
	$correct = $correct ? "T" : "F";  
	
	my $R_code = <<CHISQ;
data <- data.frame( 'true' = c($con1_t, $con2_t),
                   'false' = c($con1_f, $con2_f)
				   )
result <- as.numeric(chisq.test(data, correct = $correct)[3])
CHISQ

	
	my $R = Statistics::R->new();
	$R->run($R_code);
	my $result = $R->get('result'); 
	$R->stop;
	return $result;

}



