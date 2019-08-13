#!/usr/bin/perl -w

use strict;
use 5.010;
use Getopt::Long;
use Statistics::R;

# use Data::Printer;

my $usage = <<USAGE;
Usage: chisq-test.pl [options...] [options]files
written by corephi, group:276151571
this program is used to do chisq-test one by one.

 Options:
	-t1 
	-t2 
	-f1
	-f2
	-c  correct

USAGE

#pare the arguments
my ($con1_t, $con1_f , $con2_t, $con2_f, $correct) ;
die $usage
  unless GetOptions(
    "t1:i"            => \$con1_t,
    "t2:i"            => \$con2_t,
    "f1:i"            => \$con1_f,
    "f2:i"            => \$con2_f,
    "c"            => \$correct,
  );

  
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
print $result, "\n";




