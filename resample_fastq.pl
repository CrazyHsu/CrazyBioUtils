#!/usr/bin/perl
use 5.010001;
use strict;
use warnings;
use File::Basename;
use Algorithm::Numerical::Sample qw(sample);
use Getopt::Long;

my $usage = <<USAGE;
resample_fastq V1.0, written by corephi
This program is used for resample the sample to your specified
size, and it can also trim the sequcence at the same time.
This is useful, when you try to compare between libraries. 
----------------------------------------------------------
More scripts? Join "bioinformatics*CN" QQ group: 276151571
If you have any probel or suggestions about this program, 
please mail to: hpyu\@genetics.ac.cn
----------------------------------------------------------
Usage: resample_fastq -m [left|right|both] -l 50 1.fq 2,fq


 Options:
  -m|--mode         either "left", "right", the left or right length 
                    to keep, , default "right"
  -s|--sample-size  sample size. set to 0, if you want to keep all.
  -l|--length       the length to trim, default 0;

USAGE
my $mode        = "right";
my $len         = 0;
my $sample_size = 0;
die $usage
  unless GetOptions(
    "m|mode:s"        => \$mode,
    "s|sample-size:i" => \$sample_size,
    "l|length:i"      => \$len,
  );
die "You must specify one of length or sample size.\n"
  if ( $sample_size == 0 && $len == 0 );

my @fastq_files = @ARGV;
die "You must specify at least one fastq file.\n"
  unless @fastq_files;
my @suffixlist = qw(fq fastq);

foreach my $fastq_file (@fastq_files) {
    my ( $name, $path, $suffix ) = fileparse( $fastq_file, @suffixlist );
    my $out_file = "${name}${mode}.$len.$suffix";

    #get the total number of fastq
    my $total_num = 0;
    open my $in_fh, "<", $fastq_file
      or die "cannot open file:$!";
    while (<$in_fh>) {
        $total_num++;
    }
    close $in_fh;
    $total_num /= 4;

    #sample;
    my %sample;
    if ( $sample_size == 0 || $sample_size >= $total_num ) {
        $sample{$_} = 1 foreach 0 .. $total_num - 1;
    }
    else {
        $sample{$_} = 1 foreach sample(
            -set         => [ 0 .. $total_num - 1 ],
            -sample_size => $sample_size
        );
    }

    #openfiles
    open $in_fh, "<", $fastq_file;
    open my $out_fh, ">", $out_file
      or die "cannot create file:$!";
    while (<$in_fh>) {
        chomp;
        my $line = $_;
        my $seq_num = int( ( $. - 1 ) / 4 );
        if ( exists $sample{$seq_num} ) {
            given ($mode) {
                when ('left') {
                    if ( $. % 2 == 0 ) {
                        if ( length($line) > $len ) {
                            my $left = substr( $line, 0, $len );
                            print $out_fh $left, "\n";
                        }
                        else {
                            print $out_fh "$line\n";
                        }
                    }
                    else {
                        print $out_fh "$line\n";
                    }
                }

                when ('right') {
                    if ( $. % 2 == 0 ) {
                        if ( length($line) > $len ) {
                            my $right = substr( $line, -$len );
                            print $out_fh $right, "\n";
                        }
                        else {
                            print $out_fh "$line\n";
                        }
                    }
                    else {
                        print $out_fh "$line\n";
                    }
                }

                default { die $usage; }

            }
        }
        else {

        }

    }
    close $in_fh;
    close $out_fh;

}

