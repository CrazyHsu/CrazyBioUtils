#!/usr/bin/perl -w

use strict;
use 5.010;
use Bio::SeqIO;
use Getopt::Long;
use File::Basename;


my $usage = <<USAGE;
SYSNOPSIS
split_fa_by_len.pl -l 300 -o ./ file|glob

 Options:
   -l|--length    cutoff length, default 70 
   -o|--output    output folder for removed reads 
USAGE
my $out_folder = dirname './';
my $cutoff_len = 300;

die $usage
  unless GetOptions(
    "o|output=s"   => \$out_folder,
	"l|cut-length:i"       =>  \$cutoff_len,
  );
my $glob = shift;
die "You must spcified at least one fasta\n" unless $glob;
my @fas  = glob "$glob";
@fas = grep {-e $_} @fas;
die "At least one file specified\n" if @fas < 1;
$out_folder =~ s/[\/|\|]+$//;
mkdir $out_folder unless -e $out_folder;



my @suffx  = qw (.fa .fas .fasta);
foreach my $in_fasta(@fas) {
	my $basename =  basename( $in_fasta, @suffx );
	my $out_long_file = "$out_folder/$basename.long.$cutoff_len.fa";
	my $out_short_file = "$out_folder/$basename.short.$cutoff_len.fa";

	my $in  = Bio::SeqIO->new(-file => $in_fasta ,
						 -format => 'Fasta');
	my $out_long = Bio::SeqIO->new(-file => ">$out_long_file" ,
						 -format => 'Fasta');
	my $out_short = Bio::SeqIO->new(-file => ">$out_short_file" ,
						 -format => 'Fasta');
						 
	while ( my $seq = $in->next_seq() ) {
		my $len = $seq->length();
		if ($len >= $cutoff_len) {
			$out_long->write_seq($seq);
		}else{
			$out_short->write_seq($seq);			
		}
	}	
}

