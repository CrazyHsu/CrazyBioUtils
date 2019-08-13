#!/usr/bin/perl -w

use strict;
use 5.010;
use Getopt::Long;
use Bio::DB::Fasta;
use Bio::Seq;
use Bio::SeqIO;

#usage information
my $usage = <<USAGE;
dump_fasta_seq V1.0, written by corephi
This program is used to plot the reads coverage of each chromosome
More scripts? Please join "bioinformatics*CN" QQ group: 276151571
If you have any probel or suggestions about this program, please
mail to: hpyu\@genetics.ac.cn
----------------------------------------------------------
Usage: dump_fasta_seq -l loci.txt -f genome.fa -o result.fa

 Options:
  -l|loci     a table delimited file contains GID and loci, such as:
              AT1G01010	Chr1:3630-5899
			  and loci must have either of following format:
			  \$id:\$start,\$stop
              \$id:\$start..\$stop
              \$id:\$start-\$stop
              \$id:\$start,\$stop/\$strand
              \$id:\$start..\$stop/\$strand
              \$id:\$start-\$stop/\$strand
              \$id/\$strand")
  -f|fasta    a fasta file or direcotry of the genome
  -o|out      output fasta file name;    

Note: the loci file should not contain a header
USAGE

#parse parameter
my $fasta_file = '';
my $loci_file  = '';
my $out_file   = 'out.fa';
die $usage
  unless GetOptions(
    "l|loci=s"  => \$loci_file,
    "f|fasta=s" => \$fasta_file,
    "o|out=s"   => \$out_file,
  );

die "fasta file:$fasta_file does not exists\n" unless -e $fasta_file;
die "loci file:$loci_file does not exists\n"   unless -e $loci_file;

my $fa_db = Bio::DB::Fasta->new($fasta_file);
my $seq_io = Bio::SeqIO->new( '-format' => 'Fasta', -file => ">$out_file" );

open my $loci_fh, "<", $loci_file or die "cannot open locifile:$!\n";
while (<$loci_fh>) {
    chomp;
    my ( $gid, $loci ) = split /\t/;
    my $seqstr  = $fa_db->seq($loci);
    my $seq_obj = Bio::Seq->new(
        -display_id => $gid,
        -seq        => $seqstr
    );
    $seq_io->write_seq($seq_obj);
}
close $loci_fh;
