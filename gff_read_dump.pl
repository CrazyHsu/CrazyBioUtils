#!/usr/bin/perl -w

use strict;
use 5.010001;
use Getopt::Long;
use YAML;
use Bio::SeqIO;
use Bio::DB::SeqFeature::Store;

#usage information
my $usage = <<USAGE;
NMD_signature V1.0, written by corephi
Given a gtf and fasta file, this program will get intron sequcence 
in fasta format
More scripts? Please join "bioinformatics*CN" QQ group: 276151571
If you have any probel or suggestions about this program, please
mail to: hpyu\@genetics.ac.cn
----------------------------------------------------------
Usage: NMD_signature -f genome.fa -g gene_models.gtf -o result.fa

 Options:
  -u      uORF.fa, sORF Finder's result fasta
  -g      gene_models.gff
  -o      output fasta files, default out.fa;1

USAGE

#parse parameter
my $uORF_fa_file = '';
my $gff_file     = "";
my $out_file     = "out.fa";
die $usage
  unless GetOptions(
    "g=s" => \$gff_file,
    "u=s" => \$uORF_fa_file,
    "o=s" => \$out_file,
  );
die "uORF fasta file:$uORF_fa_file doesn't exists\n$usage" unless $uORF_fa_file;
die "gene_models file:$gff_file doesn't exists\n$usage"    unless $gff_file;

my $uORF_fa = Bio::SeqIO->new(
    -file   => $uORF_fa_file,
    -format => 'Fasta'
);
my %data;
while ( my $seq = $uORF_fa->next_seq() ) {
    my ( $tid, $start, $strand, $score ) = split /#/, $seq->id,
      my $len = $seq->length;
    my $end   = $start + $len - 1;
    $data{$tid}{uORF} = {
        start  => $start,
        end    => $end,
        len    => $len,
        strand => $strand,
        score  => $score,
    };
}

my $gene_model_db = Bio::DB::SeqFeature::Store->new(
    -adaptor           => 'memory',
    -dsn               => $gff_file,
    -autoindex         => 1,
    -index_subfeatures => 1,
);

my $iterator = $gene_model_db->get_seq_stream( -type => 'mRNA', );

while ( my $feature = $iterator->next_seq ) {
    my $id = $feature->name;
    print Dump $feature;

# $id =~ s/\.\.\./\.t\./;
# my ($type, $gene, $ref_num, $class_code, $t_count, $orm_num) = split /\./, $id;
# my $tid = join ".", ($gene, $ref_num, $class_code, $t_count);
# my @utr5 = feat2array($feature, 'five_prime_utr');
# my @utr3 = feat2array($feature, 'three_prime_utr');
# my @cds = feat2array($feature, 'CDS');

}

sub feat2array {
    my ( $feature, $type ) = @_;
    my @sub_feats = $feature->get_SeqFeatures($type);

    my %loci;
    foreach my $feat (@sub_feats) {    ##get the exon location;
        my $start = $feat->start;
        my $end   = $feat->end;
        $loci{$start} = $end;
    }
    my @exons;                         ##restore the exon info
    my $i = 0;
    foreach my $key ( sort { $a <=> $b } keys %loci ) {
        my $start = $key;
        my $end   = $loci{$key};
        my $len   = $end - $start + 1;
        $exons[$i][0] = $start;
        $exons[$i][1] = $end;
        $exons[$i][2] = $len;
        $i++;
    }
    return @exons;
}
