#!/usr/bin/perl -w

use strict;
use 5.010001;
use Getopt::Long;
use Bio::SeqIO;
use Bio::Seq;
use Bio::DB::Fasta;

#usage information
my $usage = <<USAGE;
gff_read_dump V1.0, written by corephi
Given a gtf and fasta file, this program will get intron sequcence 
in fasta format
More scripts? Please join "bioinformatics*CN" QQ group: 276151571
If you have any probel or suggestions about this program, please
mail to: hpyu\@genetics.ac.cn
----------------------------------------------------------
Usage: get_intron_from_gtf -f genome.fa -g gene_models.gtf -o result.fa

 Options:
  -f      a fasta file contains the genome sequence
  -g      gene_models.gtf
  -m      either 'transcipts', "cds", "exon", 'intron'
  -o      output fasta files, default out.fa;

USAGE

#parse parameter
my $fa_file  = '/home/hpyu/tair_data/Col_tair10.fa';
my $gtf_file = "/home/hpyu/tair_data/TAIR10_MOD.gff";
my $out_file = "out.fa";
my $mode     = '';
die $usage
  unless GetOptions(
    "f=s" => \$gtf_file,
    "g=s" => \$gtf_file,
    "o=s" => \$out_file,
    "m=s" => \$mode,
  );
die "$fa_file doesn't exists\n"  unless $fa_file;
die "$gtf_file doesn't exists\n" unless $gtf_file;
my %modes = map { $_ => 1 } qw(transcipts exon intron intergenic);
$mode = lc $mode;
die "mode $mode is not supported\n" unless exists $modes{$mode};

my $fa_db  = Bio::DB::Fasta->new($fa_file);
my $fa_out = Bio::SeqIO->new(
    -file   => ">$out_file",
    -format => 'Fasta'
);
my $gene_models_rf = read_gtf( $gtf_file, 'gene' );
my %gene_models = %$gene_models_rf;

given ($mode) {
    when ("transcipts") { }
    when ("cds")        { }
    when ("exon")       { }
    when ("intron")     { }
    default {
    }
}

sub ()

  foreach my $gid ( keys %gene_models ){ foreach my $tid () }

#===  FUNCTION  ================================================================
#         NAME: read_gtf
#      PURPOSE: given a gtf file, retreive the information into a hash
#   PARAMETERS: $gtf_file: string
#   			$mode: only one of 'gene' or 'transcript'
#      RETURNS: $hash_rf: a hash reference stored the information
#  DESCRIPTION: There are two mode of this function, one is 'gene', and the other is 'transcript'. The diffierence bewteen the two mode is the hash structure.
#			'gene' mode: $trasncripts{$gene_id}{$trasncript_id}....
# 			'transcript' mode: $transcripts{$transcript_id}, but have one more key than gene mode: $transcripts{$transcript_id}{gene_id}.
#			THe common keys are:
#			..{$transcript_id}{start}
#			..{$transcript_id}{end}
#			..{$transcript_id}{chr}
#			..{$transcript_id}{strand}
#			..{$transcript_id}{exons}
#			..{$transcript_id}{introns}
#
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
  sub read_gtf {
    my ( $gtf_file, $mode ) = @_;
    die "unknon mode:$mode, it must be one of 'gene' or 'transcript'\n"
      unless ( $mode eq 'gene' or $mode eq 'transcript' );
    my %transcripts;
    open my $gtf_fh, "<", $gtf_file or die "Cannot open gtf $gtf_file:$!\n";

    #read the gtf info
    while (<$gtf_fh>) {
        next if /^#/;
        chomp;
        my @lines     = split /\t/;
        my $attribute = pop @lines;
        my $type      = lc $lines[2];

        # say $attribute;
        my @temp = split ';', $attribute;
        my %attribttes;
        foreach my $attr (@temp) {
            $attr =~ s/^\s+//;
            next unless $attr;
            my ( $key, $value );
            if ( $attr =~ /(\S+) "(.*)"/ ) {
                ( $key, $value ) = ( $1, $2 );

                # say "$key\t\t\t$value";
                $attribttes{$key} = $value;
            }
            else {
                warn "$attr:is not well formatted\n";
            }
        }
        my $gene_id       = $attribttes{gene_id};
        my $transcript_id = $attribttes{transcript_id};
        my $chr           = $lines[0];
        my $start         = $lines[3];
        my $end           = $lines[4];
        my $strand        = $lines[6];
        if ( $type eq 'transcript' ) {

            if ( $mode eq 'gene' ) {
                $transcripts{$gene_id}{$transcript_id}{start} = $start;
                $transcripts{$gene_id}{$transcript_id}{end}   = $end;
            }
            else {
                $transcripts{$transcript_id}{start}   = $start;
                $transcripts{$transcript_id}{end}     = $end;
                $transcripts{$transcript_id}{gene_id} = $gene_id;
            }
        }
        elsif ( $type eq 'exon' ) {
            if ( $mode eq 'gene' ) {

                #check and store chr
                die
"same transcript with different chr, pleae check $transcript_id\n"
                  if ( exists $transcripts{$gene_id}{$transcript_id}{chr}
                    && $transcripts{$gene_id}{$transcript_id}{chr} ne $chr );
                $transcripts{$gene_id}{$transcript_id}{chr} = $chr;

                #check and store strand
                die
"same transcript with different strand, pleae check $transcript_id\n"
                  if ( exists $transcripts{$gene_id}{$transcript_id}{strand}
                    && $transcripts{$gene_id}{$transcript_id}{strand} ne
                    $strand );
                $transcripts{$gene_id}{$transcript_id}{strand} = $strand;

                #store exons as hash
                $transcripts{$gene_id}{$transcript_id}{exons}{$start} = $end;
            }
            else {
                #check and store chr
                die
"same transcript with different chr, pleae check $transcript_id\n"
                  if ( exists $transcripts{$transcript_id}{chr}
                    && $transcripts{$transcript_id}{chr} ne $chr );
                $transcripts{$transcript_id}{chr} = $chr;

                #check and store strand
                die
"same transcript with different strand, pleae check $transcript_id\n"
                  if ( exists $transcripts{$transcript_id}{strand}
                    && $transcripts{$transcript_id}{strand} ne $strand );
                $transcripts{$transcript_id}{strand} = $strand;

                #store exons as hash
                $transcripts{$transcript_id}{exons}{$start} = $end;
            }

       # print "$gene_id\t$transcript_id\t$type\t$chr\t$start\t$end\t$strand\n";

        }
        elsif ( $type eq 'cds' ) {
            if ( $mode eq 'gene' ) {

                #check and store chr
                die
"same transcript with different chr, pleae check $transcript_id\n"
                  if ( exists $transcripts{$gene_id}{$transcript_id}{chr}
                    && $transcripts{$gene_id}{$transcript_id}{chr} ne $chr );
                $transcripts{$gene_id}{$transcript_id}{chr} = $chr;

                #check and store strand
                die
"same transcript with different strand, pleae check $transcript_id\n"
                  if ( exists $transcripts{$gene_id}{$transcript_id}{strand}
                    && $transcripts{$gene_id}{$transcript_id}{strand} ne
                    $strand );
                $transcripts{$gene_id}{$transcript_id}{strand} = $strand;

                #store exons as hash
                $transcripts{$gene_id}{$transcript_id}{cds}{$start} = $end;
            }
            else {
                #check and store chr
                die
"same transcript with different chr, pleae check $transcript_id\n"
                  if ( exists $transcripts{$transcript_id}{chr}
                    && $transcripts{$transcript_id}{chr} ne $chr );
                $transcripts{$transcript_id}{chr} = $chr;

                #check and store strand
                die
"same transcript with different strand, pleae check $transcript_id\n"
                  if ( exists $transcripts{$transcript_id}{strand}
                    && $transcripts{$transcript_id}{strand} ne $strand );
                $transcripts{$transcript_id}{strand} = $strand;

                #store exons as hash
                $transcripts{$transcript_id}{cds}{$start} = $end;
            }

       # print "$gene_id\t$transcript_id\t$type\t$chr\t$start\t$end\t$strand\n";

        }
        else {
            next;
        }
    }
    close $gtf_fh;

    #get the intron information
    if ( $mode eq 'gene' ) {
        foreach my $gene_id ( sort keys %transcripts ) {
            foreach my $transcript_id ( keys %{ $transcripts{$gene_id} } ) {
                my $strand = $transcripts{$gene_id}{$transcript_id}{strand};
                my %exon_loci =
                  %{ $transcripts{$gene_id}{$transcript_id}{exons} };

                my ( $exon_rf, $intron_rf ) =
                  get_intron( \%exon_loci, $strand );
                $transcripts{$gene_id}{$transcript_id}{exons}   = $exon_rf;
                $transcripts{$gene_id}{$transcript_id}{introns} = $intron_rf;

                # print Dump $exon_rf;
                # print Dump $intron_rf;

            }
        }
    }
    else {
        foreach my $transcript_id ( sort keys %transcripts ) {
            my $strand    = $transcripts{$transcript_id}{strand};
            my %exon_loci = %{ $transcripts{$transcript_id}{exons} };

            my ( $exon_rf, $intron_rf ) = get_intron( \%exon_loci, $strand );
            $transcripts{$transcript_id}{exons}   = $exon_rf;
            $transcripts{$transcript_id}{introns} = $intron_rf;

            # print Dump $exon_rf;
            # print Dump $intron_rf;
        }
    }
    return \%transcripts;
}

#===  FUNCTION  ================================================================
#         NAME: get_intron
#      PURPOSE: given a hash of exon loci, return a 2d array of exon and intron
#   PARAMETERS: $exon_loci: a hash reference with a format of $hash{$start} = $end
#   			$strand: '+' or '-';
#      RETURNS: @array: [\@exons, \@introns]
#  DESCRIPTION: given a hash of exon loci, return a 2d array of exon and intron.
#		  Note: the order of exon and intron is equal to the strand. that mens, if it is '-', the loci in exon is descreased.
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub get_intron {
    my ( $loci_rf, $strand ) = @_;
    my %loci = %$loci_rf;
    my @exons;
    my @introns;

    # print Dump \%loci;

    my $i = 0;
    foreach my $key ( sort { $a <=> $b } keys %loci ) {
        my $start = $key;
        my $end   = $loci{$key};
        $exons[$i][0] = $start;
        $exons[$i][1] = $end;
        $i++;
    }
    if ( @exons > 1 ) {
        for ( my $i = 0 ; $i < $#exons ; $i++ ) {
            my $start = $exons[$i][1] + 1;
            my $end   = $exons[ $i + 1 ][0] - 1;
            $introns[$i][0] = $start;
            $introns[$i][1] = $end;
        }
    }

    # print Dump \@exons;
    @exons   = reverse @exons   if $strand eq '-';
    @introns = reverse @introns if $strand eq '-';
    return ( \@exons, \@introns );
}

