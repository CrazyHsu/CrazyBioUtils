#!/usr/bin/perl -w

use strict;
use 5.010;
use YAML;
use Getopt::Long;
use Statistics::R;

my $usage = <<USAGE;
SYSNOPSIS
cmp_edgeR_result.pl [options] -o out.tiff -L label1,label2 edgeR1.xls edgeR2.xls

 Options:
   -m --mode          'PValue', 'FDR' , default 'PValue';
   -t --threshold      significant threshold, default '0.05';
   -o --outfile       VennDiagram fileanme
   -L --label         Label for each transcriptome
Note: this file must contain header lines, and header lines contain "GID", "pvalue" and "FDR",
and the "GID" must be the first col;
USAGE

my $mode      = 'PValue';
my $out_file  = 'VennDiagram';
my $threshold = '0.05';
my $label;
die $usage
  unless GetOptions(
    "m|mode:s"      => \$mode,
    "o|outfile:s"   => \$out_file,
    "l|label=s"     => \$label,
    "t|threshold=s" => \$threshold,
  );
$out_file = $out_file . "_$mode" . "_$threshold.tiff";
my @edgeRs = @ARGV;
die $usage unless @edgeRs;
my @labels = split /,/, $label;

if ( @edgeRs == 1 ) {
    die "Only one edgeR.xls file\n";
}
elsif ( @edgeRs == 2 ) {
    my ( $edgeR1_file, $edgeR2_file ) = @edgeRs;
    my $edgeR1_rf      = read_xls($edgeR1_file);
    my $edgeR2_rf      = read_xls($edgeR2_file);
    my $edgeR1_id_file = $labels[0] . ".id.txt";
    my $edgeR2_id_file = $labels[1] . ".id.txt";
    write_id( $edgeR1_rf, $edgeR1_id_file );
    write_id( $edgeR2_rf, $edgeR2_id_file );
    plot2venn( $labels[0], $edgeR1_id_file, $labels[1], $edgeR2_id_file,
        $out_file );

}
else {
    print "This program only support two set Venndiagram";
}

#===  FUNCTION  ================================================================
#         NAME: write_id
#      PURPOSE: given a hash of transcriptome, write down the id the to a file;
#   PARAMETERS: $data_rf: a hash reference;
#   			$filename: the out put filename;
#      RETURNS: 1:0
#  DESCRIPTION: given a hash of transcriptome, write down the id to a file;
#		  Note:
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub write_id {
    my $data_rf  = shift;
    my $filename = shift;
    open my $out_fh, ">", $filename or die "Cannot create the file:$!";
    print $out_fh "id\n";
    foreach my $id ( sort keys %{$data_rf} ) {
        if ( $data_rf->{$id}{$mode} <= $threshold ) {
            print $out_fh "$id\n"
              if exists $data_rf->{$id}{logFC}
              && abs( $data_rf->{$id}{logFC} ) >= 1;
        }
    }
    close $out_fh;
}

#===  FUNCTION  ================================================================
#         NAME: read_xls
#      PURPOSE: given a xls file, retreive the information into a hash
#   PARAMETERS: $xls_file: string
#      RETURNS: $hash_rf: a hash reference stored the information
#  DESCRIPTION: The result hash constructor are:
#			..{$id}{$col2}
#			..{$id}{$col3}
#			..{$id}{$col..}
#
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub read_xls {
    my ($xls_file) = @_;
    my %data;
    open my $xls_fh, "<", $xls_file or die "Cannot open gtf $xls_file:$!\n";

    #read the gtf info
    my $header = readline $xls_fh;
    my @colnames = split /\t/, $header;

    #check the ID colums
    my $id_col = shift @colnames;
    while (<$xls_fh>) {
        chomp;
        my @colums = split /\t/;
        my $id     = shift @colums;
        foreach my $i ( 0 .. $#colnames ) {
            $data{$id}{ $colnames[$i] } = $colums[$i];
        }
    }
    return \%data;
}

#===  FUNCTION  ================================================================
#         NAME: plot2venn
#      PURPOSE: given a hash of transcriptome, write down the id the to a file;
#   PARAMETERS: $transcriptome_rf: a hash reference;
#   			$filename: the out put filename;
#      RETURNS: 1:0
#  DESCRIPTION: given a hash of transcriptome, write down the id to a file;
#		  Note:
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub plot2venn {
    use Statistics::R;
    my ( $label1, $file_1, $label2, $file_2, $out_file ) = @_;
    open my $R_fh, ">", "VennDiagram.R";
    my $R_code = <<RCODE;
library(VennDiagram)

$label1 <- read.table(file="$file_1", header=T)
$label2 <- read.table(file="$file_2", header=T)

data <- list('$label1' = $label1\$id,
  '$label2' = $label2\$id
  )
venn.plot <- venn.diagram(
	x = data,
	filename = "$out_file",
	output = TRUE,
	height = 3000,
	width = 3000,
	resolution = 300,
	compression = 'lzw',
	units = 'px',
	lwd = 6,
	fill = c("cornflowerblue", "darkorchid1"),
	alpha = 0.75,
	label.col = "white",
	cex = 3.5,
	fontfamily = "serif",
	fontface = "bold",
	cat.col = c("cornflowerblue", "darkorchid1"),
	cat.cex = 3,
	cat.fontfamily = "serif",
	cat.fontface = "bold",
	cat.dist = c(0.03, 0.03),
	cat.pos = c(-20, 14)
	);
RCODE
    print $R_fh $R_code;
    close $R_fh;
    my $R = Statistics::R->new();
    $R->run($R_code);
    my $error = $R->error;
    warn $error if $error;
    $R->stop;

}

