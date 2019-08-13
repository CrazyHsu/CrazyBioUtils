#!/usr/bin/perl -w

use strict;
use 5.010;
use YAML;
use Getopt::Long;
use Bio::DB::Sam;
use Statistics::R;

my $usage= <<USAGE;
SYSNOPSIS
reads_venn_from_bam.pl [options] -o out.tiff -L label1,label2 1.bam 2.bam

 Options:
   -o --outfile       VennDiagram fileanme
   -L --label         Label for each bam
USAGE

my $out_file = 'VennDiagram.tiff';
my $label;
die $usage 
	unless GetOptions (
			"o|outfile:s" => \$out_file,
			"l|label=s" => \$label,
			);
my @bams = @ARGV;
die $usage unless @bams;
my @labels = split /,/, $label;
die "labels number is not equal to bam number\n"
	if @labels != @bams;
	
for(my $i = 0; $i < @labels; $i++) {
	my $bam_file = $bams[$i];
	die "$bam_file does not exits\n" unless -e $bam_file;
}	

my @plot2venn_para = ();
for(my $i = 0; $i < @labels; $i++) {
	my $bam_file = $bams[$i];
	die "$bam_file does not exits\n" unless -e $bam_file;
	my $bam_rf = read_bam($bam_file);
	my $bam_id_file = $labels[$i].".id.txt";
	write_id($bam_rf, $bam_id_file);
	push @plot2venn_para, ($labels[$i], $bam_id_file)
}
	
if (@bams == 1){
	die "Only one bam.gtf file\n";
}elsif (@bams == 2){

	plot2venn(@plot2venn_para, $out_file);
}elsif (@bams == 3){
	plot3venn(@plot2venn_para, $out_file);
}elsif (@bams == 4) {	
	plot4venn(@plot2venn_para, $out_file);
}else{
	print "This program only support two set Venndiagram";
}



#===  FUNCTION  ================================================================
#         NAME: read_bam
#      PURPOSE: given a bam file, retreive read name into a hash
#   PARAMETERS: $bam_file: string 
#      RETURNS: $hash_rf: a hash reference stored the information
#  		
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub read_bam {
	my $bam_file = shift;
	my %reads;

	my $bam = Bio::DB::Bam->open($bam_file,"r");
	my $header = $bam->header;

	die $usage unless $bam;

	while (my $align = $bam->read1) {
		my $read_name =  $align->qname;
		$reads{$read_name} = 1;
	}
	return \%reads;
	
}

#===  FUNCTION  ================================================================
#         NAME: write_id
#      PURPOSE: given a hash of bam, write down the id the to a file;
#   PARAMETERS: $bam_rf: a hash reference;
#   			$filename: the out put filename;
#      RETURNS: 1:0
#  DESCRIPTION: given a hash of bam, write down the id to a file;
#		  Note: 
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub write_id{
	my $bam_rf = shift;
	my $filename = shift;
	open my $out_fh, ">", $filename or die "Cannot create the file:$!";
	print $out_fh "id\n";
	foreach my $id (sort keys %{$bam_rf}){
		print $out_fh "$id\n";
	}
	close $out_fh;
}

#===  FUNCTION  ================================================================
#         NAME: plot2venn
#      PURPOSE: given a hash of bam, write down the id the to a file;
#   PARAMETERS: $bam_rf: a hash reference;
#   			$filename: the out put filename;
#      RETURNS: 1:0
#  DESCRIPTION: given a hash of bam, write down the id to a file;
#		  Note: 
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub plot2venn {
	use Statistics::R;
	my ($label1, $file_1, $label2, $file_2, $out_file) = @_;
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
	my $R = Statistics::R->new();
	$R->run($R_code);
	$R->stop;
	
}
#===  FUNCTION  ================================================================
#         NAME: plot3venn
#      PURPOSE: given a hash of bam, write down the id the to a file;
#   PARAMETERS: $bam_rf: a hash reference;
#   			$filename: the out put filename;
#      RETURNS: 1:0
#  DESCRIPTION: given a hash of bam, write down the id to a file;
#		  Note: 
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub plot3venn {
	use Statistics::R;
	my ($label1, $file_1, $label2, $file_2, $label3, $file_3, $out_file) = @_;
	my $R_code = <<RCODE;
library(VennDiagram)

$label1 <- read.table(file="$file_1", header=T)
$label2 <- read.table(file="$file_2", header=T)
$label3 <- read.table(file="$file_3", header=T)

data <- list('$label1' = $label1\$id,
  '$label2' = $label2\$id,
  '$label3' = $label3\$id,
  )
venn.plot <- venn.diagram(
	x = data,
	filename = "$out_file",
  # category.names = c(
    # expression( bold('A'['1: subscript']) ),
    # expression( bold('B'^'2: going up') ),
    # expression( paste(bold('C'^'3'), bold('X'['i' <= 'r'^'2']^'2') ) )
  # ),
  output = TRUE,
  height = 3000,
  width = 3000,
  resolution = 300,
  compression = 'lzw',
  units = 'px',
  lwd = 6,
  lty = 'blank',
  fill = c('yellow', 'purple', 'green'),
  cex = 3.5,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 3,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
);

RCODE
	my $R = Statistics::R->new();
	$R->run($R_code);
	$R->stop;
	
}

#===  FUNCTION  ================================================================
#         NAME: plot4venn
#      PURPOSE: given a hash of bam, write down the id the to a file;
#   PARAMETERS: $bam_rf: a hash reference;
#   			$filename: the out put filename;
#      RETURNS: 1:0
#  DESCRIPTION: given a hash of bam, write down the id to a file;
#		  Note: 
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub plot4venn {
	use Statistics::R;
	my ($label1, $file_1, $label2, $file_2, $label3, $file_3, $label4, $file_4, $out_file) = @_;
	my $R_code = <<RCODE;
library(VennDiagram)

$label1 <- read.table(file="$file_1", header=T)
$label2 <- read.table(file="$file_2", header=T)
$label3 <- read.table(file="$file_3", header=T)
$label4 <- read.table(file="$file_4", header=T)

data <- list('$label1' = $label1\$id,
  '$label2' = $label2\$id,
  '$label3' = $label3\$id,
  '$label4' = $label4\$id
  )
venn.plot <- venn.diagram(
	x = data,
	filename = "$out_file",
	col = "transparent",
	fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
	alpha = 0.50,
	label.col = c("orange", "white", "darkorchid4", "white", 
	"white", "white", "white", "white", "darkblue", "white", 
	"white", "white", "white", "darkgreen", "white"),
	cex = 1.5,
	fontfamily = "serif",
	fontface = "bold",
	cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
	cat.cex = 1.5,
	cat.pos = 0,
	cat.dist = 0.07,
	cat.fontfamily = "serif",
	#rotation.degree = 270,
	margin = 0.2
	);

RCODE
	my $R = Statistics::R->new();
	$R->run($R_code);
	$R->stop;
	
}


#===  FUNCTION  ================================================================
#         NAME: plot4venn
#      PURPOSE: given a hash of bam, write down the id the to a file;
#   PARAMETERS: $bam_rf: a hash reference;
#   			$filename: the out put filename;
#      RETURNS: 1:0
#  DESCRIPTION: given a hash of bam, write down the id to a file;
#		  Note: 
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub plot5venn {
	use Statistics::R;
	my ($label1, $file_1, $label2, $file_2, $label3, $file_3, $label4, $file_4, $label5, $file_5,$out_file) = @_;
	my $R_code = <<RCODE;
library(VennDiagram)

$label1 <- read.table(file="$file_1", header=T)
$label2 <- read.table(file="$file_2", header=T)
$label3 <- read.table(file="$file_3", header=T)
$label4 <- read.table(file="$file_4", header=T)
$label5 <- read.table(file="$file_5", header=T)

data <- list('$label1' = $label1\$id,
  '$label2' = $label2\$id,
  '$label3' = $label3\$id,
  '$label4' = $label4\$id
  '$label5' = $label5\$id
  )
venn.plot <- venn.diagram(
	x = data,
	filename = "$out_file",
	col = "transparent",
	fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
	alpha = 0.50,
	label.col = c("orange", "white", "darkorchid4", "white", 
	"white", "white", "white", "white", "darkblue", "white", 
	"white", "white", "white", "darkgreen", "white"),
	cex = 1.5,
	fontfamily = "serif",
	fontface = "bold",
	cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
	cat.cex = 1.5,
	cat.pos = 0,
	cat.dist = 0.07,
	cat.fontfamily = "serif",
	#rotation.degree = 270,
	margin = 0.2
	);

RCODE
	my $R = Statistics::R->new();
	$R->run($R_code);
	$R->stop;
	
}