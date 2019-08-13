#!/usr/bin/perl -w

use strict;
use 5.010;
use Getopt::Long;
use File::Basename;
use Statistics::R;

my $usage= <<USAGE;
SYSNOPSIS
venn_diagram.pl [options] -o out.tiff -L label1,label2 1.txt 2.txt

 Options:
   -a --header        the file have header, default FALSE
   -o --outfile       VennDiagram fileanme
   -L --label         Label for each transcriptome
   -f --format        'pdf', 'bmp', 'jpeg', 'png', default, guess from output 
                      filename specified, if not 'pdf'
   -w|--width         width of output figure, default 10
   -h|--height        heigth of output figure, default 10
USAGE

my $header_flag = 0;
my $out_file = 'VennDiagram.tiff';
my $label;
my $width = 10;
my $height = 10;
my $format = '';
die $usage 
	unless GetOptions (
			"a|header" => \$header_flag,
			"o|outfile:s" => \$out_file,
			"l|label=s" => \$label,
			"f|format=s" => \$format,
			"w|width=f"    => \$width,
			"h|heigh=f"    => \$height,
			);
my @files = @ARGV;
die $usage unless @files;
my @labels = split /,/, $label;
die "labels number is not equal to transcriptome number\n"
	if @labels != @files;
#get the venndiagram format
if ($format) {
	if ($format =~ /(pdf)|(bmp)|(jpeg)|(png)/) {
		$format = $1;
		$out_file .= ".$format";
	}else{
		die "unsupported format:$format\n";
	}
}else{
	my @suffixlist = qw(.pdf .bmp .jpeg .png);
	my ($name,$path,$suffix) = fileparse($out_file,@suffixlist);
	if ($suffix) {
		$suffix =~ s/\.//;
		$format = $suffix;
	}else{
		warn 'cannot guess from output file, using pdf...\n';
		$format = 'pdf';
		$out_file .= ".$format";
	}
}
$header_flag = $header_flag ? "T" : "F";

my @plot_venn_para = ();
for(my $i = 0; $i < @labels; $i++) {
	push @plot_venn_para, ($labels[$i], $files[$i]);
}
	
if (@files == 1){
	die "Only one id file\n";
}elsif (@files == 2){
	plot2venn(@plot_venn_para, $header_flag, $out_file);
}elsif (@files ==3) {
	plot3venn(@plot_venn_para, $header_flag, $out_file);	
}elsif (@files == 4) {	
	plot4venn(@plot_venn_para, $header_flag, $out_file);
}elsif (@files == 5) {	
	plot5venn(@plot_venn_para, $header_flag, $out_file);
}
else{
	print "This program only support two set Venndiagram";
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
	my ($label1, $file_1, $label2, $file_2, $header_flag, $out_file) = @_;
	my $R_code = <<RCODE;
library(VennDiagram)

a <- read.table(file="$file_1", header=$header_flag)
b <- read.table(file="$file_2", header=$header_flag)

data <- list('$label1' = a[, 1],
  '$label2' = b[, 1]
  )
venn.plot <- venn.diagram(
	x = data,
	filename = NULL,
	output = TRUE,
	height = 3000,
	width = 3000,
	resolution = 300,
	compression = 'lzw',
	units = 'px',
	lwd = 6,
	col = c( '#a73727', '#346c95'),
	cat.col = c('#a73727', '#346c95'),
	label.col = "#404040",
	alpha = 0.75,
	cex = 3.5,
	fontfamily = "serif",
	fontface = "bold",
	cat.cex = 3,
	cat.fontfamily = "serif",
	cat.fontface = "bold",
	cat.dist = c(0.03, 0.03),
	cat.pos = c(-20, 14)
	)
$format(file='$out_file', width=$width, height=$height)
grid.draw(venn.plot)
dev.off()
RCODE
	my $R = Statistics::R->new();
	my $out = $R->run($R_code);
	print $out,"\n";
	$R->stop;
	
}

#===  FUNCTION  ================================================================
#         NAME: plot3venn
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
sub plot3venn {
	use Statistics::R;
	my ($label1, $file_1, $label2, $file_2, $label3, $file_3, $header_flag, $out_file) = @_;
	my $R_code = <<RCODE;
library(VennDiagram)

a <- read.table(file="$file_1", header=$header_flag)
b <- read.table(file="$file_2", header=$header_flag)
c <- read.table(file="$file_3", header=$header_flag)

data <- list('$label1' = a[, 1],
  '$label2' = b[, 1],
  '$label3' = c[, 1]
  )
venn.plot <- venn.diagram(
	x = data,
	filename = NULL,
	output = TRUE,
	height = 3000,
	width = 3000,
	resolution = 300,
	compression = 'lzw',
	units = 'px',
	# lty = 'blank',
	lwd = 6,
	col = c('#2d9f7d' , '#a73727', '#346c95'),
	# fill = c('#346c95', '#a73727', '#f9c924'),
	cat.col = c('#2d9f7d' , '#a73727', '#346c95'),
	label.col = "#404040",
	cex = 3.5,
	fontfamily = "serif",
	fontface = "bold",
	cat.cex = 3,
	cat.fontfamily = "serif",
	cat.default.pos = "outer",
	cat.fontface = "bold",
	cat.dist = c(0.055, 0.055, 0.085),
	rotation = 1,
	cat.pos = c(-27, 27, 135)
	)
$format(file='$out_file', width=$width, height=$height)
grid.draw(venn.plot)
dev.off()
RCODE
	my $R = Statistics::R->new();
	my $out = $R->run($R_code);
	print $out,"\n";
	$R->stop;
	
}
#===  FUNCTION  ================================================================
#         NAME: plot4venn
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
sub plot4venn {
	use Statistics::R;
	my ($label1, $file_1, $label2, $file_2, $label3, $file_3, $label4, $file_4, $header_flag, $out_file) = @_;
	my $R_code = <<RCODE;
library(VennDiagram)

a <- read.table(file="$file_1", header=$header_flag)
b <- read.table(file="$file_2", header=$header_flag)
c <- read.table(file="$file_3", header=$header_flag)
d <- read.table(file="$file_4", header=$header_flag)

data <- list('$label1' = a[, 1],
  '$label2' = b[, 1],
  '$label3' = c[, 1],
  '$label4' = d[, 1]
  )
venn.plot <- venn.diagram(
	x = data,
	filename = NULL,
	#rotation.degree = 270,
	# lty = "dotted",
	lwd = 4,
	col = c('#4472C4' , '#FFA032', '#BF3EFF', '#008000'),
	# cat.col = c('#4472C4' , '#FFA032', '#BF3EFF', '#008000'),
	label.col = "#404040",
	cex = 2.5,
	fontfamily = "serif",
	fontface = "bold",
	cat.cex = 2.5,
	cat.fontfamily = "serif"
	)
$format(file='$out_file', width=$width, height=$height)
grid.draw(venn.plot)
dev.off()
RCODE
	my $R = Statistics::R->new();
	my $out = $R->run($R_code);
	print $out,"\n";
	$R->stop;
	
}

#===  FUNCTION  ================================================================
#         NAME: plot5venn
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
sub plot5venn {
	use Statistics::R;
	my ($label1, $file_1, $label2, $file_2, $label3, $file_3, $label4, $file_4, $label5, $file_5, $header_flag, $out_file) = @_;
	my $R_code = <<RCODE;
library(VennDiagram)

a <- read.table(file="$file_1", header=$header_flag)
b <- read.table(file="$file_2", header=$header_flag)
c <- read.table(file="$file_3", header=$header_flag)
d <- read.table(file="$file_4", header=$header_flag)
e <- read.table(file="$file_5", header=$header_flag)

data <- list('$label1' = a[, 1],
  '$label2' = b[, 1],
  '$label3' = c[, 1],
  '$label4' = d[, 1],
  '$label5' = e[, 1]
  )
venn.plot <- venn.diagram(
	x = data,
	filename = NULL,
	#rotation.degree = 270,
	# lty = "dotted",
	lwd = 4,
	col = c('#4472C4' , '#FFA032', '#BF3EFF', '#008000', '808080'),
	# cat.col = c('#4472C4' , '#FFA032', '#BF3EFF', '#008000'),
	label.col = "#404040",
	cex = 2.5,
	fontfamily = "serif",
	fontface = "bold",
	cat.cex = 2.5,
	cat.fontfamily = "serif"
	)
$format(file='$out_file', width=$width, height=$height)
grid.draw(venn.plot)
dev.off()
RCODE
	my $R = Statistics::R->new();
	my $out = $R->run($R_code);
	print $out,"\n";
	$R->stop;
	
}