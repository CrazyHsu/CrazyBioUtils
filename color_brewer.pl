#!/usr/bin/perl -w

use strict;
use 5.010;
use Getopt::Long;
use File::Basename;
use Statistics::R;

my $usage= <<USAGE;
SYSNOPSIS
color_brewer.pl [options] -o out.tiff 1.txt

 Options:
   -a --header        the file have header, default FALSE
   -o --outfile       VennDiagram fileanme
   -m --max           use the readjust the max and min value
   -c --color         'RdYlBl', 'RdYlGn' or 'YlOrRd', default 'RdYlBl'
   -w|--width         width of output figure, default 10
   -h|--height        heigth of output figure, default 10
USAGE

my $header_flag = 0;
my $out_file = 'color.pdf';
my $max = 0;
my $width = 6;
my $height = 6;
my $color = 'RdYlBu';
die $usage 
	unless GetOptions (
			"a|header" => \$header_flag,
			"o|outfile:s" => \$out_file,
			"m|max=f" => \$max,
			"c|color=s" => \$color,
			"w|width=f"    => \$width,
			"h|heigh=f"    => \$height,
			);
$header_flag = $header_flag ? "T" : "F";
my $in_file = shift;
die "$in_file does not exists\n" unless -e $in_file;

open my $in_fh, "<", $in_file 
	or die "cannot open $in_file\n";
my $temp_file =  $in_file.".tmp";
open my $out_fh, ">", $temp_file
	or die "cannot create $temp_file\n";
if ($header_flag) {
	my $header = readline $in_fh;
	print $out_fh $header;
}
while (<$in_fh>) {
	chomp;
	my @lines = split /\t/;
	for (0 ..$#lines) {
		if ($max) {
			$lines[$_] = $max if  $lines[$_] > $max;
			$lines[$_] = -$max if  $lines[$_] < -$max;
		}		
	}
	print $out_fh join "\t", @lines, "\n";
}
close $in_fh;
close $out_fh;


my $R_code = <<RCODE;
library(pheatmap)
library(RColorBrewer)
data <- read.table(file="$temp_file", header=$header_flag)
my_breaks = function(x, n, center = F){
        if(center){
                m = max(abs(c(min(x, na.rm = T), max(x, na.rm = T))))
                res = seq(-m, m, length.out = n + 1)
        }
        else{
                res = seq(min(x, na.rm = T), max(x, na.rm = T), length.out = n + 1)
        }
        
        return(res)
}
#data.break <- my_breaks(data, 100, center=T)
pheatmap(data,
  color = colorRampPalette(brewer.pal(n = 9, name = "$color"))(100),
  cluster_rows =F ,cluster_cols = F,
  show_colnames = F, show_rowname = T,
#  breaks = data.break,
  filename = "$out_file",
  width = $width,
  height = $height)
RCODE
my $R = Statistics::R->new();
my $out = $R->run($R_code);
print $out,"\n";
$R->stop;
unlink $temp_file;
unlink "Rplots.pdf";



