#!/usr/bin/perl -w

use strict;
use 5.010;
use Bio::SeqIO;
use Getopt::Long;
use File::Basename;
use Statistics::R;

my $usage = <<USAGE;
SYSNOPSIS
fasta_len_statistcs.pl -L annotated,assembled annotated.fa assembled.fa

 Options:
   -L| --label    Label for each transcriptome
   -o|--output    output folder for removed reads 
USAGE
my $out_folder = dirname './';
my $label;
die $usage
  unless GetOptions(
    "o|output=s"   => \$out_folder,
    "L|label=s"   => \$label,
  );
my @files = @ARGV;
foreach my $file(@files) {
	warn "$file does not exits" unless -e $file
}
die $usage unless @files;
my @suffx  = qw (.fa .fas .fasta);

my @labels;
if ($label) {
	@labels = split /,/, $label;
	die "labels number is not equal to transcriptome number\n"
		if @labels != @files;  	
}else{
	foreach my $i(0 .. $#files) {
		$labels[$i] = basename( $files[$i], @suffx )
	}
}



$out_folder =~ s/[\/|\|]+$//;
mkdir $out_folder unless -e $out_folder;

open my $plot_fh, ">", "$out_folder/tatal.len.txt" 
	or die "cannot create $out_folder/tatal.txt file:$!\n";
print $plot_fh "lib\tlen\n";
foreach my $i (0 .. $#files) {
	my $label = $labels[$i];
	my $in_fasta = $files[$i];
	my $out_file = "$out_folder/$label.len.txt";
	open my $out_fh, ">", $out_file or die "cannot create files:\n";
	print $out_fh "id\tlen\n";
	my $in  = Bio::SeqIO->new(-file => $in_fasta ,
						 -format => 'Fasta');
						 
	while ( my $seq = $in->next_seq() ) {
		my $id = $seq->id();
		my $len = $seq->length();
		print $out_fh "$id\t$len\n";
		print $plot_fh "$label\t$len\n";
	}
	close $out_fh;
}
close $plot_fh;

my $R_code = plot_CODE();
my $R = Statistics::R->new();
my $out = $R->run($R_code);
print $out,"\n";
$R->stop;
	


sub plot_CODE {
    return <<MARKER;
data <- read.table(file="$out_folder/tatal.len.txt", header=T)
head(data)
library(ggplot2)
library(plyr)
median <- ddply(data, .(lib), summarise, median = median(len))
head(median)
pdf(file="$out_folder/boxplot.pdf", width=4, height=3)
ggplot(data = data, aes(y = len, x = lib, fill = lib)) +
  stat_boxplot(geom ='errorbar')+ 
  geom_boxplot(outlier.shape = NA)+
  coord_cartesian(ylim = quantile(data\$len, c(0.05, 0.95)))+
  geom_text(data = median, aes(x = lib, y = median, label = median), 
             size = 3, vjust = -1, color = "red") +
  ylab("Length of introns") + 
  xlab(NULL)+
  theme_bw() +
  theme(
    legend.position = "bottom",
    plot.background = element_blank(),
    # axis.line = element_line(color = 'black'),
    panel.border = element_rect(linetype = "solid", color="black"),
    # panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
dev.off()
MARKER

}
