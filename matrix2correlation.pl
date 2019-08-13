#!/usr/bin/perl -w

use strict;
use 5.010;
use Cwd 'abs_path';
use Getopt::Long;
use File::Basename;
use Statistics::R;

my $usage = <<USAGE;

matrix2correlation V1.0, written by corephi
----------------------------------------------------------
Usage: matrix2correlation.pl [options] rpkm.matrix1 rpkm.matrix2

Options:
  -t|top            top number gene used for plot, 0 means all. defualt 0
  -m|--marker-list  marker gene list, it must contain two colums, hand
                    has a header of "GID	SymbolName"
  -e|--extremum     fpkm interval, if set to 6, the minimal fpkm is trimmed
                    2^-6, and the maximum fpkm is trimmed to 2^6, default 10
USAGE

#########################################################################
##setting the parameters
#########################################################################
my $top              = 0;
my $marker_list_file = "";
my $extremum         = 10;
die $usage
  unless GetOptions(
    "m|marker-list:s" => \$marker_list_file,
    "e|extremum:i"    => \$extremum,
    "t|top:i"    => \$top,
  );

#setting the scale
my $fpkm_min = 2**-$extremum;
my $fpkm_max = 2**$extremum;

#check the infile
my @files = @ARGV;
foreach my $file (@files) {
	die "$file does not exists\n" unless -e $file;
}

#check maker_list_file
if ($marker_list_file) {
	$marker_list_file = abs_path($marker_list_file);
	die "$marker_list_file does not exists\n" unless -e $marker_list_file;
}

#normalize
warn "Normalizing...\n";
foreach my $i (0 .. $#files) {
	my $file = $files[$i];
	
	open my $in_fh, "<", $file or die "cannot open $file:$!\n";
	my $out_file = basename($file, ".txt", ".total" ).".normalized.dat";
	open my $out_fh, ">", $out_file  or die "cannot create file:$!\n";
	my $header = readline $in_fh;
	print $out_fh $header;
	while (<$in_fh>) {
		chomp;
		my @lines = split /\t/; 
		my $id = shift @lines;
		foreach my $i(0 .. $#lines) {
			$lines[$i] = $fpkm_max if $lines[$i] > $fpkm_max;
			$lines[$i] = $fpkm_min if $lines[$i] < $fpkm_min;
		}
		print $out_fh join "\t", $id, @lines;
		print $out_fh "\n";
	}
	$files[$i] = $out_file;	
	close $in_fh;
	close $out_fh;
}	

#find the top number genes
if ($top) {
	warn "Find top number entity...\n";
	foreach my $i (0 .. $#files) {
		my $file = $files[$i];
		
		#read the average and store to hash %rpkm
		open my $in_fh, "<", $file or die "cannot open $file:$!\n";
		readline $in_fh;
		my %rpkms;
		while (<$in_fh>) {
			chomp;
			my @lines = split /\t/;
			my $id = shift @lines;
			my $total = 0;
			$total += $_ foreach @lines;
			$rpkms{$id} = $total;
		}
		
		#find the top xx genes;
		my %topn = map {$_ => $rpkms{$_} } ( sort { $rpkms{$b} <=> $rpkms{$a} } keys %rpkms ) [0..$top-1];
		
		
		#
		open $in_fh, "<", $file or die "cannot open $file:$!\n";
		my $out_file = $file.".top$top.dat";
		open my $out_fh, ">", $out_file or die "cannot create file:$!\n";
		my $header = readline $in_fh;
		print $out_fh $header;
		while (<$in_fh>) {
			my $line = $_;
			my @lines = split /\t/; 
			my $id = shift @lines;
			print $out_fh $_ if exists $topn{$id};
		}
		$files[$i] = $out_file;	
		close $in_fh;
		close $out_fh;
	}	
}else{
	
}

#############################################################################
#output the marker gene
#############################################################################
if ($marker_list_file) {
	
	#read the marker id;
    open my $marker_fh, "<", $marker_list_file
      or die "cannot open file $marker_list_file:$!";
	my $header = readline  $marker_fh;
	die "maker does not had a header\n" unless $header =~ /GID\tSymbolName/;
	my %markers = ();
	while ( <$marker_fh> ) {
		chomp;
		my ($gid, $alias) = split /\t/;
		$markers{$gid} = $alias;
	}
	
	foreach my $file (@files) {
		open my $in_fh, "<", $file or die "cannot open $file:$!\n";
		my $out_file = $file.".markers.dat";
		open my $out_fh, ">",  or die "cannot create file:$!\n";
		my $header = readline $in_fh;
		print $out_fh $header;
		while (<$in_fh>) {
			my @lines = split /\t/; 
			my $id = shift @lines;
			if (exists $markers{$id}) {
				my $alias = $markers{$id};
				unshift @lines, $alias;
				print $out_fh join "\t", @lines;
				print $out_fh "\n";
			}
		}
		close $in_fh;
		close $out_fh;			
	}
	close $marker_fh;
}



##########################################################################
#PLot Correlation
##########################################################################

foreach my $file (@files) {
	my $R_corelation = correlation_plot($file);
	my $R_Marker     = markergene_plot($file.".markers.dat");

	my $R = Statistics::R->new();
	print "Start Ploting $file...\n";
	my $R_out = $R->run($R_corelation);
	$R_out .= $R->run($R_Marker) if $marker_list_file;
	print $R_out, "\n";
	$R->stop;
}


sub correlation_plot {
	my $file = shift;
    return <<CORRELATION
require(corrplot)||{install.packages("corrplot");require(corrplot)}
require(RColorBrewer)||{install.packages("RColorBrewer");require(RColorBrewer)}
require(pheatmap)||{install.packages("pheatmap");require(pheatmap)}

print("Reading data files:$file...")
rpkm <- read.delim(file="$file",stringsAsFactors=F)
rownames(rpkm) = rpkm\$gene_id
colnum <- length(colnames(rpkm))
rpkm <- rpkm[, 2:colnum]
print("Filtering low coverage genes...")
keep <- rowSums(rpkm > 0.5) > 1
print("Data samples are:")
rpkm <- rpkm[keep,]

cor1 <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")), interpolate = "spline")
cor2 <- colorRampPalette(brewer.pal(9, "OrRd"), interpolate = "spline")
log2rpkm <- log2(rpkm+1)
head(log2rpkm)
print("Plotting correlation...")
pdf(file="$file.correlation.logfpkm.pdf")
corrplot(cor(log2rpkm),col=cor1(512),method="color",
  order="hclust",
  # cl.lim=c(-1,1),
  cl.align="c",
  tl.pos="d", addCoef.col="grey")
dev.off()
print("Done!")

log2rpkm.t <- t(log2rpkm)
print("Plotting cluster heatmap...")
pdf(file="$file.cluster.pdf", width = 8, height = 3.5)
pheatmap(log2rpkm.t,
  #color =  colorRampPalette(c("black", "red"))(512),
  color = cor2($extremum + 1),
  legend_breaks = -$extremum:$extremum,
  show_colnames = F,
  cluster_rows = T,
  cluster_cols = F)
dev.off()
print("Done!")
CORRELATION
}

sub markergene_plot {
	my $file = shift;
    return <<MARKER;
print("Plotting MarkerGenes...")
marker <- read.delim(file="$file",stringsAsFactors=F)
#marker <- log2(marker)

# my_breaks = function(x, n, center = F){
        # if(center){
                # m = max(abs(c(min(x, na.rm = T), max(x, na.rm = T))))
                # res = seq(-m, m, length.out = n + 1)
        # }
        # else{
                # res = seq(min(x, na.rm = T), max(x, na.rm = T), length.out = n + 1)
        # }
        
        # return(res)
# }
# marker.break <- my_breaks(marker, 100, center=T)

pdf(file="$file.MarkerHeatmap.pdf", width = 10, height = 10)
pheatmap(marker, 
  #breaks = marker.break,
  #legend_breaks = -$extremum:$extremum,
  scale = "column"
)
dev.off()
print("Done!")
MARKER

}
