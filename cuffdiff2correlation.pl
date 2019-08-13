#!/usr/bin/perl -w

use strict;
use 5.010;
use Cwd 'abs_path';
use Getopt::Long;
use Statistics::R;

my $usage = <<USAGE;

cuffdiff2correlation V1.0, written by corephi
This program is used to fetch the read count counted by cuffdiff,and store
in a txt file that edgeR can be recognised.
----------------------------------------------------------
More scripts? Join "bioinformatics*CN" QQ group: 276151571
If you have any probel or suggestions about this program, 
please mail to: hpyu\@genetics.ac.cn
----------------------------------------------------------

Usage: cuffdiff2correlation.pl [options] [-i read_group_tracking] [GroupName1 GroupName2 ...]

Options:
  -i|--in-folder    cuffdiff's output folder, default "cuffdiff"
  -l|--level        'gene', 'isoform', default gene
  -s|--size         output the 'gene' or 'isoform' size(length) information
  -o|--out-foler    output directory, defalut 'correlation'
  -p|--plot         plot correlation and heatmap, depends on perl module
                    "Statistics::R" and "R" in your system environment,
                    so linux is recommand. If in window, you must and R
                    to your path:
  -m|--marker-list  marker gene list, it must contain two colums, hand
                    has a header of "GID	SymbolName"
  -min              the minimal averge fpkm value, default 1 is 'plot' is 
                    setted on, 0 if 'plot' is seted off.
  -max              the maxmum fpkm value, defalut 10 (2**10)
					

GroupName List is equal to cuffdiff's '-L' option, if ommited cuffdiff2edgeR will
automatically calculated one.
Note:If the cuffdiff stats is "OK", then ommited, but store it in outlier.gene
This program also discard genes tagged as "mRNA_TE_gene" in TAIR10.gff, which 
cannot be correctly recognised by cuffdiff 
USAGE

#########################################################################
##setting the parameters
#########################################################################
my $cuffdff_foler = "cuffdiff";
my $gff_file      = '';
my $level         = 'gene';
my $size_mode     = 0;
my $out_folder    = 'correlation';

my $plot_flag        = 0;
my $marker_list_file = "";
my $max         = 10;
my $fpkm_min = 1;
die $usage
  unless GetOptions(
    "i|in-foler:s"    => \$cuffdff_foler,
    "s|size"          => \$size_mode,
    "l|level:s"       => \$level,
    "o|out-folder:s"  => \$out_folder,
    "p|plot"          => \$plot_flag,
    "m|marker-list:s" => \$marker_list_file,
    "max:i"    => \$max,
    "min:f"    => \$fpkm_min,
  );

die "size mode and plot mode can NOT be setted on together\n"
  if $size_mode && $plot_flag;
$out_folder =~ s/\/|\\$//;
#setting the scale
$fpkm_min = 0 unless $plot_flag;
my $fpkm_max = 2**$max;


#cuffdiff folder
$cuffdff_foler =~ s/[\\|\/]$//igx;
die "cuffdiff directory: $cuffdff_foler does not exists!\n"
  unless -e $cuffdff_foler;

#level
my $read_group_tracking;
my $fpkm_tracking;
if ( $level eq "isoform" ) {
    $read_group_tracking = "isoforms.read_group_tracking";
    $fpkm_tracking       = "isoforms.fpkm_tracking";
}
elsif ( $level eq "gene" ) {
    $read_group_tracking = "genes.read_group_tracking";
    $fpkm_tracking       = "genes.fpkm_tracking";
}
else {
    warn "level $level is canot recogonised, use gene instead\n";
    $read_group_tracking = "genes.read_group_tracking";
    $fpkm_tracking       = "genes.fpkm_tracking";
}
$read_group_tracking = abs_path( $cuffdff_foler . "/" . $read_group_tracking );
$fpkm_tracking       = abs_path( $cuffdff_foler . "/" . $fpkm_tracking );

#groups
my $read_groups_info_file = abs_path( $cuffdff_foler . "/read_groups.info" );
print "Staring get Group information...\n";
my @groups = @ARGV;
if (@groups) {
    print "You specifed your group";
}
else {
    print "You haven't specifed your group\n";
    my %data;
    open my $read_groups_fh, "<", $read_groups_info_file
      or die "cannot open read groups information file:$!";
    readline $read_groups_fh;
    while (<$read_groups_fh>) {
        chomp;
        next if /^\s/;
        my ( $file, $condition, $replicate_num ) = split /\t/;
        $data{$condition} = 1;
    }
    close $read_groups_fh;
    @groups = sort keys %data;
    print "Automatically calculated one: @groups\n";

}

$marker_list_file = abs_path($marker_list_file) if $marker_list_file;

###########################################################################
#read size information
###########################################################################
my %length;
if ($size_mode) {
    print "getting lengh information...\n";
    open my $length_fh, "<", $fpkm_tracking;
    readline $length_fh;
    while (<$length_fh>) {
        my (
            $tracking_id, $class_code,      $nearest_ref_id,
            $gene_id,     $gene_short_name, $tss_id,
            $locus,       $length,          $coverage
        ) = split /\t/;
        if ( $length eq '-' ) {
            my @temp = split /:/, $locus;
            my ( $start, $end ) = split /-/, $temp[1];
            $length = $end - $start + 1;
        }
        $length{$tracking_id} = $length;
    }
    close $length_fh;
}

##########################################################################
#read cuffdiff gene read counts and store the groups infomation
###########################################################################
my %outliers;
my %gene_counts;

print "reading fpkm and count infomation...\n";
open my $read_group_tracking_fh, "<", $read_group_tracking
  or die "cannot open file:$read_group_tracking";
readline $read_group_tracking_fh;
my %fpkm_summary;
while (<$read_group_tracking_fh>) {
    chomp;
    next if /^\s/;
    my ( $gene_id, $groups, $replicate, $raw_flags, $internal_scaled_frags,
        $external_scaled_frags, $fpkm, $effective_length, $status )
      = split /\t/;
    $fpkm = $fpkm_max if $fpkm > $fpkm_max;
    if ( $status eq 'OK' ) {
        $gene_counts{$groups}{$replicate}{$gene_id}{count} = sprintf "%.0f",
          $raw_flags;
        $gene_counts{$groups}{$replicate}{$gene_id}{fpkm} = $fpkm;
		$fpkm_summary{$gene_id}{"${groups}_R$replicate"} = $fpkm;
    }
    else {
        $outliers{$gene_id} = $status;
        $gene_counts{$groups}{$replicate}{$gene_id}{count} = sprintf "%.0f",
          $raw_flags;
        $gene_counts{$groups}{$replicate}{$gene_id}{fpkm} = $fpkm;
    }
}
close $read_group_tracking_fh;

##########################################################################
#get the average value
##########################################################################
my %average;
foreach my $gid (keys %fpkm_summary) {
	my @samples = keys %{$fpkm_summary{$gid}};
	my $total = 0;
	foreach my $sample (@samples) {
		my $rpkm = $fpkm_summary{$gid}{$sample};
		$total += $rpkm;
	}
	my $average = $total / @samples;
	$average{$gid} = $average;
} 
##########################################################################

##########################################################################
#output the results
##########################################################################
#make edgeR directory

if ( -e $out_folder ) {
    die "directory $out_folder is not writable " unless -w "$out_folder";
}
else {
    mkdir $out_folder, 0775 or die "cannot make logs directory:$!";
}

chdir $out_folder or die "cannot change file to $out_folder";

#############################################################################
#output the marker gene
#############################################################################
if ($marker_list_file) {
    open my $marker_fh, "<", $marker_list_file
      or die "cannot open file $marker_list_file:$!";
    open my $marker_fpkm_fh, ">", "MarkerGeneFPKM.txt"
      or die "cannot open file MarkerGeneFPKM.txt:$!";

    #construct header
    my @marker_headers = ();
    foreach my $group (@groups) {
        foreach my $replicate ( sort keys %{ $gene_counts{$group} } ) {
            $replicate++;
            push @marker_headers, "${group}_R$replicate";
        }
    }
    my $header = join "\t", @marker_headers;
    print $marker_fpkm_fh $header, "\n";

    #output the contents
    #read markers
    readline $marker_fh;
    while (<$marker_fh>) {
        chomp;
        my ( $id, $symbol ) = split /\t/;
        my @contents = ($symbol);

        #get the counts infomation
        foreach my $group (@groups) {
            foreach my $replicate ( sort keys %{ $gene_counts{$group} } ) {
                my $fpkm = $gene_counts{$group}{$replicate}{$id}{fpkm};
                $fpkm = '0' unless $fpkm;
                $fpkm = $fpkm_max if $fpkm > $fpkm_max;
                push @contents, $fpkm;
            }
        }
        my $content = join "\t", @contents;
        print $marker_fpkm_fh $content, "\n";
    }
    close $marker_fh;
    close $marker_fpkm_fh;

}

#print the header
#output files
print "outputing fpkm...\n";
output_file("fpkm");
print "outputing read count...\n";
output_file("count");

chdir "../" or die "cannot change file to main direcotry";

##########################################################################
#PLot Correlation
##########################################################################
my $group_num    = @groups;
my $R_corelation = correlation_plot();
my $R_Marker     = markergene_plot();

if ($plot_flag) {
    my $R = Statistics::R->new();
    print "Start Ploting...\n";
    my $R_out = $R->run($R_corelation);
    $R_out .= $R->run($R_Marker) if $marker_list_file;
    print $R_out, "\n";
    $R->stop;
}

sub correlation_plot {
    return <<CORRELATION
require(corrplot)||{install.packages("corrplot");require(corrplot)}
require(RColorBrewer)||{install.packages("RColorBrewer");require(RColorBrewer)}
require(pheatmap)||{install.packages("pheatmap");require(pheatmap)}

print("Reading data files:$out_folder/total.fpkm...")
rpkm <- read.delim(file="$out_folder/total.fpkm",stringsAsFactors=F)
rownames(rpkm) = rpkm\$gene_id
colnum <- length(colnames(rpkm))
rpkm <- rpkm[, 2:colnum]

cor1 <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")), interpolate = "spline")
cor2 <- colorRampPalette(brewer.pal(9, "OrRd"), interpolate = "spline")
log2rpkm <- log2(rpkm+1)
head(log2rpkm)
print("Plotting correlation...")
pdf(file="$out_folder/correlation_logfpkm.pdf")
corrplot(cor(log2rpkm),col=cor1(512),method="color",
  order="hclust",addrect=$group_num,rect.col='green',
  cl.lim=c(0,1), cl.align="c",
  tl.pos="d", addCoef.col="grey")
dev.off()
print("Done!")

log2rpkm.t <- t(log2rpkm)
print("Plotting cluster heatmap...")
pdf(file="$out_folder/cluster.pdf", width = 6, height = 6)
pheatmap(log2rpkm.t,
  #color =  colorRampPalette(c("black", "red"))(512),
  color = cor1($max),
  legend_breaks = -$max:$max,
  show_colnames = F,
  cluster_rows = T,
  cluster_cols = T,
  treeheight_col= 0)
dev.off()
print("Done!")
CORRELATION
}

sub markergene_plot {
    return <<MARKER;
print("Plotting MarkerGenes...")
marker <- read.delim(file="$out_folder/MarkerGeneFPKM.txt",stringsAsFactors=F)
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

pdf(file="MarkerHeatmap.pdf", width = 10, height = 10)
pheatmap(marker, 
  #breaks = marker.break,
  #legend_breaks = -$max:$max,
  scale = "column"
)
dev.off()
print("Done!")
MARKER

}

sub output_file {
    my $mode = shift;
    open my $out_fh, ">", "total.$mode";
    open my $outlier_fh, ">", "outlier.txt";
    my @headers = qw (ID);
    my @ids;
    foreach my $group (@groups) {
        foreach my $replicate ( sort keys %{ $gene_counts{$group} } ) {
            @ids = sort keys %{ $gene_counts{$group}{$replicate} };
            $replicate++;
            push @headers, "${group}_R$replicate";
        }
    }
    push @headers, "length" if $size_mode;
    my $header = join "\t", @headers;
    print $out_fh $header, "\n";
	
	#print outliers
    foreach my $gene ( keys %outliers ) {
        print $outlier_fh "$gene\t$outliers{$gene}\n";
    }


    #print the content
    for my $id (@ids) {
        my @contents = ($id);

        #get the counts infomation
        foreach my $group (@groups) {
            foreach my $replicate ( sort keys %{ $gene_counts{$group} } ) {
                push @contents, $gene_counts{$group}{$replicate}{$id}{$mode};
            }
        }

        #get the length information
        if ($size_mode) {
            push @contents, $length{$id};
        }

        my $content = join "\t", @contents;
		if ($average{$id} >= $fpkm_min) {
			print $out_fh $content, "\n";
		}else{
			 print $outlier_fh "$id\n";
		}
    }
    close $out_fh;
    close $outlier_fh;	
}
