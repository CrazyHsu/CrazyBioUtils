#!/usr/bin/perl -w

use strict;
use 5.010;
use Cwd 'abs_path';
use Getopt::Long;
use Statistics::R;

my $usage = <<USAGE;

cuffdiff2DESeq V1.0, written by corephi
This program is used to fetch the read count counted by cuffdiff,and store
in a txt file that edgeR can be recognised.
----------------------------------------------------------
Usage: cuffdiff2DESeq.pl [options] [-i read_group_tracking] [GroupName1 GroupName2 ...]

Options:
  -i|--in-folder    cuffdiff's output folder, default "cuffdiff"
  -l|--level        'gene', 'isoform', default gene
  -o|--out-foler    output directory, defalut 'correlation'
  -p|--plot         plot correlation and heatmap, depends on perl module
                    "Statistics::R" and "R" in your system environment,
                    so linux is recommand. If in window, you must and R
                    to your path:
  -m|--marker-list  marker gene list, it must contain two colums, hand
                    has a header of "GID	SymbolName"
  -e|--extremum     fpkm interval, if set to 6, the minimal fpkm is trimmed
                    2^-6, and the maximum fpkm is trimmed to 2^6, default 6
					

GroupName List is equal to cuffdiff's '-L' option, if ommited cuffdiff2DESeq will
automatically calculated one.
If all the group have only one replicate, this program will restore swich to 
DESeq's 'blind' and 'fit-only' way to estimate dispersions
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
my $out_folder    = 'DESeq';

my $plot_flag        = 0;
my $marker_list_file = "";
my $extremum         = 6;
die $usage
  unless GetOptions(
    "i|in-foler:s"    => \$cuffdff_foler,
    "l|level:s"       => \$level,
    "o|out-folder:s"  => \$out_folder,
    "p|plot"          => \$plot_flag,
    "m|marker-list:s" => \$marker_list_file,
    "e|extremum:i"    => \$extremum,
  );

#setting the scale
my $fpkm_min = 2**-$extremum;
my $fpkm_max = 2**$extremum;

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

##########################################################################
#read cuffdiff gene read counts and store the groups infomation
###########################################################################
my %outliers;
my %gene_counts;

print "reading fpkm and count infomation...\n";
open my $read_group_tracking_fh, "<", $read_group_tracking
  or die "cannot open file:$read_group_tracking";
readline $read_group_tracking_fh;
while (<$read_group_tracking_fh>) {
    chomp;
    next if /^\s/;
    my ( $gene_id, $groups, $replicate, $raw_flags, $internal_scaled_frags,
        $external_scaled_frags, $fpkm, $effective_length, $status )
      = split /\t/;
    $fpkm = $fpkm_min if $fpkm < $fpkm_min;
    $fpkm = $fpkm_max if $fpkm > $fpkm_max;
    if ( $status eq 'OK' ) {
        $gene_counts{$groups}{$replicate}{$gene_id}{count} = sprintf "%.0f",
          $raw_flags;
        $gene_counts{$groups}{$replicate}{$gene_id}{fpkm} = $fpkm;
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
        foreach
          my $replicate ( sort { $a <=> $b } keys %{ $gene_counts{$group} } )
        {
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
            foreach my $replicate ( sort { $a <=> $b }
                keys %{ $gene_counts{$group} } )
            {
                my $fpkm = $gene_counts{$group}{$replicate}{$id}{fpkm};
                $fpkm = $fpkm_min unless $fpkm;
                $fpkm = $fpkm_min if $fpkm < $fpkm_min;
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
#output total files
print "outputing total fpkm...\n";
output_counts( "fpkm", @groups );
print "outputing total read count...\n";
output_counts( "count", @groups );
output_design(@groups);

#output comparation
print "Starting DESeq...\n";
for ( my $i = 0 ; $i < $#groups ; $i++ ) {
    for ( my $j = 1 ; $j <= $#groups ; $j++ ) {
        next if $i == $j;
        my ( $group1, $group2 ) = @groups[ $i, $j ];
        print "\t${group2}vs${group1}...\n";
        my $subfolder = "${group2}_vs_${group1}";
        if ( -e $subfolder ) {
            warn "directory $subfolder is exists, please check the group\n";
        }
        else {
            mkdir $subfolder, 0775 or die "cannot make logs directory:$!";
        }
        chdir $subfolder or die "cannot change file to $subfolder";

        output_counts( "count", $group1, $group2 );
        output_design( $group1, $group2 );

        my $group1_reps = keys %{ $gene_counts{$group1} };
        my $group2_reps = keys %{ $gene_counts{$group2} };
        if ( $group1_reps == 1 && $group2_reps == 1 ) {
            DESeq_without_replicate();
        }
        else {
            DESeq_with_replicates();
        }
        chdir "../" or die "cannot change file to output direcotry";
    }
}

chdir "../" or die "cannot change file to main direcotry";
print "Done!";

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
print("Filtering low coverage genes...")
keep <- rowSums(rpkm > 0.5) > 1
print("Data samples are:")
rpkm <- rpkm[keep,]

cor1 <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")), interpolate = "spline")
cor2 <- colorRampPalette(brewer.pal(9, "OrRd"), interpolate = "spline")
log2rpkm <- log2(rpkm+1)
head(log2rpkm)
print("Plotting correlation...")
pdf(file="correlation_logfpkm.pdf")
corrplot(cor(log2rpkm),col=cor1(512),method="color",
  order="hclust",addrect=$group_num,rect.col='green',
  cl.lim=c(0,1), cl.align="c",
  tl.pos="d", addCoef.col="grey")
dev.off()
print("Done!")

log2rpkm.t <- t(log2rpkm)
print("Plotting cluster heatmap...")
pdf(file="cluster.pdf", width = 8, height = 3.5)
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
  #legend_breaks = -$extremum:$extremum,
  scale = "column"
)
dev.off()
print("Done!")
MARKER

}

sub output_counts {
    my $mode   = shift;
    my @groups = @_;
    open my $out_fh, ">", "total.$mode";
    my @headers = qw (ID);
    my @ids;
    foreach my $group (@groups) {
        foreach
          my $replicate ( sort { $a <=> $b } keys %{ $gene_counts{$group} } )
        {
            @ids = sort keys %{ $gene_counts{$group}{$replicate} };
            $replicate++;
            push @headers, "${group}_R$replicate";
        }
    }
    my $header = join "\t", @headers;
    print $out_fh $header, "\n";

    #print the content
    for my $id (@ids) {
        my @contents = ($id);

        #get the counts infomation
        foreach my $group (@groups) {
            foreach my $replicate ( sort { $a <=> $b }
                keys %{ $gene_counts{$group} } )
            {
                push @contents, $gene_counts{$group}{$replicate}{$id}{$mode};
            }
        }

        my $content = join "\t", @contents;
        print $out_fh $content, "\n";
    }
    close $out_fh;

    open my $outlier_fh, ">", "outlier.txt";
    foreach my $gene ( keys %outliers ) {
        print $outlier_fh "$gene\t$outliers{$gene}\n";
    }
    close $outlier_fh;
}

sub output_design {
    my @groups = @_;
    open my $out_fh, ">", "design.txt";
    print $out_fh "sample\tgroup\n";
    foreach my $group (@groups) {
        foreach
          my $replicate ( sort { $a <=> $b } keys %{ $gene_counts{$group} } )
        {
            $replicate++;
            print $out_fh "${group}_R$replicate\t${group}\n";
        }
    }
    close $out_fh;
}

sub DESeq_without_replicate {
    my $Rocode = <<DESEQ;
	library("DESeq")
	print("R version:")
	version
	print("DEseq version")
	package.version("DESeq")
	CountTable <- read.table(file="total.count",header=TRUE,row.names = 1)
	Design <- read.table(file="design.txt",header=TRUE,row.names = 1)
	conds.level <- levels(Design[, 1])
	conds <- Design[, 1]
	cds <- newCountDataSet(CountTable, conds)
	cds <- estimateSizeFactors(cds)
	pData(cds)
	cds <- estimateDispersions(cds, method="blind", sharingMode='fit-only', fitType="local")
	pdf(file="dispersion.pdf")
	plotDispEsts(cds)
	dev.off()
	condA <- conds.level[1]
	print(paste("condA:", condA))
	condB <- conds.level[2]
	print(paste("condB:", condB))
	res = nbinomTest(cds, condA, condB)
	pdf(file="MAPlot.pdf")
	plotMA(res)
	dev.off()
	write.table(res,file="diff.txt",quote=FALSE,row.names=FALSE, sep="\t")
DESEQ

    open my $log_fh, ">", "DESeq.log" or die "cannot create data file\n";
    my $R     = Statistics::R->new();
    my $R_out = $R->run($Rocode);
    print $log_fh $R_out, "\n";
    $R->stop;
}

sub DESeq_with_replicates {

    my $Rocode = <<DESEQ;
	library("DESeq")
	print("R version:")
	version
	print("DEseq version")
	package.version("DESeq")
	CountTable <- read.table(file="total.count",header=TRUE,row.names = 1)
	Design <- read.table(file="design.txt",header=TRUE,row.names = 1)
	conds.level <- levels(Design[, 1])
	conds <- Design[, 1]
	cds <- newCountDataSet(CountTable, conds)
	cds <- estimateSizeFactors(cds)
	pData(cds)
	cds<-estimateDispersions(cds,method="pooled",sharingMode="maximum",fitType="local")
	pdf(file="dispersion.pdf")
	plotDispEsts(cds)
	dev.off()
	condA <- conds.level[1]
	print(paste("condA:", condA))
	condB <- conds.level[2]
	print(paste("condB:", condB))
	res = nbinomTest(cds, condA, condB)
	pdf(file="MAPlot.pdf")
	plotMA(res)
	dev.off()
	write.table(res,file="diff.txt",quote=FALSE,row.names=FALSE, sep="\t")
DESEQ

    open my $log_fh, ">", "DESeq.log" or die "cannot create data file\n";
    my $R     = Statistics::R->new();
    my $R_out = $R->run($Rocode);
    print $log_fh $R_out, "\n";
    close $log_fh;
    $R->stop;
}

