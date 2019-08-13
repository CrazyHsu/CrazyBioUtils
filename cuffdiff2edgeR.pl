#!/usr/bin/perl -w

use strict;
use 5.010;
use Cwd 'abs_path';
use Getopt::Long;
use Statistics::R;

my $usage = <<USAGE;

cuffdiff2edgeR V1.0, written by corephi
This program is used to fetch the read count counted by cuffdiff,and store
in a txt file that edgeR can be recognised.
----------------------------------------------------------
Usage: cuffdiff2edgeR.pl [options] [-i read_group_tracking] [GroupName1 GroupName2 ...]

Options:
  -i|--in-folder    cuffdiff's output folder, default "cuffdiff"
  -l|--level        'gene', 'isoform', default gene
  -o|--out-foler    output directory, defalut 'edgeR'
  -c|--cutoff       the minimal averge rpkm for edgeR to cutoff
  --fdr             FDR threadshold, default 0.01
  --fc              Fold change threadshold, default 2

GroupName List is equal to cuffdiff's '-L' option, if ommited cuffdiff2DESeq will
automatically calculated one.
If all the group have only one replicate, this program will restore swich to 
edgeR's 'blind' and 'fit-only' way to estimate dispersions
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
my $out_folder    = 'edgeR';
my $max         = 10;
my $fc = 2;
my $fdr = 0.01;
my $fpkm_min = 1;
die $usage
  unless GetOptions(
    "i|in-foler:s"    => \$cuffdff_foler,
    "l|level:s"       => \$level,
    "o|out-folder:s"  => \$out_folder,
    "max:i"    => \$max,
    "min:f"    => \$fpkm_min,
	"fdr:f"    => \$fdr,
	"fc:f"         => \$fc,
  );

#setting the scale
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

###########################################################################
#read size information
###########################################################################
my %length;
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

#print the header
#output total files
print "outputing total fpkm...\n";
output_counts( "fpkm", @groups );
print "outputing total read count...\n";
output_counts( "count", @groups );
my $group_num    = @groups;

#output comparation
print "Starting edgeR...\n";
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

        my $group1_reps = keys %{ $gene_counts{$group1} };
        my $group2_reps = keys %{ $gene_counts{$group2} };
        if ( $group1_reps == 1 && $group2_reps == 1 ) {
			warn "This experiment has no replicate, edgeR is not reconmened";
            edgeR_without_replicate();
        }
        else {
            edgeR_with_replicates($group1, $group1_reps, $group2, $group2_reps);
        }
        chdir "../" or die "cannot change file to output direcotry";
    }
}

chdir "../" or die "cannot change file to main direcotry";
print "Done!";

sub output_counts {
    my $mode   = shift;
    my @groups = @_;
    open my $out_fh, ">", "total.$mode";
    my @headers = qw (GID);
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
	push @headers, "Length";
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
		push @contents, $length{$id};
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


sub edgeR_without_replicate {
	my ($group1_name, $group1_reps, $group2_name, $group2_reps) = @_;
	my $group_num = $group1_reps + $group2_reps;
	my $len_col = $group_num + 1;
    my $Rocode = <<edgeR;
library(edgeR)

edgeR

    open my $log_fh, ">", "edgeR.log" or die "cannot create data file\n";
    my $R     = Statistics::R->new();
    my $R_out = $R->run($Rocode);
    print $log_fh $R_out, "\n";
    $R->stop;
}

sub edgeR_with_replicates {
	my ($group1_name, $group1_reps, $group2_name, $group2_reps) = @_;
	my $group_num = $group1_reps + $group2_reps;
	my $len_col = $group_num + 1;
    my $Rocode = <<edgeR;
library("edgeR")
library(ggplot2)
library(VennDiagram)
sessionInfo()
#read file
rawdata <- read.delim("total.count",row.names=1, stringsAsFactors=FALSE)
head(rawdata)
group <- c(rep ("$group1_name", $group1_reps), rep("$group2_name",$group2_reps))
data <- DGEList(counts=rawdata[,1:$group_num], group = group, genes = data.frame(Length=rawdata[,$len_col], row.names = rownames(rawdata)
))
head(data\$genes) 

#filter by rpkm
keep <- rowSums(rpkm(data, data\$genes\$Length, normalized.lib.sizes=TRUE ) > $fpkm_min) > 2
data <- data[keep,]

#recalculate
data\$samples\$lib.size <- colSums(data\$counts)


#Normalizing
data <- calcNormFactors(data)
data\$samples

#calculat rpkm
rpkm <- rpkm(data, data\$genes\$Length, normalized.lib.sizes=TRUE )
head(rpkm)

#PlotMDS
pdf(file='MDS.pdf')
plotMDS(data, col=c(rep("#00BFC4",$group1_reps), rep("#F8766D",$group2_reps)  ))
dev.off()

data <- estimateCommonDisp(data, verbose=TRUE)
data <- estimateTagwiseDisp(data)
plotBCV(data)

#WT, 
et <- exactTest(data, pair = c("$group1_name","$group2_name"))
head(et)

summary(de <- decideTestsDGE(et, p=0.01, adjust.method="BH")) 
tp <- topTags(et,n=200000, adjust.method="BH", sort.by="logFC")

#Plot shreshold
gene_list = as.data.frame(tp)
head(gene_list)
multiple = $fc
pvalue = $fdr
adj.pvalue = $fdr
gene_list\$pvalue_sgnt = as.factor(abs(gene_list\$logFC) > log2(multiple ) & gene_list\$PValue <= pvalue)
gene_list\$adj_pvalue_sgnt = as.factor(abs(gene_list\$logFC) > log2(multiple ) & gene_list\$FDR <= adj.pvalue)

#maplot
pdf(file='MAPlot_pvalue.pdf')
g = ggplot(data=gene_list, aes(x=logCPM, y=logFC, shape=pvalue_sgnt, colour=pvalue_sgnt)) +
  geom_point(alpha=0.6, size=1.2) +
  scale_colour_manual(values=c("#00BFC4", "#F8766D")) + 
  #geom_hline(yintercept=-log10(c(0.05,0.01,0.005,0.001)), color="purple") +
  #geom_vline(xintercept=c(-2,-1.5,-1,1,1.5,2), color="blue") +
  xlim(c(-4, 14)) + ylim(c(-5, 5)) +
  labs(title = "MAPlot of $group2_name vs $group1_name",
  x = "log2 average cpm",
  y= "log2 fold change ($group2_name / $group1_name)") +
  theme_bw() +
  theme(legend.position = "none")
g 
dev.off()

#volcanoplot
pdf(file='VolcanoPlot_pvalue.pdf')
g = ggplot(data=gene_list, aes(x=logFC, y=-log10(PValue), shape=pvalue_sgnt, colour=pvalue_sgnt)) +
  geom_point(alpha=0.6, size=1.2) +
  scale_colour_manual(values=c("#00BFC4", "#F8766D")) + 
  #geom_hline(yintercept=-log10(c(0.05,0.01,0.005,0.001)), color="purple") +
  #geom_vline(xintercept=c(-2,-1.5,-1,1,1.5,2), color="blue") +
  xlim(c(-5, 5)) + ylim(c(0, 7)) +
  labs(title = "Volcano Plot of $group2_name vs $group1_name",
  x = "log2 fold change ($group2_name / $group1_name)",
  y= "-log10 Pvalue") +
  theme_bw() +
  theme(legend.position = "none")
g 
dev.off()

#maplot
pdf(file='MAPlot_adjpvalue.pdf')
g = ggplot(data=gene_list, aes(x=logCPM, y=logFC, shape=adj_pvalue_sgnt, colour=adj_pvalue_sgnt)) +
  geom_point(alpha=0.6, size=1.2) +
  scale_colour_manual(values=c("#00BFC4", "#F8766D")) + 
  #geom_hline(yintercept=-log10(c(0.05,0.01,0.005,0.001)), color="purple") +
  #geom_vline(xintercept=c(-2,-1.5,-1,1,1.5,2), color="blue") +
  xlim(c(-4, 14)) + ylim(c(-5, 5)) +
  labs(title = "MAPlot of $group2_name vs $group1_name",
  x = "log2 average cpm",
  y= "log2 fold change ($group2_name / $group1_name)") +
  theme_bw() +
  theme(legend.position = "none")
g 
dev.off()

#volcanoplot
pdf(file='VolcanoPlot_adjpvalue.pdf')
g = ggplot(data=gene_list, aes(x=logFC, y=-log10(FDR), shape=adj_pvalue_sgnt, colour=adj_pvalue_sgnt)) +
  geom_point(alpha=0.6, size=1.2) +
  scale_colour_manual(values=c("#00BFC4", "#F8766D")) + 
  #geom_hline(yintercept=-log10(c(0.05,0.01,0.005,0.001)), color="purple") +
  #geom_vline(xintercept=c(-2,-1.5,-1,1,1.5,2), color="blue") +
  xlim(c(-5, 5)) + ylim(c(0, 5)) +
  labs(title = "Volcano Plot of $group2_name vs $group1_name",
  x = "log2 fold change ($group2_name / $group1_name)",
  y= "-log10 Pvalue") +
  theme_bw() +
  theme(legend.position = "none")
g 
dev.off()

#write to file
write.table(gene_list,file='Diff.xls',sep="\t", quote=F, row.names = TRUE, col.names = TRUE)
write.table(rpkm,file='Diff.rpkm.xls',sep="\t", quote=F, row.names = TRUE, col.names = TRUE)
edgeR

    open my $log_fh, ">", "edgeR.log" or die "cannot create data file\n";
    my $R     = Statistics::R->new();
    my $R_out = $R->run($Rocode);
    print $log_fh $Rocode, "\n";
    print $log_fh $R_out, "\n";
    close $log_fh;
    $R->stop;
}

