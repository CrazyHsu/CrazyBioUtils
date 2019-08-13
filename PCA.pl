#!/usr/bin/perl -w
use strict;
use 5.010;
# use Data::Printer;
use File::Basename;
use IPC::Cmd qw[can_run run];
use Getopt::Long;
use Cwd qw(abs_path chdir getcwd);
use File::Temp;


my $usage= <<USAGE;
SYSNOPSIS
PCA.pl [options] -o out.tiff matrix.txt 

 Options:
   -o --outfile       VennDiagram fileanme
   -f --format        'pdata', 'bmp', 'jpeg', 'png', default, guess from output 
                      filename specified, if not 'pdata'
   -w|--width         width of output figure, default 10
   -h|--height        heigth of output figure, default 10
USAGE

#############################################################################
#Parameters
#############################################################################
my $out_file = 'VennDiagram.tiff';
my $width = 10;
my $height = 10;
my $format = '';
die $usage 
	unless GetOptions (
			"o|outfile:s" => \$out_file,
			"l|label=s" => \$label,
			"f|format=s" => \$format,
			"w|width=f"    => \$width,
			"h|heigh=f"    => \$height,
			);
			
my $matrix_file = shift;
if ($matrix_file) {
	die "input file does not exists " unless -e $matrix_file;
}else{
	die "Not input file"
}			
can_run('Rscript') or die 'R is not installed!';
$matrix_file = abs_path($matrix_file);
  
#############################################################################
#Getting the experimental design
#This is sutiable for differentical test
#Material_Treatment_ReplicatNumber
############################################################################# 
open my $group_fh, "<", $matrix_file
	or die "Canot open $matrix_file\n$!";
my $header = readline $group_fh;
close $group_fh;
my @samples = split /\t/, $header;
shift @samples;

open my $design_fh, ">", ""

my %groups;
foreach my $sample ( sort @samples ) {
	
	my $basename = $sample;
	
	if ( $basename =~ /^(.+)_R(\d+)/ ) {
		my $group = $1;
		my $replicate = $2;
		$groups{$group}{$replicate} = $sample;
	}
	else {
		warn "Not supported group format:$basename\nTeat as samples without replicates\n";
		my $group = $basename;
		my $replicate = 1;
		$groups{$group}{$replicate} = $sample;
	}
	
}

#############################################################################
#Construct design files
############################################################################# 
my $design_fh = File::Temp->new(TEMPLATE => "temp.design.XXXXX",
                        DIR => './',
                        SUFFIX => '.txt');
my $D_fname = $design_fh->filename;
print $design_fh "Sample\tCondition\n";
foreach my $group (sort keys %groups) {
	foreach my $sample (sort keys %{$groups{$group}) {
		print $design_fh "$sample\t$group";
	}
}
close $design_fh;


#############################################################################
#Construct format
#############################################################################

if ($format) {
	if ($format =~ /(pdata)|(bmp)|(jpeg)|(png)/) {
		$format = $1;
		$out_file .= ".$format";
	}else{
		die "unsupported format:$format\n";
	}
}else{
	my @suffixlist = qw(.pdata .bmp .jpeg .png);
	my ($name,$path,$suffix) = fileparse($out_file,@suffixlist);
	if ($suffix) {
		$suffix =~ s/\.//;
		$format = $suffix;
	}else{
		warn 'cannot guess from output file, using pdata...\n';
		$format = 'pdata';
		$out_file .= ".$format";
	}
}




##############################
#Construct R code
##############################
my $R_fh = File::Temp->new(TEMPLATE => "temp.pca.XXXXX",
                        DIR => './',
                        SUFFIX => '.R');
my $R_fname = $R_fh->filename;
print $R_fh pca();
close $R_fh;



##############################
#Run R code
##############################

my $command = "Rscript $fname $file $mode $basename";

my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
  run( command => $command, verbose => 0 );
if ($success) {
}
else {
	my $stderr = join "\n", @$stderr_buf;
	warn "Something went wrong:\n$stderr";
}		

chdir $ori_folder;



sub pca {
	return  <<'SCORE';
#!/usr/bin/Rscript


args <- commandArgs(TRUE)
infile <- args[1]
design.file <- args[2]
outfile <- args[3]
basename <- args[3]
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(scales)
library(grid)
library(FactoMineR)

data<-read.table(file=infile,header=T,row.names=1)
head(data)

design <- read.table(file=design.file, header =T, row.names=1)


#convert userinput data and condition list for PCA analysis
dataForPCAinitialize<-function(data,conditionlist){
  data<-t(data)
  data<-data.frame(condition=conditionlist,data)
  return(data)
}

#get ggplot2 output result
getPCAplot <- function(data,conditionlist,isText=FALSE){
    a<-dataForPCAinitialize(data,conditionlist)
    pca <-PCA(a[,2:ncol(a)], scale.unit=T, graph=F)
    ctri<-pca$eig[,2][1:3]
       names(ctri)<-c("PC1","PC2","PC3")
    xlabtemp=paste("PC1 (",round(ctri["PC1"],2),"%)",sep="")
    ylabtemp=paste("PC2 (",round(ctri["PC2"],2),"%)",sep="")
    colnames(pca$ind$coord)[1:3]<-c("PC1","PC2","PC3")
    PC1 <- pca$ind$coord[,"PC1"]
    PC2 <- pca$ind$coord[,"PC2"]
    maxX=max(PC1)*1.5
    maxY=max(PC2)*1.5
    plotdata <- data.frame(Condition=a[,1],PC1,PC2) 
    plotdata$Condition <- factor(plotdata$Condition)
    plot <- ggplot(plotdata, aes(PC1,PC2),environment = environment()) + 
      geom_point(aes(colour = Condition),size = 5) + 
      theme(panel.border = element_rect(linetype = "dashed")) + 
      theme_bw() +
      ylim(-maxY,maxY)+
      xlim(-maxX,maxX)+
      scale_y_continuous(ylabtemp)+ scale_x_continuous(xlabtemp)+
      theme(legend.text = element_text(colour="blue", size = 16, face = "bold")) + 
      theme(legend.justification=c(1,0),legend.position="top")+
      theme(legend.title = element_text(colour="black", size=16, face="bold"))+
      scale_fill_brewer(palette="Spectral")
    if(isText){
      plot<-plot+geom_text(aes(colour = Condition, label=rownames(plotdata)), size=5, hjust=0.5, vjust=-0.5)
    }
    print(plot)
}


$format(file = outfile, width = $width, height = $height)
getPCAplot(data,design,isText=FALSE)
dev.off()

$format(file = paste(outfile, 'text', "$format", sep= "."), width = $width, height = $height)
getPCAplot(data,design,isText=TRUE)
dev.off()

SCORE
}
