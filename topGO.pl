#!/usr/bin/perl -w
use strict;
use 5.010;
# use Data::Printer;
use File::Basename;
use IPC::Cmd qw[can_run run];
use Getopt::Long;
use Cwd qw(abs_path chdir getcwd);
use File::Temp;

my $usage = <<USAGE;
SYSNOPSIS
topGO.pl [options] -r topGO.map -b exp.gid -q diffExp.gid

 Options:
   -m|mode            count or score mode
   -r|--map-file      topGO mapping file
   -b|--background    topGO all genes
   -q|--query         topGO interesting genes
   -o|--output        output folder for removed reads
   -a|--header        gene list file contains header   
   -n|--node-size     topGO node size, default 5
   -p|--pvalue        pvalue, default 0.05

Note:
In count mode:
map file and intersting genes (query) file are mandatory, if no all genes
(backround)is speficied, all genes in the mapping file will be used instead.

In score mode:
map file and all genes files are  mandatory, if no query gene list file is 
specified, a gene selection by by socre will automatically calucated. 
If query gene list is specified, the gene secection by query list will added   

U'd better use go_terms_formation.pl to format ur database before u run it.

USAGE
my $map_file       = 'NA';
my $mode       = 'count';
my $backgournd_file = 'NA';
my $query_file = 'NA';
my $out_folder = dirname './';
my $node_size = 5;
my $pvalue   = 0.05;
my $header_flag = 0;
die $usage
  unless GetOptions(
    "m|mode:s"       => \$mode,
    "b|background:s"       => \$backgournd_file,
    "r|map-file:s"       => \$map_file,
    "q|query:s"       => \$query_file,
    "n|node-size:i"       => \$node_size,
    "o|output:s"   => \$out_folder,
    "a|header"   => \$header_flag,
    "p|pvalue=f" => \$pvalue,
  );

can_run('Rscript') or die 'R is not installed!';
die "$map_file does not exits\n" unless -s $map_file;
$map_file = abs_path($map_file);
$backgournd_file = abs_path($backgournd_file);
$query_file = abs_path($query_file);
$header_flag = $header_flag ? "TRUE" : "FALSE";

$out_folder =~ s/[\/|\|]+$//;
mkdir $out_folder unless -d $out_folder;
my $ori_folder = getcwd;
chdir $out_folder;
##############################
#Construct R code
##############################
my $fh = File::Temp->new(TEMPLATE => 'temp.XXXXX',
                        DIR => './',
                        SUFFIX => '.R');
my $fname = $fh->filename;
print $fh header();
if ($mode eq 'count') {
	print $fh count_mode();
}else{
	print $fh score_mode();	
}
print $fh footer();
close $fh;

##############################
#Run R code
##############################

my $command = "Rscript $fname $map_file $query_file $backgournd_file $node_size $pvalue $header_flag > topGO.stdout.log 2> topGO.stderr.log";

my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
  run( command => $command, verbose => 0 );
if ($success) {
	warn "\tDone!\n";
}
else {
	my $stderr = join "\n", @$stderr_buf;
	warn "Something went wrong:\n$stderr";
}		

chdir $ori_folder;




sub header {
	return <<'HEADER';
#!/usr/bin/env Rscript
##############################################################################
##Note:
##The first backgound list file must contain the score
##If no query gene list file is specified, a gene selection by by socre will 
##automatically calucated. If query gene list is specified, the gene secection
##by query list will added
##############################################################################

library(ggplot2)
library(topGO)
library(Rgraphviz)
sessionInfo()

#variable
args <- commandArgs(TRUE)
annoMapFile = args[1]
queryFile = args[2]
backgourndFile = args[3]
nodeSize = args[4]
pvalue.threshold = args[5]
header_flag = as.logical(args[6])

###################################################################
#basic function
#get range
###################################################################
getRange <- function(num) {
  
  if (num > 1000 ) {
    range = c(5, 10, 20, 30, 50, 100, 200, 300, 500, 1000)
  }else if (num > 500) {
    range = c(5, 10, 20, 30, 50, 100, 200, 300, 500, num)
  }else if (num > 300) {
    range = c(5, 10, 20, 30, 50, 100, 200, 300, num)
  }else if (num > 200) {
    range = c(5, 10, 20, 30, 50, 100, 200, num)
  }else if (num > 100) {
    range = c(5, 10, 20, 30, 50, 100, num)
  }else if (num > 50) {
    range = c(5, 10, 20, 30, 50, num)
  }else if (num > 30) {
    range = c(5, 10, 20, 30, num)
  }else if (num > 20) {
    range = c(5, 10, 20, num)
  }else if (num > 10) {
    range = c(5, 10, num)
  }else if (num > 5) {
    range = c(5, num)
  }else if (num > 0) {
    range = c(num)
  }else if (num == 0) {
    warning("No nodes")
    range = (0)
  }
  range
}

#Reading custom annotations
if (file.exists(annoMapFile)) {
  geneID2GO <- readMappings(file= annoMapFile)
  str(head(geneID2GO))
}else{
  warning("cosutm annotatoin list file does not exists")
  stop()
}

HEADER
}

sub footer {
	return <<'FOOTER';

tGO("BP")
tGO("MF")
tGO("CC")

#remove the genes not in querylist
GO2geneID <- inverseList(geneID2GO)
for (i in 1:length(GO2geneID)) {
  tmp <-unlist(unname(GO2geneID[i]))
  names(tmp) <- tmp
  GO2geneID[i] <- list(unname(tmp[tmp %in% queryGenes]))
}
#remove emptoy go terms
GO2geneID <- GO2geneID[lapply(GO2geneID,length)>0]
#change list to datafrmae
GO2geneID <- lapply(GO2geneID, paste, collapse=" ")
go.terms <- names(GO2geneID)
go.geneID <- unlist(unname(GO2geneID))
GO2geneID<- data.frame(terms = go.terms,
                       ids = go.geneID
                       )
write.table(GO2geneID, file="GO2geneID.txt", sep="\t", quote=F, row.names=F)



FOOTER
}


sub count_mode {
	return <<'COUNT';

#Reading the querylist
if (file.exists(queryFile)) {
  queryGenes <- read.table(file=queryFile, stringsAsFactor=F, header = header_flag)[, 1]
  queryGenes <- unlist(as.vector(queryGenes))
  head(queryGenes)
}else{
  warning("query list file does not exists")
  stop()
}


#Reading the background list
if (file.exists(backgourndFile)) {
  backgourndGenes <- read.table(file=backgourndFile, stringsAsFactor=F, header = header_flag)[, 1]
  backgourndGenes <- unlist(as.vector(backgourndGenes))
  head(backgourndGenes)
}else{
  warning("Back ground annotatoin list doesnot exists, using the total genes")
  backgourndGenes <- names(geneID2GO)
}


#Construct the gene list
geneList <- factor(as.integer(backgourndGenes %in% queryGenes ))
names(geneList) <- backgourndGenes
head(geneList)
str(geneList)


#declear the function
tGO  <- function(ontology) {
  
  #construct GOdata
  GOdata <- new("topGOdata", 
                description = "GO analysis",
                ontology=ontology, 
                allGenes=geneList, 
                nodeSize = nodeSize,
                annot=annFUN.gene2GO, 
                gene2GO = geneID2GO)
  
  
  #Test
  resultFisher.classic <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  resultFisher.elim <- runTest(GOdata, algorithm = "elim", statistic = "fisher")

  topNodesNum <- max(length(score(resultFisher.classic)),
                     length(score(resultFisher.elim))
                     )  
  
  #Summary the results
  allRes <- GenTable(GOdata, 
                     classicFisher = resultFisher.classic, 
                     elimFisher = resultFisher.elim,
                     orderBy = "classicFisher",
                     ranksOf = "classicFisher",
                     topNodes = topNodesNum
					 )
  
  write.table(allRes, file=paste(ontology, ".details.txt", sep = ""), sep="\t", quote=F, row.names=F)
 
  #Got the significant number and Plot
  resultFisher.classic.sig <- length(allRes$classicFisher[allRes$classicFisher < pvalue.threshold])
  resultFisher.elim.sig <- length(allRes$elimFisher[allRes$elimFisher < pvalue.threshold])
  
  range.fisher.classic <- getRange(resultFisher.classic.sig)
  for (i in range.fisher.classic) {
    printGraph(GOdata, resultFisher.classic, firstSigNodes = i, fn.prefix= paste(ontology, "Fisher", sep = ".") , useInfo = "all", pdfSW = TRUE)
  }  
  
  range.fisher.elim <- getRange(resultFisher.elim.sig)
  for (i in range.fisher.elim) {
    printGraph(GOdata, resultFisher.elim, firstSigNodes = i, fn.prefix= paste(ontology, "Fisher", sep = ".") , useInfo = "all", pdfSW = TRUE)    
  }    
   
}

COUNT
}



sub score_mode {
	return  <<'SCORE';
#!/usr/bin/env Rscript

library(ggplot2)
library(topGO)
library(Rgraphviz)
sessionInfo()

#variable
args <- commandArgs(TRUE)
annoMapFile = args[1]
queryFile = args[2]
backgourndFile = args[3]
nodeSize = args[4]
pvalue.threshold = args[5]
header_flag = as.logical(args[6])





#Reading the querylist
if (file.exists(queryFile)) {
  queryGenes.raw <- read.table(file=queryFile, stringsAsFactor=F, header = header_flag)
  if (ncol(queryGenes.raw) == 2) {
    queryGenes <- queryGenes.raw[, 2]
    names(queryGenes) <- queryGenes.raw[, 1]  
  }else{
    queryGenes <- rep(0, nrow(queryGenes.raw))
    names(queryGenes) <- queryGenes.raw[, 1]
  }
  head(queryGenes) 
}else{
  warning("query list file does not exists\n Calculate from background list")
}


#Reading the background list
if (file.exists(backgourndFile)) {  
  backgourndGenes.raw <- read.table(file=backgourndFile, stringsAsFactor=F, header = header_flag)
  if (ncol(backgourndGenes.raw) >= 2) {
    backgourndGenes <- backgourndGenes.raw[, 2]
    names(backgourndGenes) <- backgourndGenes.raw[, 1]  
    head(backgourndGenes)   
  }else{
    warning("back groud list have no score")
    stop()
  }

}else{
  warning("Back ground annotatoin list doesnot exists, using the total genes")
  stop()
}


#Construct the gene list
geneList <- backgourndGenes

if (exists("queryGenes")) {
  diffGenes <- function(allScore) {
    return(allScore %in% queryGenes)
  }  
}else{
  diffGenes <- function(allScore) {
    return(allScore < pvalue.threshold)
  }  
}

#declear the function
tGO  <- function(ontology) {
  
  #construct GOdata
  GOdata <- new("topGOdata", 
                description = "GO analysis",
                ontology=ontology, 
                allGenes=geneList, 
                geneSel = diffGenes,
                nodeSize = nodeSize,
                annot=annFUN.gene2GO, 
                gene2GO = geneID2GO)
  
  
  #Test
  resultFisher.classic <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  resultFisher.elim <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
  resultKS.classic <- runTest(GOdata, algorithm = "classic", statistic = "ks")
  resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")

  topNodesNum <- max(length(score(resultFisher.classic)),
					length(score(resultFisher.elim)),
					length(score(resultKS.classic)),
					length(score(resultKS.elim))
                     )   
  #Summary the results
  allRes <- GenTable(GOdata, 
                     classicFisher = resultFisher.classic, 
                     classicKS = resultKS.classic,
                     elimFisher = resultFisher.elim,
                     elimKS = resultKS.elim,
                     orderBy = "classicFisher",
                     ranksOf = "classicFisher",
                     topNodes = topNodesNum
					 )
  write.table(allRes, file=paste(ontology, ".details.txt", sep = ""), sep="\t", quote=F, row.names=F)
					 
  
  #Got the significant number and Plot
  resultFisher.classic.sig <- length(allRes$classicFisher[allRes$classicFisher < pvalue.threshold])
  resultFisher.elim.sig <- length(allRes$elimFisher[allRes$elimFisher < pvalue.threshold])
  resultKS.classic.sig <- length(allRes$classicKS[allRes$classicKS < pvalue.threshold])
  resultKS.elim.sig <- length(allRes$elimKS[allRes$elimKS < pvalue.threshold])
  
  range.fisher.classic <- getRange(resultFisher.classic.sig)
  for (i in range.fisher.classic) {
    printGraph(GOdata, resultFisher.classic, firstSigNodes = i, fn.prefix= paste(ontology, "Fisher", sep = ".") , useInfo = "all", pdfSW = TRUE, .NO.CHAR = 30, useFullNames = T)
  }  
  
  range.fisher.elim <- getRange(resultFisher.elim.sig)
  for (i in range.fisher.elim) {
    printGraph(GOdata, resultFisher.elim, firstSigNodes = i, fn.prefix= paste(ontology, "Fisher", sep = ".") , useInfo = "all", pdfSW = TRUE, .NO.CHAR = 30, useFullNames = T)   
  }   
  
  range.ks.classic <- getRange(resultKS.classic.sig)
  for (i in range.ks.classic) {
    printGraph(GOdata, resultKS.classic, firstSigNodes = i, fn.prefix= paste(ontology, "KS", sep = ".") , useInfo = "all", pdfSW = TRUE, .NO.CHAR = 30, useFullNames = T)   
  }   
  
  range.ks.elim <- getRange(resultFisher.elim.sig)
  for (i in range.ks.elim) {
    printGraph(GOdata, resultKS.elim, firstSigNodes = i, fn.prefix= paste(ontology, "KS", sep = ".") , useInfo = "all", pdfSW = TRUE, .NO.CHAR = 30, useFullNames = T)  
  }   
}

SCORE
}
