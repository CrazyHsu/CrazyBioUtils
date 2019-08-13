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
normalty_test.pl [options] infile.txt

 Options:
   -m|mode            'col' or 'row', default 'col'
   -o|--output        output folder, default './'

USAGE
my $mode       = 'col';
my $out_folder = dirname './';

die $usage
  unless GetOptions(
    "m|mode:s"       => \$mode,
    "o|output:s"   => \$out_folder,

  );
  
my $file = shift;

can_run('Rscript') or die 'R is not installed!';
die "$file does not exits\n" unless -s $file;
$file = abs_path($file);
my $basename = basename($file, qw(.txt .csv .dat .xls));

$out_folder =~ s/[\/|\|]+$//;
mkdir $out_folder unless -d $out_folder;
my $ori_folder = getcwd;
chdir $out_folder;
##############################
#Construct R code
##############################
my $fh = File::Temp->new(TEMPLATE => "temp.$basename.XXXXX",
                        DIR => './',
                        SUFFIX => '.R');
my $fname = $fh->filename;
print $fh normalty_test();
close $fh;

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



sub normalty_test {
	return  <<'SCORE';
#!/usr/bin/Rscript

args <- commandArgs(TRUE)
file <- args[1]
mode <- args[2]
basename <- args[3]
write("Reading file...\n", stderr())
data <- read.table(file = file, header = T, row.names = 1, sep="\t")
if (mode == 'col') {
  pvalue <- rep(1, dim(data)[2])
  names(pvalue) <- colnames(data)
  for (i in 1 : (dim(data)[2]) ){
	pvalue[i] <- as.numeric(shapiro.test( as.numeric(data[, i]))[2]) 
  }
}else{
  pvalue <- rep(1, dim(data)[1])
  names(pvalue) <- rownames(data)
  for (i in 1 : (dim(data)[1]) ){
	pvalue[i] <- as.numeric(shapiro.test( as.numeric(data[i, ]))[2]) 
  }  
}

pvalue <- -log10(pvalue)

pvalue <- as.data.frame(pvalue)
write("Plotting the distribution...\n", stderr())
pdf(file = paste(basename, "pdf", sep = "."), width = 4, height = 3)
library(ggplot2)
ggplot(pvalue, aes(pvalue)) +
  geom_histogram(binwidth=2, colour="black", fill="white") +
  xlab("-log10Pvalue") +
  ylab("Counts") +
  theme_bw()
dev.off()
write("Writing Results...\n", stderr())
write.table(pvalue, file = paste(basename, "dat", sep = "."), sep = "\t", quote = F)

SCORE
}
