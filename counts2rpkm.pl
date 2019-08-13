#!/usr/bin/perl -w

use strict;
use 5.010;
use Cwd 'abs_path';
use Statistics::R;
use Data::Printer;

my $usage = <<USAGE;

counts2rpkm V1.0.0, written by corephi
This program is used to transform read count to RPKM
---------------------------------------------------------------------------------------
Usage: counts2rpkm.pl count.txt

Note:
The first column must be GID, and the last one must be length

Samples:
GID	sample1_R1	sample1_R2	sample2_R1	sample2_R2 Length

USAGE


#check input file 
my $infile       = shift @ARGV ;
die "infile:$infile does not exists\n" unless -e $infile;
my $in_fh;
open $in_fh, "<", $infile or die "cannot open file $infile:$!\n";


###############################################################################
#read group info
###############################################################################
my $header = readline $in_fh;
chomp($header);
my @samples = split /\t/, $header; 
shift @samples;
pop @samples;
close $in_fh;

my %groups_info = ();
foreach my $sample (@samples) {
	if ( $sample =~ /^(.+)_R(\d+)/ ) {
		my $group = $1;
		my $replicate = $2;
		$groups_info{$sample}{group} = $group;
		$groups_info{$sample}{replicate} = $replicate;
	}
	else {
		warn "Not supported group format:$sample\nTeat as samples without replicates\n";
		my $group = $sample;
		$groups_info{$sample}{group} = $group;
		$groups_info{$sample}{replicate} = 1;
	}	
}


my @groups = map {$groups_info{$_}{group}} @samples;
my $group_string_in_r = join ", ", map {"\"$_\""} @groups;
$group_string_in_r = "c($group_string_in_r)";
my $sample_num = @samples;
my $len_col = $sample_num + 1;


###############################################################################
#Calculate RPKM in sample
###############################################################################


my $Rocode = <<edgeR;
	
library("edgeR")
sessionInfo()
# read file
rawdata <- read.delim("$infile",row.names=1, stringsAsFactors=FALSE, header = TRUE)
rawdata <-  rawdata[grep("^_", row.names(rawdata), invert = TRUE), ]
head(rawdata)

data <- DGEList(
  counts=rawdata[,1:$sample_num], 
  group = $group_string_in_r, 
  genes = data.frame(Length=rawdata[,$len_col], 
                     row.names = rownames(rawdata))
  )


# Normalizing
data <- calcNormFactors(data)
write.table(data\$samples, file = "normFactor.txt", quote = F, sep = "\t") 


# output the normalized counts
normFactors <- data\$samples\$norm.factors
normCounts <- t( t(data\$counts) / normFactors )
normCounts <- round(normCounts)
colnames(normCounts) <- colnames(data\$counts)

write.table(normCounts, file = "Counts.normed.txt", quote = F, sep = "\t") 

# calculat rpkm
rpkm <- rpkm(data, data\$genes\$Length, normalized.lib.sizes=F )
head(rpkm)

write.table(rpkm, file = "RPKM.sample.txt", quote = F, sep = "\t") 
	
edgeR

open my $log_fh, ">", "rpkm.log" or die "cannot create data file\n";
my $R     = Statistics::R->new();
my $R_out = $R->run($Rocode);
print $log_fh $Rocode, "\n";
print $log_fh $R_out, "\n";
close $log_fh;
$R->stop;


###############################################################################
#Calculate RPKM in sample
###############################################################################
open my $rpkm_sample_fh, "<", "RPKM.sample.txt";
open my $rpkm_group_fh, ">", "RPKM.group.txt";
my $rpkm_first_line = readline $rpkm_sample_fh;
chomp ($rpkm_first_line);
my @samples_rpkm = map {s/^X(\d+)/$1/; $_} split /\t/, $rpkm_first_line; 

my %rpkms = ();

#read the rpkm
while (<$rpkm_sample_fh>) {
	chomp;
	my @lines = split /\t/;
	my $gid = shift @lines;
	foreach my $i (0 .. $#lines) {
		my $sample = $samples_rpkm[$i];
		my $rpkm = $lines[$i];
		my $group = $groups_info{$sample}{group};
		my $replicate = $groups_info{$sample}{replicate};
		$rpkms{$gid}{$group}{$replicate} = $rpkm;
	}	
}
close $rpkm_sample_fh;

#calculte and out put the rpkm by group
my %groups_tmp = map {$_ => 1} @groups;
my @sample_groups = sort keys %groups_tmp ; 
print $rpkm_group_fh "GID\t", join "\t", @sample_groups;
print $rpkm_group_fh "\n";
foreach my $gid (sort keys %rpkms) {
	my @lines = ($gid);
	foreach my $group (@sample_groups) {
		my $rpkm = 0;
		foreach my $replicate (keys %{$rpkms{$gid}{$group}} ) {
			my $tmp = $rpkms{$gid}{$group}{$replicate};
			next if $rpkm eq "Inf" ;
			if ($tmp eq 'NA' ){
				$rpkm += 0;
			}elsif ($tmp eq 'Inf') {
				$rpkm = "Inf"
			}else{
				$rpkm += $rpkms{$gid}{$group}{$replicate};
			}
		}
		my $replicate_num = keys %{$rpkms{$gid}{$group}};
		push @lines, $rpkm / $replicate_num ;
	}
	push @lines, "\n";
	print $rpkm_group_fh join "\t", @lines;
}
close $rpkm_group_fh;
