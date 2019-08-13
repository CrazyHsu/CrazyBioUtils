#!/usr/bin/perl -w

use strict;
use 5.010;
use Data::Printer;
use File::Basename;

my @genomes = @ARGV;

my @fqs = glob "/home/hpyu/Wheat/CleanReads/*.fq.gz";
my %groups;
foreach my $fq (@fqs) {
	my $basename = basename($fq, ".fq.gz");
	my @temp = split '_', $basename;
	pop @temp;
	my $material = join "_", @temp;
	$groups{$material} = 1;
}

#write merege files


foreach my $genome (@genomes) {
	open my $fh, ">", "assembly_list.txt" or die "cannot create files:\n";
	foreach my $material (sort keys %groups) {
		print $fh "$genome/$material/transcripts.gtf\n"
	}
	close $fh;
	mkdir "$genome/cuffmerge" unless "$genome/cuffmerge";
	warn "cuffmerge -o $genome/cuffmerge -g ../../Genome/IWGSC/ta_IWGSC_MIPSv2.2_HighConf_${genome}_2014Jul18.gtf -p 8 -s ../../Genome/IWGSC/ta_IWGSC_CSSassembly_MIPSv2REF_${genome}_cleaned_2014Jul16.fa assembly_list.txt 2> ${genome}/cuffmerge.stderr.log > ${genome}/cuffmerge.stdout.log\n";
	system "cuffmerge -o $genome/cuffmerge -g ../../Genome/IWGSC/ta_IWGSC_MIPSv2.2_HighConf_${genome}_2014Jul18.gtf -p 8 -s ../../Genome/IWGSC/ta_IWGSC_CSSassembly_MIPSv2REF_${genome}_cleaned_2014Jul16.fa assembly_list.txt 2> ${genome}/cuffmerge.stderr.log > ${genome}/cuffmerge.stdout.log";

	system "ln $genome/cuffmerge/merged.gtf $genome.gtf";
}

