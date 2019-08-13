#!/usr/bin/perl -w

use strict;
use 5.010;
use Data::Printer;
use File::Basename;
use Cwd 'abs_path';
use Getopt::Long;

my $usage = <<USAGE;
merge_ASE.pl -rmats label1,rMATS_folder1 -mats label2,rMATS_folder1 -asta label3,file1
	-asta    ASTALAVISTA results files
	-rmats   rMAT output files
	-mode    total|novel|known, default total
	-o       output prefix
USAGE

#########################################################################
##setting the parameters
#########################################################################
my @asta_files = ();
my @rmats_folders = ();
my $mode = "total";
my $out_prefix = 'merged';
die $usage
  unless GetOptions(
    "asta:s"    => \@asta_files,
    "rmats:s"       => \@rmats_folders,
    "mode:s"       => \$mode,
    "o:s"  => \$out_prefix,
  );

my %files = ();
#pre-check rMATS output files
foreach my $rmats_folder (@rmats_folders) {
	my ($label, $folder) = split /,/, $rmats_folder;
		if (-d $folder && -e $folder) {
		$files{rmats}{$label} = abs_path($folder);		
	}else{
		die "rMATs folder: $folder does not exits";	
	}
}

#pre-check ASTALAVISTA output files
foreach my $asta_file (@asta_files) {
	my ($label, $file) = split /,/, $asta_file;
	if (-e $asta_file) {
		$files{asta}{$label} = $file;		
	}else{
		die "asta file: $file does not exits";	
		
	}
}

#########################################################################
##Read and store AS events
#########################################################################
my %ases = ();
my $count = 0;
my %ases_idx = ();

open my $trackingid_fh, ">", "$out_prefix.tracking_id.txt";
print $trackingid_fh "Label\tFile\tOriginalID\tACCESSION\n";


#read the rMATS outputs
foreach my $label (sort keys %{$files{rmats}}) {
	my $folder = $files{rmats}{$label};
	warn "Proceeding $label...\n";
	opendir my $dir_fh, "$folder" 
		or die "Cannot open $folder:$!";
	while (my $filename = readdir $dir_fh) {
	
		#skip '.' and '..'
		next if $filename =~ /^\./;
		
		#skip the novel or known Events		
		if ($mode eq "novel") {
			if ($filename !~ "novel") {
				warn "Skipping events file:$filename\n";
				next;
			}
		}elsif ($mode eq "known") {
			if ($filename =~ "novel") {
				warn "Skipping events file:$filename\n";
				next;
			}
		}
		
		#read ASE file
		my %rmats_ase = read_rMATS_file("$folder/$filename");
		
		#check whether this events was stored
		foreach my $id (keys %rmats_ase) {
			my $seqid = $rmats_ase{$id}{seqid};
			my $as_type = $rmats_ase{$id}{as_type};
			my $IncForm = $rmats_ase{$id}{IncForm};
			my $SkipForm = $rmats_ase{$id}{SkipForm};
			
			my $accession = '';
			if (exists $ases_idx{$as_type}{$seqid}{$IncForm}{$SkipForm}) {
				$accession =  $ases_idx{$as_type}{$seqid}{$IncForm}{$SkipForm};
			}else{
				$count++;
				$accession = "ASE" . sprintf "%06d", $count;
				$ases{$accession} = $rmats_ase{$id};
				$ases_idx{$as_type}{$seqid}{$IncForm}{$SkipForm} = $accession;
			}
			print $trackingid_fh "$label\t$folder/$filename\t$id\t$accession\n";
		}
	}
}


#read the astalavista outputs


#output ase
open my $ase_out_fh, ">", "$out_prefix.ase.txt";
print $ase_out_fh "ACCESSION\tGID\tSymbol\tType\tSeqID\tStrand\t\tAnnotationi\tIncForm\tSkipForm\n";
foreach my $acc (sort keys %ases) {
	my $gid = $ases{$acc}{gid};
	my $geneSymbol = $ases{$acc}{geneSymbol};
	my $as_type = $ases{$acc}{as_type};
	my $seqid = $ases{$acc}{seqid};
	my $strand = $ases{$acc}{strand};
	my $anno_type = $ases{$acc}{anno_type};
	my $IncForm = $ases{$acc}{IncForm};
	my $SkipForm = $ases{$acc}{SkipForm};
	print $ase_out_fh "$acc\t$gid\t$geneSymbol\t$as_type\t$seqid\t$strand\t$anno_type\t$IncForm\t$SkipForm\n";
}

sub read_rMATS_file {
	my $file = shift;
	my $basename = basename($file, ".txt");
	my $anno_type = $basename =~ /novel/ ? "novel" : "known";
	
	my @tmps = split /\./, $basename;
	my $as_type = $tmps[-1];
	my %ase = ();
	
	open my $rmats_fh, "<", $file;
	readline $rmats_fh;
	while(<$rmats_fh>) {
		chomp;
		my @lines = split /\t/;
		my ($id, $gid, $geneSymbol, $seqid, $strand) = splice(@lines, 0, 5);
		$gid =~ s/\"//g;
		$geneSymbol =~ s/\"//g;
		my ($IncForm, $SkipForm );
		given ($as_type) {
			when (/RI/) {
				my ($riExonStart, $riExonEnd,
					$upstreamES, $upstreamEE, 
					$downstreamES, $downstreamEE) = @lines;
				$riExonStart += 1;
				$upstreamES += 1;
				$downstreamES += 1;
				
				$IncForm = "$riExonStart-$riExonEnd";
				$SkipForm = sort_exon_in_isoform("$upstreamES-$upstreamEE|$downstreamES-$downstreamEE");				
			}
			
			when (/A3SS|A5SS/) {
				my ($longExonStart, $longExonEnd,
					$shortES, $shortEE, 
					$flankingES, $flankingEE) = @lines;
				$longExonStart += 1;
				$shortES += 1;
				$flankingES += 1;
				
				$IncForm = sort_exon_in_isoform("$longExonStart-$longExonEnd|$flankingES-$flankingEE");
				$SkipForm = sort_exon_in_isoform("$shortES-$shortEE|$flankingES-$flankingEE");						
			}
			
			when (/SE/) {
				my ($exonStart, $exonEnd,
					$upstreamES, $upstreamEE, 
					$downstreamES, $downstreamEE) = @lines;	
				$exonStart += 1;
				$upstreamES += 1;
				$downstreamES += 1;	
				
				$IncForm = sort_exon_in_isoform("$upstreamES-$upstreamEE|$exonStart-$exonEnd|$downstreamES-$downstreamEE");
				$SkipForm = sort_exon_in_isoform("$upstreamES-$upstreamEE|$downstreamES-$downstreamEE");					
			}

			when (/MXE/) {
				my ($Exon1stStart, $Exon1stEnd,
					$Exon2ndStart, $Exon2ndEnd, 
					$upstreamES, $upstreamEE, 
					$downstreamES, $downstreamEE) = @lines;	
					
				$Exon1stStart += 1;
				$Exon2ndStart += 1;
				$upstreamES += 1;
				$downstreamES += 1;	
				
				$IncForm = sort_exon_in_isoform("$upstreamES-$upstreamEE|$Exon1stStart-$Exon1stEnd|$downstreamES-$downstreamEE");
				$SkipForm = sort_exon_in_isoform("$upstreamES-$upstreamEE|$Exon2ndStart-$Exon2ndEnd|$downstreamES-$downstreamEE");					
			}

			when (/AFE/) { #skip this type of AS
				next;
				
			}

			when (/ALE/) { #skip this type of AS
				next;
				
			}			
			
		}
		$ase{$id}{as_type} = $as_type;
		$ase{$id}{anno_type} = $anno_type;
		$ase{$id}{gid} = $gid;
		$ase{$id}{geneSymbol} = $geneSymbol;
		$ase{$id}{seqid} = $seqid;
		$ase{$id}{strand} = $strand;
		$ase{$id}{IncForm} = $IncForm;
		$ase{$id}{SkipForm} = $SkipForm;
	}
	return %ase;
}


sub sort_exon_in_isoform {
	my $isoform = shift;
	my @array = ();
	my $seqid;
	foreach (split /\|/, $isoform) {
		my ($start, $end) = split /-|,|:/;
		($start, $end) = sort {$a <=> $b} ($start, $end);
		push @array, [$start, $end];
	}
	@array = sort {$a->[0] <=> $b->[0]} @array;
	my @strs = ();
	foreach (@array) {
		push @strs, $_->[0].'-'.$_->[1];
	}
	return join "|", @strs;
}



