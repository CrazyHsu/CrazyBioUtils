#!/usr/bin/perl -w

use strict;
use 5.010;
use Data::Printer;
use File::Basename;
use Cwd 'abs_path';
use Getopt::Long;

my $usage = <<USAGE;
merge_rMATS_PSI.pl -rmats label1,label2,rMATS_folder1 -mats label2,label3,rMATS_folder2
	-ref     reference ase file 
	-c       read count cutoff, total reads suported the AS events
	-mode    total, known, novel
	-j       junction count only, default off
	-rmats   rMAT output files
	-o       output prefix
USAGE

#########################################################################
##setting the parameters
#########################################################################
my $ase_ref_file = ();
my @rmats_folders = ();
my $mode = 'total';
my $out_prefix = 'merged';
my $junction_count = 0;
my $count_cutoff = 3;
die $usage
  unless GetOptions(
    "ref:s"    => \$ase_ref_file,
    "rmats:s"       => \@rmats_folders,
    "mode:s"       => \$mode,
    "o:s"  => \$out_prefix,
    "j"  => \$junction_count,
    "c:i"  => \$count_cutoff,
  );
  
my %files = ();
#pre-check rMATS output files
foreach my $rmats_folder (@rmats_folders) {
	my ($label1, $label2, $folder) = split /,/, $rmats_folder;
		if (-d $folder && -e $folder && -e "$folder/ASEvents" && "$folder/MATS_output") {
		$files{rmats}{"${label1}_vs_${label2}"}{events} = abs_path("$folder/ASEvents");
		$files{rmats}{"${label1}_vs_${label2}"}{results} = abs_path("$folder/MATS_output");
	}else{
		die "rMATs folder: $folder seems not rMATs folder";	
	}
}


#########################################################################
##Reconstruct AS events
#########################################################################
if ($ase_ref_file && -e $ase_ref_file) {
	warn "Reading specified reference ASE file\n";
}else {
	warn "No refererence ASE file inputed, reconstructing one\n";

	my $count = 0;
	my %ases_idx = ();
	my %ases = ();


	open my $trackingid_fh, ">", "$out_prefix.tracking_id.txt";
	print $trackingid_fh "Label\tFile\tOriginalID\tACCESSION\n";


	#read the rMATS outputs
	foreach my $label (sort keys %{$files{rmats}}) {
		my $folder = $files{rmats}{$label}{events};
		warn "Reading Events Annotation $label...\n";
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
			my %rmats_ase = read_rMATS_events_file("$folder/$filename");
			
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
					$accession = "ASE" . sprintf "%08d", $count;
					$ases{$accession} = $rmats_ase{$id};
					$ases_idx{$as_type}{$seqid}{$IncForm}{$SkipForm} = $accession;
				}
				print $trackingid_fh "$label\t$folder/$filename\t$id\t$accession\n";
			}
		}
	}

	#output ase
	open my $ase_out_fh, ">", "$out_prefix.ase.txt";
	print $ase_out_fh "ACCESSION\tGID\tSymbol\tType\tSeqID\tStrand\tAnnotation\tIncForm\tSkipForm\n";
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
	$ase_ref_file = "$out_prefix.ase.txt";
	
}

#read reference ASE file
my %ase_idx =();
my %ase = ();

open my $ase_fh, "<", $ase_ref_file 
	or die "Cannot open $ase_ref_file:$!\n";
readline $ase_fh;
while (<$ase_fh>) {
	chomp;
	my ($acc, $gid, $geneSymbol, $as_type, $seqid, $strand, 
		$anno_type, $IncForm, $SkipForm) = split /\t/;
	
	#store to ase 	
	$ase{$acc}{as_type} = $as_type;
	$ase{$acc}{anno_type} = $anno_type;
	$ase{$acc}{gid} = $gid;
	$ase{$acc}{geneSymbol} = $geneSymbol;
	$ase{$acc}{seqid} = $seqid;
	$ase{$acc}{strand} = $strand;
	$ase{$acc}{IncForm} = $IncForm;
	$ase{$acc}{SkipForm} = $SkipForm;

	#store to ase index
	$ase_idx{$seqid}{$strand}{$IncForm}{$SkipForm} = $acc;
	
}
#########################################################
#Read the rMATS Events Percent Spliced In
#########################################################
my %samples_tmp = ();
my %labels_tmp = ();
foreach my $label (sort keys %{$files{rmats}}) {
	my $folder = $files{rmats}{$label}{results};

	my ($sample1, $sample2) = split /_vs_/, $label;
	$labels_tmp{$label} = 1;
	$samples_tmp{$sample1} = 1;
	$samples_tmp{$sample2} = 1;
	
	warn "Reading Events Annotation $label...\n";
	opendir my $dir_fh, "$folder" 
		or die "Cannot open $folder:$!";
	while (my $filename = readdir $dir_fh) {
	
		#skip '.' and '..'
		next if $filename =~ /^\./;
		if ($junction_count) {
			next if $filename !~ /JunctionCountOnly/;
		}else{
			next if $filename !~ /ReadsOnTargetAndJunctionCounts/;		
		}
		#read ASE file
		my %rmats_ase = read_rMATS_output_file($label,"$folder/$filename");
		
		foreach my $id (keys %rmats_ase) {
			my $seqid = $rmats_ase{$id}{seqid};
			my $strand = $rmats_ase{$id}{strand};
			my $IncForm = $rmats_ase{$id}{IncForm};
			my $SkipForm = $rmats_ase{$id}{SkipForm};
			
			my $acc = '';
			if (exists $ase_idx{$seqid}{$strand}{$IncForm}{$SkipForm}) {
				$acc = $ase_idx{$seqid}{$strand}{$IncForm}{$SkipForm};
				$ase{$acc}{samples}{$sample1} = $rmats_ase{$id}{samples}{$sample1};
				$ase{$acc}{samples}{$sample2} = $rmats_ase{$id}{samples}{$sample2};
				$ase{$acc}{versus}{$label} = $rmats_ase{$id}{versus}{$label};
			}else{
				warn "AS Events was not annotated\n";
			}
		}
	}
}

my @samples = sort keys %samples_tmp;
my @labels = sort keys %labels_tmp;
#########################################################
#Read the rMATS Events Percent Spliced In
#########################################################
open my $ase_out_fh, ">", "$out_prefix.PSI.txt";
print $ase_out_fh "ACCESSION\tGID\tSymbol\tType\tSeqID\tStrand\tAnnotation\tIncForm\tSkipForm\t";
print $ase_out_fh join "\t", @samples;
print $ase_out_fh "\n";

open my $ase_pvalue_fh, ">", "$out_prefix.pvalues.txt";
print $ase_pvalue_fh "ACCESSION\tGID\tSymbol\tType\tSeqID\tStrand\tAnnotation\tIncForm\tSkipForm\t";
print $ase_pvalue_fh join "\t", @labels;
print $ase_pvalue_fh "\n";

open my $ase_fdr_fh, ">", "$out_prefix.fdr.txt";
print $ase_fdr_fh "ACCESSION\tGID\tSymbol\tType\tSeqID\tStrand\tAnnotation\tIncForm\tSkipForm\t";
print $ase_fdr_fh join "\t", @labels;
print $ase_fdr_fh "\n";
foreach my $acc (sort keys %ase) {
	my $gid = $ase{$acc}{gid};
	my $geneSymbol = $ase{$acc}{geneSymbol};
	my $as_type = $ase{$acc}{as_type};
	my $seqid = $ase{$acc}{seqid};
	my $strand = $ase{$acc}{strand};
	my $anno_type = $ase{$acc}{anno_type};
	my $IncForm = $ase{$acc}{IncForm};
	my $SkipForm = $ase{$acc}{SkipForm};
	
	my @sample_PSI = ();
	foreach my $sample (@samples) {
		my $psi = exists $ase{$acc}{samples}{$sample} ? $ase{$acc}{samples}{$sample}{IncLevel} : "NA";
		push @sample_PSI, $psi;
	}
	
	print $ase_out_fh "$acc\t$gid\t$geneSymbol\t$as_type\t$seqid\t$strand\t$anno_type\t$IncForm\t$SkipForm\t";
	print $ase_out_fh join "\t", @sample_PSI;
	print $ase_out_fh "\n";

	my @label_pvalues = ();
	my @label_fdrs = ();
	foreach my $label (@labels) {
		my ($pvalue, $fdr) = ("NA")x2;
		if (exists $ase{$acc}{versus}{$label}) {
			$pvalue =  $ase{$acc}{versus}{$label}{pvalue} ;
			$fdr =  $ase{$acc}{versus}{$label}{fdr} ;
		}
		push @label_pvalues, $pvalue;
		push @label_fdrs, $fdr;
	}		
	print $ase_pvalue_fh "$acc\t$gid\t$geneSymbol\t$as_type\t$seqid\t$strand\t$anno_type\t$IncForm\t$SkipForm\t";
	print $ase_pvalue_fh join "\t", @label_pvalues;
	print $ase_pvalue_fh "\n";

	print $ase_fdr_fh "$acc\t$gid\t$geneSymbol\t$as_type\t$seqid\t$strand\t$anno_type\t$IncForm\t$SkipForm\t";
	print $ase_fdr_fh join "\t", @label_fdrs;
	print $ase_fdr_fh "\n";	
	
}



sub read_rMATS_output_file {
	my ($label, $file) = @_;
	
	my ($sample1, $sample2) = split /_vs_/, $label;
	
	my $basename = basename($file, ".txt");
	
	my @tmps = split /\./, $basename;
	my $as_type = $tmps[0];

	my %ase;
	open my $rmats_fh, "<", $file;
	readline $rmats_fh;
	while(<$rmats_fh>) {
		chomp;
		my @lines = split /\t/;
		my ($id, $gid, $geneSymbol, $seqid, $strand) = splice(@lines, 0, 5);
		$gid =~ s/\"//g;
		$geneSymbol =~ s/\"//g;
		
		$ase{$id}{as_type} = $as_type;
		$ase{$id}{gid} = $gid;
		$ase{$id}{geneSymbol} = $geneSymbol;
		$ase{$id}{seqid} = $seqid;
		$ase{$id}{strand} = $strand;		
		
		my ($IncForm, $SkipForm );
		given ($as_type) {
			when (/RI/) {
				my ($riExonStart, $riExonEnd,
					$upstreamES, $upstreamEE, 
					$downstreamES, $downstreamEE) = splice(@lines, 0, 6);
				$riExonStart += 1;
				$upstreamES += 1;
				$downstreamES += 1;
				
				$IncForm = "$riExonStart-$riExonEnd";
				$SkipForm = sort_exon_in_isoform("$upstreamES-$upstreamEE|$downstreamES-$downstreamEE");				
			}
			
			when (/A3SS|A5SS/) {
				my ($longExonStart, $longExonEnd,
					$shortES, $shortEE, 
					$flankingES, $flankingEE) = splice(@lines, 0, 6);
				$longExonStart += 1;
				$shortES += 1;
				$flankingES += 1;
				
				$IncForm = sort_exon_in_isoform("$longExonStart-$longExonEnd|$flankingES-$flankingEE");
				$SkipForm = sort_exon_in_isoform("$shortES-$shortEE|$flankingES-$flankingEE");						
			}
			
			when (/SE/) {
				my ($exonStart, $exonEnd,
					$upstreamES, $upstreamEE, 
					$downstreamES, $downstreamEE) = splice(@lines, 0, 6);
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
					$downstreamES, $downstreamEE) = splice(@lines, 0, 8);	
					
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
		$ase{$id}{IncForm} = $IncForm;
		$ase{$id}{SkipForm} = $SkipForm;
	
	
		my ($id2, $ijc_sample1, $sjc_sample1, $ijc_sample2, $sjc_sample2, 
		$IncFormLen, $SkipFormLen, $pvalue, $fdr, 
		$IncLevel1, $IncLevel2, $IncLevelDifference) = @lines;
		
		$ase{$id}{samples}{$sample1}{IncFormCount} = $ijc_sample1;
		$ase{$id}{samples}{$sample1}{SkipFormCount} = $sjc_sample1;
		$IncLevel1 = ($ijc_sample1 + $sjc_sample1) >= $count_cutoff ?  $IncLevel1 : "NA";
		$ase{$id}{samples}{$sample1}{IncLevel} = $IncLevel1;
		
		$ase{$id}{samples}{$sample2}{IncFormCount} = $ijc_sample2;
		$ase{$id}{samples}{$sample2}{SkipFormCount} = $sjc_sample2;
		$IncLevel2 = ($ijc_sample2 + $sjc_sample2) >= $count_cutoff ?  $IncLevel2 : "NA";
		$ase{$id}{samples}{$sample2}{IncLevel} = $IncLevel2;
		
		$ase{$id}{versus}{$label}{pvalue} = $pvalue;
		$ase{$id}{versus}{$label}{fdr} = $fdr;
		$ase{$id}{versus}{$label}{IncLevelDifference} = $IncLevelDifference;

	}
	return %ase;	
}



sub read_rMATS_events_file {
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



