#!/usr/bin/perl -w

use strict;
use 5.010;
use YAML;
use Getopt::Long;
use Data::Printer;
use File::Basename;

my $usage = <<USAGE;
SYSNOPSIS
MapmanOntology.pl -r annotation.gmt -b background.txt -q query_list

 Options:
   -r     Mapman annotation file with gmt format
   -b     background gene list
   -q     query gene list
   -o     output file
   -s     minimal gene number wthin one Mapman Germs, default 5.

USAGE

my $mapman_anno_file      = '';
my $backgourd_file        = '';
my $query_file            = '';
my $node_size = 5; 
my $output_file = 'MapmanOntology.xls';
die $usage
  unless GetOptions(
    "r=s"      => \$mapman_anno_file,
    "s:i"      => \$node_size,
    "b:s" => \$backgourd_file,
    "q:s" => \$query_file,
    "o:s" => \$output_file,

  );
die "Mapman annotation file:$mapman_anno_file does not exits\n" 
	unless -e  $mapman_anno_file;
die "Query gene list file:$query_file does not exits\n" 
	unless -e  $query_file;
	

#read annoation
my %anno = read_gmt_annotation($mapman_anno_file);
#read querylist
my %query_list = read_list($query_file, $anno{gids});
#read backgroudlist
my %background;
if ($backgourd_file) {
	my %background_list = read_list($backgourd_file, $anno{gids});
	foreach my $bin (keys %{$anno{bins}}) {
		my $name = $anno{bins}{$bin}{name};
		my @tmp_gids1 = grep {exists $background_list{$_} } keys %{$anno{bins}{$bin}{gids}};
		my %tmp_gids2 = map {$_ => 1 } @tmp_gids1 ;
		if (@tmp_gids1) {
			$background{bins}{$bin}{name} = $name;
			$background{bins}{$bin}{gids} = \%tmp_gids2;
		}
	}
	foreach my $gid (keys %{$anno{gids}}) {
		$background{gids}{$gid} = $anno{gids}{$gid};
	}
}else{
	%background = %anno;
}

#enrichment analysis
my $bgtotal = keys %{$background{gids}};
my $querytotal = keys %query_list;
my @bins = keys %{$background{bins}};
my @pvalues = ();
my @valid_bins  = ();
foreach my $bin (@bins) {
	my %gids_in_bin =  %{$background{bins}{$bin}{gids}};
	my $bgitem = keys %gids_in_bin;
	my @genes_query_bin = grep {exists $gids_in_bin{$_}} keys %query_list;
	my $queryitem = @genes_query_bin;
	$background{bins}{$bin}{queryitem} = $queryitem;
	$background{bins}{$bin}{querytotal} = $querytotal;
	$background{bins}{$bin}{bgitem} = $bgitem;
	$background{bins}{$bin}{bgtotal} = $bgtotal;
	$background{bins}{$bin}{query_entries} = \@genes_query_bin;
	push @pvalues, hypergeometric_test('phyper', $queryitem, $querytotal , $bgitem, $bgtotal);
}
# p qvalue(\@pvalues);
my @qvalues = @{qvalue(\@pvalues)};

open my $out_fh, ">", $output_file
	or die "cannnot open $output_file:$!\n";

print $out_fh "Mapman_acc\tTerm\tqueryitem\tquerytotal\tbgitem\tbgtotal\tpvalue\tFDR\tentries\n";
foreach my $i (0 .. $#bins) {
	my $bin = $bins[$i];
	my $name = $background{bins}{$bin}{name};
	my $queryitem = $background{bins}{$bin}{queryitem};
	my $querytotal = $background{bins}{$bin}{querytotal};
	my $bgitem = $background{bins}{$bin}{bgitem};
	my $bgtotal = $background{bins}{$bin}{bgtotal};
	my $pvalue = $pvalues[$i];
	my $qvalue = $qvalues[$i];
	my $entries = join "|", @{$background{bins}{$bin}{query_entries}};
	print $out_fh "$bin\t$name\t$queryitem\t$querytotal\t$bgitem\t$bgtotal\t$pvalue\t$qvalue\t$entries\n";
}


sub read_list{
	my ($list_file, $gids_rf) = @_;
	my %gids = %$gids_rf;
	open my $fh, "<", $list_file or die "Cannot open file:$list_file$!\n";
	my %list = ();
	while (<$fh>) {
		chomp;
		$list{$_} = 1 foreach grep {$gids{$_}} map {lc} split /\s+/;
	}
	return %list;
}

sub read_gmt_annotation {
	my $annot_file = shift;
	open my $fh, "<", $annot_file or die "Cannot open file:$!\n";
	readline $fh;
	my %data;
	while (<$fh>) {
		chomp;
		my ($bin_name, @gids) = split /\t/;
		next if @gids < $node_size;
		my @tmp = split /\s+/, $bin_name;
		my $bin = shift @tmp;
		my $name = join " ", @tmp;
		$data{bins}{$bin}{name} = $name;
		my %gids_hash = map {$_ => 1} @gids;
		$data{bins}{$bin}{gids} = \%gids_hash;
		foreach my $gid(@gids) {
			$data{gids}{$gid}{$bin} = 1;
		}
	}
	return %data;
}
	


sub hypergeometric_test {
	use Statistics::R;
	my ($mode, $sample_true, $sample_total , $urn_true, $urn_total, $lower_tail) = @_;
	my $urn_false = $urn_total - $urn_true;
	if ($lower_tail && ($lower_tail eq "T" || $lower_tail eq '1') )  {
		$lower_tail = "T"; 
	}else{
		$lower_tail = "F"; 
	}	
	my $R_code;	
	if ($mode eq 'phyper') {
		$R_code = phyper_code($sample_true, $urn_true , $urn_false, $sample_total, $lower_tail);
	}else{
		$R_code = dhyper_code($sample_true, $urn_true , $urn_false, $sample_total);
	}
	
	my $R = Statistics::R->new();
	$R->run($R_code);
	my $result = $R->get('result'); 
	$R->stop;
	return $result;

}

sub phyper_code($$$$$) {
	my ($sample_true, $urn_true , $urn_false, $sample_total, $lower_tail) = @_;
	return <<DHYPER;
result <- phyper($sample_true, $urn_true , $urn_false, $sample_total, lower.tail=$lower_tail)
DHYPER
}

sub dhyper_code($$$$) {
	my ($sample_true, $sample_total , $urn_true, $urn_total) = @_;
	return <<PHYPER;
result <- dhyper($sample_true, $sample_total , $urn_true, $urn_total)
	
PHYPER
}

sub qvalue {
	use Statistics::R;
	my ($pvalues_rf, $methods) = @_;
	$methods = 'BH' unless $methods;
	my @pavlues = @$pvalues_rf;
	my $num = @pavlues;
	my $p_str = join ", ", @pavlues;
	my $R_code = <<BIN;
	p <- c($p_str)
result <- p.adjust(p, method = "$methods", n = $num)
BIN
	my $R = Statistics::R->new();
	$R->run($R_code);
	my $result = $R->get('result');
	$R->stop;
	return $result;
}
