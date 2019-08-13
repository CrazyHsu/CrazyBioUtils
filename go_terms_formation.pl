#!/usr/bin/perl -w

use strict;
use 5.010;
use Getopt::Long;
use File::Basename;

my $usage = <<USAGE;
SYSNOPSIS:
topGO_mapping_formation.pl -i agriGO.map > topGO.map
written by corephi, group:276151571
this program is reformat gene to GO term mapping file

 Options:
  -m   agriGO, WEGO, topGO, default topGO
  -i   agriGO.map, default STDIN
  -o   topGO.map, default STDOUT

File formats:
agriGO.map:
	GID<tab>GO.Terms.1
	GID<tab>GO.Terms.2

topGO.map:
	GID<tab>GO.Terms.1, GO.Terms.2
	
weGO.map:
	GID<tab>GO.Terms.1<tab>GO.Terms.2<tab>GO.Terms.3

  
Examples:
topGO_mapping_formation.pl agriGO.map > topGO.map
topGO_mapping_formation.pl < agriGO.map > topGO.map
topGO_mapping_formation.pl -i agriGO.map -o topGO.map
cat agriGO.map | topGO_mapping_formation.pl >  topGO.map

USAGE

my $in_map = '-';
my $out_map = '-';
my $mode = 'topGO';
my $threashold = 0;
die $usage
  unless GetOptions(
    "i:s" => \$in_map,
    "m:s" => \$mode,
    "o:s" => \$out_map,
  );

my @in_fhs;
if ($in_map && $in_map ne '-') {
    open my $in_fh, "<", $in_map or die "cannot open file $in_map:$!\n";
	push @in_fhs, $in_fh;
}elsif (@ARGV) {
	foreach my $in_file (@ARGV) {
		if (-s $in_file) {
			open my $in_fh, "<", $in_file or die "cannot open file $in_file:$!\n";
			push @in_fhs, $in_fh;		
		}
	}
}else{
	my $in_fh = *STDIN;
	push @in_fhs, $in_fh;
}
 
my $out_fh;
if ($out_map ne '-') {
    open $out_fh, ">", $out_map or die "cannot open file $out_map:$!\n";
}
else {
    $out_fh = *STDOUT;
}

my %gid2go;
#read the map info
foreach my $in_fh (@in_fhs) {
	while (<$in_fh>) {
		s/\r?\n//;
		my ($gid, $go_term_srt) = split /\t/;
		foreach my $go_term(split /,\s*|\t+|\|/, $go_term_srt) {
			$gid2go{$gid}{$go_term} = 1;
		}
	}
	close $in_fh;
}
#output the map info
foreach my $gid (sort keys %gid2go) {
	my @terms = sort keys %{$gid2go{$gid}};
	my $dlim = ", ";
	my $term = '';
	if ($mode eq "WEGO") {
		$dlim = " ";
		$term = join ", ", @terms;
		print $out_fh "$gid\t$term\n";	
	}elsif ($mode eq "topGO"){
		$dlim = ", ";
		$term = join ", ", @terms;
		print $out_fh "$gid\t$term\n";	
	}elsif ($mode eq 'agriGO'){
		@terms = map {"$gid\t$_"} @terms;
		$term = join "\n", @terms;
		print $out_fh "$term\n";	
	}
}

close $out_fh;

