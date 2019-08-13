#!/usr/bin/perl -w

use strict;
use 5.010;
use Getopt::Long;
use File::Basename;

my $usage = <<USAGE;
SYSNOPSIS
fix_cuffmere_gtf.pl -i transcript.gtf > transcripts.fixed.gtf
written by corephi, group:276151571
this program is used to split bamfile by multi hits, splice junction etc.

 Options:
  -i   in.gtf, default STDIN
  -k   keep original attributes, defaut off
  -o   out.gtf, default STDOUT
  
Examples:
fix_cuffmere_gtf.pl cufflinks.gtf > fixed.gtf
fix_cuffmere_gtf.pl < cufflinks.gtf > fixed.gtf
fix_cuffmere_gtf.pl -i cufflinks.gtf -o fixed.gtf
cat in.gtf | fix_cuffmere_gtf.pl >  fixed.gtf

USAGE

my $in_gtf = '-';
my $out_gtf = '-';
my $keep_attribute = 0;

die $usage
  unless GetOptions(
    "i:s" => \$in_gtf,
    "o:s" => \$out_gtf,
    "k" => \$keep_attribute,
  );

my @in_fhs;
if ($in_gtf && $in_gtf ne '-') {
    open my $in_fh, "<", $in_gtf or die "cannot open file $in_gtf:$!\n";
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
if ($out_gtf ne '-') {
    open $out_fh, ">", $out_gtf or die "cannot open file $out_gtf:$!\n";
}
else {
    $out_fh = *STDOUT;
}

#read the gtf info
my %genes;
my %transcripts;
my $gene_idx = 1;
foreach my $in_fh (@in_fhs) {
	while (<$in_fh>) {
        next if /^#/;
        chomp;
		my @lines = split /\t/;
        my $attribute = pop @lines;
		
        my %attribttes;
		while ($attribute =~ /(\S+) "(\S+)";?/g) {
			my ( $key, $value ) = ( $1, $2 );
			# say $1,'|', $2;
            $attribttes{$key} = $value;
		}
        my $gene_id       = $attribttes{gene_id};
        my $transcript_id = $attribttes{transcript_id};
		
		my $new_gid = '';
		if (exists $genes{$gene_id}) {
			$new_gid = $genes{$gene_id}{name};
		}else{
			$new_gid = sprintf "YHP%08d", $gene_idx;	
			$genes{$gene_id}{name} = $new_gid;		
			$genes{$gene_id}{nex_isfm} = 1;		
			$gene_idx++
		}
		
		my $new_tid = '';
		if (exists $transcripts{$transcript_id}) {
			$new_tid = $transcripts{$transcript_id}
		}else{
			$new_tid = "$new_gid.". $genes{$gene_id}{nex_isfm};
			$transcripts{$transcript_id} = $new_tid;
			$genes{$gene_id}{nex_isfm}++		
		}
		my $new_attr = qq{gene_id "$new_gid"; transcript_id "$new_tid";};	
		if ( exists $attribttes{exon_number} ) {
			my $exon_number = $attribttes{exon_number};
			$new_attr .= qq{ exon_number "$exon_number";};
		}
		
		my $line = join "\t", @lines, $new_attr;
		print $out_fh $line, "\n";
	}	
}

close $out_fh;


