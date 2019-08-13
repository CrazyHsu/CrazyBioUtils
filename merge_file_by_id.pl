#!/usr/bin/perl -w

use strict;
use 5.010;
use Getopt::Long;

#usage information
my $usage = <<USAGE;
merge_file_by_id V1.0.2, written by corephi
This program is used to merge two file by id. It set a file
as templete, if bfile has the same line with a start of id, 
it append the bfile content to the end of current line; if 
not, the same number of bfile colums '-' will be replaced
to keep the output file have the same columns.
----------------------------------------------------------
More scripts? Join "bioinformatics*CN" QQ group: 276151571
If you have any probel or suggestions about this program, 
please mail to: hpyu\@genetics.ac.cn
----------------------------------------------------------
Usage: merge_file_by_id -g -a afile -b bfile -o outfile

Options:
  -a    afile name
  -b    bfile name
  -g    gene alias file name
  -m    missing value, default "-"
  -o    output filename, default STDOUT
  
Note: 
  make sure your first column is the identifer, and both the
  a file and b file's header has same identifer such as "ID",
  "GID" etc, so that this program can handle the header as
  well.

USAGE
my $afile           = '';
my $bfile           = '';
my $gene_alias_file = '';
my $missing_value   = "-";
my $outfile         = '-';
die $usage
  unless GetOptions(
    "a=s" => \$afile,
    "g:s" => \$gene_alias_file,
    "m:s" => \$missing_value,
    "b=s" => \$bfile,
    "o:s" => \$outfile,
  );
die $usage unless $afile;

die "a and b can not get from STDIN at the same time\n" 
	if $afile && $bfile && $afile eq '-' && $bfile eq '-';
die "a and b can not get from STDIN at the same time\n" 
	unless $afile || $bfile;

my $bfile_fh;
if ($bfile && $bfile ne '-') {
    die "b file does not exists\n" unless -e $bfile;
    open $bfile_fh, "<", $bfile or die "cannot open file $bfile:$!\n";
}else {
    $bfile_fh = *STDIN;
}

#store bfile in to hash
my %content;
my $bcolnum = 0;    #max column number in bfile
while (<$bfile_fh>) {
    next if /^#/;
    s/\r\n/\n/;
    chomp;
    my @lines = split /\t/;

    #check colnum
    if ($bcolnum) {
        warn "$bfile doesn't have same collum number\n" if $bcolnum != @lines;
    }
    else {
        $bcolnum = @lines;
    }

    my $id = shift @lines;
    if ( exists $content{$id} ) {
        foreach my $i ( 0 .. $#lines ) {
            $content{$id}->[$i] .= "|$lines[$i]";
        }
    }
    else {
        $content{$id} = \@lines;
    }
}
warn "Reading $bfile Done!\n";

#read gene alias
my %gene_alias;
if ($gene_alias_file) {
    open my $gene_alias_fh, "<", $gene_alias_file
      or die "cannot openfile $gene_alias_file:$!";
    readline $gene_alias_fh;
    while (<$gene_alias_fh>) {
        chomp;
        my ( $gid, $symbol ) = split /\t/;
        if ( exists $gene_alias{$gid} ) {
            $gene_alias{$gid} .= "|$symbol";
        }
        else {
            $gene_alias{$gid} = $symbol;
        }

    }
    close $gene_alias_fh;
}

#append content to a file
my $afile_fh;
if ($afile && $afile ne '-') {
    die "a file does not exists\n" unless -e $afile;
    open $afile_fh, "<", $afile or die "cannot open file $afile:$!\n";
}else {
    $afile_fh = *STDIN;
}

my $outfile_fh;
if ($outfile && $outfile ne '-') {
	if (-e $outfile) {
		warn "$outfile exists, appended...\n";
		open $outfile_fh, ">>", $outfile or die "cannot open file $outfile:$!\n";
	}else{
	    open $outfile_fh, ">", $outfile or die "cannot open file $outfile:$!\n";

	}
}else {
    $outfile_fh = *STDOUT;
}

my $acolnum = 0;
while (<$afile_fh>) {
    next if /^#/;
    s/\r\n/\n/;
    chomp;
    my @lines = split /\t/;

    #check colnum
    if ($acolnum) {
        die "$afile doesn't have same collum number\n" if $acolnum != @lines;
    }
    else {
        $acolnum = @lines;
    }

    my $id = $lines[0];

    #add gene alias
    if ($gene_alias_file) {
        if ( $id eq "GID" ) {
            push @lines, 'symbol';
        }
        else {
            if ( exists $gene_alias{$id} ) {
                push @lines, $gene_alias{$id};
            }
            else {
                push @lines, $missing_value;
            }
        }
    }

    #add bname
    if ( exists $content{$id} ) {
        my @cols = @{ $content{$id} };
        if ( @cols < $bcolnum - 1 ) {
            push @cols, split '', ($missing_value) x ( $bcolnum - 1 - @cols );
        }

        push @lines, @cols;
    }
    else {
        push @lines, ($missing_value) x ( $bcolnum - 1 );
        warn "No $id information in file $bfile\n";
    }

    #output
    my $line = join "\t", @lines;
    print $outfile_fh $line . "\n";

    # print $outfile_fh join "\t", @lines;
}
warn "Finished\n";

