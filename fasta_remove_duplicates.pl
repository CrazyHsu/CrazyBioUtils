#!/usr/bin/env perl
# https://github.com/shenwei356/bio_scripts

use strict;
use File::Basename;
use Getopt::Long;
use Digest::MD5 'md5_hex';
use BioUtil::Seq;
use BioUtil::Util;

local $| = 1;
$0 = basename($0);
my $usage = <<USAGE;
Remove duplicated fasta records
Usage: $0 [options] [fastafiles...]
Options:
    
   -n   Comparing by header.
   -s   Comparing by sequence.
   -i   Ignore case.
   -l   Output line length. [70]
   -h   Show this help information.
Examples:
    
    fasta_remove_duplicates.pl -s -i seq1.fa seq2.fa > uniq.fa
    fasta_remove_duplicates.pl -n seq*.fa > uniq.fa
    zcat seq.fa.gz | fasta_remove_duplicates.pl -s -i > uniq.fa
    # remove records same header and seqs
    fasta_remove_duplicates.pl -s -n -i seq1.fa > uniq.fa
https://github.com/shenwei356/bio_scripts
USAGE

my $help        = 0;
my $by_head     = 0;
my $by_seq      = 0;
my $ignore_case = 0;
my $linelength  = 70;
GetOptions(
    'help|h' => \$help,
    "n"      => \$by_head,
    "s"      => \$by_seq,
    "i"      => \$ignore_case,
    'l=i'    => \$linelength,
) or die $usage;

die $usage if $help;

# get the file list
my @files = file_list_from_argv(@ARGV);

my $md5s = {};
my ( $sum, $n ) = ( 0, 0 );
my ( $file, $next_seq, $fa, $header, $seq, $target, $md5 ) = (undef) x 7;
for $file (@files) {
    $next_seq = FastaReader($file);
    while ( $fa = &$next_seq() ) {
        ( $header, $seq ) = @$fa;

        if ($by_seq) {    # comparing by seq
            $target = $seq;
            if ($by_head) {    # comparing by head and seq
                $target = $header . $seq;
            }
        }
        elsif ($by_head) {     # comparing by head
            $target = $header;
        }

        $target = lc $target if $ignore_case;
        $md5 = md5_hex($target);

        if ( $$md5s{$md5} == 1 ) {    # duplicates
            $n++;
        }
        else {
            $$md5s{$md5} = 1;
            $sum++;
            print ">$header\n", format_seq( $seq, $linelength );
        }
        print STDERR "\rremove: $n; remain: $sum";
    }
}

print STDERR "\n";