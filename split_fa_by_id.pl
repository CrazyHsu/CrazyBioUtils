#!/usr/bin/perl -w

use strict;
use 5.010;
my $usage = "UASGE:\n  split_fa_by_id.pl file1.fa file2.fa ...\n";
die $usage unless @ARGV;
foreach my $file (@ARGV) {
    open my $file_fh, "<", $file
      or die "cannot open file:$!\n";
    my $out;
    my $out_fh;
    while (<$file_fh>) {
        if (/^>(\S+)\b/) {
            $out = $1 . '.fa';
            open $out_fh, ">", $out or die $!;
            print $out_fh $_;
        }
        else {
            print $out_fh $_;
        }
    }
    close $out_fh;
}
