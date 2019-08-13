#!/usr/bin/perl -w

use strict;
use 5.010;
use YAML;
use Getopt::Long;
use File::Copy;
use File::Basename;
use IPC::System::Simple qw(system capture);

my $usage = <<USAGE;
SYSNOPSIS
blast2go.pl [options] [options] -i example.fa -o example

 Options:
   -i --in_fa       sequences for blast2go in fa format
   -s --go-slim     converts the annotations into the GoSlim version
                    it may be wrong to in convert
   -r --iprscan     <path/to/adirectory/with/ipr_xml_files>
                    Note: not completed
   -p --progress    
   -o --out-folder  output the results to folder,default ./
   -b --b2g4pipe    <path/to/b2g4pipe/>, 
                       default:'~/blast/blast2go/b2g4pipe'
USAGE

my $in_fa      = '';
my $out_folder = '';
my $go_slim    = 0;
my $iprscan    = '';
my $progress   = 1;
my $b2g4pipe   = '/home/hpyu/blast/blast2go/b2g4pipe';
die $usage
  unless GetOptions(
    "i|in_fa=s"      => \$in_fa,
    "o|out-folder:s" => \$out_folder,
    "r|iprscan:s"    => \$iprscan,
    "s|go-slim"      => \$go_slim,
    "p|progress:i"   => \$progress,
    "b|b2g4pipe:s"   => \$b2g4pipe,
  );
print "Start blast2go...\n";
die $usage unless $in_fa;
$b2g4pipe =~ s/\/$//;
$out_folder =~ s/\/$//;


if ($out_folder) {
	if ($out_folder eq $in_fa) {
		$out_folder = $in_fa."_results";
		warn "output folder is equal to fasta filename, redirecting to $out_folder \n" 
	}
}else{
	if ($in_fa =~ /.fasta$|.fa$/) {
		$out_folder = basename ($in_fa, ".fasta", ".fa");
	}else{
		$out_folder = $in_fa."_results";
	}
}

#check directory
if (-e $out_folder && -d $out_folder) {
    die "directory $out_folder is not writable " unless -w "$out_folder";
}
else {
    mkdir $out_folder, 0775 or die "cannot make logs directory:$!";
}

my $command;
#####################################################
#blast
#####################################################
warn "Start blasting";
$command =  "blastx -db nr -query $in_fa -out blastResult.xml -evalue 0.00001 -max_target_seqs 5 -num_threads $progress -outfmt 5 2>blast.log >blast.txt";
warn $command."\n";
eval {
system $command;
};
if ($@){
die  "Something went wrong - $@\n";
}

#####################################################
#blast
#####################################################
warn "Start b2g4pipe...\n";
$command =
"java -Xmx15000m -cp $b2g4pipe/*:$b2g4pipe/ext/*: es.blast2go.prog.B2GAnnotPipe -in blastResult.xml -out $out_folder/results -prop $b2g4pipe/b2gPipe.properties -annot -dat -img -annex -wiki $b2g4pipe/html_template.html 2>b2g4pipe.log >b2g4pipe.txt";
$command .= " -goslim"       if $go_slim;
$command .= " -ips $iprscan" if $iprscan;
warn $command. "\n";
eval { system $command; };
if ($@) {
    die "Something went wrong - $@\n";
}
move("blastResult.xml","$out_folder/blastResult.xml");
warn "Done\n";

