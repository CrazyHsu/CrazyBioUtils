#!/usr/bin/perl -w

use strict;
use 5.010;

use Bio::DB::Sam;
use Getopt::Long;

my $usage = <<USAGE;
Usage: bam_splitter.pl [options...] bamfile
written by corephi, group:276151571
this program is used to split bamfile by multi hits, splice junction etc.
-------------------------------------------------------------------------

 Options:
  -m|--mode             working mode: 'tophat', 'bwa', 'general', default 
                        'tophat'.
  -s|--skip-multi-hits  discard the multiple aligned records, optional,
                        default off.
  -j|--split-junction   split the bam file in two files: one contains the
                        exon records, the the other contains splice juntion
                        records. optional, default off.
  -o|--out-prefix       prefix of output files. optional, default "output".

Note: except the single -s or -j mode, and you can combine -s and -j together
to get the UNIQUE aligned junctions and exons files at the same time.
 
USAGE

#pare the arguments
my $prefix          = "output";
my $skip_multi_hits = 0;
my $split_junction  = 0;
my $mode = 'tophat';
die $usage
  unless GetOptions(
    "m|mode:s"            => \$mode,
    "s|skip-multi-hits" => \$skip_multi_hits,
    "o|out-prefix:s"    => \$prefix,
    "j|split-junction"  => \$split_junction,
  );
die "-s -j are both setted as off, cancelled!\n"
  unless $skip_multi_hits or $split_junction;
if ($mode =~ /tophat|bwa|general/) {
}else{
	warn "unknown mode:$mode, use 'general' instead\n";
	$mode = 'general'
}
my $inbam_file = shift or die $usage;

#open the in bam file
my $inbam = Bio::DB::Bam->open( $inbam_file, "r" );
die $usage unless $inbam;
my $inheader     = $inbam->header;
my $target_names = $inheader->target_name;

#open the out bam file
my ( $multibam, $uniqbam, $sjbam, $exonbam );

if ($skip_multi_hits) {
    $multibam = Bio::DB::Bam->open( $prefix . ".multi.bam", "w" );
    $multibam->header_write($inheader);
    $uniqbam = Bio::DB::Bam->open( $prefix . ".uniq.bam", "w" );
    $uniqbam->header_write($inheader);
}
else {

}

if ($split_junction) {
    $sjbam = Bio::DB::Bam->open( $prefix . ".junction.bam", "w" );
    $sjbam->header_write($inheader);

    $exonbam = Bio::DB::Bam->open( $prefix . ".exon.bam", "w" );
    $exonbam->header_write($inheader);
}
else {

}

my $total_record_num = 0;
my $multi_record_num = 0;
my $sj_record_num = 0;
my $uniq_record_num = 0;
my $exon_record_num = 0;

if ($mode =~ /tophat|bwa/) {
	while ( my $align = $inbam->read1 ) {
		my $seqid = $target_names->[ $align->tid ];
		my $start = $align->pos + 1;
		my $end   = $align->calend;
		my $cigar = $align->cigar_str;

		$total_record_num++;

		if ($skip_multi_hits) {
			my $unique_flag = 0;
			if ($mode eq 'tophat') {
				my $aln_num = 0;
				$aln_num = $align->aux_get("NH");
				if ( $aln_num == 0 ) {
					die "Cannot find any AUX tag NH in the bam, please check\n";
				}
				elsif ( $aln_num == 1 ) {
					$unique_flag = 1;
				}
				else {
					$unique_flag = 0;
				}			
			}
			if ($mode eq 'bwa') {
				my $aln_flag = 0;
				$aln_flag = $align->aux_get("XT");
				if ( $aln_flag eq '0' ) {
					die "Cannot find any AUX tag XT in the bam, please check\n";
				}
				elsif ( $aln_flag eq 'U' ) {
					$unique_flag = 1;
				}
				else {
					$unique_flag = 0;
				}			
			}
			
			if ( $unique_flag == 1 ) {
				$uniqbam->write1($align);
				$uniq_record_num++;
			}
			else {
				$multi_record_num++;
				$multibam->write1($align);
			}

		}

		if ($split_junction) {
			if ( $cigar =~ /(\d+)N/ ) {
				$sjbam->write1($align);
				$sj_record_num++;
			}
			else {
				$exonbam->write1($align);
				$exon_record_num++;
			}
		}
		else {

		}
	}

}else{
	#sort the bam
	warn "Sorting the bam by read name\n";
	Bio::DB::Bam->sort_core(1,$inbam_file,'sorted');
	my $inbam = Bio::DB::Bam->open( "sorted.bam", "r" );
	my $inheader     = $inbam->header;
	my $target_names = $inheader->target_name;
	my @buffer;
	while ( my $align = $inbam->read1 ) {
		my $cigar = $align->cigar_str;
		my $seqname = $align->query->name;
		$total_record_num++;

		#check buffer
		if ($skip_multi_hits) {
			my $unique_flag = 0;
			if (@buffer && $buffer[-1]->query->name ne $seqname) {
				#check multiple alignment;
				my $pair_num = 0;
				my $expected_aligned_num = 1;
				if (@buffer == 1) {
					$unique_flag = 1;
				}elsif (@buffer == 2) {
					my $pair_number = 0;
					my $paired_number = 0;
					foreach my $aln (@buffer) {
						$pair_number += $aln->paired;
					}
					if ($pair_number == 2 && ($buffer[0]->qseq ne $buffer[1]->qseq || $buffer[0]->_qscore ne $buffer[1]->_qscore)) {
						$unique_flag = 1; 
					}else{
						say $buffer[-1]->query->name;
						$unique_flag = 0;
					}
				}else{
					$unique_flag = 0;
				}
				
				if ( $unique_flag) {
						$uniqbam->write1($_) foreach @buffer;
						$uniq_record_num +=	@buffer;
						@buffer = ();
				}else {
					$multi_record_num += @buffer;
					$multibam->write1($_) foreach @buffer;
					@buffer = ();
				}
			}
		}
		
		push @buffer, $align;

		if ($split_junction) {
			if ( $cigar =~ /(\d+)N/ ) {
				$sjbam->write1($align);
				$sj_record_num++;
			}
			else {
				$exonbam->write1($align);
				$exon_record_num++;
			}
		}
		else {

		}
	}
}


print "Total records: $total_record_num\n";
if ($split_junction) {
    print "Splice Junction: $sj_record_num\n";
    print "Exon: $exon_record_num\n";
}
if ($skip_multi_hits) {
    print "Unique records:$uniq_record_num\n";
    print "Multi-hits records:$multi_record_num\n";
}

