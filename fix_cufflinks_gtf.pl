#!/usr/bin/perl -w

use strict;
use 5.010;
use Getopt::Long;
use File::Basename;

my $usage = <<USAGE;
SYSNOPSIS
fix_cufflinks_gtf.pl -i transcript.gtf > transcripts.fixed.gtf
written by corephi, group:276151571
this program is used to split bamfile by multi hits, splice junction etc.

 Options:
  -i   in.gtf, default STDIN
  -k   keep original attributes, defaut off
  -o   out.gtf, default STDOUT
  
Examples:
fix_cufflinks_gtf.pl cufflinks.gtf > fixed.gtf
fix_cufflinks_gtf.pl < cufflinks.gtf > fixed.gtf
fix_cufflinks_gtf.pl -i cufflinks.gtf -o fixed.gtf
cat in.gtf | fix_cufflinks_gtf.pl >  fixed.gtf

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

 
foreach my $in_fh (@in_fhs) {
my %genes = read_gtf($in_fh, 'transcript');
write_transcripts_gtf(\%genes, 'transcript', $out_fh, $keep_attribute);
	
}
close $out_fh;



#===  FUNCTION  ================================================================
#         NAME: read_gtf
#      PURPOSE: given a gtf file, retreive the information into a hash
#   PARAMETERS: $gtf_file: string
#   			$mode: only one of 'gene' or 'transcript'
#      RETURNS: $hash_rf: a hash reference stored the information
#  DESCRIPTION: There are two mode of this function, one is 'gene', and the other is 'transcript'. The diffierence bewteen the two mode is the hash structure.
#			'gene' mode: $trasncripts{$gene_id}{$trasncript_id}....
# 			'transcript' mode: $transcripts{$transcript_id}, but have one more key than gene mode: $transcripts{$transcript_id}{gene_id}.
#			THe common keys are:
#			..{$transcript_id}{seqname}
#			..{$transcript_id}{source}
#			..{$transcript_id}{feature}
#			..{$transcript_id}{start}
#			..{$transcript_id}{end}
#			..{$transcript_id}{score}
#			..{$transcript_id}{strand}
#			..{$transcript_id}{frame}
#			..{$transcript_id}{attribute}
#----------------------------------------------------------
#			..{$transcript_id}{gene_id}
#			..{$transcript_id}{transcript_id}
#			..{$transcript_id}{exons}  #2-D array
#			..{$transcript_id}{introns}  #2-D array
#
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub read_gtf {
    my ( $gtf_fh, $mode ) = @_;
    die "unknon mode:$mode, it must be one of 'gene' or 'transcript'\n"
        unless ( $mode eq 'gene' or $mode eq 'transcript' );
    my %transcripts;

    #read the gtf info
    while (<$gtf_fh>) {
        next if /^#/;
        chomp;
        my ($seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attribute)     = split /\t/;
		my %transcript_info = (seqname => $seqname,
			source => $source,
			feature => $feature,
			start => $start,
			end => $end,
			score => $score,
			strand => $strand,
			frame => $frame,
			attribute => $attribute,
			);
		# say $attribute;
		
        my %attribttes;
		while ($attribute =~ /(\S+) "(\S+)";/ig) {
			my ( $key, $value ) = ( $1, $2 );
			# say $1,'|', $2;
            $attribttes{$key} = $value;
		}
        my $gene_id       = $attribttes{gene_id};
        my $transcript_id = $attribttes{transcript_id};
		
			
        if ( $feature eq 'transcript' ) {

            if ( $mode eq 'gene' ) {
				$transcripts{$gene_id}{$transcript_id} = \%transcript_info;
                $transcripts{$gene_id}{$transcript_id}{transcript_id}   = $transcript_id;
                $transcripts{$gene_id}{$transcript_id}{gene_id}   = $gene_id;
            }
            else {
                $transcripts{$transcript_id} = \%transcript_info;				
                $transcripts{$transcript_id}{transcript_id} = $transcript_id;				
                $transcripts{$transcript_id}{gene_id} = $gene_id;
            }
        }
        elsif ( $feature eq 'exon' ) {
            if ( $mode eq 'gene' ) {


				#check and store transcript_id
				if ( exists $transcripts{$gene_id}{$transcript_id}{transcript_id} ) {
					die "same transcript with different transcript_id, pleae check $transcript_id\n" 
						if $transcripts{$gene_id}{$transcript_id}{transcript_id} ne $transcript_id;
				}else {
					$transcripts{$gene_id}{$transcript_id}{transcript_id} = $transcript_id;
				}
				#check and store gene_id
				if ( exists $transcripts{$gene_id}{$transcript_id}{gene_id} ) {
					die "same transcript with different gene_id, pleae check $gene_id\n" 
						if $transcripts{$gene_id}{$transcript_id}{gene_id} ne $gene_id;
				}else {
					$transcripts{$gene_id}{$transcript_id}{gene_id} = $gene_id;
				}				

				#check attributes
				my @tabs = qw(seqname source strand);
				foreach my $tab (@tabs) {
					if ( exists $transcripts{$gene_id}{$transcript_id}{$tab} ) {
						warn "same transcript with different $tab, pleae check $transcript_id\n" 
							if $transcripts{$gene_id}{$transcript_id}{$tab} ne $transcript_info{$tab};
					}else {
						$transcripts{$gene_id}{$transcript_id}{$tab} = $transcript_info{$tab};
					}					
				}
			

                #store exons as hash
                $transcripts{$gene_id}{$transcript_id}{exons}{$start} = $end;
            }
            else {
				#check and store transcript_id
				if ( exists $transcripts{$transcript_id}{transcript_id} ) {
					die "same transcript with different transcript_id, pleae check $transcript_id\n" 
						if $transcripts{$transcript_id}{transcript_id} ne $transcript_id;
				}else {
					$transcripts{$transcript_id}{transcript_id} = $transcript_id;
				}
				#check and store gene_id
				if ( exists $transcripts{$transcript_id}{gene_id} ) {
					die "same transcript with different gene_id, pleae check $gene_id\n" 
						if $transcripts{$transcript_id}{gene_id} ne $gene_id;
				}else {
					$transcripts{$transcript_id}{gene_id} = $gene_id;
				}				

				#check attributes
				my @tabs = qw(seqname source strand);
				foreach my $tab (@tabs) {
					if ( exists $transcripts{$transcript_id}{$tab} ) {
						warn "same transcript with different $tab, pleae check $transcript_id\n" 
							if $transcripts{$transcript_id}{$tab} ne $transcript_info{$tab};
					}else {
						$transcripts{$transcript_id}{$tab} = $transcript_info{$tab};
					}					
				}

                #store exons as hash
                $transcripts{$transcript_id}{exons}{$start} = $end;
            }

     # print "$gene_id\t$transcript_id\t$feature\t$seqname\t$start\t$end\t$strand\n";

        }
        else {
            next;
        }
    }
    close $gtf_fh;

    #get the intron information
    if ( $mode eq 'gene' ) {
        foreach my $gene_id ( sort keys %transcripts ) {
            foreach my $transcript_id ( keys %{ $transcripts{$gene_id} } ) {			
                my $strand = $transcripts{$gene_id}{$transcript_id}{strand};
                my %loci = %{ $transcripts{$gene_id}{$transcript_id}{exons} };

                my ( $exon_rf, $intron_rf ) = get_intron( \%loci, $strand );
                $transcripts{$gene_id}{$transcript_id}{exons}   = $exon_rf;
                $transcripts{$gene_id}{$transcript_id}{introns} = $intron_rf;

				#get the transcript start or end				
				my @temp = sort {$a <=> $b} (@{$exon_rf->[0]}, @{$exon_rf->[-1]});
				my ($t_start, $t_end) = ($temp[0], $temp[-1]);
				
				my $tid = $transcripts{$gene_id}{$transcript_id}{transcript_id};
				my $gid = $transcripts{$gene_id}{$transcript_id}{gene_id};
				my $attribute = "transcript_id \"$tid\"; gene_id \"$gid\";";	
				my %transcript_info = (
					feature => 'transcript',
					start   => $t_start, 
					end     => $t_end,
					score   => '.',
					frame   => '.',
					attribute   => $attribute,
					);

				my @tabs = qw(feature start end score frame attribute);
				foreach my $tab (@tabs) {
					if ( exists $transcripts{$gene_id}{$transcript_id}{$tab} ) {

					}else {
						$transcripts{$gene_id}{$transcript_id}{$tab} = $transcript_info{$tab};
					}					
				}
				

            }
        }
    }
    else {
        foreach my $transcript_id ( sort keys %transcripts ) {
            my $strand = $transcripts{$transcript_id}{strand};
            my %loci   = %{ $transcripts{$transcript_id}{exons} };

            my ( $exon_rf, $intron_rf ) = get_intron( \%loci, $strand );
            $transcripts{$transcript_id}{exons}   = $exon_rf;
            $transcripts{$transcript_id}{introns} = $intron_rf;

			#get the transcript start or end
			my @temp = sort {$a <=> $b} (@{$exon_rf->[0]}, @{$exon_rf->[-1]});
			my ($t_start, $t_end) = ($temp[0], $temp[-1]);
			
			my $tid = $transcripts{$transcript_id}{transcript_id};
            my $gid = $transcripts{$transcript_id}{gene_id};
			warn "no gid in :$tid\n" unless $gid;
			my $attribute = "transcript_id \"$tid\"; gene_id \"$gid\";";
			
			my %transcript_info = (
				feature => 'transcript',
				start   => $t_start, 
				end     => $t_end,
				score   => '.',
				frame   => '.',
				attribute   => $attribute,
				);

			my @tabs = qw(feature start end score frame attribute);
			foreach my $tab (@tabs) {
				if ( exists $transcripts{$transcript_id}{$tab} ) {
				}else {
					$transcripts{$transcript_id}{$tab} = $transcript_info{$tab};
				}					
			}
        }
    }
    return %transcripts;
}

#===  FUNCTION  ================================================================
#         NAME: get_intron
#      PURPOSE: given a hash of exon loci, return a 2d array of exon and intron
#   PARAMETERS: $exon_loci: a hash reference with a format of $hash{$start} = $end
#   			$strand: '+' or '-';
#      RETURNS: @array: [\@exons, \@introns]
#  DESCRIPTION: given a hash of exon loci, return a 2d array of exon and intron.
#		  Note: the order of exon and intron is equal to the strand. that mens, if it is '-', the loci in exon is descreased.
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub get_intron {
    my ( $loci_rf, $strand ) = @_;
	$strand = '+' unless $strand;
    my %loci = %$loci_rf;
    my @exons;
    my @introns;

    # print Dump \%loci;

    my $i = 0;
    foreach my $key ( sort { $a <=> $b } keys %loci ) {
        my $start = $key;
        my $end   = $loci{$key};
        $exons[$i][0] = $start;
        $exons[$i][1] = $end;
        $i++;
    }
    if ( @exons > 1 ) {
        for ( my $i = 0; $i < $#exons; $i++ ) {
            my $start = $exons[$i][1] + 1;
            my $end   = $exons[ $i + 1 ][0] - 1;
            $introns[$i][0] = $start;
            $introns[$i][1] = $end;
        }
    }

    # print Dump \@exons;
    @exons   = reverse @exons   if $strand eq '-';
    @introns = reverse @introns if $strand eq '-';
    return ( \@exons, \@introns );
}
#===  FUNCTION  ================================================================
#         NAME: write_transcripts_gtf
#      PURPOSE: given a transcript or gene hash, write it to gtf files
#   PARAMETERS: $hash_rf: hashreference
#				$file_handle: file handle
#               $verbose: 1 or 0, default 0
#   			$mode: only one of 'gene' or 'transcript'
#  DESCRIPTION: There are two mode of this function, one is 'gene', and the other is 'transcript'. The diffierence bewteen the two mode is the hash structure.
#			'gene' mode: $trasncripts{$gene_id}{$trasncript_id}....
# 			'transcript' mode: $transcripts{$transcript_id}, but have one more key than gene mode: $transcripts{$transcript_id}{gene_id}.
#			THe common keys are:
#			..{$transcript_id}{seqname}
#			..{$transcript_id}{source}
#			..{$transcript_id}{feature}
#			..{$transcript_id}{start}
#			..{$transcript_id}{end}
#			..{$transcript_id}{score}
#			..{$transcript_id}{strand}
#			..{$transcript_id}{frame}
#			..{$transcript_id}{attribute}
#----------------------------------------------------------
#			..{$transcript_id}{gene_id}
#			..{$transcript_id}{transcript_id}
#			..{$transcript_id}{exons}   #2-D array
#			..{$transcript_id}{introns} #2-D array
#
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub write_transcripts_gtf {
	my ($hash_rf, $mode, $fh, $keep_attribute) = @_;
	my %hash = %$hash_rf;
	$mode = 'gene' unless $mode;
	$fh = *STDOUT unless $fh;
	if ($mode eq 'gene') {
		foreach my $gid (sort keys %hash) {
			foreach my $tid (sort keys %{$hash{$gid}}) {
				write_transcript_gtf($hash{$gid}{$tid}, $fh, $keep_attribute);
			}
		}
	}else{
		foreach my $tid (sort keys %hash) {
			write_transcript_gtf($hash{$tid}, $fh, $keep_attribute);
		}		
	}
}
#===  FUNCTION  ================================================================
#         NAME: write_transcript_gtf
#      PURPOSE: given a transcript or gene hash, write it to gtf files
#   PARAMETERS: $hash_rf: hashreference
#				$file_handle: file handle
#  DESCRIPTION: There are two mode of this function, one is 'gene', and the other is 'transcript'. The diffierence bewteen the two mode is the hash structure.
#			'gene' mode: $trasncripts{$gene_id}{$trasncript_id}....
# 			'transcript' mode: $transcripts{$transcript_id}, but have one more key than gene mode: $transcripts{$transcript_id}{gene_id}.
#			THe common keys are:
#			..{$transcript_id}{seqname}
#			..{$transcript_id}{source}
#			..{$transcript_id}{feature}
#			..{$transcript_id}{start}
#			..{$transcript_id}{end}
#			..{$transcript_id}{score}
#			..{$transcript_id}{strand}
#			..{$transcript_id}{frame}
#			..{$transcript_id}{attribute}
#----------------------------------------------------------
#			..{$transcript_id}{gene_id}
#			..{$transcript_id}{transcript_id}
#			..{$transcript_id}{exons}   #2-D array
#			..{$transcript_id}{introns} #2-D array
#
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub write_transcript_gtf{
	my ($transcrtipt_rf, $fh, $keep_attribute) = @_;
	$fh = *STDOUT unless $fh;
	my %transcrtipt = %$transcrtipt_rf;
		
	#output transcripts
	my $seqname = $transcrtipt{seqname};
	my $source = $transcrtipt{source};
	my $feature = $transcrtipt{feature};
	my $start = $transcrtipt{start};
	my $end = $transcrtipt{end};
	my $score = $transcrtipt{score};
	my $strand = $transcrtipt{strand};
	my $frame = $transcrtipt{frame};
	my $gene_id = $transcrtipt{gene_id};
	my $transcript_id = $transcrtipt{transcript_id};
	my $exon_rf = $transcrtipt{exons};
	my $intron_rf = $transcrtipt{introns};
	
	my $attribute = '';
	if ($keep_attribute) {
		$attribute = $transcrtipt{attribute};
	}else{
		$attribute = qq{gene_id "$gene_id"; transcript_id "$transcript_id";};
	}
	#output the transcript
	my $line = join "\t", ($seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attribute);
	print $fh $line, "\n";
	
	#output the exons
	my @exons = @$exon_rf;
	foreach my $num (0.. $#exons) {
		my $exon_start = $exons[$num][0];
		my $exon_end = $exons[$num][1];
		my $exon_number = $num + 1;
		my $exon_attr = $attribute.qq{ exon_number "$exon_number";};
		$line = join "\t", ($seqname, $source, 'exon', $exon_start, $exon_end, $score, $strand, $frame, $exon_attr);
		print $fh $line, "\n";
	}

}


