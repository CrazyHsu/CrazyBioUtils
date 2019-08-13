#!/usr/bin/perl -w

use strict;
use 5.010;
use Getopt::Long;
use Data::Printer;

#usage information
my $usage = <<USAGE;
get_intron_from_gff3 V1.0, written by corephi
Given a gtf and fasta file, this program will get intron sequcence 
in fasta format, it can also output the count, length and strand of
intron in text.
More scripts? Please join "bioinformatics*CN" QQ group: 276151571
If you have any probel or suggestions about this program, please
mail to: hpyu\@genetics.ac.cn
----------------------------------------------------------
Usage: get_intron_from_gff3.pl -f genome.fa -g gene_models.gtf -o intron.fa

 Options:
  -m      'sequcence' or 'locus'
  -f      a fasta file contains the genome sequence
  -g      gene_models.gff
  -o      output fasta files, default out.fa;
  -h 

USAGE

#parse parameter
my $fa_file  = '';
my $gff_file = "";
my $out_file = "introns.fa";
my $mode = 'sequence';
my $help = 0;
die $usage
  unless GetOptions(
    "f:s" => \$fa_file,
    "m:s" => \$mode,
    "g=s" => \$gff_file,
    "o=s" => \$out_file,
    "h" => \$help,
  );
#check mode
die $usage if $help;
die "Mode must setted to 'sequence' or 'locus'"  
	unless $mode =~ /sequence|locus/;
  
if ($mode eq 'sequence') {
	die "$fa_file doesn't exists\n"  unless $fa_file;	
}
die "$gff_file doesn't exists\n" unless $gff_file;


#find the non-redudant introns
my %introns = read_annotation($gff_file, 'intron');

if ($mode eq 'sequence') {
	use Bio::SeqIO;
	use Bio::Seq;
	use Bio::DB::Fasta;
	my $fa_db   = Bio::DB::Fasta->new($fa_file);
	my $fa_out = Bio::SeqIO->new(
		-file   => ">$out_file",
		-format => 'Fasta'
	);
	foreach my $locus (keys %introns) {
		my ($seqid, $start, $end) = split /:|-|\.\.|,/, $locus;
		my $seqstr   = $fa_db->seq($seqid, $start => $end);
		# my $seqstr   = $fa_db->seq($locus);
		
		my $seq_obj = Bio::Seq->new( -display_id => $locus,
                         -seq => $seqstr);
		$fa_out->write_seq($seq_obj)
	}			
}elsif ($mode eq 'locus') {
	open my $out_fh, ">", $out_file
		or die "Canont create file $out_file:$!\n";
	print $out_fh "Intron\tLength\tStrand\tCount\n";
	foreach my $locus (keys %introns) {
		my $length = $introns{$locus}{length};
		my $count = $introns{$locus}{count};
		my $strand = $introns{$locus}{strand};
		print $out_fh "$locus\t$length\t$strand\t$count\n";
	}	
}else{
	warn "Not supporting mode\n";
}











#===  FUNCTION  ================================================================
#         NAME: get_gene_len_from_gtf
#      PURPOSE: given gtf file, return a hash of gene length
#   PARAMETERS: $gtf_file
#      RETURNS: %lengh: $length{$gid} = $len
#   DESCRIPTON: length of gene is defined as the total length of non-redudant 
#               exon annotated in gtf.
#===============================================================================
sub get_gene_len_from_gtf {
	my $gtf_file = shift;
	my %genes = read_annotation($gtf_file, 'gene');
	my %length = ();
	foreach my $gid (keys %genes) {
		my %exon_loci = ();
		my @exons = ();
		foreach my $tid (keys %{$genes{$gid}}) {
			push @exons,  @{$genes{$gid}{$tid}{exons}};
		}
		foreach my $i (0..$#exons) {
			map {$exon_loci{$_} = 1} $exons[$i][0] .. $exons[$i][1];
		}

		$length{$gid} = scalar keys %exon_loci;  #length
	}
	return %length;
}

#===  FUNCTION  ================================================================
#         NAME: file_check
#      PURPOSE: check wheather the file exists
#   PARAMETERS: $file: string
#      RETURNS: NA
#===============================================================================
sub file_check {
	my $file = shift;
die "File $file does not exists\n" 
	unless -e  $file;
die "File:$file is empty\n" 
	unless -s $file;	
}


#===  FUNCTION  ================================================================
#         NAME: read_annotation
#      PURPOSE: given a file name, return a hash reference
#   PARAMETERS: $filename: a bed or gtf format filename
#				$mode: 'gene', 'transcript', or 'intron'
#      RETURNS: $reference
#  DESCRIPTION: given a file name, return a hash reference
#				if $mode is "gene", return like followings
#					$hash{$gid}{$t_id}{start}
#					$hash{$gid}{$t_id}{end}
#					$hash{$gid}{$t_id}{seqname}
#					$hash{$gid}{$t_id}{strand}
#					$hash{$gid}{$t_id}{exons}
#					$hash{$gid}{$t_id}{introns}
#					$hash{$gid}{$t_id}{start}
#					$hash{$gid}{$t_id}{start}
#					$hash{$gid}{$t_id}{start}
#				if $mode is "gene", return like followings
#					$hash{$t_id}{start}
#					$hash{$t_id}{end}
#					$hash{$t_id}{seqname}
#					$hash{$t_id}{strand}
#					$hash{$t_id}{exons}
#					$hash{$t_id}{introns}
#					$hash{$t_id}{start}
#					$hash{$t_id}{start}
#					$hash{$t_id}{start}
#				if $mode is "intron", return like followings
#			        $hash{$loci}{count}
#			        $hash{$loci}{length}
#			        $hash{$loci}{strand}
#		  Note: none
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================

sub read_annotation {
	my $in_file = shift;
	my $mode = shift;
	my @suffixlist = qw(.bed .gtf);
	use File::Basename;
	my ($name,$path,$suffix) = fileparse($in_file,@suffixlist);
	
	my %reference;

	if ($mode eq "gene") {
		if ($suffix eq '.gtf') {
			%reference = read_gtf($in_file,"gene");
		}elsif ($suffix eq '.bed'){
			%reference = read_bed($in_file,"gene");
		}else{
			die "unsupported format:$suffix\n";
		}
	}elsif ($mode eq 'transcript') {
		if ($suffix eq '.gtf') {
			%reference = read_gtf($in_file,"transcript");
		}elsif ($suffix eq '.bed'){
			%reference = read_bed($in_file,"transcript");
		}else{
			die "unsupported format:$suffix\n";
		}		
	}elsif ($mode eq 'intron') {
		my %transcripts;
		if ($suffix eq '.gtf') {
			%transcripts = read_gtf($in_file,"transcript");
		}elsif ($suffix eq '.bed'){
			%transcripts = read_bed($in_file,"transcript");
		}else{
			die "unsupported format:$suffix\n";
		}
		%reference = intron_db(\%transcripts);
	}else{
		die "$mode mode does not suppoted\n"; 
	}

	#read junction to hash
	return %reference;
	
}


#===  FUNCTION  ================================================================
#         NAME: intron_db
#      PURPOSE: given a hash ref of transcriptome, return a reference of intron hash
#   PARAMETERS: $trancritome_rf: a hash ref of transcriptome
#      RETURNS: \%array
#  DESCRIPTION: given a hash ref of transcriptome, return a reference of intron hash
#			    THe common keys are:
#			        ..introns{$loci}{count}
#			        ..introns{$loci}{length}
#			        ..introns{$loci}{strand}
#		  Note: none
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub intron_db {
	my $trancritome_rf = shift;
	my %trancritome = %$trancritome_rf ;
	my %introns;
	foreach my $tid (keys %trancritome) {
		my $seqname = $trancritome{$tid}{seqname};
		my $strand = $trancritome{$tid}{strand};
		my $intron_rf = $trancritome{$tid}{introns};
		my @intron = @$intron_rf;
		next unless @intron;
		for (my $i = 0; $i < @intron; $i++) {
			my $loci = "$seqname:$intron[$i][0]-$intron[$i][1]";
			if (exists $introns{$loci}) {
				$introns{$loci}{count}++;
			}else{
				$introns{$loci}{count} = 1;
				my $len = $intron[$i][1]- $intron[$i][0] + 1;
				$introns{$loci}{length} = $len;
				$introns{$loci}{strand} = $strand;
			}
		}
	}
	
	return %introns;
}


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
    my ( $gtf_file, $mode ) = @_;
    die "unknon mode:$mode, it must be one of 'gene' or 'transcript'\n"
        unless ( $mode eq 'gene' or $mode eq 'transcript' );
    my %transcripts;
    open my $gtf_fh, "<", $gtf_file or die "Cannot open gtf $gtf_file:$!\n";

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
						die "same transcript with different $tab, pleae check $transcript_id\n" 
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
						die "same transcript with different $tab, pleae check $transcript_id\n" 
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

				my @tabs = qw(feature start end);
				foreach my $tab (@tabs) {
					if ( exists $transcripts{$gene_id}{$transcript_id}{$tab} ) {
						die "same transcript with different $tab, pleae check $transcript_id\n" 
							if $transcripts{$gene_id}{$transcript_id}{$tab} ne $transcript_info{$tab};
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
			my $attribute = "transcript_id \"$tid\"; gene_id \"$gid\";";
			
			my %transcript_info = (
				feature => 'transcript',
				start   => $t_start, 
				end     => $t_end,
				score   => '.',
				frame   => '.',
				attribute   => $attribute,
				);

			my @tabs = qw(feature start end score frame);
			foreach my $tab (@tabs) {
				if ( exists $transcripts{$transcript_id}{$tab} ) {
					warn "same transcript with different $tab, pleae check $transcript_id\n" 
						if $transcripts{$transcript_id}{$tab} ne $transcript_info{$tab};
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
	my ($hash_rf, $mode, $fh) = @_;
	my %hash = %$hash_rf;
	$mode = 'gene' unless $mode;
	$fh = *STDOUT unless $fh;
	if ($mode eq 'gene') {
		foreach my $gid (sort keys %hash) {
			foreach my $tid (sort keys %{$hash{$gid}}) {
				write_transcript_gtf($hash{$gid}{$tid}, $fh);
			}
		}
	}else{
		foreach my $tid (sort keys %hash) {
			write_transcript_gtf($hash{$tid}, $fh);
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
	my ($transcrtipt_rf, $fh) = @_;
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
	my $attribute = $transcrtipt{attribute};
	my $gene_id = $transcrtipt{gene_id};
	my $transcript_id = $transcrtipt{transcript_id};
	my $exon_rf = $transcrtipt{exons};
	my $intron_rf = $transcrtipt{introns};

	#output the transcript
	my $line = join "\t", ($seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attribute);
	print $fh $line, "\n";
	
	#output the exons
	my @exons = @$exon_rf;
	foreach my $num (0.. $#exons) {
		my $exon_start = $exons[$num][0];
		my $exon_end = $exons[$num][1];
		$line = join "\t", ($seqname, $source, 'exon', $exon_start, $exon_end, $score, $strand, $frame, $attribute);
		print $fh $line, "\n";
	}

}

#===  FUNCTION  ================================================================
#         NAME: read_bed
#      PURPOSE: given a gtf file, retreive the information into a hash
#   PARAMETERS: $gtf_file: string 
#   			$mode: only one of 'gene' or 'transcript', default 'transcript';
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
#			..{$transcript_id}{exons}   #2-D array
#			..{$transcript_id}{introns} #2-D array
#  		
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub read_bed {
	my ($bed_file, $mode) = @_;
	$mode = 'transcript' unless $mode;
	die "unknon mode:$mode, it must be one of 'gene' or 'transcript'\n"
		unless ($mode eq 'gene' or $mode eq 'transcript');
	warn "bed does not suportd gene read mode, guess the gid by tid\n" if $mode eq 'gene';
	open my $bed_fh, "<", $bed_file or die "Cannot open file: $bed_file\n";
	my %transcripts;
	while (<$bed_fh>) {
		chomp;
		next if /^track/;
		my ($chrom,$start,$end,$name,$score,$strand,$thickStart,$thickEnd,$itemRgb,$blockCount,$blockSizes,$blockStarts) = split /\t/;
		my @blockSize = split ',', $blockSizes;
		my @blockStart = split ',', $blockStarts;
		die "uneque number of blockSise and blockStarts\n" unless @blockSize == @blockStart;
		my %loci;
		for (my $i = 0; $i < @blockSize; $i++) {
			my $exon_start = $start + $blockStart[$i] + 1;  #1-based
			my $exon_end = $start + $blockStart[$i] + $blockSize[$i];
			$loci{$exon_start} = $exon_end;
		}
		my ($exon_rf, $intron_rf) = get_intron(\%loci, $strand);
		my $gene_id = get_gid($name);
		my $attribute = "transcript_id \"$name\"; gene_id \"$gene_id\";";
		$transcripts{$name}{seqname} = $chrom;
		$transcripts{$name}{source} = 'bed';
		$transcripts{$name}{feature} = 'transcript';
		$transcripts{$name}{start} = $start + 1;
		$transcripts{$name}{end} = $end;
		$transcripts{$name}{score} = $score;
		$transcripts{$name}{strand} = $strand;
		$transcripts{$name}{frame} = '.';
		$transcripts{$name}{attribute} = $attribute;
		$transcripts{$name}{gene_id} = $gene_id;
		$transcripts{$name}{transcript_id} = $name;
		$transcripts{$name}{exons} = $exon_rf;
		$transcripts{$name}{introns} = $intron_rf;

	}
	close $bed_fh;
	return %transcripts;
}



#===  FUNCTION  ================================================================
#         NAME: get_gid
#      PURPOSE: give a transcirpt id return it's gene id and isoform number
#   PARAMETERS: $tid: string,transcirpt id 
#      RETURNS: ($gid, $num): an array 
#  DESCRIPTION: 
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub get_gid {
	my $tid = shift;
	warn "un assigned $tid" unless $tid;
	my @ids = split /\./, $tid;
	my $gid = shift @ids;
	return $gid;
}



