#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use Getopt::Long;
use IPC::Cmd qw[can_run run];

my $usage = <<USAGE;
SYSNOPSIS

cuffcmp_rename.pl -r anno.gtf -o out cond.1.gtf cond.2.gtf ...

This program will convert cufflinks gene and transcript ids to human 
readable ids based on the results of cuffcompare. The transcript id
will renamed by its original (from alternative splicing) transcripts, 
and labeled by class code defined by cuffcompare 
This can also used for tracking transcripts 
in different gtf files. 
The same transcript across different gtf files will labeled by the 
same transcript id.

 Options:
   -r|--reference     a gtf reference use for recalculate
   -d|--discard-ref   same to cuffcompre's -R, discared annoation
   -k|--keep-contain  same to cuffcompre's -C, keep contained in
   -s|--seq-dir       genome sequence directory, optional 
   -c|--class-code    class_code used for keep,default all
   -o|--out-prefix    output prefix of the reslut

   
Note: 
If two single exon transcripts overlapped with each other, it will be 
mereged as one. The transcripts with multilpe exons will be treated as 
one, only if the complete match of intron chain.

Pre-requrement:
Before you run this program, it is is highly recommoned to fixed the 
gtf annoation by fix_cufflinks_gtf.pl.
If your annoataion file is gff3, please transform it to gtf by gffread 
imbedded in cufflinks.

----------------------------------------------------------
More scripts? Join "bioinformatics*CN" QQ group: 276151571
If you have any probel or suggestions about this program, 
please mail to: hpyu\@genetics.ac.cn
----------------------------------------------------------

Examplies:

	#reference
	gffread TAIR10.gff3 -T -o TAIR10.gtf	
	fix_cufflinks_gtf.pl < TAIR10.gtf > TAIR10.fix.gtf
	
	#assembled or merged gtf files
	fix_cuffmere_gtf.pl cond.1.gtf | fix_cufflinks_gtf.pl > cond.1.fix.gtf
	fix_cuffmere_gtf.pl cond.2.gtf | fix_cufflinks_gtf.pl > cond.2.fix.gtf
	fix_cuffmere_gtf.pl cond.3.gtf | fix_cufflinks_gtf.pl > cond.3.fix.gtf

	#run cuffcmp_reanme
	cuffcmp_rename.pl -r TAIR10.fix.gtf -o conds cond.1.fix.gtf
	cond.2.fix.gtf cond.3.fix.gtf

Class_code transform:
for the compicity of gff3 formated gene model, two class codes were 
changed:
'=' -> 'y'
'.' -> 'd'

Some transcripts, such as reference annoataion that not present in RNA-seq
data, or transcripts contained in another transcript, will filterd by 
cuffcompare by default, so this will not output in this program. If u want 
to report these transcripts, please check parameter '-k' and '-d'.

USAGE

my $reference_file  = '/home/hpyu/tair_data/TAIR10_MOD.gtf';
my $prefix          = 'out';
my $seq_dir         = '';
my $strand_mode     = 0;
my $discard_ref     = 0;
my $keep_contained  = 0;
my $keep_class_code = '=cjeiopruxs.';
die $usage
  unless GetOptions(
    "r|reference:s"    => \$reference_file,
    "s|seq-dir:s"      => \$seq_dir,
    "d|discard-ref"  => \$discard_ref,
    "k|keep-contain" => \$keep_contained,
    "o|out-prefix:s"   => \$prefix,
    "c|class-code:s"   => \$keep_class_code,
  );

can_run('cuffcompare') or die 'cuffcompare is not installed!';

my @transcriptomes = grep { -s $_ } @ARGV;
die "No gtf files, plase check their existance\n" . $usage
  unless @transcriptomes;
die "$reference_file dose not exists" unless -s $reference_file;
if ($seq_dir) {
    die "$seq_dir dose not exists" unless -d $seq_dir || -s $seq_dir;
}

#########################################################################################
##Prerequest
#########################################################################################
warn "Starting cuffcompare...\n";
my $cuffcompare_param = "-r $reference_file -V ";
$cuffcompare_param .= "-s $seq_dir " if $seq_dir;
$cuffcompare_param .= "-R "          if $discard_ref;
$cuffcompare_param .= "-C "          if $keep_contained;

my $gtfs = join " ", @transcriptomes;
my $command =
"cuffcompare $cuffcompare_param $gtfs >cuffcompare.stdout.txt 2> cuffcompare.stderr.txt";

warn "\t$command\n";
my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
  run( command => $command, verbose => 0 );
if ($success) {
    warn "\tCuffcompare Done!\n";
}
else {
    my @stderrs = @$stderr_buf;
    warn "Something went wrong:\n@stderrs";
}

my $tracking_file = 'cuffcmp.tracking';

#########################################################################################
##Step 1 read gtf and construct loci mapping information
#########################################################################################

my %transcripts_database;
my %transcripts_database_idx;

warn "Reading Reference gtf file:$reference_file\n";
my $transcritpome_ref_db = read_gtf( $reference_file, 'transcript' );
$transcripts_database{reference} = $transcritpome_ref_db;
$transcripts_database_idx{reference} =
  construct_transcriptome_idx($transcritpome_ref_db);

#read gene-level annatioantion files, to check whether gene is anno or novel
my $refernce_genes_db = read_gtf( $reference_file, 'gene' );
my %ref_gene_ids = %$refernce_genes_db;
$ref_gene_ids{$_} = 1 foreach keys %ref_gene_ids;

#read transcript information;
foreach my $transcriptome (@transcriptomes) {
    warn "Reading gtf file:$transcriptome\n";
    my $transcritpome_db = read_gtf( $transcriptome, 'transcript' );
    $transcripts_database{$transcriptome} = $transcritpome_db;
    $transcripts_database_idx{$transcriptome} =
      construct_transcriptome_idx($transcritpome_db);
}

# print Dump \%transcripts_database;
#compare the trasncript information
warn "Constructing hash unitranscripts and id_refs\n";
my %unitranscripts;        #trancript hash same like %transcripts_database
my %unitranscripts_idx;    #trancript hash same like %transcripts_database
my %id_ref;  #id mpping list, locus_id REF_ID trancriptome1_id transcriptome2_id
my $loci_num = 0;
foreach my $transcriptome_a (@transcriptomes) {
    my %transcripts_a = %{ $transcripts_database{$transcriptome_a} };
    foreach my $tid_a ( keys %transcripts_a ) {

        #check if it is new transcirpts;
        my $transcripts_a_rf = $transcripts_a{$tid_a};

        # print Dump $transcripts_a_rf if $tid_a =~ "AT3G02430";
        if ( exists_in_db( $transcripts_a_rf, \%unitranscripts_idx ) ) {
            next;
        }
        else {
            #construct new loci
            $loci_num++;
            my $locus_id = sprintf "LOCI%08d", $loci_num;

            #write to database
            $id_ref{$locus_id}{$transcriptome_a} = $tid_a;
            $unitranscripts{$locus_id} = $transcripts_a_rf;
            add_idx_to_transcriptome( $transcripts_a_rf, \%unitranscripts_idx );

            #compare with reference
            my @ref_ids = exists_in_db( $transcripts_a_rf,
                $transcripts_database_idx{reference} );
            if (@ref_ids) {
                warn "$tid_a has duplicatd refernce id:@ref_ids"
                  if @ref_ids > 1;
                my $ref_id = $ref_ids[0];
                $id_ref{$locus_id}{reference} = $ref_ids[0];
                $unitranscripts{$locus_id} =
                  $transcripts_database{reference}{ $ref_ids[0] };
            }
            else {
                $id_ref{$locus_id}{reference} = '-';
            }

            #comapre with other trancsripts
            foreach my $transcriptome_b (@transcriptomes) {
                next if $transcriptome_a eq $transcriptome_b;
                my @tids = exists_in_db( $transcripts_a_rf,
                    $transcripts_database_idx{$transcriptome_b} );
                if (@tids) {
                    warn "$tid_a duplicatd in $transcriptome_b, id: @ref_ids"
                      if @ref_ids > 1;
                    $id_ref{$locus_id}{$transcriptome_b} = $tids[0];
                }
                else {
                    $id_ref{$locus_id}{$transcriptome_b} = '-';
                }
            }

        }

    }

}

# print Dump \%unitranscripts;
# print Dump \%id_ref;

################################################################################
#Step 2. reads cuffcompare id mapping system
################################################################################

open my $refmap_fh, "<", $tracking_file
  or die "cannot open file $tracking_file:$!\n";
my %transcript_ref;
my %gene_ref;
warn "Reading cuffcompare track file...\n";
while (<$refmap_fh>) {
    chomp;
    my @lines        = split /\t/;
    my $trans_id     = shift @lines;
    my $locus_id     = shift @lines;
    my $reference_id = shift @lines;
    my $class_code   = shift @lines;

    my ( $ref_id, $ref_gene );
    if ( $reference_id eq '-' ) {
        $ref_id = '-';
    }
    else {
        my ( $gene_temp, $transcript_temp ) = split /\|/, $reference_id;
        $ref_id   = $transcript_temp;
        $ref_gene = $gene_temp;
    }

    # say $reference_id;
    foreach my $transcriptome (@transcriptomes) {
        my $information = shift @lines;
        my ( $cuff_id, $cuff_gene );
        if ( $information eq '-' ) {
            $cuff_id   = '-';
            $cuff_gene = '-';
            next;
        }
        else {
            my ( $name, $info ) = split ':', $information;
            my ( $cuff_gene_temp, $cuff_id_temp, ) = split /\|/, $info;
            $cuff_id   = $cuff_id_temp;
            $cuff_gene = $cuff_gene_temp;
        }
        $transcript_ref{$transcriptome}{$cuff_id}{ref_id}     = $ref_id;
        $gene_ref{$transcriptome}{$cuff_gene}{ref_id}         = $ref_gene;
        $transcript_ref{$transcriptome}{$cuff_id}{class_code} = $class_code;
    }
}
close $refmap_fh;

# print Dump \%transcript_ref;

########################################################################################
##Step3 merge the id mappings
########################################################################################
my ( %uni_renamed_id, %isoforms, );
my $novel_gene_num = 1;
foreach my $locus_id ( sort keys %id_ref ) {
    my $ref_id = $id_ref{$locus_id}{reference};
    warn "$locus_id has no ref_id\n" unless $ref_id;

    ##find the unique gene in cufflinks, and get the transcript_ref and classcode;
    my %cuff_ref_ids     = ();
    my %cuff_class_codes = ();
    foreach my $transcriptome (@transcriptomes) {
        my $cuff_id = $id_ref{$locus_id}{$transcriptome};
        if ( $cuff_id && $cuff_id ne '-' ) {
            if ( exists $transcript_ref{$transcriptome}{$cuff_id} ) {
                my $cuff_ref_id =
                  $transcript_ref{$transcriptome}{$cuff_id}{ref_id};

# warn "$cuff_id in $transcriptome has no reference id information\n" unless $cuff_ref_id;
                my $cuff_class_code =
                  $transcript_ref{$transcriptome}{$cuff_id}{class_code};

# warn "$cuff_id in $transcriptome has no class_codeinformation\n" unless $cuff_class_code;
                $cuff_ref_ids{$cuff_ref_id}++;
                $cuff_class_codes{$cuff_class_code}++;
            }
            else {
             # warn "$cuff_id in $transcriptome has no reference information\n";
            }

        }
    }
    my @cuff_ref_ids     = keys %cuff_ref_ids;
    my @cuff_class_codes = keys %cuff_class_codes;
    my ( $cuff_ref_id, $class_code );
    warn "$locus_id has more than one class_codes:cuff_class_codes\n"
      if @cuff_class_codes > 1;
    if (@cuff_ref_ids) {
        warn "$locus_id has more than one ref_genes:@cuff_ref_ids\n"
          if @cuff_ref_ids > 1;
        $cuff_ref_id = shift @cuff_ref_ids;
    }
    else {
        warn "$locus_id has been filtered by cuffcompare, Skipping\n";
        next;
        $cuff_ref_id = '-';
    }
    if (@cuff_class_codes) {
        warn "$locus_id has more than one class_code:@cuff_class_codes\n"
          if @cuff_class_codes > 1;
        $class_code = shift @cuff_class_codes;
    }
    else {
        warn "$locus_id has no class_code\n";
        $class_code = 'n';
    }

    # print "$locus_id\t$ref_id\t$cuff_ref_id\t$class_code";
    #merge the cuff_ref_id and ref_id
    my ( $ref_gene, $isoform_id );
    #################ref_id####################################
    ##find the correct $ref_id
    if ( $ref_id eq '-' ) {
        if ( $cuff_ref_id eq '-' ) {

            # both of them has no annotaion
            my %gene_refs = ();

            #check wheather it is
            foreach my $transcriptome (@transcriptomes) {
                my $cuff_id = $id_ref{$locus_id}{$transcriptome};
                my ( $cuff_gene, ) = split_tid($cuff_id);

                # print "\t$cuff_gene";
                if ( exists $gene_ref{$transcriptome}{$cuff_gene} ) {
                    my $ref_gene_temp =
                      $gene_ref{$transcriptome}{$cuff_gene}{ref_id};
                    $gene_refs{$ref_gene_temp}++
                      if $ref_gene_temp && $ref_gene_temp ne '-';
                }
                else {
# warn "$cuff_id in $transcriptome has no gene mapping information information\n";
                }
            }

            # print Dump \%gene_refs;
            my @gene_refs = keys %gene_refs;

            # print "\t",@gene_refs;
            if ( @gene_refs == 0 ) {
                $ref_gene = sprintf "NOVEL%06d", $novel_gene_num;
                $novel_gene_num++;
            }
            elsif ( @gene_refs == 1 ) {
                $ref_gene = shift @gene_refs;
            }
            else {
                warn "$locus_id has more than one ref_gene\n";
                $ref_gene = shift @gene_refs;
            }

            #if it is new gene, so it must be a new transcript, renmae it.
            $ref_id = "$ref_gene.1";

        }
        else {
            # cuff_id has annotation
            $ref_id = $cuff_ref_id;
        }
    }
    else {
        if ( $cuff_ref_id eq '-' ) {
            warn "locus\t$locus_id has no cuff_ref_id\n";
        }
        else {
            warn
"locus(${locus_id})'s refernce id is not the same:$ref_id\t$cuff_ref_id\n"
              if ( $ref_id ne $cuff_ref_id );
            $ref_id = $cuff_ref_id;
        }
    }

    ################ref_gene####################################
    ##given the ref_id, got the ref_gene
    ( $ref_gene, ) = split_tid($ref_id);

    #write $ref_id and $ref_gene back to database;
    foreach my $transcriptome (@transcriptomes) {
        my $cuff_id = $id_ref{$locus_id}{$transcriptome};

        # say $cuff_id;
        if ( $cuff_id && $cuff_id ne '-' ) {
            my ( $cuff_gene, ) = split_tid($cuff_id);
            $transcript_ref{$transcriptome}{$cuff_id}{ref_id}     = $ref_id;
            $transcript_ref{$transcriptome}{$cuff_id}{class_code} = $class_code;
            $gene_ref{$transcriptome}{$cuff_gene}{ref_id}         = $ref_gene;
        }
        else {
            next;
        }

    }
    $id_ref{$locus_id}{reference} = $ref_id;

    #################isoform_id####################################
    #the followling start to rename the transcript_id;
    if ( $ref_id =~ /CUFF/ ) {

        #check wheather it is new novel gene;
        warn "error ref_id:$ref_id\n";
    }
    else {
        # get the  $ref_gene;
        ( $ref_gene, ) = split_tid($ref_id);

        #compatible to old version of cufflinks

        #get the isoform_id
        my $num = 1;
        $isoform_id = "$ref_id.$class_code.$num";

        # print "reanmed transcript_id:$transcript_id\n";
        while ( exists $uni_renamed_id{$isoform_id} ) {
            $num++;
            $isoform_id = "$ref_id.$class_code.$num";
        }
        $uni_renamed_id{$isoform_id} = $locus_id;
    }

    foreach my $transcriptome (@transcriptomes) {
        my $cuff_id = $id_ref{$locus_id}{$transcriptome};
        if ( $cuff_id && $cuff_id ne '-' ) {
            my ( $cuff_gene, ) = split_tid($cuff_id);

            #store isoform
            $isoforms{$transcriptome}{$cuff_id}{ref_id}     = $ref_id;
            $isoforms{$transcriptome}{$cuff_id}{ref_gene}   = $ref_gene;
            $isoforms{$transcriptome}{$cuff_id}{class_code} = $class_code;
            $isoforms{$transcriptome}{$cuff_id}{cuff_gene}  = $cuff_gene;
            $isoforms{$transcriptome}{$cuff_id}{isoform_id} = $isoform_id;
        }
        else {
            next;
        }
    }

    # print "\t$ref_gene\n";
}

# print Dump \%isoforms;

########################################################################################
##Step4 write to new file
########################################################################################

foreach my $transcriptome (@transcriptomes) {
    warn "Renaming $transcriptome...\n";
    open my $gtf_fh, "<", $transcriptome
      or die "cannot open file $transcriptome:$!\n";
    my %genes;
    open my $anno_fh,   ">", "${prefix}.anno.$transcriptome";
    open my $filter_fh, ">", "${prefix}.filtered.$transcriptome";
    open my $novel_fh,  ">", "${prefix}.novel.$transcriptome";
    my %novel;
    my %anno;

    while (<$gtf_fh>) {
        chomp;
        my $line = $_;
        $line =~ s/;$//;
        my @lines     = split /\t/, $line;
        my $attribute = pop @lines;
        my $type      = $lines[2];
        $lines[1] = "YHP";
        my @temp = grep { $_ }
          map { s/^\s+//; $_ }
          split ';', $attribute;
        my %attributes;

        foreach my $attr (@temp) {
            my ( $key, $value ) = split /\s+/, $attr;
            $value =~ s/"//g;

            # say $key,"\t","$value";
            $attributes{$key} = $value;
        }

        # print Dump \%attributes;
        my $gene_id       = $attributes{gene_id};
        my $transcript_id = $attributes{transcript_id};
        my $exon_number   = '';
        $exon_number = $attributes{exon_number}
          if exists $attributes{exon_number};

    #reanme gene_id ============================================================
        if ( exists $isoforms{$transcriptome}{$transcript_id} ) {

            #rename gene_id
            $gene_id = $isoforms{$transcriptome}{$transcript_id}{ref_gene};

            #reanme transcript id;
            $transcript_id =
              $isoforms{$transcriptome}{$transcript_id}{isoform_id};
            my ( $gid, $num, $class_code, $gnum ) = split /\./, $transcript_id;
            $class_code = '\.' unless $class_code;

            # say $keep_class_code;
            if ( $keep_class_code =~ /$class_code/ ) {

            }
            else {
                warn "Skipping $transcript_id, class_code:$class_code\n"
                  if $type eq "transcript";
                next;
            }

        }
        else {
            my @tids = exists_in_db(
                $transcripts_database{$transcriptome}{$transcript_id},
                $transcripts_database_idx{$transcriptome}
            );
            if ( @tids > 1 ) {
                warn "Skiping redundant transcript:$transcript_id\n"
                  if $type eq "transcript";
                say $filter_fh $line . ';Filter "Redundant"';
                next;
            }
            else {
                warn "Skipping Cuffcompare filtered transcript:$transcript_id\n"
                  if $type eq "transcript";
                say $filter_fh $line . ';Filter "Cuffcompare"';
                next;
            }
        }

    #reanme gene_id ========================done================================

        #output to new file
        if ( exists $ref_gene_ids{$gene_id} ) {    #check if it is novel gene
            if ( $type eq "transcript" ) {         #check if it is transcript
                if ( exists $genes{$gene_id} ) {
                    $genes{$gene_id}++;
                }
                else {
                    $genes{$gene_id} = 1;
                }
            }
            else {

            }

            #renmae attibute
            my $modifid_transcript_id = reformat_tid($transcript_id);
            $attribute =
              "gene_id \"$gene_id\"; transcript_id \"$modifid_transcript_id\";";
            $attribute .= " exon_number \"$exon_number\";"
              if $type eq "exon" && $exon_number;

            #output
            push @lines, $attribute;
            say $anno_fh join "\t", @lines;

            #statistics
            if ( $type eq 'transcript' ) {
                if (   exists $anno{$gene_id}
                    && exists $anno{$gene_id}{$transcript_id} )
                {
                    die "Duplicated $gene_id:$transcript_id\n";
                }
                else {
                    $anno{$gene_id}{$transcript_id} = 0;
                }
            }
            else {
                $anno{$gene_id}{$transcript_id}++;
            }

        }
        else {
            #output
            my $modifid_transcript_id = reformat_tid($transcript_id);
            $attribute =
              "gene_id \"$gene_id\"; transcript_id \"$modifid_transcript_id\";";
            push @lines, $attribute;
            say $novel_fh join "\t", @lines;

            # say $novel_fh $_;

            if ( $type eq 'transcript' ) {
                if (   exists $novel{$gene_id}
                    && exists $novel{$gene_id}{$transcript_id} )
                {
                    die "Duplicated $gene_id:$transcript_id\n";
                }
                else {
                    $novel{$gene_id}{$transcript_id} = 0;
                }
            }
            else {
                $novel{$gene_id}{$transcript_id}++;
            }
        }
    }
    close $gtf_fh;
    close $filter_fh;
    close $anno_fh;
    close $novel_fh;

    my $novel_gene_total       = keys %novel;
    my $novel_transcript_total = 0;
    for my $gene ( sort keys %novel ) {
        my $trans_count = keys %{ $novel{$gene} };
        $novel_transcript_total += $trans_count;
    }

    print "$transcriptome:\n";
    print "Novel gene number:$novel_gene_total\n";
    print "Novel transcript number:$novel_transcript_total\n";

    # print Dump \%anno;

    my $anno_gene_total = keys %anno;
    my $anno_transcript_total;
    my $as_gene_count;
    my $intron_gene_count;
    for my $gene ( sort keys %anno ) {
        my $trans_count = keys %{ $anno{$gene} };
        $anno_transcript_total += $trans_count;

        #check wheather it is alternitvie splied gene

        $as_gene_count++ if $trans_count > 1;

        #check whether it is intron-cotained gene
        my $intron_flag = 0;
        my $exon_count;
        foreach my $trans ( sort keys %{ $anno{$gene} } ) {
            $intron_flag = 1 if $anno{$gene}{$trans} > 1;
            $exon_count += $anno{$gene}{$trans};
        }
        $intron_gene_count++ if $intron_flag;
    }
    print "Annotated gene number:$anno_gene_total\n";
    print
"Annotated gene's transcript number(including novel):$anno_transcript_total\n";
    print "Winthin the Annotated Genes\n";
    print "Intron-containing Genes:$intron_gene_count\n";
    print "Alternative Spliced Genes:$as_gene_count\n";
    printf "Ratio:%2.2f%%\n", $as_gene_count / $intron_gene_count * 100;

    # print Dump \%isoforms;
}

warn "Done!\n";

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
        my (
            $seqname, $source, $feature, $start, $end,
            $score,   $strand, $frame,   $attribute
        ) = split /\t/;
        my %transcript_info = (
            seqname   => $seqname,
            source    => $source,
            feature   => $feature,
            start     => $start,
            end       => $end,
            score     => $score,
            strand    => $strand,
            frame     => $frame,
            attribute => $attribute,
        );

        # say $attribute;

        my %attribttes;
        while ( $attribute =~ /(\S+) "(\S+)";?/g ) {
            my ( $key, $value ) = ( $1, $2 );

            # say $1,'|', $2;
            $attribttes{$key} = $value;
        }
        my $gene_id       = $attribttes{gene_id};
        my $transcript_id = $attribttes{transcript_id};

        if ( $feature eq 'transcript' ) {

            if ( $mode eq 'gene' ) {
                $transcripts{$gene_id}{$transcript_id} = \%transcript_info;
                $transcripts{$gene_id}{$transcript_id}{transcript_id} =
                  $transcript_id;
                $transcripts{$gene_id}{$transcript_id}{gene_id} = $gene_id;
            }
            else {
                $transcripts{$transcript_id}                = \%transcript_info;
                $transcripts{$transcript_id}{transcript_id} = $transcript_id;
                $transcripts{$transcript_id}{gene_id}       = $gene_id;
            }
        }
        elsif ( $feature eq 'exon' ) {
            if ( $mode eq 'gene' ) {

                #check and store transcript_id
                if (
                    exists $transcripts{$gene_id}{$transcript_id}{transcript_id}
                  )
                {
                    die
"same transcript with different transcript_id, pleae check $transcript_id\n"
                      if $transcripts{$gene_id}{$transcript_id}{transcript_id}
                      ne $transcript_id;
                }
                else {
                    $transcripts{$gene_id}{$transcript_id}{transcript_id} =
                      $transcript_id;
                }

                #check and store gene_id
                if ( exists $transcripts{$gene_id}{$transcript_id}{gene_id} ) {
                    die
"same transcript with different gene_id, pleae check $gene_id\n"
                      if $transcripts{$gene_id}{$transcript_id}{gene_id} ne
                      $gene_id;
                }
                else {
                    $transcripts{$gene_id}{$transcript_id}{gene_id} = $gene_id;
                }

                #check attributes
                my @tabs = qw(seqname source strand);
                foreach my $tab (@tabs) {
                    if ( exists $transcripts{$gene_id}{$transcript_id}{$tab} ) {
                        warn
"same transcript with different $tab, pleae check $transcript_id\n"
                          if $transcripts{$gene_id}{$transcript_id}{$tab} ne
                          $transcript_info{$tab};
                    }
                    else {
                        $transcripts{$gene_id}{$transcript_id}{$tab} =
                          $transcript_info{$tab};
                    }
                }

                #store exons as hash
                $transcripts{$gene_id}{$transcript_id}{exons}{$start} = $end;
            }
            else {
                #check and store transcript_id
                if ( exists $transcripts{$transcript_id}{transcript_id} ) {
                    die
"same transcript with different transcript_id, pleae check $transcript_id\n"
                      if $transcripts{$transcript_id}{transcript_id} ne
                      $transcript_id;
                }
                else {
                    $transcripts{$transcript_id}{transcript_id} =
                      $transcript_id;
                }

                #check and store gene_id
                if ( exists $transcripts{$transcript_id}{gene_id} ) {
                    die
"same transcript with different gene_id, pleae check $gene_id\n"
                      if $transcripts{$transcript_id}{gene_id} ne $gene_id;
                }
                else {
                    $transcripts{$transcript_id}{gene_id} = $gene_id;
                }

                #check attributes
                my @tabs = qw(seqname source strand);
                foreach my $tab (@tabs) {
                    if ( exists $transcripts{$transcript_id}{$tab} ) {
                        warn
"same transcript with different $tab, pleae check $transcript_id\n"
                          if $transcripts{$transcript_id}{$tab} ne
                          $transcript_info{$tab};
                    }
                    else {
                        $transcripts{$transcript_id}{$tab} =
                          $transcript_info{$tab};
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
                my %loci   = %{ $transcripts{$gene_id}{$transcript_id}{exons} };

                my ( $exon_rf, $intron_rf ) = get_intron( \%loci, $strand );
                $transcripts{$gene_id}{$transcript_id}{exons}   = $exon_rf;
                $transcripts{$gene_id}{$transcript_id}{introns} = $intron_rf;

                #get the transcript start or end
                my @temp = sort { $a <=> $b }
                  ( @{ $exon_rf->[0] }, @{ $exon_rf->[-1] } );
                my ( $t_start, $t_end ) = ( $temp[0], $temp[-1] );

                my $tid = $transcripts{$gene_id}{$transcript_id}{transcript_id};
                my $gid = $transcripts{$gene_id}{$transcript_id}{gene_id};
                my $attribute = "transcript_id \"$tid\"; gene_id \"$gid\";";

                my %transcript_info = (
                    feature   => 'transcript',
                    start     => $t_start,
                    end       => $t_end,
                    score     => '.',
                    frame     => '.',
                    attribute => $attribute,
                );

                my @tabs = qw(feature start end score frame attribute);
                foreach my $tab (@tabs) {
                    if ( exists $transcripts{$gene_id}{$transcript_id}{$tab} ) {

                    }
                    else {
                        $transcripts{$gene_id}{$transcript_id}{$tab} =
                          $transcript_info{$tab};
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
            my @temp =
              sort { $a <=> $b } ( @{ $exon_rf->[0] }, @{ $exon_rf->[-1] } );
            my ( $t_start, $t_end ) = ( $temp[0], $temp[-1] );

            my $tid       = $transcripts{$transcript_id}{transcript_id};
            my $gid       = $transcripts{$transcript_id}{gene_id};
            my $attribute = "transcript_id \"$tid\"; gene_id \"$gid\";";

            my %transcript_info = (
                feature   => 'transcript',
                start     => $t_start,
                end       => $t_end,
                score     => '.',
                frame     => '.',
                attribute => $attribute,
            );

            my @tabs = qw(feature start end score frame attribute);
            foreach my $tab (@tabs) {
                if ( exists $transcripts{$transcript_id}{$tab} ) {

                }
                else {
                    $transcripts{$transcript_id}{$tab} = $transcript_info{$tab};
                }
            }
        }
    }
    return \%transcripts;
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
        for ( my $i = 0 ; $i < $#exons ; $i++ ) {
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
#         NAME: split_tid
#      PURPOSE: give a transcirpt id return it's gene id and isoform number
#   PARAMETERS: $tid: string,transcirpt id
#      RETURNS: ($gid, $num): an array
#  DESCRIPTION:
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub split_tid {
    my $tid = shift;
    warn "un assigned $tid" unless $tid;
    my @ids = split /\./, $tid;
    my $num = pop @ids;
    my $gid = join ".", @ids;
    return ( $gid, $num );
}

#===  FUNCTION  ================================================================
#         NAME: construct_transcript_idx{
#      PURPOSE: give a transcript hash return a hashref of transcript index
#               as: $transcriptome_data{$intron_chain} = $transcript_rf
#   PARAMETERS: $transcirpt  hash of transcript
#      RETURNS: $intron_chain, string
#  DESCRIPTION:
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub construct_transcript_idx {
    my $transcript_rf = shift;
    my %transcript_idx;

    my $tid    = $transcript_rf->{transcript_id};
    my $chr    = $transcript_rf->{seqname};
    my $strand = $transcript_rf->{strand};

    my @exons   = @{ $transcript_rf->{exons} };
    my @introns = @{ $transcript_rf->{introns} };
    if (@introns) {
        my @intron_strs = ();
        foreach my $i ( 0 .. $#introns ) {
            push @intron_strs, join "-", @{ $introns[$i] };
        }
        my $intron_string = join ",", @intron_strs;
        my $intron_chain = "$chr:" . $intron_string . "/$strand";
        $transcript_idx{introns}{$intron_chain} = $transcript_rf;
    }
    else {
        foreach my $pos (@exons) {
            my $locus = "$chr:$pos/$strand";
            $transcript_idx{exons}{$locus} = $transcript_rf;
        }
    }
    return \%transcript_idx;
}

#===  FUNCTION  ================================================================
#         NAME: add_idx_to_transcriptome{
#      PURPOSE:
#   PARAMETERS: $transcirpt  hash of transcript
#      RETURNS: $intron_chain, string
#  DESCRIPTION:
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub add_idx_to_transcriptome {
    my $transcript_rf    = shift;
    my $transcriptome_rf = shift;

    my $transcript_idx_rf = construct_transcript_idx($transcript_rf);

    my @introns = $transcript_rf->{introns};
    if (@introns) {
        foreach my $intron_chain ( keys %{ $transcript_idx_rf->{introns} } ) {
            warn
"intron:($intron_chain) already added in the database index, it may due to two of transcripts had same intron-chain in the referece\n"
              if exists $transcriptome_rf->{introns}{$intron_chain};
            $transcriptome_rf->{introns}{$intron_chain} = $transcript_rf;
        }
    }
    else {
        foreach my $locus ( keys %{ $transcript_idx_rf->{exons} } ) {
            warn
"exon:()$locus) already added in the database index, t may due to two of single exon transcripts had overlapped with each other in the referece\n"
              if exists $transcriptome_rf->{exons}{$locus};
            $transcriptome_rf->{exons}{$locus} = $transcript_rf;
        }
    }
}

#===  FUNCTION  ================================================================
#         NAME: construct_transcriptome_idx{
#      PURPOSE: give a transcriptome hash return a hashref of transcriptome
#               as: $transcriptome_data{$intron_chain} = $transcript_rf
#   PARAMETERS: $transcirpt  hash of transcript
#      RETURNS: $intron_chain, string
#  DESCRIPTION:
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub construct_transcriptome_idx {
    my $transcriptome_rf = shift;
    my %transcriptome_database;
    foreach my $tid ( keys %$transcriptome_rf ) {
        my $transcript_rf = $transcriptome_rf->{$tid};
        add_idx_to_transcriptome( $transcript_rf, \%transcriptome_database );
    }
    return \%transcriptome_database;
}

#===  FUNCTION  ================================================================
#         NAME: exists_in_db
#      PURPOSE: check if transcript exists in transcript database by mismatch
#   PARAMETERS: $transcript_rf: a hash reference of transcript_b
#   			$database_idx_rf: a hash reference of transcript_a
#				$mismatch:	      allowed mismatch
#      RETURNS: 1 or 0
#  DESCRIPTION:
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub exists_in_db {
    my $transcript_rf   = shift;
    my $database_idx_rf = shift;

    my $transcript_idx_rf = construct_transcript_idx($transcript_rf);

    my @tids = ();

    foreach my $intron_chain ( keys %{ $transcript_idx_rf->{introns} } ) {
        if ( exists $database_idx_rf->{introns}{$intron_chain} ) {
            my $tid = $database_idx_rf->{introns}{$intron_chain}{transcript_id};
            push @tids, $tid;
        }
        else {

        }
    }

    foreach my $locus ( keys %{ $transcript_idx_rf->{exons} } ) {
        if ( exists $database_idx_rf->{exons}{$locus} ) {
            my $tid = $database_idx_rf->{exons}{$locus}{transcript_id};
            push @tids, $tid;
        }
        else {

        }
    }

    return @tids;
}

sub reformat_tid {
    my $transcript_id = shift;
    my @tmp = split /\./, $transcript_id;
    if ( @tmp == 4 ) {
        $tmp[2] = "y" if $tmp[2] eq '=';
        $transcript_id = join ".", @tmp;
    }
    elsif ( @tmp == 5 ) {
        $transcript_id = join ".", ( $tmp[0], $tmp[1], "d", $tmp[-1] );
    }
    else {
        warn "$transcript_id is not well formated\n";
    }
    return $transcript_id;
}

