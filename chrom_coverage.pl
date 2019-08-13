#!/usr/bin/perl -w

use strict;
use 5.010;
use Getopt::Long;

# use Cwd 'abs_path';
use Statistics::R;
use IPC::Cmd qw[can_run];
use Bio::DB::Sam;

#usage information
my $usage = <<USAGE;
chrom_coverage V1.2, written by corephi and XiaoMa
This program is used to plot the reads coverage of each chromosome
More scripts? Please join "bioinformatics*CN" QQ group: 276151571
If you have any probel or suggestions about this program, please
mail to: hpyu\@genetics.ac.cn
----------------------------------------------------------
Usage: chrom_coverage [options...] bamfile 

 Options:
  -r|--ref-ids       the reference sequence ids calculated for plot, often 
                     they are chromosome or scanfold names, sepetated are by
                     comma.For example,'-r Chr1,Chr2,Chr3,Chr4,Chr5', this
                     program only plot Chr1-5. default, all the reference
                     sequence ids in bamfile.
  -w|--window-size   the size of window or bin set for calculating read
                     coverage, unit kb, if you haven't specified one, it 
                     will setted to:
                     windows-size = lengh(longetst chromsome) / 2000.
  -o|--out-prefix    Prefix of output pdf file, optional
  --pdf              output the figure in pdf format, default on 10 6, 
                     expample: '--pdf 10 6', output a pdf with width 10 and 
                     height 6
  --png              output the figure in png format, default off.
                     '--png 1366 768', output a png with width 1366 and
                     height 768
  -p|--plot-file     specify the coverage data file calculated by this 
                     program. when this is on, this program will direct use
                     it to plot. 
  -f|--filter-data   discard the outlier define by x * Hinge Spread(HS), 
                     0 means no filter, optional, default 0
  -y|--yaxis-scale   two float number D U, Lower limit is D, Upper limit is
                     U,if not specified, Y axes limit is Upper outlier,only
                     works when --filter-data is set to 0,	
  -b|--yaxis-break   the number of breaks on y-axis, default 3. if y aixs is
                     limited to 0-4, -b is setted as 4, the number on y-axis
                     will be 0,1,2,3,4
  -v|--verbose
  
  Note:
  The default parameters are optimized for arabidopsis.
  -pdf and -png can be setted at the same time, and output two format figures.
  In practice, you can use -f parameter to find a most appropriate y-axis limitation,
  and then use '-y -p -b' three papameters to polish your figure
  
  It may cost a lot of time depends on your genome size, please be patient 

USAGE

#parse parameter
my $windowsize  = 0;
my $prefix      = '';
my $data_file   = '';
my $filter_data = 0;
my $ref_seq_ids = '';
my $y_break     = 3;
my @png         = ();
my @pdf         = ();
my @yaxis_scale = ();
my $verbose     = '';
die $usage
  unless GetOptions(
    "w|window-size:f"    => \$windowsize,
    "r|ref-ids:s"        => \$ref_seq_ids,
    "o|out-prefix:s"     => \$prefix,
    "pdf:f{2}"           => \@pdf,
    "png:f{2}"           => \@png,
    "p|plot-file:s"      => \$data_file,
    "f|filter-data:f"    => \$filter_data,
    "y|yaxis-scale:f{2}" => \@yaxis_scale,
    "b|yaxis-break:i"    => \$y_break,
    "v|verbose!"         => \$verbose,
  );

# $data_file = abs_path($data_file);
my $plot_file_name = $prefix ? $prefix . "_coverage" : "coverage";

#check run environment
my $r_path = can_run('Rscript')
  or die "R is not installed!\nPlease ensure command 'Rscript' is in your PATH";

@pdf = ( 10, 6 ) if ( @pdf == 0 && @png == 0 );

my $bam_file = shift or die $usage;
die "Cannot find $bam_file, specify the correct one.\n"
  unless -e $bam_file;

#calculate chromose coverage
if ($data_file) {
    warn "You specified plot file,skip calucalate coverage...\n" if $verbose;
    die "Cannot find $data_file, specify the correct one.\n"
      unless -e $data_file;
}
else {
    $data_file = $prefix ? $prefix . "_coverage.dat" : "coverage.dat"
      unless $data_file;

    my $windowsize_coverage;
    if ($windowsize) {
        warn "you have specifed window-size:$windowsize kb\n" if $verbose;
        $windowsize_coverage = 1000 * $windowsize;
    }
    else {
        warn "you havene't specifed window-size, calculating one...\n"
          if $verbose;
        my $sam = Bio::DB::Sam->new(
            -bam       => $bam_file,
            -autoindex => 1
        );
        my $max_len = 0;
        foreach my $chr ( $sam->seq_ids ) {
            my $len = $sam->length($chr);
            $max_len = $len > $max_len ? $len : $max_len;
        }
        $windowsize_coverage = int( $max_len / 2000 );
        $windowsize          = $windowsize_coverage / 1000;
        warn "window-size is setted to $windowsize kb\n" if $verbose;
    }

    calc_coverage( $data_file, $bam_file, $windowsize_coverage );
}

##construct R code file
warn "Constructing R code...\n" if $verbose;

#construct header;
my $header = header();
if ($ref_seq_ids) {
	my @chr_temp = split /,/, $ref_seq_ids;
	@chr_temp = map {"'$_'"} @chr_temp;
	my $chr_string = join ', ', @chr_temp;
	$header .= <<CHROM;
	xc\$chrom <- ordered(xc\$chrom, levels = c($chr_string))
CHROM
}

#construct body by parameter
my $body   = '';
my $footer = '';
if ($filter_data) {
    $plot_file_name .= "_with_filter";
    $body = filter_body($filter_data);

    if (@pdf) {
        my ( $width, $height ) = @pdf;
        $footer .= auto_footer( $plot_file_name, "pdf", $width, $height );
    }

    if (@png) {
        my ( $width, $height ) = @png;
        $footer .= auto_footer( $plot_file_name, "png", $width, $height );
    }

}
elsif (@yaxis_scale) {
    $plot_file_name .= "_" . join "_", @yaxis_scale;
    $body = yaxis_body(@yaxis_scale);
    if (@pdf) {
        my ( $width, $height ) = @pdf;
        $footer .= scale_footer( $plot_file_name, "pdf", $width, $height );
    }

    if (@png) {
        my ( $width, $height ) = @png;
        $footer .= scale_footer( $plot_file_name, "png", $width, $height );
    }

}
else {
    $plot_file_name .= "_auto_yaxis";
    $body = auto_body();
    if (@pdf) {
        my ( $width, $height ) = @pdf;
        $footer .= auto_footer( $plot_file_name, "pdf", $width, $height );
    }

    if (@png) {
        my ( $width, $height ) = @png;
        $footer .= auto_footer( $plot_file_name, "png", $width, $height );
    }
}

##Plotting
my $R = Statistics::R->new();
warn "Start Ploting...\n" if $verbose;

my $R_code = $header . $body . $footer;
my $R_out  = $R->run($R_code);
warn $R_out, "\n" if $verbose;
$R->stop;

sub calc_coverage {

    my ( $data_file, $bam_file, $windowsize_coverage ) = @_;

    open my $data_fh, ">", $data_file
      or die "Cannot create file:$data_file\n";
    print $data_fh "chrom\twindow_num\tcoverage\n";

    my $sam = Bio::DB::Sam->new(
        -bam       => $bam_file,
        -autoindex => 1
    );
    my @chrs;
    my %seq_ids = map { $_ => 1 } $sam->seq_ids;
    if ($ref_seq_ids) {
        warn "You have specified reference ids, checking first...\n"
          if $verbose;
        my @chr_temp = split /,/, $ref_seq_ids;
        foreach my $chr (@chr_temp) {
            if ( exists $seq_ids{$chr} ) {
                push @chrs, $chr;
            }
            else {
                warn "$chr you specified is not contained in bamfile\n"
                  . "Please use the referece sequence id same to your gtf/gff\n";
            }
        }
        die "no referece sequence id to calculate, please check\n" unless @chrs;
        warn "After checking with bam header @chrs are left.\n" if $verbose;
    }
    else {
        warn
"You haven't specifed any referece sequcence id, automally calculating...\n"
          if $verbose;
        @chrs = $sam->seq_ids;
        warn "referece sequcence ids are:@chrs\n" if $verbose;
    }
    warn "Start Calulate coverage...\n" if $verbose;
    foreach my $chr (@chrs) {
        warn "Proceeding $chr...\n" if $verbose;
        my $length     = $sam->length($chr);
        my $bin        = int $length / $windowsize_coverage;
        my ($coverage) = $sam->features(
            -type   => "coverage:$bin",
            -seq_id => $chr,
        );
        my @data  = $coverage->coverage;
        my $index = 1;
        for my $density (@data) {
            print $data_fh "$chr\t$index\t$density\n";
            $index++;
        }
    }
    close $data_fh;
    warn "Coverage Data has been stored in $data_file\n" if $verbose;

}

sub header {
    return <<EOF;
require(ggplot2)||{install.packages("ggplot2");require(ggplot2)}
require(plyr)||{install.packages("plyr");require(plyr)}
xc<-read.table(file="$data_file",header=T,stringsAsFactors=F)
xc\$chrom <- as.factor(xc\$chrom)
EOF
}

sub filter_body {
    my $filter_data = shift;
    return <<EOF;
xc_l<-split(xc,xc[,1])
pro_xc<-function(x){
	HS<-diff(range(summary(x[,3])[2],summary(x[,3])[5]))
	U<-as.numeric(summary(x[,3])[5]+HS*$filter_data)
	D<-as.numeric(summary(x[,3])[2]-HS*$filter_data)
	x<-x[x[,3]<U,]
	x<-x[x[,3]>D,]
}
xc_nl<-lapply(xc_l,pro_xc) 
xc_nd<-ldply(xc_nl,data.frame)
xc_nd<-xc_nd[,-1]
xc_nd[,3]<-log2(xc_nd[,3])
xc_nd[xc_nd[,3]<=0,3]<-0.0001
xc <- xc_nd
plot_xc<-function(xc){
	names(xc)<-c("chrom","window_num","coverage")
	p<-ggplot(xc,aes(x=window_num,y=coverage,fill=chrom))+
	geom_bar(stat="identity")+facet_grid(chrom ~.)+
	scale_fill_manual(values=rep("blue", length(xc\$chrom)))+
	xlab("Chromosome position ($windowsize kb)")+
	ylab("Reads coverage (log2)")+
	theme_bw()+theme(legend.position="none")
	print(p)
}
EOF
}

sub auto_body {
    return <<EOF;
xc[,3]<-log2(xc[,3])
xc[xc[,3]<=0,3]<-0.0001 
plot_xc<-function(xc){
	names(xc)<-c("chrom","window_num","coverage")
	HS<-diff(range(summary(xc[,3])[2],summary(xc[,3])[5]))
	U<-as.numeric(summary(xc[,3])[5]+HS*1.5)
	p<-ggplot(xc,aes(x=window_num,y=coverage,fill=chrom))+
	geom_bar(stat="identity")+facet_grid(chrom ~.)+
	scale_fill_manual(values=rep("blue", length(xc\$chrom)))+
	xlab("Chromosome position ($windowsize kb)")+
	ylab("Reads coverage (log2)")+
	scale_y_continuous(breaks=round(seq(0.0001,ceiling(U),U/$y_break)))+
	theme_bw()+theme(legend.position="none")
	print(p)
}
EOF
}

sub yaxis_body {
    my ( $d, $u ) = @_;
    return <<EOF;
D<-as.numeric($d)
U<-as.numeric($u)
xc[,3]<-log2(xc[,3])
xc[xc[,3]<=0,3]<-0.0001
plot_xc<-function(xc,D,U){
	names(xc)<-c("chrom","window_num","coverage")
	p<-ggplot(xc,aes(x=window_num,y=coverage,fill=chrom))+
	geom_bar(stat="identity")+facet_grid(chrom ~.)+
	scale_fill_manual(values=rep("blue", length(xc\$chrom)))+
	xlab("Chromosome position ($windowsize kb)")+
	ylab("Reads coverage (log2)")+
	coord_cartesian(ylim = c(D,U)) +
	scale_y_continuous(breaks=round(seq(floor(D),ceiling(U),U/$y_break))) + 
	theme_bw()+theme(legend.position="none")
	print(p)
}	
EOF
}

sub scale_footer {
    my ( $plot_file_name, $dev, $width, $height ) = @_;
    return <<EOF;
$dev("$plot_file_name.$dev",width=$width,height=$height)
plot_xc(xc,D,U)
dev.off()
EOF
}

sub auto_footer {
    my ( $plot_file_name, $dev, $width, $height ) = @_;
    return <<EOF;
$dev("$plot_file_name.$dev",width=$width,height=$height)
plot_xc(xc)
dev.off()
EOF
}

