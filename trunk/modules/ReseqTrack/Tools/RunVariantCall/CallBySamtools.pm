package ReseqTrack::Tools::RunVariantCall::CallBySamtools;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename qw(basename fileparse);

use base qw(ReseqTrack::Tools::RunVariantCall);

=head2 new

  Arg [-bcftools_path]	:
  		string, optional, path to bcftools if it cannot be worked out from samtools path
  Arg [-vcfutils_path]	:
  		string, optional, path to vcfutils.pl if it cannot be worked out from samtools path 		
  Arg [-mpileup]   :
      string, command line options to use with "samtools mpileup"
      Here is a list of options: samtools mpileup [-EBug] [-C capQcoef]  [-l list] [-M capMapQ] [-Q minBaseQ] [-q minMapQ] 
  Arg [-bcfview]   :
      string, command line options to use with "bcftools view"
      Here is a list of options: bcftools view [-AbFGNQSucgv] [-D seqDict] [-l listLoci] [-i gapSNPratio] [-t mutRate] [-p varThres] [-P prior] [-1 nGroup1] [-d minFrac] [-U nPerm] [-X permThres] [-T trioType] 
  Arg [-vcfutils]   :
      string, command line options to use with "vcfutils.pl"

  + Arguments for ReseqTrack::Tools::RunProgram parent class

  Function  : Creates a new ReseqTrack::Tools::RunVariantCall::CallBySamtools object.
  Returntype: ReseqTrack::Tools::RunVariantCall::CallBySamtools
  Exceptions: 
  Example   : my $varCall_bySamtools = ReseqTrack::Tools::RunVariantCall::CallBySamtools->new(
                -input_files 			=> ['/path/sam1', '/path/sam2'],
                -program 				=> "/path/to/samtools",
                -working_dir 			=> '/path/to/dir/',
                -bcftools_path 			=> '/path/to/bcftools/',
                -vcfutils_path 			=> '/path/to/vcfutils.pl/',
                -reference				=> '/path/to/ref/',
                -mpileup 				=> '-ug',
                -bcfview 				=> '-bvcg',
                -vcfutils 				=> '-D 100',
                -chrom					=> '1',
                -region					=> '1-1000000',
                -output_to_working_dir 	=> 1 );

=cut

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my ( 	$bcftools_path, 
  		$vcfutils_path, 
  		$mpileup, 
  		$bcfview, 
  		$vcfutils)
    = rearrange( [
         qw( 	BCFTOOLS_PATH 
         		VCFUTILS_PATH 
         		MPILEUP 
         		BCFVIEW 
         		VCFUTILS)
		], @args);	
  
  ## Set defaults
  $self->program("/nfs/1000g-work/G1K/work/bin/samtools/samtools") if (! $self->program);
  $self->reference("/nfs/1000g-work/G1K/scratch/zheng/reference_genomes/human_g1k_v37.fasta") if (! $self->reference);
  
  $self->{'options'}->{'mpileup'} = '-ug' unless ($mpileup);
  $self->{'options'}->{'bcfview'} = '-bvcg' unless ($bcfview);
  $self->{'options'}->{'vcfutils'} = '-D 100' unless ($vcfutils);
  
  my $samtools_path = $self->program; 
  my ($tool, $dir) = fileparse($samtools_path);	
  if ( !$bcftools_path ) {
	my $derived_bcftools_path = $dir . "/bcftools/bcftools";
	if ( -e  $derived_bcftools_path ) {
		$self->bcftools_path($derived_bcftools_path);
  	}
  	else {
  		throw("Cannot find bcftools in the samtools installation. Please provide path to the program.");
  	}        
  }	
  else {
  	$self->bcftools_path($bcftools_path);
  }
  
  if ( !$vcfutils_path ) {
  	my $derived_vcfutils_path = $dir . "/bcftools/vcfutils.pl";
  	if (-e $derived_vcfutils_path) {
	  	$self->vcfutils_path($derived_vcfutils_path);
  	}
  	else {
  	    throw("Cannot find bcftools/vcfutils.pl in the samtools installation. Please provide path to the program.");
  	}    
  }
  else {
	  $self->vcfutils_path($vcfutils_path);
  }

  print "input mpileup option is $mpileup\n"; 	
  $self->options('mpileup', $mpileup);
  $self->options('bcfview', $bcfview);
  $self->options('vcfutils', $vcfutils); 
  
  if (! $self->job_name) {
      $self->generate_job_name;
  }

  return $self;
}


=head2 run_variant_calling

  Arg [1]   : ReseqTrack::Tools::RunVariantCall::CallBySamtools
  Arg [2]   : string, paths of bam files
  Function  : uses samtools mpileup and bcftools view to call variants
  Returntype: string, path of raw vcf file before filtering
  Exceptions: 
  Example   : my $raw_vcf = $self->run_variant_calling($bams);

=cut

sub run_variant_calling {
    my ($self, $input_bams, $output_raw_bcf) = @_;  # $input_bams can be an array ref

    my $cmd = $self->program . " mpileup";

    if ($self->options('mpileup')) {
        $cmd .= " " . $self->options('mpileup');
        print "mpileup option is " . $self->options('mpileup') . "\n";
    }           

	$cmd .= " -f " . $self->reference;
	
	foreach my $bam ( @$input_bams ) {  ## FIXME: use -b FILE to take a list of bam files. one bam per line
		$cmd .= " " . $bam;
	}	
	
	if ($self->chrom) {
        $cmd .= " -r " . $self->chrom;
    	if ($self->region) {
    	    $cmd .= ":" . $self->region;
    	}
    }	
	
	$cmd .= " | " . $self->bcftools_path . " view";
	
	if ($self->options('bcfview')) {
        $cmd .= " " . $self->options('bcfview');
        print "bcfview option is " . $self->options('bcfview') . "\n";
    }  

	$cmd .= " - > $output_raw_bcf";
	print "Running command.................................\n";
	my $exit;
	eval{
		$exit = $self->execute_command_line($cmd);
	};
	if($exit > 0){
		throw("Failed to run command\n$cmd\n". @_  . " exit code $exit");
	}

	return $self->intermediate_output_file($output_raw_bcf);
}

=head2 run_variant_filtering

  Arg [1]   : ReseqTrack::Tools::RunVariantCall::CallBySamtools
  Arg [2]   : string, paths of raw bcf file generated by run_variant_calling
  Arg [3]	: string, path to reference genome
  Function  : uses samtools mpileup and bcftools view to call variants
  Returntype: string, path of raw vcf file before filtering
  Exceptions: 
  Example   : my $filtered_vcf = $self->run_variant_filtering($raw_bcf);

=cut

#bcftools view var.raw.bcf | vcfutils.pl varFilter -D 100 > var.flt.vcf

sub run_variant_filtering {
    my ($self, $input_raw_bcf, $output_filtered_vcf) = @_;  


    my $cmd = $self->bcftools_path . " view $input_raw_bcf | ";
	$cmd .= $self->vcfutils_path . " varFilter";

    if ($self->options('vcfutils')) {
        $cmd .= " " . $self->options('vcfutils');
    }   
    
    $cmd .= " > $output_filtered_vcf ";
	
	print "Running command.................................\n";
	my $exit;
	eval{
		$exit = $self->execute_command_line($cmd);
	};
	if( $exit > 0 ){
		throw("Failed to run command\n$cmd\n". @_  . " exit code $exit");
	} 

	return $self;
}
    

=head2 run

  Arg [1]   : ReseqTrack::Tools::RunVariantCall::CallBySamtools
  Function  : uses samtools to call variants in $self->input_files.
  Output is files are stored in $self->output_files
  Returntype: 
  Exceptions: 
  Example   : $self->run();

=cut

sub run {
    my ($self) = @_;
	
	my $dir = $self->working_dir;
	my $input_bams = $self->input_files;
	my $chr = $self->chrom;
	my $r = $self->region;
	
	my $region;
	if ($r) {
		$region = "chr" . $chr . "_" . $r; 
	}
	elsif ($chr) {
		$region = "chr" . $chr;
	}	
		
	my $flt_out = $self->derive_output_file_name->[0];
	
	my $raw_out = $flt_out;
	$raw_out =~ s/flt/raw/;
	$raw_out =~ s/vcf/bcf/;
	
	$self->run_variant_calling($input_bams, $raw_out);
	$self->run_variant_filtering($raw_out, $flt_out);
	
	$self->files_to_delete($self->intermediate_output_file);

    return $self;

}

=head2 bcftools_path

  Arg [1]   : ReseqTrack::Tools::RunVariantCall::CallBySamtools
  Arg [2]   : string, optional, path of bcftools program
  Function  : accessor method for bcftools
  Returntype: string
  Exceptions: n/a
  Example   : my $bcftools_path = $self->bcftools_path;

=cut


sub bcftools_path {
  my ($self, $bcftools_path) = @_;
  if ($bcftools_path) {
    $self->{'bcftools_path'} = $bcftools_path;
  }
  return $self->{'bcftools_path'};
}

=head2 vcfutils_path

  Arg [1]   : ReseqTrack::Tools::RunVariantCall::CallBySamtools
  Arg [2]   : string, optional, path of vcfutils program
  Function  : accessor method for vcfutils
  Returntype: string
  Exceptions: n/a
  Example   : my $vcfutils_path = $self->vcfutils_path;

=cut


sub vcfutils_path {
  my ($self, $vcfutils_path) = @_;
  if ($vcfutils_path) {
    $self->{'vcfutils_path'} = $vcfutils_path;
  }
  return $self->{'vcfutils_path'};
}


=head2 derive_output_file_name 

  Arg [1]   : ReseqTrack::Tools::RunVariantCall::CallBySamtools object
  Function  : create an output VCF name based on input file information
  Returntype: file path
  Exceptions: 
  Example   : my $output_file = $self->derive_output_file_name->[0];

=cut


sub derive_output_file_name {  
	my ( $self ) = @_;		
	my $first_bam = basename($self->input_files->[0]);
	my @tmp = split(/\./, $first_bam);
	my $first_sample = $tmp[0];
	my $sample_cnt = @{$self->input_files} - 1;
	
	my $output_file;


	if ($self->region) {
		$output_file = $self->working_dir . "/$first_sample" . "_and_" . $sample_cnt . "_others.chr" . $self->chrom . "_" . $self->region . ".samtools.flt.vcf";
	}
	else {
		$output_file = $self->working_dir . "/$first_sample" . "_and_" . $sample_cnt . "_others.samtools.flt.vcf";
	}

	return $self->output_files($output_file);
}	

1;

=pod

=head1 NAME

ReseqTrack::Tools::RunVariantCall::CallBySamtools

=head1 SYNOPSIS

This is a class for running mpileup and BCFtools to call variants
It is a sub class of a ReseqTrack::Tools::RunVariantCall.

It generates and run command like to call SNPs and indels for multiple diploid individuals:

>samtools mpileup -P ILLUMINA -ugf ref.fa *.bam | bcftools view -bcvg - > var.raw.bcf 
>bcftools view var.raw.bcf | vcfutils.pl varFilter -D 2000 > var.flt.vcf

Petr's parameter when doing the 1kg runs:

samtools mpileup -EDS -C50 -m2 -F0.0005 -d 10000 -P ILLUMINA
bcftools view -p 0.99

and with ploidy 1 for chrY and non-pseudosomal regions of male chrX
(1-60000 and 2699521-154931043). Ploidy can be set as an additional
column to -s option of bcftools view, the accepted values are 1 or 2.

=cut



