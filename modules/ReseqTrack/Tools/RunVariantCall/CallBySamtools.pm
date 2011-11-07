package ReseqTrack::Tools::RunVariantCall::CallBySamtools;
#package RunVariantCall::CallBySamtools;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename qw(basename fileparse);

use base qw(ReseqTrack::Tools::RunVariantCall);

#use lib '/homes/zheng/reseq-personal/zheng/lib';
#use base qw(RunVariantCall);
use ReseqTrack::Tools::RunVariantCallUtils qw(derive_output_file_name);

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
  Arg [-replace_files]   :
      boolean, default 1, any input / intermediate file will be deleted after output is created.
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
                -replace_files 			=> 1,
                -output_to_working_dir 	=> 1 );
=cut

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my ( 	$bcftools_path, 
  		$vcfutils_path, 
  		$mpileup, 
  		$bcfview, 
  		$vcfutils, 
        $replace_files)
    = rearrange( [
         qw( 	BCFTOOLS_PATH 
         		VCFUTILS_PATH 
         		MPILEUP 
         		BCFVIEW 
         		VCFUTILS
                REPLACE_FILES)
		], @args);	
  
  ## Set defaults
  $self->program("/nfs/1000g-work/G1K/work/bin/samtools/samtools") if (! $self->program);
  
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
  
  #$self->output_to_working_dir($output_to_working_dir);

	
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
	
    $self->execute_command_line($cmd);

    return $output_raw_bcf;
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
    
    $self->execute_command_line($cmd);

    return $output_filtered_vcf;
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
		
	my $out1 = derive_output_file_name($dir, $input_bams, $region, "samtools");
	my $out2 = $out1;
	$out2 =~ s/raw/flt/;
	$out2 =~ s/bcf/vcf/;
	
	my $intermediate_bcf = $self->run_variant_calling($input_bams, $out1);

	my $flt_vcf = $self->run_variant_filtering($intermediate_bcf, $out2);

    $self->output_files($intermediate_bcf );
    $self->output_files($flt_vcf);

    return;

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


=head2 replace_files   ### FIXME: can delete this, it is not much use

  Arg [1]   : ReseqTrack::Tools::RunSamtools
  Arg [2]   : boolean, optional, value of replace_files flag
  Function  : accessor method for replace_files flag
  Returntype: boolean, replace_files flag
  Exceptions: n/a
  Example   : $self->replace_files(1);

=cut

sub replace_files {
    my $self = shift;
    
    if (@_) {
        $self ->{'replace_files'} = (shift) ? 1 : 0;
    }
    return $self->{'replace_files'};
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

=cut



