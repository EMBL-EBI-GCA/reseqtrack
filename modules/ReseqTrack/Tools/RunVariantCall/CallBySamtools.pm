package ReseqTrack::Tools::RunVariantCall::CallBySamtools;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename qw(basename fileparse);

use base qw(ReseqTrack::Tools::RunVariantCall);
use File::Path;

=head2 new

  Arg [-bcftools_path]    :
          string, optional, path to bcftools if it cannot be worked out from samtools path
  Arg [-vcfutils_path]    :
          string, optional, path to vcfutils.pl if it cannot be worked out from samtools path         
  Arg [-parameters]   :
      hashref, command line options to use with "samtools mpileup", "bcftools view" and "vcfutils.pl"
      Here is a list of options: samtools mpileup [-EBug] [-C capQcoef]  [-l list] [-M capMapQ] [-Q minBaseQ] [-q minMapQ] 
      Here is a list of options: bcftools view [-AbFGNQSucgv] [-D seqDict] [-l listLoci] [-i gapSNPratio] [-t mutRate] [-p varThres] [-P prior] [-1 nGroup1] [-d minFrac] [-U nPerm] [-X permThres] [-T trioType] 
  Arg [-super_pop_name]    :
      string, default is "unknownPop", used in output file name and collection name
      
  + Arguments for ReseqTrack::Tools::RunProgram parent class

  Function  : Creates a new ReseqTrack::Tools::RunVariantCall::CallBySamtools object.
  Returntype: ReseqTrack::Tools::RunVariantCall::CallBySamtools
  Exceptions: 
  Example   : my $varCall_bySamtools = ReseqTrack::Tools::RunVariantCall::CallBySamtools->new(
                -input_files             => ['/path/sam1', '/path/sam2'],
                -program                 => "/path/to/samtools",
                -working_dir             => '/path/to/dir/',  ## this is working dir and output dir
                -bcftools_path           => '/path/to/bcftools/',
                -vcfutils_path           => '/path/to/vcfutils.pl/',
                -reference               => '/path/to/ref/',
                -parameters              => {	'mpileup'=>'-ug', 
                								'bcfview'=>'-bvcg', 
                								'vcfutils'=>'-D 100' }
                -chrom                   => '1',
                -region                  => '1-1000000',
                -output_name_prefix      => PHASE1,
                -super_pop_name          => EUR (or all)
                );

=cut

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my (     $bcftools_path, 
          $vcfutils_path, 
		  $parameters,
          $super_pop_name)
    = rearrange( [
         qw(     BCFTOOLS_PATH 
                 VCFUTILS_PATH 
				 PARAMETERS
                 SUPER_POP_NAME )
        ], @args);    
  
  ## Set defaults
  $self->program("/nfs/1000g-work/G1K/work/bin/samtools_latest/samtools") if (! $self->program);
  $self->reference("/nfs/1000g-work/G1K/scratch/zheng/reference_genomes/human_g1k_v37.fasta") if (! $self->reference);
  
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

  $self->options($parameters);
  print "option hash keys are:\n";
  print join("\n", keys %{$self->options}) . "\n";
  
  $self->super_pop_name($super_pop_name);
  
  return $self;
}

=head
sub DEFAULT_OPTIONS { return {
        'mpileup' => '-ug',
        'bcfview' => '-bvcg',
        'vcfutils' => '-d 2',
        };
}
=cut

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
	
	if ($self->options->{'mpileup'} ) {
		$cmd .= " " . $self->options->{'mpileup'};
		print "mpileup option is " . $self->options->{'mpileup'} . "\n";	
	}			
	else {
		throw("Please provide parameters for running mpileup using tag -parameters");
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
	
	if ($self->options->{'bcfview'} ) {
		$cmd .= " " . $self->options->{'bcfview'};
		print "bcfview option is " . $self->options->{'bcfview'} . "\n";	
	}			
	else {
		throw("Please provide parameters for running bcfview using tag -parameters");
	}   

    $cmd .= " - > $output_raw_bcf";

    $self->created_files($output_raw_bcf);
    $self->execute_command_line($cmd);

#    return $self->intermediate_output_file($output_raw_bcf);
	return $self;	
}

=head2 run_variant_filtering

  Arg [1]   : ReseqTrack::Tools::RunVariantCall::CallBySamtools
  Arg [2]   : string, paths of raw bcf file generated by run_variant_calling
  Arg [3]    : string, path to reference genome
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
		
	if ($self->options->{'vcfutils'} ) {
		$cmd .= " " . $self->options->{'vcfutils'};
		print "vcfutils option is " . $self->options->{'vcfutils'} . "\n";	
	}		
	else {
		throw("Please provide parameters for running vcfutils using tag -parameters");
	}
    
    $cmd .= " > $output_filtered_vcf ";
    $self->execute_command_line($cmd);
    

    return $self;
}
    

=head2 run_program

  Arg [1]   : ReseqTrack::Tools::RunVariantCall::CallBySamtools
  Function  : uses samtools to call variants in $self->input_files.
              Output is files are stored in $self->output_files
              This method is called by the RunProgram run method
  Returntype: 
  Exceptions: 
  Example   : $self->run();

=cut

sub run_program {
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


=head2 super_pop_name
  Arg [1]   : ReseqTrack::Tools::RunVariantCall::CallBySamtools
  Arg [2]   : string, required, super_pop_name used for calling the variants, can be things like EUR, ALL (for all pop) and unknownPop, it will be 
            used in output file names
  Function  : accessor method for super_pop_name
  Returntype: string
  Exceptions: n/a
  Example   : my $super_pop_name = $self->super_pop_name;

=cut

sub super_pop_name {
  my ($self, $super_pop_name) = @_;
  if ($super_pop_name) {
    $self->{'super_pop_name'} = $super_pop_name;
  }
  return $self->{'super_pop_name'};
}


=head2 derive_output_file_name 

  Arg [1]   : ReseqTrack::Tools::RunVariantCall::CallBySamtools object
  Function  : create an output VCF name based on input file information; If chromosome is specified, 
  put the output file in a sub-directory called "chrN"; to reduce number of files in the output directory
  Returntype: file path
  Exceptions: 
  Example   : my $output_file = $self->derive_output_file_name->[0];

=cut


sub derive_output_file_name {  
    my ( $self ) = @_;        

    my $sample_cnt = @{$self->input_files};
    my $output_file;

    my $output_dir_by_chr;
    my $out_dir = $self->working_dir;
    $out_dir =~ s/\/$//;
    
    if ( $self->chrom ) {
        $output_dir_by_chr = $out_dir . "/chr" . $self->chrom;
    }
    else{
        $output_dir_by_chr = $out_dir;
    }    
    
    mkpath($output_dir_by_chr) unless (-e $output_dir_by_chr);
    
    if ($self->region) {
        $output_file = $output_dir_by_chr . "/" . $self->output_name_prefix . "_" . $self->super_pop_name . "_of_" . $sample_cnt . "bams.chr" . $self->chrom . "_" . $self->region . ".samtools.flt.vcf";

    }
    else {
        $output_file = $output_dir_by_chr . "/" . $self->output_name_prefix . "_" . $self->super_pop_name . "_of_" . $sample_cnt . "bams.samtools.flt.vcf";
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

For exome studies, extract on-target sites at the end using a tabix command (see mpileup web site use case example)



