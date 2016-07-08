package ReseqTrack::Tools::RunVariantCall::CallBySamtools;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists check_executable get_lines_from_file);
use POSIX qw(ceil);

use base qw(ReseqTrack::Tools::RunVariantCall);

=head2 new

  Arg [-bcftools_path]    :
          string, optional, path to bcftools if it cannot be worked out from samtools path      
  Arg [-parameters]   :
      hashref, command line options to use with "samtools mpileup", "bcftools call"
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
                -reference               => '/path/to/ref/',
                -parameters              => {	'mpileup'=>'-ug', 
                								'bcfcall'=>'-mv', 
                							 }
                -chrom                   => '1',
                -region                  => '1-1000000',
                -output_name_prefix      => PHASE1,
                -super_pop_name          => EUR (or all)
                );

=cut

sub DEFAULT_OPTIONS { return {
        mpileup => '-Eug -t DP -t SP -t AD -P ILLUMINA -pm3 -F0.2 -C50', 
        bcfcall => '-mvA',
        depth => 250,
        ploidy => 'GRCh38',
    };
}

## the above mpileup parameters are from 1KG phase3 paper, supplementary section, use -d 500 for exome data, the default -d 250 for LC data

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my ( $bcftools, $bgzip)
    = rearrange([ qw(
    BCFTOOLS BGZIP
    )], @args);    
  
  ## Set defaults
  $self->program('samtools') if (! $self->program);
  $self->bcftools($bcftools|| 'bcftools');
  $self->bgzip($bgzip|| 'bgzip');

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

    throw("do not have have a reference") if !$self->reference;
    check_file_exists($self->reference);
    check_executable($self->bcftools);
    check_executable($self->bgzip);

    throw("Please provide parameters for running mpileup") if !$self->options('mpileup');
    throw("Please provide parameters for running bcfcall") if !$self->options('bcfcall');

    my $output_vcf = $self->working_dir .'/'. $self->job_name . '.vcf.gz';
    $output_vcf =~ s{//}{/};
    $self->output_files($output_vcf);
    
    my @cmd_words;
    push(@cmd_words, $self->program, 'mpileup');
    push(@cmd_words, $self->options->{'mpileup'});

    if ($self->options->{'depth'}) {
      push(@cmd_words, '-d', $self->options->{'depth'}); # if depth is not provided, then the default of 250 will be used for -d
    }
	
    push(@cmd_words, '-f', $self->reference);

    if (my $region = $self->chrom) {
      if (defined $self->region_start && defined $self->region_end) {
        $region.= ':' . $self->region_start . '-' . $self->region_end;
      }
      push(@cmd_words, '-r', $region);
    }    
    push(@cmd_words, @{$self->input_files});

    push(@cmd_words, '|', $self->bcftools, 'call');
    push(@cmd_words, $self->options->{'bcfcall'});
    
    my $subset_ped;
    if ($self->options->{'sample_ped'}) {
        $subset_ped = $self->subset_ped_by_input_bam;
        #print "subset ped is $subset_ped\n";
    }    
    if ($self->options->{'ploidy_file'}) {
	    push(@cmd_words, "--ploidy-file", $self->options->{'ploidy_file'});
	    push(@cmd_words, "-S" ,  $subset_ped);
    }
    elsif  ($self->options->{'ploidy'})   {
        push(@cmd_words, "--ploidy", $self->options->{'ploidy'});
        push(@cmd_words, "-S" ,  $subset_ped);
    }

    push(@cmd_words, '|', $self->bgzip, '-c');
    push(@cmd_words, '>', $output_vcf);

    my $cmd = join(' ', @cmd_words);

    $self->execute_command_line ($cmd);

	unlink $subset_ped;
    return $self;
}

## as each input transposed bams will contain a subset of all samples, and each sample list are likely to be different between chrom regions
## 'bcftools call' -S demands a ped file that matches the samples in bam files  
## this sub is to create a ped file containing exactly the same samples as in the input bam
sub subset_ped_by_input_bam {
	my ($self) = @_;
	
	my $main_ped_lines = get_lines_from_file($self->options->{'sample_ped'});
	
	my $tmp_sample_list = $self->output_files->[0] . ".samples";
	$tmp_sample_list =~ s/.vcf.gz//;
	
	my $subset_ped_file = $self->output_files->[0] . ".ped";
	$subset_ped_file =~ s/.vcf.gz//;
	
	my $del_cmd = "rm $tmp_sample_list";
	$self->execute_command_line ($del_cmd) if (-e $tmp_sample_list);
	
	foreach my $bam ( @{$self->input_files})  {
		my $cmd = $self->program . " view " . $bam . " -H " . "| grep \@RG >> $tmp_sample_list";
		$self->execute_command_line ($cmd); 
	}
	
	my $samples = get_lines_from_file($tmp_sample_list);

	my %sample_hash;
	foreach my $line (@$samples) {
        chomp $line;
        my @tmps = split(/\t/, $line);
        my $sample;
        foreach my $tmp (@tmps) {
            $sample = $tmp if $tmp =~ /^SM:/;
        }    
        $sample =~ s/^SM://;
        $sample_hash{$sample} = 1;
	}
	
	open(OUT, ">", $subset_ped_file) || throw("Cannot open output file $subset_ped_file");
	foreach my $main_ped_line (@$main_ped_lines) {
        chomp $main_ped_line;
        next if ($main_ped_line =~ /Family/);       
        my @tmp = split(/\t/, $main_ped_line);
        $tmp[2] = 0;
        $tmp[3] = 0;
        if ($sample_hash{$tmp[1]}) {
        	print OUT join("\t", @tmp) . "\n";
        }
        
	}
	unlink $tmp_sample_list;
	#$self->{'subset_ped_file'} = $subset_ped_file;    
	return  $subset_ped_file;  
}
	

sub bcftools{
  my ($self, $bcftools) = @_;
  if ($bcftools) {
    $self->{'bcftools'} = $bcftools;
  }
  return $self->{'bcftools'};
}

sub bgzip{
  my ($self, $bgzip) = @_;
  if ($bgzip) {
    $self->{'bgzip'} = $bgzip;
  }
  return $self->{'bgzip'};
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
>bcftools call var.raw.bcf > var.flt.vcf

Petr's parameter when doing the 1kg runs:

samtools mpileup -EDS -C50 -m2 -F0.0005 -d 10000 -P ILLUMINA
bcftools view -p 0.99

and with ploidy 1 for chrY and non-pseudosomal regions of male chrX
(1-60000 and 2699521-154931043). Ploidy can be set as an additional
column to -s option of bcftools view, the accepted values are 1 or 2.

For exome studies, extract on-target sites at the end using a tabix command (see mpileup web site use case example)



