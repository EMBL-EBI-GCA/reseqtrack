package ReseqTrack::Tools::RunVariantCall::CallByBcftools;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists check_executable get_lines_from_file);
use POSIX qw(ceil);

use base qw(ReseqTrack::Tools::RunVariantCall);

=head2 new

  Arg [-parameters]   :
      hashref, command line options to use with "bcftools mpileup", "bcftools call"
      Here is a list of options: bcftools mpileup [-E] [-P platform STR] [-p] [-m min-ireads INT] [-F gap-frac FLOAT] [-a annotate LIST] [ -C adjust-MQ INT] 
      Here is a list of options: bcftools call [-m] [-v] [-A] [--ploidy <assembly>] [-o output <file>] [-O output-type] [-S samples-file <file>]

  Arg [-super_pop_name]    :
      string, default is "unknownPop", used in output file name and collection name
      
  + Arguments for ReseqTrack::Tools::RunProgram parent class

  Function  : Creates a new ReseqTrack::Tools::RunVariantCall::CallByBcftools object.
  Returntype: ReseqTrack::Tools::RunVariantCall::CallByBcftools
  Exceptions: 
  Example   : my $varCall_byBcftools = ReseqTrack::Tools::RunVariantCall::CallByBcftools->new(
                -input_files             => ['/path/sam1', '/path/sam2'],
                -program                 => "/path/to/bcftools",
                -samtools                 => '/path/to/samtools/',
                -working_dir             => '/path/to/dir/',  ## this is working dir and output dir
                -reference               => '/path/to/ref/',
                -options              => {	'mpileup'=>'-ug', 
                                                'bcfcall'=>'-mv', 
                			    }
                -chrom                   => '1',
                -region                  => '1-1000000',
                -output_name_prefix      => PHASE1,
                -super_pop_name          => EUR (or all)
                );

=cut

sub DEFAULT_OPTIONS { return {
        mpileup => '-E -a DP -a SP -a AD -P ILLUMINA -pm3 -F0.2 -C50', 
        bcfcall => '-mv -O z',
        depth => 700000,
        ploidy => 'GRCh38'
    };
}

## the above mpileup parameters are from 1KG phase3 paper, supplementary section, use -d "500 X number of samples" for exome data, use -d "250 x number of samples"  for LC data

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my ( $bcftools, $samtools)
    = rearrange([ qw(
    BCFTOOLS SAMTOOLS)], @args);    
  
  ## Set defaults
  $self->program('bcftools') if (! $self->program);
  $self->samtools($samtools || 'samtools');

  return $self;
}


=head2 run_program

  Arg [1]   : ReseqTrack::Tools::RunVariantCall::CallByBcftools
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
    check_executable($self->samtools);

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

    push(@cmd_words, '|', $self->program, 'call');
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

    push(@cmd_words, "-o", $output_vcf);

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
		my $cmd = $self->samtools . " view " . $bam . " -H " . "| grep \@RG >> $tmp_sample_list";
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

sub samtools{
    my ($self, $samtools) = @_;
    if ($samtools) {
	$self->{'samtools'} = $samtools;
    }
    return $self->{'samtools'};
}

	
1;

=pod

=head1 NAME

ReseqTrack::Tools::RunVariantCall::CallByBcftools

=head1 SYNOPSIS

This is a class for running BCFtools mpileup and BCFtools call to identify the variants
It is a sub class of a ReseqTrack::Tools::RunVariantCall.

It generates and run command like to call SNPs and indels for multiple diploid individuals:

>bcftools mpileup -P ILLUMINA -E ref.fa *.bam | bcftools call -o new.vcf.gz -O z




