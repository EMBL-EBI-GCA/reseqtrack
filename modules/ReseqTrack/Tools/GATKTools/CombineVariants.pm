=pod

=head1 NAME

ReseqTrack::Tools::GATKTools::CombineVariants

=head1 SYNOPSIS

Object to create a combined vcf file  using GATK CombineVariants

example

my $COMBINE_VARIANTS = ReseqTrack::Tools::GATKTools::CombineVariants->new(
     -variants_file     => {gatk=>'/path/to/gatk.vcf',samtools=>'/path/to/samtools.vcf',freebayes=>'/path/to/freebayes.vcf',},
     -java_exe        =>"/usr/bin/java" ,
     -jvm_args        =>"-Xmx4g",
     -options         => {'filteredRecordsMergeType' => ‘KEEP_UNCONDITIONAL‘,  }
     -job_name => 'snps',
     -GATK_PATH       =>$GATK/GenomeAnalysisTK/",
     -working_dir     => '/path/to/dir/',
     -reference       => '/path/to/reference',
     -minimal_vcf     => 1,
);


=cut

package ReseqTrack::Tools::GATKTools::CombineVariants;

use strict;
use warnings;

use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists check_executable);

use base qw(ReseqTrack::Tools::GATKTools);

sub DEFAULT_OPTIONS { return {

        
        };
}

sub new {
	my ( $class, @args ) = @_;
	my $self = $class->SUPER::new(@args);
	
    my ( $variants, $minimal_vcf, $suppress_commandline_header, $filter_mbq, $assume_identical_samples, $filtered_are_uncalled)
        = rearrange( [ qw(  VARIANTS 
                            MINIMAL_VCF 
                            SUPPRESS_COMMANDLINE_HEADER 
                            FILTER_MBQ 
                            ASSUME_IDENTICAL_SAMPLES 
                            FILTERED_ARE_UNCALLED
                         )], @args);

        $self->variants($variants);
        $self->minimal_vcf($minimal_vcf);
        $self->suppress_commandline_header($suppress_commandline_header);
        $self->filter_mbq($filter_mbq);
        $self->assume_identical_samples($assume_identical_samples);
        $self->filtered_are_uncalled($filtered_are_uncalled);
        
    return $self;
}

sub run_program {
	my $self = shift;

        throw "expecting variant list" if(!$self->variants);
               
        $self->check_jar_file_exists;
        check_file_exists($self->reference);
        check_executable($self->java_exe);

        my $combined_vcf = $self->working_dir . '/'
        . $self->job_name. '.combined.vcf.gz';
        
        $combined_vcf =~ s{//}{/}g;
        
        my @cmd_words = ($self->java_exe, $self->jvm_args, '-jar');
        push(@cmd_words, $self->gatk_path . '/' . $self->jar_file);
        push(@cmd_words, '-T', 'CombineVariants');
        
        throw "expecting hash ref" unless (ref($self->variants) eq 'HASH');
       	throw "More than one variants files are required to combine" unless( (keys %{$self->variants}) > 1 );
 
      if ( $self->variants ) {
      VARIANT:
      while (my ($tag, $value) = each %{$self->variants}) {
        next VARIANT if !defined $value;
      }
        foreach my $tag ( keys ( %{$self->variants} ) ) {
            my $tag_vcf = $self->variants->{$tag};
            check_file_exists($tag_vcf);
            push(@cmd_words, "--variant:$tag", $tag_vcf);
        }
    }
    
      if ( $self->options ) {
      OPTION:
      while (my ($tag, $value) = each %{$self->options}) {
        next OPTION if !defined $value;
              }
        foreach my $tag ( keys ( %{$self->options} ) ) {
            my $value = $self->options->{$tag};
                push(@cmd_words, "-$tag", $value);
         }
    }
        push(@cmd_words,'-minimalVCF') if($self->minimal_vcf);
        
        push(@cmd_words,'-suppressCommandLineHeader') if($self->suppress_commandline_header);
        
        push(@cmd_words,'-filterMBQ') if($self->filter_mbq);
        
        push(@cmd_words,'-assumeIdenticalSamples') if($self->assume_identical_samples);
        
        push(@cmd_words,'-filteredAreUncalled') if($self->filtered_are_uncalled);
        
    
        push(@cmd_words, '-R', $self->reference);
        push(@cmd_words, '-o', $combined_vcf);

        my $cmd = join(' ', @cmd_words);

        $self->output_files($combined_vcf);
        $self->execute_command_line ($cmd);
	    return;
}

sub variants {

  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{variants} = $arg;
  }
  return $self->{variants};
}

sub minimal_vcf {

  my $self = shift;
    if (@_) {
        $self->{'minimal_vcf'} = (shift) ? 1 : 0;
    }
    return $self->{'minimal_vcf'};
}

sub suppress_commandline_header {

  my $self = shift;
    if (@_) {
        $self->{'suppress_commandline_header'} = (shift) ? 1 : 0;
    }
    return $self->{'suppress_commandline_header'};
}

sub filter_mbq {

  my $self = shift;
    if (@_) {
        $self->{'filter_mbq'} = (shift) ? 1 : 0;
    }
    return $self->{'filter_mbq'};
}

sub assume_identical_samples {

  my $self = shift;
    if (@_) {
        $self->{'assume_identical_samples'} = (shift) ? 1 : 0;
    }
    return $self->{'assume_identical_samples'};
}

sub filtered_are_uncalled {

  my $self = shift;
    if (@_) {
        $self->{'filtered_are_uncalled'} = (shift) ? 1 : 0;
    }
    return $self->{'filtered_are_uncalled'};
}

1;
