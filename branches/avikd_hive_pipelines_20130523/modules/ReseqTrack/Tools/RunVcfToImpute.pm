#!/usr/bin/env perl

=pod

=head1 NAME

ReseqTrack::Tools::RunVcfToImpute

=head1 SYNOPSIS

This is a class for running vcf2impute_gen.pl
It is a sub class of a ReseqTrack::Tools::RunProgram

It generates and run command to run vep 

>perl /path/vcf2impute_gen.pl -vcf input.vcf -gen /path/output_dir/output.gen -samp_snptest -chr chrom -start region_start -end region_end


=cut


package ReseqTrack::Tools::RunVcfToImpute;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::FileSystemUtils qw (check_executable);
use File::Basename qw(fileparse);

use base qw(ReseqTrack::Tools::RunProgram);

=head2 new

  Arg [-chrom]   :
      string, chromosome to include in output files, in (chr)[1-22,X]
  Arg [-region_start]   :
      string, chromosomal coordinate
  Arg [-region_end]   :
      string, chromosomal coordinate

    + Arguments for ReseqTrack::Tools::RunProgram parent class

  Function  : Creates a new ReseqTrack::Tools::RunVcfToImpute object.
  Returntype: ReseqTrack::Tools::RunVcfToImpute
  Exceptions: 
  Example   : my $run_vcf2impute = ReseqTrack::Tools::RunVcfToImpute->new(
                -input_files => '/path/vcf', 
                -program => ‘/path/vcf2impute_gen.pl’,
                -chrom => 1,                
                -working_dir => '/path/output_dir/',
                -region_start => 100000,
                -region_end => 200000,
                );
            
=cut

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);
  
  my ( $chrom, $region_start, $region_end ) = rearrange( [ qw( CHROM REGION_START REGION_END ) ], @args);
  
  $self->chrom($chrom);
  $self->region_start($region_start);
  $self->region_end($region_end);
          
  return $self;
}

sub run_program {
	my $self = shift;
	
	throw "Can't accept multiple input files" if(scalar @{$self->input_files} >1);   ### ?
	
	my @cmd_words = ("perl ".$self->program);
	
	push(@cmd_words, "-vcf", ${$self->input_files}[0]);
	
	push(@cmd_words, "-chr", $self->chrom) if($self->chrom);
    
    push(@cmd_words, "-start", $self->region_start) if($self->region_start);
    
    push(@cmd_words, "-end", $self->region_end) if($self->region_end);
	
	my $output_dir = $self->output_dir;
	
	my $filename = fileparse(${$self->input_files}[0], qw( .vcf .vcf.gz ));
	
    my $output_gen = "$output_dir/$filename";
    
    $output_gen.="_".$self->chrom if($self->chrom);
    
    $output_gen.="_".$self->region_start if($self->region_start);
    
    $output_gen.="_".$self->region_end if($self->region_end); 
    
    $output_gen.=".gen";  
    
    push(@cmd_words, "-gen", $output_gen );
    
    my $output_samples=$output_gen.".samples",
    $output_gen.=".gz";
    
    $self->output_files($output_gen);
    
    $self->output_files($output_samples);
    
    push(@cmd_words, "-samp_snptest");
    
    
    if ( $self->options ) 
    {
      OPTION:
      while (my ($tag, $value) = each %{$self->options}) 
      {
        next OPTION if !defined $value;
      }
       foreach my $tag ( keys ( %{$self->options} ) ) 
       {
         if(defined $self->options->{$tag})
         {
           my $value = $self->options->{$tag};
           push(@cmd_words,"-".$tag." ".$value);
         }
         else
         {
           push(@cmd_words, "-". $tag);
         }
       }
    }
    
    my $cmd = join(' ', @cmd_words);
    $self->execute_command_line ($cmd);
    
    return $self;
    
        
}

sub chrom {
  my ( $self, $arg ) = @_;

  if ($arg) {
    $self->{chrom} = $arg;
  }
  return $self->{chrom};
}

sub region_start {
  my ( $self, $arg ) = @_;

  if ($arg) {
    $self->{region_start} = $arg;
  }
  return $self->{region_start};
}

sub region_end {
  my ( $self, $arg ) = @_;

  if ($arg) {
    $self->{region_end} = $arg;
  }
  return $self->{region_end};
}

sub output_gen_files {
    my $self = shift;
    my @files = grep { /\.gen.gz$/ } @{ $self->output_files };
    return \@files;
}

sub output_samples_files {
    my $self = shift;
    my @files = grep { /\.gen.samples$/ } @{ $self->output_files };
    return \@files;
}

1;