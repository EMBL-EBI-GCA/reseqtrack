#!/usr/bin/env perl

=pod

=head1 NAME

ReseqTrack::Tools::RunVep

=head1 SYNOPSIS

This is a class for running Variant Effect Predictor (VEP)

It generates and run command to run vep 

>perl /path/variant_effect_predictor.pl -i input.vcf -o /path/output_dir/output.vcf

or, if a chromosome region and tabix is provided

> /path/tabix -h 1:10000000-20000000 | perl /path/variant_effect_predictor.pl -o /path/output_dir/output.vcf

=cut


package ReseqTrack::Tools::RunVep;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::FileSystemUtils qw (check_executable);

use base qw(ReseqTrack::Tools::RunProgram);

=head2 new

  Arg [-tabix]   :
      string, path of the tabix executable (not needed if tabix is in $PATH)
  Arg [-region]   :
      string, chromosomal coordinate
    + Arguments for ReseqTrack::Tools::RunProgram parent class

  Function  : Creates a new ReseqTrack::Tools::RunVep object.
  Returntype: ReseqTrack::Tools::RunVep
  Exceptions: 
  Example   : my $run_vep = ReseqTrack::Tools::RunVep->new(
                -input_files => '/path/vcf', 
                -program => ‘/path/variant_effect_predictor.pl’,
                -options => { ‘plugin’=> ‘AncestralAlleles,/nfs/1000g-work/G1K/work/avikd/test_hive/ancesteral_alleles_for_vep/’,
                              ‘regulatory’ => ““,
                              ‘offline’ => ““,
                              ‘gmaf’ => ““,
                            },
                -working_dir => '/path/output_dir/',
                -tabix => ‘/path/tabix’,
                -region => ‘1:10000000-20000000’, );
            
=cut


sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);
  
  my ( $tabix, $region ) = rearrange( [ qw( TABIX REGION ) ], @args);
  
  $self->tabix($tabix);
  $self->region($region);
          
  return $self;
}

=head2 run_program

  Arg [1]   : ReseqTrack::Tools::RunVep
  Function  : uses variant_effect_predictor.pl  to evaluate functional consequence of 
              variants in in $self->input_files.
              Output is files are stored in $self->output_files
              This method is called by the RunProgram run method
  Returntype: 
  Exceptions: 
  Example   : $self->run();

=cut

sub run_program{
    my ($self) = @_;
    
    throw "Can't accept multiple input files" if(scalar @{$self->input_files} >1);
    print $self->region;
    
    my @cmd_words;
    my $region_tag;
    
    if($self->region)
    {
      throw("could not extract regions from VCF without tabix")  if(!$self->tabix);
      check_executable($self->tabix);
      
      push(@cmd_words, $self->tabix ,"-h","-pvcf");
      push(@cmd_words, ${$self->input_files}[0]);
      push(@cmd_words, $self->region);
      push(@cmd_words, "|" );
      
      $region_tag=$self->region;
      $region_tag=~ s{:}{_}g;
      $region_tag=~ s{-}{_}g;
                
    }
    
    push(@cmd_words, "perl", $self->program,);
    
    push(@cmd_words, '--i',${$self->input_files}[0]) unless($self->region);
        
    my $output = $self->working_dir .'/';
    
    if($region_tag)
    {
      $output.= $region_tag;
    }
    else
    {    
      $output.=$self->job_name;
    }
        
    $output.='.vcf' if(exists $self->options->{vcf});
    
    $output =~ s{//}{/};
    
    $self->output_files($output);
    
    push(@cmd_words, '-o', $output);
    
    
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
           push(@cmd_words,"--".$tag." ".$value);
         }
         else
         {
           push(@cmd_words, "--". $tag);
         }
       }
    }

    my $cmd = join(' ', @cmd_words);
    $self->execute_command_line ($cmd);
    
    return $self;    
}

sub region {
    my ($self, $arg) = @_;
    if ($arg) {
        $self->{'region'} = $arg;
    }
    return $self->{'region'};
}

sub tabix {
    my ($self, $arg) = @_;
    if ($arg) {
        $self->{'tabix'} = $arg;
    }
    return $self->{'tabix'};
}

1;
