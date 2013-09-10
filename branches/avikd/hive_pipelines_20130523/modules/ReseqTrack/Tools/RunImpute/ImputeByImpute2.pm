package ReseqTrack::Tools::RunImpute::ImputeByImpute2;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists check_executable);

use base qw(ReseqTrack::Tools::RunImpute);

=head2 new

  Arg [-parameters]   :
      hashref, command line options to use with “Impute2
            
  + Arguments for ReseqTrack::Tools::RunProgram parent class

  Function  : Creates a new ReseqTrack::Tools::RunImpute::ImputeByImpute2 object.
  Returntype: ReseqTrack::Tools::RunImpute::ImputeByImpute2
  Exceptions: 
  Example   : my $impute = ReseqTrack::Tools::RunImpute::ImputeByImpute2->new(
                -program        =>"/path/impute" ,
                -input_files            => ‘/path/to/gen|hap‘,
                -input_samples             => ‘/path/to/samples‘,
                -working_dir             => '/path/to/dir/',  ## this is working dir and output dir
                 -options              => {'Ne'=>20000 },
                -reference_config               => '/path/to/ref_config', ###format:chrom\tmap\thaps\tlegend\tsamples
                -chrom                   => '1',
                -region                  => '1-1000000',
                -use_phased              => 1,
                -phase_result            => 1
                );

=cut

sub DEFAULT_OPTIONS { return {
       'Ne'=>20000
        };
}

sub new {
	my ( $class, @args ) = @_;
	my $self = $class->SUPER::new(@args);
	
    my ( $input_samples, $use_phased, $phase_result )= rearrange([ qw( INPUT_SAMPLES USE_PHASED PHASE_RESULT )], @args); ### ?
    
    
    $self->input_samples($input_samples);
    $self->use_phased( $use_phased);
    $self->phase_result($phase_result);
        
    return $self;
}

=head2 run_program

  Arg [1]   : ReseqTrack::Tools::RunImpute::ImputeByImpute2
  Function  : uses Impute2 to impute in $self->input_files.
              Output is files are stored in $self->output_files
              This method is called by the RunProgram run method
  Returntype: 
  Exceptions: 
  Example   : $self->run();

=cut

sub run_program {
	my $self = shift;
    
    check_executable($self->program);
    check_file_exists($self->reference_config);
    check_file_exists($self->input_samples);
    
    throw "Impute2 doesn't accept multiple input files" if(scalar @{$self->input_files} >1);   ### ?
    
    
    my $output = $self->working_dir .'/'. $self->job_name;
    $output =~ s{//}{/};
    $self->output_files($output);
    
     my @cmd_words = ($self->program);
     
     if($self->use_phased)
     {
        push(@cmd_words, "-known_haps_g", ${$self->input_files}[0]);	
     }
     else
     {
            throw("No support for unphased inputs"); 
     }
     
     push(@cmd_words,"-sample_g", $self->input_samples) if($self->input_samples);
     
     if ( $self->options ) {
      OPTION:
      while (my ($tag, $value) = each %{$self->options}) {
        next OPTION if !defined $value;
      }
        foreach my $tag ( keys ( %{$self->options} ) ) {
            my $value = $self->options->{$tag};
            push(@cmd_words, "-".$tag." ".$value);
        }
    }
    
     $self->get_ref;
     
     push(@cmd_words,"-m",$self->ref_map);
     
     push(@cmd_words,"-h",$self->ref_hap);
     
     push(@cmd_words,"-l",$self->ref_legend);
     
     push(@cmd_words,"-sample_g_ref",$self->ref_samples) if($self->ref_samples); 
     
     push(@cmd_words,"-o",$output);
     
     my $region=$self->region;
     $region=~ s/-/ /; 
     push(@cmd_words,"-int",$region);
     
     if( $self->phase_result)
     {
         push(@cmd_words,"-phase");
     }
     else
     {
        throw("GEN to VCF conversion require phased output");
     }

     my $cmd = join(' ', @cmd_words);
     $self->execute_command_line ($cmd);

     return $self;
}

sub get_ref {
    my ( $self, $arg ) = @_;
    my $reference_config = $self->reference_config; ###format:chrom\tmap\thaps\tlegend\tsamples
    my $chrom= $self->chrom;
   
    if($chrom eq 'X' || $chrom eq 'chrX')
    {
     	throw("Support for ChrX coming soon");
    }
    else
    {
    	open(FH,$reference_config)||throw("can't read $reference_config\n");
    	while(my $ref_line=<FH>)
    	{
        	chomp($ref_line);
        	next if($ref_line=~ /^#/);
        
        	my @field=split(/\t/,$ref_line);
        
        	if($field[0]==$chrom)
        	{
           		$self->ref_map($field[1]);
           		$self->ref_hap($field[2]);
           		$self->ref_legend($field[3]);
           		#$self->ref_samples($field[4]) if($field[4]); ### samples column is not mandatory 
           		
        	}
	}
    }
    return $self->{get_ref};
}

sub ref_map {
   my ( $self, $arg ) = @_;

  if ($arg) {
    $self->{ref_map} = $arg;
  }
  return $self->{ref_map};
}

sub ref_hap {
   my ( $self, $arg ) = @_;

  if ($arg) {
    $self->{ref_hap} = $arg;
  }
  return $self->{ref_hap};
}

sub ref_legend {
   my ( $self, $arg ) = @_;

  if ($arg) {
    $self->{ref_legend} = $arg;
  }
  return $self->{ref_legend};
}

sub ref_samples {
   my ( $self, $arg ) = @_;

  if ($arg) {
    $self->{ref_samples} = $arg;
  }
  return $self->{ref_samples};
}

sub use_phased {
	my $self = shift;
    if (@_) {
        $self->{'use_phased'} = (shift) ? 1 : 0;
    }
    return $self->{'use_phased'};
}

sub phase_result {
	my $self = shift;
    if (@_) {
        $self->{'phase_result'} = (shift) ? 1 : 0;
    }
    return $self->{'phase_result'};
}

sub input_samples {
  my ( $self, $arg ) = @_;

  if ($arg) {
    $self->{input_samples} = $arg;
  }
  return $self->{input_samples};
}

1;

=pod

=head1 NAME

ReseqTrack::Tools::RunImpute::ImputeByImpute2

=head1 SYNOPSIS

This is a class for running Impute2 imputation
It is a sub class of a ReseqTrack::Tools::RunImpute

It generates and run command to impute variants from multiple diploid individuals

Reference file format: chrom<tab>map<tab>haps<tab>legend<tab>samples (samples is optional)

>impute2 -known_haps_g /path/phased/22_impute.haps -sample_g /path/phased/22_impute.samples -Ne 20000 -m /path/ref_panel/genetic_map_chr22_combined_b37.txt -h /path/ref_panel/chr22_impute.hap.gz -l /path/ref_panel/chr22_impute.legend.gz -o /path/run/22.31050408-36050407 -int 31050408 36050407 -phase
