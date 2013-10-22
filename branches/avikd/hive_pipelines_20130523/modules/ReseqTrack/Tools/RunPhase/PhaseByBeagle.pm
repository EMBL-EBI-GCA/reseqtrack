=pod

=head1 NAME

ReseqTrack::Tools::RunPhase::PhaseByBeagle

=head1 SYNOPSIS

This is a class for running Beagle4 Phasing (and Imputation)
It is a sub class of a ReseqTrack::Tools::RunPhase

It generates and run command to phase variants from multiple diploid individuals

>java -Xmx4000m -jar /path/b4.jar gtgl=in.vcf.gz chrom=22:16051249-16344141 impute=false out=/path/22.16051249-16344141  

Reference file format: chrom<tab>vcf<tab>map (map is optional)

=cut

package ReseqTrack::Tools::RunPhase::PhaseByBeagle;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists check_executable);

use base qw(ReseqTrack::Tools::RunPhase);

=head2 new

  Arg [-parameters]   :
      hashref, command line options to use with “beagle
            
  + Arguments for ReseqTrack::Tools::RunProgram parent class

  Function  : Creates a new ReseqTrack::Tools::RunPhase::PhaseByBeagle object.
  Returntype: ReseqTrack::Tools::RunPhase::PhaseByBeagle
  Exceptions: 
  Example   : my $varPhase_byBeagle = ReseqTrack::Tools::RunPhase::PhaseByBeagle->new(
                -java_exe        =>"/usr/bin/java" ,
                -jvm_args        =>"-Xmx4000m ",
                -input_files             => ‘/path/to/vcf‘
                -program                => ‘/path/to/beagle.jar’,
                -working_dir             => '/path/to/dir/',  ## this is working dir and output dir
                -reference_config               => '/path/to/ref_config',
                -options              => {	‘window‘=>24000, 
                								‘overlap‘=>3000, 
                								‘gprobs‘=>’true’ }
                -chrom                   => '1',
                -region                  => '1-1000000',
                -gl                      => 1
                );

=cut

# This is to set default parameters
sub DEFAULT_OPTIONS { return {
       'impute'=>'false'
        };
}


sub new {
	my ( $class, @args ) = @_;
	my $self = $class->SUPER::new(@args);

	my ( $java_exe, $jvm_args , $gl )= rearrange([ qw( JAVA_EXE JVM_ARGS GL )], @args); ### ?
    
     $self->java_exe($java_exe  || 'java' );
     $self->jvm_args( defined $jvm_args ? $jvm_args : '-Xmx4g' );
     $self->gl($gl);
     return $self;
}

=head2 run_program

  Arg [1]   : ReseqTrack::Tools::RunPhase::PhaseByBeagle
  Function  : uses Beagle4 to phase variants in $self->input_files.
              Output is files are stored in $self->output_files
              This method is called by the RunProgram run method
  Returntype: 
  Exceptions: 
  Example   : $self->run();

=cut

sub run_program {
    my ($self) = @_;
    
    check_file_exists($self->reference_config) if $self->reference_config; # ref is optional
    check_executable($self->java_exe);

    my $output_vcf = $self->working_dir .'/'. $self->job_name;
    $output_vcf =~ s{//}{/};
    
    my @cmd_words = ($self->java_exe, $self->jvm_args, '-jar');
    push(@cmd_words,$self->program);
	   	
    throw "Can't accept multiple input files" if(scalar @{$self->input_files} >1);   ### ?
		
    if($self->gl == 1)
    {	
	        push(@cmd_words, 'gl='.${$self->input_files}[0]);
    }
    else
    {
                push(@cmd_words, 'gtgl='.${$self->input_files}[0]);
    }
    
    if (my $region = $self->chrom) {
      $region .= ":".$self->region if ($self->region);
      push(@cmd_words, 'chrom='.$region);
    }
    
    if ( $self->options ) {
      OPTION:
      while (my ($tag, $value) = each %{$self->options}) {
        next OPTION if !defined $value;
      }
        foreach my $tag ( keys ( %{$self->options} ) ) {
            my $value = $self->options->{$tag};
            push(@cmd_words, $tag."=".$value);
        }
    }
 
    push(@cmd_words, 'out='.$output_vcf);
    $output_vcf.= '.vcf.gz';
   
    if($self->reference_config)
    {
	    $self->get_ref;
	
        push(@cmd_words,"ref=".$self->ref_vcf);
	
	    push(@cmd_words,"map=".$self->ref_map) if($self->ref_map);

    }    

    my $cmd = join(' ', @cmd_words);
    $self->output_files($output_vcf);
    
    #print $cmd,"\n";
    $self->execute_command_line ($cmd);
    
    return $self;
}
sub java_exe {
    my ($self, $arg) = @_;
    if ($arg) {
        $self->{'java_exe'} = $arg;
    }
    return $self->{'java_exe'};
}
sub jvm_args {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->{jvm_args} = $arg;
	}
	return $self->{jvm_args};
}
sub gl {
    my $self = shift;
    if (@_) {
        $self->{'gl'} = (shift) ? 1 : 0;
    }
    return $self->{'gl'};
}

sub get_ref {
    my ( $self, $arg ) = @_;
    my $reference_config = $self->reference_config; ###format:chrom\tvcf\tmap  OR format:chrom\tvcf
    my $chrom= $self->chrom;

    open(FH,$reference_config)||throw("can't read $reference_config\n");
    while(my $ref_line=<FH>)
    {
        chomp($ref_line);
        next if($ref_line=~ /^#/);

        my @field=split(/\t/,$ref_line);

        if($field[0] eq $chrom)
        {
           $self->ref_vcf($field[1]);
           $self->ref_map($field[2]) if($field[2]);
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

sub ref_vcf {
   my ( $self, $arg ) = @_;

  if ($arg) {
    $self->{ref_vcf} = $arg;
  }
  return $self->{ref_vcf};
}

1;

