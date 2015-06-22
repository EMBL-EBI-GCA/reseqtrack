=pod

=head1 NAME

ReseqTrack::Tools::RunFreebayes

=head1 SYNOPSIS

This is a class for running Freebayes
It is a sub class of a ReseqTrack::Tools::RunProgram.

=head1 Example


=cut

package ReseqTrack::Tools::RunFreebayes;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::FileSystemUtils qw (check_executable);
use List::Util qw (first);
use Env qw( @PATH );

use base qw(ReseqTrack::Tools::RunProgram);


=head2 DEFAULT_OPTIONS

  Function  : Called by the RunProgram parent object in constructor
  Returntype: hashref
  Example   : my %options = %{&ReseqTrack::Tools:RunFreebayes::DEFAULT_OPTIONS};

=cut

sub DEFAULT_OPTIONS { return {
        
        };
}

sub CMD_MAPPINGS { return {
    'vcf_allelic_primitives'    => \&run_vcf_allelic_primitives,
    }
}

=head2 new

  Arg [-freebayes_dir]   :
      string,  freebayes home directory
  Arg [-bgzip]   :
      string, path of the bgzip executable (not needed if bgzip is in $PATH)
      
  + Arguments for ReseqTrack::Tools::RunProgram parent class

  Function  : Creates a new ReseqTrack::Tools::RunFreebayes object.
  Returntype: ReseqTrack::Tools::RunFreebayes
  Exceptions: 
  Example   : my $run_freebayes = ReseqTrack::Tools::RunFreebayes->new(
                -freebayes_dir => ‘/path/freebayes/‘,
                -input_files => ['/path/vcf1', '/path/vcf2'],
                -working_dir => '/path/to/dir/');

=cut

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my ( $freebayes_dir, $bgzip )    = rearrange( [ qw( FREEBAYES_DIR BGZIP ) ], @args);

  #setting defaults
  $self->freebayes_dir($freebayes_dir);
  $self->bgzip($bgzip || 'bgzip');
  return $self;
}

=head2 run_program

  Arg [1]   : ReseqTrack::Tools::RunFreebayes
  Arg [2]   : string, command
  Function  : uses Freebayes to process the files in $self->input_files.
              Output is files are stored in $self->output_files
              This method is called by the RunProgram run method
  Returntype: 
  Exceptions: Throws if the command is not recognised.
  Example   : $self->run(‘command‘);

=cut

sub run_program {
    my ( $self, $command ) = @_;

    throw("need a freebayes_dir") if !$self->freebayes_dir;
    throw("freebayes_dir is not a directory") if ( !-d $self->freebayes_dir );
    
    my $subs = CMD_MAPPINGS();

    throw("Did not recognise command $command")
      if ( !defined $subs->{$command} );

    my @returned_values = &{ $subs->{$command} }($self);
    return @returned_values;
}

sub run_vcf_allelic_primitives {
    
    my $self = shift;
    
    my $executable = $self->freebayes_dir . '/vcflib/vcfallelicprimitives';
    check_executable($executable);
    
    throw("expecting single vcf") if (@{$self->input_files} >1);
    
    my $output_vcf = $self->working_dir . '/' . $self->job_name . '.allelicprimitives.vcf.gz';
    
    my @cmd_words = ($executable);
    
    push(@cmd_words,'-t','wasHaplo');
    
    push(@cmd_words, @{$self->input_files});
    
    push(@cmd_words, '|', $self->bgzip, '-c');
    push(@cmd_words, '>', $output_vcf);
    
    my $cmd = join(' ', @cmd_words);

    $self->output_files($output_vcf);
    $self->execute_command_line($cmd);
    
}

sub freebayes_dir {
    my $self = shift;
    return $self->program(@_);
}

sub bgzip {
    my ($self, $arg) = @_;
    if ($arg) {
        $self->{'bgzip'} = $arg;
    }
    return $self->{'bgzip'};
}

1;
