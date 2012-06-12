=pod

=head1 NAME

ReseqTrack::Tools::RunValidateBam

=head1 SYNOPSIS

This is a class for running validate_bam in reseqtrack/c_code
It is a sub class of a ReseqTrack::Tools::RunProgram.

=head1 Example


=cut

package ReseqTrack::Tools::RunValidateBam;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw);

use base qw(ReseqTrack::Tools::RunProgram);


sub DEFAULT_OPTIONS { return {
        };
}

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  return $self;
}



sub run_program{
  my ($self) = @_;

  foreach my $input_bam (@{$self->input_files}){
    my $output_bas = "$input_bam.bas";
    $output_bas =~ s{//}{/}g;

    my $cmd = $self->program . " -o $output_bas $input_bam";

    $self->output_files($output_bas);

    $self->execute_command_line($cmd);
  }
  return;
}


1;

