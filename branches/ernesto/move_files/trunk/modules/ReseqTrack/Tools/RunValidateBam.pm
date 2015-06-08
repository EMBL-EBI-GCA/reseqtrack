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
use ReseqTrack::Tools::Argument qw(rearrange);

use base qw(ReseqTrack::Tools::RunProgram);


sub DEFAULT_OPTIONS { return {
        };
}

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my ( $md5_file_hash )
    = rearrange( [ qw( MD5_FILE_HASH ) ], @args);

  while (my ($file, $md5) = each %$md5_file_hash) {
    $self->md5($file, $md5);
  }

  return $self;
}



sub run_program{
  my ($self) = @_;

  foreach my $input_bam (@{$self->input_files}){
    my $output_bas = "$input_bam.bas";
    $output_bas =~ s{//}{/}g;

    my @cmd_words = ($self->program, '-o', $output_bas);
    my $md5 = $self->md5($input_bam);
    if ($md5) {
      push(@cmd_words, '-m', $md5);
    }
    push(@cmd_words, $input_bam);

    my $cmd = join(' ', @cmd_words);

    $self->output_files($output_bas);
    $self->execute_command_line($cmd);
  }
  return;
}

sub md5 {
    my ($self, $file, $arg) = @_;
    if ($arg) {
        $self->{'md5'}->{$file} = $arg;
    }
    return $self->{'md5'}->{$file};
}



1;

