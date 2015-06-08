=pod

=head1 NAME

ReseqTrack::Tools::RunBamSqueeze

=head1 SYNOPSIS

This is a class for running the squeeze function of bamUtil
It is a sub class of a ReseqTrack::Tools::RunProgram.

=head1 Example


=cut

package ReseqTrack::Tools::RunBamSqueeze;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename qw(fileparse);

use base qw(ReseqTrack::Tools::RunProgram);


sub DEFAULT_OPTIONS { return {
        'keepOQ' => 1,
        'keepDups' => 1,
        };
}


sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my ( $rm_tag_types)
    = rearrange( [
         qw( RM_TAG_TYPES )
          ], @args);

  $self->rm_tag_types($rm_tag_types);

  return $self;
}



sub run_program{
  my ($self) = @_;

  foreach my $input_bam (@{$self->input_files}){
    my $prefix = fileparse($input_bam, qr/\.bam/ );
    my $output_bam = $self->working_dir . "/$prefix.squeezed.bam";
    $output_bam =~ s{//}{/}g;

    my @cmd_words = ($self->program, 'squeeze');
    push(@cmd_words, '--in', $input_bam);
    push(@cmd_words, '--out', $output_bam);

    my $tag_types = $self->rm_tag_types;
    if (@$tag_types) {
      push(@cmd_words, '--rmTags');
      push(@cmd_words, "'" . join(';', @$tag_types) . "'");
    }

    if ($self->options('keepOQ')) {
      push(@cmd_words, '--keepOQ');
    }
    if ($self->options('keepDups')) {
      push(@cmd_words, '--keepDups');
    }

    my $cmd = join(' ', @cmd_words);

    $self->output_files($output_bam);

    $self->execute_command_line($cmd);
  }
  return;
}


sub rm_tag_types {
  my ($self, $arg) = @_;

  $self->{'rm_tag_types'} ||= [];
  if ($arg) {
    foreach my $tag_type (@{ref($arg) eq 'ARRAY' ? $arg : [$arg]}) {
      throw("rm_tag_type does not have a valid format") if ($tag_type !~ /[A-Za-z][A-Za-z0-9]:[AifZHB]/);
      push(@{$self->{'rm_tag_types'}}, $tag_type);
    }
  }
  return $self->{'rm_tag_types'};
}


1;

