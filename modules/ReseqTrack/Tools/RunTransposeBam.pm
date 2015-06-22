=pod

=head1 NAME

ReseqTrack::Tools::RunTransposeBam

=head1 SYNOPSIS

This is a class for running transpose_bam in reseqtrack/c_code
It is a sub class of a ReseqTrack::Tools::RunProgram.

=head1 Example

my $run_transpose_bam = $Reseqtrack::Tools::RunTransposeBam->new(
     -input_files     => ['/path/to/bam1', '/path/to/bam2'],
     -working_dir     => '/path/to/dir/',
     -region          => 'chr4:10000-20000',
     -program         => '/path/to/executable',
);


=cut

package ReseqTrack::Tools::RunTransposeBam;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);

use base qw(ReseqTrack::Tools::RunProgram);


sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my ( $region,
    )
    = rearrange( [
         qw( REGION )
          ], @args);

  $self->region($region);

  return $self;
}


sub run_program{
    my ($self) = @_;

    my $region = $self->region;
    throw("do not have a region") if !$region;

    my $output_bam = $self->working_dir . '/' . $self->job_name . '.transposed.bam';
    my @cmd_words = ($self->program);
    push(@cmd_words, '-r', $region);
    push(@cmd_words, '-o', $output_bam);
    push(@cmd_words, @{$self->input_files});

    my $cmd = join(' ', @cmd_words);
    $self->output_files($output_bam);
    $self->execute_command_line($cmd);

    return;
}

sub region {
    my ($self, $arg) = @_;
    if ($arg) {
        $self->{'region'} = $arg;
    }
    return $self->{'region'};
}


1;

