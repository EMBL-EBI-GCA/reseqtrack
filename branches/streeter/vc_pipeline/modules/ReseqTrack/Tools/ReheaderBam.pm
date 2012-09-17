=pod

=head1 NAME

ReseqTrack::Tools::ReheaderBam

=head1 SYNOPSIS

This is a class for running samtools to reheader a bam file
It is a sub class of a ReseqTrack::Tools::RunProgram.

=head1 Example


=cut

package ReseqTrack::Tools::ReheaderBam;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename qw(fileparse);
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists check_file_does_not_exist);

use base qw(ReseqTrack::Tools::RunProgram);


=head2 DEFAULT_OPTIONS

  Function  : Called by the RunProgram parent object in constructor
  Returntype: hashref
  Example   : my %options = %{&ReseqTrack::Tools:RunSamtools::DEFAULT_OPTIONS};

  option 'reuse_old_header': boolean, flag to use the header of the input bam as the basis for the output bam. Default 0.
  option 'replace_PG': boolean, flag to remove the old @PG lines from the input bam header. Default 0.
  option 'replace_CO': boolean, flag to remove the old @CO lines from the input bam header. Default 0.

=cut

sub DEFAULT_OPTIONS { return {
        'replace_PG' => 0,
        'replace_CO' => 0,
        'reuse_old_header' => 0,
        };
}

=head2 new

  Arg [-header_lines_file]   :
      string, path of a file containing new header lines
  Arg [-extra_header_lines]   :
      string or arrayref of strings.  Each string is a header line to be appended to the other header lines
  Arg [-SQ_fields_hash]   :
      hashref. Used to modify/update the existing SQ lines. e.g. {'SP' => 'human'}
  + Arguments for ReseqTrack::Tools::RunProgram parent class

  Function  : Creates a new ReseqTrack::Tools::ReheaderBam object.
  Returntype: ReseqTrack::Tools::ReheaderBam
  Exceptions: 
  Example   : my $run_reheader_bam = ReseqTrack::Tools::ReheaderBam->new(
                -input_files => '/path/to/bam',
                -header_lines_file => '/path/to/header_lines',
                -extra_header_lines => '@CO    a comment',
                -working_dir => '/path/to/dir/',
                -SQ_fields_hash => {'SP' => 'human'},
                -options => {'reuse_old_header' => 1, 'replace_PG' => 1});

=cut

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my ( $header_lines_file, $extra_header_lines, $SQ_fields_hash,
        )
    = rearrange( [
         qw( HEADER_LINES_FILE EXTRA_HEADER_LINES SQ_FIELDS_HASH
                ) ], @args);

  #setting defaults
  if (!$self->program) {
    if ($ENV{SAMTOOLS}) {
      $self->program($ENV{SAMTOOLS} . '/samtools');
    }
    else {
      $self->program('samtools');
    }
  }
  $self->extra_header_lines($extra_header_lines);
  $self->header_lines_file($header_lines_file);

  if (ref($SQ_fields_hash) eq 'HASH') {
    while (my ($tag, $value) = each %$SQ_fields_hash) {
      $self->SQ_fields_hash($tag, $value);
    }
  }

  return $self;
}

sub get_old_header{
  my $self = shift;
  my $bam = $self->input_files->[0];

  my $header_file = $self->working_dir . '/' . $self->job_name. '.bam.oldheader';
  $header_file =~ s{//}{/}g;

  my @cmd_words = ($self->program, 'view', '-H');
  push(@cmd_words, $bam);
  push(@cmd_words, '>', $header_file);
  my $cmd = join(' ', @cmd_words);

  $self->created_files($header_file);
  $self->execute_command_line($cmd);

  open my $FH, '<', $header_file
      or die "cannot open $header_file";
  my @header_lines = <$FH>;
  close $FH;

  return \@header_lines;
}

sub make_header_file {
  my $self = shift;

  my $header_lines;
  if ($self->options('reuse_old_header')) {
    $header_lines = $self->get_old_header;
  }

  if ($self->options('replace_PG')) {
    my @new_header_lines = grep {! /^\@PG/} @$header_lines;
    $header_lines = \@new_header_lines;
  }
  if ($self->options('replace_CO')) {
    my @new_header_lines = grep {! /^\@CO/} @$header_lines;
    $header_lines = \@new_header_lines;
  }

  if (my $header_lines_file = $self->header_lines_file) {
    open my $FH, '<', $header_lines_file
        or throw "cannot open $header_lines_file: $!";
    my @new_header_lines = <$FH>;
    close $FH;
    push(@$header_lines, @new_header_lines);
  }

  push(@$header_lines, @{$self->extra_header_lines});

  foreach my $SQ_line (grep {/^\@SQ/} @$header_lines) {
    chomp $SQ_line;
    FIELD:
    while (my ($tag, $value) = each %{$self->SQ_fields_hash}) {
      next FIELD if (!defined $value);
      $SQ_line =~ s/\t$tag:[^\t]*//g;
      $SQ_line .= "\t$tag:$value";
    }
  }

  foreach my $line (@$header_lines) {
    $line =~ s/\n*$/\n/;
  }

  my $header_file = $self->working_dir . '/' . $self->job_name. '.bam.newheader';
  $header_file =~ s{//}{/}g;

  $self->created_files($header_file);

  open my $FH, '>', $header_file
      or die "cannot open $header_file";
  print $FH @$header_lines;
  close $FH;

  return $header_file;

}

sub reheader_bam {
  my ($self, $header) = @_;
  my $input_bam = $self->input_files->[0];

  check_file_exists($header);

  my $output_bam = $self->working_dir . '/' . $self->job_name. '.rehead.bam';
  $output_bam =~ s{//+}{/}g;
  check_file_does_not_exist($output_bam);

  my @cmd_words = ($self->program, 'reheader', $header, $input_bam);
  push(@cmd_words, '>', $output_bam);
  my $cmd = join(' ', @cmd_words);

  $self->output_files($output_bam);
  $self->execute_command_line($cmd);
}

sub run_program {
    my ($self) = @_;

    my $new_header_file;
    if ($self->options('reuse_old_header') && ! @{$self->extra_header_lines} && ! keys %{$self->SQ_fields_hash}) {
      $new_header_file = $self->header_lines_file;
    }
    else {
      $new_header_file = $self->make_header_file;
    }

    $self->reheader_bam($new_header_file);

    return;
}

sub header_lines_file {
  my ($self, $header_lines_file) = @_;
  if ($header_lines_file) {
    $self->{'header_lines_file'} = $header_lines_file;
  }
  return $self->{'header_lines_file'};
}

sub extra_header_lines {
  my ($self, $arg) = @_;
  $self->{'extra_header_lines'} ||= [];
  if ($arg) {
    push(@{$self->{'extra_header_lines'}}, ref($arg) eq 'ARRAY' ? @$arg : $arg);
  }
  return $self->{'extra_header_lines'};
}

sub SQ_fields_hash {
  my ($self, $tag, $value) = @_;
  $self->{'SQ_fields_hash'} ||= {};
  if ($tag) {
    $self->{'SQ_fields_hash'}->{$tag} = $value;
  }
  return $self->{'SQ_fields_hash'};
}



1;

