=pod

=head1 NAME

ReseqTrack::Tools::ReheaderBam

=head1 SYNOPSIS

This is a class for running samtools
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

=cut

sub DEFAULT_OPTIONS { return {
        'replace_PG' => 0,
        'replace_CO' => 0,
        };
}

=head2 new

  Arg [-reference_index]   :
      string, path of the reference genome .fai file
  Arg [-reference]   :
      string, path of the reference
  + Arguments for ReseqTrack::Tools::RunProgram parent class

  Function  : Creates a new ReseqTrack::Tools::RunSamtools object.
  Returntype: ReseqTrack::Tools::RunSamtools
  Exceptions: 
  Example   : my $run_alignment = ReseqTrack::Tools::RunSamtools->new(
                -input_files => ['/path/sam1', '/path/sam2'],
                -program => "samtools",
                -working_dir => '/path/to/dir/');

=cut

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my ( $header_file, $PG_fields_strings, $CO_fields_strings, $SQ_fields_hash,
        )
    = rearrange( [
         qw( HEADER_FILE PG_FIELDS_STRINGS CO_FIELDS_STRINGS SQ_FIELDS_HASH
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
  $self->PG_fields_strings($PG_fields_strings);
  $self->CO_fields_strings($CO_fields_strings);

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
  my $header_lines = $self->get_old_header;

  if ($self->options('replace_PG')) {
    my @new_header_lines = grep {! /^\@PG/} @$header_lines;
    $header_lines = \@new_header_lines;
  }
  if ($self->options('replace_CO')) {
    my @new_header_lines = grep {! /^\@CO/} @$header_lines;
    $header_lines = \@new_header_lines;
  }

  foreach my $SQ_line (grep {/^\@SQ/} @$header_lines) {
    chomp $SQ_line;
    FIELD:
    while (my ($tag, $value) = each %{$self->SQ_fields_hash}) {
      next FIELD if (!defined $value);
      $SQ_line =~ s/\t$tag:[^\t]//g;
      $SQ_line .= "\t$tag:$value";
    }
    $SQ_line .= "\n";
  }

  foreach my $PG_field(@{$self->PG_fields_strings}) {
    chomp $PG_field;
    my $PG_string = "\@PG\t$PG_field\n";
    push(@$header_lines, $PG_string);
  }

  foreach my $CO_field (@{$self->CO_fields_strings}) {
    chomp $CO_field;
    my $CO_string = "\@CO\t$CO_field\n";
    push(@$header_lines, $CO_string);
  }

  my $header_file = $self->working_dir . '/' . $self->job_name. '.bam.newheader';
  $header_file =~ s{//}{/}g;

  $self->created_files($header_file);
  $self->header_file($header_file);

  open my $FH, '>', $header_file
      or die "cannot open $header_file";
  print $FH @$header_lines;
  close $FH;

}

sub reheader_bam {
  my $self = shift;
  my $input_bam = $self->input_files->[0];
  my $header = $self->header_file;

  check_file_exists($header);

  my $output_bam = $self->working_dir . '/' . $self->job_name. '.rehead.bam';
  $output_bam =~ s{//}{/}g;
  check_file_does_not_exist($output_bam);

  my @cmd_words = ($self->program, 'reheader', $header, $input_bam);
  push(@cmd_words, '>', $output_bam);
  my $cmd = join(' ', @cmd_words);

  $self->output_files($output_bam);
  $self->execute_command_line($cmd);
}

=head2 run_program

  Arg [1]   : ReseqTrack::Tools::RunSamtools
  Arg [2]   : string, command, must be one of the following:
              merge, sort, index, fix_and_calmd, calmd, sam_to_bam
  Function  : uses samtools to process the files in $self->input_files.
              Output is files are stored in $self->output_files
              This method is called by the RunProgram run method
  Returntype: 
  Exceptions: Throws if the command is not recognised.
  Example   : $self->run();

=cut

sub run_program {
    my ($self) = @_;

    if (! $self->header_file) {
      $self->make_header_file;
    }

    $self->reheader_bam;

    return;
}

sub header_file {
  my ($self, $header_file) = @_;
  if ($header_file) {
    $self->{'header_file'} = $header_file;
  }
  return $self->{'header_file'};
}

sub PG_fields_strings {
  my ($self, $arg) = @_;
  $self->{'PG_fields_strings'} ||= [];
  if ($arg) {
    push(@{$self->{'PG_fields_strings'}}, ref($arg) eq 'ARRAY' ? @$arg : $arg);
  }
  return $self->{'PG_fields_strings'};
}

sub CO_fields_strings {
  my ($self, $arg) = @_;
  $self->{'CO_fields_strings'} ||= [];
  if ($arg) {
    push(@{$self->{'CO_fields_strings'}}, ref($arg) eq 'ARRAY' ? @$arg : $arg);
  }
  return $self->{'CO_fields_strings'};
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

