=pod

=head1 NAME

ReseqTrack::Tools::RunSplit

=head1 SYNOPSIS

This is a class for running split in reseqtrack/c_code
It is a sub class of a ReseqTrack::Tools::RunProgram.

=head1 Example


=cut

package ReseqTrack::Tools::RunSplit;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename qw(fileparse);
use File::Copy qw (move);

use base qw(ReseqTrack::Tools::RunProgram);

=head2 new

  Arg [-program]   :
      string, the 'split' executable
  Arg [-line_count_hash]   :
      hashref, max number of lines in the output files; keys are input file paths 
  Arg [-line_count]   :
      integer, max number of lines in all output files (overrides line_count_hash)
  Arg [-label_length]   :
      integer, number of characters used in output file name
  + Arguments for ReseqTrack::Tools::RunProgram parent class

  Function  : Creates a new ReseqTrack::Tools::RunSplit object.
  Returntype: ReseqTrack::Tools::RunSplit
  Exceptions: 
  Example   : my $run_split = ReseqTrack::Tools::RunSplit->new(
                -program => '/path/to/split',
                -working_dir => '/path/to/dir/',
                -line_count_hash => {'/path/file1' => 10, '/path/file2' => 20})

=cut


sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my ( $label_length, $line_count_hash, $line_count)
    = rearrange( [
         qw( LABEL_LENGTH LINE_COUNT_HASH LINE_COUNT )
          ], @args);

  $self->label_length($label_length);

  foreach my $file (keys %{$line_count_hash}) {
    $self->line_count($file, $line_count_hash->{$file});
    if (! scalar grep {$_ eq $file} @{$self->input_files}) {
      $self->input_files($file);
    }
  }
  if ($line_count) {
    foreach my $file (@{$self->input_files}) {
      $self->line_count($file, $line_count);
    }
  }

  return $self;
}


=head2 run

  Arg [1]   : ReseqTrack::Tools::RunSplit
  Function  : uses split to process the files in $self->input_files.
  Output files are stored in $self->output_files and $self->output_file_hash
  Returntype: 
  Exceptions: 
  Example   : $self->run();

=cut

sub run_program{
    my ($self) = @_;

    foreach my $file(@{$self->input_files}){
      throw ("must have a non-zero line count for $file") if (!$self->line_count($file));
    }

    foreach my $file(@{$self->input_files}){
      my ($basename, $input_dir, $suffix) = fileparse($file, qr/\.[^\/]+/);
      $suffix =~ s/\.gz$//;

      my $tmp_dir = $self->get_temp_dir;
      my $output_prefix = $tmp_dir . '/' . $basename . '.';
      $output_prefix =~ s{//}{/}g;

      my $cmd_line = $self->program();
      $cmd_line .= " -n " . $self->line_count($file);
      $cmd_line .= " -p $output_prefix ";
      $cmd_line .= "-s $suffix ";
      $cmd_line .= "-a ".$self->label_length." " if ($self->label_length);
      $cmd_line .= $file;

      $self->execute_command_line($cmd_line);

      opendir(my $DIR, $tmp_dir) or die "can't open $tmp_dir: $!";
      my @dir_files = readdir($DIR);
      closedir $DIR;
      DIR_FILE:
      foreach my $dir_file (@dir_files) {
        next DIR_FILE if (!($dir_file =~ /$basename\.(\d+)$suffix\.gz/));
        my $label = $1;
        my $dir_filepath = $tmp_dir . '/' . $dir_file;
        my $target_filepath = $self->working_dir . '/' . $dir_file;
        $dir_filepath =~ s{//}{/}g;
        $target_filepath =~ s{//}{/}g;

        $self->output_file_hash($file, $label, $target_filepath);
        $self->output_files($target_filepath);
        print "moving $dir_filepath to $target_filepath\n";
        move($dir_filepath, $target_filepath)
              or throw "move failed: $!";
      }
        
    }
}


sub line_count {
    my ($self, $file, $arg) = @_;
    if ($arg) {
        $self->{'line_count'}->{$file} = $arg;
    }
    return $self->{'line_count'}->{$file};
}

sub label_length {
    my ($self, $arg) = @_;
    if ($arg) {
        $self->{'label_length'} = $arg;
    }
    return $self->{'label_length'};
}

sub output_file_hash {
  my ( $self, $input_file, $label, $output_file ) = @_;

  if ($input_file && $label && $output_file) {
    $self->{'output_file_hash'}->{$input_file}->{$label} = $output_file;
  }
  return $self->{'output_file_hash'};
}



1;

