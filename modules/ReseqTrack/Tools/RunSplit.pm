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
  Output files are stored in $self->output_files and $self->grouped_output_files
  Returntype: 
  Exceptions: 
  Example   : $self->run();

=cut

sub run{
    my ($self) = @_;

    $self->change_dir();
    my $dir = $self->working_dir();

    foreach my $file(@{$self->input_files}){
      my ($basename, $input_dir, $suffix) = fileparse($file, qr/\.[^\/]+/);
      my $output_prefix = $dir . '/' . $basename . '.';
      $output_prefix =~ s{//}{/};
      $suffix =~ s/\.gz$//;

      my $line_count = $self->line_count($file);

      if ($line_count) {
        my $cmd_line = $self->program();
        $cmd_line .= " -n $line_count ";
        $cmd_line .= "-p $output_prefix ";
        $cmd_line .= "-s $suffix ";
        $cmd_line .= "-a ".$self->label_length." " if ($self->label_length);
        $cmd_line .= $file;

        $self->execute_command_line($cmd_line);

        opendir(my $DIR, $dir) or die "can't open $dir: $!";
        my @output_files = grep {/$basename\.\d+$suffix\.gz/} readdir($DIR);
        foreach my $output_file (@output_files) {
          $output_file = $dir . '/' . $output_file;
          $output_file =~ s{//}{/};
        }
        closedir $DIR;

        $self->output_files(\@output_files);
        $self->grouped_output_files($file, \@output_files);

      }
      else {
        my $output_file = $output_prefix . "1" . $suffix . ".gz";
        my $cmd_line = ($file =~ /\.gz$/)
              ? "cp $file $output_file"
              : "gzip -c $file > $output_file";
        $self->execute_command_line($cmd_line);
        $self->output_files($output_file);
        $self->grouped_output_files($file, $output_file);
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

sub grouped_output_files {
  my ( $self, $file, $arg ) = @_;

  if (! $self->{'grouped_output_files'}->{$file}) {
      $self->{'grouped_output_files'}->{$file} = [] ;
  }

  if ($arg) {
    if ( ref($arg) eq 'ARRAY' ) {
      push( @{ $self->{'grouped_output_files'}->{$file} }, @$arg );
    } else {
      push( @{ $self->{'grouped_output_files'}->{$file} }, $arg );
    }
  }
  return $self->{'grouped_output_files'}->{$file};
}



1;

