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
     -regions          => ['chr4:10000-20000', 'chr3'],
     -program         => '/path/to/executable',
);


=cut

package ReseqTrack::Tools::RunTransposeBam;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);

use base qw(ReseqTrack::Tools::RunProgram);

sub DEFAULT_OPTIONS { return {
        'build_index' => 0,
        'uniquify_rg' => 0,
        'shorten_input_names' => 0, # should be used when the number of input files is very large
        };
}

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my ( $regions,
    )
    = rearrange( [
         qw( REGIONS )
          ], @args);

  $self->regions($regions);

  return $self;
}


sub run_program{
    my ($self) = @_;

    my $regions = $self->regions;
    throw("do not have regions") if !$regions || !@$regions;
    my $input_files = $self->input_files;
    throw("do not have any input files") if !@$input_files;
    if ( $self->options('shorten_input_names')) {
       $input_files = [values %{$self->get_short_input_names($self->options('shorten_input_names'))}];
    }

    my $output_bam = $self->working_dir . '/' . $self->job_name . '.transposed.bam';
    my $tmp_file = $self->working_dir . '/' . $self->job_name . '.tmp';

    open (OUT, ">", $tmp_file) || throw("Cannot open tmp file $tmp_file");
    foreach (sort @$input_files) {
        print OUT $_ . " ";
    }
    
    my @cmd_words = ('cat');	 
    push(@cmd_words, $tmp_file);
    #push @cmd_words,  ' | xargs -s 300000';
    push @cmd_words,  ' | xargs';
    push @cmd_words, $self->program;
    foreach my $region (@$regions) {
      push(@cmd_words, '-r', $region);
    }
    push(@cmd_words, '-i') if $self->options('build_index');
    push(@cmd_words, '-u') if $self->options('uniquify_rg');
    push(@cmd_words, '-o', $output_bam);

=head
    my @cmd_words = ($self->program);
    foreach my $region (@$regions) {
      push(@cmd_words, '-r', $region);
    }
    push(@cmd_words, '-i') if $self->options('build_index');
    push(@cmd_words, '-u') if $self->options('uniquify_rg');
    push(@cmd_words, '-o', $output_bam);
    push(@cmd_words, sort @$input_files);
=cut
    my $cmd = join(' ', @cmd_words);

    $self->output_files($output_bam);
    if ($self->options('build_index')) {
      $self->output_files("$output_bam.bai");
    }
    
    $self->execute_command_line($cmd);
    unlink $tmp_file;

    return;
}

sub regions {
    my ($self, $arg) = @_;
    if ($arg) {
        $self->{'regions'} = $arg;
    }
    return $self->{'regions'};
}

sub output_bai_file {
    my $self = shift;
    return (grep { /\.bai$/ } @{$self->output_files})[0];
}

sub output_bam_file {
    my $self = shift;
    return (grep { /\.bam$/ } @{$self->output_files})[0];
}



1;

