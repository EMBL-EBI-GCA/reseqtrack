
package AccessibleGenome::RunBamUtil;

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

  my ( $SQ, $region_start, $region_end )
    = rearrange( [ qw( SQ REGION_START REGION_END) ], @args);

  $self->SQ($SQ);
  $self->region_start($region_start);
  $self->region_end($region_end);

  return $self;
}



sub run_program{
  my ($self) = @_;

  my $input_bams = $self->input_files;
  my $SQ = $self->SQ;

  throw("expecting one file only") if @$input_bams !=1;
  throw("SQ not defined") if !defined $SQ;

  my $output_stats = $self->working_dir . '/' . $self->job_name . '.stats.gz';
  my $regions_list = $self->working_dir . '/' . $self->job_name . '.regions_list';

  $self->created_files($regions_list);
  open my $fh, '>', $regions_list or throw("could not open $regions_list $!");
  print $fh join("\t", $SQ, $self->region_start -1, $self->region_end -1), "\n";
  close $fh;

  my @cmd_words = ($self->program, 'stats');
  push(@cmd_words, '--in', $input_bams->[0]);
  push(@cmd_words, '--regionList', $regions_list);
  push(@cmd_words, '--cBaseQC', '/dev/stdout');
  push(@cmd_words, '2> /dev/null');
  push(@cmd_words, '| cut -f1,2,10,14,15');
  push(@cmd_words, '| gzip -c');
  push(@cmd_words, '>', $output_stats);

  my $cmd = join(' ', @cmd_words);

  $self->output_files($output_stats);
  print "SQ is $SQ\n";
  return if $SQ =~ /HLA|decoy/;
    
  $self->execute_command_line($cmd);

  return;
}

sub SQ {
    my ($self, $SQ) = @_;
    if (@_ > 1) {
        $self->{'SQ'} = $SQ;
    }
    return $self->{'SQ'};
}

sub region_start {
    my ($self, $region_start) = @_;
    if (@_ > 1) {
        $self->{'region_start'} = $region_start;
    }
    return $self->{'region_start'} // -1;
}

sub region_end {
    my ($self, $region_end) = @_;
    if (@_ > 1) {
        $self->{'region_end'} = $region_end;
    }
    return $self->{'region_end'} // -1;
}




1;

