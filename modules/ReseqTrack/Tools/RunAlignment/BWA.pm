=pod

=head1 NAME

=head1 SYNOPSIS

=head1 Example

=cut

package ReseqTrack::Tools::RunAlignment::BWA;

use strict;
use warnings;
use vars qw(@ISA);


use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::RunAlignment;

@ISA = qw(ReseqTrack::Tools::RunAlignment);



=head2 new

  Arg [1]   : ReseqTrack:Tools::RunAlignment::BWA
  Arg [2]   : string, options for the single ended bwa mod
  Arg [3]   : string, options for paired end bwa mod
  Function  : create a BWA object, defaults program name to bwa and
  defaults to using ReseqTrack::Tools::RunAlignment::options for sampe options 
  if sampe options aren't defined and options are
  Returntype: ReseqTrack::Tools::RunAlignment::BWA
  Exceptions: n/a
  Example   : my $bwa = ReseqTrack::Tools::RunAlignment::BWA(
                      -program => "bwa",
                      -input => '/path/to/file'
                      -reference => '/path/to/reference',
                      -samtools => '/path/to/samtools',
                      );
=cut


sub new{
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  my ($samse_options, $sampe_options, $aln_options) =
      rearrange([qw(SAMSE_OPTIONS SAMPE_OPTIONS ALN_OPTIONS)], @args);

  $self->samse_options($samse_options);
  $self->sampe_options($sampe_options);
  $self->aln_options($aln_options);

  #setting defaults
  $self->program('bwa') unless($self->program);
  $self->sampe_options("-a 2000");
  $self->aln_options("-q 20 -l 32");
  #

  return $self;
}


=head2 accessor methods

  Arg [1]   : ReseqTrack::Tools::RunAlignment::BWA
  Arg [2]   : string, options string
  Function  : samse_options and sampe_options are specific bwa commandline
  options and these are the accessor methods for those variables
  Returntype: string
  Exceptions: n/a
  Example   : my $options = $self->samse_options;

=cut


sub samse_options{
  my ($self, $arg) = @_;
  unless($self->{samse_options}){
    $self->{samse_options} = '';
  }
  if($arg){
    $self->{samse_options} = $arg;
  }
  return $self->{samse_options};
}

sub sampe_options{
  my ($self, $arg) = @_;
  if($arg){
    $self->{sampe_options} = $arg;
  }
  return $self->{sampe_options};
}

sub aln_options{
  my ($self, $arg) = @_;
  if($arg){
    $self->{aln_options} = $arg;
  }
  return $self->{aln_options};
}


sub run{
  my ($self) = @_;
  $self->change_dir();
  my ($single_ended_bam, $paired_end_bam);
  if($self->fragment_file){
    $single_ended_bam  = $self->run_samse_alignment();
  }
  if($self->mate1_file && $self->mate2_file){
    $paired_end_bam = $self->run_sampe_alignment();
  }
  $self->output_files($single_ended_bam) if($single_ended_bam);
  $self->output_files($paired_end_bam) if($paired_end_bam);
  $self->delete_files;
}


sub run_sampe_alignment{
  my ($self) = @_;
  my $mate1_sai =  $self->run_aln_mode($self->mate1_file, $self->aln_options, "mate1");
  my $mate2_sai =  $self->run_aln_mode($self->mate2_file, $self->aln_options, "mate2");
  my $output_sam = $mate1_sai;
  $output_sam =~ s/\.mate1\.sai/\_pe\.sam/;
  my $sampe_cmd = $self->program." sampe ".$self->sampe_options ." ".
      $self->reference." ".$mate1_sai." ".$mate2_sai." ".$self->mate1_file." ".
      $self->mate2_file. " > ".$output_sam;
  print $sampe_cmd."\n";
  eval{
    my $exit = system($sampe_cmd);
    if($exit && $exit >= 1){
      throw("Failed to run ".$sampe_cmd);
    }
  };
  if($@){
    throw("Failed to run bwa sampe alignment $@");
  }
  $self->files_to_delete($output_sam);
  $self->files_to_delete($mate1_sai);
  $self->files_to_delete($mate2_sai);
  my $bam_file = $self->create_bam_from_sam($output_sam);
  return $bam_file;
}

sub run_samse_alignment{
  my ($self) = @_;
  my $sai_file = $self->run_aln_mode($self->fragment_file, $self->aln_options, "frag");
  $self->files_to_delete($sai_file);
  my $output_sam = $sai_file;
  $output_sam =~ s/\.frag\.sai/\_se\.sam/;
  my $samse_cmd = $self->program." samse ".$self->samse_options." ".
      $self->reference." ".$sai_file." ".$self->fragment_file." > ".$output_sam;
  print $samse_cmd."\n";
  eval{
    my $exit = system($samse_cmd);
    if($exit && $exit >= 1){
      throw("Failed to run ".$samse_cmd);
    }
  };
  if($@){
    throw("Failed to run bwa samse alignment $@");
  }
  $self->files_to_delete($output_sam);
  my $bam_file = $self->create_bam_from_sam($output_sam);
  return $bam_file;
}



sub run_aln_mode{
  my ($self, $input_file, $options, $file_ext) = @_;
  $options = $self->aln_options unless($options);
  my $output_file = $self->working_dir."/".$self->name;
  $output_file .= ".".$file_ext if($file_ext);
  $output_file .= ".sai";
  my $aln_command = $self->program." aln ".$options."  ".$self->reference." ".$input_file." > ".$output_file;
  print $aln_command."\n";
  eval{
    my $exit = system($aln_command);
    if($exit && $exit >= 1){
      throw("Failed to run ".$aln_command);
    }
  };
  if($@){
    throw("Failed to run bwa aln alignment $@");
  }
  return $output_file;
}

