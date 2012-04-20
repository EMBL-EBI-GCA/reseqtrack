
=pod

=head1 NAME

ReseqTrack::Tools::RunAlignment

=head1 SYNOPSIS

This is a base class for RunAlignment objects.
It is a sub class of a ReseqTrack::Tools::RunProgram.
It provides methods that are common to many RunAlignment child classes.
Child classes should wrap specific alignment algorithms.

=head1 Example


=cut

package ReseqTrack::Tools::RunAlignment;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::SequenceIndexUtils qw(assign_files);
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists);
use File::Basename;
use ReseqTrack::Tools::RunSamtools;

use base qw(ReseqTrack::Tools::RunProgram);


=head2 new

  Arg [-reference]   :
      string, path to the reference genome
  Arg [-ref_samtools_index]   :
      string, path to the reference genome index file for samtools
  Arg [-samtools]   :
      string, the samtools executable (if in $PATH) or path to samtools
  Arg [-mate1_file]   :
      string, optional, path to the mate1 file.
  Arg [-mate2_file]   :
      string, optional, path to the mate2 file.
  Arg [-fragment_file]   :
      string, optional, path to the fragment file.
  Arg [-convert_sam_to_bam]   :
      boolean, default 1, flag to convert output to bam
  Arg [-sort_bams]   :
      boolean, default 1, flag to sort output bams
  Arg [-merge_bams]   :
      boolean, default 1, flag to merge all output to a single bam
  Arg [-index_bams]   :
      boolean, default 1, flag to index output bam(s)
  + Arguments for ReseqTrack::Tools::RunProgram parent class

  Function  : Creates a new ReseqTrack::Tools::RunAlignment object.
  If mate1_file, mate2_file, fragment_file are not specified, they are assigned
  from the list of input_files
  Returntype: ReseqTrack::Tools::RunAlignment
  Exceptions: 
  Example   : my $run_alignment = ReseqTrack::Tools::RunAlignment->new(
                -program => "alignmentprogram",
                -working_dir => '/path/to/dir/',
                -reference => '/path/to/ref.fa',
                -samtools => '/path/to/samtools',
                -mate1_file => '/path/to/mate1',
                -mate2_file => '/path/to/mate2',
                -fragment_file => '/path/to/fragment' );

=cut

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my ( $reference, $ref_samtools_index, $samtools,
        $mate1_file, $mate2_file, $fragment_file, $paired_length,
        $read_group_fields,
        $merge_bams, $sort_bams, $index_bams, $convert_sam_to_bam)
        = rearrange( [
             qw( REFERENCE REF_SAMTOOLS_INDEX SAMTOOLS
             MATE1_FILE MATE2_FILE FRAGMENT_FILE PAIRED_LENGTH
             READ_GROUP_FIELDS
             MERGE_BAMS SORT_BAMS INDEX_BAMS CONVERT_SAM_TO_BAM)
                    ], @args);

  $self->mate1_file($mate1_file);
  $self->mate2_file($mate2_file);
  $self->fragment_file($fragment_file);
  $self->reference($reference);
  $self->ref_samtools_index($ref_samtools_index);
  $self->samtools($samtools);
  $self->paired_length($paired_length);
  $self->read_group_fields($read_group_fields);
  $self->flag_sam_to_bam($convert_sam_to_bam);
  $self->flag_merge_bams($merge_bams);
  $self->flag_sort_bams($sort_bams);
  $self->flag_index_bams($index_bams);

  # This is a property of the alignment module.  Default value is 1.
  $self->sams_have_header(1);

  return $self;
}

sub run_program {
    my ($self) = @_;

    if ($self->flag_sam_to_bam) {
      $self->setup_samtools_object;
    }
    if (!$self->fragment_file && !$self->mate1_file && !$self->mate2_file) {
      $self->assign_fastq_files;
    }

    $self->run_alignment;

    if (!$self->flag_sam_to_bam) {
      $self->output_files($self->sam_files);
    }
    else {
      $self->samtools_object->input_files( $self->sam_files );
      $self->samtools_object->run();
      $self->output_files( $self->samtools_object->output_files );
    }
}

=head2 run_samtools

  Arg [1]   : ReseqTrack::Tools::RunAlignment
  Function  : uses samtools to process files listed in $self->sam_files
  Returntype: 
  Exceptions: 
  Example   : $self->run_samtools;

=cut

sub setup_samtools_object {
  my $self = shift;

  my $sort_bams = $self->flag_sort_bams || $self->flag_merge_bams || $self->index_bams;
  my $samtools_object = ReseqTrack::Tools::RunSamtools->new(
                      -program                 => $self->samtools,
                      -job_name                => $self->job_name,
                      -reference_index         => $self->ref_samtools_index,
                      -reference               => $self->reference,
                      -flag_merge              => $self->flag_merge_bams,
                      -flag_sort               => $sort_bams,
                      -flag_index              => $self->flag_index_bams,
                      -flag_sam_to_bam         => 1,
                      -flag_use_header         => $self->sams_have_header,
                      );
  check_file_exists($samtools_object->program);
  if (!$self->sams_have_header) {
    $self->samtools_object->find_reference_index;
  }
  $self->samtools_object($samtools_object);
}


=head2 generate_job_name

  Arg [1]   : ReseqTrack::Tools::RunAlignment
  Function  : overrides the generate_job_name of the parent class
  Returntype: string
  Exceptions: 
  Example   : $self->generate_job_name();

=cut

sub generate_job_name {
    my $self = shift;

    my $file = ${$self->input_files}[0];
    my $job_name = basename($file);
    $job_name =~ /\w+/;
    $job_name = $& . $$;

    $self->job_name($job_name);
    return $job_name;
}



=head2 sam_files

  Arg [1]   : ReseqTrack::Tools::RunAlignent
  Arg [2]   : string or arrayref of strings
  Function  : accessor method for sam files.
  Returntype: arrayref of strings
  Exceptions: 
  Example   : $self->sam_files('path/to/file');

=cut

sub sam_files {
  my ( $self, $arg ) = @_;

  $self->{'sam_files'} ||= [];
  if ($arg) {
    push( @{ $self->{'sam_files'} }, ref($arg) eq 'ARRAY' ? @$arg : $arg);
    $self->created_files($arg, !$self->flag_sam_to_bam);
  }
  return $self->{'sam_files'};
}

=head2 run

  Arg [1]   : ReseqTrack::Tools::RunAlignment
  Function  : each child object should implement a run method
  Returntype: n/a
  Exceptions: throws as this method should be implement in the child class
  Example   : 

=cut

sub run_alignment {
  my ($self) = @_;
  throw(  $self
          . " must implement a run_alignment method as ReseqTrack::Tools::RunAlignment "
          . "does not provide one" );
}


=head2 assign_fastq_files

  Arg [1]   : ReseqTrack::Tools::RunAlignment
  Function  : assign file names into mate1, mate2 and frag
  Returntype: array of filepaths
  Exceptions: if the path doesn't match any of the regexs
  Example   : $self->assign_files;

=cut

sub assign_fastq_files {
  my ($self) = @_;

  my @files = @{$self->input_files};
  my ( $mate1, $mate2, $frag ) ;

  if ((scalar @files) ==1) {
    $frag = $files[0];
  }
  else {
    ($mate1, $mate2, $frag) =
        ReseqTrack::Tools::SequenceIndexUtils::assign_files(\@files);
  }
  $self->fragment_file($frag);
  $self->mate1_file($mate1);
  $self->mate2_file($mate2);

  return ( $mate1, $mate2, $frag );
}

=head2 output_bam_files

  Arg [1]   : ReseqTrack::Tools::RunAlignment
  Function  : gets bam files from the list of output files
  Returntype: arrayref of filepaths
  Exceptions: 
  Example   : my $bams = $self->output_bam_files;

=cut

sub output_bam_files {
  my $self = shift;
  my @files = grep {/\.[bs]am$/} @{$self->output_files};
  return \@files;
}

=head2 output_bai_files

  Arg [1]   : ReseqTrack::Tools::RunAlignment
  Function  : gets bai files from the list of output files
  Returntype: arrayref of filepaths
  Exceptions: 
  Example   : my $bai_files = $self->output_bai_files;

=cut

sub output_bai_files {
  my $self = shift;
  my @files = grep {/\.bai$/} @{$self->output_files};
  return \@files;
}


sub read_group_fields {
  my ($self, $hash) = @_;
  $self->{'read_group_fields'} ||= {};
  while (my ($tag, $value) = each %$hash) {
    $self->{'read_group_fields'}->{$tag} = $value;
  }
  return $self->{'read_group_fields'};
}



=head2 accessor methods

  Arg [1]   : ReseqTrack::Tools::RunAlignment
  Arg [2]   : variable
  Function  : store variable in object
  Returntype: variable
  Exceptions: n/a
  Example   :

=cut

sub samtools {
    my ($self, $samtools) = @_;
    if ($samtools) {
        $self->{'samtools'} = $samtools;
    }
    return $self->{'samtools'};
}

sub samtools_object {
    my ($self, $samtools_object) = @_;
    if ($samtools_object) {
        $self->{'samtools_object'} = $samtools_object;
    }
    if (!$self->{'samtools_object'}) {
      $self->setup_samtools_object;
    }
    return $self->{'samtools_object'};
}

sub ref_samtools_index {
    my ($self, $ref_samtools_index) = @_;
    if ($ref_samtools_index) {
        $self->{'ref_samtools_index'} = $ref_samtools_index;
    }
    return $self->{'ref_samtools_index'};
}

sub reference {
    my ($self, $reference) = @_;
    if ($reference) {
        $self->{'reference'} = $reference;
    }
    return $self->{'reference'};
}

sub fragment_file {
  my ($self, $fragment_file) = @_;
  if ($fragment_file) {
    $self->{'fragment_file'} = $fragment_file;
    $self->input_files($fragment_file);
  }
  return $self->{'fragment_file'};
}

sub mate1_file {
  my ($self, $mate1_file) = @_;
  if ($mate1_file) {
    $self->{'mate1_file'} = $mate1_file;
    $self->input_files($mate1_file);
  }
  return $self->{'mate1_file'};
}

sub mate2_file {
  my ($self, $mate2_file) = @_;
  if ($mate2_file) {
    $self->{'mate2_file'} = $mate2_file;
    $self->input_files($mate2_file);
  }
  return $self->{'mate2_file'};
}

sub paired_length {
    my $self = shift;

    if (@_) {
        $self->{paired_length} = shift;
    }
    return $self->{paired_length};
}

sub flag_merge_bams {
    my $self = shift;
    if (@_) {
        $self->{flag_merge_bams} = shift;
    }
    return $self->{flag_merge_bams};
}
sub flag_sort_bams {
    my $self = shift;
    if (@_) {
        $self->{flag_sort_bams} = shift;
    }
    return $self->{flag_sort_bams};
}
sub flag_index_bams {
    my $self = shift;
    if (@_) {
        $self->{flag_index_bams} = shift;
    }
    return $self->{flag_index_bams};
}
sub flag_sam_to_bam {
    my $self = shift;
    if (@_) {
        $self->{flag_sam_to_bam} = shift;
    }
    return $self->{flag_sam_to_bam};
}

sub sams_have_header {
    my $self = shift;
    if (@_) {
        $self->{sams_have_header} = shift;
    }
    return $self->{sams_have_header};
}

1;
