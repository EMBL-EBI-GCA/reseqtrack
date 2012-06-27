
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
use File::Basename qw(fileparse);
use Env qw( @PATH );

use base qw(ReseqTrack::Tools::RunProgram);


=head2 new

  Arg [-reference]   :
      string, path to the reference genome
  Arg [-samtools]   :
      string, path to samtools executable, used if output_format is 'BAM'
  Arg [-output_format]   :
      string, either 'SAM' (default) or 'BAM'.
  Arg [-mate1_file]   :
      string, optional, path to the mate1 file.
  Arg [-mate2_file]   :
      string, optional, path to the mate2 file.
  Arg [-fragment_file]   :
      string, optional, path to the fragment file.
  Arg [-paired_length]   :
      integer, insert size used by alignment program
  Arg [-first_read]   :
      integer, start alignment from this read in the fastq file
  Arg [-last_read]   :
      integer, end alignment at this read in the fastq file
  Arg [-read_group_fields]   :
      hashref, values to go in the RG tag e.g. {'ID' => 1, 'LB' => 'my_library'}
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
                -mate1_file => '/path/to/mate1',
                -mate2_file => '/path/to/mate2',
                -fragment_file => '/path/to/fragment',
                -read_group_fields => {'ID' => 1, 'LB' => 'my_library'},
                -paired_length => 3000);

=cut

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my ( $reference, $samtools, $output_format,
        $mate1_file, $mate2_file, $fragment_file, $paired_length,
        $read_group_fields, $first_read, $last_read,
        )
        = rearrange( [
             qw( REFERENCE SAMTOOLS OUTPUT_FORMAT
             MATE1_FILE MATE2_FILE FRAGMENT_FILE PAIRED_LENGTH
             READ_GROUP_FIELDS FIRST_READ LAST_READ
             )
                    ], @args);

  $self->mate1_file($mate1_file);
  $self->mate2_file($mate2_file);
  $self->fragment_file($fragment_file);
  $self->reference($reference);
  $self->paired_length($paired_length);
  $self->read_group_fields($read_group_fields);
  $self->first_read($first_read);
  $self->last_read($last_read);
  $self->output_format($output_format || 'SAM');

  if (!$samtools) {
    if ($ENV{SAMTOOLS}) {
      $samtools = $ENV{SAMTOOLS} . '/samtools';
    }
  }
  $self->samtools($samtools || 'samtools');

  return $self;
}

=head2 run_program

  Arg [1]   : ReseqTrack::Tools::RunAlignment
  Function  : Calls the run_alignment method, which must be implemented by child class
              This method is called by the RunProgram run method
              Does various checks that are needed for running any alignment program
  Returntype: 
  Exceptions: 
  Example   : $self->run();

=cut

sub run_program {
    my ($self) = @_;

    if (!$self->fragment_file && !$self->mate1_file && !$self->mate2_file) {
      $self->assign_fastq_files;
    }

    my $output_format = $self->output_format;
    throw($output_format . " is not a valid output format")
      if ($output_format ne 'BAM' && $output_format ne 'SAM');

    if ($output_format eq 'BAM' && ! -x $self->samtools) {
      if ($self->samtools =~ m{/} || ! first {-x $_} map {$_.'/'.$self->samtools} @PATH) {
        throw "cannot find samtools executable ".$self->samtools;
      }
    }

    $self->run_alignment;

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
    $job_name = $&;
    if ($self->first_read && $self->last_read) {
      $job_name .= '_'.$self->first_read .'-'.$self->last_read;
    }
    $job_name .= $$;

    $self->job_name($job_name);
    return $job_name;
}




=head2 run_alignment

  Arg [1]   : ReseqTrack::Tools::RunAlignment
  Function  : each child object should implement a run_alignment method
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

=head2 get_fastq_cmd_string

  Arg [1]   : ReseqTrack::Tools::RunAlignment
  Arg [2]   : string, fastq_type, either 'frag', 'mate1' or 'mate2'
  Function  : makes the string that should be used on the command line by child class
              for the input fastq file for the alignment program.
              This might be just a file name e.g. '/path/to/file.fq'
              or something more complex e.g. '<(sed -n 1,1000p file.fq)'
  Returntype: string
  Exceptions: if there is no fastq file of the requested type
  Example   : $self->assign_files;

=cut


sub get_fastq_cmd_string {
  my ($self, $fastq_type) = @_;

  my $fastq = $fastq_type eq 'frag' ? $self->fragment_file
            : $fastq_type eq 'mate1' ? $self->mate1_file
            : $fastq_type eq 'mate2' ? $self->mate2_file
            : '';
  throw("no file for fastq_type $fastq_type") if (!$fastq);

  my $first_read = $self->first_read;
  my $last_read = $self->last_read;
  return $fastq if (!$first_read && !$last_read);

  my $first_line = $first_read ? $first_read * 4 - 3 : 1;
  my $last_line = $last_read ? $last_read * 4 : '$';
  my $cmd = ($fastq =~ /\.gz(ip)?$/)
        ? "<(gunzip -c $fastq | sed -n $first_line,${last_line}p)"
        : "<(sed -n $first_line,${last_line}p $fastq)";
  return $cmd;
}

=head2 get_static_fastq

  Arg [1]   : ReseqTrack::Tools::RunAlignment
  Arg [2]   : string, fastq_type, either 'frag', 'mate1' or 'mate2'
  Function  : This is an alternative to get_fastq_cmd_string
              For use by alignment programs that won't accept <(....) as an input
              Creates a temporary fastq file (if necessary) and adds it to created_files
  Returntype: string
  Exceptions: if there is no fastq file of the requested type
  Example   : $self->assign_files;

=cut

sub get_static_fastq {
  my ($self, $fastq_type) = @_;
  my $fastq = $fastq_type eq 'frag' ? $self->fragment_file
            : $fastq_type eq 'mate1' ? $self->mate1_file
            : $fastq_type eq 'mate2' ? $self->mate2_file
            : '';
  throw("no file for fastq_type $fastq_type") if (!$fastq);

  my $first_read = $self->first_read;
  my $last_read = $self->last_read;
  return $fastq if (!$first_read && !$last_read);

  my $first_line = $first_read ? $first_read * 4 - 3 : 1;
  my $last_line = $last_read ? $last_read * 4 : '$';

  my $temp_fastq = $self->working_dir . "/" . $self->job_name;
  $temp_fastq .= "_1" if ($fastq_type eq 'mate1');
  $temp_fastq .= "_2" if ($fastq_type eq 'mate2');
  $temp_fastq =~ s{//}{/}g;

  my $cmd = ($fastq =~ /\.gz(ip)?$/)
        ? "gunzip -c $fastq | sed -n $first_line,${last_line}p"
        : "sed -n $first_line,${last_line}p $fastq";
  $cmd .= " > $temp_fastq";

  $self->created_files($temp_fastq);
  $self->execute_cmd($cmd);

  return $temp_fastq;
}




=head2 accessor methods

  Arg [1]   : ReseqTrack::Tools::RunAlignment
  Arg [2]   : variable
  Function  : store variable in object
  Returntype: variable
  Exceptions: n/a
  Example   :

=cut


sub read_group_fields {
  my ($self, $hash) = @_;
  $self->{'read_group_fields'} ||= {};
  while (my ($tag, $value) = each %$hash) {
    $self->{'read_group_fields'}->{$tag} = $value;
  }
  return $self->{'read_group_fields'};
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

sub first_read {
    my $self = shift;
    if (@_) {
        $self->{first_read} = shift;
    }
    return $self->{first_read};
}

sub last_read {
    my $self = shift;
    if (@_) {
        $self->{last_read} = shift;
    }
    return $self->{last_read};
}

sub samtools {
    my $self = shift;
    if (@_) {
        $self->{samtools} = shift;
    }
    return $self->{samtools};
}

sub output_format {
    my $self = shift;
    if (@_) {
        $self->{output_format} = shift;
    }
    return $self->{output_format};
}

1;
