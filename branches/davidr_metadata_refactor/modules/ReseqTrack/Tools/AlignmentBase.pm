package ReseqTrack::Tools::AlignmentBase;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::SequenceIndexUtils;
use Data::Dumper;
use File::Temp qw/ tempfile tempdir /;
use File::Basename;

sub new {
  my ($class, @args) = @_;

  my $self ={};

  bless $self,$class;
 
  my ($reference,
      $samtools,
      $bam,
      $sam,
      ) =

 rearrange([qw(
              REFERENCE
              SAMTOOLS
              BAM
              SAM
              )],
              @args);

  $self->reference($reference);
  $self->samtools($samtools);
  $self->bam($bam);
  $self->sam($sam);

  return $self;
}

sub reference{
  my ($self, $arg) = @_;
  if ($arg) {
    $self->{'reference'} = $arg;
  }
  return $self->{'reference'};
}

sub sam{
  my ($self, $arg) = @_;
  if ($arg) {
    $self->{'sam'} = $arg;
  }
  return $self->{'sam'};
}

sub bam{
  my ($self, $arg) = @_;
  if ($arg) {
    $self->{'bam'} = $arg;
  }
  return $self->{'bam'};
}

sub samtools{
  my ($self, $arg) = @_;
  if ($arg) {
    $self->{'samtools'} = $arg;
  }
  return $self->{'samtools'};
}

sub files_to_delete {
  my ( $self, $file ) = @_;
  if ($file) {
    if ( ref($file) eq 'ARRAY' ) {
      foreach my $path (@$file) {
        $self->{'files_to_delete'}->{$path} = 1;
      }
    } else {
      $self->{'files_to_delete'}->{$file} = 1;
    }
  }
  my @keys = keys( %{ $self->{'files_to_delete'} } );
  return \@keys;
}




=head2 delete_files

  Arg [1]   : ReseqTrack::Tools::RunAlignment
  Arg [2]   : arrayref of string of paths to delete
  Function  : remove given files
  Returntype: n/a
  Exceptions: throws if a file doesn't exist
  Example   :

=cut


sub delete_files{
  my ($self, $files) = @_;
  $files = $self->files_to_delete unless($files);
  foreach my $file (@$files) {
    print "Deleting ".$file."\n";
    unlink $file;
  }
}


######################
sub check_bai_present{
  my $self  = shift;

  my $bam = $self->bam;
  my $bai = $bam . ".bai";

  my $cmd = $self->samtools . " index " . $bam;


  if (! -e $bai){
    print "No bai file present. Creating\n";
    eval{
      `$cmd`;
    };

    if ($@) {
      throw( "Failed: $cmd" );      
    }

  }

 return $bai;
}



1;
