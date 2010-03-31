=pod

=head1 NAME

ReseqTrack::Tools::RunAlignment

=head1 SYNOPSIS

This is a base class for RunAlignment objects and provides some standard
accessor methods and throws exceptions when vital methods aren't implemented in
the child classes. The Child classes should wrap specific alignment algorithms

=head1 Example


=cut

package ReseqTrack::Tools::RunAlignment;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::SequenceIndexUtils;


=head2 new

  Arg [1]   : ReseqTrack::Tools::RunAlignment
  Arg [2]   : string, path to reference genome file/directory/stem as is required
  for aligner
  Arg [3]   : string/arrayref/ReseqTrack::File/ReseqTrack::Collection. The
  input argument can be several things, all must provide a filename or set of 
  filenames for the aligner to run on
  Arg [4]   : string, path to alignment program
  Arg [5]   : string, commandline options for program
  Arg [6]   : path to samtools installation so output bam can be
  indexed
  Function  : This should create the ReseqTrack::Tools::RunAlignment object 
  Returntype: ReseqTrack::Tools::RunAlignment
  Exceptions: 
  Example   : 

=cut


sub new {
  my ($class, @args) = @_;
  my $self ={};
  bless $self,$class;
  #reference genome file 
  #input sequence, we will allow single file name, list of file names
  # or a collection object with associated sequence
  #program
  #commandline options
  my ($reference, $input, $program, $options, $samtools, $working_dir, $name) = 
      rearrange([qw(REFERENCE INPUT PROGRAM OPTIONS SAMTOOLS WORKING_DIR NAME)], 
                @args);

  $self->reference($reference);
  $self->input($input);
  $self->program($program);
  $self->options($options);
  $self->samtools($samtools);
  $self->working_dir($working_dir);
  $self->name($name);
  #setting defaults
  $self->working_dir("/tmp/") unless($self->working_dir);
  unless($self->name){
    my $string = $self->mate1_file;
    $string = $self->fragment_file unless($string);
    $string =~ /^(\S+)\./;
    $self->name($1);
    throw("ReseqTrack::Tools::RunAlignment, Not sure what to do failed to define ".
          "a name from ".$string." for this run") unless($self->name);
  }
  #
  throw("Have no working directory") unless($self->working_dir && -d $self->working_dir);
  return $self
}



=head2 accessor methods

  Arg [1]   : ReseqTrack::Tools::RunAlignment
  Arg [2]   : string, generally
  Function  : These are accessor methods for variables needed to run
  and alignment like reference genomes, program and options
  Returntype: string, generally
  Exceptions: n/a
  Example   : my $reference = $run_alignment->reference;

=cut


sub reference{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'reference'} = $arg;
  }
  return $self->{'reference'};
}

sub program{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'program'} = $arg;
  }
  return $self->{'program'};
}

sub options{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'options'} = $arg;
  }
  return $self->{'options'};
}

sub samtools{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'samtools'} = $arg;
  }
  return $self->{'samtools'};
}

sub fragment_file{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'fragment_file'} = $arg;
  }
  return $self->{'fragment_file'};
}
sub mate1_file{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'mate1_file'} = $arg;
  }
  return $self->{'mate1_file'};
}
sub mate2_file{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'mate2_file'} = $arg;
  }
  return $self->{'mate2_file'};
}
sub name{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'name'} = $arg;
  }
  return $self->{'name'};
}
sub files_to_delete{
  my ($self, $file) = @_;
  if($file){
    if(ref($file) eq 'ARRAY'){
      foreach my $path(@$file){
        $self->{'files_to_delete'}->{$path} = 1;
      }
    }else{
      $self->{'files_to_delete'}->{$file} = 1;
    }
  }
  my @keys = keys(%{$self->{'files_to_delete'}});
  return \@keys;
}

=head2 working_dir

  Arg [1]   : ReseqTrack::Tools::RunAlignment
  Arg [2]   : string, path of working directory
  Function  : accessor method for working directory, also creates directory
  if it doesn't exist
  Returntype: string
  Exceptions: n/a
  Example   : my $work_dir = $self->working_dir;

=cut


sub working_dir{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'working_dir'} = $arg;
  }
  if($self->{'working_dir'}){
    unless(-d $self->{'working_dir'}){
      mkdir($self->{'working_dir'}, 775);
    }
    unless(-d $self->{'working_dir'}){
      throw("ReseqTrack::Tools::RunAlignment can run when ".$self->{'working_dir'}.
            " does not exist ");
    }
  }
  return $self->{'working_dir'};
}



=head2 change_dir

  Arg [1]   : ReseqTrack::Tools::RunAlignment
  Arg [2]   : string, path of working directory
  Function  : changes to given working directory
  Returntype: string
  Exceptions: throws if it can't change to given directory
  Example   : $self->check_dir();

=cut


sub change_dir{
  my ($self, $dir) = @_;
  $dir = $self->working_dir unless($dir);
  chdir($dir) or throw("Failed to change to ".$dir.
                       " ReseqTrack::Tools::RunAlignment check_dir");
  return $dir;
}

=head2 output_files

  Arg [1]   : ReseqTrack::Tools::RunAlignment
  Arg [2]   : string or arrayref of strings
  Function  : This is to store output filepaths
  Returntype: arrayref of strings
  Exceptions: n/a
  Example   : $self->output('path/to/file');

=cut



sub output_files{
  my ($self, $arg) = @_;
  $self->{'output'} = [] unless($self->{'output'});
  if($arg){
    if(ref($arg) eq 'ARRAY'){
      push(@{$self->{'output'}}, @$arg);
    }else{
      push(@{$self->{'output'}}, $arg);
    }
  }
      
  return $self->{'output'};
}


=head2 input

  Arg [1]   : ReseqTrack::Tools::RunAlignment
  Arg [2]   : string/arrayref/ReseqTrack::File/ReseqTrack::Collection. The
  input argument can be several things, all must provide a filename or set of 
  filenames for the aligner to run on
  Function  : setup the input file paths
  Returntype: string/arrayref/ReseqTrack::File/ReseqTrack::Collection, it returns
    what is passed in
  Exceptions: 
  Example   : 

=cut



sub input{
  my ($self, $input) = @_;
  #need to work out if we have filepaths, file objects or a collection object
  if($input){
    if(-e $input){
      $self->fragment_file($input);
    }elsif(ref($input) eq 'ARRAY'){
      my ($mate1, $mate2, $frag);
      if(-e $input->[0]){
        ($mate1, $mate2, $frag) = $self->assign_fastq_files($input);
      }elsif($input->[0]->isa("ReseqTrack::File")){
        my @names;
        foreach my $file(@$input){
          push(@names, $file->name);
        }
        ($mate1, $mate2, $frag) = $self->assign_fastq_files(\@names);
      }else{
        print STDERR "Not sure how to deal with the contents of ".$input."\n";
        foreach my $element(@$input){
          print STDERR $element."\n";
        }
        throw("ReseqTrack::Tools::RunAlignment::input Failed to process ".$input);
      }
      $self->fragment_file($frag);
      $self->mate1_file($mate1);
      $self->mate2_file($mate2);
    }elsif($input->isa("ReseqTrack::File")){
      $self->fragment_file($input->name);
    }elsif($input->isa('ReseqTrack::Collection')){
      throw("ReseqTrack::Tools::RunAlignment::input Can only handle file ".
            "collections not ".$input->table_name." collections")
          unless($input->table_name eq 'file');
      my $others = $input->others;
      my @names;
      foreach my $other(@$others){
        push(@names, $other->name);
      }
      my ($mate1, $mate2, $frag) = $self->assign_fastq_files(\@names);
      $self->fragment_file($frag);
      $self->mate1_file($mate1);
      $self->mate2_file($mate2);
    }else{
      throw("ReseqTrack::Tools::RunAlignment::input not sure how to handle input ".
            $input." what type is it?");
    }
    $self->{'input'} = $input;
  }
  return $self->{'input'};
}



=head2 run

  Arg [1]   : ReseqTrack::Tools::RunAlignment
  Function  : each child object should implement a run method
  Returntype: n/a
  Exceptions: throws as this method should be implement in the child class
  Example   : 

=cut


sub run{
  my ($self) = @_;
  throw($self." must implement a run method as ReseqTrack::Tools::RunAlignment ".
        "does not provide one");
}


=head2 create_bam_from_sam

  Arg [1]   : ReseqTrack::Tools::RunAlignment
  Arg [2]   : string, path to sam file
  Arg [3]   : binary flag if set sam files marked for deletion, this is on by
  default but it can be turned off if delete_sam is set to 0
  Function  : create a bam file from a sam file
  Returntype: string, path to bam file
  Exceptions: n/a
  Example   : 

=cut



sub create_bam_from_sam{
  my ($self, $sam, $delete_sam) = @_;
  unless(defined($delete_sam)){
    $delete_sam = 1;
  }
  my $bam = $sam;
  $bam =~ s/sam/bam/;
  my $cmd = $self->samtools." import ".$self->reference." ".$sam." ".$bam;
  print $cmd."\n";
  eval{
    my $exit = system($cmd);
    if($exit && $exit >= 1){
      throw("Failed to run ".$cmd);
    }
  };
  if($@){
    throw("Failed to run samtools import $@");
  }
  if($delete_sam){
    $self->files_to_delete($sam);
  }
  return $bam;
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
  foreach my $file(@$files){
    print "Deleting ".$file."\n";
    unlink $file;
  }
}


=head2 assign_fastq_files

  Arg [1]   : ReseqTrack::Tools::RunAlignment
  Arg [2]   : arrayref of filepaths
  Function  : assign file names into mate1, mate2 and frag
  Returntype: array of filepaths
  Exceptions: if the path doesn't match any of the regexs
  Example   : my ($mate1, $mate2, $frag) = $self->assign_files

=cut


sub assign_fastq_files{
  my ($self, $files) = @_;
  return assign_files($files);
}


1;

