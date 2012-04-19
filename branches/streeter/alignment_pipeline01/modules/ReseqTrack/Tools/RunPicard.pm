=pod

=head1 NAME

ReseqTrack::Tools::RunPicard

=head1 SYNOPSIS

This is a class for running picard
It is a sub class of a ReseqTrack::Tools::RunProgram.

=head1 Example


=cut

package ReseqTrack::Tools::RunPicard;
use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::FileSystemUtils qw(check_file_does_not_exist make_directory);
use File::Basename qw(fileparse);
use List::Util qw (first);
use Env qw( @PATH );
use File::Copy qw (copy);

use base qw(ReseqTrack::Tools::RunProgram);

=head2 new

  Arg [-picard_dir]   :
      string, directory containing picard jar files
  Arg [-program]   :
      string, the java executable
  Arg [-jvm_options]   :
      string, options for java, default is '-Xmx2g'
  Arg [-options]   :
      string, options for the picard jar file
  Arg [-tmp_dir]   :
      string, path of tmp_dir for picard, default is working_dir
  Arg [-index_bai]   :
      boolean, flag to produce an index file with suffix .bai
  Arg [-index_bam_bai]   :
      boolean, flag to produce an index file with suffix .bam.bai
  + Arguments for ReseqTrack::Tools::RunProgram parent class

  Function  : Creates a new ReseqTrack::Tools::RunPicard object.
  Returntype: ReseqTrack::Tools::RunPicard
  Exceptions: 
  Example   : my $run_picard = ReseqTrack::Tools::RunPicard->new(
                -input_files => ['/path/sam1', '/path/sam2'],
                -picard_dir => '/path/to/picard/',
                -working_dir => '/path/to/dir/',
                -index_bam_bai => 1);

=cut


sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my ( $java_exe, $picard_dir, $jvm_options,
        $index_bai, $index_bam_bai,
        $options,
        )
    = rearrange( [
         qw( JAVA_EXE PICARD_DIR JVM_OPTIONS
                INDEX_BAI INDEX_BAM_BAI
                OPTIONS)
        ], @args);

  if (!$self->program && !$java_exe) {
    $java_exe = first {-x $_} map {"$_/java"} @PATH;
  }
  $self->java_exe($java_exe);

  $self->picard_dir($picard_dir || $ENV{PICARD});
  $self->jvm_options( defined $jvm_options ? $jvm_options : '-Xmx2g' );

  $self->options($options);

  if ($index_bai && $index_bam_bai) {
      throw( "use -index_bai OR -index_bam_bai, not both" );
  }
  $self->flags('index_bam_bai', $index_bam_bai);
  $self->flags('index_bai', $index_bai);

  return $self;
}


=head2 run

  Arg [1]   : ReseqTrack::Tools::RunPicard
  Arg [2]   : string, either 'remove_duplicates', 'merge', 'sort' or 'mark_duplicates' 
  Function  : uses picard to process the files in $self->input_files.
  Output is files are stored in $self->output_files
  Returntype: 
  Exceptions: 
  Example   : $self->run('merge');

=cut

sub run_program{
    my ($self, $command) = @_;

    throw("need a picard_dir") if ! $self->picard_dir;
    throw("Please specify a command") if (! $command);

    my $tmp_dir = $self->working_dir()
          .'/'.$self->job_name.'.'.$$.'.tmp/';
    check_file_does_not_exist($tmp_dir);
    make_directory($tmp_dir);
    $self->created_files($tmp_dir);
    $self->tmp_dir($tmp_dir);

    if ($command eq 'remove_duplicates') {
        $self->run_remove_duplicates;
    }
    elsif ($command eq 'mark_duplicates') {
        $self->run_mark_duplicates;
    }
    elsif ($command eq 'merge') {
        $self->run_merge;
    }
    elsif ($command eq 'sort') {
        $self->run_sort;
    }
    elsif ($command eq 'alignment_metrics'){
        $self->run_alignment_metrics;
    }
   else {
        throw("Did not recognise command $command");
    }
}


=head2 run_remove_duplicates

  Arg [1]   : ReseqTrack::Tools::RunPicard
  Function  : uses MarkDuplicates.jar to remove duplicates
  Returntype: 
  Exceptions: 
  Example   : $self->run_remove_duplicates

=cut

sub run_remove_duplicates{
    my ($self) = @_;
    return $self->_process_duplicates(1);

}

=head2 run_mark_duplicates

  Arg [1]   : ReseqTrack::Tools::RunPicard
  Function  : uses MarkDuplicates.jar to mark duplicates
  Returntype: 
  Exceptions: 
  Example   : $self->run_mark_duplicates

=cut
sub run_mark_duplicates {
    my ($self) = @_;
    return $self->_process_duplicates(0);
}

sub _process_duplicates {
    my ($self,$remove_duplicates) = @_;

    foreach my $input (@{$self->input_files}) {
        my $suffix;
        if ($remove_duplicates) {
            $suffix = '.rmdup';
        }
        else {
            $suffix = '.mrkdup';
        }

        my $name = fileparse($input, qr/\.[sb]am/);
        my $prefix = $self->working_dir . '/' . $name . $suffix;
        $prefix =~ s{//}{/}g;

        my $bam = $prefix . '.bam';
        my $metrics = $prefix . '.metrics';

        my $jar = $self->_jar_path('MarkDuplicates.jar');

        my $cmd = $self->java_exe;
        $cmd .= ' ' . $self->jvm_options if ($self->jvm_options);
        $cmd .= ' -jar ' . $jar;
        $cmd .= ' ' . $self->options if ($self->options);
        $cmd .= ' TMP_DIR=' . $self->tmp_dir;
        $cmd .= ' INPUT=' . $input;
        $cmd .= ' OUTPUT=' . $bam;
        $cmd .= ' METRICS_FILE=' . $metrics;
        if ($remove_duplicates) {
            $cmd .= ' REMOVE_DUPLICATES=true';
        }
        else {
            $cmd .= ' REMOVE_DUPLICATES=false';
        }

        if ($self->flags('index_bai') || $self->flags('index_bam_bai')) {
            $cmd .= ' CREATE_INDEX=true';
            $self->created_files("$prefix.bai");
        }

        $self->output_files($bam);
        $self->created_files($metrics);

        $self->execute_command_line($cmd);


        if ($self->flags('index_bai') || $self->flags('index_bam_bai')) {
            my $bai = "$prefix.bai";
            if ($self->flags('index_bam_bai')) {
                my $bam_bai = "$prefix.bam.bai";
                $self->execute_command_line("mv $bai $bam_bai");
                $bai = $bam_bai;
            }
            $self->output_files($bai);
        }

    }
    return;
}

=head2 run_alignment_metrics

  Arg [1]   : ReseqTrack::Tools::RunPicard
  Function  : uses CollectAlignmentSummaryMetrics.jar to generate alignment metrics file. Reads metrics. 
  Returntype: Collection of hashrefs. Keys described at http://picard.sourceforge.net/picard-metric-definitions.shtml#AlignmentSummaryMetrics
  Exceptions: 
  Example   : my $statistics = $self->run_alignment_metrics
  
=cut
sub run_alignment_metrics {
    my ($self) = @_;

    my ($input) = @{$self->input_files};
    my ($name,$dir,$suffix) = fileparse($input, qr/\.[sb]am/);
    my $base_name = $dir.'/'. $name;
    
    my $output = $base_name.'.metrics';
    my $jar = $self->_jar_path('CollectAlignmentSummaryMetrics.jar');

    my $cmd = $self->java_exe;
    $cmd .= ' ' . $self->jvm_options if ($self->jvm_options);
    $cmd .= ' -jar ' . $jar;
    $cmd .= ' ' . $self->options if ($self->options);
    $cmd .= ' INPUT=' . $input;
    $cmd .= ' OUTPUT=' . $output;

    $self->output_files($output);
    $self->execute_command_line($cmd);


    open(my $fh, "<", $output) or die "cannot open < $output: $!";

    my (@titles, @statistics);
    
    while (my $line = <$fh>) {
        next if $line =~ m/#/;
        chomp $line;
        my @vals = split /\s/, $line;
        next unless @vals;
        
        if (@titles){
            my %metrics;
            push @statistics, \%metrics;
            for my $i (0..scalar(@titles)){                
                my ($key, $value) = ($titles[$i],$vals[$i]);
                next unless ($key && $value);                            
                $metrics{$key} = $value;
            }
        }
        else {
            @titles = @vals;
        }
    }        
    close $fh;
    
    return \@statistics;
}

=head2 run_merge

  Arg [1]   : ReseqTrack::Tools::RunPicard
  Function  : uses MergeSamFiles.jar to merge, sort and index
  Returntype: 
  Exceptions: 
  Example   : $self->run_merge

=cut


sub run_merge{
    my ($self) = @_;

    my $prefix = $self->working_dir . '/' . $self->job_name;
    $prefix =~ s{//}{/}g;

    my $input_bam_list = $self->input_files;
    my $bam = "$prefix.merge.bam";
    $self->output_files($bam);

    my $jar = $self->_jar_path('MergeSamFiles.jar');

    if (@$input_bam_list >1) {
      my $cmd = $self->java_exe;
      $cmd .= ' ' . $self->jvm_options if ($self->jvm_options);
      $cmd .= ' -jar ' . $jar;
      $cmd .= ' ' . $self->options if ($self->options);
      $cmd .= ' TMP_DIR=' . $self->tmp_dir;
      foreach my $input (@$input_bam_list) {
          $cmd .= ' INPUT=' . $input;
      }
      $cmd .= ' OUTPUT=' . $bam;

      if ($self->flags('index_bai') || $self->flags('index_bam_bai')) {
          $cmd .= ' CREATE_INDEX=true';
          $self->created_files("$prefix.bai");
      }

      $self->execute_command_line($cmd);
    }
    elsif(@$input_bam_list ==1) {
      print "copying ".$input_bam_list->[0]." to $bam\n";
      copy($input_bam_list->[0], $bam)
          or throw "copy failed: $!";
    }

    if ($self->flags('index_bai') || $self->flags('index_bam_bai')) {
        my $bai = "$prefix.bai";
        if ($self->flags('index_bam_bai')) {
            my $bam_bai = "$prefix.bam.bai";
            $self->execute_command_line("mv $bai $bam_bai");
            $bai = $bam_bai;
        }
        $self->output_files($bai);
    }


    return;
}

=head2 run_sort

  Arg [1]   : ReseqTrack::Tools::RunPicard
  Arg [2]   : String, sort order, defaults to coordinate
  Function  : uses SortSam.jar to sort
  Returntype: 
  Exceptions: 
  Example   : $self->run_sort

=cut

sub run_sort{
  my ($self,$sort_order) = @_;
    
  $sort_order ||= 'coordinate';
          
  my $jar = $self->_jar_path('SortSam.jar');

  foreach my $input (@{$self->input_files}) {
    my $name = fileparse($input, qr/\.[sb]am/);
    my $prefix = $self->working_dir . '/' . $name . '.sorted';
    $prefix =~ s{//}{/}g;
    my $output = $prefix . '.bam';
            
    my $cmd = $self->java_exe;
    $cmd .= ' ' . $self->jvm_options if ($self->jvm_options);
    $cmd .= ' -jar ' . $jar;
    $cmd .= ' ' . $self->options if ($self->options);
    $cmd .= ' SORT_ORDER='. $sort_order;
    $cmd .= ' INPUT=' . $input;
    $cmd .= ' OUTPUT=' . $output;
  
    $self->output_files($output);
    $self->execute_command_line($cmd);
  
  }

    return;   
}

sub _jar_path{
    my ($self,$jar_file) = @_;
    my $jar = $self->picard_dir . '/' . $jar_file;
    $jar =~ s{//}{/}g;

    return $jar;
}


sub picard_dir {
    my ($self, $arg) = @_;
    if ($arg) {
        $self->{'picard_dir'} = $arg;
    }
    return $self->{'picard_dir'};
}

sub jvm_options {
    my ($self, $arg) = @_;
    if ($arg) {
        $self->{'jvm_options'} = $arg;
    }
    return $self->{'jvm_options'};
}

sub options {
    my ($self, $arg) = @_;
    if ($arg) {
        $self->{'options'} = $arg;
    }
    return $self->{'options'};
}

sub tmp_dir {
    my ($self, $arg) = @_;
    if ($arg) {
        $self->{'tmp_dir'} = $arg;
    }
    return $self->{'tmp_dir'};
}

sub java_exe {
  my $self = shift;
  return $self->program(@_);
}




=head2 flags

  Arg [1]   : ReseqTrack::Tools::RunPicard
  Arg [2]   : string, name of key e.g. "index_bai" or "index_bam_bai"
  Arg [3]   : boolean, optional, turn flag on or off
  Function  : accessor method for flags to use the various samtools utilities.
  Returntype: boolean, flag
  Exceptions: n/a
  Example   : my $run_sort_flag = $self->flags{'index_bai'};

=cut

sub flags {
    my ($self, $flag_name, $flag_value) = @_;

    if (! $self->{'flags'}) {
        $self->{'flags'} = {};
    }

    throw( "flag_name not specified")
        if (! $flag_name);

    if ($flag_value) {
        $self->{'flags'}->{$flag_name} = $flag_value ? 1 : 0;
    }

    return $self->{'flags'}->{$flag_name};
}


sub output_bai_files {
  my $self = shift;
  my @files = grep {/\.bai$/} @{$self->output_files};
  return \@files;
}

sub output_bam_files {
  my $self = shift;
  my @files = grep {/\.bam$/} @{$self->output_files};
  return \@files;
}


1;

