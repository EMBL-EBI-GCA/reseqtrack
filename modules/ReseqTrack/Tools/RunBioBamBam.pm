
=pod

=head1 NAME

ReseqTrack::Tools::RunBioBamBam

=head1 SYNOPSIS

This is a class for running biobambam
It is a sub class of a ReseqTrack::Tools::RunProgram.

=head1 Example


=cut

package ReseqTrack::Tools::RunBioBamBam;
use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::FileSystemUtils qw(check_executable check_file_exists);
use File::Basename qw(fileparse);
use File::Copy qw(move);

use base qw(ReseqTrack::Tools::RunProgram);

=head2 DEFAULT_OPTIONS

  Function  : Called by the RunProgram parent object in constructor
  Returntype: hashref
  Example   : my %options = %{&ReseqTrack::Tools:RunBioBamBam::DEFAULT_OPTIONS};

=cut


sub DEFAULT_OPTIONS {
    return {
    };
}

=head2 new

  Arg [-biobambam_dir]   :
      string, directory containing biobambam executables
  Arg [-command]   :
      string, functionality of biobambam, can be "bammarkduplicates", "bammarkduplicates2", "bammerge", "bamsort", "bamtofastq"      
  Arg [-create_index]   :
      boolean, flag to create index files for all outputs
  Arg [-keep_metrics]    :
      boolean, flag to keep metrics files when they aren't the primary product (e.g. mark duplicates )
  + Arguments for ReseqTrack::Tools::RunProgram parent class

  Function  : Creates a new ReseqTrack::Tools::RunBioBamBam object.
  Returntype: ReseqTrack::Tools::RunBioBamBam
  Exceptions: 
  Example   : my $run_biobambam = ReseqTrack::Tools::RunBioBamBam->new(
                -input_files		=> ['/path/sam1', '/path/sam2'],
                -biobambam_dir		=> '/path/to/biobambam/',
                -command			=> 'bammarkduplicates',
                -keep_metrics		=> 1,
                -create_index		=> 1,
                -remove_duplicates	=> 0,
                -working_dir		=> '/path/to/dir/',
                -options			=> {"verbose" => 1, "flag2" => 0} );
=cut

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my ( $biobambam_dir, $create_index, $remove_duplicates, $keep_metrics, $command, $sort_order ) =
    rearrange(
    [
      qw( BIOBAMBAM_DIR CREATE_INDEX REMOVE_DUPLICATES KEEP_METRICS COMMAND SORT_ORDER
        )
    ],
    @args
    );
	
  $self->biobambam_dir( $biobambam_dir || $self->program || $ENV{BIOBAMBAM} );
  $self->keep_metrics( $keep_metrics || 0);
  $self->create_index($create_index || 0);
  $self->remove_duplicates($remove_duplicates || 0);
  $self->sort_order($sort_order || "coordinate");
  $self->command($command);

  return $self;
}

=head2 run_program

  Arg [1]   : ReseqTrack::Tools::RunBioBamBam
  Arg [2]   : string, command, must be one of the following:
              'bammarkduplicates', 'bammarkduplicates2', 'bammerge', 'bamsort', 'bamtofastq'
  Function  : uses biobambam to process the files in $self->input_files.
              Output is files are stored in $self->output_files
              This method is called by the RunProgram run method
  Returntype: 
  Exceptions: Throws if it cannot find the biobambam_dir. Throws if command is not recognised.
  Example   : $self->run();
  
=cut

sub run_program {
  my ( $self) = @_;

  throw("need a biobambam_dir") if !$self->biobambam_dir;
  throw("biobambam_dir is not a directory") if ( !-d $self->biobambam_dir );
  
  my $valid_cmds = valid_commands();
  throw("Did not recognise command " . $self->command )
    if ( !defined $valid_cmds->{$self->command} );
    
  my $bam_name_length = 0;
  foreach my $bam_name ( @{$self->input_files} ) {
    $bam_name_length += length($bam_name);
  }
  if ($bam_name_length > 120000) {
      throw("Please use full path for input files") if ($self->input_files->[0] !~ /\//);
      $self->options->{'shorten_input_names'}=1;
  }  

  my @returned_values;
  if ( $self->command eq "bammarkduplicates" || $self->command eq "bammarkduplicates2" ) {
	  @returned_values = $self->run_bam_mark_duplicates;
  }  
  elsif ( $self->command eq "bammerge") {
      $self->run_bam_merge;
  }   
  elsif ( $self->command eq "bamsort") {
      $self->run_bam_sort;
  }   
  else {
      print "Sorry, the BioBamBam functionality $self->command has not been implemented in the wrapper\n";
  }      
  return @returned_values;
}


=head2 run_sort

  Arg [1]   : ReseqTrack::Tools::RunBioBamBam
  Arg [2]   : String, sort order, defaults to coordinate
  Function  : uses bin/bamsort to sort
  Returntype: 
  Exceptions: 
  Example   : $self->run_bam_sort

=cut

sub run_bam_sort {
  my ($self) = @_;

  foreach my $input ( @{ $self->input_files } ) {
    my $name = fileparse( $input, qr/\.[sb]am/ );
    my $prefix = $self->working_dir . '/' . $name . '.sorted';
    $prefix =~ s{//}{/}g;
    my $output = $prefix . '.bam';
	my $output_index = $output . ".bai";
	  
    my @cmd_words = ( $self->biobambam_dir . "/" . $self->command );
  
    push( @cmd_words, 'SO=' . $self->sort_order );
    push( @cmd_words, 'I=' . $input );
    push( @cmd_words, 'O=' . $output );
    
    push( @cmd_words,
      'index=' . ( $self->create_index ? 1 : 0 ) );
    push( @cmd_words,
      'indexfilename=' . ( $self->create_index ? $output_index : "") );
        
	push( @cmd_words, 'reference=' . $self->reference) if ( $self->options->{'inputformat'} && 
															$self->options->{'inputformat'} =="cram" &&
															$self->options->{'outputformat'} && 
															$self->options->{'outputformat'} =="cram");
	 
    if ( $self->options ) {
      OPTION:
      while (my ($tag, $value) = each %{$self->options}) {
        next OPTION if !defined $value;    
        push(@cmd_words, "$tag=". $value);
      }
    }
	
    my $cmd = join( ' ', @cmd_words );

    $self->output_files($output);
    $self->created_files($output_index) if $self->create_index;

    $self->execute_command_line($cmd);

    if ( $self->create_index ) {
      #my $bai       = "$prefix.bai";
      my $index_ext = $self->options('index_ext');
      if ( $index_ext && $index_ext ne '.bai' ) {
        my $corrected_index = "$output" . $index_ext;
        move( $output_index, $corrected_index )
          or throw "move failed: $!";
        $output_index = $corrected_index;
      }
      $self->output_files($output_index);
    }

  }

  return;
}


=head2 run_bam_merge

  Arg [1]   : ReseqTrack::Tools::RunBioBamBam
  Function  : uses bammerge to merge bam files
  Returntype: 
  Exceptions: 
  Example   : $self->run_bam_merge

=cut

sub run_bam_merge {
  my ($self) = @_;

  my $prefix = $self->working_dir . '/' . $self->job_name . ".merge";
  $prefix =~ s{//}{/}g;
  my $bam = "$prefix.bam";
  my $bai = $bam . ".bai";

  my @cmd_words = ( $self->biobambam_dir . "/" . $self->command );

  my $input_files = $self->options('shorten_input_names') ? [values %{$self->get_short_input_names}]
                    : $self->input_files;
                                        
  push( @cmd_words, map { "I=$_" } @$input_files );
  
  push( @cmd_words,
    'index=' . ( $self->create_index ? 1 : 0 ) );
    
  push ( @cmd_words,
        'indexfilename=' . ( $self->create_index ? $bai : "" ));       
    
  push( @cmd_words, 'SO=' . $self->sort_order );
 
 
 if ( $self->options ) {
      OPTION:
      while (my ($tag, $value) = each %{$self->options}) {
        next OPTION if !defined $value;   
        push(@cmd_words, "$tag=". $value);
      }
  }
 
  push( @cmd_words, ' > ' . $bam );

  my $cmd = join( ' ', @cmd_words );

  $self->output_files($bam);
  $self->created_files($bai) if $self->create_index;

  $self->execute_command_line($cmd);

  if ( $self->create_index ) {
    my $index_ext = $self->options('index_ext');
    if ( $index_ext && $index_ext ne '.bai' ) {
      my $corrected_bai = "$bam" . $index_ext;
      move( $bai, $corrected_bai )
        or throw "move failed: $!";
      $bai = $corrected_bai;
    }
    $self->output_files($bai);
  }


  return;
}

=head2 run_bam_mark_duplicates

  Arg [1]   : ReseqTrack::Tools::RunBioBamBam
  Function  : uses bammarkduplicates or bammarkduplicates2 to mark duplicates; bammarkduplicates works on sorted or not sorted BAMs; bammarkduplicates2 only 
works on sorted BAMs 
  Returntype: 
  Exceptions: 
  Example   : $self->run_bam_mark_duplicates

=cut

sub run_bam_mark_duplicates {
  my ($self) = @_;

  my @metrics_data;
  my $suffix = $self->remove_duplicates ? 'rmdup' : 'mrkdup';
  my $prefix = $self->working_dir . '/' . $self->job_name . ".$suffix";
  $prefix =~ s{//}{/}g;
  my $bam = $prefix . '.bam';
  my $metrics = $prefix . '.dup_metrics' if ($self->keep_metrics);
  my $removed_dup_file = $prefix . '.removed_dups.bam'  if ( $self->remove_duplicates );

  my @cmd_words = ( $self->biobambam_dir . "/" . $self->command );

  my $input_files = $self->options('shorten_input_names') ? [values %{$self->get_short_input_names}]
                    : $self->input_files;
                                        
  push( @cmd_words, map { "I=$_" } @$input_files );
  push( @cmd_words, 'O=' . $bam );
  push( @cmd_words, 'M=' . $metrics ) if ($self->keep_metrics);
  
  push( @cmd_words,
        'rmdup=' . ( $self->remove_duplicates ? 1 : 0  ));
  push( @cmd_words,
        'D=' . $removed_dup_file ) if ( $self->remove_duplicates ); 
  push ( @cmd_words,
        'index=' . ( $self->create_index ? 1 : 0 ));
      
  if ( $self->options ) {
      OPTION:
      while (my ($tag, $value) = each %{$self->options}) {
        next OPTION if !defined $value;    
        push(@cmd_words, "$tag=". $value);
      }
  }

  my $cmd = join( ' ', @cmd_words );
  $self->output_files($bam);
  if ( $self->keep_metrics ) {
    $self->output_files($metrics);
  }
  else {
    $self->created_files($metrics);
  }

  $self->execute_command_line($cmd);
  push @metrics_data, $self->parse_metrics_file($metrics) if ($self->keep_metrics);

  if ( $self->create_index ) {
      my $bai       = "$bam.bai";
      my $index_ext = $self->options('index_ext');
      if ( $index_ext && $index_ext ne '.bai' ) {
        my $corrected_bai = "$bam" . $index_ext;
        move( $bai, $corrected_bai )
          or throw "move failed: $!";
        $bai = $corrected_bai;
      }
      $self->output_files($bai);
  }
  return ( \@metrics_data );
}

sub parse_metrics_file {
  my ( $self, $metrics_file ) = @_;

  open( my $fh, '<', $metrics_file )
    or throw("Could not open metrics file: $metrics_file");

  my @column_headers;
  my @rows;

  while (<$fh>) {
    chomp;
    last if m/## HISTOGRAM/;    # not a good fit for the statistics module
    next if m/^#/;              # comment lines
    next if ( !$_ );            # blank lines

    my @vals = split /\t/;

    if ( scalar(@column_headers) < 1 ) {    # no headers read yet
      @column_headers = @vals;
    }
    else {
      my %row_val;

      for ( my $i = 0 ; $i < scalar(@column_headers) ; $i++ ) {
        undef( $vals[$i] ) if ( defined $vals[$i] && $vals[$i] eq '' );
        $row_val{ $column_headers[$i] } = $vals[$i];
      }

      push @rows, \%row_val;
    }
  }
  return @rows;
}

sub valid_commands { return {
    'bammarkduplicates'   	=> 1,
    'bammarkduplicates2'  	=> 1,
    'bamtofastq'			=> 1,			
    'bamsort'				=> 1,
    'bammerge'				=> 1,
    };    
}

sub biobambam_dir {
  my $self = shift;
  return $self->program(@_);
}

sub command {
  my ($self, $command) = @_;
  if ($command) {
        $self->{'command'} = $command;
  }
  return $self->{'command'};
}

sub sort_order {
  my ($self, $sort_order) = @_;
  if ($sort_order) {
        $self->{'sort_order'} = $sort_order;
  }
  return $self->{'sort_order'};
}

sub keep_metrics {
  my $self = shift;
  if (@_) {
    $self->{'keep_metrics'} = (shift) ? 1 : 0;
  }
  return $self->{'keep_metrics'};
}

sub create_index {
  my $self = shift;
  if (@_) {
    $self->{'create_index'} = (shift) ? 1 : 0;
  }
  return $self->{'create_index'};
}

sub remove_duplicates {
  my $self = shift;
  if (@_) {
    $self->{'remove_duplicates'} = (shift) ? 1 : 0;
  }
  return $self->{'remove_duplicates'};
}

sub output_bai_files {
  my $self = shift;
  my @files = grep { /\.bai$/ } @{ $self->output_files };
  return \@files;
}

sub output_bam_files {
  my $self = shift;
  my @files = grep { /\.bam$/ } @{ $self->output_files };
  return \@files;
}

sub output_metrics_files {
  my $self = shift;
  my @files = grep { /metrics$/ } @{ $self->output_files };
  return \@files;
}

1;

