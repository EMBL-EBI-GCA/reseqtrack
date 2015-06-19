
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
use ReseqTrack::Tools::FileSystemUtils qw(check_executable check_file_exists);
use File::Basename qw(fileparse);
use File::Copy qw(move);

use base qw(ReseqTrack::Tools::RunProgram);

=head2 DEFAULT_OPTIONS

  Function  : Called by the RunProgram parent object in constructor
  Returntype: hashref
  Example   : my %options = %{&ReseqTrack::Tools:RunPicard::DEFAULT_OPTIONS};

=cut

sub DEFAULT_OPTIONS {
    return {
        'assume_sorted' =>
          undef,  # assume input files are sorted, even if header says otherwise
        'sort_order' => 'coordinate',    # can be 'coordinate' or 'queryname'
        'index_ext'  => '.bam.bai',      #can be '.bai' or '.bam.bai'
        'verbosity'  => 'INFO',
        'quiet'      => 'false',
        'max_records_in_ram'    => 500000,
        'remove_duplicates'     => 0,        # used by run_mark_duplicates
        'use_threading'         => 0,        # used by merge
        'validation_stringency' => undef,
        'max_file_handles' => 1000, # should be slightly less than 'ulimit -n'
        'shorten_input_names' => 0, # should be used by merge and mark_duplicates when the number of input files is very large
        'ref_flat' =>
          undef
        , # gene annotations file in ref flat format, used by CollectRnaSeqMetrics
        'ribosomal_intervals' =>
          undef
        , # Location of rRNA sequences in genome, in interval_list format, used by CollectRnaSeqMetrics
        'reference_sequence' =>
          undef
        ,    # Location of reference fasta file, used by CollectRnaSeqMetrics and ReorderSam
        'strand_specificity' => 'NONE'
        ,    #For strand-specific library prep,  used by CollectRnaSeqMetrics
        'read_group_fields' => {}, # used by AddOrReplaceReadGroups
        'keep_mapped_regex' => 'concordant_uniq|concordant_mult',
        'metrics_programs'  => [
          'CollectAlignmentSummaryMetrics', 'CollectInsertSizeMetrics',
          'QualityScoreDistribution',       'MeanQualityByCycle',
          ],
    };
}


sub CMD_MAPPINGS { return {
    'mark_duplicates'   => \&run_mark_duplicates,
    'fix_mate'          => \&run_fix_mate,
    'merge'             => \&run_merge,
    'sort'              => \&run_sort,
    'alignment_metrics' => \&run_alignment_metrics,
    'rna_seq_metrics'   => \&run_rna_alignment_metrics,
    'add_or_replace_read_groups'   => \&run_add_or_replace_read_groups,
    'demap_merge'         => \&run_demap_merge,
    'insert_size_metrics' => \&run_insert_size_metrics,
    'multiple_metrics'    => \&run_multiple_metrics,
    'reorder'             => \&run_reorder,
    };    
}

=head2 new

  Arg [-picard_dir]   :
      string, directory containing picard jar files
  Arg [-java_exe]   :
      string, the java executable, default is 'java'
  Arg [-jvm_options]   :
      string, options for java, default is '-Xmx4g'
  Arg [-create_index]   :
      boolean, flag to create index files for all outputs
  Arg [-keep_metrics]    :
      boolean, flag to keep metrics files when they aren't the primary product (e.g. mark duplicates )
  + Arguments for ReseqTrack::Tools::RunProgram parent class

  Function  : Creates a new ReseqTrack::Tools::RunPicard object.
  Returntype: ReseqTrack::Tools::RunPicard
  Exceptions: 
  Example   : my $run_picard = ReseqTrack::Tools::RunPicard->new(
                -input_files => ['/path/sam1', '/path/sam2'],
                -picard_dir => '/path/to/picard/',
                -working_dir => '/path/to/dir/',
                -create_index => 1);

=cut

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my ( $java_exe, $picard_dir, $jvm_options, $create_index, $keep_metrics, ) =
    rearrange(
    [
      qw( JAVA_EXE PICARD_DIR JVM_OPTIONS CREATE_INDEX KEEP_METRICS
        )
    ],
    @args
    );

  $self->java_exe( $java_exe         || 'java' );
  $self->picard_dir( $picard_dir     || $self->program || $ENV{PICARD} );
  $self->keep_metrics( $keep_metrics || 0 );

  $self->jvm_options( defined $jvm_options ? $jvm_options : '-Xmx4g' );
  $self->create_index($create_index);

  return $self;
}

=head2 run_program

  Arg [1]   : ReseqTrack::Tools::RunPicard
  Arg [2]   : string, command, must be one of the following:
              'mark_duplicates', 'merge', 'sort', 'fix_mate' or 'alignment_metrics' 
  Function  : uses picard to process the files in $self->input_files.
              Output is files are stored in $self->output_files
              This method is called by the RunProgram run method
  Returntype: 
  Exceptions: Throws if it cannot find the picard_dir. Throws if command is not recognised.
  Example   : $self->run('merge');

=cut

sub run_program {
  my ( $self, $command ) = @_;

  throw("need a picard_dir") if !$self->picard_dir;
  throw("picard_dir is not a directory") if ( !-d $self->picard_dir );
  check_executable( $self->java_exe );

  my $subs = CMD_MAPPINGS();

  throw("Did not recognise command $command")
    if ( !defined $subs->{$command} );

  my @returned_values = &{ $subs->{$command} }($self);

  return @returned_values;
}

sub get_valid_commands {
  my $subs = CMD_MAPPINGS();
  return keys %$subs;
}

=head2 run_mark_duplicates

  Arg [1]   : ReseqTrack::Tools::RunPicard
  Function  : uses MarkDuplicates.jar to mark duplicates
  Returntype: 
  Exceptions: 
  Example   : $self->run_mark_duplicates

=cut

sub run_mark_duplicates {
  my $self = shift;

    my $jar = $self->_jar_path('MarkDuplicates.jar');
    my @metrics_data;
    my $suffix = $self->options('remove_duplicates') ? 'rmdup' : 'mrkdup';
    my $prefix = $self->working_dir . '/' . $self->job_name . ".$suffix";
    $prefix =~ s{//}{/}g;
    my $bam     = $prefix . '.bam';
    my $metrics = $prefix . '.dup_metrics';

    my @cmd_words = ( $self->java_exe );
    push( @cmd_words, $self->jvm_options ) if ( $self->jvm_options );
    push( @cmd_words, '-jar', $jar );
    push( @cmd_words, $self->_get_standard_options );
    my $input_files = $self->options('shorten_input_names') ? [values %{$self->get_short_input_names}]
                    : $self->input_files;
    push( @cmd_words, map { "INPUT=$_" } @$input_files );
    push( @cmd_words, 'OUTPUT=' . $bam );
    push( @cmd_words, 'METRICS_FILE=' . $metrics );
    push( @cmd_words,
        'REMOVE_DUPLICATES='
          . ( $self->options('remove_duplicates') ? 'true' : 'false' ) );
    push( @cmd_words,
        'ASSUME_SORTED='
          . ( $self->options('assume_sorted') ? 'true' : 'false' ) )
      if defined $self->options('assume_sorted');
    push( @cmd_words,
        'CREATE_INDEX=' . ( $self->create_index ? 'true' : 'false' ) );
    push( @cmd_words, 'MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=' . $self->options('max_file_handles'))
      if $self->options('max_file_handles');
    my $cmd = join( ' ', @cmd_words );
    $self->output_files($bam);
    if ( $self->keep_metrics ) {
      $self->output_files($metrics);
    }
    else {
      $self->created_files($metrics);
    }

    $self->created_files("$prefix.bai") if $self->create_index;
    $self->execute_command_line($cmd);
    push @metrics_data, $self->parse_metrics_file($metrics);

    if ( $self->create_index ) {
      my $bai       = "$prefix.bai";
      my $index_ext = $self->options('index_ext');
      if ( $index_ext && $index_ext ne '.bai' ) {
        my $corrected_bai = "$prefix" . $index_ext;
        move( $bai, $corrected_bai )
          or throw "move failed: $!";
        $bai = $corrected_bai;
      }
      $self->output_files($bai);
    }

  return ( \@metrics_data );
}

=head2 run_alignment_metrics

  Arg [1]   : ReseqTrack::Tools::RunPicard
  Function  : uses CollectAlignmentSummaryMetrics.jar to generate alignment metrics file. Reads metrics. 
  Returntype: Collection of hashrefs. Keys described at http://picard.sourceforge.net/picard-metric-definitions.shtml#AlignmentSummaryMetrics
  Exceptions: 
  
=cut

sub run_alignment_metrics {
  my ($self) = @_;

  my @metrics;

  throw("Reference sequence is required for meaningful metrics")
    unless $self->options('reference_sequence');

  my $jar = $self->_jar_path('CollectAlignmentSummaryMetrics.jar');
  foreach my $input ( @{ $self->input_files } ) {
    my ( $name, $dir ) = fileparse( $input, qr/\.[sb]am/ );
    my $base_name = $dir . '/' . $name;

    my $output = $base_name . '.metrics';

    my @cmd_words = ( $self->java_exe );
    push( @cmd_words, $self->jvm_options ) if ( $self->jvm_options );
    push( @cmd_words, '-jar', $jar );
    push( @cmd_words, $self->_get_standard_options );
    push( @cmd_words, 'INPUT=' . $input );
    push( @cmd_words, 'OUTPUT=' . $output );
    push( @cmd_words,
      'ASSUME_SORTED='
        . ( $self->options('assume_sorted') ? 'true' : 'false' ) )
      if defined $self->options('assume_sorted');
    push( @cmd_words,
      'REFERENCE_SEQUENCE=' . $self->options('reference_sequence') )
      if $self->options('reference_sequence');
    my $cmd = join( ' ', @cmd_words );

    $self->output_files($output);
    $self->execute_command_line($cmd);

    push @metrics, $self->parse_metrics_file($output);
  }
  return ( \@metrics );
}

=head2 run_insert_size_metrics

  Arg [1]   : ReseqTrack::Tools::RunPicard
  Function  : uses CollectInsertSizeMetrics.jar to generate alignment metrics file. Reads metrics. 
  Returntype: Collection of hashrefs. Keys described at http://picard.sourceforge.net/picard-metric-definitions.shtml#AlignmentSummaryMetrics
  Exceptions: 
  
=cut

sub run_insert_size_metrics {
  my ($self) = @_;

  my @metrics;

  throw("Reference sequence is required for meaningful metrics")
    unless $self->options('reference_sequence');

  my $jar = $self->_jar_path('CollectInsertSizeMetrics.jar');
  foreach my $input ( @{ $self->input_files } ) {

    my ( $name, $dir ) = fileparse( $input, qr/\.[sb]am/ );
    my $base_name = $dir . '/' . $name;

    my $output    = $base_name . '.insert_size_metrics';
    my $histogram = $base_name . '.insert_size.pdf';

    my @cmd_words = ( $self->java_exe );
    push( @cmd_words, $self->jvm_options ) if ( $self->jvm_options );
    push( @cmd_words, '-jar', $jar );

    push( @cmd_words, $self->_get_standard_options );

    push( @cmd_words, 'INPUT=' . $input );
    push( @cmd_words, 'OUTPUT=' . $output );
    push( @cmd_words, 'HISTOGRAM_FILE=' . $histogram );
    push( @cmd_words,
      'ASSUME_SORTED='
        . ( $self->options('assume_sorted') ? 'true' : 'false' ) )
      if defined $self->options('assume_sorted');
    push( @cmd_words,
      'REFERENCE_SEQUENCE=' . $self->options('reference_sequence') )
      if $self->options('reference_sequence');
    my $cmd = join( ' ', @cmd_words );

    $self->output_files($output);
    $self->output_files($histogram);
    $self->execute_command_line($cmd);

    push @metrics, $self->parse_metrics_file($output);

  }
  return ( \@metrics );
}

sub run_multiple_metrics {
  my ($self) = @_;

  my @metrics;

  my %file_suffixes = (
    'CollectAlignmentSummaryMetrics' => ['.alignment_summary_metrics'],
    'CollectInsertSizeMetrics' =>
      [ '.insert_size_metrics', '.insert_size_histogram.pdf' ],
    'QualityScoreDistribution' =>
      [ '.quality_distribution_metrics', '.quality_distribution.pdf' ],
    'MeanQualityByCycle' =>
      [ '.quality_by_cycle_metrics', '.quality_by_cycle.pdf' ],
  );

  throw("Reference sequence is required for meaningful metrics")
    unless $self->options('reference_sequence');

  my $jar = $self->_jar_path('CollectMultipleMetrics.jar');
  foreach my $input ( @{ $self->input_files } ) {

    my ( $name, $dir ) = fileparse( $input, qr/\.[sb]am/ );
    my $output = $dir . '/' . $name;

    my @cmd_words = ( $self->java_exe );
    push( @cmd_words, $self->jvm_options ) if ( $self->jvm_options );
    push( @cmd_words, '-jar', $jar );

    push( @cmd_words, $self->_get_standard_options );

    push( @cmd_words, 'INPUT=' . $input );
    push( @cmd_words, 'OUTPUT=' . $output );
    push( @cmd_words,
      'ASSUME_SORTED='
        . ( $self->options('assume_sorted') ? 'true' : 'false' ) )
      if defined $self->options('assume_sorted');
    push( @cmd_words,
      'REFERENCE_SEQUENCE=' . $self->options('reference_sequence') )
      if $self->options('reference_sequence');

    for my $program ( @{ $self->options('metrics_programs') } ) {
      push( @cmd_words, 'PROGRAM=' . $program );
    }
    
    my $cmd = join( ' ', @cmd_words );
    $self->execute_command_line($cmd);

    for my $program ( @{ $self->options('metrics_programs') } ) {
      my $suffixes = $file_suffixes{$program};
      
      for my $suffix (@$suffixes) {
        my $output_file = $output . $suffix;
        $self->output_files( $output_file ) if (-e $output_file);
      }
    }

    for my $metrics_file ( @{ $self->output_metrics_files } ) {
      push @metrics, $self->parse_metrics_file($metrics_file);
    }

  }
  return ( \@metrics );
}

=head2 run_rna_alignment_metrics

  Arg [1]   : ReseqTrack::Tools::RunPicard
  Function  : uses CollectRnaSeqMetrics.jar to generate alignment metrics file. Reads metrics. 
  Returntype: Collection of hashrefs. Keys described at http://picard.sourceforge.net/picard-metric-definitions.shtml#AlignmentSummaryMetrics
  Exceptions: 
  
=cut

sub run_rna_alignment_metrics {
  my ($self) = @_;

  my @metrics;

  throw('A REF_FLAT file must be specified for RNA alignment metrics')
    unless $self->options('ref_flat');

  my $jar = $self->_jar_path('CollectRnaSeqMetrics.jar');
  foreach my $input ( @{ $self->input_files } ) {
    my ($input) = @{ $self->input_files };
    my ( $name, $dir ) = fileparse( $input, qr/\.[sb]am/ );
    my $base_name = $dir . '/' . $name;

    my $output    = $base_name . '.rnaseq_metrics';
    my @cmd_words = ( $self->java_exe );
    push( @cmd_words, $self->jvm_options ) if ( $self->jvm_options );
    push( @cmd_words, '-jar', $jar );
    push( @cmd_words, $self->_get_standard_options );
    push( @cmd_words, 'INPUT=' . $input );
    push( @cmd_words, 'OUTPUT=' . $output );
    push( @cmd_words,
      'ASSUME_SORTED='
        . ( $self->options('assume_sorted') ? 'true' : 'false' ) )
      if defined $self->options('assume_sorted');

    push( @cmd_words, 'REF_FLAT=' . $self->options('ref_flat') );
    push( @cmd_words,
      'RIBOSOMAL_INTERVALS=' . $self->options('ribosomal_intervals') )
      if $self->options('ribosomal_intervals');
    push( @cmd_words,
      'REFERENCE_SEQUENCE=' . $self->options('reference_sequence') )
      if $self->options('reference_sequence');

    push( @cmd_words,
      'STRAND_SPECIFICITY=' . $self->options('strand_specificity') );

    my $cmd = join( ' ', @cmd_words );

    $self->output_files($output);
    $self->execute_command_line($cmd);

    push @metrics, $self->parse_metrics_file($output);
  }
  return ( \@metrics );
}

=head2 run_merge

  Arg [1]   : ReseqTrack::Tools::RunPicard
  Function  : uses MergeSamFiles.jar to merge, sort and index
  Returntype: 
  Exceptions: 
  Example   : $self->run_merge

=cut

sub run_merge {
  my ($self) = @_;

  my $input_bam_list = $self->options('shorten_input_names') ? [values %{$self->get_short_input_names}]
                  : $self->input_files;
  throw("need more than two or more files to merge")
    if ( @$input_bam_list < 2 );

  my $prefix = $self->working_dir . '/' . $self->job_name . ".merge";
  $prefix =~ s{//}{/}g;
  my $bam = "$prefix.bam";

  my $jar = $self->_jar_path('MergeSamFiles.jar');

  my $sort_order = $self->options('sort_order');
  if ( $sort_order ne 'coordinate' && $sort_order ne 'queryname' ) {
    $sort_order = 'null';
  }

  my @cmd_words = ( $self->java_exe );
  push( @cmd_words, $self->jvm_options ) if ( $self->jvm_options );
  push( @cmd_words, '-jar', $jar );
  push( @cmd_words, $self->_get_standard_options );
  push( @cmd_words, map { "INPUT=$_" } @$input_bam_list );
  push( @cmd_words, 'OUTPUT=' . $bam );
  push( @cmd_words,
    'CREATE_INDEX=' . ( $self->create_index ? 'true' : 'false' ) );
  push( @cmd_words, 'SORT_ORDER=' . $sort_order );
  push( @cmd_words,
    'ASSUME_SORTED=' . ( $self->options('assume_sorted') ? 'true' : 'false' ) )
    if defined $self->options('assume_sorted');
  push( @cmd_words,
    'USE_THREADING=' . ( $self->options('use_threading') ? 'true' : 'false' ) );

  my $cmd = join( ' ', @cmd_words );

  $self->output_files($bam);
  $self->created_files("$prefix.bai") if $self->create_index;

  $self->execute_command_line($cmd);

  if ( $self->create_index ) {
    my $bai       = "$prefix.bai";
    my $index_ext = $self->options('index_ext');
    if ( $index_ext && $index_ext ne '.bai' ) {
      my $corrected_bai = "$prefix" . $index_ext;
      move( $bai, $corrected_bai )
        or throw "move failed: $!";
      $bai = $corrected_bai;
    }
    $self->output_files($bai);
  }

  return;
}

sub run_demap_merge {
  my ($self) = @_;

  my $input = $self->input_files;
  throw("need more than two or more files to merge") if ( @$input < 2 );

  my $prefix = $self->working_dir . '/' . $self->job_name . ".merge";
  $prefix =~ s{//}{/}g;
  my $bam = "$prefix.bam";

  my ( @demap, @keep_mapped );
  my $pattern = $self->options('keep_mapped_regex');
  my $re      = qr/$pattern/;
  for my $file (@$input) {

    if ( $file =~ $re ) {
      push @keep_mapped, $file;
    }
    else {
      push @demap, $file;
    }
  }

  my $jar = $self->_jar_path('DeMapSam.jar');

  my @cmd_words = ( $self->java_exe );
  push( @cmd_words, $self->jvm_options ) if ( $self->jvm_options );
  push( @cmd_words, '-jar', $jar );
  push( @cmd_words, $self->_get_standard_options );

  push( @cmd_words, map { "INPUT=$_" } @demap );
  push( @cmd_words, map { "IKM=$_" } @keep_mapped );

  push( @cmd_words, 'OUTPUT=' . $bam );
  push( @cmd_words,
    'USE_THREADING=' . ( $self->options('use_threading') ? 'true' : 'false' ) );

  my $cmd = join( ' ', @cmd_words );

  $self->output_files($bam);

  $self->execute_command_line($cmd);

  return;
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

=head2 run_sort

  Arg [1]   : ReseqTrack::Tools::RunPicard
  Arg [2]   : String, sort order, defaults to coordinate
  Function  : uses SortSam.jar to sort
  Returntype: 
  Exceptions: 
  Example   : $self->run_sort

=cut

sub run_sort {
  my ($self) = @_;

  my $sort_order = $self->options('sort_order');
  if ( $sort_order ne 'coordinate' && $sort_order ne 'queryname' ) {
    $sort_order = 'null';
  }

  my $jar = $self->_jar_path('SortSam.jar');

  foreach my $input ( @{ $self->input_files } ) {
    my $name = fileparse( $input, qr/\.[sb]am/ );
    my $prefix = $self->working_dir . '/' . $name . '.sorted';
    $prefix =~ s{//}{/}g;
    my $output = $prefix . '.bam';

    my @cmd_words = ( $self->java_exe );
    push( @cmd_words, $self->jvm_options ) if ( $self->jvm_options );
    push( @cmd_words, '-jar', $jar );
    push( @cmd_words, $self->_get_standard_options );
    push( @cmd_words, 'SORT_ORDER=' . $sort_order );
    push( @cmd_words, 'INPUT=' . $input );
    push( @cmd_words, 'OUTPUT=' . $output );
    push( @cmd_words,
      'CREATE_INDEX=' . ( $self->create_index ? 'true' : 'false' ) );

    my $cmd = join( ' ', @cmd_words );

    $self->output_files($output);
    $self->created_files("$prefix.bai") if $self->create_index;

    $self->execute_command_line($cmd);

    if ( $self->create_index ) {
      my $bai       = "$prefix.bai";
      my $index_ext = $self->options('index_ext');
      if ( $index_ext && $index_ext ne '.bai' ) {
        my $corrected_bai = "$prefix" . $index_ext;
        move( $bai, $corrected_bai )
          or throw "move failed: $!";
        $bai = $corrected_bai;
      }
      $self->output_files($bai);
    }

  }

  return;
}


=head2 run_reorder

  Arg [1]   : ReseqTrack::Tools::RunPicard
  Function  : uses ReorderSam.jar to reorder reads to match a different contig order
  Returntype: 
  Exceptions: 
  Example   : $self->run_sort

=cut

sub run_reorder {
  my ($self) = @_;

  my $jar = $self->_jar_path('ReorderSam.jar');

  throw("Reference sequence is required to reorder contigs")
    unless $self->options('reference_sequence');
  check_file_exists $self->options('reference_sequence');

  foreach my $input ( @{ $self->input_files } ) {
    my $name = fileparse( $input, qr/\.[sb]am/ );
    my $prefix = $self->working_dir . '/' . $name . '.reordered';
    $prefix =~ s{//}{/}g;
    my $output = $prefix . '.bam';

    my @cmd_words = ( $self->java_exe );
    push( @cmd_words, $self->jvm_options ) if ( $self->jvm_options );
    push( @cmd_words, '-jar', $jar );
    push( @cmd_words, $self->_get_standard_options );
    push( @cmd_words, 'INPUT=' . $input );
    push( @cmd_words, 'OUTPUT=' . $output );
    push( @cmd_words, 'REFERENCE=' . $self->options('reference_sequence'));
    push( @cmd_words,
      'CREATE_INDEX=' . ( $self->create_index ? 'true' : 'false' ) );

    my $cmd = join( ' ', @cmd_words );

    $self->output_files($output);
    $self->created_files("$prefix.bai") if $self->create_index;

    $self->execute_command_line($cmd);

    if ( $self->create_index ) {
      my $bai       = "$prefix.bai";
      my $index_ext = $self->options('index_ext');
      if ( $index_ext && $index_ext ne '.bai' ) {
        my $corrected_bai = "$prefix" . $index_ext;
        move( $bai, $corrected_bai )
          or throw "move failed: $!";
        $bai = $corrected_bai;
      }
      $self->output_files($bai);
    }

  }

  return;
}

=head2 run_fix_mate

  Arg [1]   : ReseqTrack::Tools::RunPicard
  Function  : uses FixMateInformation.jar to ensure mate-pair info is in sync
  Returntype: 
  Exceptions: 
  Example   : $self->run_fix_mate

=cut

sub run_fix_mate {
  my ($self) = @_;

  my $sort_order = $self->options('sort_order');
  if ( $sort_order ne 'coordinate' && $sort_order ne 'queryname' ) {
    $sort_order = 'null';
  }

  my $jar = $self->_jar_path('FixMateInformation.jar');

  foreach my $input ( @{ $self->input_files } ) {
    my $name = fileparse( $input, qr/\.[sb]am/ );
    my $prefix = $self->working_dir . '/' . $name . '.fixed';
    $prefix =~ s{//}{/}g;
    my $output = $prefix . '.bam';

    my @cmd_words = ( $self->java_exe );
    push( @cmd_words, $self->jvm_options ) if ( $self->jvm_options );
    push( @cmd_words, '-jar', $jar );
    push( @cmd_words, $self->_get_standard_options );
    push( @cmd_words, 'SORT_ORDER=' . $sort_order );
    push( @cmd_words, 'INPUT=' . $input );
    push( @cmd_words, 'OUTPUT=' . $output );
    push( @cmd_words,
      'CREATE_INDEX=' . ( $self->create_index ? 'true' : 'false' ) );

    my $cmd = join( ' ', @cmd_words );

    $self->output_files($output);
    $self->created_files("$prefix.bai") if $self->create_index;

    $self->execute_command_line($cmd);

    if ( $self->create_index ) {
      my $bai       = "$prefix.bai";
      my $index_ext = $self->options('index_ext');
      if ( $index_ext && $index_ext ne '.bai' ) {
        my $corrected_bai = "$prefix" . $index_ext;
        move( $bai, $corrected_bai )
          or throw "move failed: $!";
        $bai = $corrected_bai;
      }
      $self->output_files($bai);
    }

  }

  return;
}

sub run_add_or_replace_read_groups{
  my ($self) = @_;
    
  my $sort_order = $self->options('sort_order');
  if ($sort_order ne 'coordinate' && $sort_order ne 'queryname') {
      $sort_order = 'null';
  }
          
  my $jar = $self->_jar_path('AddOrReplaceReadGroups.jar');

  foreach my $input (@{$self->input_files}) {
    my $name = fileparse($input, qr/\.[sb]am/);
    my $prefix = $self->working_dir . '/' . $name . '.fixed';
    $prefix =~ s{//}{/}g;
    my $output = $prefix . '.bam';

    my @cmd_words = ($self->java_exe);
    push(@cmd_words, $self->jvm_options) if ($self->jvm_options);
    push(@cmd_words, '-jar', $jar);
    push(@cmd_words, $self->_get_standard_options);
    push(@cmd_words, 'SORT_ORDER=' . $sort_order);
    push(@cmd_words, 'INPUT=' . $input);
    push(@cmd_words, 'OUTPUT=' . $output);
    push(@cmd_words, 'CREATE_INDEX='
            . ($self->create_index ? 'true' : 'false'));

    my $rg_fields = $self->options('read_group_fields');
    push(@cmd_words, "RGID='" . $rg_fields->{'ID'} . "'");
    push(@cmd_words, "RGLB='" . $rg_fields->{'LB'} . "'");
    push(@cmd_words, "RGPL='" . $rg_fields->{'PL'} . "'");
    push(@cmd_words, "RGPU='" . $rg_fields->{'PU'} . "'");
    push(@cmd_words, "RGSM='" . $rg_fields->{'SM'} . "'");
    push(@cmd_words, "RGCN='" . $rg_fields->{'CN'} . "'") if $rg_fields->{'CN'};
    push(@cmd_words, "RGDS='" . $rg_fields->{'DS'} . "'") if $rg_fields->{'DS'};
    push(@cmd_words, "RGDT='" . $rg_fields->{'DT'} . "'") if $rg_fields->{'DT'};
    push(@cmd_words, "RGPI='" . $rg_fields->{'PI'} . "'") if $rg_fields->{'PI'};

    my $cmd = join(' ', @cmd_words);
  
    $self->output_files($output);
    $self->created_files("$prefix.bai") if $self->create_index;

    $self->execute_command_line($cmd);

    if ($self->create_index) {
        my $bai = "$prefix.bai";
        my $index_ext = $self->options('index_ext');
        if ($index_ext && $index_ext ne '.bai') {
            my $corrected_bai = "$prefix" . $index_ext;
            move($bai, $corrected_bai)
              or throw "move failed: $!";
            $bai = $corrected_bai;
        }
        $self->output_files($bai);
    }
  
  }

    return;   
}


sub _get_standard_options {
  my $self = shift;

  my @option_strings;
  push( @option_strings, 'TMP_DIR=' . $self->get_temp_dir );
  push( @option_strings, 'VERBOSITY=' . $self->options('verbosity') )
    if ( $self->options('verbosity') );
  push( @option_strings,
    'MAX_RECORDS_IN_RAM=' . $self->options('max_records_in_ram') )
    if ( $self->options('max_records_in_ram') );
  push( @option_strings,
    'VALIDATION_STRINGENCY=' . $self->options('validation_stringency') )
    if ( $self->options('validation_stringency') );
  push( @option_strings, 'QUIET=' . $self->options('quiet') )
    if ( $self->options('quiet') );
    push( @option_strings, 'CREATE_MD5_FILE='.$self->options('create_md5_file')) if $self->options('create_md5_file');

  return join( ' ', @option_strings );
}

sub _jar_path {
  my ( $self, $jar_file ) = @_;
  my $jar = $self->picard_dir . '/' . $jar_file;
  $jar =~ s{//}{/}g;

  return $jar;
}

sub jvm_options {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'jvm_options'} = $arg;
  }
  return $self->{'jvm_options'};
}

sub picard_dir {
  my $self = shift;
  return $self->program(@_);
}

sub java_exe {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'java_exe'} = $arg;
  }
  return $self->{'java_exe'};
}

sub create_index {
  my $self = shift;
  if (@_) {
    $self->{'create_index'} = (shift) ? 1 : 0;
  }
  return $self->{'create_index'};
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

sub keep_metrics {
  my $self = shift;
  if (@_) {
    $self->{'keep_metrics'} = (shift) ? 1 : 0;
  }
  return $self->{'keep_metrics'};
}

1;

