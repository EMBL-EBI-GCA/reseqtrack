
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
use ReseqTrack::Tools::FileSystemUtils qw(check_executable);
use File::Basename qw(fileparse);
use File::Copy qw (move);
use IPC::System::Simple qw( capture );
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
        'ref_flat' =>
          undef
        , # gene annotations file in ref flat format, used by CollectRnaSeqMetrics
        'ribosomal_intervals' =>
          undef
        , # Location of rRNA sequences in genome, in interval_list format, used by CollectRnaSeqMetrics
        'reference_sequence' =>
          undef
        ,    # Location of reference fasta file, used by CollectRnaSeqMetrics
        'strand_specificity' => 'NONE'
        ,    #For strand-specific library prep,  used by CollectRnaSeqMetrics
	'restore_quality' => 1, # restore original quality value, for revert_sam
	'remove_duplicate_info' => 1, # remove markduplicate info, for revert_sam
	'remove_alignment_info' => 1, # remove alignment info, for revert_sam
	'paired_end' => 1, #sam_to_fastq, treat reads in BAM as paired end
	'output_per_rg' => 1, # samtofastq, always output fastq per @RG
    };
}


sub CMD_MAPPINGS { return {
    'mark_duplicates'   => \&run_mark_duplicates,
    'fix_mate'          => \&run_fix_mate,
    'merge'             => \&run_merge,
    'sort'              => \&run_sort,
    'alignment_metrics' => \&run_alignment_metrics,
    'rna_seq_metrics'   => \&run_rna_alignment_metrics,
    'revert_sam'	=> \&run_revertsam,
    'sam_to_fastq'	=> \&run_samtofastq
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

    my ( $java_exe, $picard_dir, $jvm_options, $create_index, $keep_metrics, )
      = rearrange(
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
=head2 run_samtofastq

  Arg [1]   : ReseqTrack::Tools::RunPicard
  Arg [2]   : String, sort order, defaults to coordinate
  Function  : uses SamToFastq.jar to convert BAM to Fastq
  Returntype: 
  Exceptions: 
  Example   : $self->run_samtofastq

=cut


sub run_samtofastq {
    my $self = shift;
    my @fastq_file;
    
    my $jar = $self->_jar_path('SamToFastq.jar');
    my $paired_end=$self->options('paired_end'); ## method to get paired_end status from db
    my $output_dir = '';
    my $output_per_rg=$self->options('output_per_rg'); 

    foreach my $input ( @{ $self->input_files } ) {
        my $name = fileparse( $input, qr/\.[sb]am/ );
        $name =~ s{#}{_}g;
        #my $prefix_dir = $self->working_dir . '/' . $name .'/';
        my $prefix_dir = $self->working_dir .'/';
        $prefix_dir =~ s{//}{/}g;
                
        if($output_per_rg)
        {
            $output_dir = $prefix_dir;
        }
	system("mkdir -p $output_dir"); ## method to create dir
               
        my @cmd_words = ( $self->java_exe );
        push( @cmd_words, $self->jvm_options ) if ( $self->jvm_options );
        push( @cmd_words, '-jar', $jar );
        push( @cmd_words, $self->_get_standard_options );
        push( @cmd_words, 'INPUT=' . $input );
        push( @cmd_words,
            'OUTPUT_PER_RG='. ( $self->options('output_per_rg') ? 'true' : 'false' ) );
              
        if($output_per_rg)
        {
            push( @cmd_words, 'OUTPUT_DIR='.$output_dir );
            
            my $rg_id_array=get_rg_name($input);
            
            
            foreach my $rg_id(@{$rg_id_array})
            {
                if($paired_end)
                {
                    push(@fastq_file,$output_dir.$rg_id."_1.fastq");
                    push(@fastq_file,$output_dir.$rg_id."_2.fastq");
                }
                else
                {
                     push(@fastq_file,$output_dir.$rg_id.".fastq");               
                }
            }
        }
        
        my $cmd = join( ' ', @cmd_words );
        $self->execute_command_line($cmd);
           
      }
	my @gz_fastq_file;

	foreach my $fastq(@fastq_file)
	{
		my $gz_fastq=compress_file($fastq);	
		push @gz_fastq_file, $gz_fastq;
	}  
    	return \@gz_fastq_file;
}

sub get_rg_name {
        my $bam=shift;
        my $samtools_cmd="samtools view -H $bam|grep ^\@RG";
        my @rg_line=capture($samtools_cmd);
        my $rg_id=get_rg_id(\@rg_line);

        return $rg_id;
}

sub get_rg_id {
        my $rg_array=shift;
        my @rg_id;

        foreach my $rg_line(@{$rg_array})
        {
                chomp($rg_line);
                if($rg_line=~ /PU\:(\S+)\s/)
                {
                        push @rg_id,$1;
                        next;
                }
                elsif($rg_line=~ /ID:(\S+)\s/)
                {
                        push @rg_id,$1;
                }
        }
	@rg_id=map {s/#/_/g;$_} @rg_id;
        return \@rg_id;
}

sub compress_file{
  my $file = shift;
  my $cmd = "gzip ".$file;
  my $new_file = $file.".gz";
  if(-e $new_file){
    unlink $new_file;
  }
  my $exit = system($cmd);
  if($exit >= 1){
    throw("There is a problem with ".$cmd." non zero exit code ".$exit);
  }
  unless(-e $new_file){
    throw("Failed to produce ".$new_file." from ".$file." using ".$cmd);
  }
  return $new_file;
}

=head2 run_revertsam

  Arg [1]   : ReseqTrack::Tools::RunPicard
  Arg [2]   : String, sort order, defaults to coordinate
  Function  : uses RevertSam.jar to remove all mapping information in BAM
  Returntype: 
  Exceptions: 
  Example   : $self->run_revertsam

=cut

sub run_revertsam {
    my $self = shift;
    
    my $sort_order = $self->options('sort_order');
    my $restore_qual = $self->options('restore_quality');
    my $remove_dup_info = $self->options('remove_duplicate_info');
    my $remove_aln_info = $self->options('remove_alignment_info');
    
    
    if ( $sort_order ne 'coordinate' && $sort_order ne 'queryname' ) {
        $sort_order = 'null';
    }
    
    my $jar = $self->_jar_path('RevertSam.jar');
    
   my $revert_bam;
 
    foreach my $input ( @{ $self->input_files } ) {
        my $name = fileparse( $input, qr/\.[sb]am/ );
        my $prefix = $self->working_dir . '/' . $name . '_revert';
        $prefix =~ s{//}{/}g;
        my $output = $prefix . '.bam';
   	$revert_bam=$output; 
    
    my @cmd_words = ( $self->java_exe );
        push( @cmd_words, $self->jvm_options ) if ( $self->jvm_options );
        push( @cmd_words, '-jar', $jar );
        push( @cmd_words, $self->_get_standard_options );
        push( @cmd_words, 'SORT_ORDER=' . $sort_order );
        push( @cmd_words, 'INPUT=' . $input );
        push( @cmd_words, 'OUTPUT=' . $output );
	push( @cmd_words,
            'RESTORE_ORIGINAL_QUALITIES='
              . ( $self->options('restore_quality') ? 'true' : 'false' ) );
        push( @cmd_words,
            'REMOVE_DUPLICATE_INFORMATION='
              . ( $self->options('remove_duplicate_info') ? 'true' : 'false' ) );
        push( @cmd_words,
            'REMOVE_ALIGNMENT_INFORMATION='
              . ( $self->options('remove_alignment_info') ? 'true' : 'false' ) );
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
	return $revert_bam;
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
    foreach my $input ( @{ $self->input_files } ) {
        my $suffix = $self->options('remove_duplicates') ? '.rmdup' : '.mrkdup';
        my $name = fileparse( $input, qr/\.[sb]am/ );
        my $prefix = $self->working_dir . '/' . $name . $suffix;
        $prefix =~ s{//}{/}g;
        my $bam     = $prefix . '.bam';
        my $metrics = $prefix . '.dup_metrics';

        my @cmd_words = ( $self->java_exe );
        push( @cmd_words, $self->jvm_options ) if ( $self->jvm_options );
        push( @cmd_words, '-jar', $jar );
        push( @cmd_words, $self->_get_standard_options );
        push( @cmd_words, 'INPUT=' . $input );
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
    
    throw("Reference sequence is required for meaningful metrics") unless $self->options('reference_sequence');
    
    my $jar = $self->_jar_path('CollectAlignmentSummaryMetrics.jar');
    foreach my $input ( @{ $self->input_files } ) {
        my ($input) = @{ $self->input_files };
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

        my $output = $base_name . '.rnaseq_metrics';
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

    my $input_bam_list = $self->input_files;
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
        'ASSUME_SORTED='
          . ( $self->options('assume_sorted') ? 'true' : 'false' ) )
      if defined $self->options('assume_sorted');
    push( @cmd_words,
        'USE_THREADING='
          . ( $self->options('use_threading') ? 'true' : 'false' ) );

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
                undef( $vals[$i] ) if (defined $vals[$i] && $vals[$i] eq '' );
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

=head2 run_fix_mate

  Arg [1]   : ReseqTrack::Tools::RunPicard
  Function  : uses FixMateInformation.jar to ensure mate-pair info is in sync
  Returntype: 
  Exceptions: 
  Example   : $self->run_fix_mate

=cut

sub run_fix_mate{
  my ($self) = @_;
    
  my $sort_order = $self->options('sort_order');
  if ($sort_order ne 'coordinate' && $sort_order ne 'queryname') {
      $sort_order = 'null';
  }
          
  my $jar = $self->_jar_path('FixMateInformation.jar');

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

