
=pod

=head1 NAME

ReseqTrack::Tools::RunSamtools

=head1 SYNOPSIS

This is a class for running samtools
It is a sub class of a ReseqTrack::Tools::RunProgram.

=head1 Example


=cut

package ReseqTrack::Tools::RunSamtools;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename qw(fileparse basename);
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists);

use base qw(ReseqTrack::Tools::RunProgram);

=head2 DEFAULT_OPTIONS

  Function  : Called by the RunProgram parent object in constructor
  Returntype: hashref
  Example   : my %options = %{&ReseqTrack::Tools:RunSamtools::DEFAULT_OPTIONS};

=cut

sub DEFAULT_OPTIONS {
    return {
        'use_reference_index' => 0,    # use the -t option when importing a sam
        'input_sort_status' =>
          undef,    # can be 'c' for coordinate or 'n' for name
        'output_sort_status' =>
          undef,    # can be 'c' for coordinate or 'n' for name
        'uncompressed_output' => 0,
        'force_overwrite'     => 1,    # used by samtools merge
        'attach_RG_tag'       => 0,    # used by samtools merge
        'compute_BQ_tag'      => 1,    # used by samtools calmd
        'ext_BAQ_calc'        => 1,    # used by samtools calmd
        'flag_value'          => undef, # used by samtools view (filter)
        'keep_flag'           => 0,     # used by samtools view (filter)
        'remove_flag'         => 0,     # used by samtools view (filter)
        'mapq_threshold'      => 0,     # used by samtools view (filter)
    };
}

=head2 new

  Arg [-reference_index]   :
      string, path of the reference genome .fai file
  Arg [-reference]   :
      string, path of the reference
  + Arguments for ReseqTrack::Tools::RunProgram parent class

  Function  : Creates a new ReseqTrack::Tools::RunSamtools object.
  Returntype: ReseqTrack::Tools::RunSamtools
  Exceptions: 
  Example   : my $run_alignment = ReseqTrack::Tools::RunSamtools->new(
                -input_files => ['/path/sam1', '/path/sam2'],
                -program => "samtools",
                -working_dir => '/path/to/dir/');

=cut

sub new {
    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);

    my ( $reference_index, $reference, ) = rearrange(
        [
            qw( REFERENCE_INDEX REFERENCE
              )
        ],
        @args
    );

    #setting defaults
    if ( !$self->program ) {
        if ( $ENV{SAMTOOLS} ) {
            $self->program( $ENV{SAMTOOLS} . '/samtools' );
        }
        else {
            $self->program('samtools');
        }
    }

    $self->reference_index($reference_index);
    $self->reference($reference);

    return $self;
}

sub run_fix_and_calmd {
    my ($self) = @_;

    foreach my $input ( @{ $self->input_files } ) {
        my $prefix = fileparse( $input, qr/\.[sb]am/ );
        my $output_bam = $self->working_dir . "/$prefix.fixed.md.bam";
        $output_bam =~ s{//}{/}g;

        my @cmds;
        push( @cmds, $self->_get_file_to_sorted_bam_cmd( $input, 1, 1 ) );
        push( @cmds, $self->_get_fixmate_cmd );
        push( @cmds, $self->_get_sort_cmd );
        push( @cmds,
            $self->_get_calmd_cmd( $self->options('uncompessed_output') ) );

        my $cmd = join( ' | ', @cmds ) . " > $output_bam";

        $self->output_files($output_bam);
        $self->execute_command_line($cmd);
    }
}

sub run_calmd {
    my ($self) = @_;

    foreach my $input ( @{ $self->input_files } ) {
        my $prefix = fileparse( $input, qr/\.[sb]am/ );
        my $output_bam = $self->working_dir . "/$prefix.md.bam";
        $output_bam =~ s{//}{/}g;

        my @cmds;
        push( @cmds, $self->_get_file_to_sorted_bam_cmd( $input, 0, 1 ) );
        my $output_sort_status = $self->options('output_sort_status');
        if ( $output_sort_status && $output_sort_status eq 'n' ) {
            push( @cmds, $self->_get_calmd_cmd(1) );
            push( @cmds, $self->_get_sort_cmd(1) );
        }
        else {
            push( @cmds,
                $self->_get_calmd_cmd( $self->options('uncompressed_output') )
            );
        }

        my $cmd = join( ' | ', @cmds ) . " > $output_bam";

        $self->output_files($output_bam);
        $self->execute_command_line($cmd);
    }
}

sub run_bam_filter {
    my ($self) = @_;
    foreach my $input ( @{ $self->input_files } ) {
       my $prefix     = fileparse( $input, qr/\.[sb]am/ );
       my $output_bam = $self->working_dir . "/$prefix.filtered.bam";
       $output_bam =~ s{//}{/}g;
       
       my $keep_flag   = $self->options('keep_flag');
       my $remove_flag = $self->options('remove_flag');
       my $flag_value  = $self->options('flag_value');
       my $uncompressed = $self->options('uncompressed_output');
 
       throw('required option keep_flag or remove_flag when flag_value is provided ') 
            if ($flag_value && !$keep_flag && !$remove_flag);
            
       throw('mutually exclusive options keep_flag and remove_flag')  
            if ($keep_flag && $remove_flag);
       
       my $mapq_threshold = $self->options('mapq_threshold');
       
       my @cmd_words = ( $self->program, 'view' );
       push( @cmd_words, $uncompressed ? '-u' : '-b' );
       push( @cmd_words, '-o', $output_bam ); 
       push( @cmd_words, '-f', $flag_value ) if ( $flag_value && $keep_flag  );    
       push( @cmd_words, '-F', $flag_value ) if ( $flag_value && $remove_flag  );
       push( @cmd_words, '-q', $mapq_threshold ) if $mapq_threshold;
       
       if ( $input =~ /\.sam$/ ){
         push( @cmd_words, '-S', $input )  
       }
       else {
         push( @cmd_words, $input );
       }
 
       my $cmd = join( ' ', @cmd_words );

       $self->output_files($output_bam);
       $self->execute_command_line($cmd);
    }
}

sub run_fixmate {
    my ($self) = @_;

    foreach my $input ( @{ $self->input_files } ) {
        my $prefix = fileparse( $input, qr/\.[sb]am/ );
        my $output_bam = $self->working_dir . "/$prefix.fixed.bam";
        $output_bam =~ s{//}{/}g;

        my @cmds;
        push( @cmds, $self->_get_file_to_sorted_bam_cmd( $input, 1, 1 ) );
        push( @cmds, $self->_get_fixmate_cmd );
        if ( $self->options('output_sort_status') eq 'c' ) {
            push( @cmds, $self->_get_sort_cmd() );
        }

        my $cmd = join( ' | ', @cmds ) . " > $output_bam";

        $self->output_files($output_bam);
        $self->execute_command_line($cmd);
    }
}

sub run_view {
    my ($self) = @;;
}

sub find_reference_index {
    my ($self) = @_;

    if ( !$self->reference_index ) {
        my $reference_index = $self->reference;
        $reference_index =~ s/\.gz$//;
        $reference_index .= '.fai';
        $self->reference_index($reference_index);
    }

    check_file_exists( $self->reference_index );
    return $self->reference_index;
}

sub run_sam_to_bam {
    my ($self) = @_;

    foreach my $input ( @{ $self->input_files } ) {
        my $prefix = fileparse( $input, qr/\.[sb]am/ );

        my $bam = $self->working_dir . "/$prefix.bam";
        $bam =~ s{//}{/}g;

        my $output_sort_status = $self->options('output_sort_status');

        my $cmd;
        if ( $output_sort_status eq 'c' ) {
            $cmd =
              $self->_get_file_to_sorted_bam_cmd( $input, 0,
                $self->options('uncompressed_output') );
        }
        elsif ( $output_sort_status eq 'n' ) {
            $cmd =
              $self->_get_file_to_sorted_bam_cmd( $input, 1,
                $self->options('uncompressed_output') );
        }
        else {
            $cmd =
              $self->_get_sam_to_bam_cmd( $input,
                $self->options('uncompressed_output') );
        }
        $cmd .= " > $bam";

        $self->created_files($bam);
        $self->execute_command_line($cmd);
    }
}

sub run_bam_to_cram {
    my ($self) = @_;

    foreach my $input ( @{ $self->input_files } ) {
        my $cram = $self->working_dir . "/" . basename($input) . ".cram";
        $cram =~ s{//}{/}g;	
        
        my @cmd_words = ($self->program, 'view', '-h', '-C');
        push @cmd_words, '-T', $self->reference;
        push @cmd_words, $input;
        push @cmd_words, ">", $cram;
        
        my $cmd = join(" ", @cmd_words);
        
        $self->output_files($cram);
        $self->execute_command_line($cmd);
    }
}        
	

sub run_sort {
    my ($self) = @_;

    foreach my $input ( @{ $self->input_files } ) {
        my $prefix = fileparse( $input, qr/\.[sb]am/ );
        my $output_bam = $self->working_dir . "/$prefix.sorted.bam";
        $output_bam =~ s{//}{/}g;

        my $cmd =
            $self->options('output_sort_status') eq 'n'
          ? $self->_get_file_to_sorted_bam_cmd( $input, 1 )
          : $self->_get_file_to_sorted_bam_cmd($input);
        $cmd .= " > $output_bam";

        $self->output_files($output_bam);
        $self->execute_command_line($cmd);
    }
}

sub run_index {
    my ($self) = @_;

    foreach my $file ( @{ $self->input_files } ) {
        my $index;
        if ($file =~ /\.cram/) {
            $index = $file . ".crai";
        }
        elsif ($file =~ /\.bam/ ) {
        	$index = $file . ".bai";
    	}
    	else {
    		throw("Cannot decide what index file to create");
    	}
    	
        my $cmd = $self->program . " index " . $file;

        $self->output_files($index);
        $self->execute_command_line($cmd);
    }
}

sub run_flagstat {
    my ($self) = @_;

    my @metrics;
    foreach my $file ( @{ $self->input_files } ) {
        my $output_file = $file . ".flagstat";
        my $cmd = $self->program . " flagstat " . $file . " > $output_file";

        $self->output_files($output_file);

        $self->execute_command_line($cmd);

        open( my $fh, '<', $output_file )
          or throw("Could not open metrics file: $output_file");
        my %data;
        while (<$fh>) {
            chomp;
            if (/^(\d+).*in total/) {
                $data{total_reads} = $1;
            }
            elsif (/^(\d+).*QC failure/) {
                $data{qc_failures} = $1;
            }
            elsif (/^(\d+).*duplicates/) {
                $data{duplicates} = $1;
            }
            elsif (/^(\d+).*mapped \(/) {
                $data{mapped_reads} = $1;
            }
            elsif (/^(\d+).*paired in sequencing/) {
                $data{paired_reads} = $1;
            }
            elsif (/^(\d+).*read1/) {
                $data{read1_reads} = $1;
            }
            elsif (/^(\d+).*read2/) {
                $data{read2_reads} = $1;
            }
            elsif (/^(\d+).*properly paired/) {
                $data{mapped_proper_paired_reads} = $1;
            }
            elsif (/^(\d+).*with itself and mate mapped/) {
                $data{mapped_paired_reads} = $1;
            }
            elsif (/^(\d+).*singletons/) {
                $data{singletons} = $1;
            }
            elsif (/^(\d+).*with mate mapped to a different chr/) {
                $data{mate_mapped_to_other_chr} = $1;
            }
            elsif (/^(\d+).*with mate mapped to a different chr \(/) {
                $data{mate_mapped_to_other_chr_q5} = $1;
            }
        }
        close($fh);
        push @metrics, \%data;

    }
    $self->output_metrics_object( \@metrics );
    return ( \@metrics );
}

sub run_merge {
    my ($self) = @_;

    throw("need more than two or more files to merge")
      if ( @{ $self->input_files } < 2 );

    my $output_bam = $self->working_dir . '/' . $self->job_name . '.merged.bam';
    $output_bam =~ s{//}{/}g;
    throw("file already exists and force_overwrite is not set")
      if ( -e $output_bam && !$self->options('force_overwrite') );

    my $output_sort_status = $self->options('output_sort_status');
    my $input_sort_status  = $self->options('input_sort_status');

    my @cmd_words;
    push( @cmd_words, $self->program, 'merge' );
    push( @cmd_words, '-f' ) if ( $self->options('force_overwrite') );
    push( @cmd_words, '-r' ) if ( $self->options('attach_RG_tag') );
    push( @cmd_words, '-n' ) if ( $output_sort_status eq 'n' );

    push( @cmd_words, $output_bam );

    my $needs_sorting =
      ( $output_sort_status eq 'n' && $input_sort_status ne 'n' );
    $needs_sorting ||=
      ( $output_sort_status ne 'n' && $input_sort_status ne 'c' );

    foreach my $input ( @{ $self->input_files } ) {
        if ( $needs_sorting || $input =~ /\.sam$/ ) {
            my $file_to_sorted_bam_cmd =
              $self->_get_file_to_sorted_bam_cmd( $input,
                ( $output_sort_status eq 'n' ), 1 );
            push( @cmd_words, '<(', $file_to_sorted_bam_cmd, ')' );
        }
        else {
            push( @cmd_words, $input );
        }
    }
    my $cmd = join( ' ', @cmd_words );

    $self->output_files($output_bam);
    $self->execute_command_line($cmd);

}

=head2 run_program

  Arg [1]   : ReseqTrack::Tools::RunSamtools
  Arg [2]   : string, command, must be one of the following:
              merge, sort, index, fix_and_calmd, calmd, sam_to_bam
  Function  : uses samtools to process the files in $self->input_files.
              Output is files are stored in $self->output_files
              This method is called by the RunProgram run method
  Returntype: 
  Exceptions: Throws if the command is not recognised.
  Example   : $self->run();

=cut

sub run_program {
    my ( $self, $command ) = @_;

    my %subs = (
        'merge'         => \&run_merge,
        'sort'          => \&run_sort,
        'index'         => \&run_index,
        'fix_and_calmd' => \&run_fix_and_calmd,
        'calmd'         => \&run_calmd,
        'sam_to_bam'    => \&run_sam_to_bam,
        'bam_to_cram'	=> \&run_bam_to_cram,
        'flagstat'      => \&run_flagstat,
        'filter'        => \&run_bam_filter,
    );

    throw("Did not recognise command $command") if ( !defined $subs{$command} );

    return &{ $subs{$command} }($self);
  }

sub _get_fixmate_cmd {
    my ( $self, $bam ) = @_;
    my @cmd_words = ( $self->program, 'fixmate' );
    push( @cmd_words, $bam || '/dev/stdin' );
    push( @cmd_words, '/dev/stdout' );
    return join( ' ', @cmd_words );
}

sub _get_calmd_cmd {
    my ( $self, $uncompressed, $input ) = @_;
    throw("do not have a reference") if !( $self->reference );
    my @cmd_words = ( $self->program, 'calmd' );
    push( @cmd_words, $uncompressed ? '-u' : '-b' );
    push( @cmd_words, '-r' ) if ( $self->options('compute_BQ_tag') );
    push( @cmd_words, '-E' ) if ( $self->options('ext_BAQ_calc') );

    if ($input) {
        push( @cmd_words, '-S' ) if ( $input =~ /\.sam$/ );
        push( @cmd_words, $input );
    }
    else {
        push( @cmd_words, '-' );
    }
    push( @cmd_words, $self->reference );
    return join( ' ', @cmd_words );
}

sub _get_file_to_sorted_bam_cmd {
    my ( $self, $file, $name_sort, $uncompressed ) = @_;

    return $self->_get_sam_to_sorted_bam_cmd( $file, $name_sort, $uncompressed )
      if ( $file =~ /\.sam$/ );

    my $input_sort_status = $self->options('input_sort_status');
    return $self->_get_sort_cmd( 1, $file )
      if ( $name_sort && ( !$input_sort_status || $input_sort_status ne 'n' ) );

    return $self->_get_sort_cmd( 0, $file )
      if ( !$input_sort_status || $input_sort_status ne 'c' );

    my @cmd_words = ( $self->program, 'view', '-hb' );
    push( @cmd_words, '-u' ) if $uncompressed;
    push( @cmd_words, $file );
    return join( ' ', @cmd_words );
}

sub _get_sam_to_sorted_bam_cmd {
    my ( $self, $sam, $name_sort, $uncompressed ) = @_;

    my $pipe_to_sort = $name_sort && $self->options('input_sort_status') ne 'n';
    $pipe_to_sort ||= !$name_sort && $self->options('input_sort_status') ne 'c';

    my @cmds =
      ( $self->_get_sam_to_bam_cmd( $sam, $uncompressed || $pipe_to_sort ) );
    if ($pipe_to_sort) {
        push( @cmds, $self->_get_sort_cmd($name_sort) );
    }
    return join( ' | ', @cmds );
}

sub _get_sam_to_bam_cmd {
    my ( $self, $sam, $uncompressed ) = @_;

    my @cmd_words = ( $self->program, 'view', '-bS' );
    push( @cmd_words, '-u' ) if ($uncompressed);
    push( @cmd_words, '-t', $self->find_reference_index )
      if ( $self->options('use_reference_index') );
    push( @cmd_words, $sam || '-' );

    return join( ' ', @cmd_words );
}

sub _get_sort_cmd {
    my ( $self, $name_sort, $bam ) = @_;

    my $prefix = $self->get_temp_dir . '/' . $self->job_name;
    $prefix .= '.' . ( $name_sort ? 'nsort' : 'csort' );

    my @cmd_words = ( $self->program, 'sort', '-o' );
    push( @cmd_words, '-n' ) if ($name_sort);
    push( @cmd_words, $bam || '-' );
    push( @cmd_words, $prefix );

    return join( ' ', @cmd_words );
}

=head2 reference

  Arg [1]   : ReseqTrack::Tools::RunSamtools
  Arg [2]   : string, optional, path of genome reference
  Function  : accessor method for reference
  Returntype: string
  Exceptions: n/a
  Example   : my $reference = $self->reference;

=cut

sub reference {
    my ( $self, $reference ) = @_;
    if ($reference) {
        $self->{'reference'} = $reference;
    }
    return $self->{'reference'};
}

=head2 reference_index

  Arg [1]   : ReseqTrack::Tools::RunSamtools
  Arg [2]   : string, optional, path of reference fai file
  Function  : accessor method for reference_index
  Returntype: string
  Exceptions: n/a
  Example   : my $reference_index = $self->reference_index;

=cut

sub reference_index {
    my ( $self, $reference_index ) = @_;
    if ($reference_index) {
        $self->{'reference_index'} = $reference_index;
    }
    return $self->{'reference_index'};
}

=head2 output_metrics_object

  Arg [1]   : ReseqTrack::Tools::RunSamtools
  Arg [2]   : Array, optional, metrics array
  Function  : accessor method for metrics array
  Returntype: Array
  Exceptions: n/a
  Example   : my $reference_index = $self->output_metrics_object;

=cut

sub output_metrics_object {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'output_metrics_object'} = $arg;;
  }
  return $self->{'output_metrics_object'};
}

1;

