=pod

=head1 NAME

ReseqTrack::Tools::RunCramtools

=head1 SYNOPSIS

This is a class for running Cramtools
It is a sub class of a ReseqTrack::Tools::RunProgram.

=head1 Example


=cut

package ReseqTrack::Tools::RunCramtools;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename qw(fileparse);
use File::Copy qw(move);
use ReseqTrack::Tools::FileSystemUtils qw(check_executable check_file_exists);

use base qw(ReseqTrack::Tools::RunProgram);

=head2 DEFAULT_OPTIONS

  Function  : Called by the RunProgram parent object in constructor
  Returntype: hashref
  Example   : my %options = %{&ReseqTrack::Tools:RunCramtools::DEFAULT_OPTIONS};

=cut

sub DEFAULT_OPTIONS {
    return {
     'enumerate' => 1,
     'gzip'      => 1,
     'ignore-md5-mismatch' => 0,
     'skip-md5-check'  => 0,
     'reverse' => 0,
     'reference-fasta-file' => '',
     'prefix' => '',
    };
}    

sub CMD_MAPPINGS { 
    return {
     'fastq' => \&run_convert_fastq,
    };
}    

=head2 new

  Arg [-java_exe]   :
      string, the java executable, default is 'java'
  Arg [-jvm_options]   :
      string, options for java, default is '-Xmx4g'
  Arg [-library_layout]  :
      string, sequencing library type
  
  Function  : Creates a new ReseqTrack::Tools::RunCramtools object.
  Returntype: ReseqTrack::Tools::RunCramtools
  Exceptions: 
  Example   : my $run_cramtools = ReseqTrack::Tools::RunCramtools->new(
                -input_files => ['/path/cram1', '/path/cram2'],
                -program => '/path/to/cramtools.jar',
                -working_dir => '/path/to/dir/',
                -library_layout => 'PAIRED',
               );

=cut

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my ( $java_exe,  $jvm_options , $library_layout ) =
    rearrange(
    [
      qw( JAVA_EXE JVM_OPTIONS LIBRARY_LAYOUT )
    ],
    @args
    );

  $self->java_exe( $java_exe         || 'java' );
  $self->jvm_options( defined $jvm_options ? $jvm_options : '-Xmx4g' );
  $self->library_layout( $library_layout );
  
  return $self;
}


=head2 run_program

  Arg [1]   : ReseqTrack::Tools::RunCramtools
  Arg [2]   : string, command, must be one of the following:
              'fastq' 
  Function  : uses Cramtools to process the files in $self->input_files.
              Output is files are stored in $self->output_files
              This method is called by the RunProgram run method
  Returntype: 
  Exceptions: Throws if it cannot find the picard_dir. Throws if command is not recognised.
  Example   : $self->run('fastq');

=cut

sub run_program {
  my ( $self, $command ) = @_;

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


=head2 run_convert_fastq

  Arg [1]   : ReseqTrack::Tools::RunCramtools
  Function  : uses Cramtools.jar to convert cram files to fastq.gz
  Returntype: 
  Exceptions: 
  Example   : $self->run_convert_fastq

=cut


sub run_convert_fastq {
  my $self = shift;
 
  throw( "expecting single input" ) unless scalar @{ $self->input_files } == 1;
  my $input = ${ $self->input_files }[0];
 
  my $job_name = $self->job_name;
 
  my $prefix_name = $self->options('prefix');
  my $prefix;
 
  if ( $prefix_name ) {
    $prefix = $self->working_dir . '/' . $prefix_name;
  }
  else {
    $prefix = $self->working_dir . '/' . $job_name ;
  }

  $prefix =~ s{//}{/}g;

  my $temp_prefix = $self->get_temp_dir;  
  $temp_prefix = $temp_prefix . '/' . $job_name;
  $temp_prefix =~ s{//}{/}g;
  
  my @cmd_words = ( $self->java_exe );
  push( @cmd_words, $self->jvm_options ) if ( $self->jvm_options );
  push( @cmd_words, '-jar', $self->program , 'fastq' );
 
  push( @cmd_words, '-I' , $input );
  push( @cmd_words, '-F', $temp_prefix );

  push( @cmd_words, '--enumerate'  ) if $self->options('enumerate');
  push( @cmd_words, '-z' ) if $self->options('gzip');
  push( @cmd_words, '--skip-md5-check' ) if $self->options('skip-md5-check');
  push( @cmd_words, '--reverse' ) if $self->options('reverse');
  
  push( @cmd_words, '-R',  $self->options('reference-fasta-file') ) 
                       if $self->options('reference-fasta-file');
                       
  push( @cmd_words, '--read-name-prefix',  $self->options('read-name-prefix') ) 
                       if $self->options('read-name-prefix');
  
  
  my $cmd = join( ' ', @cmd_words );
  
  $self->execute_command_line($cmd);

  
  my $library_layout = $self->library_layout;
  $library_layout = uc($library_layout);
  
  my $from_output;
  my $to_output;
  
  if( $library_layout eq 'SINGLE' ){
    $from_output = $temp_prefix.'.fastq';
    $from_output = $from_output. '.gz' if $self->options('gzip');
 
    $to_output = $prefix.'.fastq';
    $to_output = $to_output . '.gz' if $self->options('gzip');

    move( $from_output, $to_output ) or throw("could not move $from_output to $to_output $!");
  }
  elsif ( $library_layout eq 'PAIRED' ){
    my $from_output_1 = $temp_prefix.'_1.fastq';
    my $from_output_2 = $temp_prefix.'_2.fastq';

    my $to_output_1 = $prefix.'_1.fastq';
    my $to_output_2 = $prefix.'_2.fastq';
    
    if ( $self->options('gzip') ){
      $from_output_1 = $from_output_1 .'.gz';
      $from_output_2 = $from_output_2 .'.gz';
    
      $to_output_1 = $to_output_1 .'.gz';
      $to_output_2 = $to_output_2 .'.gz';
      
      move( $from_output_1, $to_output_1 ) or throw("could not move $from_output_1 to $to_output_1 $!");
      move( $from_output_2, $to_output_2 ) or throw("could not move $from_output_2 to $to_output_2 $!");
      
      push @{$to_output}, $to_output_1, $to_output_2;
    }
    else {
     throw( "library_layout $library_layout not recognized" );
    }
  }
  
  $self->output_files($to_output);
  
  return;
}


sub java_exe {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'java_exe'} = $arg;
  }
  return $self->{'java_exe'};
}

sub jvm_options {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'jvm_options'} = $arg;
  }
  return $self->{'jvm_options'};
}

sub library_layout {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'library_layout'} = $arg;
  }
  return $self->{'library_layout'};
}

1;
