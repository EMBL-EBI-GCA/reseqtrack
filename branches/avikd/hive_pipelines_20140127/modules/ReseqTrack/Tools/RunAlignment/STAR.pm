package ReseqTrack::Tools::RunAlignment::STAR;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::RunAlignment;
use Data::Dumper;
use File::Basename;
use File::Copy;

@ISA = qw(ReseqTrack::Tools::RunAlignment);

=pod

=head1 NAME

ReseqTrack::Tools::RunAlignment::STAR

=head1 SYNOPSIS

class for running STAR. Child class of ReseqTrack::Tools::RunAlignment

=head1 Example

my $star = ReseqTrack::Tools::RunAlignment::STAR(
                      -input_files => '/path/to/file',
                      -reference => '/path/to/reference',
                      -options => { ‘outFilterScoreMin’ => 0, ‘outFilterMatchNmin’ => 0 },
                      -read_group_fields => {'ID' => 1, 'LB' => 'my_lib'},
                      -output_format => 'BAM',
                      -samtools => '/path/to/samtools',
                      );
$star->run;
my $output_file_list = $star->output_files;

=cut

sub DEFAULT_OPTIONS {
  return {
    'runMode'                 => 'alignReads',
    'runThreadN'              => 4,
    'genomeLoad'              => 'NoSharedMemory',
    'outStd'                  => 'SAM',
    'outSAMstrandField'       => 'intronMotif',
    'outSAMattributes'        => 'Standard',
    'outSJfilterReads'        => 'All',
    'alignIntronMin'          => 21,
    'alignSJDBoverhangMin'    => 3, 
    'outFilterMultimapNmax'   => 10,  
    'outFilterMismatchNmax'   => 10,
  };
}

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  #setting defaults
  if ( !$self->program ) {
    if ( $ENV{STAR} ) {
      $self->program( $ENV{STAR} . '/STAR' );
    }
    else {
      $self->program('STAR');
    }
  }

  return $self;
}

=head2 run_alignment

  Arg [1]   : ReseqTrack::Tools::RunAlignment::STAR
  Function  : uses STAR to align files in $self->input_files.
              Output is files are stored in $self->output_files
              This method is called by the RunProgram run method
  Returntype: n/a
  Exceptions: 
  Example   : $self->run();

=cut

sub run_alignment {
  my ($self) = @_;

  my $output_file = $self->working_dir() . '/' . $self->job_name;
  $output_file =~ s{//}{/};
  
  throw( "reference directory does not exist" ) unless ( -d $self->reference ); #reference is a directory for STAR
  throw( "runMode is required" ) unless( $self->options('runMode') );
  throw( "output format not supported, use \'outStd=SAM\'" ) unless ( $self->options('outStd') eq 'SAM' );
  throw("unrecognised runMode:" . $self->options('runMode')) unless ( $self->options('runMode') eq 'alignReads' ) ;
  
  
  if ( $self->fragment_file ) {
    $self->_do_alignment( $output_file . '_se', $self->fragment_file );
  }
  if ( $self->mate1_file && $self->mate2_file ) {
    $self->_do_alignment( $output_file . '_pe',
    $self->mate1_file, $self->mate2_file );
  }
}


sub _do_alignment {
  my ( $self, $output_file, @input_files ) = @_;

  #STAR doesn't output BAM, but we can pipe the output through samtools
  my $do_bam_conversion = 0;
  my $output_option     = $self->output_format;
  $do_bam_conversion = 1 if ( lc($output_option) eq 'bam' );
  
  my @restricted_options = qw /genomeDir outSAMattrRGline readFilesIn readFilesCommand/;
  my %restricted_option_list = map { $_ => 1 } @restricted_options;
  
  my @cmd_words;
  push( @cmd_words, $self->program );
  
  #genome parameters
  push( @cmd_words, '--genomeDir', $self->reference );

  
  if ( $self->read_group_fields->{'ID'} ) {
     my $rg_string = q(ID:) . $self->read_group_fields->{'ID'};
    RG:
    while (my ($tag, $value) = each %{$self->read_group_fields}) {
        next RG if ($tag eq 'ID');
        next RG if (!$value);
        if($value =~ /.*\s+/){
          $rg_string .= "\t\"" . $tag . ":" . $value . "\""; # RG format for STAR,  ID:xxx CN:yy "DS:z z z"
        }
        else {
          $rg_string .= "\t" . $tag . ":" . $value;
        }       
    }
    push( @cmd_words, '--outSAMattrRGline', $rg_string );
  }
  
  #inputs
  push( @cmd_words, '--readFilesIn', @input_files );
  push( @cmd_words, '--readFilesCommand', "zcat" )
    if ( $input_files[0] =~ /.gz$/ );

  #other options
  if ( $self->options ) {
    OPTION:
    while (my ($tag, $value) = each %{$self->options}) {
        next OPTION if !defined $value;
        next OPTION if (exists($restricted_option_list{$tag}));
    }
    
    foreach my $tag ( keys ( %{$self->options} ) ) {
      if (defined $self->options->{$tag}) {
        my $value = $self->options->{$tag};
        push (@cmd_words, "--".$tag." ".$value);
      }
    }
  }
   

    if ( $do_bam_conversion ) {
      push( @cmd_words, '|', $self->samtools, 'view -bhS -' );
    }
    $output_file .= '.'.lc($output_option);
    push( @cmd_words, '>', $output_file );
    
    $self->output_files($self->working_dir() . '/' . 'Log.out');
    $self->output_files($self->working_dir() . '/' . 'Log.progress.out');
    $self->output_files($self->working_dir() . '/' . 'Log.std.out');


  
  $self->output_files($output_file);
  my $cmd = join( ' ', @cmd_words );

  $self->execute_command_line($cmd);
  
}

sub output_bam_files {
  my $self = shift;
  my @files = grep { /\.bam$/ } @{ $self->output_files };
  return \@files;
}

1;
