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
  $self->program('STAR') unless $self->program ;
  
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
  throw( "unrecognised runMode:" . $self->options('runMode')) unless ( $self->options('runMode') eq 'alignReads' ) ;
  
  throw( "fragment and mate-pair mixed read not allowed for RNA-seq" ) if( $self->fragment_file && $self->mate1_file && $self->mate2_file );
  
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
  #my $do_bam_conversion = 0;
  #my $output_option     = $self->output_format;
  #$do_bam_conversion = 1 if ( lc($self->output_format) eq 'bam' );
  
  my $do_bam_conversion =  lc($self->output_format) eq 'bam' ? 1 : 0 ;
  
  my @restricted_options = qw /genomeDir outSAMattrRGline readFilesIn readFilesCommand outStd outFileNamePrefix/;
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
    
 #output
 my $output_tmp = $self->get_temp_dir() . '/' . $self->job_name;
 
 push( @cmd_words, '--outStd', "SAM" );
 push( @cmd_words, '--outFileNamePrefix', $output_tmp );
 
  #other options
  if ( $self->options ) {
    OPTION:
    while (my ($tag, $value) = each %{$self->options}) {
        next OPTION if !defined $value;
        throw ( "$tag not permitted" ) if (exists($restricted_option_list{$tag}));

        push (@cmd_words, "--".$tag." ".$value);
    }
  }
  
    if ( $do_bam_conversion ) {
      push( @cmd_words, '|', $self->samtools, 'view -bhS -' );
    }
    
    my $output_tmp_file = $output_tmp .'.' . lc($self->output_format);
    my $sj_tab_tmp = $output_tmp . 'SJ.out.tab';
    
    push( @cmd_words, '>', $output_tmp_file );
  
    my $cmd = join( ' ', @cmd_words );

    $self->execute_command_line($cmd);
    
    my $sj_tab =  $output_file.'.SJ.out.tab';
    
    $output_file .= '.'.lc($self->output_format);
    move( $output_tmp_file, $output_file );
    move($sj_tab_tmp, $sj_tab);
     
    $self->output_files($output_file);
    $self->output_files($sj_tab);
  
}

sub output_bam_files {
  my $self = shift;
  my @files = grep { /\.bam$/ } @{ $self->output_files };
  return \@files;
}

sub output_sam_files {
  my $self = shift;
  my @files = grep { /\.sam$/ } @{ $self->output_files };
  return \@files;
}

sub output_sj_tab_files {
  my $self = shift;
  my @files = grep { /\.SJ.out.tab$/ } @{ $self->output_files };
  return \@files;
}

1;
