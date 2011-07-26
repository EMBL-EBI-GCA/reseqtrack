package ReseqTrack::Tools::RunVerifyBamID;
 
use strict;
use warnings;
 
use Data::Dumper;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename;
use ReseqTrack::Tools::RunProgram;

use vars qw(@ISA);

@ISA = qw(ReseqTrack::Tools::RunProgram);

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my ( $reference, $claimed_sample, $bimp, $selfonly, $out_prefix,
     $maxdepth) =
    rearrange(
	      [
	       qw(
			  REFERENCE
			  CLAIMED_SAMPLE
			  BIMP
			  SELFONLY
			  OUT_PREFIX
                          MAXDEPTH
			  )
	      ],
	      @args
	     );


  $self->reference($reference);
  $self->bimp($bimp);
  $self->selfonly($selfonly);
  $self->out_prefix($out_prefix);
  $self->maxdepth($maxdepth);

 
  return $self;
}

sub run {
  
  my $self = shift;

  $self->construct_run_cmd();

  $self->change_dir();

  $self->add_outfiles();
 
  $self->execute_command_line($self->command_line);
 


  return;

}

######################

sub construct_run_cmd {
  my $self = shift;
  my $cmd;
  my $out_prefix;

  my $files =  $self->input_files;
  my $bam   = @{$files}[0];

  $cmd .= $self->program . " ";

  $cmd .= "--reference " . $self->reference . " ";

  $cmd .= "--in " . $bam . " ";

  $cmd .= "--bfile " . $self->bimp . " " if ( $self->bimp );

  $cmd .= "--selfonly "  if ( $self->selfonly );

  $cmd .= "--precise -d " . $self->maxdepth  . " "  if ($self->maxdepth);

  if ( !defined( $self->out_prefix ) ) {
    $out_prefix = $self->working_dir . "\/" . basename( $bam );
    $out_prefix =~ s/\.bam//; 
    $out_prefix =~ s/\/\//\//g;
  } else {
    $out_prefix = $self->out_prefix;
  }

  $cmd .= "--out " . $out_prefix . " --verbose ";

  #print $cmd, "\n";
  
  $self->command_line($cmd);

  return;
}

######################
sub command_line {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{command_line} = $arg;
  }
  return $self->{command_line};
}

sub add_outfiles {

  my $self = shift;
  my @output_files;
  my @extensions;

  my $work_dir = $self->working_dir;
  $work_dir =~ s/\/$//;

  my $bam = basename( @{$self->input_files}[0] );
  $bam =~ s /\.bam$//;
  my $file_base = $work_dir . '/' . $bam;

  if ( $self->selfonly ) {
    @extensions = qw (selfSM selfRG );
  } else {
    @extensions = qw (selfSM selfRG bestSM bestRG );
  }

  foreach my $ext (@extensions) {
   
    my $x = $ext . "\_file";
    my $filename =  $file_base . ".$ext";
    $self->$x($filename);
    print $self->$x,"\n";
    push( @output_files, $filename );
  }




  $self->output_files( \@output_files );
 
  return;
}

#sub bam {
#  my ( $self, $arg ) = @_;
#  if ($arg) {
#    $self->{bam} = $arg;
#  }
#  return $self->{bam};
#}

sub reference {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'reference'} = $arg;
  }
  return $self->{'reference'};
}

sub bimp {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{bimp} = $arg;
  }
  return $self->{bimp};
}

sub out_prefix {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{out_prefix} = $arg;
  }
  return $self->{out_prefix};
}

sub selfonly {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{selfonly} = $arg;
  }
  return $self->{selfonly};
}

sub selfSM_file {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{selfSM_file} = $arg;
  }
  return $self->{selfSM_file};
}

sub selfRG_file {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{selfRG_file} = $arg;
  }
  return $self->{selfRG_file};
}

sub bestRG_file {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{bestRG_file} = $arg;
  }
  return $self->{bestRG_file};
}
sub bestSM_file {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{bestSM_file} = $arg;
  }
  return $self->{bestSM_file};
}

sub maxdepth {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{maxdepth} = $arg;
  }
  return $self->{maxdepth};
}


1;
