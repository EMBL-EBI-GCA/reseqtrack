package ReseqTrack::Tools::RunAlignment::BFAST;

use strict;
use warnings;
use vars qw(@ISA);
use Data::Dumper;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::RunAlignment;
use ReseqTrack::Tools::FileSystemUtils qw(delete_directory);

@ISA = qw(ReseqTrack::Tools::RunAlignment);

sub new {
 my ( $class, @args ) = @_;
 my $self = $class->SUPER::new(@args);

 my $bfast = "/nfs/1000g-work/G1K/work/bin/bfast-0.6.4d/bin/bfast";
 my $prebfast =
"/nfs/1000g-work/G1K/work/bin/bfast-0.6.4d/bin/1kg2bfastfastq/1kg_2_bfast_fastq";

 my ( $build, $space, $aln_options, $preprocess_exe, $read_length ) =
   rearrange( [qw(BUILD SPACE  ALN_OPTIONS PREPROCESS_EXE READ_LENGTH)],
              @args );

 $self->build("grc37");
 $self->program($bfast) if ( !$self->program );
 $self->preprocess_exe($prebfast) if ( !$preprocess_exe );
 $self->build($build) if ($build);    # only have GRC37 indexes at moment
 $self->read_length($read_length);
 $self->create_tmp_process_dir();

 return $self;
}

sub create_cmds {
 my $self = shift;
 $self->set_options();
 $self->set_read_length() if ( !$self->read_length );
 $self->set_indexes();
 $self->create_fastq_cmd_line();
 $self->create_match_cmd_line();
 $self->create_localalign_cmd_line();
 $self->create_postprocess_cmd_line();
}

sub run {
 my ($self) = @_;
 my $exit;

 $self->create_cmds();

 $self->change_dir();

 print "Working in ", $self->working_dir(), "\n";

 # 'match procedure, incorporates preprocess pipe
 $exit = system( $self->match_cmd );
 if ( $exit && $exit >= 1 ) {
  throw( "Failed to run " . $self->match_cmd );
 }

 #localalign
 eval {
  $exit = system( $self->localalign_cmd );
  if ( $exit && $exit >= 1 ) {
   throw( "Failed to run " . $self->localalign_cmd );
  }
 };
 if ($@) {
  throw("Failed to run bwa sampe alignment $@");
 }

 #postprocess
 eval {
  $exit = system( $self->postprocess_cmd );
  if ( $exit && $exit >= 1 ) {
   throw( "Failed to run " . $self->postprocess_cmd );
  }
 };
 if ($@) {
  throw("Failed to run bwa sampe alignment $@");
 }

 print "REALLY ??????? I made it to here . Wow";

 print "Delete tmp dir\n";
 delete_directory( $self->working_dir() );

}




sub set_read_length {
 my $self = shift;
 my $r_length;

 if ( $self->read_length ) {
  print "Read Length already set. Skipping";
  return;
 }

 # Assuming all reads are about the same length.

 # $self->fragment_file("/home/smithre/bobo.gz");
 #print  $self->fragment_file,"\n";

 if ( $self->fragment_file ) {

  my $frag = $self->fragment_file();

  open( my $FQ, "zcat $frag |" ) || die "Could not open: $frag";
  my $line_id   = <$FQ>;
  my $line_seq  = <$FQ>;
  my $line_sep  = <$FQ>;
  my $line_qual = <$FQ>;
  close($FQ);
  $r_length = length($line_seq);
  print "read length = $r_length \n";
  $self->read_length($r_length);
  return;
 }

 if ( $self->mate1_file() ) {
  my $frag = $self->mate1_file();
  open( my $FQ, "zcat $frag |" ) || die "Could not open: $frag";
  my $line_id   = <$FQ>;
  my $line_seq  = <$FQ>;
  my $line_sep  = <$FQ>;
  my $line_qual = <$FQ>;
  close($FQ);
  $r_length = length($line_seq);
  print "read length = $r_length \n";
  $self->read_length($r_length);
  return;
 }

 if ( $self->mate2_file() ) {
  my $frag = $self->mate2_file();
  open( my $FQ, "zcat $frag |" ) || die "Could not open: $frag";
  my $line_id   = <$FQ>;
  my $line_seq  = <$FQ>;
  my $line_sep  = <$FQ>;
  my $line_qual = <$FQ>;
  close($FQ);
  $r_length = length($line_seq);
  print "read length = $r_length \n";
  $self->read_length($r_length);
  return;
 }

 throw "Could not set read_length" if ( !$self->read_length() );

}

sub create_postprocess_cmd_line {
 my $self = shift;
 my $cmd_line;
 my $tmp_dir = $self->working_dir() . '/';

 $cmd_line = $self->program() . " postprocess ";
 $cmd_line .= " -f " . $self->reference() . " ";
 $cmd_line .= " -i  ${tmp_dir}$$.baf  ";
 $cmd_line .= " -n 4 -Q 1000 -t  > ${tmp_dir}$$.sam ";
 print $cmd_line, "\n\n";
 $self->postprocess_cmd($cmd_line);
 $self->sam("$$.sam");
}

sub create_localalign_cmd_line {
 my $self = shift;
 my $cmd_line;

 my $tmp_dir = $self->working_dir() . '/';

 $cmd_line = $self->program() . " localalign ";
 $cmd_line .= " -f " . $self->reference() . " ";

 $cmd_line .= " -m  ${tmp_dir}$$.bmf ";

 $cmd_line .= " -A 1 -o 20 -n 4 -t > $tmp_dir$$.baf ";
 print $cmd_line, "\n\n";
 $self->localalign_cmd($cmd_line);

}

sub create_match_cmd_line {
 my $self = shift;
 my $cmd_line;

 my $tmp_dir = $self->working_dir() . '/';

 $cmd_line = $self->fastq_cmd();
 $cmd_line .= $self->program() . " match ";
 $cmd_line .= " -f " . $self->reference() . " ";
 $cmd_line .= $self->ali_options() . " ";
 $cmd_line .= ">  ${tmp_dir}$$.bmf";

 #print "\n\n", $cmd_line, "\n";
 print $cmd_line, "\n\n";
 $self->match_cmd($cmd_line);

}

sub create_fastq_cmd_line {
 my $self = shift;
 my $cmd_line;

 my $tmp_dir = $self->working_dir() . '/';

 $cmd_line = $self->preprocess_exe() . " ";
 $cmd_line .= $self->mate1_file() . " "    if ( defined $self->mate1_file );
 $cmd_line .= $self->mate2_file() . " "    if ( defined $self->mate2_file );
 $cmd_line .= $self->fragment_file() . " " if ( defined $self->fragment_file );

 $cmd_line .= " \| ";

 #print $cmd_line, "\n\n";
 #print"\n\n", $cmd_line, "\n";
 $self->fastq_cmd($cmd_line);
}

sub set_options {
 my $self = shift;

 #throw "No read length given" if ( !defined( $self->read_length() ) );
 $self->ali_options('-A 1 -K 8 -M 384 -n 4 -Q 25000 ');
}

sub set_indexes {
 my $self = shift;

 throw "No read length given" if ( !defined( $self->read_length() ) );

 if ( $self->build() eq "grc37" ) {

  if ( $self->read_length() < 40 ) {
   $self->bfast_indexes('/nfs/1000g-work/G1K/work/reference/BFAST/short');
   $self->reference(
             '/nfs/1000g-work/G1K/work/reference/BFAST/short/human_g1k_v37.fa');
  }
  else {
   $self->bfast_indexes('/nfs/1000g-work/G1K/work/reference/BFAST/long/');
   $self->reference(
              '/nfs/1000g-work/G1K/work/reference/BFAST/long/human_g1k_v37.fa');
  }
 }
 else {
  throw "Do not have color space ref file for ncbi36";
 }

}

sub bfast_indexes {
 my ( $self, $arg ) = @_;

 if ($arg) {
  $self->{bfast_indexes} = $arg;
 }
 return $self->{bfast_indexes};

}

sub ali_options {
 my ( $self, $arg ) = @_;

 if ($arg) {
  $self->{ali_options} = $arg;
 }
 return $self->{ali_options};
}

sub preprocess_exe {
 my ( $self, $arg ) = @_;

 if ($arg) {
  $self->{preprocess_exe} = $arg;
 }
 return $self->{preprocess_exe};
}

sub read_length {
 my ( $self, $arg ) = @_;

 if ($arg) {
  $self->{read_length} = $arg;
 }
 return $self->{read_length};
}

sub build {
 my ( $self, $arg ) = @_;

 if ($arg) {
  $self->{build} = $arg;
 }
 return $self->{build};
}

sub fastq_cmd {
 my ( $self, $arg ) = @_;
 if ($arg) {
  $self->{fastq_cmd} = $arg;
 }
 return $self->{fastq_cmd};
}

sub match_cmd {
 my ( $self, $arg ) = @_;

 if ($arg) {
  $self->{match_cmd} = $arg;
 }
 return $self->{match_cmd};
}

sub localalign_cmd {
 my ( $self, $arg ) = @_;

 if ($arg) {
  $self->{localalign_cmd} = $arg;
 }
 return $self->{localalign_cmd};
}

sub postprocess_cmd {
 my ( $self, $arg ) = @_;

 if ($arg) {
  $self->{postprocess_cmd} = $arg;
 }
 return $self->{postprocess_cmd};
}

sub space {
 my ( $self, $arg ) = @_;

 if ($arg) {
  $self->{space} = $arg;
 }
 return $self->{space};
}

sub sam {
 my ( $self, $arg ) = @_;

 if ($arg) {
  $self->{sam} = $arg;
 }
 return $self->{sam};
}

sub bam {
 my ( $self, $arg ) = @_;

 if ($arg) {
  $self->{bam} = $arg;
 }
 return $self->{bam};
}

1;

