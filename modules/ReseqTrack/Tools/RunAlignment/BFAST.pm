package ReseqTrack::Tools::RunAlignment::BFAST;

use strict;
use warnings;
use vars qw(@ISA);
use Data::Dumper;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::RunAlignment;
use ReseqTrack::Tools::FileSystemUtils qw(delete_directory check_files_exists);

@ISA = qw(ReseqTrack::Tools::RunAlignment);

sub new {
 my ( $class, @args ) = @_;
 my $self = $class->SUPER::new(@args);

 my ( $read_length, $short_reference, $long_reference,$samtools,
      $program, $preprocess_exe, )
   = rearrange(
  [
   qw(
     READ_LENGTH
     SHORT_REFERENCE
     LONG_REFERENCE
     SAMTOOLS
     PROGRAM
     PREPROCESS_EXE
     )
  ],
  @args
   );

 $self->program($program) if ( !$self->program );
 $self->preprocess_exe($preprocess_exe);
 $self->read_length($read_length);
 $self->create_tmp_process_dir();
 $self->short_reference($short_reference);
 $self->long_reference($long_reference);
 $self->snpbin($snpbin);
 $self->samtools($samtools);

 return $self;
}

sub create_cmds {
 my $self = shift;
 $self->set_options();
 $self->set_read_length() if ( !$self->read_length );
 $self->set_reference();
 $self->create_fastq_cmd_line();
 $self->create_match_cmd_line();
 $self->create_localalign_cmd_line();
 $self->create_postprocess_cmd_line();
}

sub run {
 my ($self) = @_;
 my $exit;

 $self->create_cmds();
 $self->create_bam();
 $self->change_dir();

 print "Working in ", $self->working_dir(), "\n";

 # 'match procedure, incorporates preprocess pipe
 eval {
  $exit = system( $self->match_cmd );
  if ( $exit && $exit >= 1 ) {
   throw( "Failed to run " . $self->match_cmd );
  }
 };
 if ($@) {
  throw("Failed to run bfast match: $@");
 }

 #localalign
 eval {
  $exit = system( $self->localalign_cmd );
  if ( $exit && $exit >= 1 ) {
   print STDERR ( "Failed to run " . $self->localalign_cmd );
  }
 };
 if ($@ ||  $exit >= 1) {
  throw("Failed to run bfast localalign: $@");
 }

 #postprocess
 eval {
  $exit = system( $self->postprocess_cmd );
  if ( $exit && $exit >= 1 ) {
   print STDERR( "Failed to run " . $self->postprocess_cmd );
  }
 };
 if ($@ || $exit > 1) {
  throw("Failed to bfast postprocess:  $@");
 }

$self->create_bam();

 print "REALLY ??????? I made it to here . Wow";

 print "Delete tmp dir\n";
 #delete_directory( $self->working_dir() );

}


############
sub create_bam {
 my $self = shift;
 
 my $sam = $self->sam();
 
 my $make_bam = $self->samtools . "  import " . $self->reference  . " $sam  $$.bam ";
 
 print "$make_bam\n";
 $self->bam("$$.bam");
 print $self->bam,"\n";
}



sub set_read_length {
 my $self = shift;
 my $r_length;

 if ( $self->read_length ) {
  print "Read Length already set. Skipping";
  return;
 }

 # Assuming all reads are about the same length.
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
######################

sub set_reference {
 my $self = shift;

 throw "No read length given" if ( !defined( $self->read_length() ) );

 if ( $self->read_length() < 40 ) {
  my $ref = $self->short_reference;
  $self->reference($ref);
 }
 else {
  my $ref = $self->long_reference;
  $self->reference($ref);
 }
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

sub reference {
 my ( $self, $arg ) = @_;
 if ($arg) {
  $self->{reference} = $arg;
 }
 return $self->{reference};
}

sub short_reference {
 my ( $self, $arg ) = @_;
 if ($arg) {
  $self->{short_reference} = $arg;
 }
 return $self->{short_reference};
}

sub long_reference {
 my ( $self, $arg ) = @_;
 if ($arg) {
  $self->{long_reference} = $arg;
 }
 return $self->{long_reference};
}

sub samtools {
 my ( $self, $arg ) = @_;
 if ($arg) {
  $self->{samtools} = $arg;
 }
 return $self->{samtools};
}

sub snpbin {
 my ( $self, $arg ) = @_;
 if ($arg) {
  $self->{snpbin} = $arg;
 }
 return $self->{snpbin};
}

sub program {
 my ( $self, $arg ) = @_;
 if ($arg) {
  $self->{program} = $arg;
 }
 return $self->{program};
}

1;

1;
