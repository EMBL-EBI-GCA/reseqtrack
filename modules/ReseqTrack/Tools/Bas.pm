=pod

=head1 NAME

ReseqTrack::Tools::Bas

=head1 SYNOPSIS

A pretty much self-contained  package that should generate a bas file
for a bam.

=cut

package ReseqTrack::Tools::Bas;

use strict;
use warnings;
use Data::Dumper;
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::RunAlignment;
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::SequenceIndexUtils;
use ReseqTrack::Tools::FileSystemUtils qw (create_tmp_process_dir run_md5 delete_directory );
use ReseqTrack::Tools::BamUtils qw (CHECK_AND_PARSE_FILE_NAME);
use File::Temp qw/ tempfile tempdir /;
use File::Basename;
use File::Copy;
 
use vars qw(@ISA);

@ISA = qw(ReseqTrack::Tools::RunAlignment);

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my (
      $need_tags,
      $sequence_index, 
      $md5,
      $verbose,
      $in_parent,
      $output_dir,
     )
    =
      rearrange(
		[
		 qw(
              NEED_TAGS
              SEQUENCE_INDEX
              MD5
              VERBOSE
              IN_PARENT
              OUTPUT_DIR
      )
		],
		@args
	       );

  $self->bam ( @{$self->input_files}[0]) if ( defined @{$self->input_files}[0]) ;



  $self->perl_exe ( '/nfs/1000g-work/G1K/work/bin/local-perl/bin/perl');
  $self->sequence_index($sequence_index);
  $self->need_tags($need_tags);
  $self->bam_md5($md5);
  $self->in_parent($in_parent);
  $self->output_dir($output_dir);

  if ( ! defined ($self->samtools)){

    if ( defined $ENV{SAMTOOLS}){
      my $exe = $ENV{SAMTOOLS} . '/'. "samtools";
      $exe =~ s/\/\//\//g;
      print "Setting path to samtools from ENV\n";
      $self->samtools($exe) if ( -e $exe);
      if (! -e $exe){
	throw "Path to samtools not set. Required for processing\n";
      }
    }
  }

  $self->set_required_vars;
  return $self;
}



sub run {

  my $self = shift;

  if ( ! $self->bam_md5 ) {
    $self->bam_md5 ( run_md5 ($self->bam) );
  }

  $self->get_bam_run_id_info();

  $self->tmp_process_dir;

  if ( $self->need_tags ) {
    $self->add_tags;
  } else {
    $self->bam_to_process ($self->bam);	
  }
   

  $self->create_bas;

  $self->correct_bas_file_convention;

  if ($self->output_dir){
      my $new_loc =$self->output_dir . '/'. basename($self->bam) . '.bas';
       $new_loc =~ s/\/\//\//g;
     
      my $old_loc = $self->tmp_dir . '/' . $self->tmp_bas;
      print "Copy $old_loc to $new_loc\n" if $self->verbose; 
       $old_loc =~ s/\/\//\//g;
      my  $fred = copy ( $old_loc, $new_loc);
      print $fred,"\n" if $self->verbose;
      return;
  }


  if ($self->in_parent){
    my $new_loc = $self->bam . '.bas';
    print "Moving ", $self->tmp_bas, " ", $new_loc,"\n" if $self->verbose;
    move ( $self->tmp_bas, $new_loc);
    my $sanger_inline = $self->tmp_dir . '/'. "\_Inline";
    $sanger_inline =~ s/\/\//\//g;
    print  $sanger_inline,"\n" if $self->verbose;
    `chmod -R 775 $sanger_inline` if (-e $sanger_inline);
    $self->files_to_delete ($self->tmp_dir);
  }


  return;
}

 
sub create_bas {

  my $self = shift;
  my $bam  = $self->bam_to_process;

  if ( !-e $bam) {
    die "Bam file: $bam does not exist\n";
  }

  my $go_here =  $self->tmp_dir;
  print "Changing dir to $go_here\n" if $self->verbose;
  chdir ( $go_here);

  my $perl = $self->perl_exe;


  my $bas =  basename ($bam) .'.bas';

  $bas = $self->tmp_dir . '/' . $bas;
  $bas =~ s/\/\//\//g;
  print "Making $bas\n" if $self->verbose;


  $self->tmp_bas ($bas);

  my $perl_mods =  '-MVRPipe::File -MVRPipe::Steps::bam_stats -e "VRPipe::Steps::bam_stats';

  my $make_bas_file_cmd ="$perl  ${perl_mods}->bas('$bam',  '$bas')\" ";
  print $make_bas_file_cmd,"\n" if $self->verbose;

  eval{
    `$make_bas_file_cmd`;
  };
  die "bas creaton failed: $@" if $@;

  return;
}



=head2 add_tags

  Function  : uses samtools add correct NM tags to bam
  Example   : my $bam = $self->run_sam_to_bam($sam);

=cut

sub add_tags {
  my $self = shift;


  my $bam     = $self->bam;
  my $ref     = $self->reference;

  my $tmp_bam = $self->tmp_dir . '/' . basename($bam) . '_'. $$;



  my $sam_exe = $self->samtools;

  print "Adding tags\n" if $self->verbose;
  print  "$sam_exe fillmd -b  $bam $ref > $tmp_bam\n" if $self->verbose;

  eval{  
    `$sam_exe fillmd -b  $bam $ref  > $tmp_bam `;
  };
  die "Adding tags failed: $@" if $@;

  $self->bam_to_process ($tmp_bam);
  print  $self->bam_to_process,"\n" if $self->verbose;
  return;
}


sub tmp_process_dir {
  my $self = shift;
  my $tmp;

  $tmp =  create_tmp_process_dir ($self->working_dir); 
  $self->tmp_dir ( $tmp);

  return;
}


=head2 correct_bas_file_convention
  Arg [1]   : ReseqTrack::Tools::RunSamtools
  Arg [2]   : string, path of sam file
  Function  : Picard processed bams can have RG info in PU field of header.
              Look in PU field SRR/ERR strings
  Returntype: None
=cut

sub correct_bas_file_convention{
  my $self = shift;

  print "RG check\n";

  my $read_group_info = $self->rg_info;

  my $bam_name  = basename ($self->bam_to_process);
  my $bas_file  = $self->tmp_bas;
  $bas_file =~ s /\/\//\//g;
  my $bam_md5   = $self->bam_md5;


  my @data;
  my @hold;
  my $newline;

  $bam_name =~   s/\.bam//g;	# goes in col 1 of bas file

  open (my $IN, '<' ,$bas_file) || die "Failed to open $bas_file";


  while (<$IN>) {

    @data = split /\t/; 


    if ($data[0] eq "bam_filename") {
      $newline = join ("\t",@data);
      push (@hold, $newline);
      next;
    }

    $data[0] = $bam_name;	# original bam file name
    $data[1] = $bam_md5 ;	# not tmp bam md5 if you added tags 


    if (! ( $data[6] =~ /^[ERR|SRR]/)) {
      print "Invalid run id identifier : $data[6] :: ";

      if (defined  $$read_group_info{$data[6]}) {
	print "Picard processed bam ??? Checking PU designation .. ";	
	my $tmp_rg =  $$read_group_info{$data[6]}{PU};

	if ( $tmp_rg =~ /^[ERR|SRR]/){
	  print "Replacing $data[6] with " , $$read_group_info{$data[6]}{PU},"\n";
	  $data[6] = $$read_group_info{$data[6]}{PU};
	}
      }
      else{
	print "Do not know what to replace bad identifier with. Ignoring\n";
      }
    }
    else{
      print "$data[6] looks OK\n";
    }

    my $newline = join ("\t", @data);

    push (@hold, $newline);
    
  }
  close ($IN);

  my $tmp_file = $bas_file . ".org";
  move ($bas_file, $tmp_file);



  open (my $OUT, '>' ,"$bas_file") || die "Failed to open $bas_file for rewrite";
  foreach my $i (@hold) {
    print $OUT $i;
  }
  close ($OUT);

  my $new_location = $bas_file;
  $new_location =~ s/\.org//;

  move ( $bas_file, $new_location);

  print "See $new_location for new bas file\n" if $self->verbose;

  chmod(0775, $bas_file);

  print "Done RG check\n";

  return;
}




sub set_required_vars {
  my $self = shift;

 $ENV{'PERL5LIB'} = '/nfs/1000g-work/G1K/work/bin/vr-pipe/modules/:/nfs/1000g-work/G1K/work/bin/local-perl/local-lib/lib/perl5/x86_64-linux:/nfs/1000g-work/G1K/work/bin/local-perl/local-lib/lib/perl5';


  if ( ! (defined $ENV{'PICARD'} )) {
    $ENV{'PICARD'} = '/nfs/1000g-work/G1K/work/bin/picard/picard_143/picard-tools-1.43/';
  }

  if ( ! (defined $ENV{'SAMTOOLS'} )) {
    $ENV{'SAMTOOLS'} = '/nfs/1000g-work/G1K/work/bin/samtools';
  }

  if ( ! (defined $ENV{'PERL_INLINE_DIRECTORY'})) {
    $ENV{'PERL_INLINE_DIRECTORY'} = '/homes/rseqpipe/.Inline';
  }

  if (! ($ENV{'PATH'} =~ /samtools/i) ) {
    $ENV{'PATH'}= $ENV{'PATH'} . ':' . $ENV{'SAMTOOLS'};
  }

  return;
}


#####
#Generic get/set methods

sub tmp_dir {
  my ($self, $arg) = @_;
  if ($arg) {
    $self->{'tmp_dir'} = $arg;
  }
  return $self->{'tmp_dir'};
}

sub bam_to_process {

  my ($self, $arg) = @_;
  if ($arg) {
    $self->{'bam_to_process'} = $arg;
  }
  return $self->{'bam_to_process'};

}


sub need_tags{
  my ($self, $arg) = @_;
  if (defined $arg) {
    $self->{'need_tags'} = $arg;
  }
  return $self->{'need_tags'};
}

sub change_dir {
  my ( $self, $dir ) = @_;
  $dir = $self->working_dir unless ($dir);
  chdir($dir)
    or throw( "Failed to change to " . $dir );
  return $dir;
}


sub sequence_index{
  my ($self, $arg) = @_;
  if ($arg) {
    $self->{'sequence_index'} = $arg;
  }
  return $self->{'sequence_index'};
}

sub tmp_bas{
  my ($self, $arg) = @_;
  if ($arg) {
    $self->{'tmp_bas'} = $arg;
  }
  return $self->{'tmp_bas'};
}

sub bas{
  my ($self, $arg) = @_;
  if ($arg) {
    $self->{'bas'} = $arg;
  }
  return $self->{'bas'};
}

sub bam_md5 {
  my ($self, $arg) = @_;
  if ($arg) {
    $self->{'bam_md5'} = $arg;
  }
  return $self->{'bam_md5'};

}

sub perl_exe {
  my ($self, $arg) = @_;
  if ($arg) {
    $self->{'perl_exe'} = $arg;
  }
  return $self->{'perl_exe'};
}

sub verbose {
  my ($self, $arg) = @_;
  if ($arg) {
    $self->{'verbose'} = $arg;
  }
  return $self->{'verbose'};
}

sub in_parent {
  my ($self, $arg) = @_;
  if (defined $arg) {
    $self->{'in_parent'} = $arg;
  }
  return $self->{'in_parent'};
}

sub output_dir {
  my ($self, $arg) = @_;
  if (defined $arg) {
    $self->{'output_dir'} = $arg;
  }
  return $self->{'output_dir'};
}


sub bam {
  my ($self, $arg) = @_;
  if (defined $arg) {
    $self->{'bam'} = $arg;
  }
  return $self->{'bam'};
}



sub get_bam_run_id_info{

#BCM SOLID bams has RG in PU field of bam header. Must catch.
#@RG     ID:1    PL:SOLiD        PU:SRR097880    LB:ANG_TG.HG00142-1_1sA SM:HG00142      CN:BCM
#@RG     ID:2    PL:SOLiD        PU:SRR097881    LB:ANG_TG.HG00142-1_1sA SM:HG00142      CN:BCM
#
# Extract RG head info and change when in sub  correct_bas_file_convention


  my $self = shift;
  my $bam  = $self->bam;
  my $samtools  = $self->samtools;
  my $key = "-";
  my %rg_hash;

  my @rg_info =  `$samtools view $bam -H | grep "^\@RG"`;


  die "No read group info found in $bam\n" if ( scalar  @rg_info < 1);
  print "Found RG info for ",  scalar  @rg_info ," runs\n";

  foreach my $line (@rg_info) {
 #   print $line;
    my @aa = split /\s+/,$line;
    shift @aa;

    my $key = "-";

    foreach my $x (@aa) {
      my ($p1, $p2) = split /:/,$x;
      ($key = $p2) if ( $p1 eq "ID");
    }
    die "No ID in RG row \n" if ( $key eq "-");


    foreach my $x (@aa) {
      my ($p1, $p2) = split /:/,$x;
      next if ($p2 eq $key);
      $rg_hash{$key}{$p1} =  $p2;
    }
  }

  $self->rg_info(\%rg_hash);

  return;
}

sub rg_info {
  my ($self, $arg) = @_;
  if (defined $arg) {
    $self->{'rg_info'} = $arg;
  }
  return $self->{'rg_info'};
}



1;



