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
use ReseqTrack::Tools::AlignmentBase;
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::SequenceIndexUtils;
use ReseqTrack::Tools::FileSystemUtils qw (create_tmp_process_dir run_md5 delete_directory );
use ReseqTrack::Tools::BamUtils qw (CHECK_AND_PARSE_FILE_NAME);
use File::Temp qw/ tempfile tempdir /;
use File::Basename;
use File::Copy;

use vars qw(@ISA);

@ISA = qw(ReseqTrack::Tools::AlignmentBase);

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);
 
  my (
      $need_tags,     
      $working_dir,
      $sequence_index, 
      $release_date,
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
              WORKING_DIR
              SEQUENCE_INDEX
              RELEASE_DATE
              MD5
              VERBOSE
              IN_PARENT
              OUTPUT_DIR
)
		],
		@args
	       );
  
  $self->perl_exe ( '/nfs/1000g-work/G1K/work/bin/local-perl/bin/perl');

  $self->working_dir($working_dir);
  $self->sequence_index($sequence_index);
  $self->need_tags($need_tags);
  $self->release_date($release_date);
  $self->bam_md5($md5);
  $self->parse_study_name;
  $self->in_parent($in_parent);
  $self->output_dir($output_dir);
  $self->set_required_vars;

  if( ! $self->samtools){
    $self->samtools('/nfs/1000g-work/G1K/work/bin/samtools/samtools');
  }


  return $self;
}

sub run {

  my $self = shift;

  if ( ! $self->bam_md5 ) {
    $self->bam_md5 ( run_md5 ($self->bam) );
  }

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
      exit;
  }


  if ($self->in_parent){
    my $new_loc = $self->bam . '.bas';
    print "Moving ", $self->tmp_bas, " ", $new_loc,"\n" if $self->verbose;
    move ( $self->tmp_bas, $new_loc);
    my $sanger_inline = $self->tmp_dir . '/'. "\_Inline";
    $sanger_inline =~ s/\/\//\//g;
    print  $sanger_inline,"\n" if $self->verbose;
    `chmod -R 775 $sanger_inline` if (-e $sanger_inline);

    delete_directory ($self->tmp_dir);
  }
#  delete_directory ($self->tmp_dir);

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
  my $release_date = $self->release_date;

  my $bas =  basename ($bam) .'.bas';
  
  print "Making $bas\n" if $self->verbose;

  $bas =~ s/\/\//\//g;
  $self->tmp_bas ($bas);
  
  my $make_bas_file_cmd ="$perl -MVertRes::Utils::Sam -e \"VertRes::Utils::Sam->new->bas('$bam', '$release_date', '$bas')\" ";
  print $make_bas_file_cmd,"\n" if $self->verbose;

  eval{
    `$make_bas_file_cmd`;
  };
  die "bas creaton failed: $@" if $@;
 
  return;
}




sub parse_study_name{
  my $self = shift;
  my @aa = split /\./, $self->bam;

  my ($sample, $platform, $algorithm, $project, $analysis_grp, $chr, $date);

  $analysis_grp = "UNKNOWN_GROUP";
 
  ($sample, $platform, $algorithm, $project, $analysis_grp, $chr,$date) = 
    CHECK_AND_PARSE_FILE_NAME ($self->bam);

  $self->study_name ( $analysis_grp);
  $self->release_date ( $date);

  print "release date= $date\n" if $self->verbose;
  print "study name = ",  $self->study_name,"\n" if $self->verbose;

  return;
}



sub add_tags {
  my $self = shift;

  my $sam_exe = $self->samtools;
  my $bam     = $self->bam;
  my $ref     = $self->reference;

  my $tmp_bam = $self->tmp_dir . '/' . basename($bam) . '_'. $$;
 
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



sub correct_bas_file_convention{

  my $self = shift;
  
  my $bam_name  = basename ($self->bam_to_process);
  my $bas_file  = $self->tmp_dir . '/' . $self->tmp_bas;
  $bas_file =~ s /\/\//\//g;
  my $bam_md5   = $self->bam_md5;
  my $study_name = $self->study_name;

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
    $data[2] = $study_name;
 
    my $newline = join ( "\t", @data);

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

  return; 
}


sub set_required_vars {
  my $self = shift;
  

$ENV{'PERL5LIB'} = '/nfs/1000g-work/G1K/work/bin/local-perl/local-lib/lib/perl5/x86_64-linux:/nfs/1000g-work/G1K/work/bin/local-perl/local-lib/lib/perl5:/homes/smithre/OneKGenomes/reseqtrack/modules:/nfs/1000g-work/G1K/work/bin/vr-codebase/modules/';

#  print  $ENV{'PERL5LIB'},"\n";


  if ( ! (defined $ENV{'PICARD'} )) {
    $ENV{'PICARD'} = '/nfs/1000g-work/G1K/work/bin/picard/picard_143/picard-tools-1.43/';
  }

  if ( ! (defined $ENV{'SAMTOOLS'} )) {
    $ENV{'SAMTOOLS'} = '/nfs/1000g-work/G1K/work/bin/samtools';
  }

  if ( ! defined ($ENV{'PERL_INLINE_DIRECTORY'})) {
    $ENV{'PERL_INLINE_DIRECTORY'} = '/homes/rseqpipe/.Inline';
  }

  if (! ($ENV{'PATH'} =~ /samtools/) ) {
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

sub study_name {
  my ($self, $arg) = @_;
  if ($arg) {
    $self->{'study_name'} = $arg;
  }
  return $self->{'study_name'};
}

sub release_date {
  
  my ($self, $arg) = @_;
  if ($arg) {
    $self->{'release_date'} = $arg;
  }
  return $self->{'release_date'};

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

sub working_dir{
  my ($self, $arg) = @_;
  if ($arg) {
    $self->{'working_dir'} = $arg;
  }
  return $self->{'working_dir'};
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


1;



