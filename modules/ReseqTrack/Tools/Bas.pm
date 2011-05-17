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
  

  $self->set_required_vars;


  if ( ! $self->bam_md5 ){
   $self->bam_md5 ( run_md5 ($self->bam) );
  }

  return $self;
}

sub run {
 my $self = shift;

 $self->tmp_process_dir;

 if ($self->need_tags){
   $self->add_tags;
 }
   
 $self->create_bas;
 $self->correct_bas_file_convention;
 delete_directory ($self->tmp_dir);

 return;
}

=head2 create_bas

  Arg [1]   : ReseqTrack::Bas
  Function  : Run Sanger module that creates bas file. 
  Returntype: None

=cut

sub create_bas {

  my $self = shift;
  my $bam  = $self->bam_to_process;

  my $perl = $self->perl_exe;
  my $release_date = $self->release_date;
  my $bas = $bam .'.bas';

  $self->tmp_bas ($bas);
 
  my $make_bas_file_cmd ="$perl -MVertRes::Utils::Sam -e \"VertRes::Utils::Sam->new->bas('$bam', '$release_date', '$bas')\" ";
  print $make_bas_file_cmd,"\n" if $self->verbose;

   eval{
     `$make_bas_file_cmd`;
   };
   die "Adding tags failed: $@" if $@;

  return;
}


=head2 parse_study_name

  Arg [1]   : ReseqTrack::Bas
  Function  : Column 3 of bas file needs correct study type
              This extracts from bam name. Assuming correctly formated.
              
  Returntype: None

=cut

sub parse_study_name{
  my $self = shift;
  my @aa = split /\./, $self->bam;
  print "@aa\n" if $self->verbose;

  $self->study_name ( $aa[-3]);
  print "study name = ",  $self->study_name,"\n" if $self->verbose;
  return;
}

=head2 add_tags

  Arg [1]   : ReseqTrack::Bas
  Function  : For correct stats bam needs to have correct tags added.
              Have to run through samtools 'fillmd'
  Returntype: None

=cut

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

=head2   tmp_process_dir

  Arg [1]   : ReseqTrack::Bas
  Function  : Need tmp dir for processing
  Returntype: None
 
=cut

sub tmp_process_dir {
 my $self = shift;
 my $tmp;

 $tmp =  create_tmp_process_dir ($self->working_dir); 
 $self->tmp_dir ( $tmp);

 return;
}

=head2  correct_bas_file_convention

  Arg [1]   : ReseqTrack::Bas
  Function  : Add correct bam file name and study name to tmp bas file
  Returntype: None
  Exceptions: none
 
=cut

sub correct_bas_file_convention{

  my $self = shift;
  
  my $bam_name  = basename ($self->bam_to_process);
  my $bas_file  = $self->tmp_bas;
  my $bam_md5   = $self->bam_md5;
  my $study_name = $self->study_name;

  my @data;
  my @hold;
  my $newline;

  $bam_name =~   s/\.bam//g; # goes in col 1 of bas file


  open (IN, '<' ,$bas_file) || die "Failed to open $bas_file";


  while (<IN>) {
   
    @data = split /\t/; 
  

    if ($data[0] eq "bam_filename"){
      $newline = join ("\t",@data);
      push (@hold, $newline);
      next;
    }

    $data[0] = $bam_name;           # original bam file name
    $data[1] = $bam_md5 ;      # not tmp bam md5 if you added tags 
    $data[2] = $study_name;
 
    my $newline = join ( "\t", @data);

    push (@hold, $newline);
    
  }
  close (IN);

   my $tmp_file = $bas_file . ".org";
   move ($bas_file, $tmp_file);

  open (OUT, '>' ,"$bas_file") || die "Failed to open $bas_file for rewrite";
  foreach my $i (@hold) {
    print OUT $i;
  }
  close (OUT);

  
  my $new_location = $self->bam . '.bas';

  move ( $bas_file, $new_location);

  print "See $new_location for new bas file\n" ;

  chmod(0775, $bas_file);

  return; 
}

=head2  set_required_vars

  Function  : Set ENV variables + perl path so code runs. Only run with
  local perl installation at moment
  Returntype:None
 
=cut

sub set_required_vars {
 my $self = shift;
  
  $ENV{'PERL5LIB'} = '/nfs/1000g-work/G1K/work/bin/local-perl/local-lib/lib/perl5/x86_64-linux:/nfs/1000g-work/G1K/work/bin/local-perl/local-lib/lib/perl5:/homes/smithre/OneKGenomes/reseqtrack/modules:/nfs/1000g-work/G1K/work/bin/VertebrateResequencing-vr-codebase-3ddb4db/modules/';

# /nfs/1000g-work/G1K/work/bin/VertebrateResequencing-vr-codebase-3ddb4db/modules/VertRes

  $self->samtools('/nfs/1000g-work/G1K/work/bin/samtools_latest/samtools/samtools'); 
  

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

1;



