package ReseqTrack::Tools::GetFastq::ENAFuse;

use strict;
use warnings;

use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::ERAUtils qw ();
use Digest::SHA qw(sha512_hex);
use Net::FTP;

use base qw(ReseqTrack::Tools::GetFastq);

=pod

=head1 NAME

ReseqTrack::Tools::GetFastq::ENAFuse

=head1 SYNOPSIS

class for getting fastq files through the ena fuse layer. Child class of ReseqTrack::Tools::GetFastq

=head1 Example

my $fastq_getter = ReseqTrack::Tools::GetFastq::ENAFuse(
                      -output_dir => '/path/to/dir'
                      -run_info => $my_run,
                      -db => $era_db,
                      -source_root_dir => '/mount/ena/dir',
                      -fuse_user => 'username',
                      -fuse_password => 'mypassword',
                      -fuse_user => 'username',
                      -fuse_ftphost => 'fusehost.ebi.ac.uk',
                      );
$fastq_getter->run;
my $output_file_list = $fastq_getter->output_files;

=cut


sub new {
    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);

    my ( $fuse_password, $fuse_user, $fuse_ftphost, $fuse_mount_dir, $source_root_dir)
        = rearrange( [ qw( FUSE_PASSWORD FUSE_USER FUSE_FTPHOST FUSE_MOUNT_DIR SOURCE_ROOT_DIR)], @args);

    $self->fuse_password($fuse_password);
    $self->fuse_user($fuse_user);
    $self->fuse_ftphost($fuse_ftphost);
    $self->fuse_mount_dir($fuse_mount_dir || '/protected');

    # This overrides the default set in ReseqTrack::Tools::GetFastq
    $self->source_root_dir($source_root_dir || '/era-priv');

    return $self;
}

sub check_status {
  my $self = shift;
  my $era_rmia = $self->db->get_ERARunMetaInfoAdaptor;
  my $status = $era_rmia->get_status($self->run_info->source_id);
  if ($status ne 'private' && $status ne 'public') {
    print "do not recognise status $status for ".$self->run_info->source_id."\n";
    return 0;
  }
  return 1;
}

sub get_fastq_details {
  my $self = shift;
  throw("do not have a fuse user name") if !$self->fuse_user;
  throw("do not have a fuse password") if !$self->fuse_password;

	my $study_id = $self->run_info->experiment->study->source_id;

  my $string = '/' . $self->fuse_user . '/' . $study_id
            . '/' . $self->fuse_password . '/';
  my $digest = sha512_hex($string);

  ###### Alternative version if sha512_hex does not work
  #my $digest = `echo -n '$string' |  sha512sum`;
  #chomp $digest;
  #$digest =~ s/\s.*//;
  ######

  my $root_dir .= join('/', $self->fuse_mount_dir,
                        $self->fuse_user,
                        $study_id,
                        $digest,
                        $self->source_root_dir);
  $root_dir =~ s{//}{/}g;

  my ($db_md5_hash, $db_size_hash, $name_hash) =
      ReseqTrack::Tools::ERAUtils::get_fastq_details($study_id, $self->db, $root_dir);
  $self->db_md5_hash($db_md5_hash);
  $self->db_size_hash($db_size_hash);
  $self->name_hash($name_hash);
}

sub get_files {
  my $self = shift;
  if (! $self->fuse_ftphost) {
    warn("no fuse_ftphost, so will get files using 'copy'");
    return $self->SUPER::get_files;
  }

  my $ftp = Net::FTP->new($self->fuse_ftphost)
    or throw('cannot connect to '.$self->fuse_ftphost.": $@");
  $ftp->login('anonymous')
    or throw('cannot login to '.$self->fuse_ftphost.' as anonymous: '.$ftp->message);
  $ftp->binary;

  my $output_hash = $self->output_hash;
  while (my($source_path, $output_path) = each %$output_hash) {
    $ftp->get($source_path, $output_path)
        or throw("Failed ftp get from $source_path on ".$self->fuse_ftphost." to $output_path: ".$ftp->message);
    chmod 0644, $output_path;
    $self->output_files($output_path);
  }
}


sub fuse_password{
  my ($self, $fuse_password) = @_;
  if(defined($fuse_password)){
    $self->{fuse_password} = $fuse_password;
  }
  return $self->{fuse_password};
}

sub fuse_user{
  my ($self, $fuse_user) = @_;
  if(defined($fuse_user)){
    $self->{fuse_user} = $fuse_user;
  }
  return $self->{fuse_user};
}

sub fuse_ftphost{
  my ($self, $fuse_ftphost) = @_;
  if(defined($fuse_ftphost)){
    $self->{fuse_ftphost} = $fuse_ftphost;
  }
  return $self->{fuse_ftphost};
}

sub fuse_mount_dir{
  my ($self, $fuse_mount_dir) = @_;
  if(defined($fuse_mount_dir)){
    $self->{fuse_mount_dir} = $fuse_mount_dir;
  }
  return $self->{fuse_mount_dir};
}




1;

