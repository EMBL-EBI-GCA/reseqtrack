package ReseqTrack::Tools::GetFastq::ENAFuse;

use strict;
use warnings;

use Digest::SHA qw(sha512_hex);

use base qw(ReseqTrack::Tools::GetFastq);

=pod

=head1 NAME

ReseqTrack::Tools::GetFastq::ENAFuse

=head1 SYNOPSIS

class for getting fastq files through the ena fuse layer. Child class of ReseqTrack::Tools::GetFastq

=head1 Example

my $fastq_getter = ReseqTrack::Tools::GetFastq::ENAFuse(
                      -output_dir => '/path/to/dir'
                      -run_meta_info => $my_rmi,
                      -db => $era_db,
                      -source_root_dir => '/mount/ena/dir',
                      -fuse_user => 'username',
                      -fuse_password => 'mypassword',
                      -fuse_user => 'username',
                      -fuse_scphost => 'fusehost.ebi.ac.uk',
                      );
$fastq_getter->run;
my $output_file_list = $fastq_getter->output_files;

=cut


sub new {
    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);

    my ( $fuse_password, $fuse_user, $fuse_scphost, $fuse_scpuser)
        = rearrange( [ qw( FUSE_PASSWORD FUSE_USER FUSE_SCPHOST FUSE_SCPUSER)], @args);

    $self->fuse_password($fuse_password);
    $self->fuse_user($fuse_user);
    $self->fuse_scphost($fuse_scphost);
    $self->fuse_scpuser($fuse_scpuser || 'anonymous');

    return $self;
}

sub check_status {
  my $self = shift;
  my $era_rmia = $self->db->get_ERARunMetaInfoAdaptor;
  my $status = $era_rmia->get_status($self->run_meta_info->run_id);
  if ($status ne 'private' && $status ne 'public') {
    print "do not recognise status $status for ".$self->run_meta_info->run_id."\n";
    return 0;
  }
  return 1;
}

sub get_fastq_details {
  my $self = shift;
  throw("do not have a fuse user name") if !$self->fuse_user;
  throw("do not have a fuse password") if !$self->fuse_password;

  my $string = '/' . $self->fuse_user . '/' . $self->run_meta_info->study_id
            . '/' . $self->fuse_password . '/';
  my $digest = sha512_hex($string);
  my $root_dir = $self->source_root_dir;
  $root_dir .= '/' if ($root_dir !~ /\/$/);
  $root_dir .= $self->fuse_user .'/'. $self->run_meta_info->study_id .'/'. $digest;

  my ($db_md5_hash, $db_size_hash, $name_hash) =
      ReseqTrack::Tools::ERAUtils::get_fastq_details($self->run_meta_info->run_id, $self->db, $root_dir);
  $self->db_md5_hash($db_md5_hash);
  $self->db_size_hash($db_size_hash);
  $self->name_hash($name_hash);
}

sub get_files {
  my $self = shift;
  if (! $self->fuse_scphost) {
    warn("no fuse_scphost, so will get files using 'copy'");
    return $self->SUPER::get_files;
  }

  my $scp = Net::SCP->new($self->fuse_scphost, $self->fuse_scpuser);
  my $output_hash = $self->output_hash;
  while (my($source_path, $output_path) = each %$output_hash) {
    $scp->scp($source_path, $output_path)
        or throw("Failed scp from $source_path on ".$self->fuse_scphost." to $output_path: $scp->{errstr}");
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

sub fuse_scphost{
  my ($self, $fuse_scphost) = @_;
  if(defined($fuse_scphost)){
    $self->{fuse_scphost} = $fuse_scphost;
  }
  return $self->{fuse_scphost};
}

sub fuse_scpuser{
  my ($self, $fuse_scpuser) = @_;
  if(defined($fuse_scpuser)){
    $self->{fuse_scpuser} = $fuse_scpuser;
  }
  return $self->{fuse_scpuser};
}




1;

