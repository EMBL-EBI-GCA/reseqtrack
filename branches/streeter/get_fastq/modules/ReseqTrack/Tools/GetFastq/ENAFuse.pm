package ReseqTrack::Tools::GetFastq::ENAFuse;

use strict;
use warnings;

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
                      );
$fastq_getter->run;
my $output_file_list = $fastq_getter->output_files;

=cut


sub new {
    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);

    my ( $fuse_password, $fuse_user)
        = rearrange( [ qw( FUSE_PASSWORD FUSE_USER )], @args);

    $self->fuse_password($fuse_password);
    $self->fuse_user($fuse_user);

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

sub make_output_hash {
  my $self = shift;
  throw("do not have a fuse user name") if !$self->fuse_user;
  throw("do not have a fuse password") if !$self->fuse_password;

  my $output_dir = $self->output_dir;
  $output_dir .= '/' . $self->fuse_user;
  $output_dir .= '/' . $self->run_meta_info->study_id;
  $output_dir .= '/' . $self->fuse_user . $self->run_meta_info->study_id . $self->fuse_password;
  my $name_hash = $self->name_hash;
  my %output_hash;
  while (my ($basename, $source_path) = each %$name_hash) {
    my $output_path = $output_dir . '/' . $basename;
    $output_path =~ s{//}{/}g;
    $output_hash{$source_path} = $output_path;
  }
  $self->output_hash(\%output_hash);
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


1;
