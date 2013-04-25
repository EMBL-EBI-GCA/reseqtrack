
package ReseqTrack::HiveProcess::FileReleaseChecks;

use strict;

use base ('ReseqTrack::HiveProcess::BranchableProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception qw(throw);
use Digest::MD5::File qw(file_md5_hex);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

# project-specific classes can override this.
sub get_check_subs {
  return {quick => [\&check_size, \&check_name], slow => [\&check_md5]};
}

sub check_size {
  my ($self) = @_;
  my $db_size = $self->param('db_size');
  my $file_size = -s $self->param('dropbox_filename');
  return $db_size == $file_size ? 0 : 'sizes do not match';
}

sub check_md5 {
  my ($self) = @_;
  my $db_md5 = $self->param('db_md5');
  my $file_md5 = file_md5_hex($self->param('dropbox_filename'));
  return $db_md5 eq $file_md5 ? 0 : 'md5 does not match';
}

# this is project-specific and must be implented by a child class
sub check_name {
  return 0;
}

sub run {
    my $self = shift @_;

    my $db = ReseqTrack::DBSQL::DBAdaptor->new(%{$self->param('reseqtrack_db')});

    my $log_adaptor = $db->get_RejectLogAdaptor;

    my $check_subs = $self->get_check_subs();
    my $failed_checks = 0;
    CHECK:
    foreach my $sub (@{${$self->get_check_subs}{$self->param('check_class')}}) {
      if (my $reject_message =  &$sub($self)) {
        my $reject_log = ReseqTrack::RejectLog->new(
          -file_id => $self->param('db_file_id'),
          -is_reject => 'y',
          -reject_reason => $reject_message,
          );
        $log_adaptor->store($reject_log, 1);
        $failed_checks = 0;
        last CHECK;
      }
    }
    if (! $failed_checks) {
      $self->flows_this_branch(3);
    }
    return;
}

1;

