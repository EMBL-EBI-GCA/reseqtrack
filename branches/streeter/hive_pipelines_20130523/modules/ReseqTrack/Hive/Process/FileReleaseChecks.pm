
package ReseqTrack::Hive::Process::FileReleaseChecks;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
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
  my ($self, $is_reject_ref, $reject_message_ref) = @_;
  my $db_size = $self->param('db_size');
  my $file_size = -s $self->param('dropbox_filename');
  if ($db_size != $file_size) {
    $$is_reject_ref = 'n';
    $$reject_message_ref = "sizes do not match: $db_size $file_size";
    return 0;
  }
  return 1;
}

sub check_md5 {
  my ($self, $is_reject_ref, $reject_message_ref) = @_;
  my $db_md5 = $self->param('db_md5');
  my $file_md5 = file_md5_hex($self->param('dropbox_filename'));
  if ($file_md5 ne $db_md5) {
    $$is_reject_ref = 'y';
    $$reject_message_ref = "md5 does not match: $file_md5";
    return 0;
  }
  return 1;
}

# this is project-specific and must be implented by a child class
sub check_name {
  return 1;
}

sub run {
    my $self = shift @_;

    my $db = ReseqTrack::DBSQL::DBAdaptor->new(%{$self->param('reseqtrack_db')});

    my $log_adaptor = $db->get_RejectLogAdaptor;

    my $check_subs = $self->get_check_subs();
    my ($is_reject, $reject_message);
    CHECK:
    foreach my $sub (@{${$self->get_check_subs}{$self->param('check_class')}}) {
      if (! &$sub($self, \$is_reject, \$reject_message)) {
        my $reject_log = ReseqTrack::RejectLog->new(
          -file_id => $self->param('db_file_id'),
          -is_reject => $is_reject,
          -reject_reason => $reject_message,
          );
        $log_adaptor->store($reject_log, 1);
        $self->input_job->autoflow(0);
        last CHECK;
      }
    }
    return;
}

1;

