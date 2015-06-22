
package ReseqTrack::Hive::Process::FileRelease::Checks;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception qw(throw);
use Digest::MD5::File qw(file_md5_hex);
use File::stat;


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

# project-specific classes can override this.
sub get_check_subs {
  my ($self) = @_;
  return {quick => [map {$self->can($_)} qw(check_ctime check_update_time check_size check_name)],
          slow => [map {$self->can($_)} qw(check_ctime check_update_time check_md5)]};
}

sub check_size {
  my ($self, $dropbox_path, $file_object) = @_;
  my $db_size = $file_object->size;
  my $file_size = -s $dropbox_path;
  if ($db_size == 0) {
    $self->is_reject(0);
    $self->reject_message("file size is 0 in database");
    return 0;
  }
  if ($db_size < $file_size) {
    $self->is_reject(0);
    $self->reject_message("sizes do not match: $db_size $file_size");
    return 0;
  }
  if ($db_size > $file_size) {
    $self->is_reject(1);
    $self->reject_message("file size is greater than expected: $db_size $file_size");
    return 0;
  }

  my $file_details = $self->param('file');
  if (my $index_dbID = $file_details->{'db'}->{'index_dbID'}) {
    my $index_extension = $file_details->{'dropbox'}->{'index_ext'};
    my $index_file = $file_object->adaptor->fetch_by_dbID($index_dbID);
    throw("did not find file $index_dbID in database") if !$index_file;
    my $index_db_size = $index_file->size;
    my $index_file_size = -s $dropbox_path.$index_extension;
    if (!$index_file_size || $index_file_size != $index_db_size) {
      $self->is_reject(1);
      $self->reject_message("index file size error: $index_db_size $index_file_size");
      return 0;
    }
  }

  return 1;
}

sub check_md5 {
  my ($self, $dropbox_path, $file_object) = @_;
  my $db_md5 = $file_object->md5;
  my $file_md5 = file_md5_hex($dropbox_path);
  if ($file_md5 ne $db_md5) {
    $self->is_reject(1);
    $self->reject_message("md5 does not match: $file_md5");
    return 0;
  }

  my $file_details = $self->param('file');
  if (my $index_dbID = $file_details->{'db'}->{'index_dbID'}) {
    my $index_extension = $file_details->{'dropbox'}->{'index_ext'};
    my $index_file = $file_object->adaptor->fetch_by_dbID($index_dbID);
    throw("did not find file $index_dbID in database") if !$index_file;
    my $index_db_md5 = $index_file->md5;
    my $index_file_md5 = file_md5_hex($dropbox_path.$index_extension);
    if ($index_file_md5 ne $index_db_md5) {
      $self->is_reject(1);
      $self->reject_message("index md5 does not match: $index_file_md5");
      return 0;
    }
  }
  
  return 1;
}

sub check_ctime {
  my ($self, $dropbox_path, $file_object) = @_;
  my $st = stat($dropbox_path) or throw("could not stat $dropbox_path: $!");
  my $ctime_now = $st->ctime;
  my $file_details = $self->param('file');
  my $job_ctime = $file_details->{'dropbox'}->{'ctime'};
  if ($ctime_now != $job_ctime) {
    $self->is_reject(0);
    $self->reject_message("file changed since pipeline job created");
    return 0;
  }

  if (my $index_job_ctime = $file_details->{'dropbox'}->{'index_ctime'}) {
    my $index_path = $dropbox_path . $file_details->{'dropbox'}->{'index_ext'};
    my $index_st = stat($index_path) or throw("could not stat $index_path: $!");
    my $index_ctime_now = $index_st->ctime;
    if ($index_ctime_now != $index_job_ctime) {
      $self->is_reject(0);
      $self->reject_message("index file changed since pipeline job created");
      return 0;
    }
  }

  return 1;
}

sub check_update_time {
  my ($self, $dropbox_path, $file_object) = @_;
  my $db_updated = $file_object->updated;
  my $file_details = $self->param('file');
  my $job_updated = $file_details->{'db'}->{'updated'};
  if ($db_updated ne $job_updated) {
    $self->is_reject(0);
    $self->reject_message("checking job is out of date with database");
    return 0;
  }

  if (my $index_dbID = $file_details->{'db'}->{'index_dbID'}) {
    my $index_file = $file_object->adaptor->fetch_by_dbID($index_dbID);
    throw("did not find file $index_dbID in database") if !$index_file;
    my $index_db_updated = $index_file->updated;
    my $index_job_updated = $file_details->{'db'}->{'index_updated'};
    if ($index_db_updated ne $index_job_updated) {
      $self->is_reject(0);
      $self->reject_message("checking job is out of date with database");
      return 0;
    }
  }
  return 1;
}

# this is project-specific and must be implented by a child class
sub check_name {
  return 1;
}

sub is_reject {
  my ($self, $is_reject) = @_;
  if (defined $is_reject) {
    if (!$is_reject || $is_reject =~ /^n/i) {
      $self->param('is_reject', 0);
    }
    else {
      $self->param('is_reject', 1);
    }
  }
  return $self->param('is_reject');
}

sub reject_message {
  my ($self, $reject_message) = @_;
  if (defined $reject_message) {
    $self->param('reject_message', $reject_message);
  }
  return $self->param_is_defined('reject_message') ? $self->param('reject_message') : '';
}

sub param_defaults {
  return {
    'ps_attributes' => {},
    'is_failed' => undef,
  };
}


sub run {
    my $self = shift @_;

    my $file_details = $self->param_required('file');
    my $db_params = $self->param_required('reseqtrack_db');
    my $check_class = $self->param_required('check_class');
    my $ps_id = $self->param_required('ps_id');
    my $ps_attributes = $self->param('ps_attributes');
    
    throw("input_id not correctly formulated") if ! defined $file_details->{'dropbox'};
    throw("input_id not correctly formulated") if ! defined $file_details->{'db'};
    throw("input_id not correctly formulated") if ! defined $file_details->{'dropbox'}->{'path'};
    throw("input_id not correctly formulated") if ! defined $file_details->{'dropbox'}->{'ctime'};
    throw("input_id not correctly formulated") if ! defined $file_details->{'db'}->{'dbID'};
    throw("input_id not correctly formulated") if ! defined $file_details->{'db'}->{'updated'};

    my $db = ReseqTrack::DBSQL::DBAdaptor->new(%{$db_params});
    my $psa = $db->get_PipelineSeedAdaptor;
    my $pipeline_seed = $psa->fetch_by_dbID($ps_id);
    throw("did not get pipeline_seed for $ps_id") if !$pipeline_seed;

    # sanity check:
    my $self_dbname = $self->dbc->dbname;
    my $ps_dbname = $pipeline_seed->hive_db->name;
    throw("dbnames do not match $self_dbname $ps_dbname") if $self_dbname ne $ps_dbname;

    my $file_id = $file_details->{'db'}->{'dbID'};
    my $file_object = $db->get_FileAdaptor->fetch_by_dbID($file_id);
    throw("did not find file $file_id in database") if !$file_object;

    CHECK:
    foreach my $sub (@{${$self->get_check_subs}{$check_class}}) {
      $self->is_reject(0);
      $self->reject_message('');
      my $success = &$sub($self, $file_details->{'dropbox'}->{'path'}, $file_object);
      if (! $success) {
        my $attribute_name = $self->is_reject ? 'reject message' : 'message';
        $ps_attributes->{$attribute_name} = $self->reject_message;
        $self->output_param('ps_attributes', $ps_attributes);
        $self->output_param('is_failed', 1);
        last CHECK;
      }
    }
    return;
}

1;

