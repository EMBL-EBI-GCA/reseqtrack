
package ReseqTrack::Hive::Process::RunMetaInfoFactory;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception qw(throw);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;

    my $type_branch = $self->param_required('type_branch');
    my $output_dir = $self->output_dir;

    my $allowed_status_arr = $self->param_is_defined('allowed_status') ? $self->get_param_values('allowed_status')
                            : ['public', 'private'];
    my $allowed_platform_arr = $self->param_is_defined('allowed_platform') ? $self->get_param_values('allowed_platform')
                            : ['ILLUMINA'];

    my $db = ReseqTrack::DBSQL::DBAdaptor->new(%{$self->param('reseqtrack_db')});

    if (lc($type_branch) eq 'sample') {
      my $sql = "SELECT sample_id, sample_name FROM run_meta_info";
      my @bind_values;
      if ($self->param_is_defined('study_id')) {
        $sql .= ' WHERE study_id = ?';
        push(@bind_values, $self->param('study_id'));
      }
      $sql .= " GROUP BY sample_id";
      my $sth = $db->dbc->prepare($sql) or die "could not prepare $sql: ".$db->dbc->errstr;;
      $sth->execute(@bind_values) or die "could not execute $sql: ".$sth->errstr;

      while (my $row = $sth->fetchrow_arrayref) {
        my ($sample_id, $sample_name) = @$row;
        $self->prepare_factory_output_id($sample_name, {'sample_id' => $sample_id});
      }
    }
    elsif(lc($type_branch) eq 'library') {
      my $sql = 'SELECT DISTINCT library_name FROM run_meta_info WHERE 1';
      my @bind_values;
      if ($self->param_is_defined('sample_id')) {
        $sql .= ' AND sample_id = ?';
        push(@bind_values, $self->param('sample_id'));
      }
      if ($self->param_is_defined('study_id')) {
        $sql .= ' AND study_id = ?';
        push(@bind_values, $self->param('study_id'));
      }
      my $sth = $db->dbc->prepare($sql) or die "could not prepare $sql: ".$db->dbc->errstr;;
      $sth->execute(@bind_values) or die "could not execute $sql: ".$sth->errstr;
      foreach my $library_name (map {$_->[0]} @{$sth->fetchall_arrayref()}) {
        $self->prepare_factory_output_id($library_name, {'library_name' => $library_name});
      }
    }
    elsif(lc($type_branch) eq 'run') {
      my $sql = 'SELECT run_id FROM run_meta_info where 1';
      my @bind_values;
      if ($self->param_is_defined('sample_id')) {
        $sql .= ' AND sample_id = ?';
        push(@bind_values, $self->param('sample_id'));
      }
      if ($self->param_is_defined('library_name')) {
        $sql .= ' AND library_name = ?';
        push(@bind_values, $self->param('library_name'));
      }
      if ($self->param_is_defined('study_id')) {
        $sql .= ' AND study_id = ?';
        push(@bind_values, $self->param('study_id'));
      }
      if (@$allowed_status_arr) {
        $sql .= ' AND status in (' . join(',', map {'?'} @$allowed_status_arr) . ')';
        push(@bind_values, @$allowed_status_arr);
      }
      if (@$allowed_platform_arr) {
        $sql .= ' AND instrument_platform in (' . join(',', map {'?'} @$allowed_platform_arr) . ')';
        push(@bind_values, @$allowed_platform_arr);
      }
      my $sth = $db->dbc->prepare($sql) or die "could not prepare $sql: ".$db->dbc->errstr;;
      $sth->execute(@bind_values) or die "could not execute $sql: ".$sth->errstr;
      foreach my $run_id (map {$_->[0]} @{$sth->fetchall_arrayref()}) {
        $self->prepare_factory_output_id($run_id, {'run_id' => $run_id});
      }
    }

}

1;

