
package ReseqTrack::HiveProcess::RunMetaInfoFactory;

use strict;

use base ('ReseqTrack::HiveProcess::BranchableProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception qw(throw);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;

    my $type_branch = $self->param('type_branch') || die "'type_branch' is an obligatory parameter";
    my $output_dir = $self->param('output_dir') || die "'output_dir' is an obligatory parameter";

    my $db = ReseqTrack::DBSQL::DBAdaptor->new(%{$self->param('reseqtrack_db')});

    if (lc($type_branch) eq 'sample') {
      #my $sql = "SELECT DISTINCT sample_id FROM run_meta_info";
      my $sql = "SELECT sample_id, sample_name FROM run_meta_info";
      my @bind_values;
      if (my $study_id = $self->param('branch_study_id')) {
        $sql .= ' WHERE study_id = ?';
        push(@bind_values, $study_id);
      }
      $sql .= " GROUP BY sample_id";
      my $sth = $db->dbc->prepare($sql) or die "could not prepare $sql: ".$db->dbc->errstr;;
      $sth->execute(@bind_values) or die "could not execute $sql: ".$sth->errstr;

      while (my $row = $sth->fetchrow_arrayref) {
        my ($sample_id, $sample_name) = @$row;
        $self->output_child_branches('SAMPLE_ID' => $sample_id, 'SAMPLE_NAME' => $sample_name, 'label' => $sample_name, 'output_dir' => "$output_dir/$sample_id");
      }
      #foreach my $sample_id (map {$_->[0]} @{$sth->fetchall_arrayref()}) {
        #$self->output_child_branches('SAMPLE_ID' => $sample_id, 'label' => $sample_id, 'output_dir' => "$output_dir/$sample_id");
      #}
    }
    elsif(lc($type_branch) eq 'library') {
      my $sql = 'SELECT DISTINCT library_name FROM run_meta_info WHERE 1';
      my @bind_values;
      if (my $branch_sample_id = $self->param('branch_sample_id')) {
        $sql .= ' AND sample_id = ?';
        push(@bind_values, $branch_sample_id);
      }
      if (my $study_id = $self->param('branch_study_id')) {
        $sql .= ' AND study_id = ?';
        push(@bind_values, $study_id);
      }
      my $sth = $db->dbc->prepare($sql) or die "could not prepare $sql: ".$db->dbc->errstr;;
      $sth->execute(@bind_values) or die "could not execute $sql: ".$sth->errstr;
      foreach my $library_name (map {$_->[0]} @{$sth->fetchall_arrayref()}) {
        $self->output_child_branches('LIBRARY_NAME' => $library_name, 'label' => $library_name, 'output_dir' => "$output_dir/$library_name");
      }
    }
    elsif(lc($type_branch) eq 'run') {
      my $sql = 'SELECT run_id FROM run_meta_info where 1';
      my @bind_values;
      if (my $branch_sample_id = $self->param('branch_sample_id')) {
        $sql .= ' AND sample_id = ?';
        push(@bind_values, $branch_sample_id);
      }
      if (my $library_name = $self->param('branch_library_name')) {
        $sql .= ' AND library_name = ?';
        push(@bind_values, $library_name);
      }
      if (my $study_id = $self->param('branch_study_id')) {
        $sql .= ' AND study_id = ?';
        push(@bind_values, $study_id);
      }
      my $sth = $db->dbc->prepare($sql) or die "could not prepare $sql: ".$db->dbc->errstr;;
      $sth->execute(@bind_values) or die "could not execute $sql: ".$sth->errstr;
      foreach my $run_id (map {$_->[0]} @{$sth->fetchall_arrayref()}) {
        $self->output_child_branches('RUN_ID' => $run_id, 'label' => $run_id, 'output_dir' => "$output_dir/$run_id");
      }
    }

}

1;

