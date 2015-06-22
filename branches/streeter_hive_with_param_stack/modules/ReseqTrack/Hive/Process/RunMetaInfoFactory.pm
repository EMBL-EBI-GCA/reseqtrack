
package ReseqTrack::Hive::Process::RunMetaInfoFactory;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception qw(throw);

sub param_defaults {
  return {
    allowed_platform => ['ILLUMINA'],
    allowed_status => ['public', 'private'],
  };
}


sub libraries_factory {
    my ($self) = @_;
    my $sample_id = $self->param_required('sample_id');
    my $allowed_strategies = $self->param_to_flat_array('allowed_strategy');
    my $allowed_platforms = $self->param_to_flat_array('allowed_platform');

    my $db = ReseqTrack::DBSQL::DBAdaptor->new(%{$self->param('reseqtrack_db')});
    my %library_names;
    EXPERIMENT:
    foreach my $experiment (@{$db->get_ExperimentAdaptor->fetch_by_sample_id($sample_id)}) {
      if (@$allowed_strategies) {
        next EXPERIMENT if ! grep {$experiment->library_strategy eq $_} @$allowed_strategies;
      }
      next EXPERIMENT if ! grep {$experiment->instrument_platform eq $_} @$allowed_platforms;
      $library_names{$experiment->library_name} = 1;
    }
    foreach my $library_name (keys %library_names) {
      $self->prepare_factory_output_id({'library_name' => $library_name});
    }
}

sub runs_factory {
    my ($self) = @_;
    my $library_name = $self->param_required('library_name');
    my $sample_id = $self->param_required('sample_id');
    my $allowed_statuses = $self->param_to_flat_array('allowed_status');

    my $db = ReseqTrack::DBSQL::DBAdaptor->new(%{$self->param('reseqtrack_db')});

    my $experiments = $db->get_ExperimentAdaptor->fetch_by_column_name('library_name', $library_name);
    my $ra = $db->get_RunAdaptor;
    foreach my $experiment (@{$db->get_ExperimentAdaptor->fetch_by_column_name('library_name', $library_name)}) {
      RUN:
      foreach my $run (@{$ra->fetch_by_experiment_id($experiment->dbID)}) {
        next RUN if $run->sample_id != $sample_id;
        next RUN if ! grep {$run->status eq $_} @$allowed_statuses;
        $self->prepare_factory_output_id({'run_id' => $run->dbID, 'run_source_id' => $run->source_id});
      }
    }
}

sub run {
    my $self = shift @_;

    my %factories = (
        'run' => \&runs_factory,
        'library' => \&libraries_factory,
        );

    &{$factories{$self->param_required('factory_type')}}($self);
}

1;

