
package ReseqTrack::Hive::Process::RunMetaInfoFactory;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception qw(throw);

sub param_defaults {
  return {
    output_experiment_columns => [],
    output_experiment_attributes => [],
    output_run_columns => [],
    output_run_attributes => [],
    output_study_columns => [],
    output_study_attributes => [],

    require_experiment_columns => {
          instrument_platform => ['ILLUMINA'], 
        },
    require_run_columns => {
          allowed_status => ['public', 'private'],
        },
    require_study_columns => {},

    require_run_attributes => {},
    require_experiment_attributes => {},
    require_study_attributes => {},
    exclude_run_attributes => {},
    exclude_experiment_attributes => {},
    exclude_study_attributes => {},

  };
}


sub libraries_factory {
    my ($self) = @_;
    my $sample_id = $self->param_required('sample_id');
    my $require_experiment_columns = $self->param('require_experiment_columns') || param_defaults()->{'require_experiment_columns'};
    my $require_experiment_attributes = $self->param('require_experiment_attributes') || param_defaults()->{'require_experiment_attributes'};
    my $exclude_experiment_attributes = $self->param('exclude_experiment_attributes') || param_defaults()->{'exclude_experiment_attributes'};
    my $require_study_columns = $self->param('require_study_columns') || param_defaults()->{'require_study_columns'};
    my $require_study_attributes = $self->param('require_study_attributes') || param_defaults()->{'require_study_attributes'};
    my $exclude_study_attributes = $self->param('exclude_study_attributes') || param_defaults()->{'exclude_study_attributes'};

    my $db = ReseqTrack::DBSQL::DBAdaptor->new(%{$self->param('reseqtrack_db')});
    my $sta = $db->get_StudyAdaptor;
    my $ea = $db->get_ExperimentAdaptor;
    my %library_names;
    EXP:
    foreach my $experiment (@{$ea->fetch_by_sample_id($sample_id)}) {

      foreach my $column_name (keys %$require_experiment_columns) {
        my $required = $require_experiment_columns->{$column_name};
        $required = [$required] if !ref($required);
        my $val = &{$ea->column_mappings($experiment)->{$column_name}}();
        next EXP if ! grep {$val eq $_} @$required;
      }
      foreach my $attr_name (keys %$require_experiment_attributes) {
        my ($attribute) = grep {$_->attribute_name eq $attr_name} @{$experiment->attributes};
        next EXP if !$attribute;
        my $required = $require_experiment_attributes->{$attr_name};
        $required = [$required] if !ref($required);
        next EXP if !grep {$attribute->attribute_value eq $_} @$required;
      }
      ATTR:
      foreach my $attr_name (keys %$exclude_experiment_attributes) {
        my ($attribute) = grep {$_->attribute_name eq $attr_name} @{$experiment->attributes};
        next ATTR if !$attribute;
        my $excluded = $exclude_experiment_attributes->{$attr_name};
        $excluded = [$excluded] if !ref($excluded);
        next EXP if grep {$attribute->attribute_value eq $_} @$excluded;
      }

      my $study = $sta->fetch_by_dbID($experiment->study_id);
      foreach my $column_name (keys %$require_study_columns) {
        my $required = $require_study_columns->{$column_name};
        $required = [$required] if !ref($required);
        my $val = &{$sta->column_mappings($study)->{$column_name}}();
        next EXP if ! grep {$val eq $_} @$required;
      }
      foreach my $attr_name (keys %$require_study_attributes) {
        my ($attribute) = grep {$_->attribute_name eq $attr_name} @{$study->attributes};
        next EXP if !$attribute;
        my $required = $require_study_attributes->{$attr_name};
        $required = [$required] if !ref($required);
        next EXP if !grep {$attribute->attribute_value eq $_} @$required;
      }
      ATTR:
      foreach my $attr_name (keys %$exclude_study_attributes) {
        my ($attribute) = grep {$_->attribute_name eq $attr_name} @{$study->attributes};
        next ATTR if !$attribute;
        my $excluded = $exclude_study_attributes->{$attr_name};
        $excluded = [$excluded] if !ref($excluded);
        next EXP if grep {$attribute->attribute_value eq $_} @$excluded;
      }

      $library_names{$experiment->library_name} = 1;
    }
    foreach my $library_name (keys %library_names) {
      $self->prepare_factory_output_id({'library_name' => $library_name});
    }
}



sub runs_factory {
    my ($self) = @_;
    my $library_name = $self->param('library_name');
    my $experiment_source_id = $self->param('experiment_source_id');
    throw('need library name or experiment source id')
           if ( !$library_name && !$experiment_source_id);  
           
    my $sample_id = $self->param_required('sample_id');
    my $output_run_columns = $self->param('output_run_columns') || [];
    my $output_run_attributes = $self->param('output_run_attributes') || [];
    my $output_experiment_columns = $self->param('output_experiment_columns') || [];
    my $output_experiment_attributes = $self->param('output_experiment_attributes') || [];
    my $output_study_columns = $self->param('output_study_columns') || [];
    my $output_study_attributes = $self->param('output_study_attributes') || [];

    my $require_run_columns = $self->param('require_run_columns') || param_defaults()->{'require_run_columns'};
    my $require_run_attributes = $self->param('require_run_attributes') || param_defaults()->{'require_run_attributes'};
    my $exclude_run_attributes = $self->param('exclude_run_attributes') || param_defaults()->{'exclude_run_attributes'};
    my $require_experiment_columns = $self->param('require_experiment_columns') || param_defaults()->{'require_experiment_columns'};
    my $require_experiment_attributes = $self->param('require_experiment_attributes') || param_defaults()->{'require_experiment_attributes'};
    my $exclude_experiment_attributes = $self->param('exclude_experiment_attributes') || param_defaults()->{'exclude_experiment_attributes'};
    my $require_study_columns = $self->param('require_study_columns') || param_defaults()->{'require_study_columns'};
    my $require_study_attributes = $self->param('require_study_attributes') || param_defaults()->{'require_study_attributes'};
    my $exclude_study_attributes = $self->param('exclude_study_attributes') || param_defaults()->{'exclude_study_attributes'};

    my $db = ReseqTrack::DBSQL::DBAdaptor->new(%{$self->param('reseqtrack_db')});

    my $ra = $db->get_RunAdaptor;
    my $ea = $db->get_ExperimentAdaptor;
    my $sta = $db->get_StudyAdaptor;
    
    my $experiments;
    if ( $library_name ) {
      $experiments = $ea->fetch_by_column_name('library_name', $library_name);
    }
    elsif ( $experiment_source_id && !$library_name ) {
      my $experiment = $ea->fetch_by_source_id( $experiment_source_id );  
      push @{$experiments}, $experiment;
    }
    
    EXP:
    foreach my $experiment (@{$experiments}) {
      foreach my $column_name (keys %$require_experiment_columns) {
        my $required = $require_experiment_columns->{$column_name};
        $required = [$required] if !ref($required);
        my $val = &{$ea->column_mappings($experiment)->{$column_name}}();
        next EXP if ! grep {$val eq $_} @$required;
      }
      foreach my $attr_name (keys %$require_experiment_attributes) {
        my ($attribute) = grep {$_->attribute_name eq $attr_name} @{$experiment->attributes};
        next EXP if !$attribute;
        my $required = $require_experiment_attributes->{$attr_name};
        $required = [$required] if !ref($required);
        next EXP if !grep {$attribute->attribute_value eq $_} @$required;
      }
      ATTR:
      foreach my $attr_name (keys %$exclude_experiment_attributes) {
        my ($attribute) = grep {$_->attribute_name eq $attr_name} @{$experiment->attributes};
        next ATTR if !$attribute;
        my $excluded = $exclude_experiment_attributes->{$attr_name};
        $excluded = [$excluded] if !ref($excluded);
        next EXP if grep {$attribute->attribute_value eq $_} @$excluded;
      }

      my $study = $sta->fetch_by_dbID($experiment->study_id);
      foreach my $column_name (keys %$require_study_columns) {
        my $required = $require_study_columns->{$column_name};
        $required = [$required] if !ref($required);
        my $val = &{$sta->column_mappings($study)->{$column_name}}();
        next EXP if ! grep {$val eq $_} @$required;
      }
      foreach my $attr_name (keys %$require_study_attributes) {
        my ($attribute) = grep {$_->attribute_name eq $attr_name} @{$study->attributes};
        next EXP if !$attribute;
        my $required = $require_study_attributes->{$attr_name};
        $required = [$required] if !ref($required);
        next EXP if !grep {$attribute->attribute_value eq $_} @$required;
      }
      ATTR:
      foreach my $attr_name (keys %$exclude_study_attributes) {
        my ($attribute) = grep {$_->attribute_name eq $attr_name} @{$study->attributes};
        next ATTR if !$attribute;
        my $excluded = $exclude_study_attributes->{$attr_name};
        $excluded = [$excluded] if !ref($excluded);
        next EXP if grep {$attribute->attribute_value eq $_} @$excluded;
      }



      RUN:
      foreach my $run (@{$ra->fetch_by_experiment_id($experiment->dbID)}) {
        next RUN if $run->sample_id != $sample_id;

        foreach my $column_name (keys %$require_run_columns) {
          my $required = $require_run_columns->{$column_name};
          $required = [$required] if !ref($required);
          my $val = &{$ra->column_mappings($run)->{$column_name}}();
          next RUN if ! grep {$val eq $_} @$required;
        }
        foreach my $attr_name (keys %$require_run_attributes) {
          my ($attribute) = grep {$_->attribute_name eq $attr_name} @{$run->attributes};
          next RUN if !$attribute;
          my $required = $require_run_attributes->{$attr_name};
          $required = [$required] if !ref($required);
          next RUN if !grep {$attribute->attribute_value eq $_} @$required;
        }
        ATTR:
        foreach my $attr_name (keys %$exclude_run_attributes) {
          my ($attribute) = grep {$_->attribute_name eq $attr_name} @{$run->attributes};
          next ATTR if !$attribute;
          my $excluded = $exclude_run_attributes->{$attr_name};
          $excluded = [$excluded] if !ref($excluded);
          next RUN if grep {$attribute->attribute_value eq $_} @$excluded;
        }


        my %output_hash = (run_id => $run->dbID);

        foreach my $column_name (@$output_run_columns) {
          $output_hash{$column_name} = &{$ra->column_mappings($run)->{$column_name}}();
        }
        ATTRIBUTE:
        foreach my $attr_name (@$output_run_attributes) {
          my ($attribute) = grep {$_->attribute_name eq $attr_name} @{$run->attributes};
          next ATTRIBUTE if !$attribute;
          $output_hash{$attr_name} = $attribute->attribute_value;
        }

        foreach my $column_name (@$output_experiment_columns) {
          $output_hash{$column_name} = &{$ea->column_mappings($experiment)->{$column_name}}();
        }
        ATTRIBUTE:
        foreach my $attr_name (@$output_experiment_attributes) {
          my ($attribute) = grep {$_->attribute_name eq $attr_name} @{$experiment->attributes};
          next ATTRIBUTE if !$attribute;
          $output_hash{$attr_name} = $attribute->attribute_value;
        }

        foreach my $column_name (@$output_study_columns) {
          $output_hash{$column_name} = &{$sta->column_mappings($study)->{$column_name}}();
        }
        ATTRIBUTE:
        foreach my $attr_name (@$output_study_attributes) {
          my ($attribute) = grep {$_->attribute_name eq $attr_name} @{$study->attributes};
          next ATTRIBUTE if !$attribute;
          $output_hash{$attr_name} = $attribute->attribute_value;
        }

        $self->prepare_factory_output_id(\%output_hash);
      }
    }
    $db->dbc->disconnect_when_inactive(1);
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

