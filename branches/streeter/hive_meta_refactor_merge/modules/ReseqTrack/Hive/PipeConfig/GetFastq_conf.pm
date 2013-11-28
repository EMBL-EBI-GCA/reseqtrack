
package ReseqTrack::Hive::PipeConfig::GetFastq_conf;

use strict;
use warnings;

use base ('ReseqTrack::Hive::PipeConfig::ReseqTrackGeneric_conf');


=head2 default_options

    Description : Implements default_options() interface method of Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that is used to initialize default options.
                  In addition to the standard things it defines two options, 'first_mult' and 'second_mult' that are supposed to contain the long numbers to be multiplied.

=cut

sub default_options {
    my ($self) = @_;

    return {
        %{ $self->SUPER::default_options() },               # inherit other stuff from the base class

        'pipeline_name' => 'get_fastq',                     # name used by the beekeeper to prefix job names on the farm

        seeding_module => 'ReseqTrack::Hive::PipeSeed::Run',
        seeding_options => {
            output_columns => $self->o('run_columns'),
            output_attributes => $self->o('run_attributes'),
            require_columns => $self->o('require_run_columns'),
            exclude_columns => $self->o('exclude_run_columns'),
            require_attributes => $self->o('require_run_attributes'),
            exclude_attributes => $self->o('exclude_run_attributes'),
            output_sample_columns => $self->o('sample_columns'),
            output_sample_attributes => $self->o('sample_attribute_keys'),
            output_experiment_columns => $self->o('experiment_columns'),
            output_experiment_attributes => $self->o('experiment_attribute_keys'),
          },

        get_fastq_module => undef, # module default is ReseqTrack::Tools::GetFastq
        fastq_source_root_dir => undef, # module default is /nfs/era-pub
        clobber => 0,
        fastqc_exe => '/nfs/1000g-work/G1K/work/bin/FastQC/fastqc',

        fastq_type => undef, # use file_type_rule table
        fastqc_summary_type => undef, # use file_type_rule table
        fastqc_report_type => undef, # use file_type_rule table
        fastqc_zip_type => undef, # use file_type_rule table

        'sample_attribute_keys' => [],
        'sample_columns' => ['sample_source_id', 'sample_alias'],
        'run_attribute_keys' => [],
        'run_columns' => ['run_id', 'run_source_id'],
        'study_attribute_keys' => [],
        'study_columns' => ['study_source_id'],
        'experiment_attribute_keys' => [],
        'experiment_columns' => [],
        require_run_attributes => {},
        exclude_run_attributes => {},

        'fastq_output_dir' => $self->o('root_output_dir'),
        'fastq_output_layout' => '#sample_alias#/sequence_read',
        'fastqc_output_dir' => $self->o('fastq_output_dir'),
        'fastqc_output_layout' => '#fastq_output_layout#/fastqc',

        fastq_name_file_module => 'ReseqTrack::Hive::NameFile::BaseNameFile',
        fastq_name_file_method => 'basic',
        fastq_name_file_params => { new_dir => '#fastq_output_dir#/#fastq_output_layout#' },

        fastqc_name_file_module => 'ReseqTrack::Hive::NameFile::BaseNameFile',
        fastqc_name_file_method => 'basic',
        fastqc_name_file_params => { new_dir => '#fastqc_output_dir#/#fastqc_output_layout#' },

    };
}


sub resource_classes {
    my ($self) = @_;
    return {
            %{$self->SUPER::resource_classes},  # inherit 'default' from the parent class
            '200Mb' => { 'LSF' => '-C0 -M200 -q production -R"select[mem>200] rusage[mem=200]"' },
    };
}

sub pipeline_wide_parameters {
    my ($self) = @_;
    return {
        %{$self->SUPER::pipeline_wide_parameters},

        dir_label_params => ['study_source_id', 'sample_source_id', 'run_source_id'],
    };
}


sub pipeline_analyses {
    my ($self) = @_;

    my @analyses;
    push(@analyses, {
            -logic_name    => 'get_seeds',
            -module        => 'ReseqTrack::Hive::Process::SeedFactory',
            -meadow_type => 'LOCAL',
            -parameters    => {
                seeding_module => $self->o('seeding_module'),
                seeding_options => $self->o('seeding_options'),
            },
            -analysis_capacity  =>  1,  # use per-analysis limiter
            -flow_into => {
                2 => [ 'get_fastq' ],
            },
      });
    push(@analyses, {
            -logic_name    => 'get_fastq',
            -module        => 'ReseqTrack::Hive::Process::GetFastq',
            -parameters    => {
              module => $self->o('get_fastq_module'),
              source_root_dir => $self->o('fastq_source_root_dir'),
              clobber => $self->o('clobber'),
              era_dbuser => $self->o('era_dbuser'),
              era_dbpass => $self->o('era_dbpass'),
            },
            -flow_into => {
                1 => [ 'store_fastq' ],
            },
            -rc_name => '200Mb',
            -analysis_capacity  =>  25,  # use per-analysis limiter
            -hive_capacity  =>  -1,
      });
    push(@analyses, {
            -logic_name    => 'store_fastq',
            -module        => 'ReseqTrack::Hive::Process::LoadFile',
            -meadow_type => 'LOCAL',
            -parameters    => {
              type => $self->o('fastq_type'),
              collect => 1,
              clobber => $self->o('clobber'),
              name_file_module => $self->o('fastq_name_file_module'),
              name_file_method => $self->o('fastq_name_file_method'),
              name_file_params => $self->o('fastq_name_file_params'),
              fastq_output_dir => $self->o('fastq_output_dir'),
              fastq_output_layout => $self->o('fastq_output_layout'),
              md5 => '#fastq_md5#',
              collection_name => '#run_source_id#',
              file => '#fastq#',
            },
            -flow_into => {
                1 => { 'fastq_factory' => {'fastq' => '#file#'}},
            },
      });
    push(@analyses, {
            -logic_name    => 'fastq_factory',
            -module        => 'ReseqTrack::Hive::Process::JobFactory',
            -meadow_type => 'LOCAL',
            -parameters    => {
              factory_value => '#fastq#',
            },
            -flow_into => {
                '2->A' => { 'fastqc' => {'fastq' => '#factory_value#'}},
                'A->1' => [ 'store_fastqc_summary' ],
            },
      });
    push(@analyses, {
            -logic_name    => 'fastqc',
            -module        => 'ReseqTrack::Hive::Process::RunFastqc',
            -parameters    => {
              store_attributes => 1,
              program_file => $self->o('fastqc_exe'),
            },
            -flow_into => {
                1 => [ ':////accu?fastqc_summary=[]', ':////accu?fastqc_report=[]', ':////accu?fastqc_zip=[]']
            },
            -rc_name => '200Mb',
            -analysis_capacity  =>  50,  # use per-analysis limiter
            -hive_capacity  =>  -1,
      });
    push(@analyses, {
            -logic_name    => 'store_fastqc_summary',
            -module        => 'ReseqTrack::Hive::Process::LoadFile',
            -meadow_type => 'LOCAL',
            -parameters    => {
              type => $self->o('fastqc_summary_type'),
              collect => 1,
              name_file_module => $self->o('fastqc_name_file_module'),
              name_file_method => $self->o('fastqc_name_file_method'),
              name_file_params => $self->o('fastqc_name_file_params'),
              fastqc_output_dir => $self->o('fastqc_output_dir'),
              fastq_output_layout => $self->o('fastq_output_layout'),
              fastqc_output_layout => $self->o('fastqc_output_layout'),
              clobber => $self->o('clobber'),
              collection_name => '#run_source_id#',
              file => '#fastqc_summary#',
            },
            -flow_into => {
                1 => { 'store_fastqc_report' => {fastqc_summary => '#file#'}},
            },
      });
    push(@analyses, {
            -logic_name    => 'store_fastqc_report',
            -module        => 'ReseqTrack::Hive::Process::LoadFile',
            -meadow_type => 'LOCAL',
            -parameters    => {
              type => $self->o('fastqc_report_type'),
              collect => 1,
              name_file_module => $self->o('fastqc_name_file_module'),
              name_file_method => $self->o('fastqc_name_file_method'),
              name_file_params => $self->o('fastqc_name_file_params'),
              fastqc_output_dir => $self->o('fastqc_output_dir'),
              fastq_output_layout => $self->o('fastq_output_layout'),
              fastqc_output_layout => $self->o('fastqc_output_layout'),
              clobber => $self->o('clobber'),
              collection_name => '#run_source_id#',
              file => '#fastqc_report#',
            },
            -flow_into => {
                1 => { 'store_fastqc_zip' => {fastqc_report => '#file#'}},
            },
      });
    push(@analyses, {
            -logic_name    => 'store_fastqc_zip',
            -module        => 'ReseqTrack::Hive::Process::LoadFile',
            -meadow_type => 'LOCAL',
            -parameters    => {
              type => $self->o('fastqc_zip_type'),
              collect => 1,
              name_file_module => $self->o('fastqc_name_file_module'),
              name_file_method => $self->o('fastqc_name_file_method'),
              name_file_params => $self->o('fastqc_name_file_params'),
              fastqc_output_dir => $self->o('fastqc_output_dir'),
              fastq_output_layout => $self->o('fastq_output_layout'),
              fastqc_output_layout => $self->o('fastqc_output_layout'),
              clobber => $self->o('clobber'),
              collection_name => '#run_source_id#',
              file => '#fastqc_zip#',
            },
            -flow_into => {
                1 => { 'mark_seed_complete' => {fastqc_zip => '#file#'}},
            },
      });
    push(@analyses, {
            -logic_name    => 'mark_seed_complete',
            -module        => 'ReseqTrack::Hive::Process::UpdateSeed',
            -parameters    => {
              is_complete  => 1,
            },
            -meadow_type => 'LOCAL',
      });

    return \@analyses;
}

1;

