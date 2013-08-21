
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


        get_fastq_module => undef, # module default is ReseqTrack::Tools::GetFastq
        fastq_source_root_dir => undef, # module default is /nfs/era-pub
        clobber => 0,
        fastqc_output_dir => $self->o('fastq_output_dir'),
        fastqc_exe => '/nfs/1000g-work/G1K/work/bin/FastQC/fastqc',

        fastq_type => undef, # use file_type_rule table
        fastqc_summary_type => undef, # use file_type_rule table
        fastqc_report_type => undef, # use file_type_rule table
        fastqc_zip_type => undef, # use file_type_rule table

        use_label_management => 0,

    };
}


sub resource_classes {
    my ($self) = @_;
    return {
            %{$self->SUPER::resource_classes},  # inherit 'default' from the parent class
            '200Mb' => { 'LSF' => '-C0 -M200 -q production -R"select[mem>200] rusage[mem=200]"' },
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
              seed_label => 'source_id',
              output_columns => 'source_id',
              output_sample_columns => 'sample_alias',
              output_sample_attributes => 'population',

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
              fastq_output_dir => $self->o('fastq_output_dir'),
              output_dir => q(#expr(join('/', $fastq_output_dir, $population, $sample_alias, 'sequence'))expr#),
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
              collection_name => '#source_id#',
              file => '#fastq#',
            },
            -flow_into => {
                1 => [ 'fastq_factory' ],
            },
      });
    push(@analyses, {
            -logic_name    => 'fastq_factory',
            -module        => 'ReseqTrack::Hive::Process::JobFactory',
            -meadow_type => 'LOCAL',
            -parameters    => {
              factory_value => '#fastq#',
              temp_param_sub => { 2 => [['fastq','factory_value']]}, # temporary hack pending updates to hive code
            },
            -flow_into => {
                '2->A' => [ 'fastqc' ],
                'A->1' => [ 'store_fastqc_summary' ],
            },
      });
    push(@analyses, {
            -logic_name    => 'fastqc',
            -module        => 'ReseqTrack::Hive::Process::RunFastqc',
            -parameters    => {
              store_attributes => 1,
              program_file => $self->o('fastqc_exe'),
              fastqc_output_dir => $self->o('fastqc_output_dir'),
              output_dir => q(#expr(join('/', $fastqc_output_dir, $population, $sample_alias, 'sequence', 'fastqc'))expr#),
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
              collection_name => '#source_id#',
              file => '#fastqc_summary#',
            },
            -flow_into => {
                1 => [ 'store_fastqc_report' ],
            },
      });
    push(@analyses, {
            -logic_name    => 'store_fastqc_report',
            -module        => 'ReseqTrack::Hive::Process::LoadFile',
            -meadow_type => 'LOCAL',
            -parameters    => {
              type => $self->o('fastqc_report_type'),
              collect => 1,
              collection_name => '#source_id#',
              file => '#fastqc_report#',
            },
            -flow_into => {
                1 => [ 'store_fastqc_zip' ],
            },
      });
    push(@analyses, {
            -logic_name    => 'store_fastqc_zip',
            -module        => 'ReseqTrack::Hive::Process::LoadFile',
            -meadow_type => 'LOCAL',
            -parameters    => {
              type => $self->o('fastqc_zip_type'),
              collect => 1,
              collection_name => '#source_id#',
              file => '#fastqc_zip#',
            },
            -flow_into => {
                1 => [ 'mark_seed_complete' ],
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

