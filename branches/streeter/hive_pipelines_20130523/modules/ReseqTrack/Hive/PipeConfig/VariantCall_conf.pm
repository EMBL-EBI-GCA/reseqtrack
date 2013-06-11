

package ReseqTrack::Hive::PipeConfig::VariantCall_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');  # All Hive databases configuration files should inherit from HiveGeneric, directly or indirectly


=head2 default_options

    Description : Implements default_options() interface method of Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that is used to initialize default options.
                  In addition to the standard things it defines two options, 'first_mult' and 'second_mult' that are supposed to contain the long numbers to be multiplied.

=cut

sub default_options {
    my ($self) = @_;

    return {
        %{ $self->SUPER::default_options() },               # inherit other stuff from the base class

        'pipeline_name' => 'vc',                     # name used by the beekeeper to prefix job names on the farm

        'reseqtrack_db'  => {
            -host => $self->o('host'),
            -port => 4197,
            -user => 'g1kro',
            -pass => '',
            #-dbname => 'nextgen_track', # set on the command line
        },

        'type_bam'    => 'BAM',
        'transpose_exe' => $self->o('ENV', 'RESEQTRACK').'/c_code/transpose_bam/transpose_bam',
        'samtools_exe' => '/nfs/1000g-work/G1K/work/bin/samtools/samtools',
        'bcftools_exe' => '/nfs/1000g-work/G1K/work/bin/samtools/bcftools/bcftools',
        'vcfutils_exe' => '/nfs/1000g-work/G1K/work/bin/samtools/bcftools/vcfutils.pl',
        'bgzip_exe' => '/nfs/1000g-work/G1K/work/bin/tabix/bgzip',
        'gatk_dir' => '/nfs/1000g-work/G1K/work/bin/gatk/dist/',
        'freebayes_exe' => '/nfs/1000g-work/G1K/work/bin/freebayes/bin/freebayes',

        call_by_gatk_options => {}, # use module defaults
        call_by_samtools_options => {}, # use module defaults
        call_by_freebayes_options => {}, # use module defaults

        'fai' => $self->o('reference') . '.fai',
        'target_bed_file' => undef,

        'final_label' => $self->o('pipeline_name'),
        'call_window_size' => 50000,
        'transpose_window_size' => 50000000,

        'call_by_samtools' => 1,
        'call_by_gatk' => 1,
        'call_by_freebayes' => 1,


    };
}


=head2 pipeline_create_commands

    Description : Implements pipeline_create_commands() interface method of Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that lists the commands that will create and set up the Hive database.
                  In addition to the standard creation of the database and populating it with Hive tables and procedures it also creates two pipeline-specific tables used by Runnables to communicate.

=cut

sub pipeline_create_commands {
    my ($self) = @_;

    my $sql_1 = '
    CREATE TABLE reseqtrack_file (
      file_id int(10) unsigned NOT NULL AUTO_INCREMENT,
      name VARCHAR(64000) NOT NULL,
      PRIMARY KEY (file_id)
    )';

    return [
        @{$self->SUPER::pipeline_create_commands},  # inheriting database and hive tables' creation

        $self->db_execute_command('pipeline_db', $sql_1),

    ];
}


=head2 pipeline_wide_parameters

    Description : Interface method that should return a hash of pipeline_wide_parameter_name->pipeline_wide_parameter_value pairs.
                  The value doesn't have to be a scalar, can be any Perl structure now (will be stringified and de-stringified automagically).
                  Please see existing PipeConfig modules for examples.

=cut

sub pipeline_wide_parameters {
    my ($self) = @_;
    return {
        %{$self->SUPER::pipeline_wide_parameters},          # here we inherit anything from the base class

        'root_output_dir' => $self->o('root_output_dir'),
        'reseqtrack_db' => $self->o('reseqtrack_db'),
        fai => $self->o('fai'),

        'file_param_names' => ['bam', 'bai', 'vcf'],

    };
}

sub resource_classes {
    my ($self) = @_;
    return {
            %{$self->SUPER::resource_classes},  # inherit 'default' from the parent class
            '200Mb' => { 'LSF' => '-C0 -M200 -q production -R"select[mem>200] rusage[mem=200]"' },
            '500Mb' => { 'LSF' => '-C0 -M500 -q production -R"select[mem>500] rusage[mem=500]"' },
            '1Gb'   => { 'LSF' => '-C0 -M1000 -q production -R"select[mem>1000] rusage[mem=1000]"' },
            '2Gb' => { 'LSF' => '-C0 -M2000 -q production -R"select[mem>2000] rusage[mem=2000]"' },
            '4Gb' => { 'LSF' => '-C0 -M4000 -q production -R"select[mem>4000] rusage[mem=4000]"' },
            '5Gb' => { 'LSF' => '-C0 -M5000 -q production -R"select[mem>5000] rusage[mem=5000]"' },
    };
}


=head2 pipeline_analyses

    Description : Implements pipeline_analyses() interface method of Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that defines the structure of the pipeline: analyses, jobs, rules, etc.
                  Here it defines three analyses:

                    * 'start' with two jobs (multiply 'first_mult' by 'second_mult' and vice versa - to check the commutativity of multiplivation).
                      Each job will dataflow (create more jobs) via branch #2 into 'part_multiply' and via branch #1 into 'add_together'.

                    * 'part_multiply' initially without jobs (they will flow from 'start')

                    * 'add_together' initially without jobs (they will flow from 'start').
                       All 'add_together' jobs will wait for completion of 'part_multiply' jobs before their own execution (to ensure all data is available).

    There are two control modes in this pipeline:
        A. The default mode is to use the '2' and '1' dataflow rules from 'start' analysis and a -wait_for rule in 'add_together' analysis for analysis-wide synchronization.
        B. The semaphored mode is to use '2->A' and 'A->1' semaphored dataflow rules from 'start' instead, and comment out the analysis-wide -wait_for rule, relying on semaphores.

=cut

sub pipeline_analyses {
    my ($self) = @_;

    my @analyses;
    push(@analyses, {
            -logic_name    => 'callgroups_factory',
            -module        => 'ReseqTrack::Hive::Process::JobFactory',
            -meadow_type => 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
            -parameters    => {
                factory_value => '#callgroup#',
                temp_param_sub => { 2 => [['callgroup','factory_value']]}, # temporary hack pending updates to hive code
            },
            -input_ids => [{callgroup => [split(',', $self->o('callgroup'))]}],
            -flow_into => {
                2 => [ 'find_source_bams' ],   # will create a semaphored fan of jobs
            },
      });
    push(@analyses, {
            -logic_name    => 'find_source_bams',
            -module        => 'ReseqTrack::Hive::Process::ImportCollection',
            -meadow_type => 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
            -parameters    => {
                collection_type => $self->o('type_bam'),
                collection_name=> '#callgroup#',
                output_param => 'bam',
            },
            -flow_into => {
                1 => [ 'regions_factory_1' ],
            },
      });
    push(@analyses, {
            -logic_name    => 'regions_factory_1',
            -module        => 'ReseqTrack::Hive::Process::SequenceSliceFactory',
            -meadow_type => 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
            -parameters    => {
                num_bases => $self->o('transpose_window_size'),
                max_sequences => 2000,
                bed => $self->o('target_bed_file'),
            },
            -flow_into => {
                '2->A' => [ 'transpose_bam' ],
                'A->1' => [ 'decide_mergers'],
            },
      });
    push(@analyses, {
            -logic_name    => 'transpose_bam',
            -module        => 'ReseqTrack::Hive::Process::TransposeBam',
            -parameters    => {
                program_file => $self->o('transpose_exe'),
                create_index => 1,
                bed => $self->o('target_bed_file'),
                region_overlap => 100,
            },
            -rc_name => '2Gb',
            -analysis_capacity  =>  4,  # use per-analysis limiter
            -hive_capacity  =>  -1,
            -flow_into => {
                '1->A' => [ 'regions_factory_2' ],
                'A->1' => [ 'collect_vcf' ],
            },
      });
    push(@analyses, {
            -logic_name    => 'regions_factory_2',
            -module        => 'ReseqTrack::Hive::Process::SequenceSliceFactory',
            -parameters    => {
                num_bases => $self->o('call_window_size'),
                max_sequences => 1,
                bed => $self->o('target_bed_file'),
            },
            -rc_name => '200Mb',
            -analysis_capacity  =>  4,  # use per-analysis limiter
            -hive_capacity  =>  -1,
            -flow_into => {
                '2' => [ 'decide_callers',
                          ':////accu?bp_start=[fan_index]',
                          ':////accu?bp_end=[fan_index]',
                          ],
            },
      });
    push(@analyses, {
          -logic_name => 'decide_callers',
          -module        => 'ReseqTrack::Hive::Process::FlowDecider',
          -meadow_type=> 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
          -parameters => {
              require_true => {
                  1 => $self->o('call_by_samtools'),
                  2 => $self->o('call_by_gatk'),
                  3 => $self->o('call_by_freebayes'),
              }
          },
            -flow_into => {
                '1' => [ 'call_by_samtools' ],
                '2' => [ 'call_by_gatk' ],
                '3' => [ 'call_by_freebayes' ],
            },
      });
    push(@analyses, {
            -logic_name    => 'collect_vcf',
            -module        => 'ReseqTrack::Hive::Process::FlowDecider',
            -meadow_type => 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
            -parameters    => {
                require_true => {
                  1 => $self->o('call_by_samtools'),
                  2 => $self->o('call_by_gatk'),
                  3 => $self->o('call_by_freebayes'),
                  4 => 1,
              },
              delete_param => ['bam','bai'],
            },
            -flow_into => {
                1 => [ ':////accu?samtools_vcf=[fan_index]'],
                2 => [ ':////accu?gatk_vcf=[fan_index]'],
                3 => [ ':////accu?freebayes_vcf=[fan_index]'],
                4 => [ ':////accu?bp_start=[fan_index]',
                       ':////accu?bp_end=[fan_index]',
                    ],
            },
      });
    push(@analyses, {
          -logic_name    => 'call_by_samtools',
          -module        => 'ReseqTrack::Hive::Process::RunVariantCall',
          -parameters    => {
              module_name => 'CallBySamtools',
              reference => $self->o('reference'),
              samtools => $self->o('samtools_exe'),
              bcftools => $self->o('bcftools_exe'),
              vcfutils => $self->o('vcfutils_exe'),
              bgzip => $self->o('bgzip_exe'),
              options => $self->o('call_by_samtools_options'),
              region_overlap => 100,
              temp_param_sub => { 1 => [['samtools_vcf','vcf']]}, # temporary hack pending updates to hive code
          },
          -rc_name => '500Mb',
          -analysis_capacity  =>  50,  # use per-analysis limiter
          -hive_capacity  =>  -1,
          -flow_into => {
              1 => [ ':////accu?samtools_vcf=[fan_index]' ],
          },
      });
    push(@analyses, {
          -logic_name    => 'call_by_gatk',
          -module        => 'ReseqTrack::Hive::Process::RunVariantCall',
          -parameters    => {
              module_name => 'CallByGATK',
              reference => $self->o('reference'),
              gatk_dir => $self->o('gatk_dir'),
              options => $self->o('call_by_gatk_options'),
              region_overlap => 100,
              temp_param_sub => { 1 => [['gatk_vcf','vcf']]}, # temporary hack pending updates to hive code
          },
          -rc_name => '2Gb',
          -analysis_capacity  =>  50,  # use per-analysis limiter
          -hive_capacity  =>  -1,
          -flow_into => {
              1 => [ ':////accu?gatk_vcf=[fan_index]' ],
          },
      });
    push(@analyses, {
          -logic_name    => 'call_by_freebayes',
          -module        => 'ReseqTrack::Hive::Process::RunVariantCall',
          -parameters    => {
              module_name => 'CallByFreebayes',
              reference => $self->o('reference'),
              freebayes => $self->o('freebayes_exe'),
              bgzip => $self->o('bgzip_exe'),
              options => $self->o('call_by_freebayes_options'),
              region_overlap => 100,
              temp_param_sub => { 1 => [['freebayes_vcf','vcf']]}, # temporary hack pending updates to hive code
          },
          -rc_name => '2Gb',
          -analysis_capacity  =>  50,  # use per-analysis limiter
          -hive_capacity  =>  -1,
          -flow_into => {
              1 => [ ':////accu?freebayes_vcf=[fan_index]' ],
          },
      });
    push(@analyses, {
          -logic_name => 'decide_mergers',
          -module        => 'ReseqTrack::Hive::Process::FlowDecider',
          -meadow_type=> 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
          -parameters => {
              require_true => {
                  1 => $self->o('call_by_samtools'),
                  2 => $self->o('call_by_gatk'),
                  3 => $self->o('call_by_freebayes'),
              },
              temp_param_sub => {
                1 => [['vcf','samtools_vcf'],['gatk_vcf', 'undef'],['freebayes_vcf','undef'],['caller', '"samtools"']],
                2 => [['vcf','gatk_vcf'],['samtools_vcf', 'undef'],['freebayes_vcf','undef'],['caller', '"gatk"']],
                3 => [['vcf','freebayes_vcf'],['samtools_vcf', 'undef'],['gatk_vcf','undef'],['caller', '"freebayes"']],
              }, # temporary hack pending updates to hive code
          },
            -flow_into => {
                '1' => [ 'merge_vcf' ],
                '2' => [ 'merge_vcf' ],
                '3' => [ 'merge_vcf' ],
            },
      });


    push(@analyses, {
          -logic_name    => 'merge_vcf',
          -module        => 'ReseqTrack::Hive::Process::MergeVcf',
          -parameters    => {
              analysis_label => '#expr("call_by_".$caller)expr#',
              bgzip => $self->o('bgzip_exe'),
              delete_param => ['vcf'],
          },
          -rc_name => '200Mb',
          -analysis_capacity  =>  50,  # use per-analysis limiter
          -hive_capacity  =>  -1,
      });

    return \@analyses;
}

1;

