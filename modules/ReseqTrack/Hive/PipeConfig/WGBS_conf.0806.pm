package ReseqTrack::Hive::PipeConfig::WGBS_conf;

use strict;
use warnings;

use base ('ReseqTrack::Hive::PipeConfig::ReseqTrackGeneric_conf');


sub default_options {
    my ($self) = @_;

    return {
        %{ $self->SUPER::default_options() },

        seeding_module  => 'ReseqTrack::Hive::PipeSeed::BasePipeSeed',
        seeding_options => {
            collection_type    => 'FASTQ',
	    output_columns => $self->o('sample_columns'),
            output_attributes => $self->o('sample_attributes'),
            require_columns => $self->o('require_sample_columns'),
            exclude_columns => $self->o('exclude_sample_columns'),
            require_attributes => $self->o('require_sample_attributes'),
            exclude_attributes => $self->o('exclude_sample_attributes'),
        },

        regexs     => undef,

	#Bismark mapper option
	'multicore' => 2, # Sets the number of parallel instances of Bismark to be run concurrently. This forks
	                   # the Bismark alignment step very early on so that each individual Spawn of Bismark
	                   # processes only every n-th sequence (n being set by --multicore).
 
	#Methylation extractor options
	'cutoff' => 5, # The minimum number of times a methylation state has to be seen for that nucleotide
                       #  before its methylation percentage is reported.


        'collection_name'      => '#sample_source_id#',
        'build_collection'     => 1,

	sample_columns => ['sample_id', 'sample_source_id', 'sample_alias'],
        sample_attributes => [],
	experiment_columns => ['instrument_platform', 'paired_nominal_length'],
	experiment_attributes => [],
	study_columns => ['study_source_id'],
	study_attributes => [],
	run_columns => ['run_source_id', 'center_name', 'run_alias'],
	run_attributes => [],

	require_run_attributes => {},
        require_experiment_attributes => {},
        require_study_attributes => {},
        require_sample_attributes => {},
        exclude_run_attributes => {},
        exclude_experiment_attributes => {},
        exclude_study_attributes => {},
        exclude_sample_attributes => {},
	require_experiment_columns => { instrument_platform => ['ILLUMINA'], },
        require_run_columns => { status => ['public', 'private'], },
        require_study_columns => {},
        require_sample_columns => {},
        exclude_sample_columns => {},

        'dir_label_params_list' => [
            "sample_source_id", "experiment_source_id",
            "run_source_id",    "chunk_label"
        ],

        name_file_module       => 'ReseqTrack::Hive::NameFile::BaseNameFile',
        name_file_method       => 'basic',
        unfilt_bam_file_params => {
            new_dir       => '#final_output_dir#/bismark_mapper',
            new_basename  => '#sample_source_id#.bismark.unfilt.GRCh38',
            add_datestamp => 1,
            suffix        => '.bam',
        },
	
	unfilt_bai_file_params => {
            new_dir       => '#final_output_dir#/bismark_mapper',
            new_basename  => '#sample_source_id#.bismark.unfilt.GRCh38',
            add_datestamp => 1,
            suffix        => '.bai',
        },

        unfilt_flagstat_file_params => {
            new_dir       => '#final_output_dir#/bismark_mapper',
            new_basename  => '#sample_source_id#.bismark.unfilt.GRCh38',
            add_datestamp => 1,
            suffix        => '.flagstat',
        },

	dup_flagstat_file_params => {
            new_dir       => '#final_output_dir#/bismark_mapper',
            new_basename  => '#sample_source_id#.dup.bismark.GRCh38',
            add_datestamp => 1,
            suffix        => '.flagstat',
        },

        dedup_name_file_params => {
            new_dir       => '#final_output_dir#/bismark_mapper',
            new_basename  => '#sample_source_id#.dedup.bismark.GRCh38',
            add_datestamp => 1,
            suffix        => '.bam',
        },

        dedup_flagstat_file_params => {
            new_dir       => '#final_output_dir#/bismark_mapper',
            new_basename  => '#sample_source_id#.dedup.bismark.GRCh38',
            add_datestamp => 1,
            suffix        => '.flagstat',
        },

        fastqc_name_file_params => {
            new_dir       => '#final_output_dir#/fastqc',
            add_datestamp => 1,
            suffix        => '.txt',
        },

        fastqc_zip_name_file_params => {
            new_dir       => '#final_output_dir#/fastqc',
            add_datestamp => 1,
            suffix        => '.zip',
        },

        mbiastxt_name_file_params => {
            new_dir       => '#final_output_dir#/bismark_methcall',
            add_datestamp => 1,
            suffix        => '.txt',
        },

        mbiaspng_name_file_params => {
            new_dir       => '#final_output_dir#/bismark_methcall',
            add_datestamp => 1,
            suffix        => '.png',
        },

	bedgraph_name_file_params => {
            new_dir       => '#final_output_dir#/bismark_methcall',
            add_datestamp => 1,
            suffix        => '.gz',
        },

	chhcontext_name_file_params => {
            new_dir       => '#final_output_dir#/bismark_methcall',
            add_datestamp => 1,
            suffix        => '.txt',
        },

	cpgcontext_name_file_params => {
            new_dir       => '#final_output_dir#/bismark_methcall',
            add_datestamp => 1,
            suffix        => '.txt',
        },

	chgcontext_name_file_params => {
            new_dir       => '#final_output_dir#/bismark_methcall',
            add_datestamp => 1,
            suffix        => '.txt',
        },

	splitting_name_file_params => {
            new_dir       => '#final_output_dir#/bismark_methcall',
            add_datestamp => 1,
            suffix        => '.txt',
        },

        'filter_duplicates' => 1,
        'duplicate_flag_value' =>
          '1024',    ## sam flag for PCR or optical duplicate
        'sam_dedup_filter_options' => {
            flag_value  => $self->o('duplicate_flag_value'),
            remove_flag => $self->o('filter_duplicates'),
        },
    };
}

sub resource_classes {
    my ($self) = @_;
    return {
        %{ $self->SUPER::resource_classes }
        ,    # inherit 'default' from the parent class
        '200Mb' => {
                'LSF' => '-C0 -M200 -q '
              . $self->o('lsf_queue')
              . ' -R"select[mem>200] rusage[mem=200]"'
        },
        '500Mb' => {
                'LSF' => '-C0 -M500 -q '
              . $self->o('lsf_queue')
              . ' -R"select[mem>500] rusage[mem=500]"'
        },
        '800Mb' => {
                'LSF' => '-C0 -M800 -q '
              . $self->o('lsf_queue')
              . ' -R"select[mem>800] rusage[mem=800]"'
        },
        '1Gb' => {
                'LSF' => '-C0 -M1000 -q '
              . $self->o('lsf_queue')
              . ' -R"select[mem>1000] rusage[mem=1000]"'
        },
        '2Gb' => {
                'LSF' => '-C0 -M2000 -q '
              . $self->o('lsf_queue')
              . ' -R"select[mem>2000] rusage[mem=2000]"'
        },
        '3Gb' => {
                'LSF' => '-C0 -M3000 -q '
              . $self->o('lsf_queue')
              . ' -R"select[mem>3000] rusage[mem=3000]"'
        },
        '4Gb' => {
                'LSF' => '-C0 -M4000 -q '
              . $self->o('lsf_queue')
              . ' -R"select[mem>4000] rusage[mem=4000]"'
        },
        '5Gb' => {
                'LSF' => '-n 20 -C0 -M5000 -q '
              . $self->o('lsf_queue')
              . ' -R"select[mem>5000] rusage[mem=5000]"'
        },
        '6Gb' => {
                'LSF' => '-C0 -M6000 -q '
              . $self->o('lsf_queue')
              . ' -R"select[mem>6000] rusage[mem=6000]"'
        },
        '8Gb' => {
                'LSF' => '-C0 -M8000 -q '
              . $self->o('lsf_queue')
              . ' -R"select[mem>8000] rusage[mem=8000]"'
        },
        '10Gb' => {
                'LSF' => '-C0 -M10000 -q '
              . $self->o('lsf_queue')
              . ' -R"select[mem>10000] rusage[mem=10000]"'
        },
        '15Gb' => {
                'LSF' => '-C0 -M15000 -q '
              . $self->o('lsf_queue')
              . ' -R"select[mem>15000] rusage[mem=15000]"'
        },
	'70Gb' => {
                'LSF' => '-n 20 -C0 -M70000 -q '
              . $self->o('lsf_queue')
              . ' -R"select[mem>70000] rusage[mem=70000]"'
        },
    };
}

sub pipeline_wide_parameters {
    my ($self) = @_;
    return {
        %{ $self->SUPER::pipeline_wide_parameters },
    };
}


sub hive_meta_table {
    my ($self) = @_;
    return {
        %{ $self->SUPER::hive_meta_table }
        ,    # here we inherit anything from the base class
        'hive_use_param_stack' =>
          1,    # switch on the implicit parameter propagation mechanism
    };
}

sub pipeline_analyses {
    my ($self) = @_;

    my @analyses;

    push(
        @analyses,
        {
            -logic_name  => 'get_seeds',
            -module      => 'ReseqTrack::Hive::Process::SeedFactory',
            -meadow_type => 'LOCAL',
            -parameters  => {
                seeding_module  => $self->o('seeding_module'),
                seeding_options => $self->o('seeding_options'),
            },
            -flow_into => { 2 => ['libraries_factory'], },
        }
    );

    push(@analyses, {
            -logic_name    => 'libraries_factory',
            -module        => 'ReseqTrack::Hive::Process::RunMetaInfoFactory',
            -meadow_type => 'LOCAL',
            -parameters    => {
                factory_type => 'library',
                require_experiment_columns => $self->o('require_experiment_columns'),
                require_study_columns => $self->o('require_study_columns'),
                require_experiment_attributes => $self->o('require_experiment_attributes'),
                require_study_attributes => $self->o('require_study_attributes'),
                exclude_experiment_attributes => $self->o('exclude_experiment_attributes'),
                exclude_study_attributes => $self->o('exclude_study_attributes'),
            },
            -flow_into => {
		'2->A' => [ 'runs_factory' ],
                'A->1' => [ 'store_mapper_report'  ],
            },
	 });

    push(@analyses, {
            -logic_name    => 'runs_factory',
            -module        => 'ReseqTrack::Hive::Process::RunMetaInfoFactory',
            -meadow_type => 'LOCAL',
            -parameters    => {
                factory_type => 'run',
                require_experiment_columns => $self->o('require_experiment_columns'),
                require_study_columns => $self->o('require_study_columns'),
                require_run_columns => $self->o('require_run_columns'),
                require_experiment_attributes => $self->o('require_experiment_attributes'),
                require_study_attributes => $self->o('require_study_attributes'),
                require_run_attributes => $self->o('require_run_attributes'),
                exclude_experiment_attributes => $self->o('exclude_experiment_attributes'),
                exclude_study_attributes => $self->o('exclude_study_attributes'),
                exclude_run_attributes => $self->o('exclude_run_attributes'),

                output_run_columns => $self->o('run_columns'),
                output_study_columns => $self->o('study_columns'),
                output_experiment_columns => $self->o('experiment_columns'),
                output_run_attributes => $self->o('run_attributes'),
                output_study_attributes => $self->o('study_attributes'),
                output_experiment_attributes => $self->o('experiment_attributes'),
            },
            -flow_into => {
                2 => [ 'find_source_fastqs' ],
            },
	 });

    push(@analyses, {
            -logic_name    => 'find_source_fastqs',
            -module        => 'ReseqTrack::Hive::Process::ImportCollection',
            -meadow_type => 'LOCAL',
            -parameters    => {
                collection_type => $self->o('type_fastq'),
                collection_name => '#run_source_id#',
                output_param => 'fastq',
            },
	     -flow_into =>
	     { 1 => { 'fastq_factory' => { 'fastq' => '#fastq#' } }, },
	 });
    
    push(
        @analyses,
        {
            -logic_name  => 'fastq_factory',
            -module      => 'ReseqTrack::Hive::Process::JobFactory',
            -meadow_type => 'LOCAL',
            -parameters  => { factory_value => '#fastq#', },

            -flow_into => {
                '2->A' => { 'run_fastqc' => { 'fastq' => '#factory_value#' } },
                'A->1' => ['store_fastqc_summary'],
            },
        }
	);
 
    push(
        @analyses,
	{
            -logic_name => 'run_fastqc',
            -module     => 'ReseqTrack::Hive::Process::RunFastqc',
            -parameters => {
		output_dir => $self->o('final_output_dir'),
                store_attributes => 1,
                program_file     => $self->o('fastqc_exe'),
            },
            -rc_name           => '500Mb',
            -analysis_capacity => 50,        # use per-analysis limiter
            -hive_capacity     => -1,
            -flow_into         => {
                1 => [
                    ':////accu?fastqc_summary=[]',
                    ':////accu?fastqc_report=[]',
                    ':////accu?fastqc_zip=[]'
                ]
            },
	}
	);

        push(
        @analyses,
	    {
		-logic_name  => 'store_fastqc_summary',
		-module      => 'ReseqTrack::Hive::Process::LoadFile',
		-meadow_type => 'LOCAL',
		-parameters  => {
		    type             => $self->o('fastqc_summary_type'),
		    name_file_module => $self->o('name_file_module'),
		    name_file_method => $self->o('name_file_method'),
		    name_file_params => $self->o('fastqc_name_file_params'),
		    final_output_dir => $self->o('final_output_dir'),
		    collection_name  => $self->o('collection_name'),
		    collect          => $self->o('build_collection'),
		    file             => '#fastqc_summary#',
		},
		-flow_into => {
		    1 =>
		    { 'store_fastqc_report' => { fastqc_summary => '#file#' } },
		},
	    }
	    );

        push(
        @analyses,
	    {
            -logic_name  => 'store_fastqc_report',
            -module      => 'ReseqTrack::Hive::Process::LoadFile',
            -meadow_type => 'LOCAL',
            -parameters  => {
                type               => $self->o('fastqc_report_type'),
                name_file_module   => $self->o('name_file_module'),
                name_file_method   => $self->o('name_file_method'),
                name_file_params   => $self->o('fastqc_name_file_params'),
                final_output_dir   => $self->o('final_output_dir'),
                collection_name    => $self->o('collection_name'),
                collect            => $self->o('build_collection'),
                file               => '#fastqc_report#',
                reseqtrack_options => { delete_param => ['store_fastqc_zip'], }
            },
            -flow_into =>
	    { 1 => { 'store_fastqc_zip' => { fastqc_report => '#file#' } }, },
	    }
	    );

        push(
        @analyses,
	    {
            -logic_name  => 'store_fastqc_zip',
            -module      => 'ReseqTrack::Hive::Process::LoadFile',
            -meadow_type => 'LOCAL',
            -parameters  => {
                type             => $self->o('fastqc_zip_type'),
                name_file_module => $self->o('name_file_module'),
                name_file_method => $self->o('name_file_method'),
                name_file_params => $self->o('fastqc_zip_name_file_params'),
                final_output_dir => $self->o('final_output_dir'),
                collection_name  => $self->o('collection_name'),
                collect          => $self->o('build_collection'),
                file             => '#fastqc_zip#',
            },
            -flow_into => {
                1 => {
                    'bismark_mapper' => {
                        fastq      => '#fastq#',
                        fastqc_zip => '#file#'
                    }
                },
            }
	    }
	    );

        push(
        @analyses,
	    {
            -logic_name => 'bismark_mapper',
            -module     => 'ReseqTrack::Hive::Process::RunBismark',
            -parameters => {
                base         => '#run_source_id#',
                program_file => $self->o('bismark_exe'),
                samtools     => $self->o('samtools_exe'),
                reference    => $self->o('reference'),
                multicore    => $self->o('multicore'),
                command      => 'aln',
                rg_id        => '#run_source_id#',
                rg_sample    => '#sample_source_id#',
                output_dir   => $self->o('final_output_dir').'/bismark_mapper',
            },
            -hive_capacity => 200,
            -rc_name       => '70Gb',
	    -flow_into => {
              1 => [ 
		  ':////accu?bam=[]',':////accu?mapper_report=[]'
		  ]
	    },
	    }
	    );
           
    push(
        @analyses,
	{
	    -logic_name  => 'store_mapper_report',
	    -module      => 'ReseqTrack::Hive::Process::LoadFile',
	    -meadow_type => 'LOCAL',
	    -parameters  => {
		file => '#mapper_report#',
		type             => $self->o('mapper_report_type'),
		name_file_module => $self->o('name_file_module'),
		name_file_method => $self->o('name_file_method'),
		name_file_params =>  {
		    add_datestamp => 1,
		    suffix        => '.txt',
		},
		final_output_dir => $self->o('final_output_dir')."/bismark_mapper",
		collection_name  => $self->o('collection_name'),
		collect          => $self->o('build_collection'),
	    },
	    -rc_name       => '200Mb',
	    -hive_capacity => 200,
	    -flow_into     => {
		1 => {
		    'store_run_bam' => {file => '#bam#' }},
	    },
	}
	);

    push(
        @analyses,
        {
            -logic_name => 'store_run_bam',
            -module     => 'ReseqTrack::Hive::Process::LoadFile',
            -parameters => {
		type             => $self->o('run_bam_type'),
                name_file_module => $self->o('name_file_module'),
                name_file_method => $self->o('name_file_method'),
		name_file_params =>  {
                    add_datestamp => 1,
                    suffix        => '.bam',
		},
                final_output_dir => $self->o('final_output_dir')."/bismark_mapper",
                collection_name  => $self->o('collection_name'),
                collect          => $self->o('build_collection'),
            },
            -analysis_capacity => 20,
            -flow_into =>
            { 1 => { 'bam_sort' => { 'bam' => '#file#' }, }, },
        }
        );

    push(
        @analyses,
        {
            -logic_name => 'bam_sort',
            -module     => 'ReseqTrack::Hive::Process::RunBioBamBam',
            -parameters => {
                biobambam_dir      => $self->o('biobambam_dir'),
                command            => "bamsort",
                output_dir       => $self->o('final_output_dir').'/bismark_mapper',
                options            => { fixmates => 1 },
                create_index       => 0,
            },
            -analysis_capacity => 20,
            -flow_into =>
	    { 1 => { 'decide_merge_libraries' => { 'bam' => '#bam#' }, }, },
        }
	);
	    
    push(@analyses, {
	-logic_name => 'decide_merge_libraries',
	-module        => 'ReseqTrack::Hive::Process::BaseProcess',
	-meadow_type=> 'LOCAL',
	-parameters => {
	    reseqtrack_options => {
		flows_do_count_param => 'bam',
		flows_non_factory => { 
		    1 => 1,
		    2 => 1,
		},
		flows_do_count => {
		    1 => '1',
		    2 => '2+',
		},
	    },
	},
	-flow_into => {
	    1 => [ 'unfilt_flagstat'],
	    2 => { 'merge_libraries' => {'bam' => '#bam#' }},
	},
	 });

    push(@analyses, {
          -logic_name => 'merge_libraries',
          -module        => 'ReseqTrack::Hive::Process::RunBioBamBam',
          -parameters => {
	      reseqtrack_options => {
		  delete_param => ['bam'],
	      },
              biobambam_dir => $self->o('biobambam_dir'),
              command => 'bammerge',
              create_index => 0,
	      output_dir       => $self->o('final_output_dir').'/bismark_mapper',
          },
          -rc_name => '500Mb',
          -hive_capacity  =>  200,
          -flow_into => {
              1 => [ 'unfilt_flagstat'],
          },
	 });


    push(
        @analyses,
        {
            -logic_name => 'unfilt_flagstat',
            -module     => 'ReseqTrack::Hive::Process::RunSamtools',
            -parameters => {
                program_file   => $self->o('samtools_exe'),
                command        => 'flagstat',
                add_attributes => 0,
                reference      => $self->o('reference'),
            },
            -rc_name       => '2Gb',
            -hive_capacity => 200,
            -flow_into     => {
                1 => {
		    'store_unfilt_flagstat' => {
			'file'=> '#metrics#',
			'bam'=> '#bam#'
                    }
                }
	    },
        }
	);


      push(
        @analyses,
	  {
            -logic_name  => 'store_unfilt_flagstat',
            -module      => 'ReseqTrack::Hive::Process::LoadFile',
            -meadow_type => 'LOCAL',
            -parameters  => {
                type             => $self->o('unfilt_flagstat_type'),
                name_file_module => 'ReseqTrack::Hive::NameFile::BaseNameFile',
                name_file_method => $self->o('name_file_method'),
                name_file_params => $self->o('unfilt_flagstat_file_params'),
                final_output_dir => $self->o('final_output_dir'),
                collection_name  => $self->o('collection_name'),
                collect          => $self->o('build_collection'),
            },
            -rc_name       => '200Mb',
            -hive_capacity => 200,
            -flow_into     => {
                1 => {
                    'store_unfilt_bam' => {
                        'file'            => '#bam#',
                        'unfilt_flagstat' => '#file#'
                    }
                }
            },
	  }
	  );

      push(
        @analyses,
	  {
            -logic_name  => 'store_unfilt_bam',
            -module      => 'ReseqTrack::Hive::Process::LoadFile',
            -meadow_type => 'LOCAL',
            -parameters  => {
                type             => $self->o('unfilt_bam_type'),
		final_output_dir => $self->o('final_output_dir'),
                name_file_module => $self->o('name_file_module'),
                name_file_method => $self->o('name_file_method'),
                name_file_params => $self->o('unfilt_bam_file_params'),
                collection_name  => $self->o('collection_name'),
                collect          => $self->o('build_collection'),
            },
            -rc_name       => '200Mb',
            -hive_capacity => 200,
	    -flow_into     => {
                1 => {
                    'mark_duplicates' => {
                        'bam'        => '#file#',
                        'unfilt_bam' => '#file#'
                    },
                },
            }
	  }
	  );

      push(
        @analyses,
	  {
            -logic_name => 'mark_duplicates',
            -module     => 'ReseqTrack::Hive::Process::RunBioBamBam',
            -parameters => {
                biobambam_dir => $self->o('biobambam_dir'),
                command       => 'bammarkduplicates',
                create_index  => 0,
		output_dir => $self->o('final_output_dir')."/bismark_mapper",
            },
            -rc_name       => '800Mb',
            -hive_capacity => 100,
            -flow_into     => { 1 => { 'dup_flagstat' => { 'bam' => '#bam#' } }, },
	  }
	  );

    push(
        @analyses,
	{
            -logic_name => 'dup_flagstat',
            -module     => 'ReseqTrack::Hive::Process::RunSamtools',
            -parameters => {
                program_file   => $self->o('samtools_exe'),
                command        => 'flagstat',
                add_attributes => 1,
                reference      => $self->o('reference'),
            },
            -rc_name       => '2Gb',
            -hive_capacity => 200,
            -flow_into     => {
                1 => {
                    'dup_attributes' => {
			'dup_metrics' => '#metrics#',
                        'dup_attribute_metrics' => '#attribute_metrics#',
                        'bam'  => '#bam#'
                    }
                }
            },
     }
     );

    push(@analyses, {
         -logic_name => 'dup_attributes',
         -module        => 'ReseqTrack::Hive::Process::UpdateAttribute',
         -parameters => {
             attribute_metrics => '#dup_attribute_metrics#',
             collection_type => $self->o('unfilt_bam_type'),
             collection_name => $self->o('collection_name'),
         },
         -rc_name => '200Mb',
         -hive_capacity  =>  200,
         -flow_into => {
             1 => {
                 'store_dup_flagstat' => {
                     'file'=> '#dup_metrics#',
                     'bam'=> '#bam#'
                 }
             },
         },
	 });

    push(
        @analyses,
	{
            -logic_name  => 'store_dup_flagstat',
            -module      => 'ReseqTrack::Hive::Process::LoadFile',
            -meadow_type => 'LOCAL',
            -parameters  => {
                type             => $self->o('dup_flagstat_type'),
                name_file_module => 'ReseqTrack::Hive::NameFile::BaseNameFile',
                name_file_method => 'basic',
                name_file_params => $self->o('dup_flagstat_file_params'),
                final_output_dir => $self->o('final_output_dir'),
                collection_name  => $self->o('collection_name'),
                collect          => $self->o('build_collection'),
            },
            -rc_name       => '200Mb',
            -hive_capacity => 200,
            -flow_into     => {
                1 => {
                    'dedup_bam' => {
                        'bam'            => '#bam#',
                        'dup_flagstat' => '#file#'
                    },
                }
            },
	}
	);

    push(
        @analyses,
        {
            -logic_name => 'dedup_bam',
            -module     => 'ReseqTrack::Hive::Process::RunSamtools',
            -parameters => {
		reseqtrack_options => {
		    delete_param => ['bam'],
		},
                program_file     => $self->o('samtools_exe'),
                command          => 'filter',
                samtools_options => $self->o('sam_dedup_filter_options'),
		output_dir => $self->o('final_output_dir')."/bismark_mapper",
            },
            -rc_name       => '2Gb',
            -hive_capacity => 200,
           -flow_into =>
	    { 1 => { 'dedup_flagstat' => { 'bam' => '#bam#' } }, },
        }
	);

       push(
        @analyses,
	  {
            -logic_name => 'dedup_flagstat',
            -module     => 'ReseqTrack::Hive::Process::RunSamtools',
            -parameters => {
                program_file   => $self->o('samtools_exe'),
                command        => 'flagstat',
                add_attributes => 0,
                reference      => $self->o('reference'),
            },
            -rc_name       => '2Gb',
            -hive_capacity => 200,
            -flow_into     => {
                1 => {
                    'store_dedup_flagstat' => {
                        'file' => '#metrics#',
                        'bam'  => '#bam#'
                    }
                }
            },
	  }
	  );

      push(
        @analyses,
	  {
            -logic_name  => 'store_dedup_flagstat',
            -module      => 'ReseqTrack::Hive::Process::LoadFile',
            -meadow_type => 'LOCAL',
            -parameters  => {
                type             => $self->o('dedup_flagstat_type'),
                name_file_module => 'ReseqTrack::Hive::NameFile::BaseNameFile',
                name_file_method => 'basic',
                name_file_params => $self->o('dedup_flagstat_file_params'),
                final_output_dir => $self->o('final_output_dir'),
                collection_name  => $self->o('collection_name'),
                collect          => $self->o('build_collection'),
            },
            -rc_name       => '200Mb',
            -hive_capacity => 200,
            -flow_into     => {
                1 => {
                    'store_dedup_bam' => {
                        'file'            => '#bam#',
                        'dedup_flagstat' => '#file#'
                    },
                }
            },
	  }
	  );

      push(
        @analyses,
	  {
            -logic_name  => 'store_dedup_bam',
            -module      => 'ReseqTrack::Hive::Process::LoadFile',
            -meadow_type => 'LOCAL',
            -parameters  => {
                type             => $self->o('dedup_bam_type'),
                name_file_module => $self->o('name_file_module'),
                name_file_method => $self->o('name_file_method'),
                name_file_params => $self->o('dedup_name_file_params'),
                final_output_dir => $self->o('final_output_dir'),
                collection_name  => $self->o('collection_name'),
                collect          => $self->o('build_collection'),
            },
            -rc_name       => '200Mb',
            -hive_capacity => 200,
            -flow_into =>
	    { 1 => { 'bam_factory' => { 'dedup_bam' => '#file#' }, }, }
	  }
	  );

    push(
        @analyses,
        {
            -logic_name  => 'bam_factory',
            -module      => 'ReseqTrack::Hive::Process::JobFactory',
            -meadow_type => 'LOCAL',
            -parameters  => { factory_value => '#dedup_bam#', },
            -flow_into   => {
                '2->A' =>
		{ 'bismark_methcall' => { 'bam' => '#factory_value#' } },
                'A->1' => ['mark_seed_complete'],
            },
        }
	);

     push(
        @analyses,
	 {
            -logic_name => 'bismark_methcall',
            -module     => 'ReseqTrack::Hive::Process::RunBismark',
            -parameters => {
                program_file => $self->o('bismark_methcall_exe'),
                command      => 'methext',
                'output_dir' => $self->o('final_output_dir')
                  . "/bismark_methcall",
                cutoff => $self->o('cutoff'),
		multicore => => $self->o('multicore')
            },
            -rc_name       => '5Gb',
            -hive_capacity => 200,
                    -flow_into         => {
                        1 => { 'store_mbias_txt' => { 'file' => '#mbias_txt#' },
                               'store_mbias_png' => { 'file' => '#mbias_png#' },
                               'store_bedgraph' => {'file' => '#bedgraph#' },
                               'store_chhcontext' => {'file' => '#chh_context#' },
                               'store_cpgcontext' => {'file' => '#cpg_context#' },
                               'store_chgcontext' => {'file' => '#chg_context#' },
                               'store_splitting' => {'file' => '#splitting#' },
                        },
                },
	 }
	 );

    push(
        @analyses,
        {
            -logic_name => 'store_mbias_txt',
            -module     => 'ReseqTrack::Hive::Process::LoadFile',
            -meadow_type => 'LOCAL',
            -parameters => {
                clobber => 1,
                type             => $self->o('mbiastxt_type'),
                name_file_module => $self->o('name_file_module'),
                name_file_method => $self->o('name_file_method'),
                name_file_params => $self->o('mbiastxt_name_file_params'),
                final_output_dir => $self->o('final_output_dir'),
                collection_name  => $self->o('collection_name'),
                collect          => $self->o('build_collection')
            },
            -rc_name       => '200Mb',
            -hive_capacity => 200,
        }
        );

    push(
        @analyses,
        {
            -logic_name => 'store_mbias_png',
            -module     => 'ReseqTrack::Hive::Process::LoadFile',
            -meadow_type => 'LOCAL',
            -parameters => {
                type             => $self->o('mbiaspng_type'),
                name_file_module => $self->o('name_file_module'),
                name_file_method => $self->o('name_file_method'),
                name_file_params => $self->o('mbiaspng_name_file_params'),
                final_output_dir => $self->o('final_output_dir'),
                collection_name  => $self->o('collection_name'),
                collect          => $self->o('build_collection'),
            },
            -rc_name       => '200Mb',
            -hive_capacity => 200,
        }
        );

     push(
        @analyses,
	 {
            -logic_name => 'store_bedgraph',
            -module     => 'ReseqTrack::Hive::Process::LoadFile',
            -meadow_type => 'LOCAL',
            -parameters => {
                type             => $self->o('bedgraph_type'),
                name_file_module => $self->o('name_file_module'),
                name_file_method => $self->o('name_file_method'),
                name_file_params => $self->o('bedgraph_name_file_params'),
                final_output_dir => $self->o('final_output_dir'),
                collection_name  => $self->o('collection_name'),
                collect          => $self->o('build_collection'),
            },
            -rc_name       => '200Mb',
            -hive_capacity => 200,
	 }
	 );

    push(
        @analyses,
        {
            -logic_name => 'store_chhcontext',
            -module     => 'ReseqTrack::Hive::Process::LoadFile',
            -meadow_type => 'LOCAL',
            -parameters => {
                type             => $self->o('chhcontext_type'),
                name_file_module => $self->o('name_file_module'),
                name_file_method => $self->o('name_file_method'),
                name_file_params => $self->o('chhcontext_name_file_params'),
                final_output_dir => $self->o('final_output_dir'),
                collection_name  => $self->o('collection_name'),
                collect          => $self->o('build_collection'),
            },
            -rc_name       => '200Mb',
            -hive_capacity => 200,
        }
        );

    push(
        @analyses,
        {
            -logic_name => 'store_cpgcontext',
            -module     => 'ReseqTrack::Hive::Process::LoadFile',
            -meadow_type => 'LOCAL',
            -parameters => {
                type             => $self->o('cpgcontext_type'),
                name_file_module => $self->o('name_file_module'),
                name_file_method => $self->o('name_file_method'),
                name_file_params => $self->o('cpgcontext_name_file_params'),
                final_output_dir => $self->o('final_output_dir'),
                collection_name  => $self->o('collection_name'),
                collect          => $self->o('build_collection'),
            },
            -rc_name       => '200Mb',
            -hive_capacity => 200,
        }
        );

    push(
        @analyses,
        {
            -logic_name => 'store_chgcontext',
            -module     => 'ReseqTrack::Hive::Process::LoadFile',
            -meadow_type => 'LOCAL',
            -parameters => {
                type             => $self->o('chgcontext_type'),
                name_file_module => $self->o('name_file_module'),
                name_file_method => $self->o('name_file_method'),
                name_file_params => $self->o('chgcontext_name_file_params'),
                final_output_dir => $self->o('final_output_dir'),
                collection_name  => $self->o('collection_name'),
                collect          => $self->o('build_collection'),
            },
            -rc_name       => '200Mb',
            -hive_capacity => 200,
        }
        );

    push(
        @analyses,
        {
            -logic_name => 'store_splitting',
            -module     => 'ReseqTrack::Hive::Process::LoadFile',
            -meadow_type => 'LOCAL',
            -parameters => {
                type             => $self->o('splitting_type'),
                name_file_module => $self->o('name_file_module'),
                name_file_method => $self->o('name_file_method'),
                name_file_params => $self->o('splitting_name_file_params'),
                final_output_dir => $self->o('final_output_dir'),
                collection_name  => $self->o('collection_name'),
                collect          => $self->o('build_collection'),
            },
            -rc_name       => '200Mb',
            -hive_capacity => 200,
        }
        );

    push(
        @analyses,
        {
            -logic_name  => 'mark_seed_complete',
            -module      => 'ReseqTrack::Hive::Process::UpdateSeed',
            -parameters  => { is_complete => 1, },
            -meadow_type => 'LOCAL',
        }
    );

    return \@analyses;
}

1;

