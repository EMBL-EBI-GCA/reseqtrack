
package ReseqTrack::Hive::PipeConfig::ReseqTrackGeneric_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');  # All Hive databases configuration files should inherit from HiveGeneric, directly or indirectly


=head2 default_options

    Description : Implements default_options() interface method of Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that is used to initialize default options.
                  In addition to the standard things it defines two options, 'first_mult' and 'second_mult' that are supposed to contain the long numbers to be multiplied.

=cut

sub default_options {
    my ($self) = @_;

    my $super_default_options = $self->SUPER::default_options();
    # override defaults set in parent class:
    $super_default_options->{'pipeline_db'}->{'-port'} = 4175;
    $super_default_options->{'pipeline_db'}->{'-user'} = 'g1krw';

    return {
        %{ $super_default_options },

        'host' => 'mysql-g1k',
        'reseqtrack_db'  => {
            -host => $self->o('host'),
            -port => 4175,
            -user => undef, # set on the command line
            -pass => undef, # set on the command line
            -dbname => $self->o('reseqtrack_db_name'),
        },

        root_output_dir => $self->o('ENV', 'PWD'), # Should be set to something more sensible
        dir_label_params => [],
        lsf_queue => 'production',
	lsf_resource => '',
	java_exe => '/usr/bin/java',
    }; ## when lsf_resource is null, the LSF will try to figure out the filesystem to use based on file path
}


=head2 pipeline_create_commands

    Description : Implements pipeline_create_commands() interface method of Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that lists the commands that will create and set up the Hive database.
                  In addition to the standard creation of the database and populating it with Hive tables and procedures it also creates two pipeline-specific tables used by Runnables to communicate.

=cut

sub pipeline_create_commands {
    my ($self) = @_;

   
    my $reseqtrack_file_sql =  '
    CREATE TABLE reseqtrack_file (
      file_id int(10) unsigned NOT NULL AUTO_INCREMENT,
      name VARCHAR(64000) NOT NULL,
      PRIMARY KEY (file_id)
      )';

    return [
        @{$self->SUPER::pipeline_create_commands},  # inheriting database and hive tables' creation
        'db_cmd.pl -url '.$self->pipeline_url()." -sql '$reseqtrack_file_sql'",
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

        'reseqtrack_db' => $self->o('reseqtrack_db'),
        'root_output_dir' => $self->o('root_output_dir'),
        'dir_label_params' => $self->o('dir_label_params'),

    };
}

sub hive_meta_table {
    my ($self) = @_;
    return {
        %{$self->SUPER::hive_meta_table},          # here we inherit anything from the base class
        'hive_use_param_stack' => 1,
    };
}


1;

