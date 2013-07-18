
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
            -user => 'g1kro',
            -pass => undef, # set on the command line
            -dbname => $self->o('reseqtrack_db_name'),
        },

        root_output_dir => $self->o('ENV', 'PWD'), # Should be set to something more sensible

        use_label_management => 1, # boolean.  Must be 1 if you want to have structured output directories and meaningful job names
        use_reseqtrack_file_table => 1, # boolean. File paths in job input ids will be converted to e.g. F1234F
        forbid_duplicate_file_id => $self->forbid_duplicate_file_id, # boolean. If true then an index will be built on use_reseqtrack_file_table.  Table will be searched before file name is converted.

    };
}

sub forbid_duplicate_file_id {return 0;}


=head2 pipeline_create_commands

    Description : Implements pipeline_create_commands() interface method of Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that lists the commands that will create and set up the Hive database.
                  In addition to the standard creation of the database and populating it with Hive tables and procedures it also creates two pipeline-specific tables used by Runnables to communicate.

=cut

sub pipeline_create_commands {
    my ($self) = @_;

    my @sql;
    push(@sql, '
    CREATE TABLE reseqtrack_file (
      file_id int(10) unsigned NOT NULL AUTO_INCREMENT,
      name VARCHAR(64000) NOT NULL,
      PRIMARY KEY (file_id)
      )');

    if ($self->forbid_duplicate_file_id) {
      push(@sql, 'CREATE INDEX name_idx ON reseqtrack_file (name(256))');
    }

    return [
        @{$self->SUPER::pipeline_create_commands},  # inheriting database and hive tables' creation

        map {$self->db_execute_command('pipeline_db', $_)} @sql,
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
        'use_label_management' => $self->o('use_label_management'),
        'use_reseqtrack_file_table' => $self->o('use_reseqtrack_file_table'),
        'forbid_duplicate_file_id' => $self->o('forbid_duplicate_file_id'),

    };
}


1;

