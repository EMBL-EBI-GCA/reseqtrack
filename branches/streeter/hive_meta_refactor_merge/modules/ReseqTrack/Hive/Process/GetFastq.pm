
package ReseqTrack::Hive::Process::GetFastq;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::ERAUtils qw(get_erapro_conn);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;

    my $run_id = $self->param_required('run_id');
    my $db_params = $self->param_required('reseqtrack_db');
    my $module = $self->param_is_defined('module') ? $self->param('module') : 'ReseqTrack::Tools::GetFastq';
    my $source_root_dir = $self->param_is_defined('source_root_dir') ? $self->param('source_root_dir') : undef;
    my $clobber = $self->param_is_defined('clobber') ? $self->param('clobber') : undef;
    my $era_dbuser = $self->param_required('era_dbuser');
    my $era_dbpass = $self->param_required('era_dbpass');

    eval "require $module" or throw "cannot load module $module $@";

    my $db = ReseqTrack::DBSQL::DBAdaptor->new(%{$db_params});
    my $run = $db->get_RunAdaptor->fetch_by_dbID($run_id);
    throw("did not get run with dbID $run_id") if !$run;

    my $era_db = get_erapro_conn($era_dbuser, $era_dbpass);
    
    my %constructor_hash;
    if ($self->param_is_defined('module_options')) {
      while (my ($key, $value) = each %{$self->param('module_options')}) {
        $constructor_hash{'-'.$key} = $value;
      }
    }
    $constructor_hash{-output_dir} = $self->output_dir;
    $constructor_hash{-run_meta_info} = $rmi;
    $constructor_hash{-source_root_dir} = $source_root_dir;
    $constructor_hash{-clobber} = $clobber;
    $constructor_hash{-db} = $era_db;
    my $fastq_getter = $module->new (%constructor_hash);

    $db->dbc->disconnect_when_inactive(1);

    my $num_files = $fastq_getter->run;

    if (!$num_files) {
      $self->flows_non_factory(undef);
      return;
    }

    $self->output_param('fastq', $fastq_getter->output_files);
    $self->output_param('fastq_md5', $fastq_getter->md5_hash);
}

1;

