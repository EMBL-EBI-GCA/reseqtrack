
package ReseqTrack::Hive::Process::GetFastq;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::ERAUtils qw(get_erapro_conn);


sub param_defaults {
  return {
    module => 'ReseqTrack::Tools::GetFastq',
    source_root_dir => undef,
    clobber => undef,
    module_options => {},
    era_dbname => undef,
  };
}


sub run {
    my $self = shift @_;

    my $run_id = $self->param_required('run_id');
    my $db_params = $self->param_required('reseqtrack_db');
    my $module = $self->param('module') // param_defaults()->{'module'};
    my $source_root_dir = $self->param('source_root_dir');
    my $clobber = $self->param('clobber');
    my $era_dbuser = $self->param_required('era_dbuser');
    my $era_dbpass = $self->param_required('era_dbpass');
    my $era_dbname = $self->param('era_dbname');

    eval "require $module" or throw "cannot load module $module $@";

    my $db = ReseqTrack::DBSQL::DBAdaptor->new(%{$db_params});
    my $run = $db->get_RunAdaptor->fetch_by_dbID($run_id);
    throw("did not get run with dbID $run_id") if !$run;

    my $era_db = get_erapro_conn($era_dbuser, $era_dbpass, $era_dbname);
    
    my %constructor_hash;
    while (my ($key, $value) = each %{$self->param('module_options')}) {
      $constructor_hash{'-'.$key} = $value;
    }
    $constructor_hash{-output_dir} = $self->output_dir;
    $constructor_hash{-run_info} = $run;
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

