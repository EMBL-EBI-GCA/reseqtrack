package ReseqTrack::Hive::Process::RunCramtools;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::Exception qw(throw);
use File::Basename qw(fileparse);
use ReseqTrack::Tools::RunCramtools;
use ReseqTrack::Tools::MetaDataUtils qw( fetch_metadata_object lookup_property );
use ReseqTrack::Tools::FileSystemUtils qw(check_directory_exists check_file_exists);


sub param_defaults {
  return {
    java_exe     => undef,
    program_file => undef,
    jvm_args     => undef,
    command      => undef,
    options      => {},
  };
}

sub run {
  my $self = shift @_;
  $self->param_required('cram');
  my $crams     = $self->param_as_array( 'cram' );
  my $command   = $self->param_required( 'command' );
  my $db_params = $self->param_required( 'reseqtrack_db' );
       
  my %options = %{$self->param( 'options' )};
            
  my @allowed_cmds = ReseqTrack::Tools::RunCramtools->get_valid_commands;
  throw( "Don't recognise command $command. Acceptable commands are: @allowed_cmds")
    if ( !grep { $command eq $_ } @allowed_cmds );
      
  throw('expecting single cram file') unless scalar @{$crams} == 1;
    
  if ( $command eq 'fastq' ){    ## hack for cram to fastq conversion
    my $cram_file_path = $$crams[0];
      
    check_file_exists($cram_file_path);
      
    my ($cram_prefix, $dir, $suffix) = fileparse( $cram_file_path, qr/\..*/ );
    my $prefix = $cram_prefix;
    $options{ prefix => $prefix };           ## seting prefix for fastq file name
      
    ## library layout is required for cram to fastq conversion
    my $library_layout;
          
    unless ( $library_layout ){
      my $db = ReseqTrack::DBSQL::DBAdaptor->new(%{$db_params}); 
      my $file_object = $db->get_FileAdaptor->fetch_by_name($cram_file_path);
      my $file_dbID = $file_object->{dbID};
      
      my $collection_objects = $db->get_CollectionAdaptor->fetch_by_other_id_and_table_name($file_dbID, 'file');
      throw("expecting a single collection") unless scalar @$collection_objects == 1;
        
      $collection_object = $$collection_objects[0];
      my $collection_name = $collection_object->name;
        
      my $meta_objects = fetch_metadata_object( $collection_name, $db);
      throw( "couldn't found metadata object for $collection_name" ) unless $meta_objects;    ## assuming run/experiment name as collection name
        
      $library_layout = lookup_property( $meta_objects, 'library_layout');
      throw('did not get library_layout') unless $library_layout;
        
       $db->dbc->disconnect_when_inactive(1); 
    }
      
  }    
    
  my $cramtools_objects = ReseqTrack::Tools::RunCramtools->new(
              -input_files    => $crams,
              -program        => $self->param('program_file'),
              -working_dir    => $self->output_dir,
              -library_layout => $library_layout,
              -options        =>  \%options,
             );
               
  $self->run_program($cramtools_objects, $command);
               
  my $output_fastqs = $cramtools_object->output_files;           
  $self->output_param('fastq', $output_fastqs);
}

1;