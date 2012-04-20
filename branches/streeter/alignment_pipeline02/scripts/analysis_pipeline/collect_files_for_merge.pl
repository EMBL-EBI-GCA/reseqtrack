#!/sw/arch/bin/perl -w

use strict;

use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception;
use Getopt::Long;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4175;
my $dbname;
my $type_key;
my $type_input;
my $type_collection;
my $level;
my $help;

&GetOptions( 
	    'dbhost=s'      => \$dbhost,
	    'dbname=s'      => \$dbname,
	    'dbuser=s'      => \$dbuser,
	    'dbpass=s'      => \$dbpass,
	    'dbport=s'      => \$dbport,
	    'level=s'      => \$level,
	    'type_key=s'       => \$type_key,
	    'type_input=s'       => \$type_input,
	    'type_collection=s'       => \$type_collection,
	    'help!'         => \$help,
	   );


if ($help) {
    exec('perldoc', $0);
    exit(0);
}


my $db = ReseqTrack::DBSQL::DBAdaptor->new(
					   -host => $dbhost,
					   -user => $dbuser,
					   -port => $dbport,
					   -dbname => $dbname,
					   -pass => $dbpass,
					  );

my $ca = $db->get_CollectionAdaptor;
my %inputs;
foreach my $collection (@{$ca->fetch_by_type($type_input)}) {
  $inputs{$collection->name} = $collection->others;
}


if ($level eq 'RUN') {
  my $key_collections = $ca->fetch_by_type($type_key);

  COLLECTION:
  foreach my $key_collection (@$key_collections) {
    my @others;
    foreach my $name (map {$_->name} @{$key_collection->others}) {
      next COLLECTION if (!exists $inputs{$name});
      push(@others, @{$inputs{$name}});
    }
    my $output_collection = ReseqTrack::Collection->new(
            -name => $key_collection->name, -type => $type_collection,
            -others => \@others);
    $ca->store($output_collection);
  }
}

if ($level eq 'LIBRARY') {
  my $rmia = $db->get_RunMetaInfoAdaptor;
  my %library_runs;
  foreach my $rmi (@{$rmia->fetch_all}) {
    push(@{$library_runs{$rmi->sample_id}{$rmi->library_name}}, $rmi->run_id);
  }
  foreach my $sample_id (keys %library_runs) {
    LIBRARY:
    while (my ($library_name, $run_id_list) = each %{$library_runs{$sample_id}}) {
      my @others;
      foreach my $run_id (@$run_id_list) {
        next LIBRARY if (!exists $inputs{$run_id});
        push(@others, @{$inputs{$run_id}});
      }
      my $output_collection = ReseqTrack::Collection->new(
              -name => $sample_id.'_'.$library_name, -type => $type_collection,
              -others => \@others);
      $ca->store($output_collection);
    }
  }
}

if ($level eq 'SAMPLE') {
  my $rmia = $db->get_RunMetaInfoAdaptor;
  my %sample_libraries;
  foreach my $rmi (@{$rmia->fetch_all}) {
    push(@{$sample_libraries{$rmi->sample_id}}, $rmi->library_name);
  }
  SAMPLE:
  foreach my $sample_id (keys %sample_libraries) {
    my @others;
    foreach my $library_name (@{$sample_libraries{$sample_id}}) {
      my $inputs_key = $sample_id.'_'.$library_name;
      next SAMPLE if (!exists $inputs{$inputs_key});
      push(@others, @{$inputs{$inputs_key}});
    }
    my $output_collection = ReseqTrack::Collection->new(
            -name => $sample_id, -type => $type_collection,
            -others => \@others);
    $ca->store($output_collection);
  }
}

=pod

=head1 NAME

reseqtrack/scripts/event/collect_files_for_merge.pl

=head1 SYNOPSIS

This script finds files that are ready to be merged for the next level of processing
(i.e. run level, library level or sample level)
Files must stored in collections, whose names are used to determine how they will be
grouped together in the output collections.

For merging to sample level, the input collections must be named SampleID_LibraryName.
The output collection will be named by the sample_id.

For merging to library level, the input collections must be named by the run_id.
The output collection will be named SampleID_LibraryName.

For merging to run level, a 'key_collection' will be used to determine the names of the
input collections. The name of key_collection must be run_id, and the names of
key_collection->others will be the identifiers. Usually the output collection from
split_fastq.pl should be used as the key_collection.
The output collection will be named by the run_id

=head1 OPTIONS

  database options:

    -dbhost, the name of the mysql-host
    -dbname, the name of the mysql database
    -dbuser, the name of the mysql user
    -dbpass, the database password if appropriate
    -dbport, the port the mysql instance is running on

  other options:

  -level, must be either 'RUN', 'LIBRARY' or 'SAMPLE'

  -type_input, type of the collections that will be searched for files

  -type_collection, type of the output collection that will be written 

  -type_key, type of the 'key_collection' as described above.  Only needed when merging
  to run level

  -help, flag to print this help and exit

=head1 Examples

    $DB_OPTS="-dbhost mysql-host -dbuser rw_user -dbpass **** -dbport 4197 -dbname my_database"

  perl reseqtrack/event/collect_files_for_merge.pl $DB_OPTS -level SAMPLE
    -type_input LIBRARY_BAM -type_collection LIBRARY_BAM_SET

=cut

