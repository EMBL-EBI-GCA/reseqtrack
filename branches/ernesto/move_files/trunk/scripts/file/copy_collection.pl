#!/sw/arch/bin/perl -w

use strict;

use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::Tools::GeneralUtils;
use ReseqTrack::Tools::Exception;
use File::Basename;
use Getopt::Long;
use ReseqTrack::Tools::RunMetaInfoUtils qw(create_directory_path);
use ReseqTrack::Tools::FileSystemUtils qw(run_md5);
use File::Path qw(make_path);
use Data::Dumper;
use File::Copy;

$| = 1;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport;
my $dbname;
my $run = 0;
my $name;
my $type;
my $directory_layout;
my $output_dir;
my $run_id_regex    = '[ESD]RR\d{6}';
my $sample_id_regex = '[ESD]RS\d{6}';
my $experiment_id_regex = '[ESD]RX\d{6}';
my $check_md5;
my $file_mode;

&GetOptions(
  'dbhost=s'           => \$dbhost,
  'dbname=s'           => \$dbname,
  'dbuser=s'           => \$dbuser,
  'dbpass=s'           => \$dbpass,
  'dbport=s'           => \$dbport,
  'name=s'             => \$name,
  'type=s'             => \$type,
  'output_dir=s'       => \$output_dir,
  'directory_layout=s' => \$directory_layout,
  'run!'               => \$run,
  'run_id_regex=s'     => \$run_id_regex,
  'sample_id_regex=s'  => \$sample_id_regex,
  'check_md5!'         => \$check_md5,
  'file_mode=s'             => \$file_mode,
);

$file_mode = oct($file_mode) if ($file_mode);

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $dbhost,
  -user   => $dbuser,
  -port   => $dbport,
  -dbname => $dbname,
  -pass   => $dbpass,
);

my $ca = $db->get_CollectionAdaptor;
my $fa = $db->get_FileAdaptor;
my $ra = $db->get_RunMetaInfoAdaptor;

my $c = $ca->fetch_by_name_and_type( $name, $type );
throw("Failed to find a collection for $name $type") unless ($c);
throw("Collection $name $type does not contain files, nothing to move")
  unless ( $c->table_name eq 'file' );
throw("Directory layout is required") unless ($directory_layout);
throw("Output dir is required")       unless ($output_dir);

my $run_meta_info;
if ( $name =~ /$run_id_regex/ ) {
  $run_meta_info = $ra->fetch_by_run_id($&);
}
elsif ( $name =~ /$sample_id_regex/ ) {
  my $rmi_list = $ra->fetch_by_sample_id($&);
  $run_meta_info = $rmi_list->[0] if (@$rmi_list);
}
elsif ( $name =~ /$experiment_id_regex/ ) {
  my $rmi_list = $ra->fetch_by_experiment_id($&);
  $run_meta_info = $rmi_list->[0] if (@$rmi_list);
}

if ($run_meta_info) {
  $output_dir =
    create_directory_path( $run_meta_info, $directory_layout, $output_dir );
}

my ($err, %opts);
$opts{err}     = \$err;
$opts{mask}    = $file_mode if $file_mode;
$opts{verbose} = 1;

make_path( $output_dir, \%opts );

for my $f ( @{ $c->others } ) {

  my $old_path = $f->name;
  my $new_path = $output_dir . '/' . basename( $f->name );

  my $old_size = $f->size || -s $old_path;

  if ( $new_path eq $old_path ) {
    print "Destination is the same as the source, aborting copy $old_path $/"
      if ( !$run );
  }
  elsif ($run) {
    my $cp_success = copy( $old_path, $new_path );

    throw("Copy from $old_path to $new_path failed: $!") unless ($cp_success);

    my $new_size = -s $new_path;

    throw("Files size mismatch between $old_path and $new_path, copy failed")
      unless ( $new_size == $old_size );

    if ($check_md5) {
      my $old_md5 = $f->md5 || run_md5($old_path);
      my $new_md5 = run_md5($new_path);
      throw("md5 mismatch between $old_path and $new_path, copy failed")
        unless ( $old_md5 eq $new_md5 );
    }
    
    chmod($file_mode, $new_path) if $file_mode;
  }
  else {
    print "Would copy file from $old_path to $new_path$/";
  }
}

