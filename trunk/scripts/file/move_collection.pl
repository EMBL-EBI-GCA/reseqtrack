#!/usr/bin/env perl

use strict;
use warnings;

use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::Tools::GeneralUtils;
use ReseqTrack::Tools::Exception;
use File::Basename;
use Getopt::Long;
use ReseqTrack::Tools::RunMetaInfoUtils qw(create_directory_path);
use Data::Dumper;
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
);

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

if ($run_meta_info) {
    $output_dir =
      create_directory_path( $run_meta_info, $directory_layout, $output_dir );
}

for my $f ( @{ $c->others } ) {

    my $old_path = $f->name;
    my $new_path = $output_dir . '/' . basename( $f->name );

    if ( $new_path eq $old_path ) {
        print
          "Destination is the same as the source, aborting move $old_path $/"
          if ( !$run );
    }
    elsif ($run) {
        move_file_in_db_and_dir( [$f], $output_dir, $f->type, $db );
    }
    else {
        print "Would move file from $old_path to $new_path$/" if ( !$run );
    }
}

