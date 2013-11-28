#!/usr/bin/env perl
use strict;
use warnings;
use ReseqTrack::Pipeline;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception;
use Getopt::Long;

my ($dbhost, $dbuser, $dbpass, $dbport, $dbname);
my $file;
my ($read, $write, $update, $help);

&GetOptions(
    'dbhost=s'  => \$dbhost,
    'dbname=s'  => \$dbname,
    'dbuser=s'  => \$dbuser,
    'dbpass=s'  => \$dbpass,
    'dbport=s'  => \$dbport,
    'file=s'    => \$file,
    'read!'     => \$read,
    'write!'    => \$write,
    'help!'     => \$help,
);

if ($help) {
    perldocs();
}

throw("Must specify file (-file) to read from or write to") if !$file;
throw("Must specify -read or -write") if (!$read && !$write) || ($read && $write);

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
    -host   => $dbhost,
    -user   => $dbuser,
    -port   => $dbport,
    -dbname => $dbname,
    -pass   => $dbpass,
);
my $pa = $db->get_PipelineAdaptor;

if ($read) {
  foreach my $pipeline (@{parse_file($file)}) {
    if (my $existing = $pa->fetch_by_name($pipeline->name)) {
      $pipeline->dbID($existing->dbID);
      create_history($pipeline, $existing);
      $pa->update($pipeline);
    }
    else {
      $pa->store($pipeline);
    }
  }
}
elsif ($write) {
  my $pipelines = $pa->fetch_all;
  dump_pipelines($file, $pipelines);
}

#######################################################
sub dump_pipelines {
  my ($file, $pipelines) = @_;
  open my $fh, '>', $file or throw("could not open $file $!");
  foreach my $pipeline (@$pipelines) {
    print $fh "[" , $pipeline->name , "]\n";
    print $fh "table_name=", $pipeline->table_name, "\n";
    print $fh "config_module=", $pipeline->config_module, "\n";
    print $fh "config_options=", $pipeline->config_options, "\n" if $pipeline->config_options;
    print $fh "\n";
  }
  close $fh;
}

sub parse_file {
  my ($file) = @_;
  my @pipelines;
  open my $fh, '<', $file or throw("could not open $file $!");
  PIPELINE:
  while (1) {
    my $line = <$fh>;
    last PIPELINE if !$line;
    redo PIPELINE if $line !~ /\S/; # blank line
    my ($name) = $line =~ /\[(.*)\]/;
    throw("could not find a pipeline name in $line") if !$name;
    my $pipeline = ReseqTrack::Pipeline->new( -name => $name);
    my @config_options;
    LINE:
    while (my $line = <$fh>) {
      last LINE if $line !~ /\S/; # blank line
      chomp $line;
      my ($key, $val) = split('=', $line, 2);
      throw("could not interpret $line") if !defined $val;
      if ($key eq 'table_name') {
        $pipeline->table_name($val);
      }
      elsif ($key eq 'config_module') {
        $pipeline->config_module($val);
      }
      elsif ($key eq 'config_options') {
        push(@config_options, $val);
      }
      else {
        throw("could not interpret key $key");
      }
    }
    $pipeline->config_options(join(' ', @config_options));
    push(@pipelines, $pipeline);
  }
  close $fh;
  return \@pipelines;
}

sub create_history {
  my ($new, $old) = @_;
  my @comments;
  if ($new->table_name ne $old->table_name) {
    push(@comments, "table_name changed from ".$old->table_name." to ".$new->table_name);
  }
  if ($new->config_module ne $old->config_module) {
    push(@comments, "config_module changed from ".$old->config_module." to ".$new->config_module);
  }
  if ($new->config_options ne $old->config_options) {
    push(@comments, "config_options changed from ".$old->config_options." to ".$new->config_options);
  }
  foreach my $comment (@comments) {
    my $history = ReseqTrack::History->new(
      -other_id => $old->dbID,
      -table_name => 'pipeline',
      -comment => $comment,
      );
    $new->history($history);
  }
}

sub perldocs {
    exec('perldoc', $0);
    exit(0);
}

=pod

=head1 NAME

ReseqTrack/scripts/event/load_pipeline_from_conf.pl

=head1 SYNOPSIS

    This script takes a config file in the following format and parse it to create
    and load pipeline objects into the event table

        [pipeline_name]
        table_name=file
        config_module=ReseqTrack::Hive::PipeConfig::MyModule
        config_options= (any flags you want to give to the hive init_pipeline.pl script)

    It can also dump existing events in this format


=head2 OPTIONS

        -dbhost, the name of the mysql-host
        -dbname, the name of the mysql database
        -dbuser, the name of the mysql user
        -dbpass, the database password if appropriate
        -dbport, the port the mysql instance is running on, this defaults to 4197 the 
         standard port for mysql-g1kdcc.ebi.ac.uk
        -file, name of file to read in or write out
        -read, binary flag to indicate the file specified by -file should be parsed and 
         loaded into the database
        -write, binary flag to indicate the file specified should be written in the given
          format based on the workflow objects already in the database
        -help, binary flag to indicate the help should be printed


=head1 Example:

$DB_OPTS="-dbhost mysql-host -dbuser rw_user -dbpass **** -dbport 4197 -dbname my_database"


    This dumps the existing event objects into a file of the given name:
    perl ReseqTrack/scripts/pipeline/load_pipeline_from_conf.pl  $DB_OPTS -file /path/to/pipeline.conf -write
   
    
    This will load the pipeline objects specifed in pipeline.conf:
    perl ReseqTrack/scripts/pipeline/load_pipeline_from_conf.pl  $DB_OPTS -file /path/to/pipeline.conf -read
    The write process will overwrite any existing file with the given name

=cut

