#!/sw/bin/perl -w

use strict;
use ReseqTrack::InputString;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception;
use ReseqTrack::Tools::FileSystemUtils;
use Getopt::Long;


my $dbhost;
my $dbuser;
my $dbpass;
my $dbport;
my $dbname;

my $file;
my $type;
my $run;

&GetOptions(
  'dbhost=s'  => \$dbhost,
  'dbname=s'  => \$dbname,
  'dbuser=s'  => \$dbuser,
  'dbpass=s'  => \$dbpass,
  'dbport=s'  => \$dbport,
  'type=s' => \$type,
  'file=s' => \$file,
  'run=s' => \$run,
    );

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
    -host   => $dbhost,
    -user   => $dbuser,
    -port   => $dbport,
    -dbname => $dbname,
    -pass   => $dbpass,
);

my $lines = get_lines_from_file($file);
my $isa = $db->get_InputStringAdaptor;
my @inputs;
foreach my $line(@$lines){
  my @values = split /\t/, $line;
  $type = $values[1] unless($type);
  my $input_string = create_input_string($values[0], $type);
  push(@inputs, $input_string);
}

foreach my $input(@inputs){
  $isa->store($input);
}

sub create_input_string{
  my ($name, $type) = @_;
  my $input_string = ReseqTrack::InputString->new(
    -name => $name,
    -type => $type,
      );
  return $input_string;
}
