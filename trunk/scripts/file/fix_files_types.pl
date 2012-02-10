#!/sw/arch/bin/perl -w

use strict;
use Getopt::Long;
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileUtils qw(create_history copy_file_object);
use ReseqTrack::Tools::FileSystemUtils;

$| = 1;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4197;
my $dbname;
my $file_list;
my $new_type;
my $old_type;

&GetOptions(
	    'dbhost=s'       => \$dbhost,
	    'dbname=s'       => \$dbname,
	    'dbuser=s'       => \$dbuser,
	    'dbpass=s'       => \$dbpass,
	    'dbport=s'       => \$dbport,
	    'file_list=s' => \$file_list,
	    'new_type=s' => \$new_type,
	    'old_type=s' => \$old_type,
	    );


my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host => $dbhost,
  -user => $dbuser,
  -port => $dbport,
  -dbname => $dbname,
  -pass => $dbpass,
    );


my $list = get_lines_from_file($file_list);

#my %hash;
#foreach my $file(@$list){
#  $hash{$file} = 1;
#}
my $fa = $db->get_FileAdaptor;

foreach my $path(@$list){
  my $file = $fa->fetch_by_name($path);
  unless($file){
    print "Failed to get file object from ".$path."\n";
    next;
  }
  if($old_type){
    unless($file->type eq $old_type){
      print $file->type." does not match ".$old_type." skipping\n";
      next;
    }
  }
  my $new_file = copy_file_object($file);
  $new_file->type($new_type);
  my $history = create_history($new_file, $file);
  $new_file->history($history);
  $fa->update($new_file);
}
