#!/sw/arch/bin/perl

use strict;
use Getopt::Long;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception;

$| = 1;


my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4197;
my $dbname;
my @types;
my $ftp_root = '/nfs/1000g-archive/vol1/ftp/';

&GetOptions(
	    'dbhost:s' => \$dbhost,
	    'dbuser:s' => \$dbuser,
	    'dbpass:s' => \$dbpass,
	    'dbport:s' => \$dbport,
	    'dbname:s' => \$dbname,
	    'ftp_root:s' => \$ftp_root,
	   );

my $db = ReseqTrack::DBSQL::DBAdaptor->new
  (
   -host   => $dbhost,
   -user   => $dbuser,
   -pass => $dbpass,
   -port   => $dbport,
   -dbname => $dbname,
  );

my $rmia = $db->get_RunMetaInfoAdaptor;
my $fa = $db->get_FileAdaptor;

my $rmis = $rmia->fetch_all;
my $files = $fa->fetch_by_type('FILTERED_FASTQ');
my $file_hash = run_hash_from_fastq($files);
foreach my $rmi(@$rmis){
  next if($rmi->status eq 'public');
  my $files = $file_hash->{$rmi->run_id};
  foreach my $file(@$files){
    next unless($file->name =~ /$ftp_root/);
    my $new = $file->name;
    $new =~ s/ftp/withdrawn/;
    print $file->name."\t".$new."\n";
  }
}



sub run_hash_from_fastq{
  my ($files) = @_;
  my %hash;
  foreach my $file(@$files){
    $file->filename =~ /([E|S]RR\d+)/;
    push(@{$hash{$1}}, $file);
  }
  return \%hash;
}
