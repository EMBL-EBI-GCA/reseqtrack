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

&GetOptions(
	    'dbhost:s' => \$dbhost,
	    'dbuser:s' => \$dbuser,
	    'dbpass:s' => \$dbpass,
	    'dbport:s' => \$dbport,
	    'dbname:s' => \$dbname,
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
my $filtered_fastq = $fa->fetch_by_type('FILTERED_FASTQ');
my $withdrawn_fastq = $fa->fetch_by_type('WITHDRAWN_FILTERED_FASTQ');
my $archive_fastq = $fa->fetch_by_type('ARCHIVE_FASTQ');


my $filtered_hash = run_hash_from_fastq($filtered_fastq);
my $withdrawn_hash = run_hash_from_fastq($withdrawn_fastq);
my $archive_hash = run_hash_from_fastq($archive_fastq);


foreach my $rmi(@$rmis){
  my @files;
  my $sample_name = $rmi->sample_name;
  push(@files, @{$filtered_hash->{$rmi->run_id}}) if($filtered_hash->{$rmi->run_id});
  push(@files, @{$withdrawn_hash->{$rmi->run_id}}) if($withdrawn_hash->{$rmi->run_id});
  push(@files, @{$archive_hash->{$rmi->run_id}}) if($archive_hash->{$rmi->run_id});
  my @problems;
  foreach my $file(@files){
    next if($file->name =~ /$sample_name/);
    push(@problems, $file);
  }
  foreach my $problem(@problems){
    my $name = $problem->name;
    my $new_name = $problem->name;
    $name =~ /([NA|HG]\d+)/;
    my $old_sample = $1;

    $new_name =~ s/$old_sample/$sample_name/;
    if($problem->type =~ /^WITHDRAWN/){
      $new_name =~ s/withdrawn/ftp/;
    }
    print $name."\t".$new_name."\n";
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
