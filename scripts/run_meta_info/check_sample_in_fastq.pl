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
  my $sample_name = $rmi->sample_name;
  if($sample_name =~ /K$/){
    throw($sample_name." seems to be modified");
  }
  my $run_id = $rmi->run_id;
  my $filtered_fastq = $filtered_hash->{$run_id};
  my $withdrawn_fastq = $withdrawn_hash->{$run_id};
  my $archive_fastq = $archive_hash->{$run_id};
  unless(($filtered_fastq && @$filtered_fastq >= 1) ||
	 ($withdrawn_fastq && @$withdrawn_fastq >= 1)){
    #print STDERR $run_id." has no filtered or withdrawn fastq\n";
    next;
  }
  my @bad_files;
  if($filtered_fastq && @$filtered_fastq >= 1){
    my $tmp = check_sample($filtered_fastq, $sample_name);
    push(@bad_files, @$tmp) if($tmp && @$tmp >= 1);
  }
  if($withdrawn_fastq && @$withdrawn_fastq >= 1){
    my $tmp = check_sample($withdrawn_fastq, $sample_name);
    push(@bad_files, @$tmp) if($tmp && @$tmp >= 1);
  }
  if($archive_fastq && @$archive_fastq >= 1){
    my $tmp = check_sample($archive_fastq, $sample_name);
    push(@bad_files, @$tmp) if($tmp && @$tmp >= 1);
  }
  print STDERR $run_id." ".$sample_name." has ".@bad_files.
    " files which need to be fixed\n" if(@bad_files >= 1);
  foreach my $file(@bad_files){
    my $new = $file->name;
    my $old_sample;
    $new =~ /(HG\d+)/;
    $old_sample = $1;
    unless($old_sample){
      $new =~ /(NA\d+)/;
      $old_sample = $1;
    }
    unless($old_sample  =~ /[NA|HG]\d+/){
      throw("Failed to get old sample from ".$new." have ".$old_sample." instead");
    }
   if($sample_name eq $old_sample){
      print STDERR "Something is wrong for ".$run_id." ".$sample_name." matches ".
	" ".$old_sample."\n";
    }
    if($sample_name =~ /K$/){
      throw($sample_name." seems to be modified");
    }
    if($file->type =~ /^WITHDRAWN/){
      $new =~ s/withdrawn/ftp/;
    }
    print $file->name." ".$new."\n";
  }
}


sub check_sample{
  my ($files, $name) = @_;
  my @bad_files;
  foreach my $file(@$files){
    unless($file->name =~ /$name/){
      push(@bad_files, $file);
    }
  }
  return \@bad_files;
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
