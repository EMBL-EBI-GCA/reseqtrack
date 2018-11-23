package ReseqTrack::Tools::ERAUtils;

use strict;
use warnings;
use Exporter;
use File::Basename;
use File::Copy;

use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::History;
use ReseqTrack::DBSQL::ERADBAdaptor;
use ReseqTrack::Tools::FileSystemUtils;
use ReseqTrack::Tools::RunMetaInfoUtils;
use ReseqTrack::Tools::Intersection;
use Data::Dumper;


use vars qw (@ISA  @EXPORT);

@ISA = qw(Exporter);

@EXPORT = qw(get_erapro_conn 
             get_era_fastq 
             compare_era_and_dcc_meta_info
             fix_sample_swap
             get_fastq_details
	   );


sub get_erapro_conn{
  # print 'test_me';
  my ($dbuser, $dbpass, $dbname) = @_;
  
  my $db = ReseqTrack::DBSQL::ERADBAdaptor->new(
    -user   => $dbuser,
    -dbname => $dbname // 'ERAPRO_HX',
    -pass   => $dbpass,
      );
  return $db;
  # print Dumper($db);
  # $VAR1 = bless( {
  #                '_dbc' => bless( {
  #                                   '_port' => 3306,
  #                                   '_host' => 'mysql',
  #                                   '_driver' => 'Oracle',
  #                                   '_password' => ,
  #                                   '_dbname' => 'ERAPRO',
  #                                   '_username' => 'ops$laura',
  #                                   '_timeout' => 0
  #                                 }, 'ReseqTrack::DBSQL::DBConnection' )
  #              }, 'ReseqTrack::DBSQL::ERADBAdaptor' );

}


sub create_run_id_based_ftp_path{
  my ($run_id, $ftp_server) = @_;
  #$run_id =~ /([E|S]RR)\d+/;
  #my $type = $1;
  #$type = lc($type);
  $run_id =~ /([E|S]RR\d\d\d)\d+/;
  my $first_dir = $1;
  my $ftp_path = $ftp_server;
  #$ftp_path .= "/".$type;
  $ftp_path .= "/".$first_dir."/".$run_id."/";
  return $ftp_path;
}

sub compare_era_and_dcc_meta_info{
  my ($era_meta_info, $dcc_meta_info) = @_;
  my $era_hash = get_meta_info_hash($era_meta_info);
  my $dcc_hash = get_meta_info_hash($dcc_meta_info);
  my @era_ids = keys(%$era_hash);
  my @dcc_ids = keys(%$dcc_hash);
  my $era_set = ReseqTrack::Tools::Intersection->new(
    -list => \@era_ids,
      );
  my $dcc_set = ReseqTrack::Tools::Intersection->new(
    -list => \@dcc_ids,
      );
  my $both = $era_set->and($dcc_set);
  my $era_only = $era_set->not($dcc_set);
  my $dcc_only = $dcc_set->not($era_set);
  my %diff_hash;
  my %era_only_hash;
  my %dcc_only_hash;
  foreach my $run_id(@{$both->list}){
    unless(are_run_meta_infos_identical($era_hash->{$run_id}, 
                                        $dcc_hash->{$run_id})){
      #print $run_id." is diff\n";
      $diff_hash{$run_id} = $era_hash->{$run_id} 
    }      
  }
  foreach my $run_id(@{$era_only->list}){
    $era_only_hash{$run_id} = $era_hash->{$run_id};
  }
  foreach my $run_id(@{$dcc_only->list}){
    $dcc_only_hash{$run_id} = $dcc_hash->{$run_id};
  }
  return (\%era_only_hash, \%dcc_only_hash, \%diff_hash);
}

sub fix_sample_swap{
  my ($run_id, $old_sample, $new_sample, $db, $run) = @_;
  if($old_sample eq $new_sample){
    warning($old_sample." is the same as ".$new_sample." not going to run");
    return undef;
  }
  if(!$old_sample || !$new_sample){
    warning("Can't fix ".$run_id." without both an old sample ".$old_sample.
            " and a new sample ".$new_sample);
    return undef;
  }
  my $fa = $db->get_FileAdaptor;
  my $files = $fa->fetch_all_like_name($run_id);
  if(!$files || @$files == 0){
    warning("No files found associated with ".$run_id);
    return undef;
  }
  my @fastq;
  foreach my $file(@$files){
    if($file->filename =~ /\.fastq\.gz$/){
      push(@fastq, $file);
    }else{
      warning($file->name." is associated with ".$run_id." whose sample has ".
              "swapped from ".$old_sample." to ".$new_sample);
    }
  }
  my @non_archived;
  my @archived;
  foreach my $fastq(@fastq){
    if($fastq->name =~ /$old_sample/){
      unless($fastq->name =~ /\/nfs\/1000g-archive\/vol1/){
        push(@non_archived, $fastq);
      }else{
        push(@archived, $fastq);
      }
    }else{
      warning($fastq->name." doesn't have ".$old_sample." in path so doesn't need ".
              "fixing");
    }
  }
  foreach my $file(@non_archived){
    my $new_path = $file->name;
    my $old_path = $file->name;
    $new_path =~ s/$old_sample/$new_sample/;
    #print "mv ".$old_path." ".$new_path."\n";
    if($run){
      if(-e $new_path){
      throw("Can't move ".$old_path." on top of ".$new_path);
       }
      move($old_path, $new_path);
      if(-e $new_path){
        $file->name($new_path);
        my $history = ReseqTrack::History->new(
          -other_id => $file->dbID,
          -table_name => 'file',
          -comment => "Moving from ".$old_path." to ".$new_path,
            );
        $file->history($history);
        $fa->update($file);
      }else{
        throw("Failed to move ".$old_path." to ".$new_path);
      }
    }else{
      print "Need to update ".$file->name." to ".$new_path."\n";
    }
  }
  my ($action, $location);
  if($run && @archived){
    
  }
  my $aa = $db->get_ArchiveAdaptor;
  foreach my $file(@archived){
    my $new_path = $file->name;
    my $old_path = $file->name;
    $new_path =~ s/$old_sample/$new_sample/;
    if($run){
      if(-e $new_path){
        throw("Can't move ".$old_path." on top of ".$new_path);
      }
      my $archive = create_archive_from_object($file); 
      $aa->store($archive);
    }else{
      print "Need to update ".$file->name." to ".$new_path."\n";
    }
  }
  if($run){
    if(@archived >= 1){
      my $find_all = 1;
    UPDATE:while($find_all){
      my $archives = $aa->fetch_all;
      $find_all = 0 if(!$archives || @$archives == 0);
      cleanup_archive($archives, $db);
      if(!$find_all){
        last UPDATE;
      }
      sleep(180);
    }
    }
  }
  return $run_id;
}


sub get_fastq_details{
  my ($run_id, $era_db, $ftp_root) = @_;
  
  my $sql = "select file_name, dir_name, md5, volume_name, bytes from fastq_file ".
      "where run_id = ?";
  my $sth = $era_db->dbc->prepare($sql);
  my %md5hash;
  my %sizehash;
  my %namehash;
  $sth->execute($run_id);
  while(my ($filename, $dirname, $md5, $vol, $bytes) = $sth->fetchrow){
    my $ftp_path = $ftp_root;
    $ftp_path .= "/" unless($ftp_path =~ /\/$/);
    $ftp_path .= $vol."/" if($vol);
    $ftp_path .= $dirname."/";
    $ftp_path .= $filename;
    $md5hash{$ftp_path} = $md5;
    $sizehash{$ftp_path} = $bytes;
    $namehash{$filename} = $ftp_path;
  }
  return (\%md5hash, \%sizehash, \%namehash);
}




sub convert_date{
  my ($old_date) = @_;
  #2009-NOV-16  21:39:50
  print "Given ".$old_date."\n";
  my ($date, $time) = split /\s+/, $old_date;
  my ($year, $month, $day) = split /\-/, $date;
  my $date_hash = date_hash();
  my $num_month = $date_hash->{$month};
  if(!$num_month){
    foreach my $key(keys(%$date_hash)){
      print "Month ".$key." ".$date_hash->{$key}."\n";
    }
    throw("Failed to find a month for ".$month);
  }
  my $new_date = $year."-".$num_month."-".$day;
  $new_date .= " ".$time;
  print "Returning ".$new_date."\n";
  return $new_date;
}



sub date_hash{
  my %hash;
  $hash{'JAN'} = "01";
  $hash{'FEB'} = "02";
  $hash{'MAR'} = "03";
  $hash{'APR'} = "04";
  $hash{'MAY'} = "05";
  $hash{'JUN'} = "06";
  $hash{'JUL'} = "07";
  $hash{'AUG'} = "08";
  $hash{'SEP'} = "09";
  $hash{'OCT'} = 10;
  $hash{'NOV'} = 11;
  $hash{'DEC'} = 12;
  return \%hash;
}

1;
