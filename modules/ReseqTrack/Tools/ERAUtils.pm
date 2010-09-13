package ReseqTrack::Tools::ERAUtils;

use strict;
use warnings;
use Exporter;
use Cwd;
use File::Basename;
use File::Copy;

use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::History;
use ReseqTrack::DBSQL::ERADBAdaptor;
use ReseqTrack::Tools::FileSystemUtils;
use ReseqTrack::Tools::RunMetaInfoUtils;
use ReseqTrack::Tools::Intersection;

use vars qw (@ISA  @EXPORT);

@ISA = qw(Exporter);

@EXPORT = qw(get_erapro_conn 
             get_era_fastq 
             compare_era_and_dcc_meta_info
             fix_sample_swap
             get_fastq_details
	     convert_population);


sub get_erapro_conn{
  my ($dbuser, $dbpass) = @_;
  my $dbname = 'ERAPRO';
  my $dbhost = 'ERAPRO';
  my $dbport = 1561;
  my $db = ReseqTrack::DBSQL::ERADBAdaptor->new(
    -host   => $dbhost,
    -user   => $dbuser,
    -port   => $dbport,
    -dbname => $dbname,
    -pass   => $dbpass,
      );
  return $db;
}


sub get_era_fastq{
  my ($run_id, $output_dir, $ftp_server, $clobber, $wget_program, 
      $wget_options, $era_db) = @_;
  #check for existing files
  throw("Can't get fastq files if no run_id is specified") unless($run_id);
  $wget_program = "wget" unless($wget_program);
  $wget_options = '-t 1 --ignore_length' unless($wget_options);
  $output_dir = getcwd() unless($output_dir);
  #check to see if files already exist;
  my $existing_files = find_file($run_id."*", $output_dir);
  if($existing_files && @$existing_files > 0){
    if($clobber){
      print STDERR "Deleting files associated with ".$run_id." in ".$output_dir."\n";
      foreach my $path(@$existing_files){
        print STDERR "Deleting ".$path."\n";
        unlink $path;
      }
    }else{
      print STDERR @$existing_files." files already exist starting with ".
          $run_id."\n";
      foreach my $path(@$existing_files){
        print STDERR $path."\n";
      }
      print STDERR "Can not continue with existing files\n";
      return undef;
    }
  }
  #get fastq paths/md5s/sizes
  my ($md5hash, $sizehash, $namehash) = get_fastq_details($run_id, $era_db, $ftp_server);
  my @paths = keys(%$md5hash);
  if(!@paths || @paths == 0){
    throw("Failed to find any paths for ".$run_id." in ".$era_db->dbc->dbname);
  }
  my $wget_stem = create_wget_stem($output_dir, $wget_program, $wget_options);
  foreach my $path(@paths){
    my $cmd = $wget_stem." ".$path;
    my $exit = system($cmd);
    if($exit >= 1){
      throw("Failed to run ".$cmd." exited with ".$exit);
    }
  }
  my ($list, $hash) = list_files_in_dir($output_dir, 1);
  my %md5_hash;

  foreach my $file(@$list){
    next unless($file =~ /$run_id/);
    my $name = basename($file);
    my $ftppath = $namehash->{$name};
    foreach my $name(keys(%$namehash)){
      print $name." ".$namehash->{$name}."\n";
    }
    if(!$ftppath){
      throw("Failed to find an ftp path for ".$name);
    }
    my $md5 = $md5hash->{$ftppath};
    print "Getting size using ftp path ".$ftppath."\n";
    my $size = $sizehash->{$ftppath};
    my $output_size = -s $file;
    unless($size == $output_size){
      throw("There is a problem with the output from ".$ftppath." the database ".
            "gives ".$size." but perl gives ".$output_size);
    }
    my $output_md5 = run_md5($file);
    unless($output_md5 eq $md5){
       throw("There is a problem with the output from ".$ftppath." the database ".
            "gives ".$md5." but run_md5 gives ".$output_md5);
    }
    $md5_hash{$file} = $md5;
  }
  #return md5 hash
  return \%md5_hash;
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
sub create_wget_stem{
  my ($output_dir, $wget_program, $wget_options) = @_;
  my $cmd_stem = $wget_program." ".$wget_options;
  $cmd_stem .= " -P ".$output_dir if($output_dir);
  return $cmd_stem;
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
  #print $sql."\n";
  #prunt $run_id."\n";
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

sub convert_population{
  my ($string, $run_id, $study_id) = @_;
  my $pop;
  if($string =~ /yri/i){
    $pop = 'YRI';
  }elsif($string =~ /yoruba/i){
    $pop = 'YRI';
  }elsif($string =~ /southern\s+han\s+chinese/i){
    $pop = 'CHS';
  }elsif($string =~ /CHS/i){
    $pop = 'CHS';
  }elsif($string =~ /han chinese/i){
    $pop = 'CHB';
  }elsif($string =~ /CHB/i){
    $pop = 'CHB';
  }elsif($string =~ /japan/i){
    $pop = 'JPT';
  }elsif($string =~ /JPT/i){
    $pop = 'JPT';
  }elsif($string =~ /CEU/i){
    $pop = 'CEU'; 
  }elsif($string =~ /CEPH/i){
    $pop = 'CEU';
  }elsif($string =~ /tuscan/i){
    $pop = 'TSI';
  }elsif($string =~ /TSI/i){
    $pop = 'TSI';
  }elsif($string =~ /denver/i){
    $pop = 'CHD';
  }elsif($string =~ /CHD/i){
    $pop = 'CHD';
  }elsif($string =~ /Luhya/i){
    $pop = 'LWK';
  }elsif($string =~ /LWK/i){
    $pop = 'LWK';
  }elsif($string =~ /UTAH/i){
    $pop = 'CEU';
  }elsif($string =~ /ASW/){
    $pop = 'ASW';
  }elsif($string =~ /MXL/){
    $pop = 'MXL';
  }elsif($string =~ /African-American/){
    $pop = 'ASW';
  }elsif($string =~ /Mexican-American/){
    $pop = 'MXL';
  }elsif($string =~ /UK/){
    $pop = 'GBR';	
  }elsif($string =~ /British/){
    $pop = 'GBR';
  }elsif($string =~ /British\s+\(GBR\)/){
    $pop = 'GBR';
  }elsif($string =~ /GBR/i){
    $pop = 'GBR';
  }elsif($string =~ /FIN/i){
    $pop = 'FIN';
  }elsif($string =~ /SHC/){
    $pop = 'CHS';
  }elsif($string =~ /Puerto\s+Rican/i){
    $pop = 'PUR';
  }elsif($string =~ /pur/i){
    $pop = 'PUR';
  }elsif($string =~ /Colombian/){
    $pop = 'CLM';
  }elsif($string =~ /CLM/){
    $pop = 'CLM';
  }elsif($string =~ /Gujarati/){
    $pop = 'GIH';
  }elsif($string =~ /Maasai/){
    $pop = 'MKK';
  }elsif($string =~ /Spanish/i){
    $pop = 'IBS';
  }elsif($string =~ /IBS/i){
    $pop = 'IBS';
  }else{
    throw("Failed to find pop for ".$string." ".$run_id." ".$study_id);
   }
  return $pop;
}


1;
