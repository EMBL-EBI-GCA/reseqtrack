package ReseqTrack::DBSQL::ERARunMetaInfoAdaptor;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::DBSQL::BaseAdaptor;
use ReseqTrack::RunMetaInfo;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::RunMetaInfoUtils qw (are_run_meta_infos_identical);
#use ReseqTrack::Tools::ERAUtils qw (convert_population convert_center_name);

use File::Basename;

@ISA = qw(ReseqTrack::DBSQL::BaseAdaptor);

sub new {
  my ($class, $db) = @_;

  my $self = $class->SUPER::new($db);

  return $self;
}

sub table_name{
  return "era.g1k_sequence_index";
}

sub columns{
  return "run_id, study_id, study_name, center_name, submission_id, submission_date, ".
      "sample_id, sample_name, population, experiment_id, instrument_platform, ".
      "instrument_model, library_name, run_name, run_block_name, paired_length, ".
      "library_layout, run_file_name, status, read_count, base_count,err_fastq_available";
}

sub fetch_by_dbID{
  my ($self) = @_;
  throw("Can't call fetch_by_dbID on ".$self." doesn't make sense");
}


sub fetch_by_run_id{
  my ($self, $run_id) = @_;
  my $sql = "select ".$self->columns." from ".$self->table_name. 
      " where run_id = ?";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $run_id);
  $sth->execute;
  my $rowHashref = $sth->fetchrow_hashref;
  print $sql."\n";
  if(!$rowHashref || keys(%$rowHashref) == 0){
    #print "Failed to find anything for ".$run_id."\n";
    return undef;
  }
  print "Have run id ".$rowHashref->{RUN_ID}." submission date ".$rowHashref->{SUBMISSION_DATE}."\n";
  my $object = $self->object_from_hashref($rowHashref);
  $sth->finish;
  return $object;
}

sub fetch_by_run_file_name{
  my ($self, $run_id) = @_;
  my $sql = "select ".$self->columns." from ".$self->table_name. 
      " where run_file_name = ?";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $run_id);
  $sth->execute;
  my $rowHashref = $sth->fetchrow_hashref;
  if(!$rowHashref || keys(%$rowHashref) == 0){
    #print "Failed to find anything for ".$run_id."\n";
    return undef;
  }
  my $object = $self->object_from_hashref($rowHashref);
  $sth->finish;
  return $object;
}

sub fetch_all_between_dates{
  my ($self, $start_day, $start_month, $start_year, 
      $end_day, $end_month, $end_year) = @_;
  my $sql = "select ".$self->columns." from ".$self->table_name;
  my $after_clause;
  my $before_clause;
  if($start_month && $start_year){
    $start_day = 31 if(!$start_day);
    $after_clause = "submission_date > TO_DATE('".$start_day."-".$start_month."-".$start_year."', 'DD-MM-YYYY')";
  }
  if($end_month && $end_year){
    $end_day = 1 if(!$end_day);
    $before_clause = "submission_date < TO_DATE('".$end_day."-".$end_month."-".$end_year."', 'DD-MM-YYYY')";
  }
  $sql .= " where ";
  $sql .= $after_clause." " if($after_clause);
  $sql .= "and " if($after_clause && $before_clause);
  $sql .= $before_clause if($before_clause);
  print $sql."\n";
  my $sth = $self->prepare($sql);
  print "Have statement ".$sth."\n";
  $sth->execute;
  print "Have executed statement\n";
  my @index;
  while(my ($rowHashref) = $sth->fetchrow_hashref){
    my $index = $self->object_from_hashref($rowHashref);
    push(@index, $index);
  }
  $sth->finish;
  return \@index;
}


sub object_from_hashref{
  my ($self, $hashref) = @_;
  throw("Can't create a object from an undefined hashref") 
      if(!$hashref || keys(%$hashref) == 0);
   my $new_date = $self->convert_date($hashref->{SUBMISSION_DATE}, $hashref->{RUN_ID}) if($hashref->{SUBMISSION_DATE});
  if(!$new_date && $hashref->{SUBMISSION_DATE}){
    throw("Failed to convert ".$hashref->{SUBMISSION_DATE});
  }
  my $population = convert_population($hashref->{POPULATION}, $hashref->{RUN_ID}, $hashref->{STUDY_ID});
  my $center_name = convert_center_name($hashref->{CENTER_NAME});
  my $object = ReseqTrack::RunMetaInfo->new(
    -run_id => $hashref->{RUN_ID},
    -study_id => $hashref->{STUDY_ID},
    -study_name => $hashref->{STUDY_NAME},
    -center_name => $center_name,
    -submission_id => $hashref->{SUBMISSION_ID},
    -submission_date => $new_date,
    -sample_id => $hashref->{SAMPLE_ID},
    -sample_name => $hashref->{SAMPLE_NAME},
    -population => $population,
    -experiment_id => $hashref->{EXPERIMENT_ID},
    -instrument_platform => $hashref->{INSTRUMENT_PLATFORM},
    -instrument_model => $hashref->{INSTRUMENT_MODEL},
    -run_name => $hashref->{RUN_NAME},
    -run_block_name => $hashref->{RUN_BLOCK_NAME},
    -library_name => $hashref->{LIBRARY_NAME},
    -paired_length => $hashref->{PAIRED_LENGTH},
    -library_layout => $hashref->{LIBRARY_LAYOUT},
    -run_file_name => $hashref->{RUN_FILE_NAME},
    -status => $hashref->{STATUS},
    -archive_read_count => $hashref->{READ_COUNT},
    -archive_base_count => $hashref->{BASE_COUNT},
    -err_fastq_available=> $hashref->{ERR_FASTQ_AVAILABLE},
      );  
   return $object;
}


sub convert_date{
  my ($self, $old_date, $run_id) = @_;
  #2009-NOV-16  21:39:50
  my $verbose = 0;
  $verbose = 1 if($run_id eq 'SRR015450');
  print "Converting ".$old_date."\n" if($verbose);
  return undef unless($old_date);
  my ($date, $time) = split /\s+/, $old_date;

  $time = '00:00:00' unless($time);
  my ($day, $month, $year) = split /\-/, $date;
  my $date_hash = date_hash();
  if($year =~ /\d\d/ && ( length ($year) == 2) ){
    $year = "20".$year
  }
  my $num_month = $date_hash->{$month};
  if(!$num_month){
    #foreach my $key(keys(%$date_hash)){
    #  print "Month ".$key." ".$date_hash->{$key}."\n";
    #}
    throw("Failed to find a month for ".$month);
  }
  my $new_date = $year."-".$num_month."-".$day;
  $new_date .= " ".$time if($time);
  print "Returning ".$new_date."\n" if($verbose);
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
sub convert_center_name{
  my ($name) = @_;
  if($name && $name eq 'ABHTD'){
    return 'ABI';
  }else{
    return $name;
  }
}

sub convert_population{
  my ($string, $run_id, $study_id) = @_;
  my $pop;
  if($string =~ /yri/i){
    $pop = 'YRI';
  }elsif($string =~ /yoruba/i){
    $pop = 'YRI';
  }elsif($string =~ /han chinese/i){
    $pop = 'CHB';
  }elsif($string =~ /japan/i){
    $pop = 'JPT';
  }elsif($string =~ /CEU/i){
    $pop = 'CEU'; 
  }elsif($string =~ /CEPH/i){
    $pop = 'CEU';
  }elsif($string =~ /tuscan/i){
    $pop = 'TSI';
  }elsif($string =~ /denver/i){
    $pop = 'CHD';
  }elsif($string =~ /Luhya/i){
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
  }elsif($string =~ /FIN/i){
    $pop = 'FIN';
  }elsif($string =~ /SHC/){
    $pop = 'CHS';
  }else{
    throw("Failed to find pop for ".$string." ".$run_id." ".$study_id);
   }
  return $pop;
}

sub store{
  throw("ReseqTrack::DBSQL::ERARunMetaInfoAdaptor can't store information");
}

1;
