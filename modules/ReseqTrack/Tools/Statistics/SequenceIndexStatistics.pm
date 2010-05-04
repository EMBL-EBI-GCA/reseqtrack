=pod

=head1 NAME

ReseqTrack::Tools::Statistics::SequenceIndexStatistics

=head1 SYNOPSIS

This is an object to provide some basic sequence index statistics

=head1 Example


=cut

package ReseqTrack::Tools::Statistics::SequenceIndexStatistics;

use strict;

use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::SequenceIndexUtils;
use ReseqTrack::Tools::ERAUtils;
use ReseqTrack::Tools::GeneralUtils;

sub new {
  my ($class, @args) = @_;
  my $self ={};
  bless $self,$class;
  my ($db, $new_index, $old_index, $skip_p2,
      $skip_p3) = rearrange([qw(DB NEW_INDEX OLD_INDEX SKIP_P2 SKIP_P3)], @args);

  # setting defaults
  $skip_p2 = 1 unless(defined($skip_p2));
  $skip_p3 = 1 unless(defined($skip_p3));
  #####
  
  $self->db($db);
  $self->new_index($new_index);
  $self->old_index($old_index);
  $self->skip_p2($skip_p2);
  $self->skip_p3($skip_p3);
  unless($self->new_index && $self->old_index){
    $self->fetch_index_files();
  }
  return $self;
}

#Accessor methods

sub db{
  my ($self, $db) = @_;
  if($db){
    throw("Must pass SequenceIndexStatistics:db a ReseqTrack::DBSQL::DBAdaptor not ".
	  " $db") unless($db->isa("ReseqTrack::DBSQL::DBAdaptor"));
    $self->{db} = $db;
  }
  return $self->{db};
}

sub new_index{
  my ($self, $new_index) = @_;
  if($new_index){
    throw("SequenceIndexStatistics:new_index ".$new_index." should exist")
      unless(-e $new_index);
    throw("SequenceIndexStatistcs:new_index ".$new_index." must match pattern ".
	  "YYYYMMDD.sequence.index") unless($new_index =~ /\d+\.sequence\.index/);
    $self->{new_index} = $new_index;
  }
  return $self->{new_index};
}

sub old_index{
  my ($self, $old_index) = @_;
  if($old_index){
    throw("SequenceIndexStatistics:old_index ".$old_index." should exist")
      unless(-e $old_index);
    throw("SequenceIndexStatistcs:new_index ".$old_index." must match pattern ".
	  "YYYYMMDD.sequence.index") unless($old_index =~ /\d+\.sequence\.index/);
    $self->{old_index} = $old_index;
  }
  return $self->{old_index};
}

sub skip_p2{
  my ($self, $skip_p2) = @_;
  if($skip_p2){
    $self->{skip_p2} = $skip_p2;
  }
  return $self->{skip_p2};
}

sub skip_p3{
  my ($self, $skip_p3) = @_;
  if($skip_p3){
    $self->{skip_p3} = $skip_p3;
  }
  return $self->{skip_p3};
}

sub new_per_sample{
  my ($self, $hashref) = @_;
  if($hashref){
    throw("Must pass new_per_sample a hashref not ".$hashref)
      unless(ref($hashref) eq 'HASH');
    $self->{new_per_sample} = $hashref;
  }
  return $self->{new_per_sample};
}

sub new_per_run{
  my ($self, $hashref) = @_;
  if($hashref){
    throw("Must pass new_per_run a hashref not ".$hashref)
      unless(ref($hashref) eq 'HASH');
    $self->{new_per_run} = $hashref;
  }
  return $self->{new_per_run};
}
sub new_per_population{
  my ($self, $hashref) = @_;
  if($hashref){
    throw("Must pass new_per_population a hashref not ".$hashref)
      unless(ref($hashref) eq 'HASH');
    $self->{new_per_population} = $hashref;
  }
  return $self->{new_per_population};
}

sub new_per_platform{
  my ($self, $hashref) = @_;
  if($hashref){
    throw("Must pass new_per_platform a hashref not ".$hashref)
      unless(ref($hashref) eq 'HASH');
    $self->{new_per_platform} = $hashref;
  }
  return $self->{new_per_platform};
}

sub new_per_center_population{
  my ($self, $hashref) = @_;
  if($hashref){
    throw("Must pass new_per_center_population a hashref not ".$hashref)
      unless(ref($hashref) eq 'HASH');
    $self->{new_per_center_population} = $hashref;
  }
  return $self->{new_per_center_population};
}

sub old_per_sample{
  my ($self, $hashref) = @_;
  if($hashref){
    throw("Must pass old_per_sample a hashref not ".$hashref)
      unless(ref($hashref) eq 'HASH');
    $self->{old_per_sample} = $hashref;
  }
  return $self->{old_per_sample};
}

sub old_per_run{
  my ($self, $hashref) = @_;
  if($hashref){
    throw("Must pass old_per_run a hashref not ".$hashref)
      unless(ref($hashref) eq 'HASH');
    $self->{old_per_run} = $hashref;
  }
  return $self->{old_per_run};
}

sub old_per_population{
  my ($self, $hashref) = @_;
  if($hashref){
    throw("Must pass old_per_population a hashref not ".$hashref)
      unless(ref($hashref) eq 'HASH');
    $self->{old_per_population} = $hashref;
  }
  return $self->{old_per_population};
}

sub old_per_platform{
  my ($self, $hashref) = @_;
  if($hashref){
    throw("Must pass old_per_platform a hashref not ".$hashref)
      unless(ref($hashref) eq 'HASH');
    $self->{old_per_platform} = $hashref;
  }
  return $self->{old_per_platform};
}

sub old_per_center_population{
  my ($self, $hashref) = @_;
  if($hashref){
    throw("Must pass old_per_center_population a hashref not ".$hashref)
      unless(ref($hashref) eq 'HASH');
    $self->{old_per_center_population} = $hashref;
  }
  return $self->{old_per_center_population};
}

sub stats_hash{
  my ($self, $hashref) = @_;
  if($hashref){
    throw("Must pass stats a hashref not ".$hashref)
      unless(ref($hashref) eq 'HASH');
    $self->{stats} = $hashref;
  }
  return $self->{stats_hash};
}

#Input methods

sub fetch_index_files{
  my ($self) = @_;
  my $fa = $self->db->get_FileAdaptor;
  my $files = $fa->fetch_by_type('INDEX');
  my %index_files;
  foreach my $file(@$files){
    next if($file->filename =~ /alignment/);
    next if($file->filename eq 'sequence.index');
    $file->filename =~ /(\d+)\.sequence\.index/;
    $index_files{$1} = $file->name;
  }
  my @dates = sort {$a <=> $b} keys(%index_files);
  $self->new_index =~ /(\d+)\.sequence\.index/ if($self->new_index);
  my $new_date = $1;
  my $old_date;
  $new_date = $dates[-1] unless($new_date);
  unless($self->new_index){
    $self->new_index($index_files{$new_date});
  }
  for(my $i=0; $i<@dates; $i++){
    next unless($new_date eq $dates[$i]);
    unless($i == 0){
      $old_date = $dates[$i-1];
    }
  }
  if($old_date && !$self->old_index){
    $self->old_index($index_files{$old_date});
  }
  throw("SequenceIndexStatistics:fetch_index_files Failed to define a new index ".
	"file") unless($self->new_index);
}

sub make_stats{
  my ($self) = @_;

  my ($run, $sample,$pop,$platform,$center) = $self->parse_index($self->new_index);
  $self->new_per_sample($sample);
  $self->new_per_population($pop);
  $self->new_per_platform($platform);
  $self->new_per_center_population($center);
  $self->new_per_run($run);
  if($self->old_index){
    my ($old_run, $old_sample,$old_pop,$old_platform,$old_center) = 
      $self->parse_index($self->old_index);
    $self->old_per_sample($old_sample);
    $self->old_per_population($old_pop);
    $self->old_per_platform($old_platform);
    $self->old_per_center_population($old_center);
    $self->old_per_run($old_run);
  }
}


sub parse_index{
  my ($self, $index) = @_;

  my $hash = get_index_hash($index);
  my %per_run;
  my %per_sample;
  my %per_population;
  my %per_platform;
  my %per_center_population;
  foreach my $file(keys(%{$hash})){
    my @values = split /\t/, $hash->{$file};
    next if($values[20]);
    my $pop = convert_population($values[10], $values[2], 0);
    next if($values[3] eq 'SRP000032' && $self->skip_p2);
    next if($values[3] eq 'SRP000033' && $self->skip_p3);
    unless($values[24] =~ /^\d+$/){
      throw($file." doesn't see to have a base count");
    }
    $per_run{$values[2]} += $values[24];
    $per_sample{$values[9]} += $values[24];
    $per_sample{total} += $values[24];
    $per_population{$pop} += $values[24];
    $per_population{'total'} += $values[24];
    $per_platform{$values[12]} += $values[24];
    $per_platform{'total'} += $values[24];
    $per_center_population{$values[5]} += $values[24];
    #$per_center_population{$values[5]}{total} += $values[24];
    $per_center_population{total} += $values[24];
  }
  return (\%per_run, \%per_sample, \%per_population, \%per_platform, 
	  \%per_center_population);
}

sub calculate_summary_stats{
  my ($self) = @_;
  my $new_run = $self->new_per_run;
  my $new_sample = $self->new_per_sample;
  my $new_population = $self->new_per_population;
  my $new_platform = $self->new_per_platform;
  my $new_center = $self->new_per_center_population;
  my ($old_sample, $old_population, $old_platform, $old_center, $old_run);
  if($self->old_index){
    $old_sample = $self->old_per_sample;
    $old_population = $self->old_per_population;
    $old_platform = $self->old_per_platform;
    $old_center = $self->old_per_center_population;
    $old_run = $self->old_per_run;
  }
  #date
  #number of accessions
  #number of samples
  #number of samples > 4x
  #total gb
  #gb per population
  #gb per platform
  #gb per center/population
  my %stats_hash;
  $stats_hash{new}{date} = $self->get_index_date($self->new_index);
  $stats_hash{new}{accessions} = keys(%$new_run);
  $stats_hash{new}{samples} = keys(%$new_sample);
  foreach my $sample(keys(%$new_sample)){
    $stats_hash{new}{samples_greater_than_4x}++ 
      if($new_sample->{$sample} > 12000000000);
  }
  foreach my $pop(keys(%$new_population)){
    $stats_hash{new}{population}{$pop} = convert_to_giga($new_population->{$pop});
  }
  foreach my $platform(keys(%$new_platform)){
     $stats_hash{new}{platform}{$platform} = 
       convert_to_giga($new_platform->{$platform});
  }
  foreach my $center(keys(%$new_center)){
    $stats_hash{new}{center}{$center} =
      convert_to_giga($new_center->{$center});
  }
  if($self->old_index){
    $stats_hash{old}{date} = $self->get_index_date($self->old_index);
    $stats_hash{old}{accessions} = keys(%$old_run);
    $stats_hash{old}{samples} = keys(%$old_sample);
    foreach my $sample(keys(%$old_sample)){
      $stats_hash{old}{samples_greater_than_4x}++ 
	if($old_sample->{$sample} > 12000000000);
    }
    foreach my $pop(keys(%$old_population)){
      $stats_hash{old}{population}{$pop} = convert_to_giga($old_population->{$pop});
    }
    foreach my $platform(keys(%$old_platform)){
      $stats_hash{old}{platform}{$platform} = 
	convert_to_giga($old_platform->{$platform});
    }
    foreach my $center(keys(%$old_center)){
      $stats_hash{old}{center}{$center} =
	convert_to_giga($old_center->{$center});
    }
    #diffs
    $stats_hash{diff}{date} = current_date;
    $stats_hash{diff}{accessions} =
      ($stats_hash{new}{accessions} - $stats_hash{old}{accessions});
    $stats_hash{diff}{samples} =
      ($stats_hash{new}{samples} - $stats_hash{old}{samples});
    $stats_hash{diff}{samples_greater_than_4x} =
      ($stats_hash{new}{samples_greater_than_4x} - 
       $stats_hash{old}{samples_greater_than_4x});
    foreach my $pop(keys(%$new_population)){
      my $new_pop = $stats_hash{new}{population}{$pop};
      my $old_pop = $stats_hash{old}{population}{$pop};
      $old_pop = 0 unless($old_pop);
      $stats_hash{diff}{population}{$pop} = $new_pop - $old_pop;
    }
    foreach my $plat(keys(%$new_platform)){
      my $new_plat = $stats_hash{new}->{platform}->{$plat};
      my $old_plat = $stats_hash{old}->{platform}->{$plat};
      $old_plat = 0 unless($old_plat);
      $stats_hash{diff}{platform}{$plat} = $new_plat - $old_plat;
    }
    foreach my $center(keys(%$new_center)){
	my $new_pop = $stats_hash{new}->{center}->{$center};
	my $old_pop = $stats_hash{old}->{center}->{$center};
	$old_pop = 0 unless($old_pop);
	$stats_hash{diff}{center}{$center} = $new_pop - $old_pop;
    }
  }
  $self->stats_hash(\%stats_hash);
  return \%stats_hash;
}

sub print_stats{
  my ($self, $stats_hash) = @_;

  $stats_hash = $self->stats_hash unless($stats_hash);
  $stats_hash = $self->calculate_summary_stats() unless($stats_hash);

  my @headers = ('old', 'new', 'diff');
  my @rows = ('date', 'accessions', 'samples', 'samples_greater_than_4x', 
	      'population', 'platform', 'center');
  my $type_hash = $self->row_type;
  print join(", ", "Category", @headers)."\n";;
  foreach my $row(@rows){
    my $type = $type_hash->{$row};
    if($type eq 'SCALAR'){
      print $row.", ";
      foreach my $header(@headers){
	my $value = $stats_hash->{$header}->{$row};
	throw($value. "This isn't a scalar it is a hash") if(ref($value) eq 'HASH');
	print $value.", ";
      }
      print "\n";
    }elsif($type eq 'HASH'){
      my @keys = keys(%{$stats_hash->{new}->{$row}});
      my @sorted_keys = sort {$a cmp $b} @keys;
      my %string;
      my $first = 0;
      print $row.", \n";
      #print join(", ", "\t", @keys)."\n";
      foreach my $key(@sorted_keys){
	my $string;
	$string = $key.", ";
	foreach my $header(@headers){
	  my $value = $stats_hash->{$header}->{$row}->{$key};
	  $value = 0 unless($value && $value >= 1);
	  $string .= $stats_hash->{$header}->{$row}->{$key}.", ";
	}
	$string{$key} = $string;
      }
      foreach my $key(@sorted_keys){
	print $string{$key}."\n";
      }
      print "\n";
    }else{
      throw("Not sure how to print ".$type);
    }
    print "\n";
  }
}


sub row_type{
  my ($self) = @_;
  my %hash;
  $hash{date} = 'SCALAR';
  $hash{accessions} = 'SCALAR';
  $hash{samples} = 'SCALAR';
  $hash{samples_greater_than_4x} = 'SCALAR';
  $hash{population} = 'HASH';
  $hash{platform} = 'HASH';
  $hash{center} = 'HASH';
  return \%hash;
}

sub get_index_date{
  my ($self, $index) = @_;
  $index =~ /(\d+)\.sequence\.index/;
  my $date = $1;
  return $date;
}

1;
