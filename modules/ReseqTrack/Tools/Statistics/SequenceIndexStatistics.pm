=pod

=head1 NAME

ReseqTrack::Tools::Statistics::SequenceIndexStatistics

=head1 SYNOPSIS

This is an object to provide some basic sequence index statistics

=head1 Example

my $sequence_index_stats = ReseqTrack::Tools::Statistics::SequenceIndexStatistics
  ->new(
	-db => $db,
	-new_index => $new_index_file,
	-old_index => $old_index_file,
       );

$sequence_index_stats->make_stats;
$sequence_index_stats->print_stats;


=cut

package ReseqTrack::Tools::Statistics::SequenceIndexStatistics;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::SequenceIndexUtils;
use ReseqTrack::Tools::ERAUtils;
use ReseqTrack::Tools::GeneralUtils;
use ReseqTrack::Tools::Statistics;

@ISA = qw(ReseqTrack::Tools::Statistics);

=head2 new

  Arg [1]   : ReseqTrack::Tools::Statistics:SequenceIndexStatistics
  Arg [NEW_INDEX]   : string, path to sequence.index file
  Arg [OLD_INDEX]   : string, path to sequence index file
  Function  : create ReseqTrack::Tools::Statistics::SequenceIndexStatistics object
  Returntype: ReseqTrack::Tools::Statistics:SequenceIndexStatistics
  Exceptions: n/a
  Example   : my $stats = ReseqTrack::Tools::Statistics::SequenceIndexStatistics
  ->new(
	-db => $db,
	-new_index => $new_index_file,
	-old_index => $old_index_file,
       );

=cut



sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  unless($self->new_index && $self->old_index){
    $self->fetch_index_files();
  }
  return $self;
}

#Accessor methods




=head2 new/old_index

  Arg [1]   : ReseqTrack::Tools::Statistics::SequenceIndexStatistics
  Arg [2]   : string, filepath
  Function  : accessor method for index file paths
  Returntype: string
  Exceptions: throws if file doesn't exist or if filename doesn't match the
  expected patter YYYYMMDD.sequence.index
  Example   :

=cut



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





=head2 hashref_holders

  Arg [1]   : ReseqTrack::Tools::Statistics::SequenceIndexStatistics
  Arg [2]   : hashref
  Function  : accessor methods for several hashes used to hold stats
  Returntype: hashref
  Exceptions: throws if not passed a hashref
  Example   : my $new_per_sample = $self->new_per_sample();

=cut



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


=head2 fetch_index_files

  Arg [1]   : ReseqTrack::Tools::Statistics::SequenceIndexStatistics
  Function  : fetch the most recent index file and the previous verion from
  the given database
  Returntype: n/a
  Exceptions: throw if unable to define a new_index file or if index file name
  doesn't match the expected patter
  Example   : $self->fetch_index_files

=cut



sub fetch_index_files{
  my ($self) = @_;
  my $fa = $self->db->get_FileAdaptor;
  #fetch all files of type INDEX
  my $files = $fa->fetch_by_type('INDEX');
  my %index_files;
  foreach my $file(@$files){
    #skip any alignment index
    next if($file->filename =~ /alignment/);
    #skip the root sequence index file
    next if($file->filename eq 'sequence.index');
    $file->filename =~ /(\d+)\.sequence\.index/;
    if(!$1){
      throw("Don't know how to parse ".$file->name);
    }
    #associate the date with a particular file
    $index_files{$1} = $file->name;
  }
  #sort dates numerically, the newest date should always be the largest number
  #provided the YYYYMMDD format is followed
  my @dates = sort {$a <=> $b} keys(%index_files);
  $self->new_index =~ /(\d+)\.sequence\.index/ if($self->new_index);
  my $new_date = $1;
  my $old_date;
  #The newest date is the new file unless it is already defined
  $new_date = $dates[-1] unless($new_date);
  unless($self->new_index){
    $self->new_index($index_files{$new_date});
  }
  #Iterate through the other dates, once the new date is matched take the previous
  #date to defind the old index
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


=head2 make_stats

  Arg [1]   : ReseqTrack::Tools::Statistics::SequenceIndexStatistics
  Function  : create the statistics given thes specified index files
  Returntype: n/a
  Exceptions: n/a
  Example   : $self->make_stats;

=cut



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



=head2 parse_index

  Arg [1]   : ReseqTrack::Tools::Statistics::SequenceIndexStatistics
  Arg [2]   : string, path to index file
  Function  : to parse the index files to produce base pair statistics
  Returntype: array of hashrefs
  Exceptions: throw if a run doesn't have a base count
  Example   :

=cut



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
    if($values[5] eq 'ABHTD'){
      $values[5] = 'ABI';
    }
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


=head2 calculate_summary_stats

  Arg [1]   : ReseqTrack::Tools::Statistics::SequenceIndexStatistics
  Function  : create diffs between given stats and populate a hash with the 
  statistics
  Returntype:hashref
  Exceptions: n/a
  Example   : $stats_hash = $self->calculate_summary_stats

=cut



sub calculate_summary_stats{
  my ($self) = @_;

  #get summary hashes
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
  #number of samples > 10Gb ...    was number of samples > 4x
  #total gb
  #gb per population
  #gb per platform
  #gb per center
  my %stats_hash;
  #Define new stats and convert into gigabases
  $stats_hash{new}{Date} = $self->get_index_date($self->new_index);
  $stats_hash{new}{'# Accessions'} = keys(%$new_run);
  $stats_hash{new}{'# Samples'} = keys(%$new_sample);
  foreach my $sample(keys(%$new_sample)){
    #samples which are greater than 4x have more than 12Gb of sequence
    $stats_hash{new}{'# Samples greater than 10Gb'}++ 
      if($new_sample->{$sample} > 10000000000);
  }
  foreach my $pop(keys(%$new_population)){
    $stats_hash{new}{'Population in Gb'}{$pop} = convert_to_giga($new_population->{$pop});
  }
  foreach my $platform(keys(%$new_platform)){
     $stats_hash{new}{'Platform in Gb'}{$platform} = 
       convert_to_giga($new_platform->{$platform});
  }
  foreach my $center(keys(%$new_center)){
    $stats_hash{new}{'Center in Gb'}{$center} =
      convert_to_giga($new_center->{$center});
  }
  if($self->old_index){
    $stats_hash{old}{Date} = $self->get_index_date($self->old_index);
    $stats_hash{old}{'# Accessions'} = keys(%$old_run);
    $stats_hash{old}{'# Samples'} = keys(%$old_sample);
    foreach my $sample(keys(%$old_sample)){
      $stats_hash{old}{'# Samples greater than 10Gb'}++ 
	if($old_sample->{$sample} > 10000000000);
    }
    foreach my $pop(keys(%$old_population)){
      $stats_hash{old}{'Population in Gb'}{$pop} = convert_to_giga($old_population->{$pop});
    }
    foreach my $platform(keys(%$old_platform)){
      $stats_hash{old}{'Platform in Gb'}{$platform} = 
	convert_to_giga($old_platform->{$platform});
    }
    foreach my $center(keys(%$old_center)){
      $stats_hash{old}{'Center in Gb'}{$center} =
	convert_to_giga($old_center->{$center});
    }
    #Calculate diffs
    $stats_hash{diff}{'Date'} = current_date;
    $stats_hash{diff}{'# Accessions'} =
      ($stats_hash{new}{'# Accessions'} - $stats_hash{old}{'# Accessions'});
    $stats_hash{diff}{'# Samples'} =
      ($stats_hash{new}{'# Samples'} - $stats_hash{old}{'# Samples'});
    $stats_hash{diff}{'# Samples greater than 10Gb'} =
      ($stats_hash{new}{'# Samples greater than 10Gb'} - 
       $stats_hash{old}{'# Samples greater than 10Gb'});
    foreach my $pop(keys(%$new_population)){
      my $new_pop = $stats_hash{new}{'Population in Gb'}{$pop};
      my $old_pop = $stats_hash{old}{'Population in Gb'}{$pop};
      $old_pop = 0 unless($old_pop);
      $stats_hash{diff}{'Population in Gb'}{$pop} = $new_pop - $old_pop;
    }
    foreach my $plat(keys(%$new_platform)){
      my $new_plat = $stats_hash{new}->{'Platform in Gb'}->{$plat};
      my $old_plat = $stats_hash{old}->{'Platform in Gb'}->{$plat};
      $old_plat = 0 unless($old_plat);
      $stats_hash{diff}{'Platform in Gb'}{$plat} = $new_plat - $old_plat;
    }
    foreach my $center(keys(%$new_center)){
	my $new_pop = $stats_hash{new}->{'Center in Gb'}->{$center};
	my $old_pop = $stats_hash{old}->{'Center in Gb'}->{$center};
	$old_pop = 0 unless($old_pop);
	$stats_hash{diff}{'Center in Gb'}{$center} = $new_pop - $old_pop;
    }
  }
  $self->stats_hash(\%stats_hash);
  return \%stats_hash;
}



=head2 print_stats

  Arg [1]   : ReseqTrack::Tools::Statistics::SequenceIndexStatistics
  Arg [2]   : hashref, hash of stats as produced by calculate summary stats
  Function  : to print to STDOUT a csv file with stats
  Returntype: n/a
  Exceptions: throws if doesn't not how to print a particular row of stats
  Example   : $self->print_stats();

=cut


sub print_stats{
  my ($self, $stats_hash) = @_;

  $stats_hash = $self->stats_hash unless($stats_hash);
  $stats_hash = $self->calculate_summary_stats() unless($stats_hash);

  my @headers = ('old', 'new', 'diff');
  my @rows = ('Date', '# Accessions', '# Samples', '# Samples greater than 10Gb', 
	      'Population in Gb', 'Platform in Gb', 'Center in Gb');
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
      throw("Not sure how to print ".$row." ".$type);
    }
    print "\n";
  }
}




=head2 row_type

  Arg [1]   : ReseqTrack::Tools::Statistics::SequenceIndexStatistics
  Function  : return a hash which defined data types for particular rows of
  statistics, can be either Scalar or HASH
  Returntype: hashref
  Exceptions: 
  Example   : 

=cut



sub row_type{
  my ($self) = @_;
  my %hash;
  $hash{Date} = 'SCALAR';
  $hash{'# Accessions'} = 'SCALAR';
  $hash{'# Samples'} = 'SCALAR';
  $hash{'# Samples greater than 10Gb'} = 'SCALAR';
  $hash{'Population in Gb'} = 'HASH';
  $hash{'Platform in Gb'} = 'HASH';
  $hash{'Center in Gb'} = 'HASH';
  return \%hash;
}


=head2 get_index_date

  Arg [1]   : ReseqTrack::Tools::Statistics::SequenceIndexStatistics
  Arg [2]   : string, file path
  Function  : parses date from index file name assuming format
  YYYYMMDD.sequence.index
  Returntype: string, date string
  Exceptions: n/a
  Example   : my $date = $self->get_index_date($index_path);

=cut



sub get_index_date{
  my ($self, $index) = @_;
 # print STDERR "Parsing ".$index."\n";
  $index =~ /(\d+)\.sequence\.index/;
  my $date = $1;
  # print STDERR "RETURNING ".$date."\n";
  return $date;
}

1;
