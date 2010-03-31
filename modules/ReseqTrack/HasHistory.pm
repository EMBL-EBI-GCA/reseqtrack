package ReseqTrack::HasHistory;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);

use ReseqTrack::Base;

@ISA = qw(ReseqTrack::Base);

sub new {
  my ($class, @args) = @_;
  my  $self = $class->SUPER::new(@args);
  my ($history, $statistics) = rearrange(['HISTORY', 'STATISTICS'], @args);
  $self->history($history);
  $self->statistics($statistics);
  return $self;
}

=head2 history

  Arg [1]   : ReseqTrack::HasHistory
  Arg [2]   : ReseqTrack::History or array ref of ReseqTrack::History
  Function  : add history objects to internal arrayref. If no history objects
  are defined but a dbID and adaptor is the method will try and fetch any history 
  objects
  Returntype: arrayref of ReseqTrack::History
  Exceptions: throw if not passed a ReseqTrack::History object or arrayref of said
  objects
  Example   : 

=cut



sub history{
  my ($self, $arg) = @_;
  if($arg){
    $self->populate_history;
    if(ref($arg) eq 'ARRAY'){
      throw("Must pass ReseqTrack::HasHistory::history an arrayref of history objects") if($arg->[0] && (!$arg->[0]->isa("ReseqTrack::History")));
      unless($arg->[0]){
        throw($arg." appears to contain an undefined entry");
      }
      push(@{$self->{history}}, @$arg);
    }elsif($arg->isa("ReseqTrack::History")){
      push(@{$self->{history}}, $arg);
    }else{
      throw("Must give ReqseqTrack::File::history either a ReseqTrack::History ".
            "object or an arrayref of History objects not ".$arg);
    }
    $self->{history} = $self->uniquify_histories($self->{history});
  }
  
  if(!$self->{history} || @{$self->{history}} == 0){
    $self->populate_history;
  }
  
  return $self->{history};
}



=head2 populate_history

  Arg [1]   : ReseqTrack::HasHistory
  Function  : gets history objects from the database and ensures they are all unique
  Returntype: n/a
  Exceptions: warns if non unique history objects are found
  Example   : 

=cut


sub populate_history{
  my ($self) = @_;
  my @histories;
  if($self->adaptor && $self->dbID){
    my $hist_a = $self->adaptor->db->get_HistoryAdaptor;
    my $objects = $hist_a->fetch_by_other_id_and_table_name($self->dbID, 
                                                            $self->object_table_name);
    push(@histories, @$objects) if($objects && @$objects >= 1);
  }
  push(@histories, @{$self->{history}}) if($self->{history});
  $self->{history} = $self->uniquify_histories(\@histories);
}

sub uniquify_histories{
  my ($self, $histories) = @_;
  my %hash;
 HISTORY:foreach my $history(@$histories){
   my $comment = $history->comment;
   my $time = $history->time;
   $time = '' if(!$time);
   my $unique_string = $comment."-".$time;
   if($hash{$unique_string}){
     warning($unique_string." is already defined skipping");
     next HISTORY;
   }
   $hash{$unique_string} = $history;
 }
  my @values = values(%hash);
  return \@values;
}


=head2 refresh_history

  Arg [1]   : ReseqTrack::File
  Function  : This will overwrite the existing history array and pull a new one
  from the database
  Returntype: arrayref of ReseqTrack::History
  Exceptions: warns if array current contains objects
  Example   : 

=cut


sub refresh_history{
  my ($self) = @_;
  if($self->adaptor && $self->dbID){
    my $hist_a = $self->adaptor->get_HistoryAdaptor;
    my $objects = $hist_a->fetch_by_object_id_and_table_name($self->dbID, 
                                                             $self->object_table_name);
    if(@{$self->{history}}){
      warning($self."->{history} contains objects you will lose them now");
    }
    $self->{history} = $objects;
  }
  return $self->{history};
}

sub statistics{
  my ($self, $arg) = @_;
  if($arg){
    $self->populate_statistics(1); # $now has two stats objects, read and base
    if(ref($arg) eq 'ARRAY'){
      if(@$arg >= 1){
        throw("Must pass ReseqTrack::HasHistory::statistics an arrayref of statistics objects")
            unless($arg->[0]->isa("ReseqTrack::Statistic"));
        push(@{$self->{statistics}}, @$arg);
      }
    }elsif($arg->isa("ReseqTrack::Statistic")){
      push(@{$self->{statistics}}, $arg);
    }else{
      throw("Must give ReqseqTrack::File::statistics either a ReseqTrack::History ".
            "object or an arrayref of History objects not ".$arg);
    }
    $self->{statistics} = $self->uniquify_statistics($self->{statistics});
  }
  if(!$self->{statistics} || @{$self->{statistics}} == 0){
    $self->populate_statistics;
  }
  #print "at the end of statistics block in HasHistory.pm, att value is " . $self->{statistics}->[0]->attribute_value . "\n" if ($self->{statistics}->[0]);  
  #print "at the end of statistics block in HasHistory.pm, att value is " . $self->{statistics}->[1]->attribute_value . "\n" if ($self->{statistics}->[1]);
  #print "at the end of statistics block in HasHistory.pm, att value is " . $self->{statistics}->[2]->attribute_value . "\n" if ($self->{statistics}->[2]);
  #print "at the end of statistics block in HasHistory.pm, dbID value is " . $self->{statistics}->[0]->dbID . "\n" if ($self->{statistics}->[0]);
  #print "at the end of statistics block in HasHistory.pm, dbID value is " . $self->{statistics}->[1]->dbID . "\n" if ($self->{statistics}->[1]);
  return $self->{statistics};
}

sub populate_statistics{
  my ($self, $dont_unique) = @_;
  my @statistics;
  if($self->adaptor && $self->dbID){
    my $hist_a = $self->adaptor->db->get_StatisticsAdaptor;
    my $objects = $hist_a->fetch_by_other_id_and_table_name($self->dbID, 
                                                            $self->object_table_name);
    push(@statistics, @$objects) if($objects && @$objects >= 1);
    #print "Inside populate_stats, att v is " . $$objects[0]->attribute_value . "\n" if (@statistics );
    #print "Inside populate_stats, dbID is " . $$objects[0]->dbID . "\n" if (@statistics );
  }
  push(@statistics, @{$self->{statistics}}) if($self->{statistics});
  unless($dont_unique){
    $self->{statistics} = $self->uniquify_statistics(\@statistics);
  }else{
    $self->{statistics} = \@statistics;
  }
}

sub refresh_statistics{
  my ($self) = @_;
  if($self->adaptor && $self->dbID){
    my $hist_a = $self->adaptor->get_StatisticsAdaptor;
    my $objects = $hist_a->fetch_by_object_id_and_table_name
        ($self->dbID, $self->object_table_name);
    if(@{$self->{statistics}}){
      warning($self."->{statistics} contains objects you will lose them now");
    }
    $self->{statistics} = $objects;
  }
  return $self->{statistics};
}

sub uniquify_statistics{
  my ($self, $statistics) = @_;
  my %hash;
  my %id; #key is obj att name, value is its dbID
  # This function takes an array of stats objects for the same file, 
  # merge stats with the same attribute name and update the  attribute value.
  # some stats objects have not been loaded in db so they do not have dbID, in these cases, dbID would be inherited from existing tats object
  HISTORY:foreach my $stats(@$statistics){
  if($hash{$stats->attribute_name}){
      if ($stats->dbID) {
          $id{$stats->attribute_name} = $stats->dbID;
      }
      else {
	  $stats->dbID($id{$stats->attribute_name});
      }	
      $hash{$stats->attribute_name} = $stats;#this way the same attribute name always gets the latest value
      warning($stats->attribute_name." is already defined, the attribute value is overwritten with the new value");      
      next HISTORY;
   }
   else{
      $hash{$stats->attribute_name} = $stats;
      $id{$stats->attribute_name}=$stats->dbID;
   }
  }
  my @values = values(%hash);
  #for my $i (@values) {
  #	print "after uniquify, stats obj have dbID: ". $i->dbID ."\n";
  #}
  return \@values;
}

sub replace_statistic{
  my ($self, $statistic) = @_;
  my %hash;
  foreach my $statistics(@{$self->statistics}){
    $hash{$statistics->attribute_name} = $statistics;#this way the same attribute name always gets the latest value
  }
  
  $hash{$statistic->attribute_name} = $statistic;
  my @values = values(%hash);
  return \@values;
}

sub object_table_name{
  my ($self) = @_;
  throw($self." must implement an object_table_name method to name the table ".
        "its data is stored in");
}

1;
