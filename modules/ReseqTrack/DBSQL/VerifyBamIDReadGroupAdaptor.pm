package ReseqTrack::DBSQL::VerifyBamIDReadGroupAdaptor;
  
use strict;
use warnings;

use ReseqTrack::DBSQL::BaseAdaptor;
use ReseqTrack::Base;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename;
use Data::Dumper;
use ReseqTrack::VerifyBamIDReadGroup;


use vars qw(@ISA);
@ISA = qw(ReseqTrack::DBSQL::BaseAdaptor);

sub new {
  my ( $class, $db ) = @_;

  my $self = $class->SUPER::new($db);

  return $self;
}

sub columns {
  return "
        verifybamid_readgroup.verifybamid_readgroup_id,
        verifybamid_readgroup.other_id,
        verifybamid_readgroup.run_id,
        verifybamid_readgroup.selfibd, 
        verifybamid_readgroup.selfmix, 
        verifybamid_readgroup.best_sample,
        verifybamid_readgroup.bestibd,
        verifybamid_readgroup.bestmix,
        verifybamid_readgroup.status
    ";
}

sub table_name {
  return "verifybamid_readgroup";
}

sub store {
  my ( $self, $verifybamid_readgroup ) = @_;

  throw(  "Can't store "
	  . $verifybamid_readgroup
	  . " using ReseqTrack::DBSQL::VerifyBamIDReadGroupAdaptor" )
    unless (
	    $verifybamid_readgroup->isa(
					"ReseqTrack::VerifyBamIDReadGroup")
	   );

  my $sql =
    "insert into verifybamid_readgroup (other_id, "
      . "run_id, selfibd, selfmix ,best_sample,bestibd,bestmix,status) "
	. " values( ?,?,?,?,?,?,?,?) ";

  my $sth = $self->prepare($sql);

  $sth->bind_param( 1, $verifybamid_readgroup->other_id );
  $sth->bind_param( 2, $verifybamid_readgroup->run_id );
  $sth->bind_param( 3, $verifybamid_readgroup->selfibd );
  $sth->bind_param( 4, $verifybamid_readgroup->selfmix );
  $sth->bind_param( 5, $verifybamid_readgroup->best_sample );
  $sth->bind_param( 6, $verifybamid_readgroup->bestibd );
  $sth->bind_param( 7, $verifybamid_readgroup->bestmix );
  $sth->bind_param( 8, $verifybamid_readgroup->status );

  my $rows_inserted = $sth->execute();
  my $dbID          = $sth->{'mysql_insertid'};

  $verifybamid_readgroup->dbID($dbID);
  $verifybamid_readgroup->adaptor($self);
  $sth->finish();

  return $verifybamid_readgroup;
}

sub update {
  my ( $self, $verifybamid_readgroup ) = @_;

  throw(  "Can't store "
	  . $verifybamid_readgroup
	  . " using ReseqTrack::DBSQL::VerifyBamIDReadGroupAdaptor" )
    unless (
	    $verifybamid_readgroup->isa(
					"ReseqTrack::VerifyBamIDReadGroup"
				       )
	   );

  my $sql =
    "update verifybamid_readgroup set "
      . " selfibd = ? , selfmix = ? , best_sample = ?, "
	. "bestibd  = ? ,bestmix = ?, status = ? where  other_id = ?  and run_id = ?";

 
  
  my $sth = $self->prepare($sql);

  $sth->bind_param( 1, $verifybamid_readgroup->selfibd );
  $sth->bind_param( 2, $verifybamid_readgroup->selfmix );
  $sth->bind_param( 3, $verifybamid_readgroup->best_sample );
  $sth->bind_param( 4, $verifybamid_readgroup->bestibd );
  $sth->bind_param( 5, $verifybamid_readgroup->bestmix );
  $sth->bind_param( 6, $verifybamid_readgroup->status );
  $sth->bind_param( 7, $verifybamid_readgroup->other_id );
  $sth->bind_param( 8, $verifybamid_readgroup->run_id );

  $sth->execute();
  $sth->finish();

  return $verifybamid_readgroup;
}



sub object_from_hashref {
  my ( $self, $hashref ) = @_;

  throw(
	" Can't create a ReseqTrack::VerifyBamIDReadGroup from an undefined hashref
      "
       ) if ( !$hashref );

  my $OBJ = "ReseqTrack::VerifyBamIDReadGroup";
  my $obj = $OBJ->new(

		      -adaptor               => $self,
		      -dbID                  => $hashref->{verifybamid_readgroup_id},
		      -other_id              => $hashref->{other_id},
		      -run_id                => $hashref->{run_id},
		      -selfibd               => $hashref->{selfibd},
		      -selfmix               => $hashref->{selfmix},
		      -best_sample           => $hashref->{best_sample},
		      -bestibd               => $hashref->{bestibd},
		      -bestmix               => $hashref->{bestmix},
		      -status                => $hashref->{status},
		     );
  return $obj;
}


sub fetch_by_other_id_and_run_id{
    
  my ($self, $other_id, $run_id) = @_;
  my $sql = "select ". $self->columns." from " .  $self->table_name .
    " where other_id = ? and run_id = ?";
      
   
     
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $other_id);
  $sth->bind_param(2, $run_id);
  $sth->execute;
  
  my @results;
  while (my $rowHashref = $sth->fetchrow_hashref) {
    my $result = $self->object_from_hashref($rowHashref);
    push(@results, $result);
  }
  $sth->finish;

  throw($other_id." has returned multiple objects ".@results." not sure what to do") 
    if (@results && @results >= 2);
 

  if ( @results ==0 ){
   print "No results found for  $other_id $run_id\n";
   return;

  }

  my $result = $results[0];
    
  return $result;
}


sub fetch_by_status{
    
  my ($self, $status) = @_;
  my $sql = "select run_id from " .  $self->table_name .
    " where status = ?";
     
       
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $status);
  $sth->execute;
  
  my %failed_run_ids;
  my @results;
  
  while (my $rowHashref = $sth->fetchrow_hashref) {
    my $result = $self->object_from_hashref($rowHashref);
    push(@results, $result);
    $failed_run_ids{$result->run_id}{run_status} = 1;    
  }
  $sth->finish;
 
  

#  print "Got $ctr run_ids that failed verifybamid\n";
  
  return (\@results ,\%failed_run_ids);
}


sub fetch_invalid_run_id_designation{
    
  my ($self) = @_;
  my $sql =  "select distinct other_id from " .  $self->table_name . 
    " where run_id not like \"%RR%\"";

  my $ctr = 0;    
      
  my $sth = $self->prepare($sql);
  $sth->execute;
  
  my @results;

  while (my $rowArrayref = $sth->fetchrow_arrayref) {
    $ctr++;
    my $result = @$rowArrayref[0];
    push(@results, $result);
  }
  $sth->finish;

#  print "Got $ctr run_ids with non \"*RR*\" run_ids\n";
  
  return \@results;
}


sub fetch_by_other_id {
    
  my ($self, $other_id) = @_;
  my $sql = "select ". $self->columns." from " .  $self->table_name .
    " where other_id = ? ";
      
   
     
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $other_id);
  $sth->execute;
  
  my @results;
  while (my $rowHashref = $sth->fetchrow_hashref) {
    my $result = $self->object_from_hashref($rowHashref);
    push(@results, $result);
  }
  $sth->finish;

#  print "$other_id: has ".@results." read group results \n"; 
  
  return \@results;
}

sub fetch_by_run_id {
    
  my ($self, $run_id) = @_;
  my $sql = "select ". $self->columns." from " .  $self->table_name .
    " where run_id = ? ";
      
   
     
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $run_id);
  $sth->execute;
  
  my @results;
  while (my $rowHashref = $sth->fetchrow_hashref) {
    my $result = $self->object_from_hashref($rowHashref);
    push(@results, $result);
  }
  $sth->finish;

#  print "$run_id: has ".@results." read group results \n"; 
  
  return \@results;
}


sub update_run_id {

  my ( $self, $verifybamid_readgroup ) = @_;
  
  throw(  "Can't store "
	  . $verifybamid_readgroup
	  . " using ReseqTrack::DBSQL::VerifyBamIDReadGroupAdaptor" )
    unless (
	    $verifybamid_readgroup->isa(
					"ReseqTrack::VerifyBamIDReadGroup"
				       )
	   );

  my $sql = "update verifybamid_readgroup set  run_id = ? ".
    " where verifybamid_readgroup_id  = ?";
	        
  my $sth = $self->prepare($sql);

  $sth->bind_param( 1, $verifybamid_readgroup->run_id );
  $sth->bind_param( 2, $verifybamid_readgroup->dbID );

  $sth->execute();
  $sth->finish();

  return $verifybamid_readgroup;
}

1;

