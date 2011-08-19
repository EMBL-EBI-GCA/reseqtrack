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
        verifybamid_ReadGroup.verifybamid_ReadGroup_id,
        verifybamid_ReadGroup.other_id,
        verifybamid_ReadGroup.run_id,
        verifybamid_ReadGroup.SELFIBD, 
        verifybamid_ReadGroup.SELFMIX, 
        verifybamid_ReadGroup.BEST_SM,
        verifybamid_ReadGroup.BESTIBD,
        verifybamid_ReadGroup.BESTMIX,
        verifybamid_ReadGroup.status
    ";
}

sub table_name {
  return "verifybamid_ReadGroup";
}

sub store {
  my ( $self, $verifybamid_ReadGroup ) = @_;

  throw(  "Can't store "
	  . $verifybamid_ReadGroup
	  . " using ReseqTrack::DBSQL::VerifyBamIDReadGroupAdaptor" )
    unless (
	    $verifybamid_ReadGroup->isa(
					"ReseqTrack::VerifyBamIDReadGroup")
	   );

  my $sql =
    "insert into verifybamid_ReadGroup (other_id, "
      . "run_id, SELFIBD, SELFMIX ,BEST_SM,BESTIBD,BESTMIX,status) "
	. " values( ?,?,?,?,?,?,?,?) ";

  my $sth = $self->prepare($sql);

  $sth->bind_param( 1, $verifybamid_ReadGroup->other_id );
  $sth->bind_param( 2, $verifybamid_ReadGroup->run_id );
  $sth->bind_param( 3, $verifybamid_ReadGroup->SELFIBD );
  $sth->bind_param( 4, $verifybamid_ReadGroup->SELFMIX );
  $sth->bind_param( 5, $verifybamid_ReadGroup->BEST_SM );
  $sth->bind_param( 6, $verifybamid_ReadGroup->BESTIBD );
  $sth->bind_param( 7, $verifybamid_ReadGroup->BESTMIX );
  $sth->bind_param( 8, $verifybamid_ReadGroup->status );

  my $rows_inserted = $sth->execute();
  my $dbID          = $sth->{'mysql_insertid'};

  $verifybamid_ReadGroup->dbID($dbID);
  $verifybamid_ReadGroup->adaptor($self);
  $sth->finish();

  return $verifybamid_ReadGroup;
}

sub update {
  my ( $self, $verifybamid_ReadGroup ) = @_;
  
  throw(  "Can't store "
	  . $verifybamid_ReadGroup
	  . " using ReseqTrack::DBSQL::VerifyBamIDReadGroupAdaptor" )
    unless (
	    $verifybamid_ReadGroup->isa(
					"ReseqTrack::VerifyBamIDReadGroup"
				       )
	   );

  my $sql =
    "update verifybamid_ReadGroup set "
      . " SELFIBD = ? , SELFMIX = ? , BEST_SM = ?, "
	. "BESTIBD  = ? ,BESTMIX = ?, status = ? where  other_id = ?  and run_id = ?";

 
  
  my $sth = $self->prepare($sql);

  $sth->bind_param( 1, $verifybamid_ReadGroup->SELFIBD );
  $sth->bind_param( 2, $verifybamid_ReadGroup->SELFMIX );
  $sth->bind_param( 3, $verifybamid_ReadGroup->BEST_SM );
  $sth->bind_param( 4, $verifybamid_ReadGroup->BESTIBD );
  $sth->bind_param( 5, $verifybamid_ReadGroup->BESTMIX );
  $sth->bind_param( 6, $verifybamid_ReadGroup->status );
  $sth->bind_param( 7, $verifybamid_ReadGroup->other_id );
  $sth->bind_param( 8, $verifybamid_ReadGroup->run_id );

  $sth->execute();
  $sth->finish();

  return $verifybamid_ReadGroup;
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
		      -dbID                  => $hashref->{verifybamid_ReadGroup_id},
		      -other_id              => $hashref->{other_id},
		      -run_id                => $hashref->{run_id},
		      -SELFIBD               => $hashref->{SELFIBD},
		      -SELFMIX               => $hashref->{SELFMIX},
		      -BEST_SM               => $hashref->{BEST_SM},
		      -BESTIBD               => $hashref->{BESTIBD},
		      -BESTMIX               => $hashref->{BESTMIX},
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
 
  my $result = $results[0];
    
  return \$result;
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

  my ( $self, $verifybamid_ReadGroup ) = @_;
  
  throw(  "Can't store "
	  . $verifybamid_ReadGroup
	  . " using ReseqTrack::DBSQL::VerifyBamIDReadGroupAdaptor" )
    unless (
	    $verifybamid_ReadGroup->isa(
					"ReseqTrack::VerifyBamIDReadGroup"
				       )
	   );

  my $sql = "update verifybamid_ReadGroup set  run_id = ? ".
    " where verifybamid_readgroup_id  = ?";
	        
  my $sth = $self->prepare($sql);

  $sth->bind_param( 1, $verifybamid_ReadGroup->run_id );
  $sth->bind_param( 2, $verifybamid_ReadGroup->dbID );

  $sth->execute();
  $sth->finish();

  return $verifybamid_ReadGroup;
}

1;

