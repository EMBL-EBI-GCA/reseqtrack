package ReseqTrack::DBSQL::VerifyBamIDSampleAdaptor;

use strict;
use warnings;
 
use ReseqTrack::DBSQL::BaseAdaptor;
use ReseqTrack::Base;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename;
use Data::Dumper;
use ReseqTrack::VerifyBamIDSample;
use vars qw(@ISA);
@ISA = qw(ReseqTrack::DBSQL::BaseAdaptor);

sub new {
  my ( $class, $db ) = @_;

  my $self = $class->SUPER::new($db);

  return $self;
}

#cat HG00308.selfSM  | cut -f 1,4,5,6,17,19,21,25,26,27,28
#SEQ_SM  SELFIBD SELFIBDLLK      SELFIBDLLK-     HET-A1% ALT-A1% #DP     %MIX    %HOM    BESTHOMMIXLLK   BESTHOMMIXLLK-
#HG00308 0.680   7.305e+04       2.892e+03       0.64466 0.32623 0.294   0.180   0.000   -4.718e+03      2.480e+01

sub columns {
  return "
        verifybamid_Sample.verifybamid_Sample_id,
        verifybamid_Sample.other_id,
        verifybamid_Sample.table_name,
        verifybamid_Sample.SEQ_SM,
        verifybamid_Sample.SELFIBD, 
        verifybamid_Sample.SELFIBDLLK,
        verifybamid_Sample.SELFIBDLLKdiff,
        verifybamid_Sample.HET_A1,
        verifybamid_Sample.ALT_A1,
        verifybamid_Sample.DP,
        verifybamid_Sample.MIX,
        verifybamid_Sample.HOM  ,
        verifybamid_Sample.BESTHOMMIXLLK,
        verifybamid_Sample.BESTHOMMIXLLKdiff,
        verifybamid_Sample.num_run_ids,
        verifybamid_Sample.num_low_selfIBD_run_ids,
        verifybamid_Sample.failed,
        verifybamid_Sample.status,
        verifybamid_Sample.performed    "
}

sub table_name {
  return "verifybamid_Sample";
}

sub store {
  my ( $self, $verifybamid_Sample ) = @_;


  # print Dumper $verifybamid_Sample;

  throw(  "Can't store "
	  . $verifybamid_Sample
	  . " using ReseqTrack::DBSQL::VerifyBamIDSampleAdaptor" )
    unless (
	    $verifybamid_Sample->isa(
				     "ReseqTrack::VerifyBamIDSample")
	   );

  my $sql =
    "insert into verifybamid_Sample (other_id, table_name,"
      . "SEQ_SM, SELFIBD, SELFIBDLLK, SELFIBDLLKdiff,HET_A1, ALT_A1,"
	. "DP, MIX,HOM  , BESTHOMMIXLLK, BESTHOMMIXLLKdiff,"
	  . "num_run_ids, num_low_selfIBD_run_ids, failed, status,performed ) "
	    . " values( ?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,now() ) ";

  my $sth = $self->prepare($sql);
  $sth->bind_param( 1,  $verifybamid_Sample->other_id );
  $sth->bind_param( 2,  $verifybamid_Sample->table_name );
  $sth->bind_param( 3,  $verifybamid_Sample->SEQ_SM );
  $sth->bind_param( 4,  $verifybamid_Sample->SELFIBD );
  $sth->bind_param( 5,  $verifybamid_Sample->SELFIBDLLK );
  $sth->bind_param( 6,  $verifybamid_Sample->SELFIBDLLKdiff );
  $sth->bind_param( 7,  $verifybamid_Sample->HET_A1 );
  $sth->bind_param( 8,  $verifybamid_Sample->ALT_A1 );
  $sth->bind_param( 9,  $verifybamid_Sample->DP );
  $sth->bind_param( 10, $verifybamid_Sample->MIX );
  $sth->bind_param( 11, $verifybamid_Sample->HOM );
  $sth->bind_param( 12, $verifybamid_Sample->BESTHOMMIXLLK );
  $sth->bind_param( 13, $verifybamid_Sample->BESTHOMMIXLLKdiff );
  $sth->bind_param( 14, $verifybamid_Sample->num_run_ids );
  $sth->bind_param( 15, $verifybamid_Sample->num_low_selfIBD_run_ids );
  $sth->bind_param( 16, $verifybamid_Sample->failed );
  $sth->bind_param( 17, $verifybamid_Sample->status );

  my $rows_inserted = $sth->execute();
  my $dbID          = $sth->{'mysql_insertid'};

  $verifybamid_Sample->dbID($dbID);
  $verifybamid_Sample->adaptor($self);
  $sth->finish();

  return $verifybamid_Sample;
}

sub update {

  my ( $self, $verifybamid_Sample ) = @_;
	
  throw(  "Can't update "
	  . $verifybamid_Sample
	  . " using ReseqTrack::DBSQL::VerifyBamIDSampleAdaptor" )
    unless (
	    $verifybamid_Sample->isa(
				     "ReseqTrack::VerifyBamIDSample")
	   );

  my $sql =
    " update verifybamid_Sample  set SEQ_SM = ?, SELFIBD = ? ,"
      . " SELFIBDLLK = ? , SELFIBDLLKdiff = ? ,HET_A1 = ? , ALT_A1 = ?,"
	. " DP = ? , MIX = ? ,HOM = ?   , BESTHOMMIXLLK = ? , "
	  . " BESTHOMMIXLLKdiff = ?, num_run_ids = ? , "
	    . " num_low_selfIBD_run_ids = ?, failed = ? , status = ? " 
	      . ", performed = now()  where other_id  = ? and table_name = ?";

  my $sth = $self->prepare($sql);

  $sth->bind_param( 1,  $verifybamid_Sample->SEQ_SM );
  $sth->bind_param( 2,  $verifybamid_Sample->SELFIBD );
  $sth->bind_param( 3,  $verifybamid_Sample->SELFIBDLLK );
  $sth->bind_param( 4,  $verifybamid_Sample->SELFIBDLLKdiff );
  $sth->bind_param( 5,  $verifybamid_Sample->HET_A1 );
  $sth->bind_param( 6,  $verifybamid_Sample->ALT_A1 );
  $sth->bind_param( 7,  $verifybamid_Sample->DP );
  $sth->bind_param( 8,  $verifybamid_Sample->MIX );
  $sth->bind_param( 9,  $verifybamid_Sample->HOM );
  $sth->bind_param( 10, $verifybamid_Sample->BESTHOMMIXLLK );
  $sth->bind_param( 11, $verifybamid_Sample->BESTHOMMIXLLKdiff );
  $sth->bind_param( 12, $verifybamid_Sample->num_run_ids );
  $sth->bind_param( 13, $verifybamid_Sample->num_low_selfIBD_run_ids );
  $sth->bind_param( 14, $verifybamid_Sample->failed );
  $sth->bind_param( 15, $verifybamid_Sample->status );
  $sth->bind_param( 16, $verifybamid_Sample->other_id );
  $sth->bind_param( 17, $verifybamid_Sample->table_name );

  $sth->execute();
  $sth->finish();

  return $verifybamid_Sample;
}

sub object_from_hashref {
  my ( $self, $hashref ) = @_;

  throw(
	"Can't create a ReseqTrack::VerifyBamIDSample from an undefined hashref
	  "
       ) if ( !$hashref );

  my $OBJ="ReseqTrack::VerifyBamIDSample";
  my $obj = $OBJ->new(

		      -adaptor               => $self,
		      -dbID => $hashref->{verifybamid_Sample_id},
		      -other_id              => $hashref->{other_id},
		      -table_name            => $hashref->{table_name},
		      -SEQ_SM                => $hashref->{SEQ_SM},
		      -SELFIBD               => $hashref->{SELFIBD},
		      -SELFIBDLLK            => $hashref->{SELFIBDLLK},
		      -SELFIBDLLKdiff        => $hashref->{SELFIBDLLKdiff},
		      -HET_A1                => $hashref->{HET_A1},
		      -ALT_A1                => $hashref->{ALT_A1},
		      -DP                    => $hashref->{DP},
		      -MIX                   => $hashref->{MIX},
		      -HOM                   => $hashref->{HOM},
		      -BESTHOMMIXLLK         => $hashref->{BESTHOMMIXLLK},
		      -BESTHOMMIXLLKdiff     => $hashref->{BESTHOMMIXLLKdiff},
		      -num_run_ids           => $hashref->{num_run_ids},
		      -num_low_selfIBD_run_ids =>$hashref->{num_low_selfIBD_run_ids},
		      -failed    => $hashref->{failed},
		      -status    => $hashref->{status},
		      -performed => $hashref->{performed},

		     );
  return $obj;
}


sub fetch_by_other_id{
	
  my ($self, $other_id) = @_;
  my $sql = "select ".$self->columns." from  verifybamid_Sample ".
    "where other_id = ? ";
      
     
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $other_id);
  $sth->execute;
  
  my @results;

  while (my $rowHashref = $sth->fetchrow_hashref) {
    my $result = $self->object_from_hashref($rowHashref);
    push(@results, $result);
  }
  $sth->finish;

  throw($other_id . 
	" has returned multiple objects ".@results." not sure what to do") 
    if (@results && @results >= 2);

  my $result = $results[0];
  
#  print "Got Sample result for file_id = $other_id\n";
  
  return $result;
}



1;
