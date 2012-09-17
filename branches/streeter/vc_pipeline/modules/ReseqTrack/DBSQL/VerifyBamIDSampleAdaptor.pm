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
        verifybamid_sample.verifybamid_sample_id,
        verifybamid_sample.other_id,
        verifybamid_sample.table_name,
        verifybamid_sample.sample_name,
        verifybamid_sample.selfibd, 
        verifybamid_sample.selfibdllk,
        verifybamid_sample.selfibdllkdiff,
        verifybamid_sample.het_a1,
        verifybamid_sample.alt_a1,
        verifybamid_sample.dp,
        verifybamid_sample.mix,
        verifybamid_sample.hom  ,
        verifybamid_sample.besthommixllk,
        verifybamid_sample.besthommixllkdiff,
        verifybamid_sample.num_run_ids,
        verifybamid_sample.num_low_selfibd_run_ids,
        verifybamid_sample.sequence_index,
        verifybamid_sample.analysis_group,
        verifybamid_sample.chr20,  
        verifybamid_sample.failed,
        verifybamid_sample.status,
        verifybamid_sample.performed    "
}

sub table_name {
  return "verifybamid_sample";
}

sub store {
  my ( $self, $verifybamid_sample ) = @_;

  throw(  "Can't store "
	  . $verifybamid_sample
	  . " using ReseqTrack::DBSQL::VerifyBamIDSampleAdaptor" )
    unless (
	    $verifybamid_sample->isa(
				     "ReseqTrack::VerifyBamIDSample")
	   );

  my $sql =
    "insert into verifybamid_sample (other_id, table_name,"
      . "sample_name, selfibd, selfibdllk, selfibdllkdiff,het_a1, alt_a1,"
	. "dp, mix,hom  , besthommixllk, besthommixllkdiff,"
	  . "num_run_ids, num_low_selfibd_run_ids, sequence_index,analysis_group, " 
	    . " chr20, failed, status,performed ) "
	      . " values( ?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,now() ) ";

  my $sth = $self->prepare($sql);
  $sth->bind_param( 1,  $verifybamid_sample->other_id );
  $sth->bind_param( 2,  $verifybamid_sample->table_name );
  $sth->bind_param( 3,  $verifybamid_sample->sample_name );
  $sth->bind_param( 4,  $verifybamid_sample->selfibd );
  $sth->bind_param( 5,  $verifybamid_sample->selfibdllk );
  $sth->bind_param( 6,  $verifybamid_sample->selfibdllkdiff );
  $sth->bind_param( 7,  $verifybamid_sample->het_a1 );
  $sth->bind_param( 8,  $verifybamid_sample->alt_a1 );
  $sth->bind_param( 9,  $verifybamid_sample->dp );
  $sth->bind_param( 10, $verifybamid_sample->mix );
  $sth->bind_param( 11, $verifybamid_sample->hom );
  $sth->bind_param( 12, $verifybamid_sample->besthommixllk );
  $sth->bind_param( 13, $verifybamid_sample->besthommixllkdiff );
  $sth->bind_param( 14, $verifybamid_sample->num_run_ids );
  $sth->bind_param( 15, $verifybamid_sample->num_low_selfibd_run_ids );

  $sth->bind_param( 16, $verifybamid_sample->sequence_index );
  $sth->bind_param( 17, $verifybamid_sample->analysis_group );
  $sth->bind_param( 18, $verifybamid_sample->chr20 );
  $sth->bind_param( 19, $verifybamid_sample->failed );
  $sth->bind_param( 20, $verifybamid_sample->status );

  my $rows_inserted = $sth->execute();
  my $dbID          = $sth->{'mysql_insertid'};

  $verifybamid_sample->dbID($dbID);
  $verifybamid_sample->adaptor($self);
  $sth->finish();

  return $verifybamid_sample;
}

sub update {

  my ( $self, $verifybamid_sample ) = @_;
	
  throw(  "Can't update "
	  . $verifybamid_sample
	  . " using ReseqTrack::DBSQL::VerifyBamIDSampleAdaptor" )
    unless (
	    $verifybamid_sample->isa(
				     "ReseqTrack::VerifyBamIDSample")
	   );

  my $sql =
    " update verifybamid_sample  set sample_name = ?, selfibd = ? ,"
      . " selfibdllk = ? , selfibdllkdiff = ? ,het_a1 = ? , alt_a1 = ?,"
	. " dp = ? , mix = ? ,hom = ?   , besthommixllk = ? , "
	  . " besthommixllkdiff = ?, num_run_ids = ? , "
	    . " num_low_selfibd_run_ids = ?, sequence_index = ? , analysis_group = ? , "
	      . " chr20 = ? ,failed = ? , status = ? " 
	      . ", performed = now()  where other_id  = ? and table_name = ?";

  my $sth = $self->prepare($sql);

  $sth->bind_param( 1,  $verifybamid_sample->sample_name );
  $sth->bind_param( 2,  $verifybamid_sample->selfibd );
  $sth->bind_param( 3,  $verifybamid_sample->selfibdllk );
  $sth->bind_param( 4,  $verifybamid_sample->selfibdllkdiff );
  $sth->bind_param( 5,  $verifybamid_sample->het_a1 );
  $sth->bind_param( 6,  $verifybamid_sample->alt_a1 );
  $sth->bind_param( 7,  $verifybamid_sample->dp );
  $sth->bind_param( 8,  $verifybamid_sample->mix );
  $sth->bind_param( 9,  $verifybamid_sample->hom );
  $sth->bind_param( 10, $verifybamid_sample->besthommixllk );
  $sth->bind_param( 11, $verifybamid_sample->besthommixllkdiff );
  $sth->bind_param( 12, $verifybamid_sample->num_run_ids );
  $sth->bind_param( 13, $verifybamid_sample->num_low_selfibd_run_ids );
  $sth->bind_param( 14, $verifybamid_sample->sequence_index );
  $sth->bind_param( 15, $verifybamid_sample->analysis_group );
  $sth->bind_param( 16, $verifybamid_sample->chr20 );
  $sth->bind_param( 17, $verifybamid_sample->failed );
  $sth->bind_param( 18, $verifybamid_sample->status );
  $sth->bind_param( 19, $verifybamid_sample->other_id );
  $sth->bind_param( 20, $verifybamid_sample->table_name );

  $sth->execute();
  $sth->finish();

  return $verifybamid_sample;
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
		      -dbID                  => $hashref->{verifybamid_sample_id},
		      -other_id              => $hashref->{other_id},
		      -table_name            => $hashref->{table_name},
		      -sample_name           => $hashref->{sample_name},
		      -selfibd               => $hashref->{selfibd},
		      -selfibdllk            => $hashref->{selfibdllk},
		      -selfibdllkdiff        => $hashref->{selfibdllkdiff},
		      -het_a1                => $hashref->{het_a1},
		      -alt_a1                => $hashref->{alt_a1},
		      -dp                    => $hashref->{dp},
		      -mix                   => $hashref->{mix},
		      -hom                   => $hashref->{hom},
		      -besthommixllk         => $hashref->{besthommixllk},
		      -besthommixllkdiff     => $hashref->{besthommixllkdiff},
		      -num_run_ids           => $hashref->{num_run_ids},
		      -num_low_selfibd_run_ids =>$hashref->{num_low_selfibd_run_ids},
		      -sequence_index        =>$hashref->{sequence_index},
		      -analysis_group        =>$hashref->{analysis_group},
              -chr20                  =>$hashref->{chr20},
		      -failed    => $hashref->{failed},
		      -status    => $hashref->{status},
		      -performed => $hashref->{performed},

		     );
  return $obj;
}


sub fetch_by_other_id{
	
  my ($self, $other_id) = @_;
  my $sql = "select ".$self->columns." from  verifybamid_sample ".
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
