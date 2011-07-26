package ReseqTrack::VerifyBamIDSample;

use strict;
use warnings;
use vars qw(@ISA);
use ReseqTrack::DBSQL::BaseAdaptor;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename;
use Data::Dumper;
use ReseqTrack::Tools::RunVerifyBamIDUtils qw (ALT_A1
  BESTHOMMIXLLK   BESTHOMMIXLLKdiff DP HET_A1
  HOM MIX  SEQ_SM run_id SELFIBD
  SELFIBDLLK SELFIBDLLKdiff SEQ_RG SEQ_SM performed other_id table_name);
  
@ISA = qw(ReseqTrack::Base);

sub new {

    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);

    my (
        $other_id,
        $table_name,
        $SEQ_SM,
        $SELFIBD,
        $SELFIBDLLK,
        $SELFIBDLLKdiff,
        $HET_A1,
        $ALT_A1,
        $DP,
        $MIX,
        $HOM,
        $BESTHOMMIXLLK,
        $BESTHOMMIXLLKdiff,
        $num_run_ids,
        $num_low_selfIBD_run_ids,
        $failed,
        $performed,

      ) = rearrange(
        [
            qw(
              other_id
              table_name
              SEQ_SM
              SELFIBD
              SELFIBDLLK
              SELFIBDLLKdiff
              HET_A1
              ALT_A1
              DP
              MIX
              HOM  
              BESTHOMMIXLLK
              BESTHOMMIXLLKdiff
              num_run_ids
              num_low_selfIBD_run_ids
              FAILED
              performed
              )
        ],
        @args
      );
    $self->other_id ($other_id);
    $self->table_name($table_name);
    $self->SEQ_SM($SEQ_SM);
    $self->SELFIBD($SELFIBD);
    $self->SELFIBDLLK($SELFIBDLLK);
    $self->SELFIBDLLKdiff($SELFIBDLLKdiff);
    $self->HET_A1($HET_A1);
    $self->ALT_A1($ALT_A1);
    $self->DP($DP);
    $self->MIX($MIX);
    $self->HOM($HOM);
    $self->BESTHOMMIXLLK($BESTHOMMIXLLK);
    $self->BESTHOMMIXLLKdiff($BESTHOMMIXLLKdiff);
    $self->num_run_ids($num_run_ids);
    $self->num_low_selfIBD_run_ids($num_low_selfIBD_run_ids);
    $self->failed($failed);
    $self->performed($performed);
    return $self;
}

sub num_run_ids {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{num_run_ids} = $arg;
    }
    return $self->{num_run_ids};
}
sub  num_low_selfIBD_run_ids{
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{num_low_selfIBD_run_ids} = $arg;
    }
    return $self->{num_low_selfIBD_run_ids};
}

sub failed {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{failed} = $arg;
    }
    return $self->{failed};
}
    

#cat HG00308.selfRG  | cut -f 1,4,5,6,17,19,21,25,26,27,28
#SEQ_RG  SELFIBD SELFIBDLLK      SELFIBDLLK-     HET-A1% ALT-A1% #DP     %MIX    %HOM    BESTHOMMIXLLK   BESTHOMMIXLLK-
#HG00308 0.680   7.305e+04       2.892e+03       0.64466 0.32623 0.294   0.180   0.000   -4.718e+03      2.480e+01


1;
