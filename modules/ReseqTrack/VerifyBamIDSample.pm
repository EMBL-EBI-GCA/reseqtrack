package ReseqTrack::VerifyBamIDSample;

use strict;
use warnings;
use vars qw(@ISA);
use ReseqTrack::DBSQL::BaseAdaptor;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename;
use Data::Dumper;

  
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
        $status,
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
              STATUS
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
    $self->status($status);
    $self->performed($performed);
    return $self;
}


sub status {
  my ( $self, $arg ) = @_;
  if ( defined $arg ) {
	  $self->{status} = $arg;
	}
  return $self->{status};
}

sub performed {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{performed} = $arg;
	}
	return $self->{performed};
}

#5. SELFIBDLLK/BESTIBDLLK : Log likelihood of the sequence reads
# given the MLE SELFIBD/BESTIBD with SELF_SM/BEST_SM
sub SELFIBDLLK {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{SELFIBDLLK} = $arg;
	}
	return $self->{SELFIBDLLK};
}

#6. SELFIBDLK-/BESTIBDLK- : Difference of log-likelihood between
# SELFIBDLLK/BESTIBDLLK and the likelihood of reads under no contamination (SELFIBD=1)
sub SELFIBDLLKdiff {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{SELFIBDLLKdiff} = $arg;
	}
	return $self->{SELFIBDLLKdiff};
}


#4. SELFIBD/BESTIBD
sub SELFIBD {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{SELFIBD} = $arg;
	}
	return $self->{SELFIBD};
}

#SEQ_SM : Sample ID of the sequenced sample. Obtained from @RG header / SM tag in the BAM file
sub SEQ_SM {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{SEQ_SM} = $arg;
	}
	return $self->{SEQ_SM};
}

#25. %MIX : Maximum-likelihood estimate of % of contamination based on
# two-sample mixture model from population allele frequency.
sub MIX {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{MIX} = $arg;
	}
	return $self->{MIX};
}



sub HOM {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{HOM} = $arg;
	}
	return $self->{HOM};
}


#HET-A1% : Fraction of reference bases in Heterozygous site
sub HET_A1 {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{HET_A1} = $arg;
	}
	return $self->{HET_A1};
}


#21.#DP>1 : Number of sites with depth of 2 or greater
sub DP {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{DP} = $arg;
	}
	return $self->{DP};
}

sub BESTHOMMIXLLKdiff {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{BESTHOMMIXLLKdiff} = $arg;
	}
	return $self->{BESTHOMMIXLLKdiff};
}


sub BESTHOMMIXLLK {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{BESTHOMMIXLLK} = $arg;
	}
	return $self->{BESTHOMMIXLLK};
}

#ALT-A1% : Fraction of reference bases in HomAlt site
sub ALT_A1 {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{ALT_A1} = $arg;
	}
	return $self->{ALT_A1};
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
    

sub other_id {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{other_id} = $arg;
    }
    return $self->{other_id};
}

sub table_name {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{table_name} = $arg;
    }
    return $self->{table_name};
}



#INFO extracted from results files

#cat HG00308.selfRG  | cut -f 1,4,5,6,17,19,21,25,26,27,28
#SEQ_RG  SELFIBD SELFIBDLLK      SELFIBDLLK-     HET-A1% ALT-A1% #DP     %MIX    %HOM    BESTHOMMIXLLK   BESTHOMMIXLLK-
#HG00308 0.680   7.305e+04       2.892e+03       0.64466 0.32623 0.294   0.180   0.000   -4.718e+03      2.480e+01


1;
