package ReseqTrack::Tools::RunVerifyBamIDUtils;
 
use strict;
use warnings;

use Exporter;
use File::Basename;

use vars qw (@ISA  @EXPORT);

@ISA    = qw(Exporter);
@EXPORT = qw(
  other_id
  table_name
  ALT_A1
  BESTHOMMIXLLK
  BESTHOMMIXLLKdiff
  DP
  HET_A1
  HOM MIX
  swap_candidate
  SEQ_SM
  SEQ_RG
  SELFIBD
  SELFIBDLLK
  SELFIBDLLKdiff
  run_id
  SEQ_SM
  performed
);

#ALT-A1% : Fraction of reference bases in HomAlt site
sub ALT_A1 {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{ALT_A1} = $arg;
	}
	return $self->{ALT_A1};
}

sub BESTHOMMIXLLK {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{BESTHOMMIXLLK} = $arg;
	}
	return $self->{BESTHOMMIXLLK};
}

sub BESTHOMMIXLLKdiff {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{BESTHOMMIXLLKdiff} = $arg;
	}
	return $self->{BESTHOMMIXLLKdiff};
}

#21.#DP>1 : Number of sites with depth of 2 or greater
sub DP {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{DP} = $arg;
	}
	return $self->{DP};
}

#HET-A1% : Fraction of reference bases in Heterozygous site
sub HET_A1 {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{HET_A1} = $arg;
	}
	return $self->{HET_A1};
}

sub HOM {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{HOM} = $arg;
	}
	return $self->{HOM};
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

sub performed {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{performed} = $arg;
	}
	return $self->{performed};
}

sub swap_candidate {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{swap_candidate} = $arg;
	}
	return $self->{swap_candidate};
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

sub run_id {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{run_id} = $arg;
	}
	return $self->{run_id};
}

#SEQ_SM : Sample ID of the sequenced sample. Obtained from @RG header / SM tag in the BAM file
sub SEQ_SM {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{SEQ_SM} = $arg;
	}
	return $self->{SEQ_SM};
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

1;
