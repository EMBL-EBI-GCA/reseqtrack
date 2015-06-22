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
		$sample_name,
		$selfibd,
		$selfibdllk,
		$selfibdllkdiff,
		$het_a1,
		$alt_a1,
		$dp,
		$mix,
		$hom,
		$besthommixllk,
		$besthommixllkdiff,
		$num_run_ids,
		$num_low_selfibd_run_ids,
		$sequence_index,
		$analysis_group,
		$chr20,
		$failed,
		$status,
		$performed,

	  ) = rearrange(
		[
			qw(
			  OTHER_ID
			  TABLE_NAME
			  SAMPLE_NAME
			  SELFIBD
			  SELFIBDLLK
			  SELFIBDLLKDIFF
			  HET_A1
			  ALT_A1
			  DP
			  MIX
			  HOM
			  BESTHOMMIXLLK
			  BESTHOMMIXLLKDIFF
			  NUM_RUN_IDS
			  NUM_LOW_SELFIBD_RUN_IDS
			  SEQUENCE_INDEX
			  ANALYSIS_GROUP
			  CHR20
			  FAILED
			  STATUS
			  PERFORMED
			  )
		],
		@args
	  );

	$self->other_id($other_id);
	$self->table_name($table_name);
	$self->sample_name($sample_name);
	$self->selfibd($selfibd);
	$self->selfibdllk($selfibdllk);
	$self->selfibdllkdiff($selfibdllkdiff);
	$self->het_a1($het_a1);
	$self->alt_a1($alt_a1);
	$self->dp($dp);
	$self->mix($mix);
	$self->hom($hom);
	$self->besthommixllk($besthommixllk);
	$self->besthommixllkdiff($besthommixllkdiff);
	$self->num_run_ids($num_run_ids);
	$self->num_low_selfibd_run_ids($num_low_selfibd_run_ids);
	$self->sequence_index($sequence_index);
	$self->analysis_group($analysis_group);
	$self->chr20($chr20);
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

#5. selfibdllk/BESTIBDllk : Log likelihood of the sequence reads
# given the MLE selfibd/BESTIBD with SELF_SM/BEST_SM
sub selfibdllk {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{selfibdllk} = $arg;
	}
	return $self->{selfibdllk};
}

#6. selfibdLK-/BESTIBDLK- : Difference of log-likelihood between
# selfibdllk/BESTIBDllk and the likelihood of reads under no contamination (selfibd=1)
sub selfibdllkdiff {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{selfibdllkdiff} = $arg;
	}
	return $self->{selfibdllkdiff};
}

#4. selfibd/BESTIBD
sub selfibd {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{selfibd} = $arg;
	}
	return $self->{selfibd};
}

#sample_name : Sample ID of the sequenced sample. Obtained from @RG header / SM tag in the BAM file
sub sample_name {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{sample_name} = $arg;
	}
	return $self->{sample_name};
}

#25. %mix : Maximum-likelihood estimate of % of contamination based on
# two-sample mixture model from population allele frequency.
sub mix {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{mix} = $arg;
	}
	return $self->{mix};
}

sub hom {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{hom} = $arg;
	}
	return $self->{hom};
}

#HET-A1% : Fraction of reference bases in Heterozygous site
sub het_a1 {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{het_a1} = $arg;
	}
	return $self->{het_a1};
}

#21.#dp>1 : Number of sites with depth of 2 or greater
sub dp {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{dp} = $arg;
	}
	return $self->{dp};
}

sub besthommixllkdiff {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{besthommixllkdiff} = $arg;
	}
	return $self->{besthommixllkdiff};
}

sub besthommixllk {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{besthommixllk} = $arg;
	}
	return $self->{besthommixllk};
}

#ALT-A1% : Fraction of reference bases in HomAlt site
sub alt_a1 {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{alt_a1} = $arg;
	}
	return $self->{alt_a1};
}

sub num_run_ids {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{num_run_ids} = $arg;
	}
	return $self->{num_run_ids};
}

sub num_low_selfibd_run_ids {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{num_low_selfibd_run_ids} = $arg;
	}
	return $self->{num_low_selfibd_run_ids};
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

sub sequence_index {

	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{sequence_index} = $arg;
	}
	return $self->{sequence_index};

}

sub analysis_group {

	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{analysis_group} = $arg;
	}
	return $self->{analysis_group};

}

sub chr20 {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{chr20} = $arg;
	}
	return $self->{chr20};

}

#INFO extracted from results files

#cat HG00308.selfRG  | cut -f 1,4,5,6,17,19,21,25,26,27,28
#SEQ_RG  SELFIBD SELFIBDLLK      SELFIBDLLK-     HET-A1% ALT-A1% #DP     %MIX    %HOM    besthomMIXLLK   besthomMIXLLK-
#HG00308 0.680   7.305e+04       2.892e+03       0.64466 0.32623 0.294   0.180   0.000   -4.718e+03      2.480e+01

1;
