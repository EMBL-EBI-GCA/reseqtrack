package ReseqTrack::VerifyBamIDReadGroup;
 
use strict;
use warnings;
use vars qw(@ISA);
use ReseqTrack::DBSQL::BaseAdaptor;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename;
use Data::Dumper;
use ReseqTrack::Tools::RunVerifyBamIDUtils qw (
  swap_candidate run_id SELFIBD other_id );

@ISA = qw(ReseqTrack::Base);

sub new {

	my ( $class, @args ) = @_;
	my $self = $class->SUPER::new(@args);

	my (
		$other_id,
		$run_id,
		$SELFIBD,
		$SELFMIX,
		$BEST_SM,
		$BESTIBD,
		$BESTMIX,

	  ) = rearrange(
		[
			qw(
			  other_id
			  RUN_ID
			  SELFIBD
			  SELFMIX
			  BEST_SM
			  BESTIBD
			  BESTMIX
			  )
		],
		@args
	  );
	$self->other_id($other_id);
	$self->run_id($run_id);
    
	$self->SELFIBD($SELFIBD);
	$self->SELFMIX($SELFMIX);
	$self->BEST_SM($BEST_SM);
	$self->BESTIBD($BESTIBD);
	$self->BESTMIX($BESTMIX);
	#print Dumper $self;

	return $self;
}


sub SELFMIX {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{SELFMIX} = $arg;
    }
    return $self->{SELFMIX};
}

sub BEST_SM {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{BEST_SM} = $arg;
    }
    return $self->{BEST_SM};
}

sub BESTIBD {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{BESTIBD} = $arg;
    }
    return $self->{BESTIBD};
}


sub BESTMIX {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{BESTMIX} = $arg;
    }
    return $self->{BESTMIX};
}

1;

