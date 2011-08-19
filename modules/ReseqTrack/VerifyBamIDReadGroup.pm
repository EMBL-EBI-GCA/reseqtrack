package ReseqTrack::VerifyBamIDReadGroup;
 
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
		$run_id,
		$SELFIBD,
		$SELFMIX,
		$BEST_SM,
		$BESTIBD,
		$BESTMIX,
	        $status

	  ) = rearrange(
		[
			qw(
			  OTHER_ID
			  RUN_ID
			  SELFIBD
			  SELFMIX
			  BEST_SM
			  BESTIBD
			  BESTMIX
                          STATUS
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
	$self->status($status);	


	return $self;
}
sub run_id {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{run_id} = $arg;
	}
	return $self->{run_id};
}



sub status {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{status} = $arg;
    }
    return $self->{status};
}


sub other_id {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{other_id} = $arg;
    }
    return $self->{other_id};
}


#4. SELFIBD/BESTIBD
sub SELFIBD {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{SELFIBD} = $arg;
	}
	return $self->{SELFIBD};
}



sub swap_candidate {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{swap_candidate} = $arg;
	}
	return $self->{swap_candidate};
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

