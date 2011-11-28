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
		$selfibd,
		$selfmix,
		$best_sample,
		$bestibd,
		$bestmix,
	    $status

	  ) = rearrange(
		[
			qw(
			  OTHER_ID
			  RUN_ID
			  SELFIBD
			  SELFMIX
			  BEST_SAMPLE
			  BESTIBD
			  BESTMIX
              STATUS
			  )
		],
		@args
	  );
	$self->other_id($other_id);
	$self->run_id($run_id);
    
	$self->selfibd($selfibd);
	$self->selfmix($selfmix);
	$self->best_sample($best_sample);
	$self->bestibd($bestibd);
	$self->bestmix($bestmix);
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


#4. SELFIBD/bestibd
sub selfibd {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{selfibd} = $arg;
	}
	return $self->{selfibd};
}



sub swap_candidate {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{swap_candidate} = $arg;
	}
	return $self->{swap_candidate};
}




sub selfmix {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{selfmix} = $arg;
    }
    return $self->{selfmix};
}

sub best_sample {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{best_sample} = $arg;
    }
    return $self->{best_ssample};
}

sub bestibd {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{bestibd} = $arg;
    }
    return $self->{bestibd};
}


sub bestmix {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{bestmix} = $arg;
    }
    return $self->{bestmix};
}


1;

