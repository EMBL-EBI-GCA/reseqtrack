package ReseqTrack::DBSQL::GenotypeResultsAdaptor;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::DBSQL::BaseAdaptor;
use ReseqTrack::Base;
use ReseqTrack::GenotypeResults;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::FileUtils qw(are_files_identical);
use File::Basename;

@ISA = qw(ReseqTrack::DBSQL::BaseAdaptor);

sub new {
	my ( $class, $db ) = @_;

	my $self = $class->SUPER::new($db);

	return $self;
}


sub columns {
	return " genotype_results.genotype_results_id,
        genotype_results.other_id, 
        genotype_results.table_name,
        genotype_results.reference  ,
        genotype_results.snps_bin  ,
        genotype_results.aligner,
        genotype_results.version,
        genotype_results.validation_method,
        genotype_results.max_bases,
        genotype_results.percent_mapped  , 
        genotype_results.summary,
        genotype_results.verdict    ,
        genotype_results.performed";
}

sub table_name {
	return "genotype_results";
}

sub store {
	my ( $self, $genotype_results ) = @_;
	throw(  "Can't store "
		  . $genotype_results
		  . " using ReseqTrack::DBSQL::GenotypeResultsAdaptor" )
	  unless ( $genotype_results->isa("ReseqTrack::GenotypeResults") );

	my $sql = "insert ignore into genotype_results (
                    other_id,
                    table_name,
                    reference,
                    snps_bin,
                    aligner,
                    version,
                    validation_method,
                    max_bases,
                    percent_mapped,
                    summary,
                    verdict,
                    performed  )"
	  .

	  "values(?, ?, ?, ?, ?,?, ?, ?,?, ?,?,now() )";

	my $sth = $self->prepare($sql);
	$sth->bind_param( 1,  $genotype_results->other_id );
	$sth->bind_param( 2,  $genotype_results->table_name );
	$sth->bind_param( 3,  $genotype_results->reference );
	$sth->bind_param( 4,  $genotype_results->snps_bin );
	$sth->bind_param( 5,  $genotype_results->aligner );
	$sth->bind_param( 6,  $genotype_results->version );
	$sth->bind_param( 7,  $genotype_results->validation_method );
	$sth->bind_param( 8,  $genotype_results->max_bases );
	$sth->bind_param( 9,  $genotype_results->percent_mapped );
	$sth->bind_param( 10, $genotype_results->summary );
	$sth->bind_param( 11, $genotype_results->verdict );

	my $rows_inserted = $sth->execute();
	my $dbID          = $sth->{'mysql_insertid'};



	$genotype_results->dbID($dbID);
	$genotype_results->adaptor($self);
	$sth->finish();

	return $genotype_results;
}

sub fetch_by_aligner {
	my ( $self, $type ) = @_;
	my $sql =
	  "select " . $self->columns . " from genotype_results " . "where aligner = ? ";
	my $sth = $self->prepare($sql);
	$sth->bind_param( 1, $type );
	$sth->execute;
	my @results;
	while ( my $rowHashref = $sth->fetchrow_hashref ) {
		my $result = $self->object_from_hashref($rowHashref);
		push( @results, $result );
	}
	$sth->finish;
	return \@results;
}

sub fetch_by_reference {
	my ( $self, $type ) = @_;
	my $sql =
	  "select " . $self->columns . " from genotype_results where reference = ? ";
	my $sth = $self->prepare($sql);
	$sth->bind_param( 1, $type );
	$sth->execute;
	my @results;
	while ( my $rowHashref = $sth->fetchrow_hashref ) {
		my $result = $self->object_from_hashref($rowHashref);
		push( @results, $result );
	}
	$sth->finish;
	return \@results;
}

sub fetch_by_other_id {
	my ( $self, $type ) = @_;
	my $sql =
	  "select " . $self->columns . " from genotype_results where other_id = ? ";
	my $sth = $self->prepare($sql);
	$sth->bind_param( 1, $type );
	$sth->execute;
	my @results;
	while ( my $rowHashref = $sth->fetchrow_hashref ) {
		my $result = $self->object_from_hashref($rowHashref);
		push( @results, $result );
	}
	$sth->finish;
	return \@results;
}

sub fetch_by_verdict {
	my ( $self, $verdict ) = @_;
	my $sql =
	  "select " . $self->columns . " from genotype_results where verdict = ? ";
	my $sth = $self->prepare($sql);
	$sth->bind_param( 1, $verdict );
	$sth->execute;
	my @results;
	while ( my $rowHashref = $sth->fetchrow_hashref ) {
		my $result = $self->object_from_hashref($rowHashref);
		push( @results, $result );
	}
	$sth->finish;
	return \@results;
}

sub object_from_hashref {
	my ( $self, $hashref ) = @_;
	throw("Can't create a ReseqTrack::Genotype_Results from an undefined hashref")
	  if ( !$hashref );
	my $geno_obj = ReseqTrack::GenotypeResults->new(
		-adaptor           => $self,
		-genotype_results_id     => $hashref->{genotype_results_id},
		-other_id          => $hashref->{other_id},
		-table_name        => $hashref->{table_name},
		-reference         => $hashref->{reference},
		-snps_bin          => $hashref->{snps_bin},
		-aligner           => $hashref->{aligner},
		-version           => $hashref->{version},
		-validation_method => $hashref->{validation_method},
		-max_bases         => $hashref->{max_bases},
		-percent_mapped    => $hashref->{percent_mapped},
		-summary           => $hashref->{summary}, 
		-verdict           => $hashref->{verdict},
		-performed         => $hashref->{performed}
	);
	return $geno_obj;
}

1;
