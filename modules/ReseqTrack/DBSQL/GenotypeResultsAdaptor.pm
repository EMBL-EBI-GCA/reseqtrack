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
use Data::Dumper;

@ISA = qw(ReseqTrack::DBSQL::BaseAdaptor);

sub new {
	my ( $class, $db ) = @_;

	my $self = $class->SUPER::new($db);

	return $self;
}


sub columns {
	return "
        genotype_results.genotype_results_id,
        genotype_results.table_name,
        genotype_results.other_id, 
        genotype_results.name,
        genotype_results.claimed,
        genotype_results.top_hit,
        genotype_results.second_hit,
        genotype_results.ratio_2_to_1,
        genotype_results.ratio_claimed,
        genotype_results.reference  ,
        genotype_results.snps_bin  ,
        genotype_results.aligner,
        genotype_results.version,
        genotype_results.validation_method,
        genotype_results.max_bases,
        genotype_results.percent_mapped  , 
        genotype_results.percent_reads_used,
        genotype_results.verdict    ,
        genotype_results.cfg_file,
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
       
        table_name,
        other_id, 
        name,
        claimed,
        top_hit,
        second_hit,
        ratio_2_to_1,
        ratio_claimed,
        reference  ,
        snps_bin  ,
        aligner,
        version,
        validation_method,
        max_bases,
        percent_mapped  , 
        percent_reads_used,
        verdict    ,
        cfg_file,
        performed        
       )"
	  .

	  "values(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,now() )";

	my $sth = $self->prepare($sql);

	$sth->bind_param(1 ,  $genotype_results->table_name );
	$sth->bind_param(2 ,  $genotype_results->other_id );
	$sth->bind_param(3 ,  $genotype_results->name);
	$sth->bind_param(4 ,  $genotype_results->claimed);
	$sth->bind_param(5 ,  $genotype_results->top_hit);
	$sth->bind_param(6 ,  $genotype_results->second_hit);
	$sth->bind_param(7 ,  $genotype_results->ratio_2_to_1);
	$sth->bind_param(8 ,  $genotype_results->ratio_claimed);
	$sth->bind_param(9 ,  $genotype_results->reference );
	$sth->bind_param(10,  $genotype_results->snps_bin );
	$sth->bind_param(11 ,  $genotype_results->aligner );
	$sth->bind_param(12 ,  $genotype_results->version );
	$sth->bind_param(13 ,  $genotype_results->validation_method );
	$sth->bind_param(14 ,  $genotype_results->max_bases );
	$sth->bind_param(15 ,  $genotype_results->percent_mapped );
        $sth->bind_param(16 ,  $genotype_results->percent_reads_used );
	$sth->bind_param(17 ,  $genotype_results->verdict );
        $sth->bind_param(18 ,  $genotype_results->cfg_file );  

	my $rows_inserted = $sth->execute();
	my $dbID          = $sth->{'mysql_insertid'};



	$genotype_results->dbID($dbID);
	$genotype_results->adaptor($self);
	$sth->finish();

	return $genotype_results;
}

sub object_from_hashref {
  my ( $self, $hashref ) = @_;

  throw("Can't create a ReseqTrack::GenotypeResults from an undefined hashref")
    if ( !$hashref );


  my $OBJ ="ReseqTrack::GenotypeResults";
  my $geno_obj = $OBJ->new(
			   -adaptor           => $self,
			   -genotype_results_id  => $hashref->{genotype_results_id},
			   -table_name        => $hashref->{table_name},
			   -other_id          => $hashref->{other_id},
			   -name              => $hashref->{name},
			   -claimed           => $hashref->{claimed},
			   -top_hit           => $hashref->{top_hit},
			   -second_hit        => $hashref->{second_hit},
			   -ratio_2_to_1      => $hashref->{ratio_2_to_1},
			   -ratio_claimed     => $hashref->{ratio_claimed},
			   -reference         => $hashref->{reference},
			   -snps_bin          => $hashref->{snps_bin},
			   -aligner           => $hashref->{aligner},
			   -version           => $hashref->{version},
			   -validation_method => $hashref->{validation_method},
			   -max_bases         => $hashref->{max_bases},
			   -percent_mapped    => $hashref->{percent_mapped},
			   -percent_reads_used=> $hashref->{percent_reads_used},
			   -verdict           => $hashref->{verdict},
			   -cfg_file          => $hashref->{cfg_file},
			   -performed         => $hashref->{performed},
	);
	return $geno_obj;
}

sub fetch_by_name {
	my ( $self, $name ) = @_;
	my $sql =
	  "select " . $self->columns . " from genotype_results " . "where name = ? ";
	my $sth = $self->prepare($sql);
	$sth->bind_param( 1, $name );
	$sth->execute;
	my @results;
	while ( my $rowHashref = $sth->fetchrow_hashref ) {
	  my $result = $self->object_from_hashref($rowHashref);
	  push( @results, $result );
	}
	$sth->finish;

	throw($name." has returned multiple objects ".@results." not sure what to do") 
	   if ( scalar (@results) > 1);

	my $gt_obj = $results[0];
	return $gt_obj;
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

sub update{

	my ( $self, $genotype_results ) = @_;
	throw(  "Can't update "
		  . $genotype_results
		  . " using ReseqTrack::DBSQL::GenotypeResultsAdaptor" )
	  unless ( $genotype_results->isa("ReseqTrack::GenotypeResults") );


#	print Dumper  ($genotype_results);
	
	my $sql = "update genotype_results  set table_name = ?, other_id   = ? , claimed  = ?, top_hit = ?, second_hit = ? , ratio_2_to_1 = ?, ratio_claimed = ?, reference  = ? , snps_bin     = ? ,aligner    = ?, version    = ? , validation_method = ?, max_bases     = ?, percent_mapped = ? , percent_reads_used = ?,  verdict = ? ,  cfg_file = ?,  performed = now()  where name  = ? ";
	print "$sql\n";
#	exit;
	my $sth = $self->prepare($sql);

	$sth->bind_param( 1 ,  $genotype_results->table_name );
	$sth->bind_param( 2 ,  $genotype_results->other_id );
	$sth->bind_param( 3 ,  $genotype_results->claimed);
	$sth->bind_param( 4 ,  $genotype_results->top_hit);
	$sth->bind_param( 5 ,  $genotype_results->second_hit);
	$sth->bind_param( 6 ,  $genotype_results->ratio_2_to_1);
	$sth->bind_param( 7 ,  $genotype_results->ratio_claimed);
	$sth->bind_param( 8 ,  $genotype_results->reference );
	$sth->bind_param( 9 ,  $genotype_results->snps_bin );
	$sth->bind_param(10 ,  $genotype_results->aligner );
	$sth->bind_param(11 ,  $genotype_results->version );
	$sth->bind_param(12 ,  $genotype_results->validation_method );
	$sth->bind_param(13 ,  $genotype_results->max_bases );
	$sth->bind_param(14 ,  $genotype_results->percent_mapped );
        $sth->bind_param(15 ,  $genotype_results->percent_reads_used );
	$sth->bind_param(16 ,  $genotype_results->verdict );
        $sth->bind_param(17 ,  $genotype_results->cfg_file );  
	$sth->bind_param(18 ,  $genotype_results->name);

        $sth->execute();
	$sth->finish();

	return $genotype_results;
}


1;
