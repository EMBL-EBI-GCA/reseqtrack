package ReseqTrack::DBSQL::VerifyBamIDAdaptor;

use strict;
use warnings;

use ReseqTrack::DBSQL::BaseAdaptor;
use ReseqTrack::Base;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename;
use Data::Dumper;
use ReseqTrack::VerifyBamID;

use vars qw(@ISA);



@ISA = qw(ReseqTrack::DBSQL::BaseAdaptor);

sub new {
	my ( $class, $db ) = @_;

	my $self = $class->SUPER::new($db);

	return $self;
}

sub columns {
	return
"verifybamid.file_id,
verifybamid.sample,
verifybamid.read_group,
verifybamid.chip_id,
verifybamid.snps,
verifybamid.num_reads,
verifybamid.avg_depth,
verifybamid.free_contam,
verifybamid.free_mlogl_est_contam,
verifybamid.free_mlogl_zero_contam,
verifybamid.free_ref_bias_ref_het,
verifybamid.free_ref_bias_refhomalt,
verifybamid.chip_contam,
verifybamid.chip_mlogl_est_contam,
verifybamid.chip_mlogl_zero_contam,
verifybamid.chip_ref_bias_ref_het,
verifybamid.chip_ref_bias_refhomalt,
verifybamid.depth_homref_site,
verifybamid.rel_depth_het_site,
verifybamid.rel_depth_homalt_site,
verifybamid.run_mode              ,
verifybamid.used_genotypes,
verifybamid.target_region         ,
verifybamid.vcf,
verifybamid.verdict,
verifybamid.performed";
}

sub table_name {
	return "verifybamid";
}

sub store {
	my ( $self, $verifybamid ) = @_;

	throw(  "Can't store "
		  . $verifybamid
		  . " using ReseqTrack::DBSQL::VerifyBamID" )
	  unless ( $verifybamid->isa("ReseqTrack::VerifyBamID") );


	my @columns = qw( file_id sample read_group chip_id
snps num_reads avg_depth free_contam free_mlogl_est_contam
free_mlogl_zero_contam free_ref_bias_ref_het free_ref_bias_refhomalt
chip_contam chip_mlogl_est_contam chip_mlogl_zero_contam
chip_ref_bias_ref_het chip_ref_bias_refhomalt depth_homref_site
rel_depth_het_site rel_depth_homalt_site run_mode
used_genotypes target_region vcf verdict);



	my $inset_cols = join( "\,", @columns );

	my $needed_q_marks = "?," x scalar @columns;
	chop $needed_q_marks;

	my $sql = "insert into verifybamid ( $inset_cols, performed) values( $needed_q_marks,now() )";



	my $sth = $self->prepare($sql);
	$sth->bind_param( 1, $verifybamid->file_id);
	$sth->bind_param( 2, $verifybamid->sample);
	$sth->bind_param( 3, $verifybamid->read_group);
	$sth->bind_param( 4, $verifybamid->chip_id);
	$sth->bind_param( 5, $verifybamid->snps);
	$sth->bind_param( 6, $verifybamid->num_reads);
	$sth->bind_param( 7, $verifybamid->avg_depth);
	$sth->bind_param( 8, $verifybamid->free_contam);
	$sth->bind_param( 9, $verifybamid->free_mlogl_est_contam);
	$sth->bind_param( 10, $verifybamid->free_mlogl_zero_contam);
	$sth->bind_param( 11, $verifybamid->free_ref_bias_ref_het);
	$sth->bind_param( 12, $verifybamid->free_ref_bias_refhomalt);
	$sth->bind_param( 13, $verifybamid->chip_contam);
	$sth->bind_param( 14, $verifybamid->chip_mlogl_est_contam);
	$sth->bind_param( 15, $verifybamid->chip_mlogl_zero_contam);
	$sth->bind_param( 16, $verifybamid->chip_ref_bias_ref_het);
	$sth->bind_param( 17, $verifybamid->chip_ref_bias_refhomalt);
	$sth->bind_param( 18, $verifybamid->depth_homref_site);
	$sth->bind_param( 19, $verifybamid->rel_depth_het_site);
	$sth->bind_param( 20, $verifybamid->rel_depth_homalt_site);
	$sth->bind_param( 21, $verifybamid->run_mode);
	$sth->bind_param( 22, $verifybamid->used_genotypes);
	$sth->bind_param( 23, $verifybamid->target_region);
	$sth->bind_param( 24, $verifybamid->vcf);
	$sth->bind_param( 25, $verifybamid->verdict);

	my $rows_inserted = $sth->execute();
	my $dbID          = $sth->{'mysql_insertid'};

	$verifybamid->dbID($dbID);
	$verifybamid->adaptor($self);
	$sth->finish();

	return $verifybamid;

}

sub update {
	my ( $self, $verifybamid ) = @_;



	throw(  "Can't update "
		  . $verifybamid
		  . " using ReseqTrack::DBSQL::VerifyBamID" )
	  unless ( $verifybamid->isa("ReseqTrack::VerifyBamID") );

	my @cols = qw(  sample  chip_id snps num_reads avg_depth free_contam free_mlogl_est_contam
free_mlogl_zero_contam free_ref_bias_ref_het free_ref_bias_refhomalt
chip_contam chip_mlogl_est_contam chip_mlogl_zero_contam
chip_ref_bias_ref_het chip_ref_bias_refhomalt depth_homref_site
rel_depth_het_site rel_depth_homalt_site run_mode
used_genotypes target_region verdict);

	my $update_cols = join (" = ? ," , @cols);
        $update_cols .= " = ?"; 

	my $sql = "update verifybamid set $update_cols , performed = now() where ( file_id = ? and read_group = ? and vcf = ?) ";

	my $sth = $self->prepare($sql);
	my $ctr = 0;

	$sth->bind_param( 1, $verifybamid->sample);
	$sth->bind_param( 2, $verifybamid->chip_id);
	$sth->bind_param( 3, $verifybamid->snps);
	$sth->bind_param( 4, $verifybamid->num_reads);
	$sth->bind_param( 5, $verifybamid->avg_depth);
	$sth->bind_param( 6, $verifybamid->free_contam);
	$sth->bind_param( 7, $verifybamid->free_mlogl_est_contam);
	$sth->bind_param( 8, $verifybamid->free_mlogl_zero_contam);
	$sth->bind_param( 9, $verifybamid->free_ref_bias_ref_het);
	$sth->bind_param( 10, $verifybamid->free_ref_bias_refhomalt);
	$sth->bind_param( 11, $verifybamid->chip_contam);
	$sth->bind_param( 12, $verifybamid->chip_mlogl_est_contam);
	$sth->bind_param( 13, $verifybamid->chip_mlogl_zero_contam);
	$sth->bind_param( 14, $verifybamid->chip_ref_bias_ref_het);
	$sth->bind_param( 15, $verifybamid->chip_ref_bias_refhomalt);
	$sth->bind_param( 16, $verifybamid->depth_homref_site);
	$sth->bind_param( 17, $verifybamid->rel_depth_het_site);
	$sth->bind_param( 18, $verifybamid->rel_depth_homalt_site);
	$sth->bind_param( 19, $verifybamid->run_mode);
	$sth->bind_param( 20, $verifybamid->used_genotypes);
	$sth->bind_param( 21, $verifybamid->target_region);
	$sth->bind_param( 22, $verifybamid->verdict);
	$sth->bind_param( 23, $verifybamid->file_id );
	$sth->bind_param( 24, $verifybamid->read_group );
	$sth->bind_param( 25, $verifybamid->vcf );
	
	$sth->execute();
	$sth->finish();

	return $verifybamid;

}

sub object_from_hashref {
	my ( $self, $hashref ) = @_;

	throw(
		"Can't create a ReseqTrack::VerifyBamID from an undefined hashref
	  "
	) if ( !$hashref );

	my $OBJ = "ReseqTrack::VerifyBamID";
	my $obj = $OBJ->new(
			    -adaptor                => $self,
			    -dbID                   => $hashref->{verifybamid_id},
			    -file_id => $hashref->{file_id},
			    -sample => $hashref->{sample},
			    -read_group => $hashref->{read_group},
			    -chip_id => $hashref->{chip_id},
			    -snps => $hashref->{snps},
			    -num_reads => $hashref->{num_reads},
			    -avg_depth => $hashref->{avg_depth},
			    -free_contam => $hashref->{free_contam},
			    -free_mlogl_est_contam => $hashref->{free_mlogl_est_contam},
			    -free_mlogl_zero_contam => $hashref->{free_mlogl_zero_contam},
			    -free_ref_bias_ref_het => $hashref->{free_ref_bias_ref_het},
			    -free_ref_bias_refhomalt => $hashref->{free_ref_bias_refhomalt},
			    -chip_contam => $hashref->{chip_contam},
			    -chip_mlogl_est_contam => $hashref->{chip_mlogl_est_contam},
			    -chip_mlogl_zero_contam => $hashref->{chip_mlogl_zero_contam},
			    -chip_ref_bias_ref_het => $hashref->{chip_ref_bias_ref_het},
			    -chip_ref_bias_refhomalt => $hashref->{chip_ref_bias_refhomalt},
			    -depth_homref_site => $hashref->{depth_homref_site},
			    -rel_depth_het_site => $hashref->{rel_depth_het_site},
			    -rel_depth_homalt_site => $hashref->{rel_depth_homalt_site},
			    -run_mode => $hashref->{run_mode},
			    -used_genotypes => $hashref->{used_genotypes},
			    -target_region => $hashref->{target_region},
			    -vcf => $hashref->{vcf},
			    -verdict => $hashref->{verdict},
			    -performed => $hashref->{performed},
	);

	return $obj;
}


sub fetch_by_vcf_file_id_and_readgroup {
    
  my ($self, $vcf, $file_id, $read_group) = @_;
  my $sql = "select ". $self->columns." from " .  $self->table_name .
    " where vcf=? and file_id = ? and read_group = ?";
          
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $vcf);
  $sth->bind_param(2, $file_id);
  $sth->bind_param(3, $read_group);
  $sth->execute;
  
  my @results;
  while (my $rowHashref = $sth->fetchrow_hashref) {
    my $result = $self->object_from_hashref($rowHashref);
    push(@results, $result);
  }
  $sth->finish;

  throw("file_id ",$file_id." has returned multiple objects ".@results." not sure what to do") 
    if (@results && @results >= 2);

  if ( @results ==0 ){
#   print "No results found for file_id=  $file_id  read_group $read_group\n";
   return;

  }

  my $result = $results[0];
    
  return $result;
}


=head
sub fetch_by_file_id_and_readgroup{
    
  my ($self, $file_id, $read_group) = @_;
  my $sql = "select ". $self->columns." from " .  $self->table_name .
    " where file_id = ? and read_group = ?";
          
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $file_id);
  $sth->bind_param(2, $read_group);
  $sth->execute;
  
  my @results;
  while (my $rowHashref = $sth->fetchrow_hashref) {
    my $result = $self->object_from_hashref($rowHashref);
    push(@results, $result);
  }
  $sth->finish;

  throw("file_id ",$file_id." has returned multiple objects ".@results." not sure what to do") 
    if (@results && @results >= 2);

  if ( @results ==0 ){
#   print "No results found for file_id=  $file_id  read_group $read_group\n";
   return;

  }

  my $result = $results[0];
    
  return $result;
}

=cut

return 1;

