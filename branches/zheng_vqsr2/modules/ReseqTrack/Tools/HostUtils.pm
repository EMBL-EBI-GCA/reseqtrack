=pod

=head1 NAME

ReseqTrack::Tools::HostUtils;

=head1 SYNOPSIS

A method collection class which contains useful utility methods for ReseqTrack::Host
object handling.

You can either import all the contained method into your name space

use ReseqTrack::Tools::HostUtils;

or just specific methods

use ReseqTrack::Tools::HostUtils qw(create_host_object)

=cut

package ReseqTrack::Tools::HostUtils;

use strict;
use warnings;
use Exporter;
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Host;
use ReseqTrack::Tools::GeneralUtils;
use vars qw (@ISA  @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(get_host_object check_name);



=head2 get_host_object

  Arg [1]   : string, name of host 
  Arg [2]   : ReseqTrack::DBSQL::DBAdaptor
  Arg [3]   : int, 0/1 to indicate if the host is remote
  Arg [4]   : int, 0/1 to indicate that the host name does not needed to be checked  
  Function  : create a new host object or fetch an existing one from the db
  Returntype: ReseqTrack::Host
  Exceptions: throws if name is not of defined type
  Example   : my $host = get_host_object("1000genomes.ebi.ac.uk", $db);

=cut



sub get_host_object{
  my ($name, $db, $remote, $no_check) = @_;
  unless($no_check){
    unless(check_name($name)){
      throw($name." is not a valid name for a host in ".$db->dbc->dbname);
    }
  }
  my $ha = $db->get_HostAdaptor;
  my $host = $ha->fetch_by_name($name);
  unless($host){
    $host = ReseqTrack::Host->new
      (
       -name => $name,
       -remote => $remote
      );
  }
  return $host;
}



=head2 check_name

  Arg [1]   : String, name
  Function  : confirms name is part of allowed set
  Returntype: 0/1 depending if it is correct
  Exceptions: 
  Example   : 

=cut



sub check_name{
  my ($name) = @_;
  my %valid_hash;
  $valid_hash{"1000genomes.ebi.ac.uk"} = 1;
  $valid_hash{"sanger"} = 1;
  $valid_hash{"tgen"} = 1;
  $valid_hash{"ncbi"} = 1;
  $valid_hash{"broad"} = 1;
  $valid_hash{"baylor"} = 1;
  $valid_hash{"boston_college"} = 1;
  return $valid_hash{$name};
}

1;
