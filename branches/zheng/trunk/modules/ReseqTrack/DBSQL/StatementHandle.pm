#This is a complete crib from Bio::EnsEMBL::DBSQL::StatementHandle;

package ReseqTrack::DBSQL::StatementHandle;

use vars qw(@ISA);
use strict;

use ReseqTrack::Tools::Exception qw(warning stack_trace_dump);

use DBD::mysql;
use DBI;

#use Time::HiRes qw(time);

@ISA = qw(DBI::st);


# As DBD::mysql::st is a tied hash can't store things in it,
# so have to have parallel hash
my %dbchash;
my %dbc_sql_hash;


sub dbc {
  my $self = shift;

  if (@_) {
    my $dbc = shift;
    if(!defined($dbc)) {
      # without delete key space would grow indefinitely causing mem-leak
      delete($dbchash{$self});
    } else {
      $dbchash{$self} = $dbc;
    }
  }

  return $dbchash{$self};
}

sub sql {
  my $self = shift;

  if (@_) {
    my $sql = shift;
    if(!defined($sql)) {
      # without delete key space would grow indefinitely causing mem-leak
      delete($dbc_sql_hash{$self});
    } else {
      $dbc_sql_hash{$self} = $sql;
    }
  }

  return $dbc_sql_hash{$self};
}

sub DESTROY {
  my ($self) = @_;

  my $dbc = $self->dbc;
  $self->dbc(undef);
  my $sql = $self->sql;
  $self->sql(undef);

  # Re-bless into DBI::st so that superclass destroy method is called if
  # it exists (it does not exist in all DBI versions).
  bless( $self, 'DBI::st' );

  # The count for the number of kids is decremented only after this
  # function is complete. Disconnect if there is 1 kid (this one)
  # remaining.
  if (    $dbc
       && $dbc->disconnect_when_inactive()
       && $dbc->connected
       && ( $dbc->db_handle->{Kids} == 1 ) )
  {
    if ( $dbc->disconnect_if_idle() ) {
      warn("Problem disconnect $self around sql = $sql\n");
    }
  }
} ## end sub DESTROY

1;
