=pod

=head1 NAME

ReseqTrack::Tools::GeneralUtils;

=head1 SYNOPSIS

This is a collection of generally useful methods for the ReseqTrack code base

=head1 Examples

my $time = current_time();

=cut

package ReseqTrack::Tools::GeneralUtils;

use strict;
use warnings;
use Exporter;
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use Sys::Hostname;
use Socket;

use vars qw (@ISA  @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(current_time parse_movelist get_input_arg create_lock_string
             delete_lock_string is_locked useage convert_to_giga current_date 
	     create_filename calculate_coverage trim_spaces);



=head2 current_time

  Arg [1]   : n/a
  Function  : returns the current time in the from yyyy-mm-dd hh:mm:ss
  Returntype: see above
  Exceptions: none
  Example   : my $time = current_time

=cut



sub current_time{
  my @months = qw(01 02 03 04 05 06 07 08 09 10 11 12);
  my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
  my $year = 1900 + $yearOffset;
  $dayOfMonth = sprintf("%02.d", $dayOfMonth);
  $second = sprintf("%02.d", $second);
  $minute = sprintf("%02.d", $minute);
  $hour = sprintf("%02.d", $hour);
  my $theTime = "$year-$months[$month]-$dayOfMonth $hour:$minute:$second";
  return $theTime;
}

=head2 current_date

  Arg [1]   : n/a
  Function  : returns the current time in the from yyyymmdd
  Returntype: see above
  Exceptions: none
  Example   : my $time = current_time

=cut



sub current_date{
  my @months = qw(01 02 03 04 05 06 07 08 09 10 11 12);
  my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
  my $year = 1900 + $yearOffset;
  $dayOfMonth = sprintf("%02.d", $dayOfMonth);
  my $theTime = $year.$months[$month].$dayOfMonth;
  return $theTime;
}


=head2 parse_movelist

  Arg [1]   : string, file_path
  Function  : parses file and returns hash putting column 1 as the keys and column 2 as the values 
  Returntype: hashref
  Exceptions: throws on failure to open file
  Example   : my %move_list = %{parse_movelist($file)};

=cut



sub parse_movelist{
  my ($move_list) = @_;
  open(FH, $move_list) or throw("FAiled to open ".$move_list." $!");
  my %hash;
  while(<FH>){
    chomp;
    my @values = split /\t/, $_;
    #print "Adding old ".$values[0]." new ".$values[1]." to the list\n";
    $hash{$values[0]} = $values[1];
  }
  return \%hash;
}


=head2 get_input_arg

  Function  : waits for input from STDIN and returns '1' if input =~m/y/i
              and '0' if input matches /n/i.
  Returntype: 1 or 0
  Exceptions: none

=cut

sub get_input_arg {
  while (defined (my $line=<STDIN>)){
   chomp($line) ;
   if ( $line=~m/y/i){
      return 1 ;
   }elsif( $line =~m/n/i){
     return 0 ;
   }
   print "Wrong input - only answer 'y' or 'n'\n" ;
  }
}

sub create_lock_string{
  my ($name, $ma) = @_;
  my $host = qualify_hostname(hostname());
  my $user = scalar getpwuid($<);
  my $lock_str = join ":", "$user\@$host", $$, time();

  if($ma){
    $ma->store($name, $lock_str);
  }

  return $lock_str;

}

sub delete_lock_string{
  my ($name, $ma) = @_;
  #print "Removing lock for ".$name."\n";
  $ma->remove_by_meta_key($name);
}

sub is_locked{
  my ($name, $ma) = @_;
  my $string = $ma->fetch_meta_value_by_meta_key($name);
  if($string){
    #laura@ebi-231.ebi.ac.uk:27736:1267806066 
    my($user, $host, $pid, $started) = $string =~ /(\w+)@(.*):(\d+):(\d+)/;
    $started = scalar localtime $started;
    my $dbname = $ma->dbc->dbname;
    my $dbhost = $ma->dbc->host;
    my $error_str = ("Error: $name in place\n\n".
                     "\tdb       $dbname\@$dbhost\n".
                     "\tpid      $pid on ".
                     "host $host\n\tstarted  $started\n\n".
                     "The process above must be terminated before this ".
                     "script can be run.\nIf the process does not exist, ".
                     "remove the lock by removing the lock from the ".
                     "database:\n\ndelete from meta where ".
                     "meta_key = '$name';\n\n\n\n" );
    print STDERR $error_str;
    throw("Can't run there may be another process running");
  }
}

=head2 qualify_hostname

  Arg [1]   : string, output of Sys::Hostname::hostname
  Function  : produces fully qualified host name
  Returntype: string
  Exceptions: none
  Example   : my $host = $self->qualify_hostname(hostname());

=cut


sub qualify_hostname{
  my ($hostname) = @_;

  my $addr = gethostbyname($hostname);
  my $host = gethostbyaddr($addr, AF_INET);
  return $host;
}


sub useage{
  exec('perldoc', $0);
  exit(0);
}

sub convert_to_giga{
  my ($base_count) = @_;
  $base_count = 0 unless($base_count);
 # my $gigabase = $base_count/1000000000;
  my $gigabase = sprintf("%.2f", $base_count/1000000000);
  return $gigabase;
}

sub calculate_coverage{
  my ($base_count, $genome_size) = @_;
  $genome_size = 2850000000 unless($genome_size);
  return $base_count unless($base_count && $base_count >= 1);
  return $base_count/$genome_size;
}

sub create_filename{
  my ($dir, $stem, $ext) = @_;
  my $rand = int(rand(10000));
  my $name = $stem.".".$$.".".$rand;
  $name .= $ext if($ext);
  my $path = $dir."/".$name;
  while(-e $path){
    $rand = int(rand(10000));
    my $name = $stem.".".$$.".".$rand;
    $name .= $ext if($ext);
    $path = $dir."/".$name;
  }
  $path =~ s/\/\//\//g;
  return $path;
}

sub trim_spaces{
  my $string = shift;
  $string =~ s/^\s+//;
  $string =~ s/\s+$//;
  return $string;
}


1;
