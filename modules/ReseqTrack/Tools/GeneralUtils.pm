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
	     create_filename calculate_coverage trim_spaces execute_system_command
	     execute_pipe_system_command get_params get_open_file_handle get_time_stamps);




=head2 useage

  Arg [1]   : N/A
  Function  : executes perldocs for given script
  Returntype: N/A
  Exceptions: 
  Example   : 

=cut


sub useage{
  exec('perldoc', $0);
  exit(0);
}



=head2 current_time

  Arg [1]   : n/a
  Function  : returns the current time in the from yyyy-mm-dd hh:mm:ss
  Returntype: see above
  Exceptions: none
  Example   : my $time = current_time

=cut



sub current_time{
  my ($second, $minute, $hour, $day, $month, $year) = localtime();
  return sprintf("%04d-%02d-%02d %02d:%02d:%02d", $year+1900, $month+1, $day, $hour, $minute, $second);
}

=head2 get_time_stamps

  Arg [1]   : n/a
  Function  : returns three time stamps, yyyy-mm-dd.hh.mm, yyyy-mm, or yyyy-mm-dd
  Returntype: see above
  Exceptions: none
  Example   : my ($time_stamp, $month_stamp, $day_stamp) = get_time_stamps

=cut

sub get_time_stamps {
        my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime();
        $mon++;
        $year = $year + 1900;
        my $time = $year . "-" . $mon . "-" . $mday . "." . $hour . "." . $min;
        my $month = $year . "-" . $mon ;
        my $day = $year . "-" . $mon . "-" . $mday;
        return ($time, $month, $day);
}


=head2 current_date

  Arg [1]   : n/a
  Function  : returns the current time in the from yyyymmdd
  Returntype: see above
  Exceptions: none
  Example   : my $time = current_time

=cut



sub current_date{
  my (%options) = @_;
  my ($year, $month, $day) = (localtime())[5,4,3];
  if ($options{hyphens}) {
    return sprintf("%04d-%02d-%02d", $year+1900, $month+1, $day);
  }
  return sprintf("%04d%02d%02d", $year+1900, $month+1, $day);
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
    $values[0] =~ s/\s+$//;
    $values[1] =~ s/\s+$//;
    $values[0] =~ s/^\s+//;
    $values[1] =~ s/^\s+//;
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
    print STDERR "$name in place\n\n".
                     "\tdb       $dbname\@$dbhost\n".
                     "\tpid      $pid on ".
                     "host $host\n\tstarted  $started\n\n".
                     "There may be another process running\n".
                     "If the process does not exist, ".
                     "remove the lock by removing the lock from the ".
                     "database:\n\ndelete from meta where ".
                     "meta_key = '$name';\n\n\n\n" ;
    my $error_str =  "There may be another process running; " .
                     "it must be terminated before this script can be run.\n" .
                     "If the process does not exist, ".
                     "remove the lock by removing the lock from the ".
                     "database:\n\ndelete from meta where ".
                     "meta_key = '$name';\n\n\n\n" ;
    throw("There may be another process running");
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

=head2 execute_system_command

  Arg [1]   : string containing the command line
  Function  : executes the command line and waits for the exit code
  Returntype: exit code of system
  Exceptions: throws if command failed to execute or if exit code is not zero.
  Also throws if the command line is an empty string.
  Example   : execute_system_command('/path/to/executable -options > output')

=cut

sub execute_system_command{
    my $command_line = shift;

    throw("no command line to execute")
        if (! $command_line);

    system( $command_line );

    throw("command failed: $! $command_line")
        if ( $? == -1 );

    my $signal = $? & 127;
    throw("process died with signal $signal $command_line")
        if ($signal);

    my $exit = $? >> 8;
    throw("command exited with value $exit $command_line")
        if ($exit != 0);

    return $exit;

}

=head2 execute_pipe_system_command

  Arg [1]   : string containing the command line
  Function  : executes the command line and collects the first line of output. Best suited to simple output where creating a file is unecessary
  Returntype: scalar of the first line of output, line end removed
  Exceptions: throws if command failed to execute or if exit code is not zero.
  Also throws if the command line is an empty string.
  Example   : my $line_count = execute_pipe_system_command('cat /path/to/file | grep -v pattern_to_filter_out | wc -l')

=cut
sub execute_pipe_system_command{
	my ($command_line, $verbose) = @_;
		
	throw("no command to execute") if (! $command_line);

	my $result;
	
	open (my $fh, '-|', $command_line) or throw("command failed to open: $! $command_line");
	
	$result = <$fh>;
	chomp $result;
		
	throw("command failed: $! $command_line")
        if ( $? == -1 );

    my $signal = $? & 127;
    throw("process died with signal $signal $command_line")
        if ($signal);

    my $exit = $? >> 8;
    throw("command exited with value $exit $command_line")
        if ($exit != 0);	
	
	close $fh;
	 
	return $result;
}

=head2  get_params

   Arg [1]   : hash ref
   Function  : adds key-value pairs listed in input file to hash unless already 
               defined
   Returntype: hash ref
   Exceptions: none
   Example   : get_params ( cfg_file, \%input);

   Note_1    : If the same key appears several times in the config file only the 
               first one is used
   Note_2    : If a key has been defined before in $input hash, the value for this 
               key in the config file is skipped

=cut

sub get_params {

    my $file  = shift;
    my $input = shift;


    throw ("Could not open $file") if ( !-e $file);

    open my $IN, '<', $file || throw "No config file found";

    while (<$IN>) {
        chomp $_;
        next if ( !$_ );
        next if ( /^#/ );
        my @aa = split /=/;
        $aa[1] =~ s/\s+$//g;

        if (defined  $$input{ $aa[0] }){
          print "$aa[0] already set. Skipping cfg entry\n";
          next;
        }

        $$input{ $aa[0] } = $aa[1] ;
    }
    close $IN;

    return ( $input );
}


=head2 get_open_file_handle

  Arg [1]   : file name
  Function  : check file extension and open file handle
  Returntype: file handle
  Exceptions: none
  Example   : get_open_file_handle ( file_name);

=cut


sub get_open_file_handle {

  my $filename = shift;
  my $fh_in;

  throw "No file found: $filename" if (! -e $filename );

  if ( $filename =~ /\.gz$/i ) {
    open( $fh_in, "zcat $filename | " )
      or throw("failed zcat $filename: $!");
  }
  else {
    open( $fh_in, '<', $filename) or throw("failed open $filename: $!");
  }

  return $fh_in;
}




1;
