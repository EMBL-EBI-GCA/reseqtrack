#!/sw/arch/bin/perl
use strict;
use ReseqTrack::Event;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception;
use Getopt::Long;
use Data::Dumper;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport;
my $dbname;

my $infile;
my $outfile = "event_table.dat";    #default name of output file
my $run;

my ( $read, $write, $update, $load, $help, $verbose, );

&GetOptions(
    'dbhost=s'  => \$dbhost,
    'dbname=s'  => \$dbname,
    'dbuser=s'  => \$dbuser,
    'dbpass=s'  => \$dbpass,
    'dbport=s'  => \$dbport,
    'file=s'    => \$infile,
    'read!'     => \$read,
    'write!'    => \$write,
    'help!'     => \$help,
    'run!'      => \$run,
    'verbose!'  => \$verbose,
    'outfile=s' => \$outfile,
);

if ($help) {
  perldocs();
}

#Could allow this. Add to db then dump contents. Always dump contents
#after modifying db.
if ( $read && $write ) {
    throw "Must specify -read or -write. Not both\n";
}

if ( ! $read && ! $write ){
     	 print "Must specify -read or -write.\n";
	 perldocs();
} 

if ( $read && !$infile ) {
    throw "Must specify file (-file) with list of events when in read mode";
}

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
    -host   => $dbhost,
    -user   => $dbuser,
    -port   => $dbport,
    -dbname => $dbname,
    -pass   => $dbpass,
);

if ($read) {
    my $events_list = parse_events_list( $infile, $verbose );
    process_event_from_file( $events_list, $verbose );
}

if ($write) {
    dump_event_table_contents( $outfile, $verbose );
}

#######################################################
sub dump_event_table_contents {

    #output contents of 'event' table to fomatted text file

    my $outfile = shift;
    my $verbose = shift;

    my $ea     = $db->get_EventAdaptor;
    my $events = $ea->fetch_all;

    open( OUT, '>', "event_table.dat" ) or die "File open error: $!";

    printf "Event data entries:%d \n", scalar @$events;

    foreach my $ievent (@$events) {
        printf OUT "[$ievent->{name}]\n"          if $ievent->{name};
        printf OUT "program=$ievent->{program}\n" if $ievent->{program};
        printf OUT "program_version=$ievent->{program_version}\n"
          if $ievent->{program_version};
        printf OUT "options=$ievent->{options}\n" if $ievent->{options};
        printf OUT "input_flag=$ievent->{input_flag}\n"
          if $ievent->{input_flag};
        printf OUT "farm_options=$ievent->{farm_options}\n"
          if $ievent->{farm_options};
        printf OUT "batch_size=$ievent->{batch_size}\n"
          if $ievent->{batch_size};
        printf OUT "output_path=$ievent->{output_path}\n"
          if $ievent->{output_path};
        printf OUT "type=$ievent->{type}\n" if $ievent->{type};
        printf OUT "table_name=$ievent->{table_name}\n"
          if $ievent->{table_name};
        printf OUT "\n";
    }

    close(OUT);
    print "\nFinished. See $outfile file\n\n";
}

sub process_event_from_file {

    #loop through array of hashes and load/update event info
    # into database

    my $events  = shift;     #array of hashes with event info
    my $verbose = shift;
    my $class   = "Event";
    my $updated = 0;
    my $stored  = 0;

    printf "\nFiles to process:%d \n", scalar @$events if $verbose;

    my $ea = $db->get_EventAdaptor;

    foreach my $ievent (@$events) {

        print "Processing: " . $ievent->{name} . "\n" if $verbose;

        print Dumper ($ievent) if $verbose;

        my %params;
        foreach my $key ( keys(%$ievent) ) {
            my $new = "-" . $key;
            $params{$new} = $ievent->{$key};
        }

        my $jevent = ReseqTrack::Event->new( %params, );

        print Dumper ($jevent) if $verbose;

        my $stored_event = $ea->fetch_by_name( $jevent->name );

        if ($stored_event) {
            $jevent->dbID( $stored_event->dbID );
            $jevent->adaptor($ea);
            $ea->update($jevent);
            $updated++;
        }
        else {
            $ea->store($jevent);
            $stored++;
        }
    }

    print "Events stored : $stored\n";
    print "Events updated: $updated\n";
}

sub parse_events_list {

    # assumming properly formatted file.
    # Parse file and create an array of hashes containing file event
    # information

    my $filename = shift;    #formatted list of event data
    my $verbose  = shift;

    my $this_file = ();
    my $ctr       = 0;
    my $name      = "name";
    my @array_events;

    open( INFILE, $filename ) or die "File designation error: $!";

    while (<INFILE>) {

        $ctr++;
        chomp($_);

        #if blank line should be at end of event data
        if ( $_ =~ /^$/ ) {

            if ( $this_file && $this_file->{$name} ) {
                print "Storing $this_file->{$name}\n" if $verbose;
                push( @array_events, $this_file );
            }

            $this_file = {};
            next;    #blank line
        }

        #extract filename from between []
        if ( $_ =~ /\[(.*?)\]/ ) {
            my $header = $1;    # file name
            $this_file->{$name} = $header;
            print "\[$this_file->{$name}\]\n" if $verbose;
            next;
        }

        #should not have data with no file name associated
        if ( !( $_ =~ /^$/ ) && !$this_file->{name} ) {
            warning("Unassociated data line $ctr.Ignoring");
            next;
        }

        if ( $_ =~ /\=/ ) {
            ( my $key, my $value ) = split /\=/, $_;

            if ( !exists $this_file->{$key} ) {
                $this_file->{$key} = $value;
                print $key. "=" . $this_file->{$key} . "\n" if $verbose;
            }
            else {
                warning("Duplicate key data. line $ctr.Ignoring");
            }
        }

    }

    #store last entry
    if ( $this_file && $this_file->{$name} ) {
        print "Storing $this_file->{$name}\n" if $verbose;
        push( @array_events, $this_file );
        $this_file = ();
    }

    #if last entry is corrupt
    if ( $this_file && !$this_file->{$name} ) {

        #warning ("Residual data at fileend.Ignoring\n");
        $this_file = ();
    }

    if ( scalar @array_events == 0 ) {
        throw("No events loaded from $filename. Cannot continue.");
    }

    printf "Files to process:%d \n", scalar @array_events if $verbose;

    close(INFILE);

    return \@array_events;
}

sub perldocs {
    exec( 'perldoc', $0 );
    exit(0);
}

=pod

=head1 NAME

ReseqTrack/scripts/event/load_event_from_conf.pl

=head1 SYNOPSIS

    This script takes a config file in the following format and parse it to create
    and load event objects into the event table

        [event_name]
        program=/path/to/program
        program_version=program version
        options=command line options for version
        input_flag=command line flag to associate with input
        farm_options=options for LSF submission of event
        batch_size=number of inputs to associate in a single LSF submission
        output_path=where to direct the stdout and stderr of the file
        type=the type of the input data required by the event
        table_name=the table name needed to retrieve the input data

    It can also dump existing events in this format


=head2 OPTIONS

        -dbhost, the name of the mysql-host
        -dbname, the name of the mysql database
        -dbuser, the name of the mysql user
        -dbpass, the database password if appropriate
        -dbport, the port the mysql instance is running on, this defaults to 4197 the 
         standard port for mysql-g1kdcc.ebi.ac.uk
        -file, name of file to read in or write out
        -read, binary flag to indicate the file specified by -file should be parsed and 
         loaded into the database
        -write, binary flag to indicate the file specified should be written in the given
          format based on the workflow objects already in the database
        -help, binary flag to indicate the help should be printed


=head1 Example:


$DB_OPTS="-dbhost mysql-host -dbuser rw_user -dbpass **** -dbport 4197 -dbname my_database"


    This dumps the existing event objects into a file of the given name:
    perl ReseqTrack/scripts/event/load_event_from_conf.pl  $DB_OPTS -file /path/to/event.conf -write
   
    
    This will load the event objects specifed in event.conf:
    perl ReseqTrack/scripts/event/load_event_from_conf.pl  $DB_OPTS -file /path/to/event.conf -read
    The write process will overwrite any existing file with the given name

=head2 Other useful scripts

    ReseqTrack/workflow/load_workflow_from_conf.pl is another script which works in a
    very similar manner but for workflow objects

=cut

