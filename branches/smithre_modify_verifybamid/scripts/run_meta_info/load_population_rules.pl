#!/sw/bin/perl -w

use strict;
use ReseqTrack::PopulationRule;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception;
use ReseqTrack::Tools::FileSystemUtils qw( get_lines_from_file );
use Getopt::Long;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport;
my $dbname;

my $read;
my $write;
my $clear;
my $file;
my $help;

&GetOptions(
    'dbhost=s' => \$dbhost,
    'dbname=s' => \$dbname,
    'dbuser=s' => \$dbuser,
    'dbpass=s' => \$dbpass,
    'dbport=s' => \$dbport,
    'read!'     => \$read,
    'write!'    => \$write,
    'clear!'    => \$clear,
    'file=s'   => \$file,
    'help!'    => \$help,
);

if ($help) {
    perldocs();
}

throw("Must specify -read or -write") if (!$read && !$write);
throw("Must specify only one of -read or -write") if ($read && $write);
throw("Must specify -file to read from or write to") if (!$file);

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
    -host   => $dbhost,
    -user   => $dbuser,
    -port   => $dbport,
    -dbname => $dbname,
    -pass   => $dbpass,
);

my $pra = $db->get_PopulationRuleAdaptor;

if ($read) {
    if ($clear) {
        $pra->delete_all;
    }
    parse_file();
}

if ($write) {
    dump_population_rules();
}





sub parse_file {
    my $existing_prs = $pra->fetch_all;
    my %existing_rule_orders;
    foreach my $pr (@$existing_prs) {
        $existing_rule_orders{$pr->rule_order} = 1;
    }


    my $lines = get_lines_from_file($file);
    my %population_rules;
    foreach my $line (@$lines) {
        next $line if ($line =~ /^#/);
        my ($rule_order, $population, $match_regex) = split(/\s+/, $line, 3);

        throw("rule_order must be defined") if (!$rule_order);
        throw("population_rule already exists in database with rule_order $rule_order")
            if ($existing_rule_orders{$rule_order});
        throw("population_rule multiply defined with rule_order $rule_order in $file")
            if ($population_rules{$rule_order});

        my $population_rule = ReseqTrack::PopulationRule->new(
                                -rule_order => $rule_order,
                                -population => $population,
                                -match_regex => $match_regex);
        $population_rules{$rule_order} = $population_rule;
    }

    foreach my $rule_order (keys %population_rules) {
        $pra->store($population_rules{$rule_order});
    }
}

sub dump_population_rules {
    open my $OUT, '>', $file or die "cannot open $file $!";
    my $population_rules = $pra->fetch_all_in_order;

    foreach my $pr (@$population_rules) {
        print $OUT join("\t", $pr->rule_order, $pr->population, $pr->match_regex), "\n";
    }
    close $OUT;
}




sub perldocs {
    exec('perldoc', $0);
    exit(0);
}

=pod

=head1 NAME

ReseqTrack/scripts/files/load_population_rules.pl

=head1 SYNOPSIS

    This script can perform either of the following functions:

    1) Get existing population_rule objects from the database and write them to a config file

    2) Take a config file and parse it to create and load population_rule objects into
    the population_rule table

    Suggested use for updating the population_rule table:
    First, use -write to write existing population_rules to a config file.
    Edit the config file to add / delete / modify population_rules
    Then, use -read -clear to clear the population_rule table and to add the new population_rules


=head2 CONFIG FILE FORMAT

    Each population_rule is written on a single line
    Lines beginning with '#' are ignored
    Each population_rule line contains three columns separated by whitespace

    Column 1:   integer, the order in which the rule is considered. When assigning populations,
                each population_rule is tested in order until a match is successful.
    Column 2:   string, population code e.g. 'GBR'
    Column 3:   string, regular expressions separated by brackets and keywords NOT, AND, OR.
                When assigning populations, a descriptive string is tested against these regular expressions.

    e.g.
        # rule_order  population  match_regex
        1   YRI     /YRI/i
        2   CHB     /CHB/ OR /han chinese/
        3   TSI     /t[uo]scan/i
        4   OTHER   //

=head2 OPTIONS

        -dbhost, the name of the mysql-host
        -dbname, the name of the mysql database
        -dbuser, the name of the mysql user
        -dbpass, the database password if appropriate
        -dbport, the port the mysql instance is running on
        -file, name of file to read from or write to
        -read, flag to read file and insert into database
        -write, flag to write database file_type_rules to the file
        -clear, flag to delete all existing file_type_rules from the database. Only used with the -read flag.
        -help, flag to print this help and exit


=head1 Example:


$DB_OPTS="-dbhost mysql-host -dbuser rw_user -dbpass **** -dbport 4197 -dbname my_database"

    This will load the population_rule objects specifed in population_rule.conf:
    perl ReseqTrack/scripts/run_meta_info/load_population_rules.pl  $DB_OPTS -file /path/to/population_rule.conf -read

    This will dump the existing population_rule objects into current_population_rules.conf:
    perl ReseqTrack/scripts/run_meta_info/load_population_rules.pl  $DB_OPTS -file /path/to/current_population_rules.conf -write

=cut

