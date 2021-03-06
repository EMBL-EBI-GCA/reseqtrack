#!/usr/bin/env perl

use strict;
use warnings;
use ReseqTrack::FileTypeRule;
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
my $ftra = $db->get_FileTypeRuleAdaptor;

if ($read) {
    if ($clear) {
        $ftra->delete_all;
    }
    parse_file();
}

if ($write) {
    dump_ftrs();
}



sub parse_file {
    my $existing_ftrs = $ftra->fetch_all;
    my %existing_rule_orders;
    foreach my $ftr (@$existing_ftrs) {
        $existing_rule_orders{$ftr->rule_block_order}{$ftr->rule_order} = 1;
    }


    my $lines = get_lines_from_file($file);
    my %file_type_rules;
    foreach my $line (@$lines) {
        next $line if ($line =~ /^#/);
        my ($rule_block_order, $rule_order, $file_type, $match_regex) = split(/\s+/, $line, 4);

        throw("rule_block_order must be defined") if (!$rule_block_order);
        throw("rule_order must be defined") if (!$rule_order);
        throw("file_type_rule already exists in database with rule_block_order $rule_block_order and rule_order $rule_order")
            if ($existing_rule_orders{$rule_block_order}{$rule_order});
        throw("file_type_rule multiply defined with rule_block_order $rule_block_order and rule_order $rule_order in $file")
            if ($file_type_rules{$rule_block_order}{$rule_order});

        my $file_type_rule = ReseqTrack::FileTypeRule->new(
                                -rule_block_order => $rule_block_order,
                                -rule_order => $rule_order,
                                -file_type => $file_type,
                                -match_regex => $match_regex);
        $file_type_rules{$rule_block_order}{$rule_order} = $file_type_rule;
    }

    foreach my $rule_block_order (keys %file_type_rules) {
        foreach my $rule_order (keys %{$file_type_rules{$rule_block_order}}) {
            $ftra->store($file_type_rules{$rule_block_order}{$rule_order});
        }
    }
}

sub dump_ftrs {
    open my $OUT, '>', $file or die "cannot open $file $!";
    my $ftrs = $ftra->fetch_all_in_order;

    foreach my $block (@$ftrs) {
        foreach my $ftr (@$block) {
            print $OUT join("\t", $ftr->rule_block_order, $ftr->rule_order,
                                $ftr->file_type, $ftr->match_regex), "\n";
        }
    }
    close $OUT;
}




sub perldocs {
    exec('perldoc', $0);
    exit(0);
}

=pod

=head1 NAME

ReseqTrack/scripts/files/load_file_type_rules.pl

=head1 SYNOPSIS

    This script can perform either of the following functions:

    1) Get existing file_type_rule objects from the database and write them to a config file

    2) Take a config file and parse it to create and load file_type_rule objects into
    the file_type_rule table

    Suggested use for updating the file_type_rule table:
    First, use -write to write existing file_type_rules to a config file.
    Edit the config file to add / delete / modify file_type_rules
    Then, use -read -clear to clear the file_type_rule table and to add the new file_type_rules


=head2 CONFIG FILE FORMAT

    Each file_type_rule is written on a single line
    Lines beginning with '#' are ignored
    Each file_type_rule line contains four columns separated by whitespace

    Column 1:   integer, the number of the block to which the file_type_rule
                belongs. When assigning file types, each block is considered in order.
    Column 2:   integer, the order within the block.  When assigning file types, each
                file_type_rule in a block is tested in order until a match is successful.
    Column 3:   string, file type.  e.g. 'BAM'.  The character '*' is used to modify
                existing file types e.g. 'PILOT_*'
    Column 4:   string, regular expressions separated by brackets and keywords NOT, AND, OR.
                When assigning file types, the file path is tested against these regular expressions.

    e.g.
        # rule_block_order  rule_order  file_type  match_regex
        1   1   BAM     /\.bam$/
        1   2   FASTQ   /\.fastq\.gz$/ OR /\.fq.gz$/
        1   3   README  /README/ AND NOT /changelog/i
        1   4   OTHER   //

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

    This will load the file_type_rule objects specifed in file_type_rule.conf:
    perl ReseqTrack/scripts/event/load_file_type_rules.pl  $DB_OPTS -file /path/to/file_type_rule.conf -read

    This will dump the existing file_type_rule objects into current_file_type_rules.conf:
    perl ReseqTrack/scripts/event/load_file_type_rules.pl  $DB_OPTS -file /path/to/current_file_type_rules.conf -write

=cut

