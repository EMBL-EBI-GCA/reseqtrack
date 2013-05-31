use strict;
use ReseqTrack::Tools::Exception;
use ReseqTrack::Tools::UpdateRunMetaInfo;
use ReseqTrack::Tools::ERAUtils;
use Getopt::Long;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4175;
my $dbname;
my $era_dbuser;
my $era_dbpass;
my @modules;

