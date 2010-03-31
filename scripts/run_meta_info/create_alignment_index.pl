#!/sw/arch/bin/perl 

use strict;
use ReseqTrack::Tools::Exception;
use ReseqTrack::Tools::FileUtils;
use File::Basename;
use ReseqTrack::DBSQL::DBAdaptor;

use Getopt::Long;

$| = 1;

my $file_list;
my $staging_root = "/nfs/1000g-archive/vol1/ftp/";
my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4175;
my $dbname;
my $type = 'PILOT_BAM';

&GetOptions(
    'dbhost=s'       => \$dbhost,
    'dbname=s'       => \$dbname,
    'dbuser=s'       => \$dbuser,
    'dbpass=s'       => \$dbpass,
    'dbport=s'       => \$dbport,
    'file_list=s'    => \$file_list,
    'staging_root=s' => \$staging_root,
    'type=s'         => \$type,
);

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
    -host   => $dbhost,
    -user   => $dbuser,
    -port   => $dbport,
    -dbname => $dbname,
    -pass   => $dbpass,
);

my $fa    = $db->get_FileAdaptor;
my @infiles;
my $bamfiles = $fa->fetch_by_type("PILOT_BAM");
push (@infiles, @$bamfiles);
my $baifiles = $fa->fetch_by_type("PILOT_BAI");
push (@infiles, @$baifiles);
my $basfiles = $fa->fetch_by_type("PILOT_BAS");
push (@infiles, @$basfiles);

my $files = \@infiles;

throw "No files of type $type" if ( !$files );

#need to split files into .bam and .bai
my %hash;
my @bams;
print STDERR "Have " . @$files . " files \n";
foreach my $file (@$files) {
    next unless ( $file->full_path =~ /1000g-archive/ );
    next if ( $file->name =~ /pu$/ );

    #next unless($file->full_path =~ /data/);
    my $stem = $file->name;
    $stem =~ s/\.bai$//;
    $stem =~ s/\.bas$//;
    $stem =~ s/\.bam$//;

    #print STDERR "Looking at ".$file->name."\n";
    $hash{$stem} = {} if ( !$hash{$stem} );
    if ( $file->name =~ /\.bai/ ) {
        $hash{$stem}->{bai} = $file;
    }
    elsif ( $file->name =~ /\.bas/ ) {

        #print STDERR "HAVE BAS FILE Stem is ".$stem."\n";
        $hash{$stem}->{bas} = $file;
    }
    else {
        $hash{$stem}->{bam} = $file;
        push( @bams, $file );
    }
}

#my $extra_hash = $hash{'NA19046.454.ssaha.SRP000033.2009.07'};
#foreach my $key(keys(%$extra_hash)){
#  print $key." NA19046.454.ssaha.SRP000033.2009.07 ".$extra_hash->{$key}."\n";
#}
my @sorted = sort { $a->name cmp $b->name } @$files;
foreach my $file (@bams) {

    #next unless($file->is_current);
    #data/NA19239/alignment/NA19239.chromX.SRP000032.2009_02.bam
    my $path = $file->full_path;
    $path =~ s/$staging_root//;
    my ( $name, $dir ) = fileparse($path);
    $name =~ /(NA\d+).+(SRP\d+).+bam/;
    my $individual = $1;
    my $study      = $2;
    my $md5        = $file->md5;
    my $stem       = $file->name;
    my $bai        = undef;
    my $bai_md5    = undef;
    my $bas        = undef;
    my $bas_md5    = undef;
    $stem =~ s/\.bam//;

    #print STDERR "STEM IS ".$stem."\n";
    if ( $hash{$stem}->{bai} ) {
        $bai = $hash{$stem}->{bai}->full_path;
        $bai =~ s/$staging_root//;
        $bai_md5 = $hash{$stem}->{bai}->md5;
    }
    else {
        unless ( $path =~ /unmapped/ ) {
            throw( "Failed to find a bai file for " . $path );
        }
    }
    if ( $hash{$stem}->{bas} ) {

        #print STDERR "Have BAS FILE for ".$stem."\n";
        $bas = $hash{$stem}->{bas}->full_path;
        $bas =~ s/$staging_root//;
        $bas_md5 = $hash{$stem}->{bas}->md5;
    }
    if ( !$bas || !$bas_md5 ) {

        #print STDERR "WARNING ".$stem." has no bas file\n";
    }
    print join( "\t",
        $path, $md5, $study, $individual, $bai, $bai_md5, $bas, $bas_md5 )
      . "\n";
}

