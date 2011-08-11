#adapted from Rick's script /homes/smithre/OneKGenomes/reseq-personal/smithre/may_2011/scripts/make_bas_scrip.pl

use strict;
use warnings;


use ReseqTrack::Tools::Bas;
use ReseqTrack::Tools::AlignmentBase;
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);

use File::Basename;
use ReseqTrack::Tools::BamUtils;
use Getopt::Long;
use Data::Dumper;

my $input_bam_name;
my $working_dir = '/nfs/nobackup/resequencing_informatics/rseqpipe/create_bas';
my $need_tags   = 0;
my $in_parent   = 0;
my $samtools;
my $reference;

&GetOptions(
  'bam=s'		=> \$input_bam_name,
  'working_dir=s'       => \$working_dir,
  'need_tags!'          => \$need_tags,
  'in_parent!'          => \$in_parent,
  'samtools=s'          => \$samtools,
  'reference=s'         => \$reference,
);

my $bam_basename = basename($input_bam_name);

my ($sample2, $platform2, $algorithm2, $project2, $analysis2, $chrom2, $date2) = CHECK_AND_PARSE_FILE_NAME($bam_basename);
print "release data is $date2\ninput file is $input_bam_name\n";

my $bas =  ReseqTrack::Tools::Bas->new (
                     -reference => $reference,
                     -bam => $input_bam_name,
                     -working_dir => $working_dir,
                     -need_tags=> $need_tags,
                     -release_date => $date2,
		     -in_parent => $in_parent,
		     -samtools => $samtools,
                    );
#print Dumper ($bas);

$bas->run;

#print Dumper ($bas);
#  -output_dir => "/homes/smithre/",
