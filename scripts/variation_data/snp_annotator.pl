#!/sw/bin/arch/perl -w

use strict;
use Getopt::Long;
use File::Basename;

use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::Annotation::SNPAnnotation;

my $input_file;
my $output_file;
my $working_dir;
my $vcftools;
my @source_vcfs;
my $check_ensembl_sources;
my $registry_file;
my $source;
my $strict_gvf;
my $type = 'SNV';
my $species;

$| = 1;

&GetOptions(
	    'input_file:s' => \$input_file,
	    'output_file:s' => \$output_file,
	    'working_dir:s' => \$working_dir,
	    'vcftools:s' => \$vcftools,
	    'source_vcfs:s@' => \@source_vcfs,
	    'check_ensembl_sources!' => \$check_ensembl_sources,
	    'registry_file:s' => \$registry_file,
	    'source:s' => \$source,
	    'strict_gvf!' => \$strict_gvf,
	    'type:s' => \$type,
	    'species:s' => \$species,
	   ) or throw("Failed to get options $!");


my $annotate = ReseqTrack::Tools::Annotation::SNPAnnotation->new
  (
   -input_file => $input_file,
   -output_file => $output_file,
   -working_dir => $working_dir,
   -vcftools => $vcftools,
   -source_vcf => \@source_vcfs,
   -check_ensembl => $check_ensembl_sources,
   -registry => $registry_file,
   -source => $source,
   -strict_gvf => $strict_gvf,
   -type => $type,
   -species => $species,
  );

$annotate->process_input();
$annotate->print_output();
