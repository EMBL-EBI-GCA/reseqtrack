#!/sw/bin/arch/perl -w

=pod

=head1 NAME

reseqtrack/scripts/variation_data/snp_annotator.pl

=head1 SYNOPSIS

This script will take in a set of snps and use the ensembl variation database 
it is pointed at to calculate snp consequences. It can also give information about
which of the given snps can be found in the source in the Ensembl database or if
any of the snps are found in other supplied vcf files.

=head1 OPTIONS

-input_file, This is the input file to be processed. The file can be in 3 different
formats, VCF, pileup and a tab locational format with the 
columns, chr, start, end, allele_string (N/N), strand (-1 || 1) and name, if no name
is defined then a name is generated using the chromosome, coordinates and 
allele string, all the files can be raw or gzipped

-output_file, This is where the output is written in gvf format (http://www.sequenceontology.org/gvf.html) with the sources placed in additional columns after the 9
which belong to gvf if requested

-working_dir, This is where the tmp files produced by the process are written

-registry_file, This is the registry file needed to setup connected to the ensembl
databases. It should be in the format

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;


Bio::EnsEMBL::DBSQL::DBAdaptor->new
  ( '-species' => "Homo_sapiens",
    '-group'   => "core",
    '-port'    => 5306,
    '-host'    => 'ensembldb.ensembl.org',
    '-user'    => 'anonymous',
    '-pass'    => '',
    '-dbname'  => 'homo_sapiens_core_54_36p',
  );

Bio::EnsEMBL::Variation::DBSQL::DBAdaptor->new
  ( '-species' => "Homo_sapiens",
    '-group'   => "variation",
    '-port'    => 5306,
    '-host'    => 'ensembldb.ensembl.org',
    '-user'    => 'anonymous',
    '-pass'    => '',
    '-dbname'  => 'homo_sapiens_variation_54_36p', );


Bio::EnsEMBL::Registry->add_alias("Homo_sapiens","human");

If you don't want to supply a registry file you can also supply database arguments
and then the script uses the registry method load_registry_from_db

-dbhost, the host

-dbport, the port number

-dbuser, the username

-dbpass, the password

-db_version, the version of ensembl you want to query e.g 54

-species, the species to use to get databases from the registry file

-type, the type if variant calls, this defaults to SNV

-source, the source of the given variant calls, this is used in one of the gvf
columns, this defaults to the input filename

-check_ensembl_sources, a binary flag to indicate that you want any sources for 
a particular variant location from the ensembl database to be indicated

-vcftools, the location of your vcftools install

-source_vcfs, paths to vcf files your snps should be identified in using the vcftools
positions filter. These must be raw vcf files and cannot be compressed

-strict_gvf, this means that only the first 9 columns will be printed regardless of
other options specified

-help, This makes the script print out its options

=head1 Examples

This will identify what sources in ensembl/a vcf files the given snps occur in

perl reseqtrack/scripts/variation_data/snp_annotator.pl 
-input_file JPT.SRP000033.2010_03.sites.vcf.gz -registry ensembl.registry 
-source Pilot3 -output JPT.SRP000033.2010_03.snp_annotator -check_ensembl 
-source_vcf CEU.SRP000031.2010_03.genotypes.vcf -vcftools bin/vcftools/vcftools  
-working_dir /tmp

=cut

use strict;
use Getopt::Long;
use File::Basename;

use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::Annotation::SNPAnnotation;
use ReseqTrack::Tools::GeneralUtils qw(useage);

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
my $help;
my $dbhost;
my $dbpass;
my $dbport;
my $dbuser;
my $ensembl_version;


$| = 1;

&GetOptions(
	    'input_file:s' => \$input_file,
	    'output_file:s' => \$output_file,
	    'working_dir:s' => \$working_dir,
	    'vcftools:s' => \$vcftools,
	    'source_vcfs:s@' => \@source_vcfs,
	    'check_ensembl_sources!' => \$check_ensembl_sources,
	    'registry_file:s' => \$registry_file,
	    'dbhost:s' => \$dbhost,
	    'dbuser:s' => \$dbuser,
	    'dbpass:s' => \$dbpass,
	    'dbport:s' => \$dbport,
	    'ensembl_version:s' => \$ensembl_version,
	    'source:s' => \$source,
	    'strict_gvf!' => \$strict_gvf,
	    'type:s' => \$type,
	    'species:s' => \$species,
	    'help!' => \$help,
	   ) or throw("Failed to get options $!");

if($help){
  useage();
}

my $registry = 'Bio::EnsEMBL::Registry';
if($registry_file){
  $registry->load_all($registry_file);
}else{
  if($dbhost && $dbuser){
    $registry->load_registry_from_db( -host => $dbhost,
				      -user => $dbuser,
				      -pass => $dbpass,
				      -port => $dbport,
				      -db_version => $ensembl_version
				    );
  }else{
    throw("Failed to provide sufficient info to make a database connection");
  }
}


#create object
my $annotate = ReseqTrack::Tools::Annotation::SNPAnnotation->new
  (
   -input_file => $input_file,
   -output_file => $output_file,
   -working_dir => $working_dir,
   -vcftools => $vcftools,
   -source_vcf => \@source_vcfs,
   -check_ensembl => $check_ensembl_sources,
   -registry => $registry,
   -source => $source,
   -strict_gvf => $strict_gvf,
   -type => $type,
   -species => $species,
  );

#parse file
$annotate->process_input();
#calculate sources/consequences
$annotate->produce_output();
#print gvf out
$annotate->print_output();
#cleanup any tmp files
$annotate->cleanup_output();
