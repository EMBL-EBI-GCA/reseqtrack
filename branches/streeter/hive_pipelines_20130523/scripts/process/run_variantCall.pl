#!/usr/bin/env perl

use strict;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileUtils qw(create_object_from_path);
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists);
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use ReseqTrack::Tools::GeneralUtils qw(get_params);
use Getopt::Long;
use ReseqTrack::Tools::Loader::File;

$| = 1;

my @allowed_algorithms = qw(samtools umake gatk freebayes);

my %input;
GetOptions(
        \%input,    
        'dbhost=s',   
        'dbname=s',      
        'dbuser=s',
        'dbpass=s', 
        'dbport=s',  
        'cfg_file=s',
        'output_dir=s',
        'reference=s',
        'collection_name=s',
        'input_type=s',
        'host=s',
        'output_file_type=s',
        'algorithm=s',
        'program=s',
        'options=s%',
        'chrom=s',
        'region_start=s',
        'region_end=s',
        'store!',
        'update!',

        # specific for gatk:
        'java_exe=s',
        'jvm_args=s',
        'gatk_path=s',

        # specific for samtools:
        'bcftools=s',
        'vcfutils=s',
        'bgzip=s',

        );

if ( defined $input{cfg_file} ) {
  get_params( $input{cfg_file}, \%input );
}


throw("Must specify an output directory") if (!$input{output_dir});
throw("Must specify a reference") if (!$input{reference});
throw("Must specify a collection name") if (!$input{collection_name});
throw("Must specify an input_type") if (!$input{input_type});
$input{update} //= 0;
$input{host} ||= '1000genomes.ebi.ac.uk';
$input{output_file_type} ||= "VCF";

if ( !$input{algorithm} || ! grep {$input{algorithm} eq $_} @allowed_algorithms) {
  throw("must specify an algorithm from " . join(", ", @allowed_algorithms));
}
my $caller_module = load_caller_module($input{algorithm});

my @allowed_options = keys %{eval '&'."$caller_module".'::DEFAULT_OPTIONS'};
foreach my $option (keys %{$input{options}}) {
  throw("Don't recognise option $option. Acceptable options are: @allowed_options")
    if (! grep {$option eq $_ } @allowed_options);
}

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $input{dbhost},
  -user   => $input{dbuser},
  -port   => $input{dbport},
  -dbname => $input{dbname},
  -pass   => $input{dbpass},
    );

my $ca = $db->get_CollectionAdaptor;
my $collection = $ca->fetch_by_name_and_type($input{collection_name}, $input{input_type});
throw("Failed to find a collection for ".$input{collection_name}." ".$input{input_type}." from ".$input{dbname}) 
    if(!$collection);

my $bam_objects = $collection->others;
throw "no bams in collection" if !@$bam_objects;
my @bam_names;
foreach my $bam (@$bam_objects) {
  my $bam_name = $bam->name;
  throw("not a bam $bam_name") if $bam_name !~ /\.bam$/;
  throw("bam does not exist") if ! -f $bam_name;
  throw("bai does not exist") if ! -f "$bam_name.bai";
  push(@bam_names, $bam_name);
}


$db->dbc->disconnect_when_inactive(1);

my $caller_object = $caller_module->new(
                  -input_files             => \@bam_names,
                  -working_dir             => $input{output_dir},
                  -reference               => $input{reference},
                  -program                 => $input{program},
                  -job_name                => $input{collection_name},
                  -options                 => $input{options},
                  -chrom                   => $input{chrom},
                  -region_start            => $input{region_start},
                  -region_end              => $input{region_end},
                  );
if ($input{algorithm} eq 'gatk') {
  $caller_object->java_exe($input{java_exe});
  $caller_object->jvm_args($input{jvm_args});
  $caller_object->gatk_path($input{gatk_path});
}
if ($input{algorithm} eq 'samtools') {
  $caller_object->bcftools($input{bcftools});
  $caller_object->vcfutils($input{vcfutils});
  $caller_object->bgzip($input{bgzip});
}

$caller_object->run();
$db->dbc->disconnect_when_inactive(0);

if ($input{store}) {
  FILE:
  foreach my $outfile (@{$caller_object->output_files}) {
    check_file_exists($outfile);
    throw("unexpected output file $outfile") if $outfile !~ /\.vcf(?:\.gz)?$/;
    my $loader = ReseqTrack::Tools::Loader::File->new(
        -file => [$outfile],
        -do_md5 => 1,
        -hostname => $input{host},
        -db => $db,
        -assign_types => 0,
        -check_types => 0,
        -type => $input{output_file_type},
        -update_existing => $input{update},
    );
    $loader->process_input();
    $loader->create_objects();
    $loader->sanity_check_objects();
    $loader->load_objects();
  }
}

sub load_caller_module {
  my $algorithm_name = shift;
  my $module_name = 'CallBy' . ucfirst($algorithm_name);
  my $file = "ReseqTrack/Tools/RunVariantCall/$module_name.pm";
  eval {
    require "$file";
  };
  if ($@) {
    throw("cannot load $file: $@")
  }
  my $module = "ReseqTrack::Tools::RunVariantCall::$module_name";
  return $module;
}
    
                
sub help_info {

 exec('perldoc', $0);

}


=head1 SYNOPSIS

This is a generic script to call variants from BAM files. Algorthms that can be used for making varaint call include: 
samtools (http://samtools.sourceforge.net/mpileup.shtml)   
gatk (http://www.broadinstitute.org/gsa/wiki/index.php/Unified_genotyper)
umake (http://genome.sph.umich.edu/wiki/UMAKE)
freebayes (https://github.com/ekg/freebayes)

The input bam files are taken from a collection in the database. The output vcf will be written to the database.

=head1 OPTIONS

  database options:

    -dbhost, the name of the mysql-host
    -dbname, the name of the mysql database
    -dbuser, the name of the mysql user
    -dbpass, the database password if appropriate
    -dbport, the port the mysql instance is running on

  other options:

  -algorithm, must be one of 'samtools', 'gatk', 'freebayes', 'umake'

  -collection_name, name of a collection containing the input bam files
  -input_type, type of the collection containing the input bam files
  -output_file_type, file type when storing output vcf file in the database

  -output_dir, directory for the vcf file

  -reference, path to the reference fasta file

  -program, the program executable required by the caller module

  -options, for constructing the options hash passed to the caller module object
  e.g. -options dcov=100

  -chrom, string, only call variants on this chromosome
  -region_start, integer, specifies a region for calling variants
  -region_end, integer, specifies a region for calling variants

  -update, boolean, for loading input files into the database
  -store, boolean, for loading input files into the database

  -host, default is '1000genomes.ebi.ac.uk', needed for storing output vcf file

  some options specific for the gatk algorithm:
  -java_exe, path to the java executable. The GATKTools class can guess at its location if it is not specified.
  -jvm_args, options of java. The GATKTools class uses default values if nothing is specified.
  -gatk_path, path to a directory containing the gatk jar files.
  Alternatively, the GATKTools class will use the $GATK environment variable if set.

  some options specific for the samtools algorithm:
  -bcftools, path to executable
  -vcfutils, path to perl script
  -bgzip, path to executable (or default will be used)


=head1 Examples

  $DB_OPTS="-dbhost mysql-host -dbuser rw_user -dbpass **** -dbport 4197 -dbname my_database"

  perl reseqtrack/process/run_variantCall.pl $DB_OPTS 
    -algorithm samtools -output_dir /path/to/dir -reference /path/to/ref.fa
    -collection_name GBR -input_type BAM
    -chrom 10 -region_start 500000 -region_end 600000
    -store
    -options mpileup='-EDS -e20 -h100 -L250 -o40 -C50 -m1 -F0.002 -d 2500 -P ILLUMINA -ug'

=cut
