#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::DBSQL::FileAdaptor;
use ReseqTrack::DBSQL::VerifyBamIDAdaptor;  
use ReseqTrack::VerifyBamID;
use ReseqTrack::Tools::RunVerifyBamID;
use ReseqTrack::Tools::GeneralUtils qw (get_params get_open_file_handle );
use ReseqTrack::Tools::BamUtils;
use ReseqTrack::Tools::FileSystemUtils
  qw (create_tmp_process_dir delete_directory get_lines_from_file);
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::BamUtils qw (get_bam_run_id_info);
use ReseqTrack::Tools::QC::GLFUtils qw ( have_snps);
use File::Basename;

my %input;
$input{echo_cmd_line}           = 0;
$input{test}                    = 0;
$input{save_files_for_deletion} = 0;
my $bam_header_info;
my $VBAM;
my $BAM;
my ( $sample,$read_group);

GetOptions(
	\%input,         'dbhost=s',	'dbname=s',      'dbuser=s',
	'dbpass=s',      'dbport=s',	'working_dir=s', 'verbose!',
	'name=s',        'cfg_file=s',	'out_prefix=s',  'echo_cmd_line!',
	'debug!',        'update!',	'test!',         'save_files_for_deletion!',
	'snps_list=s',   'vcf_root:s', 'vcf_suffix:s', 'mapping=s', 'vcf:s', 'name=s',  'collection_type=s',
	'chrom=s',       'program=s',	'not_event!',    'options=s', 'run_mode=s',
	   'bam_name=s', 'do_read_groups!',
);

if ( $input{vcf_root} && (!$input{mapping} || $input{vcf_suffix}) ) {
	throw("Please provide a sample to pop mapping file and a vcf_suffix if you want to do by population verifybam with vcf_root $input{vcf_root}");
}
if ( $input{vcf_root} && $input{vcf} ) {
	throw("Want to run by population verifybam or not? Don't provide vcf_root and vcf at the same time");
}
	 
&get_params( $input{cfg_file}, \%input )if (defined $input{cfg_file});

$input{run_mode} = "self" if (!defined $input{run_mode});

my ( $db, $fa, $ca, $va ) = get_db_adaptors( \%input );

if ( defined $input{name} ) {    #pull bam name from collection
	$input{bam_name} =
	  get_bam_name_from_collection( $input{name},
		$input{collection_type} );
}
else {
	print "Using ", $input{bam_name},"\n";
	$input{bam_name} = $input{bam_name};
}

$BAM = $fa->fetch_by_name( $input{bam_name} );

# need file_id for table
$input{file_id} = $BAM->dbID;   

#BCM Solid bams might have RG info in odd spot
$bam_header_info = &get_bam_run_id_info( $input{bam_name});

&set_analysis_mode( \%input);

&determine_if_should_run( \%input ) if ( !defined $input{no_event} );

#&set_options( \%input ) if ( !defined $input{options} );
&set_options( \%input );

$VBAM = &create_run_program_object (\%input );

#print Dumper $VBAM;


$db->dbc->disconnect_when_inactive(1);

$$VBAM->run;

if ( $input{run_mode} eq "self"){
	$sample      =  read_output_file ($$VBAM->selfSM, $input{debug});
	$read_group  =  read_output_file ($$VBAM->selfRG, $input{debug});
}

if ( $input{run_mode} eq "best"){
	$sample     =  read_output_file ($$VBAM->bestSM, $input{debug});
	$read_group =  read_output_file ($$VBAM->bestRG, $input{debug});
}

###################


my $sample_self = @{$sample}[0]; # 1 line file evaluate bam as whole

throw ( "Do not have self file for this bam .Something failed\n") if ( !$sample_self);


$sample_self->verdict("HIGH_FREE_MIX") if ($sample_self->free_contam > $input{mix_cutoff} );
$sample_self->verdict("HIGH_CHIP_MIX") if ($sample_self->chip_contam > $input{mix_cutoff} );

if ( ($sample_self->chip_contam > 0.95) && $sample_self->free_contam < 0.05){
  print "Possible sample swap\n";
  $sample_self ->verdict("POSSIBLE_SWAP"); #overwrite verdict
}

if ($read_group){

  foreach my $rg (@$read_group){

    correct_run_id_pu_mismatch ( $rg, $bam_header_info);

    $rg->verdict("HIGH_FREE_MIX") if ( $rg->free_contam > $input{mix_cutoff} );

    #if no genotypes in vcf you will not get numeric $rg->chip_contam
    if ($rg->chip_contam ne "NA"){

      $rg->verdict("HIGH_CHIP_MIX") if ( $rg->chip_contam > $input{mix_cutoff} );
  

      #checking lanes for possible sample swaps
      if ( ($rg->chip_contam > 0.95) && $rg->free_contam < 0.05){
	print "Possible sample swap\n";
	$rg->verdict("POSSIBLE_SWAP"); #overwrite verdict
      }

    }


    # You will have to do something with this is moving to FREE only mode
    # my $large_diff = 10;
    # my $diff = abs($rg->free_mlogl_est_contam - $rg->free_mlogl_zero_contam);

    # if ( $rg->free_contam < $input{mix_cutoff} && $diff > $large_diff ){
    #   print "Possible sample contamination $diff\n";
    #   $rg->verdict("SAMPLE_CONTAM");
    # }


  }
}


my $have_result = $va->fetch_by_vcf_file_id_and_readgroup( $sample_self->vcf, $sample_self->file_id, $sample_self->read_group);  ##### FIXME, need to handle multiple entries

if ($have_result ) {
	### When file_id and read_group are the same, update if same vcf was used; create new record if different vcf was used
		if ( $have_result->vcf eq $sample_self->vcf ) {
			$va->update( $sample_self ) unless ($input{test});
		}
		else {
			$va->store ( $sample_self ) unless ($input{test});
		}	
}
else {	
	$va->store(  $sample_self ) unless ($input{test});
}

foreach my $rg (@$read_group){
  my $have_result = $va->fetch_by_vcf_file_id_and_readgroup( $rg->vcf, $rg->file_id, $rg->read_group);
#  $have_result ?  $va->update ($rg) :  $va->store($rg) unless ($input{test}); 
	if ($have_result) {
		### When file_id and read_group are the same, update if same vcf was used; create new record if different vcf was used
		if ( $have_result->vcf eq $rg->vcf ) {
			$va->update( $rg ) unless ($input{test});
		}
		else {
			$va->store ( $rg ) unless ($input{test});
		}	
	}
	else {	
		$va->store(  $rg ) unless ($input{test});
	}
}


print $$VBAM->working_dir,"\n";

#$$VBAM->files_to_delete($$VBAM->working_dir);


#==============================================================================
#===========================     SUBS  ========================================
#==============================================================================
sub which_vcf {
	my ($input) = @_;
	my %sample_pop_hash;
	my $lines = get_lines_from_file($$input{mapping});
	foreach my $line (@$lines ) {
		my ($pop, $sample) = split (/\t/, $line);
		$sample_pop_hash{$sample} = $pop;
	}
	
	my $vcf;
	if ($$input{vcf_root}) {
		my ($sample) = split (/\./,  basename($$input{bam_name}));
		my $pop = $sample_pop_hash{$sample};
		if (defined $pop ) {
			my $possible_vcf = $$input{vcf_root} . "/$pop." . $$input{vcf_suffix};
			if ( -e $possible_vcf ) {
				$vcf = $possible_vcf;
			}
			else {
				warn("VCF file for population $pop does not exist, use VCF for ALL populations.");
				$vcf = $$input{vcf_root} . "/ALL." . $$input{vcf_suffix};
			}		
			
		}	
		else {
			throw("sample $sample doesn't mapped to any population!");
		}		
	}
  else {
      $vcf = $$input{vcf};
  }           
	return $vcf; 
}  
	
sub create_run_program_object {
  my ($input) = shift;

  my $tmp_dir = create_tmp_process_dir( $input{working_dir} );
	
  my %option_hash;
  if ($$input{options} ) { ### options are passed in as a string; this bit is to convert it to a hash
      my @bits = split (/--/, $$input{options} );
      foreach my $bit (@bits) {
          $bit =~ s/^\s+|\s+$//g;
          my $flag = "--" . $bit;
          next if ($flag eq "--");

          if ($flag !~ /\s+/) {
              $option_hash{$flag} = 1;
          }
          else {
          	my @bits2 = split (/\s+/, $flag);
          	$option_hash{$bits2[0]} = $bits2[1];
          }     
      }    
  }

	print "options hash is: ";
	foreach my $k ( keys%option_hash ) {
		print "$k\t$option_hash{$k}\n";
	}		
	
      
  my $vcf = which_vcf($input);
  
  my $PROG_TYPE = 'ReseqTrack::Tools::RunVerifyBamID';
  my $PROG = $PROG_TYPE->new(
			     -program                 => $$input{program},
			     -vcf                     => $vcf,
			     -input_files             => $$input{bam_name},
			     -echo_cmd_line           => $$input{echo_cmd_line},
			     -out_prefix              => $$input{out_prefix},
			     -options		     	=> \%option_hash,		
			     -working_dir             => $tmp_dir,
			     -save_files_for_deletion => $$input{save_files_for_deletion},
			     -debug                   => $$input{debug},
			     -run_mode                => $$input{run_mode},
			    );

  die "\$PROG not created\n" if ( !$PROG );

  print Dumper $PROG if (defined $$input{debug});

  return \$PROG ;
}


sub convert_to_verifybamid_object{
  my ($hash,$read_group) = @_;

  my $OBJ = 'ReseqTrack::VerifyBamID';

  my $vcf = which_vcf(\%input);

  my $VBID = $OBJ->new (

			-file_id                   => $input{file_id},
			-sample                    => $$hash{$read_group}{'#SEQ_ID'},
			-read_group                => $$hash{$read_group}{RG},
			-chip_id                   => $$hash{$read_group}{CHIP_ID},
			-snps                      => $$hash{$read_group}{"#SNPS"},
			-num_reads                 => $$hash{$read_group}->{"#READS"},
			-avg_depth                 => $$hash{$read_group}{AVG_DP},
			-free_contam               => $$hash{$read_group}->{FREEMIX},
			-free_mlogl_est_contam     => $$hash{$read_group}->{FREELK1},
			-free_mlogl_zero_contam    => $$hash{$read_group}->{FREELK0},
			-free_ref_bias_ref_het     => $$hash{$read_group}->{FREE_RH},
			-free_ref_bias_refhomalt   => $$hash{$read_group}->{FREE_RA},
			-chip_contam               => $$hash{$read_group}->{CHIPMIX},
			-chip_mlogl_est_contam     => $$hash{$read_group}->{CHIPLK0},
			-chip_mlogl_zero_contam    => $$hash{$read_group}->{CHIPLK1},
			-chip_ref_bias_ref_het     => $$hash{$read_group}->{CHIP_RH},
			-chip_ref_bias_refhomalt   => $$hash{$read_group}->{CHIP_RA},
			-depth_homref_site         => $$hash{$read_group}->{DPREF},
			-rel_depth_het_site        => $$hash{$read_group}->{RDPHET},
			-rel_depth_homalt_site     => $$hash{$read_group}->{RDPALT},
			-run_mode                  => $input{run_mode},
			-used_genotypes            => $input{used_genotypes},
			-target_region             => $input{chrom},
			-vcf                       => basename($vcf),			  
			-verdict                   => "PASS",

		       );
			#	-vcf                       => basename ($input{vcf}),
  print Dumper $VBID if ($input{debug});

  return $VBID;

}

sub read_output_file {

	my ($file,$debug) = @_;

	my @headers;
	my @objs;

	if (! -e $file){
	  print "$file does not exist.Skipping\n";
	  return;
	}
	my $IN = get_open_file_handle($file);
	while (<$IN>) {
		chomp;
		my %hash;

		if (/^\#SEQ/) {
			@headers = split /\t/;
			next;
		}

		my @values = split /\t/;
		foreach my $i ( 0 ... $#headers ) {
			$hash{ $values[1] }{ $headers[$i] } = $values[$i];
		}
		
		push ( @objs , convert_to_verifybamid_object (\%hash, $values[1]) );	
	}

	close $IN;

	return \@objs;
}


sub set_options{
	my ( $input) = @_;	

	if ( defined ($input{run_mode}) ){
	  print "Run mode specified = ", $input{run_mode},"\n";
	}
	else{
	  $input{run_mode} = "self"; 
	}

	#default to low coverage setting
	if ( $input{bam_name} =~ /exome/i )
	  {
	    $input{options} .= " --maxDepth 255 --precise "; # use to be 1000 maxDepth
	  }	
	 print "options =", $input{options},"\n" if  (defined $input{debug});
	 
	
	if ( ! defined $input{do_read_groups}) {
		$input{options} .= " --free-full --ignoreRG ";
	}
	else {
		$input{options} .= " --free-full ";
	}
	
	if (	$input{analysis_type} ne "FREE"	) {
		$input{options} .= " --chip-full ";
	}
	
		
		
	

	return;
}

sub get_bam_name_from_collection {

	  my ( $input, $coll_name, $type ) = @_;
	  print "Pulling bam name from COLLECTION:", $input{name}, "\n";

	  throw "No collection type given (option = -type)\n"
		unless ( defined $input{collection_type} );

	  throw
		"No 'chrom' or 'mapped'  given (eg -chrom chrom20 \/ -chrom mapped)\n"
		unless ( defined $input{chrom} );

	  my $coll =
		$ca->fetch_by_name_and_type( $input{name},
		  $input{collection_type} );
	  my $msg =
		  "No collection info found for "
		. $input{name} . " "
		. $input{collection_type} . "\n";

	  throw $msg if ( !$coll );

	  my $others = $coll->others;
	  my @all_coll_bams;
	  my $use_bam;
	  foreach my $file (@$others) {
		  push( @all_coll_bams, $file->name );
		  my $basename = basename( $file->name );
		  my @bits = split( /\./, $basename );
		  if ( $bits[1] eq $input{chrom} ) {
			  $use_bam = $file->name;
		  }
	  }

	  if ( !$use_bam ) {
		  my $msg = "Could not find bam fitting chrom " . $input{chrom};
		  $msg .= " in collection bams:\n";
		  my $line = join( "\n", @all_coll_bams );
		  $msg .= $line;
		  throw "$msg";

	  }
	  print "\nRunning: $use_bam\n";

	  return $use_bam;
}

sub determine_if_should_run {

	  my ($input) = @_;

	  #where is my bam?
	  if ( !-e $input{bam_name} ) {
		  my $msg = $input{bam_name} . " does not appear to exist";
		  warn $msg;
	  }

	  #'chrom=mapped', assuming properly name bam.
	  if ( defined $input{mapped}
		  && !( $input{bam_name} =~ /\.mapped\./i ) )
	  {
		  my $msg = "Flagged to run on genome-wide set of snps.\n";
		  $msg .= "Cannot find 'mapped' in file name\n";
		  $msg .= $input{bam_name} . "\n";
		  warning "$msg";
		  exit;
	  }

	  #do not run on unmapped bams
	  if ( $input{bam_name} =~ /unmapped/i ) {
		  my $msg = $input{bam_name} . "\n";
		  $msg .= "\nYou cannot run on what appears to be an unmapped bam\n";
		  $msg .= "Skipping this bam\n";
		  warning "$msg";
		  exit;
	  }

	  #do not run on any chrom11 bams
	  if ( $input{bam_name} =~ /chrom11/i ) {
		  my $msg = $input{bam_name} . "\n";
		  $msg .= "\nRunning on what appears to be an chrom11 bam\n";
		  $msg .=
"This bam probably has an overlapping type with bams that should be tested\n";
		  $msg .= "Skipping this bam\n";
		  warning "$msg";
		  exit;
	  }

	  #do not run on exome chrom20 bams
	  if (   ( $input{bam_name} =~ /chrom20/i )
		  && ( $input{bam_name} =~ /exome/i ) )
	  {
		  my $msg = "\nRunning on what appears to be an exome chrom20 bam\n";
		  $msg .=
"This bam probably has an overlapping type with bams that should be tested\n";
		  $msg .= "Skipping this bam\n";
		  warning "$msg";
		  exit;
	  }

	  my $ftp_root='/nfs/1000g-archive/vol1/ftp/';

	  if ( !$input{bam_name} =~ /$ftp_root/ ) {
		  my $msg .= $input{name};
		  $msg .= "\nis not on $ftp_root . Not running\n";
		  throw "$msg";
	  }

	  #print "Bam looks testable\n";
	  return;
}

sub set_analysis_mode {

	  my ( $input ) = @_;

	  my %hash;
	  my $cmd;
	  my ( $sample, $platform, $algorithm, $project, $analysis, $chrom, $date )
		= CHECK_AND_PARSE_FILE_NAME($$input{bam_name});

	  $$input{sample} = $sample;

#	  my $vcf = $$input{vcf};
	  my $vcf = which_vcf($input);
	  my $tbi = $vcf . '.tbi';

	  if ( -e $tbi ) {
		  $cmd = "tabix -H $vcf | grep \"^#CHROM\"";
	  }
	  elsif ( $vcf =~ /gz$/ ) {
		  $cmd = "zcat $vcf | head -300 | grep \"^#CHROM\"";
	  }
	  elsif ( $vcf =~ /vcf$/ ) {
		  $cmd = "cat $vcf | head -300 | grep \"^#CHROM\"";
	  }

	  print "$cmd\n"  if ($input{debug});
	  my $CHROM_line = `$cmd`;

	  map { $hash{$_}++ } ( split /\t/, $CHROM_line );

	  #you will have to play with this if you move to a sites only vcf
	  if ( !defined $hash{$sample} ) {
		  print "No genotypes for $sample. Using seq-only FREE analysis\n";
		  $$input{analysis_type} = "FREE";
          $$input{options} .= " --chip-none "; #turn off chip_mix calculation
		  $$input{mix_cutoff}    = 0.03;
		  $$input{used_genotypes} = 0;
	  }
	  else {
		  print "Have genotypes for $sample using CHIP analysis\n";
		  $$input{analysis_type}= "FREE"; # this is correct, but leave --chip_mix on
		  $$input{used_genotypes} = 1;
		  $$input{mix_cutoff}    = 0.02;
	  }
 
       return;
}

sub correct_run_id_pu_mismatch {

	  my ( $rg_obj, $bam_header_info ) = @_;

	  my $rg = $rg_obj->read_group;

	  if ( !( $rg =~ /^[E|S]RR/ ) ) {

	    if ( defined $$bam_header_info{ $rg_obj->read_group } ) {

	      my $tmp_rg = $$bam_header_info{ $rg_obj->read_group}{PU};

	      if ( $tmp_rg =~ /^[ERR|SRR]/ ) {
		print "Replacing ", $rg, " with $tmp_rg \n";
		$rg_obj->read_group($tmp_rg);
			  }
	      else {
		my $msg = "Baffled by: $tmp_rg (not SRR ERR). Ignoring\n";
		warn "$msg";
	      }
	    }
	  }

	  return;
}

sub get_db_adaptors {
	  my $input = shift;

	  my $db = ReseqTrack::DBSQL::DBAdaptor->new(
		  -host   => $$input{dbhost},
		  -user   => $$input{dbuser},
		  -port   => $$input{dbport},
		  -dbname => $$input{dbname},
		  -pass   => $$input{dbpass},
	  );

	  my $fa = $db->get_FileAdaptor;
	  my $ca = $db->get_CollectionAdaptor();
	  my $va = $db->get_VerifyBamIDAdaptor();

	  return ( $db, $fa, $ca,$va );
}

=pod

=head1 NAME

  run_20120626_verifybamid.pl

=head1 SYNOPSIS

  The Umich VerifybamID program to QC bam files

 More info
 http://www.ebi.ac.uk/seqdb/confluence/display/1000GEN/BAM+QC+Using+VerifyBamID
 
 for useful information.

  example command line:
  perl run_20120626_verifybamid.pl -chrom mapped -cfg_file $reseq-personal/rseqpipe/prog_conf/verifybamid_20120620_low_coverage_whole_genome.cfg -bam_name $FTP_SITE/data/HG00867/alignment/HG00867.mapped.SOLID.bfast.CDX.low_coverage.20111114.bam

 '-test', will run script but not store results

  for running on a single file
  
  -cgf_file contents can be overwritten on the command line 
  dbhost=mysql-g1kdcc-public
  dbuser=g1krw
  dbport=4197
  dbname=g1k_archive_staging_track
  dbpass=xxxx
  vcf=/nfs/1000g-work/G1K/work/REFERENCE/20120626_verifybam_vcfs/WHOLE_GENOME/omni_wg_2141_76000_filtered.vcf.gz
  working_dir=/nfs/1000g-work/G1K/scratch/rseqpipe/verifybamid/tmp/LOW_COVERAGE_WG
  program=/nfs/1000g-work/G1K/work/bin/verifybamid_20120620/verifyBamID/bin/verifyBamID
  chrom=mapped
  collection_type=BAM

 bam must be on ftp site to run as an event
 ( You can override all checks ( chrom 11, unmapped bams) using '-no_event' option

=cut


=pod

perl /nfs/1000g-work/G1K/work/zheng/reseqtrack/scripts/qc/run_verifyBam.phase3.pl -chrom mapped -cfg_file /nfs/1000g-work/G1K/work/rseqpipe/perl_code/reseq-personal/rseqpipe/prog_conf/verifybamid_20120620_exome_mapped.test3.cfg -name HG00377.ILLUMINA.bwa.exome > vbam.log &

Config files are
/nfs/1000g-work/G1K/work/rseqpipe/perl_code/reseq-personal/rseqpipe/prog_conf/verifybamid_by_pop.20130305_exome.test.cfg
/nfs/1000g-work/G1K/work/rseqpipe/perl_code/reseq-personal/rseqpipe/prog_conf/verifybamid_by_pop.20130305_wgs.test.cfg