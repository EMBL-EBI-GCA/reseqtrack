#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::DBSQL::FileAdaptor;
use ReseqTrack::VerifyBamIDSample;
use ReseqTrack::VerifyBamIDReadGroup;
use ReseqTrack::Tools::RunVerifyBamID;
use ReseqTrack::Tools::GeneralUtils qw (get_params);
use ReseqTrack::Tools::BamUtils;
use ReseqTrack::Tools::FileSystemUtils qw (create_tmp_process_dir delete_directory);
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
my %input;
my %rg_info;


$input{echo_cmd_line} = 0;
$input{test}          = 0;

my $Sample;
my $passed_step1 = 0;
my $BAD          = 0;
my $selfSM;
my $selfRG;
my $bestRG;
my $bestSM;
my $got_sample_result;
my $got_readgroup_result;
my $Sample_adaptor;
my $RG_adaptor;
my $bam_file_obj;

GetOptions(
	   \%input,    'dbhost=s',   'dbname=s',      'dbuser=s',
	   'dbpass=s', 'dbport=s',   'working_dir=s', 'verbose!',
	   'bam=s',    'cfg_file=s', 'out_prefix=s', 'echo_cmd_line!',
	   'bimp=s', 'debug!', 'selfonly!', 'update!', 'selfSM=s',
	   'selfRG=s','bestRG=s', 'test!', 'save_files_for_deletion!');


if ( defined $input{cfg_file} ) {
  get_params( $input{cfg_file}, \%input );
}

#print Dumper %input;


my ( $db, $fa )  = get_db_adaptors( \%input );
$bam_file_obj = $fa->fetch_by_name( $input{bam} );

$Sample_adaptor = $db->get_VerifyBamIDSampleAdaptor();
$RG_adaptor     = $db->get_VerifyBamIDReadGroupAdaptor();



die "Could not retrieve " . $input{bam} . " from file table\n"
  if ( !$bam_file_obj );

print "File id = ", $bam_file_obj->dbID, "\n";




$got_sample_result    = $Sample_adaptor->fetch_by_other_id( $bam_file_obj->dbID );
$got_readgroup_result = $RG_adaptor->fetch_by_other_id( $bam_file_obj->dbID );


if ( $input{selfonly} && scalar (@$got_readgroup_result) ){
  print Dumper $got_readgroup_result;
 throw ("\nYou are running with 'selfonly' option, but there are read group\n" .
       "results associated with file_id = " .  $bam_file_obj->dbID . "\n" )
}


my $update = 0;
$update = $input{update} if (defined $input{update});


if ( $got_sample_result && ! $update) {
  print "Have results for " , $input{bam}, " update option = 0. Exitting\n";
  exit;
}

if ( $got_sample_result &&  $update) {
  print "Have results for " , $input{bam}, " update option = 1. Updating\n";
} 

if ( ! $got_sample_result) {
  print "No previous results found. Processing\n";
}


my ( $sample2, $platform2, $algorithm2, $project2, $analysis2, $chrom2, $date2 )
  = CHECK_AND_PARSE_FILE_NAME( $input{bam} );


$db->dbc->disconnect_when_inactive(2);

my $tmp_dir     = create_tmp_process_dir ($input{working_dir});

my $VOBJ = 'ReseqTrack::Tools::RunVerifyBamID';

my $VBAM = $VOBJ->new(
		      -program       => $input{program},
		      -reference     => $input{reference},
		      -input_files   => $input{bam},
		      -bfile         => $input{bfile},
		      -bimp          => $input{bimp},
		      -echo_cmd_line => $input{echo_cmd_line},
		      -out_prefix    => $input{out_prefix},
		      -options       => $input{options},
		      -selfonly      => $input{selfonly},
		      -working_dir   => $tmp_dir,
		     );


die "\$VBAM not created\n" if (!$VBAM);



$VBAM->run;



exit if $input{echo_cmd_line};


$selfSM = get_verifybamid_file($VBAM->selfSM_file, 0);
$selfRG = get_verifybamid_file($VBAM->selfRG_file, 1);
$bestRG = get_verifybamid_file($VBAM->bestRG_file, 1) if (defined $VBAM->bestRG_file);
$bestSM = get_verifybamid_file($VBAM->bestSM_file, 0) if (defined $VBAM->bestSM_file);

# Step 1. Look for SELFIBD < 0.98 and %MIX > 0.02
( $Sample, $passed_step1 ) = &create_Sample_object_step1($selfSM);


#Step 2. Check high selfSM_MIX against low %MIX in selfRG file.
#High selfSM_MIX  and low selfRG %MIX is bad.
$BAD += &check_selfRG_data( $Sample, $selfRG );


#Step 3. Check all SELFIBD is selfRG file are > 0.98
#print "Step 3  Check selfRG ------------\n";
$BAD += &check_SELFIBD_in_selfRG( $selfRG, $Sample->num_run_ids );



#Step 4  Check for sample swaps
$BAD += &check_for_sample_swaps($bestRG, $Sample->SEQ_SM) if ($bestRG) ;

print $input{bam}, "  is BAD ($BAD)\n" if ($BAD);
print "\n-----------------------------------------\n";



foreach my $key ( keys %rg_info ) {
  my $x =  $rg_info{$key};
 

  my $RG_OBJ = "ReseqTrack::VerifyBamIDReadGroup";
  my $rg_obj = $RG_OBJ->new(
			    -other_id => $bam_file_obj->dbID,
			    -run_id   => $key,
			    -SELFIBD  => $rg_info{$key}{'SELFIBD'},
			    -SELFMIX  => $rg_info{$key}{'SELFMIX'},
			    -BEST_SM  => $rg_info{$key}{'BEST_SM'},
			    -BESTIBD  => $rg_info{$key}{'BESTIBD'},
			    -BESTMIX  => $rg_info{$key}{'BESTMIX'},
			   );
	
  my $got_result = $RG_adaptor->fetch_by_other_id_and_run_id( $bam_file_obj->dbID, $key );
	
  if ($got_result) {
    $RG_adaptor->update($rg_obj)  unless ($input{test});
  } else {
    $RG_adaptor->store($rg_obj)   unless ($input{test}) ;
  }
	 
}




if ($BAD) {
  $Sample->failed(1);
}

if (! $passed_step1 && (! $bestRG || ! $bestSM )) {
  print "Only have selfonly results. Failed SELFIBD\/MIX criteria\n";
  $Sample->failed(1);
}


if ( $got_sample_result) {
  $Sample_adaptor->update($Sample) unless ($input{test});

} else {
  $Sample_adaptor->store($Sample) unless ($input{test}) ;
}

if ( -e $VBAM->{working_dir}){
  delete_directory ($VBAM->{working_dir}) unless $input{debug};
}


#=====================================================
#=====================================================
#=====================================================
sub check_SELFIBD_in_selfRG {
  my ( $selfRG, $totalRG ) = @_;
  my $ctr = 0;
  my $BAD = 0;
  print "In  check_SELFIBD_in_selfRG\n";

  foreach my $key ( keys %$selfRG ) {
    next if ( $key =~ /header/ );

    my $arr_ref = $$selfRG{$key};
    my $sibd    = $$arr_ref[ $$selfRG{header}{'SELFIBD'} ];
    my $mix     = $$arr_ref[ $$selfRG{header}{'%MIX'} ];

    if ( $sibd eq "N\/A" || $sibd < 0.98 ) {
      $ctr++;
      print
	"$key (SELFIBD threshold  < 0.98) SELFIBD $sibd \%mix $mix ($ctr\/$totalRG)\n";
      $BAD++;
    }
  }


  return $BAD;
}




sub check_for_sample_swaps {

  print "Step 4  Check for swaps ------------\n";

  my ($bestRG) = shift;
  my $BAD    = 0;

  my $contam = 0;
  foreach my $key ( keys %$bestRG ) {
    next if ( $key =~ /header/ );

    my $arr_ref = $$bestRG{$key};
    my $rg      = $$arr_ref[ $$selfRG{header}{'RG'} ];
    my $claimed = $$arr_ref[ $$bestRG{header}{'SEQ_SM'} ];
    my $best    = $$arr_ref[ $$bestRG{header}{'BEST_SM'} ];
    my $bibd    = $$arr_ref[ $$bestRG{header}{'BESTIBD'} ];




    $rg_info{$rg}{BESTMIX} = $$arr_ref[ $$bestRG{header}{'%MIX'} ];
    $rg_info{$rg}{BESTIBD} = $bibd;
    $rg_info{$rg}{BEST_SM} = $best;

   

    if ( ( $claimed ne $best ) && ( $bibd > 0.98 ) ) {
      print
	"$key possible sample swap. Claimed = $claimed Best = $best BESTIBD $bibd  \n";
      $contam++;
      $BAD++;
      next;
    }
  }

  return $BAD;
}

sub check_selfRG_data {
  my ( $Sample, $selfRG ) = @_;
  my $totalRG = 0;
  my $low_mix = 0;
  my $BAD     = 0;
  my $selfSM_MIX = 0.0;
   
  $totalRG = scalar( keys %$selfRG );
  $totalRG--;
     
  #if ( !$passed_step1 ) {

  print "Step 2 Check selfRG\n";


 # print "total RG = $totalRG\n";

  my $high_mix = 0;
  my $ctr      = 0;

  foreach my $key ( keys %$selfRG ) {

    next if ( $key =~ /header/ );
    my $arr_ref  = $$selfRG{$key};
    my $rg       = $$arr_ref[ $$selfRG{header}{'RG'} ];
    my $srg_sibd = $$arr_ref[ $$selfRG{header}{'SELFIBD'} ];
    my $srg_mix  = $$arr_ref[ $$selfRG{header}{'%MIX'} ];

    $rg_info{$rg}{'SELFIBD'} = $srg_sibd;
    $rg_info{$rg}{'SELFMIX'} = $srg_mix;


    if (! (defined  $rg_info{$rg}{BESTMIX})) {
      $rg_info{$rg}{BESTMIX} = 'N/A';
    }

    if (! (defined  $rg_info{$rg}{BESTIBD})) {
      $rg_info{$rg}{BESTIBD} = 'N/A';
    }

    if (! (defined  $rg_info{$rg}{BEST_SM})) {
      $rg_info{$rg}{BEST_SM} = 'N/A';
    }

  
    if ( $srg_mix <= 0.02 && $selfSM_MIX >= 0.02 ) {
      $low_mix++;
    }
  }

  my $poss_contam = 0.0;

  if ($low_mix) {
    $poss_contam = $low_mix / $totalRG * 100.0;
    $BAD++;

  }

  print
    "\% (high SM mix\/low RG mix)  $low_mix \/ $totalRG = $poss_contam\n";

  #	}

  $Sample->num_run_ids($totalRG);
  $Sample->num_low_selfIBD_run_ids($low_mix);


  return $BAD;
}

sub create_Sample_object_step1 {

  #Load info from selfSM file
  #if SELFIBD < 0.98 or %MIX > 0.02%. Possible problem

  $selfSM = shift;
  my $passed_step1 = 0;

  foreach my $key ( keys %$selfSM ) {
    next if ( $key =~ /header/ );
    my $a_ref = $$selfSM{$key};

    my $OBJ = "ReseqTrack::VerifyBamIDSample";
    $Sample = $OBJ->new(
			-table_name        => "file",
			-other_id          => $bam_file_obj->dbID,
			-SEQ_SM            => $$a_ref[ $$selfSM{header}{'SEQ_SM'} ],
			-SELFIBD           => $$a_ref[ $$selfSM{header}{'SELFIBD'} ],
			-SELFIBDLLK        => $$a_ref[ $$selfSM{header}{'SELFIBDLLK'} ],
			-SELFIBDLLKdiff    => $$a_ref[ $$selfSM{header}{'SELFIBDLLK-'} ],
			-HET_A1            => $$a_ref[ $$selfSM{header}{'HET-A1%'} ],
			-ALT_A1            => $$a_ref[ $$selfSM{header}{'ALT-A1%'} ],
			-DP                => $$a_ref[ $$selfSM{header}{'#DP'} ],
			-MIX               => $$a_ref[ $$selfSM{header}{'%MIX'} ],
			-HOM               => $$a_ref[ $$selfSM{header}{'%HOM'} ],
			-BESTHOMMIXLLK     => $$a_ref[ $$selfSM{header}{'BESTHOMMIXLLK'} ],
			-BESTHOMMIXLLKdiff => $$a_ref[ $$selfSM{header}{'BESTHOMMIXLLK-'} ],
		       );

    my $selfibd    = $$a_ref[ $$selfSM{header}{'SELFIBD'} ];
    my $selfSM_MIX = $$a_ref[ $$selfSM{header}{'%MIX'} ];

    print "Step 1:  :";
    print "selfSM $key : SELFIBD $selfibd \%MIX $selfSM_MIX.\n";

    if ( $selfibd eq "N/A") {
      print "Possible no snps for sample\n";
      return ( $Sample, $passed_step1 );	
    }   


    if ( $selfibd < 0.98 || $selfSM_MIX > 0.02 ) {
      print "Flag as 'POSSIBLE BAD\n";
      if ( $selfibd > ( 1.0 - $selfSM_MIX ) ) {
	print "Alt measure: $selfibd > (1.0 - $selfSM_MIX ) says OK\n";
      } else {
	print "\n";
      }

    } else {
      print "Pass Step 1:  :\n";
      $passed_step1 = 1;	# still have to do Step 3
    }

  }

  return ( $Sample, $passed_step1 );
}

sub get_verifybamid_file {
  my ( $file, $key_col ) = @_;
  my $ctr = 0;

  my @header;
  my %dat;

  print "Reading $file data\n";


  if (!-e $file) {
    print "$file does not exist.  --self_only option on ??\n";
    return;
  }

  open my $FH, '<', $file, or die "Could not open $file";

  while (<$FH>) {
   
    chomp;

    if (/^SEQ_SM/) {
      my @header = split /\t/;
      my $x      = 0;
      foreach my $h (@header) {
	$dat{header}{$h} = $x++;
      }
      next;
    }

    my @aa = split /\t/;
    $dat{ $aa[$key_col] } = \@aa;

  }



  close($FH);
  return ( \%dat );
}

sub assign_output_files {

  my $file_list = shift;

  my $selfSM = "";
  my $selfRG = "";
  my $bestRG = "";
  my $bestSM = "";

  foreach my $name (@$file_list) {
    
    print $name, "\n";
    if ( $name =~ /\.selfSM/ ) {
      $selfSM = $name;
      next;
    }

    if ( $name =~ /\.selfRG/ ) {
      $selfRG = $name;
      next;
    }
    if ( $name =~ /\.bestRG/ ) {
      $bestRG = $name;
      next;
    }
    if ( $name =~ /\.bestSM/ ) {
      $bestSM = $name;
      next;
    }

  }

  return ( $selfSM, $selfRG, $bestRG, $bestSM );
  exit;
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

  return ( $db, $fa );
}



=pod

=head1 NAME

ReseqTrack/scripts/qc/run_verifybam.pl

=head1 SYNOPSIS

 This script will take runs the program VerifyBamID on a bam file and stores
 a subsection of the results in the DB table VerifyBam_Sample and VerifyBam_ReadGroup.
 Different types of bam require different command lines. These are loaded via a cfg file.

=head1 OPTIONS

Database options

This is the database the objects will be written to

 -dbhost, the name of the mysql-host
 -dbname, the name of the mysql database
 -dbuser, the name of the mysql user
 -dbpass, the database password if appropriate
 -update,    Use if you wish to update results. If results are present in DB
            and -update is not used the script exit.
 -working_dir, Base directory for processing in temp directories. These 
           should be removed after processing is complete.

 -save_files_for_deletion, RunProgram.pm option that prevents deletion of output files.   

 VerifyBamID options.
 -selfonly,  only check sample against reference snps. Quicker.
 -bimp,      bed format file containing reference snps info.
            (bed file via something ~ ......
            vcftools --plink --vcf reference_snps.vcf
            rename out ref_snps out*
            plink --file ref_snps --maf 0.01 --geno 0.05 --make-bed
            rename plink ref_snps plink*)
 -reference, location of fasta formatted reference sequence. Also present
            must be ref.umfa file. VerifyBamID will create one if missing


 Note: If there are results for a run_id in the Sample and Read_group tables you cannot 
 run the script with the option combination  '-update -selfonly' as this will update the
 VerifyBamID_Sample table, but not the associated VerifyBamID_ReadGroup entries.


example

 perl reseqtrack/scripts/qc/run_verifybamid.pl -dbhost mysql-g1kdcc-public -dbport 4197 -dbuser g1krw -dbpass **** -dbname g1k_archive_staging_track -cfg_file /nfs/1000g-work/G1K/work/rseqpipe/perl_code/reseq-personal/rseqpipe/prog_conf/verifybam_chr20.cfg -working_dir /tmp -update -bam /nfs/1000g-archive/vol1/ftp/data/NA11831/exome_alignment/NA11831.mapped.SOLID.bfast.CEU.exome.20110411.bam


 the cfg file will look similar to 

 ftp_root=/nfs/1000g-archive/vol1/ftp/
 program=/nfs/1000g-work/G1K/work/bin/verifyBamID/bin/verifyBamID
 bimp=/nfs/1000g-work/G1K/work/REFERENCE/verifybam/CHROM20/chr20
 reference=/nfs/1000g-work/G1K/work/REFERENCE/verifybam/VerifyBamID/verifyBamID-0.0.5/reference/human_g1k_v37.fa
 working_dir=/nfs/1000g-work/G1K/scratch/rseqpipe/verifybamid/tmp
 options=-m 10 -g 5e-3
 selfonly=0

 Different options are required for different types of bams ( ie low coverage, exome ).


=cut
