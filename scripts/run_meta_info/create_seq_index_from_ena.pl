#! /usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;

use Getopt::Long;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::ERAUtils qw(get_erapro_conn get_fastq_details);
use ReseqTrack::Tools::SequenceIndexUtils qw(assign_files);

my @era_params;
my @studies_to_add;
my $analysis_group;
my $clob_read_length = 66000;
my $header;
my $file_format;
my $help;
my $ftp_prefix = "ftp://ftp.sra.ebi.ac.uk/";

&GetOptions(
    'era_dbuser=s'       => \$era_params[0],
    'era_dbpass=s'       => \$era_params[1],
    'era_dbname=s'       => \$era_params[2],
    'clob_read_length=i' => \$clob_read_length,
    'new_study=s'        => \@studies_to_add,
    'analysis_group=s'   => \$analysis_group,
    'header'             => \$header,
    'file_format=s'      => \$file_format,
    'help'               => \$help,
);

if ($help) {
    usage();
}

if ($header) {
    header();
}
if (!(($file_format eq "fastq") | ($file_format eq "bam") | ($file_format eq "pacbio"))) {
    die "Must specify -file_format (fastq, pacbio or bam) for files being indexed. fastq for standard submissions of fastq, pacbio for pacbio hdf5 and bam for bam only data sets such as 10X\n";
}

#era connection
my $era_db = get_erapro_conn(@era_params);
$era_db->dbc->db_handle->{LongReadLen} = $clob_read_length;

#get adaptors
my $study_adapt = $era_db->get_StudyAdaptor();
my $sample_adapt = $era_db->get_SampleAdaptor();
my $experiment_adapt = $era_db->get_ExperimentAdaptor();
my $run_adapt = $era_db->get_RunAdaptor();

#non-fastq file info sql
my $info_sql = 'select run_id, XMLType.getStringVal(run_xml), filename, checksum from run, xmltable(\'//FILE\' PASSING run_xml COLUMNS filename varchar2(512) PATH \'@filename\',  checksum varchar2(512) PATH \'@checksum\')(+)
where run_id = ?';
my $file_sth = $era_db->dbc->prepare($info_sql);

#work through study list...
for my $study_id (@studies_to_add) {

    my $study = $study_adapt->fetch_by_study_id($study_id);
    my %samples;
    my %experiments;

    my $runs = $run_adapt->fetch_by_study_id($study_id);

    #work through runs to create sequence index entries at run level
    RUN:
    for my $run (@$runs) {
        my $sample_pop;
        next RUN if ($run->status ne 'public');

        #need experiment and sample
        my $sample;
        my $run_experiment = $experiment_adapt->fetch_by_source_id($run->experiment_id);

        my $run_sample = $sample_adapt->fetch_by_source_id($run_experiment->sample_id);
        $sample_adapt->attach_attributes($run_sample);
        my $sample_attribs = $run_sample->attributes();

        for my $attrib (@$sample_attribs) {
            if ($attrib->attribute_name() eq "POPULATION") {
                $sample_pop = $attrib->attribute_value();
            }
        }


        #run stats available for fastq and bam but not pacbio
        my $run_stats = $run_adapt->get_run_stats($run->source_id());
        my ($md5hash, $sizehash, $namehash, @run_files, $mate1, $mate2, $frag);

        #fastq
        if ($file_format eq "fastq") {
            ($md5hash, $sizehash, $namehash) = get_fastq_details($run->source_id(), $era_db, $ftp_prefix);

            #assign files to work out pairing of files
            for my $filename (keys %$namehash) {
                push @run_files, $filename;
            }
            ($mate1, $mate2, $frag) = assign_files(\@run_files);
        }
        elsif (($file_format eq "bam") | ($file_format eq "pacbio")) {
            $file_sth->execute($run->source_id());
            while (my $hashref = $file_sth->fetchrow_hashref) {
                my $filename = $hashref->{FILENAME};
                my $md5 = $hashref->{CHECKSUM};
                $$namehash{$filename} = $ftp_prefix . "vol1/" . $filename;
                $$md5hash{$ftp_prefix . "vol1/" . $filename} = $md5;
            }
        }


        #for each file returned for the run
        for my $filename (keys %$namehash) {
            my $ftp_path = $$namehash{$filename};
            my $index_line = $ftp_path;
            $index_line = $index_line . "\t" . $$md5hash{$ftp_path};
            $index_line = $index_line . "\t" . $run->source_id();
            $index_line = $index_line . "\t" . $study_id;
            $index_line = $index_line . "\t" . @$study[0]->title();
            $index_line = $index_line . "\t" . $run->center_name();
            $index_line = $index_line . "\t" . $run->submission_id();
            $index_line = $index_line . "\t" . $run->submission_date();
            $index_line = $index_line . "\t" . $run_sample->source_id();
            $index_line = $index_line . "\t" . $run_sample->sample_alias();
            $index_line = $index_line . "\t" . $sample_pop;
            $index_line = $index_line . "\t" . $run_experiment->source_id();
            $index_line = $index_line . "\t" . $run_experiment->instrument_platform();
            $index_line = $index_line . "\t" . $run_experiment->instrument_model();
            $index_line = $index_line . "\t" . $run_experiment->library_name();
            $index_line = $index_line . "\t" . $run->run_alias();
            $index_line = $index_line . "\t" . $run_experiment->paired_nominal_length();
            $index_line = $index_line . "\t" . $run_experiment->library_layout();
            if ($mate1 && $filename eq $mate1) {
                $index_line = $index_line . "\t" . $$namehash{$mate2};
            }
            elsif ($mate2 && $filename eq $mate2) {
                $index_line = $index_line . "\t" . $$namehash{$mate1};
            }
            else {
                $index_line = $index_line . "\t";
            }
            if ($file_format eq "pacbio") {
                #no base and read count run stats for pacbio
                $index_line = $index_line . "\t";
                $index_line = $index_line . "\t"
            }
            else {
                $index_line = $index_line . "\t" . $$run_stats{'READ_COUNT'};
                $index_line = $index_line . "\t" . $$run_stats{'BASE_COUNT'};
            }
            $index_line = $index_line . "\t" . $analysis_group;
            print $index_line . "\n";
        }
    }
}

sub header {
    #create index file header
    #print date line and header in VCF like format
    ##dateline
    my ($day, $month, $year) = (localtime)[3, 4, 5];
    $month = sprintf '%02d', $month + 1;
    $day = sprintf '%02d', $day;
    $year = $year + 1900;
    my $fileDate = "##FileDate=" . $year . $month . $day . "\n";
    print $fileDate;

    print "##ENA_FILE_PATH=path to ENA file on ENA ftp site\n";
    print "##MD5=md5sum of file\n";
    print "##RUN_ID=SRA/ERA run accession\n";
    print "##STUDY_ID=SRA/ERA study accession\n";
    print "##STUDY_NAME=Name of study\n";
    print "##CENTER_NAME=Submission centre name\n";
    print "##SUBMISSION_ID=SRA/ERA submission accession\n";
    print "##SUBMISSION_DATE=Date sequence submitted, YYYY-MM-DD\n";
    print "##SAMPLE_ID=SRA/ERA sample accession\n";
    print "##SAMPLE_NAME=Sample name\n";
    print "##POPULATION=Sample population. Further information may be available with the data collection.\n";
    print "##EXPERIMENT_ID=Experiment accession\n";
    print "##INSTRUMENT_PLATFORM=Type of sequencing machine\n";
    print "##INSTRUMENT_MODEL=Model of sequencing machine\n";
    print "##LIBRARY_NAME=Library name\n";
    print "##RUN_NAME=Name of machine run\n";

    #print "##RUN_BLOCK_NAME=Name of machine run sector  (This is no longer recorded so this column is entirely null, it was left in so as not to disrupt existing sequence index parsers)\n";
    print "##INSERT_SIZE=Submitter specified insert size/paired nominal length\n";

    print "##LIBRARY_LAYOUT=Library layout, this can be either PAIRED or SINGLE\n";

    print "##PAIRED_FASTQ=Name of mate pair file if exists (Runs with failed mates will have a library layout of PAIRED but no paired fastq file)\n";

    #print "##WITHDRAWN=0/1 to indicate if the file has been withdrawn, only present if a file has been withdrawn\n";
    ##print "##WITHDRAWN_DATE=This is the date the index file is generated\n";
    ##print "##COMMENT=Comment about reason for withdrawal\n";

    print "##READ_COUNT=Read count for the file\n";
    print "##BASE_COUNT=Basepair count for the file\n";

    print "##ANALYSIS_GROUP=Analysis group is used to identify groups, or sets, of data. Further information may be available with the data collection.\n";
    print join("\t", ("#ENA_FILE_PATH", "MD5SUM", "RUN_ID", "STUDY_ID", "STUDY_NAME", "CENTER_NAME", "SUBMISSION_ID", "SUBMISSION_DATE", "SAMPLE_ID", "SAMPLE_NAME", "POPULATION", "EXPERIMENT_ID", "INSTRUMENT_PLATFORM", "INSTRUMENT_MODEL", "LIBRARY_NAME", "RUN_NAME", "INSERT_SIZE", "LIBRARY_LAYOUT", "PAIRED_FASTQ", "READ_COUNT", "BASE_COUNT", "ANALYSIS_GROUP")) . "\n";

}

sub usage {
    exec('perldoc', $0);
    exit(0);
}


1;

