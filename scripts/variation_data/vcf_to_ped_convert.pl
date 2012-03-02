#!/usr/bin/env perl

=head1 LICENSE

	Copyright (c) 1999-2011 The European Bioinformatics Institute and
	Genome Research Limited.  All rights reserved.

	This software is distributed under a modified Apache license.
	For license details, please see

	http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

	Please email comments or questions to the 1000 genomes project information 
	list at <info@1000genomes.org>

=head1	AUTHOR

  	Ian Streeter (info@1000genomes.org)

=cut


use Getopt::Long;
use Net::FTP;
use Env qw( @PATH );

use strict;
use warnings;

my @populations;
my $sample_panel;
my $ftp_host;
my $vcf;
my $region;
my $tabix = 'tabix';
my $output_ped;
my $output_info;
my $output_dir;
my $help;

GetOptions('population=s' => \@populations,
            'vcf=s' => \$vcf,
            'sample_panel_file=s' => \$sample_panel,
            'region=s' => \$region,
            'tabix=s' => \$tabix,
            'output_ped=s' => \$output_ped,
            'output_info=s' => \$output_info,
            'output_dir=s' => \$output_dir,
            'help!' => \$help,
            );

if ($help) {
    exec('perldoc', $0);
}

die("required arguments: vcf, sample_panel_file, region, population") if (! $vcf || ! $sample_panel || ! $region || ! @populations);
die("$output_dir is not a directory") if ($output_dir && ! -d $output_dir);

my $is_compressed = $vcf =~ /\.g?gz(ip)?$/;
die("cannot find executable $tabix") if ($is_compressed && ! -x $tabix && ! grep {-x $_} map {"$_/$tabix"} @PATH);
die("remote vcf file must be compressed by bgzip") if (!$is_compressed && $vcf =~ /ftp:\/\//);

my ($region_chromosome, $region_start, $region_end);
if ($region =~ /^(\w+):(\d+)-(\d+)$/) {
  ($region_chromosome, $region_start, $region_end) = ($1, $2, $3);
}
else {
  die("did not recognise region $region");
}

if (! $output_ped) {
    $output_ped = "$region.ped";
    $output_ped =~ s{:}{_};
}
if (! $output_info) {
    $output_info = "$region.info";
    $output_info =~ s{:}{_};
}
if ($output_dir) {
    $output_ped = $output_dir . '/' . $output_ped;
    $output_ped =~ s{//}{/}g;

    $output_info = $output_dir . '/' . $output_info;
    $output_info =~ s{//}{/}g;
}

my %individuals;
my @markers;
my %genotypes;

get_individuals();
get_markers_genotypes();

print_info();
print_ped();

print "Created ".$output_info." and ".$output_ped."\n";



sub get_markers_genotypes {
    my %base_codes = ('A' => 1, 'C' => 2, 'G' => 3, 'T' => 4);

    my $vcf_opener = $is_compressed ? "$tabix -h $vcf $region |" : "<$vcf";
    open my $VCF, $vcf_opener
        or die("cannot open vcf $!");

    my %column_indices;

    LINE:
    while (my $line = <$VCF>) {
        next LINE if ($line =~ /^\#\#/);
        chomp $line;
        my @columns  = split(/\t/, $line);

        if ($line =~ /^\#/) {
            foreach my $i (0..$#columns) {
                $column_indices{$columns[$i]} = $i;
            }
            next LINE;
        }

        my ($chromosome, $position, $name, $ref_allele, $alt_alleles) = @columns;
        next LINE if ($chromosome ne $region_chromosome);
        next LINE if ($position < $region_start);
        last LINE if ($position > $region_end);

        my @allele_codes = map {$base_codes{$_} || 0} $ref_allele, (split(/,/, $alt_alleles));
        next LINE if ((scalar grep {$_} @allele_codes) < 2);

        my %marker_genotypes;
        my %alleles_present;
        foreach my $population (keys %individuals) {
            INDIVIDUAL:
            foreach my $individual (@{$individuals{$population}}) {
                next INDIVIDUAL if (! $column_indices{$individual});
                my $genotype_string = $columns[ $column_indices{$individual} ];
                if ($genotype_string =~ /(\d+)(?:\/|\|)(\d+)/) {
                  my @genotype_codes = ($allele_codes[$1], $allele_codes[$2]);

                  $alleles_present{$_} = 1 foreach (grep {$_} @genotype_codes);
                  $marker_genotypes{$population}{$individual} = \@genotype_codes;
                }
                else {
                  $marker_genotypes{$population}{$individual} = [0,0];
                }
            }
        }

        next LINE if ((scalar grep {$_} keys %alleles_present) < 2);

        foreach my $population (keys %marker_genotypes) {
            foreach my $individual (keys %{$marker_genotypes{$population}}) {
                push(@{$genotypes{$population}{$individual}}, $marker_genotypes{$population}{$individual});
            }
        }

        if ($name eq '.') {
            $name = "$chromosome:$position";
        }
        push(@markers, [$name,$position]);

    }

    close $VCF;

    if ($is_compressed) {
      my $exit_status = $? >>8;
      die("tabix exited with status $exit_status") if $exit_status;
    }

    return;
}

sub print_ped {

    open my $FILE, '>', $output_ped
        or die "cannot open $output_ped $!";
    foreach my $population (keys %genotypes) {
        my $pedigree_counter = 1;
        foreach my $individual (keys %{$genotypes{$population}}) {
            my $pedigree = $population . '_' . $pedigree_counter;
            print $FILE join("\t", $pedigree, $individual, 0, 0, 0, 0,);
            foreach my $genotype_codes (@{$genotypes{$population}->{$individual}}) {
                print $FILE "\t", $genotype_codes->[0], ' ', $genotype_codes->[1];
            }
            print $FILE "\n";
            $pedigree_counter ++;
        }
    }
    close $FILE;
    return;
}

sub print_info {

    open my $FILE, '>', $output_info
        or die "cannot open $output_info $!";
    foreach my $marker (@markers) {
        print $FILE join("\t", @$marker), "\n";
    }
    close $FILE;
    return;
}




sub get_individuals {

    my @sample_panel_lines;

    if ($sample_panel =~ /ftp:\/\/([\w.]+)(\/\S+)/) {
        my $ftp_host = $1;
        my $path = $2;

        my $ftp = Net::FTP->new($ftp_host);
        $ftp->login or die('Cannot login ' , $ftp->message);

        my $sample_panel_content;
        open my $PANEL, '>', \$sample_panel_content;
        $ftp->get($path, $PANEL) or die ('could not $sample_panel ' , $ftp->message);
        $ftp->quit;
        close $PANEL;

        @sample_panel_lines = split(/\n/, $sample_panel_content);
    }
    else {
        open my $FILE, '<', $sample_panel
            or die("cannot open $sample_panel $!");
        @sample_panel_lines = <$FILE>;
        close $FILE;
    }

    my %allowed_pops_hash;
    foreach my $pop (@populations) {
        $allowed_pops_hash{$pop} = 1;
        $individuals{$pop} = [];
    }

    foreach my $line (@sample_panel_lines) {
        my ($individual, $population) = split(/\s+/, $line);
        if ($allowed_pops_hash{$population}) {
            push(@{$individuals{$population}}, $individual);
        }
    }
}



##################################################################################################################################

=pod

=head1 NAME 

	vcf_to_ped_converter.pl

=head1 SYNOPSIS

        This script prepares the input files for Haploview	

=head1	VERSION

	1.0

=head1	REQUIRED ARGUMENTS

	-vcf		    Path to a locally or remotely accessible vcf file.
                            The vcf file must be compressed by bgzip and indexed by tabix if it is a remote file.
                            The vcf format is a tab format for presenting variation sites and 
			    genotypes data and is described at http://vcftools.sourceforge.net/specs.html.
                            This tool takes both vcf4.0 and vcf4.1 format files.
	-sample_panel_file  Path to a locally or remotely accessible sample panel file, listing all individuals (first column)
                            and their population (second column)
	-region		    Chromosomal region in the format of chr:start-end (1:1000000-100500).
	-population         A population name, which must appear in the second column of the sample panel file.
                            Can be specified more than once for multiple populations.

=head1	OPTIONAL ARGUMENTS

	-tabix		    Path to the tabix executable; default is to search the path for 'tabix'
	-output_ped	    Name of the output ped file (linkage pedigree file);
                            default is region.ped (e.g. 1_100000-100500.ped)
        -output_info        Name of the output info file (marker information file);
                            default is region.info (e.g. 1_1000000-100500.info)
        -output_dir         Name of a directory to place the output_ped and output_info files
	-help		    Print out help menu
			
=head1	OUTPUT FILES

        The file formats of the linkage pedigree and marker information files are described at
        http://http://www.broadinstitute.org/science/programs/medical-and-population-genetics/haploview/input-file-formats-0
        Both files are needed as input for Haploview.

=head1 EXAMPLE

perl ~/ReseqTrack/scripts/variation_data/vcf_to_ped_converter.pl -vcf ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr13.phase1_integrated_calls.20101123.snps_indels_svs.genotypes.vcf.gz -sample_panel_file ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/phase1_integrated_calls.20101123.ALL.sample_panel -region 13:32889611-32973805 -population GBR -population FIN
