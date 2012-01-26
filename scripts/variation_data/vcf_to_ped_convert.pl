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
die("cannot find executable $tabix") if (! -x $tabix && ! grep {-x "$_/$tabix"} @PATH);
die("$output_dir is not a directory") if ($output_dir && ! -d $output_dir);

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


my $individuals = get_individuals($sample_panel, \@populations);
my ($markers, $genotypes) = get_markers_genotypes($vcf, $region, $tabix, $individuals);

print_info($markers, $output_info);
print_ped($genotypes, $output_ped);








sub get_markers_genotypes {
    my ($vcf, $region, $tabix, $individuals) = @_;

    my %base_codes = ('A' => 1, 'C' => 2, 'G' => 3, 'T' => 4);

    my @markers;
    my %genotypes;

    open my $VCF, "$tabix -h $vcf $region |"
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

        my @allele_codes = map {$base_codes{$_} || 0} $ref_allele, (split(/,/, $alt_alleles));
        next LINE if ((scalar grep {$_} @allele_codes) < 2);

        my %marker_genotypes;
        my %alleles_present;
        foreach my $population (keys %$individuals) {
            INDIVIDUAL:
            foreach my $individual (@{$individuals->{$population}}) {
                next INDIVIDUAL if (! $column_indices{$individual});
                my $genotype_string = $columns[ $column_indices{$individual} ];
                $genotype_string =~ /(\d+)(?:\/|\|)(\d+)/;
                my @genotype_codes = ($allele_codes[$1], $allele_codes[$2]);

                $alleles_present{$_} = 1 foreach (@genotype_codes);
                $marker_genotypes{$population}{$individual} = \@genotype_codes;
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
    return \@markers, \%genotypes;
}

sub print_ped {
    my ($genotypes, $file) = @_;

    open my $FILE, '>', $file
        or die "cannot open $file $!";
    foreach my $population (keys %$genotypes) {
        my $pedigree_counter = 1;
        foreach my $individual (keys %{$genotypes->{$population}}) {
            my $pedigree = $population . '_' . $pedigree_counter;
            foreach my $genotype_codes (@{$genotypes->{$individual}}) {
                print $FILE "\t", $genotype_codes->[0], ' ', $genotype_codes->[1];
            }
            print $FILE "\n";
            $pedigree_counter ++;
        }
    }
    close $FILE;
}

sub print_info {
    my ($markers, $file) = @_;

    open my $FILE, '>', $file
        or die "cannot open $file $!";
    foreach my $marker (@$markers) {
        print $FILE join("\t", @$marker), "\n";
    }
    close $FILE;
    return;
}




sub get_individuals {
    my ($sample_panel, $allowed_pops) = @_;

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
    my %individuals;
    foreach my $pop (@$allowed_pops) {
        $allowed_pops_hash{$pop} = 1;
        $individuals{$pop} = [];
    }

    foreach my $line (@sample_panel_lines) {
        my ($individual, $population) = split(/\s+/, $line);
        if ($allowed_pops_hash{$population}) {
            push(@{$individuals{$population}}, $individual);
        }
    }
    return \%individuals;
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

	-vcf		    Path to a locally or remotely accessible tabix indexed vcf file.
                            The vcf file must be compressed by bgzip and indexed by tabix.
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
