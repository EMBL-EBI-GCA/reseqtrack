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

  	Ian Streeter (streeter@ebi.ac.uk)

=cut


use Getopt::Long;
use Net::FTP;
use Env qw( @PATH );

use strict;
use warnings;

my @populations;
my $panel;
my $ftp_host;
my $vcf;
my $region;
my $tabix = 'tabix';
my $ped;
my $info;
my $help;

GetOptions('population=s' => \@populations,
            'vcf=s' => \$vcf,
            'panel=s' => \$panel,
            'region=s' => \$region,
            'tabix=s' => \$tabix,
            'ped=s' => \$ped,
            'info=s' => \$info,
            'help!' => \$help,
            );

if ($help) {
    exec('perldoc', $0);
}

die("required arguments: vcf, panel, region, population") if (! $vcf || ! $panel || ! $region || ! @populations);
die("cannot find executable $tabix") if (! -x $tabix && ! grep {-x "$_/$tabix"} @PATH);

if (! $ped) {
    $ped = "$region.ped";
    $ped =~ s{:}{_};
}
if (! $info) {
    $info = "$region.info";
    $info =~ s{:}{_};
}


my $individuals = get_individuals($panel, \@populations);
my ($markers, $genotypes) = get_markers_genotypes($vcf, $region, $tabix, $individuals);

print_info($markers, $info);
print_ped($genotypes, $ped);








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
        INDIVIDUAL:
        foreach my $individual (@$individuals) {
            next INDIVIDUAL if (! $column_indices{$individual});
            my $genotype_string = $columns[ $column_indices{$individual} ];
            $genotype_string =~ /(\d+)(?:\/|\|)(\d+)/;
            my @genotype_codes = ($allele_codes[$1], $allele_codes[$2]);

            $alleles_present{$_} = 1 foreach (@genotype_codes);
            $marker_genotypes{$individual} = \@genotype_codes;
        }

        next LINE if ((scalar grep {$_} keys %alleles_present) < 2);

        foreach my $individual (keys %marker_genotypes) {
            push(@{$genotypes{$individual}}, $marker_genotypes{$individual});
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
    foreach my $individual (keys %$genotypes) {
        print $FILE join("\t", $individual, 1, 0, 0, 0, 0,);
        foreach my $genotype_codes (@{$genotypes->{$individual}}) {
            print $FILE "\t", $genotype_codes->[0], ' ', $genotype_codes->[1];
        }
        print $FILE "\n";
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
    my ($panel, $allowed_pops) = @_;

    my @panel_lines;

    if ($panel =~ /ftp:\/\/([\w.]+)(\/\S+)/) {
        my $ftp_host = $1;
        my $path = $2;

        my $ftp = Net::FTP->new($ftp_host);
        $ftp->login or die('Cannot login ' , $ftp->message);

        my $panel_content;
        open my $PANEL, '>', \$panel_content;
        $ftp->get($path, $PANEL) or die ('could not read panel ' , $ftp->message);
        $ftp->quit;
        close $PANEL;

        @panel_lines = split(/\n/, $panel_content);
    }
    else {
        open my $FILE, '<', $panel
            or die("cannot open $panel $!");
        @panel_lines = <$FILE>;
        close $FILE;
    }

    my %allowed_pops_hash;
    foreach my $pop (@$allowed_pops) {
        $allowed_pops_hash{$pop} = 1;
    }

    my @individuals;
    foreach my $line (@panel_lines) {
        my ($individual, $population) = split(/\s+/, $line);
        if ($allowed_pops_hash{$population}) {
            push(@individuals, $individual);
        }
    }
    return \@individuals;
}



##################################################################################################################################

=pod

=head1 NAME 

	vcf_to_haploview.pl

=head1 SYNOPSIS

        This script prepares the input files for Haploview	

=head1	VERSION

	1.0

=head1	REQUIRED ARGUMENTS

	-vcf		Path to a locally or remotely accessible tabix indexed vcf file.
                        The vcf file must be compressed by bgzip and indexed by tabix.
                        The vcf format is a tab format for presenting variation sites and 
			genotypes data and is described at http://vcftools.sourceforge.net/specs.html.
                        This tool takes both vcf4.0 and vcf4.1 format files.
	-panel		Path to a locally or remotely accessible panel file, listing all individuals (first column)
                        and their population (second column)
	-region		Chromosomal region in the format of chr:start-end (1:1000000-100500).
	-population     A population name, which must appear in the second column of the panel file.
                        Can be specified more than once for multiple populations.

=head1	OPTIONAL ARGUMENTS

	-tabix			Path to the tabix executable; default is to search the path for 'tabix'
	-ped			Name of the output ped file (linkage pedigree file);
                                default is region.ped (e.g. 1_100000-100500.ped)
        -info                   Name of the output info file (marker information file);
                                default is region.info (e.g. 1_1000000-100500.info)
	-help			Print out help menu
			
=head1	OUTPUT FILES

        The file formats of the linkage pedigree and marker information files are described at
        http://http://www.broadinstitute.org/science/programs/medical-and-population-genetics/haploview/input-file-formats-0
        Both files are needed as input for Haploview.

=head1 EXAMPLE

perl ~/ReseqTrack/scripts/variation_data/vcf_to_haploview.pl -vcf ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr13.phase1_integrated_calls.20101123.snps_indels_svs.genotypes.vcf.gz -panel ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/phase1_integrated_calls.20101123.ALL.panel -region 13:32889611-32973805 -population GBR -population FIN
