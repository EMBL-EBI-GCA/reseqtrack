package ReseqTrack::VerifyBamID;

use strict;
use warnings;
use vars qw(@ISA);
use ReseqTrack::DBSQL::BaseAdaptor;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename;
use Data::Dumper;

@ISA = qw(ReseqTrack::Base);

sub new {
    my ( $class, @args ) = @_;

    my $self = $class->SUPER::new(@args);

    my (
        $file_id,                 $sample,
        $read_group,              $chip_id,
        $snps,                    $num_reads,
        $avg_depth,               $free_contam,
        $free_mlogl_est_contam,   $free_mlogl_zero_contam,
        $free_ref_bias_ref_het,   $free_ref_bias_refhomalt,
        $chip_contam,             $chip_mlogl_est_contam,
        $chip_mlogl_zero_contam,  $chip_ref_bias_ref_het,
        $chip_ref_bias_refhomalt, $depth_homref_site,
        $rel_depth_het_site,      $rel_depth_homalt_site,
        $run_mode,                $used_genotypes,
        $target_region,           $vcf,
        $verdict,                 $performed,
      )
      = rearrange(
        [
            qw(
              FILE_ID
              SAMPLE
              READ_GROUP
              CHIP_ID
              SNPS
              NUM_READS
              AVG_DEPTH
              FREE_CONTAM
              FREE_MLOGL_EST_CONTAM
              FREE_MLOGL_ZERO_CONTAM
              FREE_REF_BIAS_REF_HET
              FREE_REF_BIAS_REFHOMALT
              CHIP_CONTAM
              CHIP_MLOGL_EST_CONTAM
              CHIP_MLOGL_ZERO_CONTAM
              CHIP_REF_BIAS_REF_HET
              CHIP_REF_BIAS_REFHOMALT
              DEPTH_HOMREF_SITE
              REL_DEPTH_HET_SITE
              REL_DEPTH_HOMALT_SITE
              RUN_MODE
              USED_GENOTYPES
              TARGET_REGION
              VCF
              VERDICT
              PERFORMED
              )
        ],
        @args
      );

    $self->file_id($file_id);
    $self->sample($sample);
    $self->read_group($read_group);
    $self->chip_id($chip_id);
    $self->snps($snps);
    $self->num_reads($num_reads);
    $self->avg_depth($avg_depth);
    $self->free_contam($free_contam);
    $self->free_mlogl_est_contam($free_mlogl_est_contam);
    $self->free_mlogl_zero_contam($free_mlogl_zero_contam);
    $self->free_ref_bias_ref_het($free_ref_bias_ref_het);
    $self->free_ref_bias_refhomalt($free_ref_bias_refhomalt);
    $self->chip_contam($chip_contam);
    $self->chip_mlogl_est_contam($chip_mlogl_est_contam);
    $self->chip_mlogl_zero_contam($chip_mlogl_zero_contam);
    $self->chip_ref_bias_ref_het($chip_ref_bias_ref_het);
    $self->chip_ref_bias_refhomalt($chip_ref_bias_refhomalt);
    $self->depth_homref_site($depth_homref_site);
    $self->rel_depth_het_site($rel_depth_het_site);
    $self->rel_depth_homalt_site($rel_depth_homalt_site);
    $self->run_mode($run_mode);
    $self->used_genotypes($used_genotypes);
    $self->target_region($target_region);
    $self->vcf($vcf);
    $self->verdict($verdict);
    $self->performed($performed);

    return $self;
}

sub file_id {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{file_id} = $arg;
    }
    return $self->{file_id}

}

sub sample {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{sample} = $arg;
    }
    return $self->{sample}

}

sub read_group {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{read_group} = $arg;
    }
    return $self->{read_group}

}

sub chip_id {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{chip_id} = $arg;
    }
    return $self->{chip_id}

}

sub snps {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{snps} = $arg;
    }
    return $self->{snps}

}

sub num_reads {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{num_reads} = $arg;
    }
    return $self->{num_reads}

}

sub avg_depth {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{avg_depth} = $arg;
    }
    return $self->{avg_depth}

}

sub free_contam {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{free_contam} = $arg;
    }
    return $self->{free_contam}

}

sub free_mlogl_est_contam {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{free_mlogl_est_contam} = $arg;
    }
    return $self->{free_mlogl_est_contam}

}

sub free_mlogl_zero_contam {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{free_mlogl_zero_contam} = $arg;
    }
    return $self->{free_mlogl_zero_contam}

}

sub free_ref_bias_ref_het {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{free_ref_bias_ref_het} = $arg;
    }
    return $self->{free_ref_bias_ref_het}

}

sub free_ref_bias_refhomalt {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{free_ref_bias_refhomalt} = $arg;
    }
    return $self->{free_ref_bias_refhomalt}

}

sub chip_contam {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{chip_contam} = $arg;
    }
    return $self->{chip_contam}

}

sub chip_mlogl_est_contam {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{chip_mlogl_est_contam} = $arg;
    }
    return $self->{chip_mlogl_est_contam}

}

sub chip_mlogl_zero_contam {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{chip_mlogl_zero_contam} = $arg;
    }
    return $self->{chip_mlogl_zero_contam}

}

sub chip_ref_bias_ref_het {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{chip_ref_bias_ref_het} = $arg;
    }
    return $self->{chip_ref_bias_ref_het}

}

sub chip_ref_bias_refhomalt {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{chip_ref_bias_refhomalt} = $arg;
    }
    return $self->{chip_ref_bias_refhomalt}

}

sub depth_homref_site {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{depth_homref_site} = $arg;
    }
    return $self->{depth_homref_site}

}

sub rel_depth_het_site {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{rel_depth_het_site} = $arg;
    }
    return $self->{rel_depth_het_site}

}

sub rel_depth_homalt_site {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{rel_depth_homalt_site} = $arg;
    }
    return $self->{rel_depth_homalt_site}

}

sub run_mode {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{run_mode} = $arg;
    }
    return $self->{run_mode}

}

sub used_genotypes {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{used_genotypes} = $arg;
    }
    return $self->{used_genotypes}

}

sub target_region {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{target_region} = $arg;
    }
    return $self->{target_region}

}

sub vcf {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{vcf} = $arg;
    }
    return $self->{vcf}

}

sub verdict {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{verdict} = $arg;
    }
    return $self->{verdict}

}

sub performed {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{performed} = $arg;
    }
    return $self->{performed}

}

1;

