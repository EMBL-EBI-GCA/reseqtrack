
=pod

=head1 NAME

ReseqTrack::Tools::QC::PPQT

=head1 SYNOPSIS

This is a class for running PhantomPeakQualTools to assess the 
quality of a chip seq experiment.

https://code.google.com/p/phantompeakqualtools/

This module only concerns itself with the qc functionality of ppqt,
not the IDR or peak calling functions

=head1 Example

my $ppqt = ReseqTrack::Tools::QC::PPQT->new(
  -program => '/path/to/ppqt/runspp.r',
  -rscript_path => '/path/to/Rscript'
  -input_files => ['my_bam_file.bam'],
  -job_name => 'bob'
);

my $metrics = $ppqt->run;
my $metrics_file = $ppqt->output_files->[0];

=cut

package ReseqTrack::Tools::QC::PPQT;

use strict;
use warnings;
use File::Basename;

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);

use ReseqTrack::Tools::GeneralUtils
  qw(execute_system_command execute_pipe_system_command);

use ReseqTrack::Tools::FileSystemUtils
  qw(check_file_exists check_file_does_not_exist);

use base qw(ReseqTrack::Tools::RunProgram);

sub new {
    my ( $class, @args ) = @_;

    my $self = $class->SUPER::new(@args);

    my (
        $rscript_path, $keep_metrics, $samtools_path,
        $keep_plot,    $keep_rdata,   $no_dups
      )
      = rearrange(
        [
            qw( RSCRIPT_PATH KEEP_METRICS SAMTOOLS_PATH KEEP_PLOT KEEP_RDATA NO_DUPS)
        ],
        @args
      );

    $self->rscript_path($rscript_path);
    $self->keep_metrics_file($keep_metrics);
    $self->keep_plot($keep_plot);
    $self->keep_rdata($keep_rdata);
    $self->samtools_path($samtools_path);
    $self->no_dups($no_dups);

    #setting defaults
    if ( !$self->rscript_path ) {
        $self->rscript_path('Rscript');
    }
    if ( !$self->program ) {
        $self->program('run_spp.R');
    }

    return $self;
}

sub run_ppqt {
    my ($self) = @_;

    my @cmd_args;
    my $job_name    = $self->job_name;
    my ($name,$dir) = fileparse($self->input_files->[0]);
  
    my $base_name   = $dir . '/' . $job_name;
    my $param_file  = $base_name . '.ppqt_metrics';
    my $plot_file   = $base_name . '.ppqt.pdf';
    my $r_file      = $base_name . '.ppqt.rdata';

    my $file = $self->processed_file() || $self->input_files->[0];

    check_file_exists($file);

    push @cmd_args, $self->rscript_path;
    push @cmd_args, $self->program;

    push @cmd_args, '-rf';
    push @cmd_args, '-c=' . $file;
    push @cmd_args, '-out=' . $param_file;

    if ( $self->keep_metrics_file ) {
        $self->output_files($param_file);
    }
    else {
        $self->created_files($param_file);
    }

    if ( $self->keep_plot ) {
        $self->output_files($plot_file);
        push @cmd_args, '-savp=' . $plot_file;
    }
    if ( $self->keep_rdata ) {
        $self->output_files($r_file);
        push @cmd_args, '-savd=' . $r_file;
    }
    
    my $path = $ENV{PATH};
    
    $ENV{PATH} = dirname($self->samtools_path).':'.$path;
    
    $self->execute_command_line( join( ' ', @cmd_args ) );
    
    $ENV{PATH} = $path;
    
    return $param_file;
}

sub parse_metrics {
    my ( $self, $param_file ) = @_;
    open my $fh, '<', $param_file or throw("could not open parameter file: $!");
    my $line = <$fh>;
    chomp $line;

    my %metrics;
    (
        $metrics{filename},    $metrics{numReads},
        $metrics{estFraglen},  $metrics{corr_estFragLen},
        $metrics{phantomPeak}, $metrics{corr_phantomPeak},
        $metrics{argmin_corr}, $metrics{min_corr},
        $metrics{nsc},         $metrics{rsc},
        $metrics{quality_tag}
    ) = split /\t/, $line;

    my %quality_decode = (
        '-2' => 'veryLow',
        '-1' => 'Low',
        '0'  => 'Medium',
        '1'  => 'High',
        '2'  => 'VeryHigh',
    );

    $metrics{quality_label} = $quality_decode{ $metrics{quality_tag} };

    return [ \%metrics ];
}

sub run_program {
    my ($self) = @_;

    $self->force_no_dups if ( $self->no_dups );
    my $params_file = $self->run_ppqt;
    my $metrics     = $self->parse_metrics($params_file);

    return $metrics;
}

sub force_no_dups {
    my ($self)       = @_;
    die "Need samtools path" unless ($self->samtools_path() && -e $self->samtools_path());
    my $temp_dir     = $self->get_temp_dir();
    my $base_name    = basename( $self->input_files->[0] );
    my $dedup_target = $temp_dir . '/' . $base_name;
    my @cmd_args     = (
        $self->samtools_path(), 'view', '-bF', 1024, $self->input_files->[0],
        '>', $dedup_target
    );
    $self->execute_command_line( join( ' ', @cmd_args ) );
    $self->processed_file($dedup_target);
}

sub processed_file {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{'processed_file'} = $arg;
    }
    return $self->{'processed_file'};
}

sub rscript_path {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{'rscript_path'} = $arg;
    }
    return $self->{'rscript_path'};
}

sub samtools_path {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{'samtools_path'} = $arg;
    }
    return $self->{'samtools_path'};
}

sub keep_metrics_file {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{'keep_metrics_file'} = $arg;
    }
    return $self->{'keep_metrics_file'};
}

sub keep_plot {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{'keep_plot'} = $arg;
    }
    return $self->{'keep_plot'};
}

sub keep_rdata {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{'keep_rdata'} = $arg;
    }
    return $self->{'keep_rdata'};
}

sub no_dups {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{'no_dups'} = $arg;
    }
    return $self->{'no_dups'};
}

1;
