package ReseqTrack::Tools::RunWiggler;

use strict;
use warnings;
use vars qw(@ISA);
use Carp;

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename;
use File::Path qw(make_path);
use File::Copy;
use base qw(ReseqTrack::Tools::RunProgram);

=pod

=head1 NAME

ReseqTrack::Tools::RunWiggler;

=head1 SYNOPSIS

module for converting bam into wig using align2rawsignal aka wiggler

See docs at https://code.google.com/p/align2rawsignal/ to ensure environment variables are set
=head1 Example 1

my $conv = ReseqTrack::Tools::RunWiggler(
                      -input => '/path/to/bam_file'
                      -chromosomes => '/path/to/chromosomes',
                      -output_format => 'bw'
);
$conv->run;
my ($output_file) = @{$conv->output_files};

=head1 Example 2

my $cram = ReseqTrack::Tools::RunWiggler->new (
	-input_files => $bam,
	-output_format => 'bg',
	-program => $input{program},
);

=cut

sub DEFAULT_OPTIONS {
    return {
    };
}

sub new {
    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);

    my ( $output_format, $bedGraphToBigWig_path, $chrom_sizes_file, $dedupe, $clobber, $mcr_root, $samtools ) =
      rearrange(
        [qw(OUTPUT_FORMAT BEDGRAPHTOBIGWIG_PATH CHROM_SIZES_FILE DEDUPE CLOBBER MCR_ROOT SAMTOOLS)],
        @args );

    $self->output_format($output_format);
    $self->bedGraphToBigWig_path($bedGraphToBigWig_path);
    $self->chrom_sizes_file($chrom_sizes_file);
    $self->dedupe($dedupe);
    $self->clobber($clobber);
    $self->mcr_root($mcr_root);
    $self->samtools($samtools);

    return $self;
}

sub run_program {
    my ($self) = @_;

    my @cmd = ( $self->program );

    if ( $self->output_format eq 'bw' ) {
        die("Need bedGraphToBigWig path to encode output as bw")
          unless ( $self->bedGraphToBigWig_path );

        die("Need chrom sizes file path to encode output as bw")
          unless $self->chrom_sizes_file;
    }

    my $job_name      = $self->job_name();
    my $temp_dir      = $self->get_temp_dir();
    my $final_dir     = $self->working_dir();
    my $format_suffix = $self->output_format();
    my $mcr_root      = $self->mcr_root();
    my $samtools      = $self->samtools ? $self->samtools : 'samtools';

    my $mcr_jre      = "$mcr_root/sys/java/jre/glnxa64/jre/lib/amd64";
    my $xapplres_dir = "$mcr_root/X11/app-defaults";
    my $ld_lib_path  = join( ':',
       $ENV{LD_LIBRARY_PATH},     "$mcr_root/runtime/glnxa64",
       "$mcr_root/bin/glnxa64",   "$mcr_root/sys/os/glnxa64",
       "$mcr_jre/native_threads", "$mcr_jre/server",
       $mcr_jre );

    $ENV{MCRROOT}         = $mcr_root;
    $ENV{MCRJRE}          = $mcr_jre;
    $ENV{XAPPLRESDIR}     = $xapplres_dir;
    $ENV{LD_LIBRARY_PATH} = $ld_lib_path; 
    $ENV{MCR_CACHE_ROOT} = $temp_dir . '/MCR_CACHE';
    make_path( $ENV{MCR_CACHE_ROOT} );

    my $output = "$temp_dir/$job_name.$format_suffix";
    my $dest   = "$final_dir/$job_name.a2rs.$format_suffix";
    
    throw("$dest already exists") if (-e $dest && !$self->clobber);

    my @input_files = @{ $self->input_files };

    if ( $self->dedupe() ) {
        print "DEDUPE$/";
        for ( my $i = 0 ; $i < scalar(@input_files) ; $i++ ) {
            my $file      = $input_files[$i];
            my $temp_file = "$temp_dir/temp$i.bam";
            $self->execute_command_line( "$samtools view -b -F 1024 $file > $temp_file" );
            $input_files[$i] = $temp_file;

        }
    }

    for my $input_file (@input_files) {
        push @cmd, '-i=' . $input_file;
    }

    if ( $self->options ) {
        while ( my ( $tag, $value ) = each %{ $self->options } ) {
            push( @cmd, "-$tag=$value" );
        }
    }

    if ( $self->output_format eq 'bw' ) {
        push @cmd, '-of=bg';
        my $intermediate = "$temp_dir/temp.bg";
        push @cmd, "-o=$intermediate";
        push @cmd,
          (
            ';', $self->bedGraphToBigWig_path,
            $intermediate, $self->chrom_sizes_file, $output
          );
    }
    else {
        push @cmd, '-of=' . $self->output_format if $self->output_format;
        push @cmd, '-o=' . $output;
    }

    my $cmd = join ' ', @cmd;
    $self->execute_command_line($cmd);
    move( $output, $dest );
    $self->output_files($dest);
}

sub output_format {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{'output_format'} = $arg;
    }

    return $self->{'output_format'};
}

sub bedGraphToBigWig_path {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{'bedGraphToBigWig_path'} = $arg;
    }

    return $self->{'bedGraphToBigWig_path'};
}

sub chrom_sizes_file {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{'chrom_sizes_file'} = $arg;
    }

    return $self->{'chrom_sizes_file'};
}

sub dedupe {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{'dedupe'} = $arg;
    }

    return $self->{'dedupe'};
}

sub clobber {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{'clobber'} = $arg;
    }

    return $self->{'clobber'};
}

sub mcr_root {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{'mcr_root'} = $arg;
    }

    return $self->{'mcr_root'};
}

sub samtools {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{'samtools'} = $arg;
    }

    return $self->{'samtools'};
}

1;
