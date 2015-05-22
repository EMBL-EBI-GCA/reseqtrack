package ReseqTrack::Tools::RunPeakCall::Hotspot;

use strict;
use warnings;
use vars qw(@ISA);

use File::Basename;
use File::Copy;
use IO::Compress::Gzip qw(gzip $GzipError);

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::RunProgram;
use File::chdir;

use base qw(ReseqTrack::Tools::RunPeakCall);

sub DEFAULT_OPTIONS {
  return {
    genome     => 'hg19',
    tag_length => 'K36',
    fdr        => 0.05,
    dupok      => 'T',      # allow duplicates
    chk_chr    => 'chrX',
  };
}

=pod

=head1 NAME

ReseqTrack::Tools::RunPeakCall::Hotspot
=head1 SYNOPSIS

Class for running the Hotspot peak caller. Developed with Hotspot v3 and will not work with earlier versions.

=head1 Example

my $hotspot = ReseqTrack::Tools::RunPeakCall::Hotspot(
					-input_files => '/path/to/file.bam'
					-working_dir => $output_dir,
					-options => \%options,
					-job_name => $name,
					-program = '/path/to/hotspot_dir/'
                      );
$hotspot->run;
my $output_file_list = $hotspot->output_files;

Note that the program should be the hotspot directory, not an executable. The module will look for scripts, executables and reference files under this location:
my $ref_data_dir  = "$hotspot_dir/data";
my $bin_dir       = "$hotspot_dir/hotspot-deploy/bin/";
my $hotspot_bin   = "$bin_dir/hotspot5";
my $wavePeaks_bin = "$bin_dir/wavePeaks";

=head1
Options:

Options are named to match the Hotspot documentation.

genome:   Which genome to map to. Defaults to hg19
tag_length:   Size of tags. Defaults to K36
genome and tag_length will be used ot find the chromosome file and mappability file in the hotspot reference data directory:
  my $chrom_file    = "$ref_data_dir/$genome.chromInfo.bed";
  my $mappable_file = "$ref_data_dir/$genome.$tag_length.mappable_only.bed";
dupok: Should duplicates be retained? (T or F). Defaults to T
fdr: Permitted false discovery rate. Defaults to 0.05
chr_chk: check chromsome

=cut

sub new {
  my ( $class, @args ) = @_;

  my $self = $class->SUPER::new(@args);

  return $self;
}

sub can_read_bam {
  return 1;
}

sub can_use_control {
  return 0;
}

sub control_required {
  return 0;
}

sub run_program {
  my ($self) = @_;

  $self->write_config_file();
  $self->write_runhotspot();
  $self->run_hotspot();
  $self->find_output();
}

sub write_config_file {
  my ($self) = @_;

  my $input_file = $self->input_files->[0];
  my $genome     = $self->options('genome');
  my $tag_length = $self->options('tag_length');

  my $dupok = $self->options('dupok');

  my $fdr      = $self->options('fdr');
  my $temp_dir = $self->get_temp_dir;
  my $chk_chr  = $self->options('chk_chr');

  my $hotspot_dir   = $self->program;
  my $ref_data_dir  = "$hotspot_dir/data";
  my $bin_dir       = "$hotspot_dir/hotspot-deploy/bin/";
  my $hotspot_bin   = "$bin_dir/hotspot5";
  my $wavePeaks_bin = "$bin_dir/wavePeaks";
  my $chrom_file    = "$ref_data_dir/$genome.chromInfo.bed";
  my $mappable_file = "$ref_data_dir/$genome.$tag_length.mappable_only.bed";

  $hotspot_bin   =~ s!//!/!g;
  $wavePeaks_bin =~ s!//!/!g;
  $chrom_file    =~ s!//!/!g;
  $mappable_file =~ s!//!/!g;

  my $config_file = "$temp_dir/runall.tokens.txt";

  throw("Cannot find input file: $input_file")      unless ( -e $input_file );
  throw("Cannot find chromosome file: $chrom_file") unless ( -e $chrom_file );
  throw("Cannot find mappable file: $mappable_file")
    unless ( -e $mappable_file );
  throw("Cannot find bin dir: $bin_dir")          unless ( -d $bin_dir );
  throw("Cannot execute hotspot: $hotspot_bin")   unless ( -x $hotspot_bin );
  throw("Cannot execute wavePeaks: $hotspot_bin") unless ( -x $wavePeaks_bin );

  open( my $fh, '>', $config_file )
    or throw("Could not write to $config_file: $!");

  print $fh <<END_TOKENS;
[script-tokenizer]

#######################################
## Notes:  If duplicate token definitions exist, the last definition in
##         the file will be used. Tokens can use .ini variables and declarations.
##         See http://docs.python.org/library/configparser.html
#######################################

#######################################
# Global tokens (used by most scripts)
#######################################

## Tags file in bam format (file extension .bam), or starched bed file
## (file extension .bed.starch).  If the latter, file should be in the
## location specified by _OUTDIR_.  If the former, the bed.starch file
## will be genereated, and put in _OUTDIR_.
_TAGS_: $input_file

## Genome
_GENOME_ = $genome
## Tag length
_K_ = $tag_length
## Chromosome coordinates, bed format.
_CHROM_FILE_ = $chrom_file
## Location of uniquely mappable positions in the genome for this tag length.
_MAPPABLE_FILE_ = $mappable_file

## Set DUPOK to T for DNaseI data, F for ChIP-seq data (DUPOK = T means allow duplicate reads)
_DUPOK_ = $dupok

## FDR levels. Set to N if you do not want FDR thresholding (for
## example, if you just want SPOT score computed.)
## _FDRS_ = "N"
_FDRS_ = "$fdr"

## Tag density, 150bp window, sliding every 20bp, used for
## peak-finding.  Will be generated, based on the _TAGS_ file, if it
## does not exist. Assumed to be starched bed file, extension
## bed.starch.  Can also be a directory, in which case the density file
## name will be assumed to be the name of the tags file, minus the bam
## or bed.starch extension, with the added extension
## tagdensity.bed.starch.
_DENS_: $temp_dir

## Output directories (can all be the same location).  Use full path names.
## _OUTDIR_ contains tags files in converted bed.starch and lib.txt formats (for hotspot
## program), and hotspot and peak results.
## _RANDIR_ contains generated random tags (for FDR thresholding) and hotspots called on random tags.
_OUTDIR_ = $temp_dir
_RANDIR_ = $temp_dir

## Set to T if you want scripts to skip steps that have already been done.
_CHECK_ = T
## If _CHECK_ = T, outputs are checked for completeness by searching
## for results for the following chromsome.
_CHKCHR_ = $chk_chr

## Hotspot program binary
_HOTSPOT_ = $hotspot_bin

## Peak-finding program.
_PKFIND_BIN_ = $wavePeaks_bin
## Peak-finding smoothing level. If the resolution of the input file
## is x, then the results are smoothed out to a scale of (2^level)*x.
_PKFIND_SMTH_LVL_ = 3

## Random number seed, used for generating random tags for FDR thresholding.
_SEED_=101

## Hotspot program parameters
_THRESH_ = 2
_WIN_MIN_ = 200
_WIN_MAX_ = 300
_WIN_INCR_ = 50
_BACKGRD_WIN_ = 50000
_MERGE_DIST_ = 150
_MINSIZE_ = 10
END_TOKENS

  close $fh;
  return $config_file;
}

sub write_runhotspot {
  my ($self) = @_;

  my $hotspot_dir      = $self->program;
  my $script_tokenizer = "$hotspot_dir/ScriptTokenizer/src/script-tokenizer.py";
  my $pipe_dir         = "$hotspot_dir/pipeline-scripts";

  throw("Cannot find the script tokenizer: $script_tokenizer")
    unless -e $script_tokenizer;
  throw("Cannot find pipeline script dir: $pipe_dir") unless -d $pipe_dir;

  my $temp_dir = $self->get_temp_dir;

  my @cmd_words = ($script_tokenizer);
  push @cmd_words, "--output-dir=$temp_dir";
  push @cmd_words, $temp_dir . '/runall.tokens.txt';
  push @cmd_words, map { $pipe_dir . '/' . $_ } $self->hotspot_scripts;

  my $cmd = join ' ', @cmd_words;

warn $cmd,"\n";

  $self->execute_command_line($cmd);

}

sub run_hotspot {
  my ($self) = @_;

  my $temp_dir = $self->get_temp_dir;
  
  for my $script ($self->hotspot_scripts) {
    local $CWD = $temp_dir; ## locally changing working dir to temp
    my $target = "$temp_dir/$script.tok";
   # my $target = "$script.tok";
    warn $CWD,"\t" ,$target,"\n";

    throw("Script has not been created: $target") unless -e $target;
    $self->execute_command_line($target);
  }

}

sub find_output {
  my ($self) = @_;

  my $temp_dir   = $self->get_temp_dir;
  my $input_file = $self->input_files->[0];
  my $fdr        = $self->options('fdr');

  my ( $basename, $path, $suffix ) = fileparse( $input_file, '.bed', '.bam' );

  my $output_dir = "$temp_dir/$basename-both-passes";
  throw("Cannot find output dir: $output_dir") unless ( -d $output_dir );

  my $hotspots_file = "$output_dir/$basename.hotspot.twopass.fdr$fdr.bed";
  my $peaks_file =
    "$output_dir/$basename.hotspot.twopass.fdr$fdr.merge.pks.bed";

  throw("Cannot find hotspots file: $hotspots_file")
    unless ( -e $hotspots_file );
  throw("Cannot find peak file: $peaks_file") unless ( -e $peaks_file );

## gzip the bed
  my $hotspots_file_gz = $hotspots_file. '.gz';
  my $peaks_file_gz = $peaks_file . '.gz';
  
  my $status_hotspot = gzip $hotspots_file => $hotspots_file_gz 
        or die "gzip failed: $GzipError\n";
        
  my $status_peak = gzip $peaks_file => $peaks_file_gz 
        or die "gzip failed: $GzipError\n";   

## remove bed file
  $self->created_files( $hotspots_file );
  $self->created_files( $peaks_file );

  my $hotspots_file_target =  
    $self->working_dir . '/' . $self->job_name . '.hotspots.bed.gz';

  my $peaks_file_target =
    $self->working_dir . '/' . $self->job_name . '.peaks.bed.gz';

  $self->output_files($hotspots_file_target);
  $self->output_files($peaks_file_target);

  move( $hotspots_file_gz, $hotspots_file_target );
  move( $peaks_file_gz,    $peaks_file_target );

}

sub hotspot_scripts {
  return qw(
    run_make_lib
    run_wavelet_peak_finding
    run_10kb_counts
    run_generate_random_lib
    run_pass1_hotspot
    run_pass1_merge_and_thresh_hotspots
    run_pass2_hotspot
    run_rescore_hotspot_passes
    run_spot
    run_thresh_hot.R
    run_both-passes_merge_and_thresh_hotspots
    run_add_peaks_per_hotspot
  );

}

sub output_hotspot_bed {
  my $self = shift;
  my @files = grep { /\.hotspots.bed.gz$/ } @{ $self->output_files };
  return $files[0];
}

sub output_peak_bed {
  my $self = shift;
  my @files = grep { /\.peaks.bed.gz$/ } @{ $self->output_files };
  return $files[0];
}

1;

