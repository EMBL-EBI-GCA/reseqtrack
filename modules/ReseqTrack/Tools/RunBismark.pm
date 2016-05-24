package ReseqTrack::Tools::RunBismark;

use strict;
use warnings;
use vars qw(@ISA);

use File::Basename;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use base qw(ReseqTrack::Tools::RunProgram);


=pod

=head1 NAME

ReseqTrack::Tools::RunBismark

=head1 SYNOPSIS

class for running Bismark. Child class of ReseqTrack::Tools::RunProgram

=cut

sub DEFAULT_OPTIONS { return {
        'n' => 1, # The maximum number of mismatches permitted in the "seed". Bowtie1 only option
	'algorithm' => 'bowtie', # bismark_mapper option: Should Bismark use bowtie or bowtie2
	'report_mode' => 0, #  bismark_methylation_extractor option: comprehensive is the default
	'bedgraph' => 1, #  bismark_methylation_extractor option.
	'nondirectional' => 0, #  bismark_methylation_extractor option. Set this to 1 if the library is non-directional
        };
}

sub new {
    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);	

    #setting defaults
    $self->program('bismark') unless ( $self->program );

    my ( $base, $cutoff, $fragment_file, $mate1_file, $mate2_file, $multicore, $rg_tag,$rg_id,$rg_sample,$run_mode,$reference )= rearrange( [qw(BASE CUTOFF FRAGMENT_FILE MATE1_FILE MATE2_FILE MULTICORE RG_TAG RG_ID RG_SAMPLE RUNMODE REFERENCE)], @args);

    $self->base($base);
    $self->cutoff($cutoff);
    $self->fragment_file($fragment_file);
    $self->mate1_file($mate1_file);
    $self->mate2_file($mate2_file);
    $self->multicore($multicore);
    $self->reference($reference);
    $self->rg_tag($rg_tag);     
    $self->rg_id($rg_id);
    $self->rg_sample($rg_sample);
    $self->runmode($run_mode);

    return $self;
}

=head2 run_alignment

  Arg [1]   : ReseqTrack::Tools::RunBismark
  Function  : uses Bismark to map the Bisulfite converted reads in $self->fragment_file 
              or $self->mate1_file/$self->mate2_file
              Output is files are stored in $self->output_files
              This method is called by the RunProgram run method

  Example   : my $bismark = ReseqTrack::Tools::RunBismark->new(
                      -base  => "basename",
                      -reference => '/path/to/reference',
                      -read_group_fields => {'ID' => 1, 'LB' => 'my_lib'},
                      -working_dir => '/path/to/dir/',
                      -mate1_file => '/path/to/mate1',
                      -mate2_file => '/path/to/mate2',
                      -fragment_file => '/path/to/fragment',
                      -rg_tag 1, # Write out a Read Group tag to the resulting SAM/BAM file. 
                      -rg_id SAMPLE1 # Sets the ID field in the @RG header line.
                      -rg_sample SAMPLE1 # Sets the ID field in the @RG header line.
                      );

              $bismark->run_alignment;

              my $output_file_list = $bismark->output_files; 
=cut

sub run_alignment {
    my ($self) = @_;

    my @cmd_args;

    my $algorithm = $self->options('algorithm');

    push @cmd_args, $self->program;
    push @cmd_args, "-n ".$self->options->{'n'};
    push @cmd_args, $self->reference;

    if ( $self->fragment_file ) {
        push @cmd_args, $self->fragment_file;
    }

    if ( $self->mate1_file && $self->mate2_file ) {
	push @cmd_args, "-1 ",$self->mate1_file,"-2 ",$self->mate2_file;
    }

    push @cmd_args, "--rg_tag " if $self->rg_tag();
    push @cmd_args, "--rg_id ".$self->rg_id() if $self->rg_id();
    push @cmd_args, "--rg_sample ".$self->rg_sample() if $self->rg_sample();
    push @cmd_args, "-B ".$self->base() if $self->base() && !$self->multicore(); #multicore and base can't be used simultaneously, wait until Bismark developers fix this
    push @cmd_args, "-o ".$self->working_dir();
    push @cmd_args, "--multicore ".$self->multicore() if $self->multicore();
    
    my $aln_cmd = join(' ', @cmd_args);
    $self->execute_command_line($aln_cmd);

    my @output_files;
    opendir DH, $self->working_dir;
    if ($self->base()) {
	my $base=$self->base;
	@output_files= grep {/${base}.+[\.bam|\.report.txt]/} readdir DH;
    } else {
        throw("I need a prefix name for output name. Please use Bismark's 'base' parameter");
    }
    throw("Could not retrieve the \@output_files") if scalar(@output_files)==0;

    @output_files=grep $_ !~/bismark/, @output_files;

    for (my $i=0;$i<scalar(@output_files);$i++) {
	my $res;
	$res=rename($self->working_dir."/".$output_files[$i],$self->working_dir."/".$self->base."_bismark.bam") if $output_files[$i]=~/\.bam$/;
	$res=rename($self->working_dir."/".$output_files[$i],$self->working_dir."/".$self->base."_bismark_report.txt") if $output_files[$i]=~/\.txt$/;
	throw("[ERROR] Rename failed!") if !$res;
    }
    closedir DH;
    
    $self->bam_file($self->base."_bismark.bam");
    $self->report_file($self->base."_bismark_report.txt");

    $self->output_format('BAM');
    return;
}

=head2 run_methylation_extractor

  Arg [1]   : ReseqTrack::Tools::RunBismark
  Function  : uses Bismark's bismark_methylation_extractor to do the methylation calls in BAM files passed with 
              $self->input_files(['file1.bam'])

              Output is files are stored in $self->output_files
              This method is called by the RunProgram run method

  Example   : my $bismark = ReseqTrack::Tools::RunBismark->new(
                      -input_files => ['file1.bam']
                      -working_dir => '/path/to/dir/',
                      -runmode => 'SINGLE'
                      );

              $bismark->run_methylation_extractor;

              my $output_file_list = $bismark->output_files;
=cut

sub run_methylation_extractor {
    my ($self) = @_;

    my @cmd_args;

    $self->output_format('methcall');

    $self->program('bismark_methylation_extractor');
    push @cmd_args, $self->program;

    if ($self->runmode) {
        throw("Don't recognise runmode ".$self->runmode.". Acceptable runmodes are: 'SINGLE', 'PAIRED'")
	    if ($self->runmode ne 'SINGLE' && $self->runmode ne 'PAIRED');
	push @cmd_args, "-s " if $self->runmode eq 'SINGLE';
	push @cmd_args, "-p " if $self->runmode eq 'PAIRED';
    }

    push @cmd_args, "--".$self->options->{'report_mode'}." " if $self->options->{'report_mode'};
    push @cmd_args, "--bedGraph " if $self->options->{'bedgraph'};
    push @cmd_args, "--non_directional" if $self->options->{'nondirectional'};
    push @cmd_args, "--cutoff ".$self->cutoff() if $self->cutoff;

    my $input= $self->input_files->[0];
    push @cmd_args, $input;

    my $aln_cmd = join(' ', @cmd_args);

    $self->execute_command_line($aln_cmd);
}

sub cpg_context {
    my ( $self, $arg ) = @_;

    if ( defined $arg ) {
        $self->{'cpg_context'} = $arg;
    } else {
        my $working_dir=$self->working_dir;
        my $input= $self->input_files->[0];

        my($filename, $directories, $suffix) = fileparse($input, qr/\.[^.]*/);

        opendir DH, $working_dir or throw("Cannot open $working_dir: $!");

        my @files= grep { /^CpG.+${filename}\.txt$/ } readdir DH;

        throw("unexpected number of splitting_report files: " . scalar @files) if (@files != 2 && @files !=4);

        closedir DH;

        @files=map{$working_dir."/".$_} @files;

        $self->{'cpg_context'}=\@files;

        return \@files;
    }
}

sub chh_context {
    my ( $self, $arg ) = @_;

    if ( defined $arg ) {
        $self->{'chh_context'} = $arg;
    } else {
        my $working_dir=$self->working_dir;
        my $input= $self->input_files->[0];

        my($filename, $directories, $suffix) = fileparse($input, qr/\.[^.]*/);

        opendir DH, $working_dir or throw("Cannot open $working_dir: $!");

        my @files= grep { /^CHH.+${filename}\.txt$/ } readdir DH;

        throw("unexpected number of splitting_report files: " . scalar @files) if (@files != 2 && @files !=4);

        closedir DH;

        @files=map{$working_dir."/".$_} @files;

        $self->{'chh_context'}=\@files;

        return \@files;
    }
}

sub chg_context {
    my ( $self, $arg ) = @_;

    if ( defined $arg ) {
        $self->{'chg_context'} = $arg;
    } else {
        my $working_dir=$self->working_dir;
        my $input= $self->input_files->[0];

        my($filename, $directories, $suffix) = fileparse($input, qr/\.[^.]*/);

        opendir DH, $working_dir or throw("Cannot open $working_dir: $!");

        my @files= grep { /^CHG.+${filename}\.txt$/ } readdir DH;

        throw("unexpected number of splitting_report files: " . scalar @files) if (@files != 2 && @files !=4);

        closedir DH;

	@files=map{$working_dir."/".$_} @files;

        $self->{'chg_context'}=\@files;

        return \@files;
    }
}

sub splitting {
    my ( $self, $arg ) = @_;

    if ( defined $arg ) {
        $self->{'splitting'} = $arg;
    } else {
        my $working_dir=$self->working_dir;
        my $input= $self->input_files->[0];

        my($filename, $directories, $suffix) = fileparse($input, qr/\.[^.]*/);

        opendir DH, $working_dir or throw("Cannot open $working_dir: $!");
	
	my @files= grep {/^${filename}_splitting_report\.txt$/} readdir DH;

        throw("unexpected number of splitting_report files: " . scalar @files) if @files != 1;

        closedir DH;

	$self->{'splitting'}=$working_dir."/".$files[0];

        return $working_dir."/".$files[0];
    }
}


sub mbias_txt {
    my ( $self, $arg ) = @_;

    if ( defined $arg ) {
        $self->{'mbias_txt'} = $arg;
    } else {
	my $working_dir=$self->working_dir;
	my $input= $self->input_files->[0];

	my($filename, $directories, $suffix) = fileparse($input, qr/\.[^.]*/);

	opendir DH, $working_dir or throw("Cannot open $working_dir: $!");
	my @files= grep {/^$filename.+M-bias\.txt$/} readdir DH;

	throw("unexpected number of mbias_txt files: " . scalar @files) if @files != 1;
   
	closedir DH;

	$self->{'mbias_txt'}=$working_dir."/".$files[0];

	return $working_dir."/".$files[0];
    }
}


sub mbias_png {
    my ( $self, $arg ) = @_;

    if ( defined $arg ) {
        $self->{'mbias_png'} = $arg;
    } else {
        my $working_dir=$self->working_dir;
        my $input= $self->input_files->[0];

        my($filename, $directories, $suffix) = fileparse($input, qr/\.[^.]*/);

        opendir DH, $working_dir or throw("Cannot open $working_dir: $!");

	my @files = grep { /^$filename.+M-bias.+\.png$/ } readdir DH;

        throw("unexpected number of mbias_png files: " . scalar @files) if @files != 1;

        closedir DH;
	
	$self->{'mbias_png'}=$working_dir."/".$files[0];
	
        return $working_dir."/".$files[0];
    }
}

sub bedgraph {
    my ( $self, $arg ) = @_;

    if ( defined $arg ) {
        $self->{'bedgraph'} = $arg;
    } else {
        my $working_dir=$self->working_dir;
        my $input= $self->input_files->[0];

        my($filename, $directories, $suffix) = fileparse($input, qr/\.[^.]*/);

        opendir DH, $working_dir or throw("Cannot open $working_dir: $!");

        my @files = grep { /^$filename.+bedGraph\.gz$/ } readdir DH;

        throw("unexpected number of mbias_png files: " . scalar @files) if @files != 1;

        closedir DH;

        $self->{'bedgraph'}=$working_dir."/".$files[0];

        return $working_dir."/".$files[0];
    }
}

=head2 run_program

  Arg [1]   : ReseqTrack::Tools::RunBismark
  Arg [2]   : string, command, must be one of the following:
              aln, methext
  Function  : uses Bismark to process the files in $self->input_files.
              Output is files are stored in $self->output_files
              This method is called by the RunProgram run method
  Returntype:
  Exceptions: Throws if the command is not recognised.
  Example   : $self->run();

=cut

sub run_program {
    my ( $self, $command ) = @_;

    my %subs = (
        'aln'         => \&run_alignment,
        'methext'          => \&run_methylation_extractor	
	);

    throw("Did not recognise command $command") if ( !defined $subs{$command} );

    return &{ $subs{$command} }($self);
}

sub cutoff {
    my ( $self, $arg ) = @_;

    if ($arg) {
        $self->{cutoff} = $arg;
    }
    return $self->{cutoff};
}


sub bam_file {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
	$self->{'bam_file'} = $arg;
    }
    return $self->{'bam_file'};
}

sub report_file {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
	$self->{'report_file'} = $arg;
    }
    return $self->{'report_file'};
}

sub fragment_file {
    my ( $self, $arg ) = @_;

    if ($arg) {
        $self->{fragment_file} = $arg;
    }
    return $self->{fragment_file};
}

sub mate1_file {
    my ( $self, $arg ) = @_;

    if ($arg) {
        $self->{mate1_file} = $arg;
    }
    return $self->{mate1_file};
}

sub mate2_file {
    my ( $self, $arg ) = @_;

    if ($arg) {
        $self->{mate2_file} = $arg;
    }
    return $self->{mate2_file};
}

sub multicore {
    my ( $self, $arg ) = @_;

    if ($arg) {
        $self->{multicore} = $arg;
    }
    return $self->{multicore};
}

sub runmode {
    my ( $self, $arg ) = @_;

    if ($arg) {
        $self->{runmode} = $arg;
    }
    return $self->{runmode};
}

sub reference {
    my ( $self, $arg ) = @_;

    if ($arg) {
        $self->{reference} = $arg;
    }
    return $self->{reference};
}


sub base {
    my ( $self, $arg ) = @_;

    if ($arg) {
        $self->{base} = $arg;
    }
    return $self->{base};
}

sub rg_tag {
    my ( $self, $arg ) = @_;

    if ($arg) {
        $self->{rg_tag} = $arg;
    }
    return $self->{rg_tag};
}

sub rg_id {
    my ( $self, $arg ) = @_;

    if ($arg) {
        $self->{rg_id} = $arg;
    }
    return $self->{rg_id};
}

sub rg_sample {
    my ( $self, $arg ) = @_;

    if ($arg) {
        $self->{rg_sample} = $arg;
    }
    return $self->{rg_sample};
}

sub output_format {
    my ( $self, $arg ) = @_;

    if ($arg) {
        $self->{output_format} = $arg;
    }
    return $self->{output_format};
}


1;
