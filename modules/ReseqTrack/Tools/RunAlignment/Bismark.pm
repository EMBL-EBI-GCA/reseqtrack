package ReseqTrack::Tools::RunAlignment::Bismark;

use strict;
use warnings;
use vars qw(@ISA);

use File::Basename;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::RunAlignment;

@ISA = qw(ReseqTrack::Tools::RunAlignment);

=pod

=head1 NAME

ReseqTrack::Tools::RunAlignment::Bismark

=head1 SYNOPSIS

class for running Bismark. Child class of ReseqTrack::Tools::RunAlignment

=cut

sub DEFAULT_OPTIONS { return {
        'n' => 1, # The maximum number of mismatches permitted in the "seed". Bowtie1 only option
	'algorithm' => 'bowtie', # Should Bismark use bowtie or bowtie2
	'report_mode' => 'comprehensive',
	'bedgraph' => 1
        };
}

sub new {
    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);	

    #setting defaults
    $self->program('bismark') unless ( $self->program );

    my ( $base,$rg_tag,$rg_id,$rg_sample,$run_mode )= rearrange( [qw(BASE RG_TAG RG_ID RG_SAMPLE RUNMODE)], @args);

    $self->base($base);
    $self->rg_tag($rg_tag);     
    $self->rg_id($rg_id);
    $self->rg_sample($rg_sample);
    $self->runmode($run_mode);

    return $self;
}

=head2 run_alignment

  Arg [1]   : ReseqTrack::Tools::RunAlignment::Bismark
  Function  : uses Bismark to map the Bisulfite converted reads in $self->fragment_file 
              or $self->mate1_file/$self->mate2_file
              Output is files are stored in $self->output_files
              This method is called by the RunProgram run method

  Example   : my $bismark = ReseqTrack::Tools::RunAlignment::Bismark(
                      -base  => "basename",
                      -reference => '/path/to/reference',
                      -read_group_fields => {'ID' => 1, 'LB' => 'my_lib'},
                      -working_dir => '/path/to/dir/',
                      -mate1_file => '/path/to/mate1',
                      -mate2_file => '/path/to/mate2',
                      -fragment_file => '/path/to/fragment',
                      -rg_tag 1,
                      -rg_id SAMPLE1
                      -rg_sample SAMPLE1
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
    push @cmd_args, "-B ".$self->base() if $self->base();
    push @cmd_args, "-o ".$self->working_dir();
    
    my $aln_cmd = join(' ', @cmd_args);

    $self->execute_command_line($aln_cmd);

    my $output_file;
    if ($self->base()) {
	$output_file.=$self->working_dir."/".$self->base().".bam";
    } else {
	my($filename, $directories) = fileparse($self->fragment_file);
	$output_file.=$self->working_dir."/".$filename."_bismark_bt2.bam";
    }

    $self->output_format('BAM');
    $self->output_files($output_file);
    return;
}

=head2 run_methylation_extractor

  Arg [1]   : ReseqTrack::Tools::RunAlignment::Bismark
  Function  : uses Bismark's bismark_methylation_extractor to do the methylation calls in BAM files passed with 
              $self->input_files(['file1.bam'])

              Output is files are stored in $self->output_files
              This method is called by the RunProgram run method

  Example   : my $bismark = ReseqTrack::Tools::RunAlignment::Bismark(
                      -input_files => ['file1.bam']
                      -working_dir => '/path/to/dir/',
                      -single => 1
                      );

              $bismark->run_methylation_extractor;

              my $output_file_list = $bismark->output_files;
=cut

sub run_methylation_extractor {
    my ($self) = @_;

    my @cmd_args;

    $self->program('bismark_methylation_extractor');
    push @cmd_args, $self->program;

    if ($self->runmode) {
        throw("Don't recognise runmode ".$self->runmode.". Acceptable runmodes are: 'SINGLE', 'PAIRED'")
	    if ($self->runmode ne 'SINGLE' && $self->runmode ne 'PAIRED');
	push @cmd_args, "-s " if $self->runmode eq 'SINGLE';
	push @cmd_args, "-p " if $self->runmode eq 'PAIRED';
    }

    push @cmd_args, "--".$self->options->{'report_mode'}." ";
    push @cmd_args, "--bedGraph " if $self->options->{'bedgraph'};

    my $input= $self->input_files->[0];
    push @cmd_args, $input;

    my $aln_cmd = join(' ', @cmd_args);

    $self->execute_command_line($aln_cmd);

    my($filename, $directories, $suffix) = fileparse($input, qr/\.[^.]*/);

    my @output_files;
    push @output_files, $self->working_dir."/".$filename.".M-bias_R1.png";
    push @output_files, $self->working_dir."/".$filename.".M-bias.txt";
    push @output_files, $self->working_dir."/".$filename."_splitting_report.txt";
    push @output_files, $self->working_dir."/"."CHG_context_".$filename.".txt";
    push @output_files, $self->working_dir."/"."CHH_context_".$filename.".txt";
    push @output_files, $self->working_dir."/"."CpG_context_".$filename.".txt";
    push @output_files, $self->working_dir."/".$filename.".bedGraph.gz";

    $self->output_files(\@output_files);
}

sub runmode {
    my ( $self, $arg ) = @_;

    if ($arg) {
        $self->{runmode} = $arg;
    }
    return $self->{runmode};
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

1;
