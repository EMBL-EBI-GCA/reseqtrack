package ReseqTrack::Tools::RunAlignment::Stampy;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);


use base qw(ReseqTrack::Tools::RunAlignment);

sub new {

    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);

    my ( $genome_prefix, $hash_prefix, $build_genome_flag, $se_options, $pe_options)
        = rearrange( [
            qw( GENOME_PREFIX HASH_PREFIX BUILD_GENOME_FLAG SE_OPTIONS PE_OPTIONS )
                ], @args);


    #setting defaults
    if (! $self->program) {
        $self->program('stampy');
    }

    $self->genome_prefix($genome_prefix);
    $self->hash_prefix($hash_prefix);
    $self->build_genome_flag($build_genome_flag);
    $self->se_options($se_options);
    $self->pe_options($pe_options);


    return $self;

}
#################################################################


sub run {
    my ($self) = @_;

    $self->change_dir();

    if ($self->build_genome_flag){
        $self->build_genome();
        $self->build_hash();
    }

    if ($self->fragment_file) {
        my $sam = $self->run_se_alignment();
        $self->sam_files($sam);
    }

    if ($self->mate1_file && $self->mate2_file) {
        my $sam = $self->run_pe_alignment();
        $self->sam_files($sam);
    }

    $self->run_samtools;

    return;
}

sub run_se_alignment {
    my $self = shift;

    my $sam = $self->working_dir() . '/'
        . $self->job_name
        . '_se.sam';
    $sam =~ s{//}{/};

    my $cmd_line;
    $cmd_line = $self->program();
    $cmd_line .= " -g " . $self->genome_prefix;
    $cmd_line .= " -h " . $self->hash_prefix;
    
    if ($self->se_options) {
        $cmd_line .= " " . $self->se_options;
    }

    $cmd_line .= " -o " . $sam;
    $cmd_line .= " -M " . $self->fragment_file;

    $self->execute_command_line($cmd_line);

    return $sam ;
}

sub run_pe_alignment {
    my $self = shift;

    my $sam = $self->working_dir() . '/'
        . $self->job_name
        . '_pe.sam';
    $sam =~ s{//}{/};

    my $cmd_line;
    $cmd_line = $self->program();
    $cmd_line .= " -g " . $self->genome_prefix;
    $cmd_line .= " -h " . $self->hash_prefix;
    
    if ($self->pe_options) {
        $cmd_line .= " " . $self->pe_options;
    }

    if ($self->paired_length ) {
        $cmd_line .= " --insertsize=" . $self->paired_length;
    }

    $cmd_line .= " -o " . $sam;
    $cmd_line .= " -M " . $self->mate1_file . " " . $self->mate2_file;

    $self->execute_command_line($cmd_line);

    return $sam ;
}

sub build_genome {
    my $self = shift;

    my $genome_prefix = $self->working_dir . '/'
                . $self->job_name;
    $genome_prefix =~ s{//}{/};

    my $cmd_line = $self->program;
    $cmd_line .= " -G " . $genome_prefix;
    $cmd_line .= " " . $self->reference;

    $self->execute_command_line($cmd_line);

    $self->genome_prefix($genome_prefix);

    my $genome = $genome_prefix . ".stidx";
    $self->files_to_delete($genome);

    return;

}

sub build_hash {
    my $self = shift;

    my $hash_prefix = $self->working_dir . '/'
                . $self->job_name;
    $hash_prefix =~ s{//}{/};

    my $cmd_line = $self->program;
    $cmd_line .= " -g " . $self->genome_prefix;
    $cmd_line .= " -H " . $hash_prefix;

    $self->execute_command_line($cmd_line);

    $self->hash_prefix($hash_prefix);

    my $hash = $hash_prefix . ".sthash";
    $self->files_to_delete($hash);

    return;

}



sub genome_prefix {
    my ( $self, $arg ) = @_;

    if ($arg) {
        $self->{genome_prefix} = $arg;
    }
    return $self->{genome_prefix};
}

sub hash_prefix {
    my ( $self, $arg ) = @_;

    if ($arg) {
        $self->{hash_prefix} = $arg;
    }
    return $self->{hash_prefix};
}

sub build_genome_flag {
    my ( $self, $arg ) = @_;

    if (defined $arg) {
        $self->{build_genome_flag} = $arg;
    }
    return $self->{build_genome_flag};
}

sub pe_options {
    my $pelf = shift;

    if (@_) {
        $pelf->{pe_options} = shift;
    }
    return $pelf->{pe_options};
}

sub se_options {
    my $self = shift;

    if (@_) {
        $self->{se_options} = shift;
    }
    return $self->{se_options};
}


1;

