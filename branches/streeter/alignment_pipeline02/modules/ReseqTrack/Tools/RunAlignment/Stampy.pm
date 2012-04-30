package ReseqTrack::Tools::RunAlignment::Stampy;

use strict;
use warnings;
use ReseqTrack::Tools::FileSystemUtils qw(check_file_does_not_exist);
use ReseqTrack::Tools::Argument qw(rearrange);
use List::Util qw (first);
use Env qw( @PATH );


use base qw(ReseqTrack::Tools::RunAlignment);

sub DEFAULT_OPTIONS { return {
        'gap_open_penalty' => '',
        'gap_extension_penalty' => '',
        'fast' => 0,
        'allow_overwriting' => 1,
        };
}

sub new {

    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);

    my ( $genome_prefix, $hash_prefix, $build_genome_flag,)
        = rearrange( [
            qw(
            GENOME_PREFIX HASH_PREFIX BUILD_GENOME_FLAG
                )], @args);


    #setting defaults
    if (!$self->program) {
      if ($ENV{STAMPY}) {
        $self->program($ENV{STAMPY} . '/stampy');
      }
      else {
        $self->program(first {-x $_} map {"$_/stampy"} @PATH);
      }
    }

    $self->genome_prefix($genome_prefix);
    $self->hash_prefix($hash_prefix);
    $self->build_genome_flag($build_genome_flag);

    return $self;

}
#################################################################


sub run_alignment {
    my ($self) = @_;

    if ($self->build_genome_flag){
        $self->build_genome();
        $self->build_hash();
    }

    if ($self->fragment_file) {
      $self->run_se_alignment;
    }

    if ($self->mate1_file && $self->mate2_file) {
      $self->run_pe_alignment;
    }

    return;
}

sub run_se_alignment {
    my $self = shift;
    return if !($self->fragment_file);

    my $sam = $self->working_dir() . '/'
        . $self->job_name
        . '_se.sam';
    $sam =~ s{//}{/};

    my @cmd_words = ($self->program);
    push(@cmd_words, '-g', $self->genome_prefix);
    push(@cmd_words, '-h', $self->hash_prefix);
    push(@cmd_words, '--gapopen='.$self->options('gap_open_penalty'))
            if $self->options('gap_open_penalty');
    push(@cmd_words, '--gapextend='.$self->options('gap_extension_penalty'))
            if ($self->options('gap_extension_penalty'));
    push(@cmd_words, '--fast') if ($self->options('fast'));
    push(@cmd_words, '--overwrite') if ($self->options('allow_overwriting'));

    if ($self->read_group_fields->{'ID'}) {
      my $rg_string = "--readgroup=ID:" . $self->read_group_fields->{'ID'};
      RG:
      while (my ($tag, $value) = each %{$self->read_group_fields}) {
        next RG if ($tag eq 'ID');
        next RG if (!$value);
        $rg_string .= ",$tag:$value";
      }
      push(@cmd_words, $rg_string);
    }
    push(@cmd_words, '-o', $sam);
    push(@cmd_words, '-M', $self->get_static_fastq('frag'));

    my $cmd_line = join(' ', @cmd_words);

    $self->sam_files($sam);
    $self->execute_command_line($cmd_line);

    return;
}

sub run_pe_alignment {
    my $self = shift;
    return if (!$self->mate1_file || !$self->mate2_file);

    my $sam = $self->working_dir() . '/'
        . $self->job_name
        . '_pe.sam';
    $sam =~ s{//}{/};

    my @cmd_words = ($self->program);
    push(@cmd_words, '-g', $self->genome_prefix);
    push(@cmd_words, '-h', $self->hash_prefix);
    push(@cmd_words, '--gapopen='.$self->options('gap_open_penalty'))
            if $self->options('gap_open_penalty');
    push(@cmd_words, '--gapextend='.$self->options('gap_extension_penalty'))
            if ($self->options('gap_extension_penalty'));
    push(@cmd_words, '--fast') if ($self->options('fast'));
    push(@cmd_words, '--overwrite') if ($self->options('allow_overwriting'));

    push(@cmd_words, '--insertsize='.$self->paired_length)
            if ($self->paired_length);

    if ($self->read_group_fields->{'ID'}) {
      my $rg_string = "--readgroup=ID:" . $self->read_group_fields->{'ID'};
      RG:
      while (my ($tag, $value) = each %{$self->read_group_fields}) {
        next RG if ($tag eq 'ID');
        next RG if (!$value);
        $rg_string .= ",$tag:$value";
      }
      push(@cmd_words, $rg_string);
    }
    push(@cmd_words, '-o', $sam);
    push(@cmd_words, '-M', map {$self->get_static_fastq($_)} ('mate1', 'mate2'));

    my $cmd_line = join(' ', @cmd_words);

    $self->sam_files($sam);
    $self->execute_command_line($cmd_line);

    return;
}

sub build_genome {
    my $self = shift;

    my $genome_prefix = $self->working_dir . '/'
                . $self->job_name;
    $genome_prefix =~ s{//}{/};

    my $cmd_line = $self->program;
    $cmd_line .= " -G " . $genome_prefix;
    $cmd_line .= " " . $self->reference;

    $self->genome_prefix($genome_prefix);
    my $genome = $genome_prefix . ".stidx";
    check_file_does_not_exist($genome);
    $self->created_files($genome);

    $self->execute_command_line($cmd_line);

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

    $self->hash_prefix($hash_prefix);
    my $hash = $hash_prefix . ".sthash";
    check_file_does_not_exist($hash);
    $self->created_files($hash);

    $self->execute_command_line($cmd_line);

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

