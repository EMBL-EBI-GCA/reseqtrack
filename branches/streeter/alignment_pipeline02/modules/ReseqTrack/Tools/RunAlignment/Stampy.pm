package ReseqTrack::Tools::RunAlignment::Stampy;

use strict;
use warnings;
use ReseqTrack::Tools::FileSystemUtils qw(check_file_does_not_exist);
use ReseqTrack::Tools::Argument qw(rearrange);
use List::Util qw (first);
use Env qw( @PATH );


use base qw(ReseqTrack::Tools::RunAlignment);

sub new {

    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);

    my ( $genome_prefix, $hash_prefix, $build_genome_flag, $se_options, $pe_options)
        = rearrange( [
            qw( GENOME_PREFIX HASH_PREFIX BUILD_GENOME_FLAG SE_OPTIONS PE_OPTIONS )
                ], @args);


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
    $self->se_options($se_options);
    $self->pe_options($pe_options);


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
    if (-e $sam) {
      delete_file($sam);
    }

    my $cmd_line;
    $cmd_line = $self->program();
    $cmd_line .= " -g " . $self->genome_prefix;
    $cmd_line .= " -h " . $self->hash_prefix;
    
    if ($self->se_options) {
        $cmd_line .= " " . $self->se_options;
    }

    if ($self->read_group_fields->{'ID'}) {
      $cmd_line .= " --readgroup=ID:" . $self->read_group_fields->{'ID'};
      RG:
      while (my ($tag, $value) = each %{$self->read_group_fields}) {
        next RG if ($tag eq 'ID');
        next RG if (!$value);
        $cmd_line .= ",$tag:$value";
      }
    }

    $cmd_line .= " -o " . $sam;
    $cmd_line .= " -M " . $self->fragment_file;

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
    if (-e $sam) {
      delete_file($sam);
    }

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

    if ($self->read_group_fields->{'ID'}) {
      $cmd_line .= " --readgroup=ID:" . $self->read_group_fields->{'ID'};
      RG:
      while (my ($tag, $value) = each %{$self->read_group_fields}) {
        next RG if ($tag eq 'ID');
        next RG if (!$value);
        $cmd_line .= ",$tag:$value";
      }
    }

    $cmd_line .= " -o " . $sam;
    $cmd_line .= " -M " . $self->mate1_file . " " . $self->mate2_file;

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

