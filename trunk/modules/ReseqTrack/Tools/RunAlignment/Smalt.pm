package ReseqTrack::Tools::RunAlignment::Smalt;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename qw(fileparse);


use base qw(ReseqTrack::Tools::RunAlignment);

sub new {

    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);

    my ( $index_prefix, $build_index_flag, $se_options, $pe_options)
        = rearrange( [
            qw( INDEX_PREFIX BUILD_index_FLAG SE_OPTIONS PE_OPTIONS )
                ], @args);



    $self->index_prefix($index_prefix);
    $self->build_index_flag($build_index_flag);
    $self->se_options($se_options);
    $self->pe_options($pe_options);


    return $self;

}
#################################################################


sub run {
    my ($self) = @_;

    $self->change_dir();

    if ($self->build_index_flag){
        $self->build_index();
    }

    if ($self->fragment_file) {
        my $sam = $self->run_se_alignment();
        $self->sam_files($sam);
    }

    if ($self->mate1_file && $self->mate2_file) {
        my $sam = $self->run_pe_alignment();
        $self->sam_files($sam);
    }

    $self->run_samtools(0);

    return;
}

sub run_se_alignment {
    my $self = shift;

    my $sam = $self->working_dir() . '/'
        . $self->job_name
        . '_se.sam';
    $sam =~ s{//}{/};

    my $frag_file = $self->fragment_file;
    if ($frag_file =~ /\.gz$/) {
        $frag_file = $self->decompress_file($frag_file);
    }

    my $cmd_line = $self->program . ' map';
    $cmd_line .= ' ' . $self->se_options if $self->se_options;
    $cmd_line .= ' -f sam -o ' . $sam;
    $cmd_line .= ' ' . $self->index_prefix;
    $cmd_line .= ' ' . $frag_file;

    $self->execute_command_line($cmd_line);

    return $sam ;
}

sub run_pe_alignment {
    my $self = shift;

    my $sam = $self->working_dir() . '/'
        . $self->job_name
        . '_pe.sam';
    $sam =~ s{//}{/};

    my $mate1_file = $self->mate1_file;
    if ($mate1_file =~ /\.gz$/) {
        $mate1_file = $self->decompress_file($mate1_file);
    }

    my $mate2_file = $self->mate2_file;
    if ($mate2_file =~ /\.gz$/) {
        $mate2_file = $self->decompress_file($mate2_file);
    }


    my $cmd_line = $self->program . ' map';
    $cmd_line .= ' ' . $self->se_options if $self->se_options;
    $cmd_line .= ' -f sam -o ' . $sam;
    $cmd_line .= ' ' . $self->index_prefix;
    $cmd_line .= ' ' . $mate1_file . ' ' . $mate2_file;

    $self->execute_command_line($cmd_line);

    return $sam ;
}

sub build_index {
    my $self = shift;

    my $index_prefix = $self->working_dir . '/'
                    . $self->job_name;
    $index_prefix =~ s{//}{/};

    my $reference = $self->reference;
    if ($reference =~ /\.gz$/) {
        $reference = $self->decompress_file($reference);
    }

    my $cmd_line = $self->program;
    $cmd_line .= ' index ' . $index_prefix . ' ' . $reference;

    $self->execute_command_line($cmd_line);

    $self->index_prefix($index_prefix);

    my $sma = $index_prefix . '.sma';
    my $smi = $index_prefix . '.smi';

    $self->files_to_delete([$sma, $smi]);

    return;

}

sub decompress_file {
    my ($self, $file) = @_;

    my $name =  fileparse($file, qr/\.gz/);
    my $file_gunzipped = $self->working_dir . '/' . $name;
    $file_gunzipped =~ s{//}{/}g;

    my $gunzip_cmd = 'gunzip -c ' . $file . ' > ' . $file_gunzipped;
    $self->execute_command_line($gunzip_cmd);
    
    $self->files_to_delete($file_gunzipped);
    return $file_gunzipped;
}




sub index_prefix {
    my ( $self, $arg ) = @_;

    if ($arg) {
        $self->{index_prefix} = $arg;
    }
    return $self->{index_prefix};
}


sub build_index_flag {
    my ( $self, $arg ) = @_;

    if (defined $arg) {
        $self->{build_index_flag} = $arg;
    }
    return $self->{build_index_flag};
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

