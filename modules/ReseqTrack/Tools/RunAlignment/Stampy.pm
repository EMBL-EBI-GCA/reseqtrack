package ReseqTrack::Tools::RunAlignment::Stampy;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::RunAlignment;
use ReseqTrack::Tools::FileSystemUtils qw(delete_directory check_files_exists );

use File::Basename;
use File::stat;

our @ISA = qw(ReseqTrack::Tools::RunAlignment);

sub new {

    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);

    my (
        $genome_prefix, $hash_prefix,
    )
        = rearrange(
            [
                qw(
                      GENOME_PREFIX
                      HASH_PREFIX
                )
            ],
            @args
        );

    $self->genome_prefix($genome_prefix);
    $self->hash_prefix($hash_prefix);


    return $self;

}
#################################################################


sub run {
    my ($self) = @_;
    my $exit;

    $self->change_dir();

    $self->create_cmd_line();

    eval {
        $exit = system( $self->cmd_line );
        #print $self->cmd_line, "\n"; $exit=0;
        if ( $exit && $exit >= 1 ) {
            print STDERR ( "Failed to run " . $self->cmd_line );
        }
    };
    if ( $@ || $exit >= 1 ) {
        throw("Failed to run stampy: $@");
    }

    my $sam = $self->sam();
    my $bam = $self->create_bam_from_sam($sam);

    $self->output_files($bam);

    $self->delete_files();

    return 0 ;
}



############

sub create_cmd_line {
    my $self = shift;

    my $out_sam = $self->working_dir() . '/' . $$ . '.sam';
    $out_sam =~ s{//}{/};
    $self->sam($out_sam);

    my $cmd_line;
    $cmd_line = $self->program();
    $cmd_line .= " -g " . $self->genome_prefix();
    $cmd_line .= " -h " . $self->hash_prefix();
    $cmd_line .= " -o " . $out_sam;
    $cmd_line .= " -M ";

    if (! $self->skip_fragment){
        $cmd_line .= $self->fragment_file() . " "
            if ( $self->fragment_file);
    }

    if (! $self->skip_mate_files){
        $cmd_line .= $self->mate1_file() . " " if ( $self->mate1_file );
        $cmd_line .= $self->mate2_file() . " " if ( $self->mate2_file );
    }
            

    $self->cmd_line($cmd_line);

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

sub sam {
    my ( $self, $arg ) = @_;

    if ($arg) {
        $self->{sam} = $arg;
    }
    return $self->{sam};
}

sub cmd_line {
    my ( $self, $arg ) = @_;

    if ($arg) {
        $self->{cmd_line} = $arg;
    }
    return $self->{cmd_line};
}



1;

