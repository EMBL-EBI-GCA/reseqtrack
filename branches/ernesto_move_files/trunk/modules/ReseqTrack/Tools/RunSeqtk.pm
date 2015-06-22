
=pod

=head1 NAME

ReseqTrack::Tools::RunSeqtk

=head1 SYNOPSIS

This is a class for running seqtk
Right now only the fqchk function has been implemented; other functionalities can be added easily if needed.
It is a sub class of a ReseqTrack::Tools::RunProgram.

=head1 Example


=cut

package ReseqTrack::Tools::RunSeqtk;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename qw(fileparse);
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists);

use base qw(ReseqTrack::Tools::RunProgram);


sub new {
    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);

    #setting defaults
    if ( !$self->program ) {
        if ( $ENV{SEQTK} ) {
            $self->program( $ENV{SEQTK} . '/seqtk' );
        }
        else {
            $self->program('seqtk');
        }
    }

    return $self;
}

sub run_fqchk {

    my ($self) = @_;

    foreach my $input ( @{ $self->input_files } ) {

        my $prefix = fileparse( $input);
        my $output_file = $self->working_dir . "/$prefix.tmp";
        $output_file =~ s{//}{/}g;

        my @cmds;
        push @cmds, $self->program . " fqchk ";
        push @cmds, $input;
        push @cmds, ">" . $output_file;

        my $cmd = join( ' ', @cmds );

        $self->output_files($output_file);
        $self->execute_command_line($cmd);
    }
}

1;
