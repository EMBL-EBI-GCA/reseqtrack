package ReseqTrack::Tools::GetFastq::ENAAspera;

use strict;
use warnings;

use ReseqTrack::Tools::Aspera;
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::ERAUtils qw();

use base qw(ReseqTrack::Tools::GetFastq);

=pod

=head1 NAME

ReseqTrack::Tools::GetFastq::ENAAspera

=head1 SYNOPSIS

class for getting fastq files using Aspera protocol. Child class of ReseqTrack::Tools::GetFastq

=head1 Example

my $fastq_getter = ReseqTrack::Tools::GetFastq::ENAAspera(
                      -output_dir => '/path/to/dir',
                      -run_info => $my_run,
                      -db => $era_db,
                      -source_root_dir => '/mount/ena/dir',
                      );
$fastq_getter->run();
my $output_file_list = $fastq_getter->output_files;

=cut


sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);

    my ($ena_ascphost, $ascp_param, $source_root_dir)
        = rearrange([ qw(ENA_ASCPHOST ASCP_PARAM SOURCE_ROOT_DIR) ], @args);

    $self->ena_ascphost($ena_ascphost || 'fasp.sra.ebi.ac.uk');
    $self->ascp_param($ascp_param || {'P' => 33001});

    # This overrides the default set in ReseqTrack::Tools::GetFastq
    $self->source_root_dir($source_root_dir || '/');

    return $self;
}


sub get_files {
    my $self = shift;

    my $output_hash = $self->output_hash;
    my $aspera_getter = ReseqTrack::Tools::Aspera -> new(
            -username     => "fasp",
            -aspera_url   => $self->ena_ascphost,
            -ascp_param   => $self->ascp_param
    );
    while (my ($source_path, $output_path) = each %$output_hash) {
        $aspera_getter->run_download(
            -remote_path  => $source_path,
            -local_path   => $output_path,
        );
    }

}


sub ena_ascphost {
    my ($self, $ena_ascphost) = @_;
    if (defined($ena_ascphost)) {
        $self->{ena_ascphost} = $ena_ascphost;
    }
    return $self->{ena_ascphost};
}


sub ascp_param {
    my ($self, $ascp_param) = @_;
    if (defined($ascp_param)) {
        $self->{ascp_param} = $ascp_param;
    }
    return $self->{ascp_param};
}

1;
