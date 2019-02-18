package ReseqTrack::Tools::GetFastq::ENAAspera;

use strict;
use warnings;

use ReseqTrack::Hive::Process::Aspera;
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
                      -output_dir => '/path/to/dir'
                      -run_info => $my_run,
                      -db => $era_db,
                      -source_root_dir => '/mount/ena/dir',
                      );
$fastq_getter->run;
my $output_file_list = $fastq_getter->output_files;

=cut


sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);

    my ($ena_ascphost, $source_root_dir)
        = rearrange([ qw(ENA_ASCPHOST SOURCE_ROOT_DIR) ], @args);

    $self->ena_ascphost($ena_ascphost || 'fasp.sra.ebi.ac.uk:');

    # This overrides the default set in ReseqTrack::Tools::GetFastq
    $self->source_root_dir($source_root_dir || '/');

    return $self;
}


sub get_files {
    my $self = shift;

    my $output_hash = $self->output_hash;
    while (my ($source_path, $output_path) = each %$output_hash) {
        my $aspera_getter = ReseqTrack::Hive::Process::Aspera -> new(
            # Required:
            -filename     => $source_path,
            #-ascp_exe => ,
            -username     => "fasp",
            -aspera_url   => $self->ena_ascphost,
            -work_dir     => "",
            # Optional:
            -download_dir => $output_path,
            # -trim_path    =>
            # -ascp_param => ,
        );
        $aspera_getter->run;
    }

}


sub ena_ascphost {
    my ($self, $ena_ascphost) = @_;
    if (defined($ena_ascphost)) {
        $self->{ena_ascphost} = $ena_ascphost;
    }
    return $self->{ena_ascphost};
}

1;


