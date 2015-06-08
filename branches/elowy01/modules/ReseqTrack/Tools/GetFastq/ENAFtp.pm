package ReseqTrack::Tools::GetFastq::ENAFtp;

use strict;
use warnings;

use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::ERAUtils qw ();
use Net::FTP;

use base qw(ReseqTrack::Tools::GetFastq);

=pod

=head1 NAME

ReseqTrack::Tools::GetFastq::ENAFtp

=head1 SYNOPSIS

class for getting fastq files using Ftp protocol. Child class of ReseqTrack::Tools::GetFastq

=head1 Example

my $fastq_getter = ReseqTrack::Tools::GetFastq::ENAFtp(
                      -output_dir => '/path/to/dir'
                      -run_info => $my_run,
                      -db => $era_db,
                      -source_root_dir => '/mount/ena/dir',
                      );
$fastq_getter->run;
my $output_file_list = $fastq_getter->output_files;

=cut


sub new {
    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);

    my ( $ena_ftphost, $source_root_dir)
        = rearrange( [ qw( ENA_FTPHOST SOURCE_ROOT_DIR)], @args);

    $self->ena_ftphost($ena_ftphost || 'ftp.sra.ebi.ac.uk');

    # This overrides the default set in ReseqTrack::Tools::GetFastq
    $self->source_root_dir($source_root_dir || '/');

    return $self;
}


sub get_files {
  my $self = shift;
  
  my $ftp = Net::FTP->new($self->ena_ftphost)
    or throw('cannot connect to '.$self->ena_ftphost.": $@");
  $ftp->login('anonymous')
    or throw('cannot login to '.$self->ena_ftphost.' as anonymous: '.$ftp->message);
  $ftp->binary;

  my $output_hash = $self->output_hash;
  while (my($source_path, $output_path) = each %$output_hash) {
    $ftp->get($source_path, $output_path)
        or throw("Failed ftp get from $source_path on ".$self->ena_ftphost." to $output_path: ".$ftp->message);
    chmod 0644, $output_path;
    $self->output_files($output_path);
  }
}


sub ena_ftphost{
  my ($self, $ena_ftphost) = @_;
  if(defined($ena_ftphost)){
    $self->{ena_ftphost} = $ena_ftphost;
  }
  return $self->{ena_ftphost};
}

1;
