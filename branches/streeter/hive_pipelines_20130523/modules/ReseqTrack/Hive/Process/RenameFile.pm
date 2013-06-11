
package ReseqTrack::Hive::Process::RenameFile;

use strict;
use warnings;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::FileSystemUtils qw(check_directory_exists);
use ReseqTrack::Tools::Exception qw(throw);
use File::Copy qw( move );


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;

    my $param_name = $self->param_required('file_param_name');
    $self->param_required($param_name);
    my $files = $self->get_param_values($param_name);
    throw("too many files: ".join(' ', @$files)) if @$files >1;

    my $suffix = $self->param_required('suffix');

    my $output_dir = $self->output_dir;
    my $job_name = $self->job_name;

    my $old_filename = $files->[0];
    throw("do not have a file") if !$old_filename;

    my $new_filename = "$output_dir/$job_name.$suffix";
    $new_filename =~ s{//+}{/}g;
    $old_filename =~ s{//+}{/}g;
    throw("this module shouldn't be used if new_file and old_file are the same") if $new_filename eq $old_filename;

    check_directory_exists($output_dir);
    move($old_filename, $new_filename) or throw("could not move $old_filename to $new_filename $!");

    $self->output_param($param_name, $new_filename);

}


1;

