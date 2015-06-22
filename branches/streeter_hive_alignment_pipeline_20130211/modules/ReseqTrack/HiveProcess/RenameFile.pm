
package ReseqTrack::HiveProcess::RenameFile;

use strict;
use warnings;

use base ('ReseqTrack::HiveProcess::BranchableProcess');
use ReseqTrack::Tools::FileSystemUtils qw(check_directory_exists);
use ReseqTrack::Tools::Exception qw(throw);
use File::Copy qw( move );


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;
    my $files = $self->param('old_filename') || die "'file' is an obligatory parameter";
    my $suffix = $self->param('suffix') || die "'suffix' is an obligatory parameter";
    my $output_dir = $self->output_dir;
    my $job_name = $self->job_name;

    throw("will only rename one file") if ref($files) eq 'ARRAY' && scalar @$files >1;
    my $old_filename = ref($files) eq 'ARRAY' ? $files->[0] : $files;
    throw("do not have a file") if !$old_filename;

    my $new_filename = "$output_dir/$job_name.$suffix";
    $new_filename =~ s{//+}{/}g;
    $old_filename =~ s{//+}{/}g;
    throw("this module shouldn't be used if new_file and old_file are the same") if $new_filename eq $old_filename;

    check_directory_exists($output_dir);
    move($old_filename, $new_filename) or throw("could not move $old_filename to $new_filename $!");

    $self->output_this_branch('new_filename' => $new_filename);

}


1;

