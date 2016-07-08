
package AccessibleGenome::MoveFile;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists check_directory_exists check_file_does_not_exist);
use File::Copy qw( move );
use File::Basename qw( dirname );
use Cwd qw( abs_path);


sub param_defaults {
  return {
    name_file_module => 'ReseqTrack::Hive::NameFile::BaseNameFile',
    name_file_method => 'default',
    name_file_params => {},
    clobber => 0,
  };
}


sub run {
    my $self = shift @_;

    my $current_file_path = $self->param_required('file');
    my $clobber = $self->param('clobber') ? 1 : 0;

    check_file_exists($current_file_path);

    my $name_file_module = $self->param('name_file_module') // param_defaults()->{'name_file_module'};
    eval "require $name_file_module" or throw "cannot load module $name_file_module $@";
    my $renamer = $name_file_module->new(
          -method => $self->param('name_file_method') // param_defaults()->{'name_file_method'},
          -params => $self->param('name_file_params'),
          -files => [$current_file_path],
          );
    $renamer->derive_paths;
    my $file_map = $renamer->file_map;

    my $new_path = $file_map->{$current_file_path};
    if($new_path) {;
      my $dir_name = dirname($new_path);
      check_directory_exists($dir_name);
      if (!$clobber && (abs_path($new_path) ne abs_path($current_file_path))) {
        check_file_does_not_exist($new_path);
      }
      move($current_file_path, $new_path) or throw("could not move $current_file_path to $new_path $!");
      $current_file_path = $new_path
    }

    $self->output_param('file', $current_file_path);

}

1;

