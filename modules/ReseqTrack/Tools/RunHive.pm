
=pod

=head1 NAME

ReseqTrack::Tools::RunHive

=head1 SYNOPSIS

This is a class for running hive scripts
It is a sub class of a ReseqTrack::Tools::RunProgram.

=head1 Example


=cut

package ReseqTrack::Tools::RunHive;
use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::FileSystemUtils qw(check_executable check_directory_exists);
use Cwd qw(getcwd);

use base qw(ReseqTrack::Tools::RunProgram);

=head2 DEFAULT_OPTIONS

  Function  : Called by the RunProgram parent object in constructor
  Returntype: hashref
  Example   : my %options = %{&ReseqTrack::Tools:RunPicard::DEFAULT_OPTIONS};

=cut

sub DEFAULT_OPTIONS {
    return {
      'loop' => 1,
      'use_log_dir' => 0,
    };
}


sub CMD_MAPPINGS { return {
    'run'   => \&run_beekeeper,
    'seed'    => \&run_seed_pipeline,
    'sync'    => \&run_sync_pipeline,
    'init'    => \&run_init_pipeline,
    };    
}

=head2 new

  Arg [-hive_scripts_dir]   :
      string, directory containing hive scripts files
  Arg [-url]   :
      string, the url of the hive database
  Arg [-hive_user]   :
      string, the username for access to the hive database
  Arg [-hive_password]   :
      string, the password for access to the hive database
  + Arguments for ReseqTrack::Tools::RunProgram parent class

  Function  : Creates a new ReseqTrack::Tools::RunHive object.
  Returntype: ReseqTrack::Tools::RunHive
  Exceptions: 
  Example   : my $run_hive = ReseqTrack::Tools::RunHive->new(
                -hive_scripts_dir => '/path/to/ensembl/ensembl-hive/scripts',
                -working_dir => '/path/to/dir/',
                -url => 'mysql://mysql-g1k:4175/my_db'
                -hive_user => 'me',
                -hive_password => 'my_password'
                )

=cut

sub new {
    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);

    my ( $hive_scripts_dir, $hive_host, $hive_port, $hive_dbname, $hive_user, $hive_password )
      = rearrange(
        [
            qw( HIVE_SCRIPTS_DIR HIVE_HOST HIVE_PORT HIVE_DBNAME HIVE_USER HIVE_PASSWORD
              )
        ],
        @args
      );

    $self->hive_scripts_dir( $hive_scripts_dir     || $self->program );
    $self->hive_port($hive_port);
    $self->hive_host($hive_host);
    $self->hive_dbname($hive_dbname);
    $self->hive_user($hive_user);
    $self->hive_password($hive_password);

    $self->working_dir(getcwd) if ! $self->working_dir;

    return $self;
}

=head2 run_program

  Arg [1]   : ReseqTrack::Tools::RunHive
  Arg [2]   : string, command, must be one of the following:
              'beekeeper', 'seed', 'sync', 'init'
  Function  : Runs hive scripts to manage hive pipelines
              This method is called by the RunProgram run method
  Returntype: 
  Exceptions: Throws if it cannot find the hive_scripts_dir. Throws if command is not recognised.
  Example   : $self->run('merge');

=cut

sub run_program {
    my ( $self, $command, @args ) = @_;

    throw("hive_scripts_dir is not a directory") if ( !-d $self->hive_scripts_dir );

    my $subs = CMD_MAPPINGS();

    throw("Did not recognise command $command")
      if ( !defined $subs->{$command} );

    throw("hive_user not defined") if !$self->hive_user;
    throw("hive_password not defined") if !$self->hive_password;
    throw("hive_port not defined") if !$self->hive_port;
    throw("hive_host not defined") if !$self->hive_host;
    throw("hive_dbname not defined") if !$self->hive_dbname;

    my @returned_values = &{ $subs->{$command} }($self, @args);
    return @returned_values;
}

sub get_valid_commands {
    my $subs = CMD_MAPPINGS();
    return keys %$subs;
}

sub make_url {
  my ($self) = @_;
  my $hive_user = $self->hive_user;
  my $hive_password = $self->hive_password;
  my $hive_dbname = $self->hive_dbname;
  my $hive_port = $self->hive_port;
  my $hive_host = $self->hive_host;

  return "mysql://$hive_user:$hive_password\@$hive_host:$hive_port/$hive_dbname";

}

sub run_sync_pipeline {
    my $self = shift;

    my $url = $self->make_url;

    my $script = $self->hive_scripts_dir . '/beekeeper.pl';

    my @cmd_words = ( 'perl', $script, '-url', $url, '-sync');
    my $cmd = join( ' ', @cmd_words );
    $self->execute_command_line($cmd);

}

sub run_beekeeper {
    my $self = shift;

    my $url = $self->make_url;

    my $script = $self->hive_scripts_dir . '/beekeeper.pl';

    my @cmd_words = ( 'perl', $script, '-url', $url);
    push( @cmd_words, ($self->options('loop') ? '-loop' : '-run'));

    if ($self->options('use_log_dir')) {
      my $log_dir = $self->working_dir . '/hive_log_dir';
      check_directory_exists($log_dir);
      push( @cmd_words, '-hive_log_dir', $log_dir);
    }
    if ($self->options('submission_options')){
      push( @cmd_words, '-submission_options',$self->options('submission_options'));
    }
    push(@cmd_words, '&>', '/dev/null');

    my $cmd = join( ' ', @cmd_words );
    $self->execute_command_line($cmd);

}

sub run_seed_pipeline {
    my ($self, $logic_name, $input_id) = @_;

    throw("need a logic name") if !defined $logic_name;
    throw("need an input_id") if !defined $input_id;

    my $url = $self->make_url;

    my $script = $self->hive_scripts_dir . '/seed_pipeline.pl';

    my @cmd_words = ( 'perl', $script, '-url', $url);
    push( @cmd_words, '-logic_name', $logic_name);
    push( @cmd_words, '-input_id', q(').$input_id.q('));

    my $cmd = join( ' ', @cmd_words );
    $self->execute_command_line($cmd);

}

sub run_init_pipeline {
    my ($self, $module, $init_options) = @_;

    throw("must have a module") if !$module;

    my $script = $self->hive_scripts_dir . '/init_pipeline.pl';

    my @cmd_words = ( sprintf('PATH=%s:$PATH', $self->hive_scripts_dir));
    push(@cmd_words, 'perl', $script, $module);
    push(@cmd_words, $init_options);
    push(@cmd_words, '-host', $self->hive_host);
    push(@cmd_words, '-password', $self->hive_password);
    push(@cmd_words, '-pipeline_db', '-user='.$self->hive_user);
    push(@cmd_words, '-pipeline_db', '-dbname='.$self->hive_dbname);
    push(@cmd_words, '-pipeline_db', '-port='.$self->hive_port);

    my $cmd = join( ' ', @cmd_words );
    $self->execute_command_line($cmd);

}

sub hive_scripts_dir {
    my $self = shift;
    return $self->program(@_);
}

sub hive_dbname {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{'hive_dbname'} = $arg;
    }
    return $self->{'hive_dbname'};
}

sub hive_host {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{'hive_host'} = $arg;
    }
    return $self->{'hive_host'};
}

sub hive_port {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{'hive_port'} = $arg;
    }
    return $self->{'hive_port'};
}

sub hive_user {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{'hive_user'} = $arg;
    }
    return $self->{'hive_user'};
}

sub hive_password {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{'hive_password'} = $arg;
    }
    return $self->{'hive_password'};
}

1;

