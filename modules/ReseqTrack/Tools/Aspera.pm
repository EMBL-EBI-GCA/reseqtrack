package ReseqTrack::Tools::Aspera;

use strict;
use warnings;
use autodie;

use File::Path qw(make_path remove_tree);
use File::Basename qw(dirname basename);
use File::Temp qw(tempdir);
use ReseqTrack::Tools::Argument qw(rearrange);

use base ('ReseqTrack::Tools::RunProgram');


=head2 DEFAULT_OPTIONS

  Function  : Called by the RunProgram parent object in constructor
  Returntype: hashref
  Example   : my %options = %{&ReseqTrack::Tools:RunSamtools::DEFAULT_OPTIONS};

=cut

sub DEFAULT_OPTIONS {
    return {
        'ascp_exe'     => 'ascp',
        'ascp_param'   => {
            'k' => 2,       # Resume if checksums match, else re-transfer
            'Q' => undef,   # Queue the downloads fairly (is probably deprecated)
            'T' => undef,   # Don't use encryption
            'r' => undef,   # ? (is probably deprecated)
            'i' => '~/.aspera/cli/etc/asperaweb_id_dsa.openssh', # Default location for private key file
        },
    };
}


=head2 new

  Arg [-ascp_exe]     :
      string, path to aspera CLI
  Arg [-username]     :
      string, user account on the target server
  Arg [-aspera_url]   :
      string, the root location of the aspera server
  Arg [-ascp_param] :
      hash, key-value pairs of parameters for aspera. Parmameters without values are passed in with a value of undef

  + Arguments for ReseqTrack::Tools::RunProgram parent class

  Function  : Creates a new ReseqTrack::Tools::Aspera object.
  Returntype: ReseqTrack::Tools::Aspera
  Exceptions:
  Example   : my $run_aspera = ReseqTrack::Tools::Aspera->new(
                                -aspera_url  = 'fasp.sra.ebi.ac.uk',
                                -username    = 'era-fasp',
                                -ascp_param  = {'P' => 33001, 'k' => 2}
                              );
              $run_aspera->run_download(
                                -remote_path  = $source_path,
                                -local_path   = $output_path,
                              );

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);

    my ($ascp_exe, $username, $aspera_url, $ascp_param) = rearrange([
        qw(ASCP_EXE USERNAME ASPERA_URL ASCP_PARAM ) ],
        @args
    );

    $self->ascp_exe($ascp_exe // $self->options('ascp_exe'));
    $self->check_ascp_exe;  # Check that the aspera executable can be accessed

    $self->username($username);
    $self->aspera_url($aspera_url);

    $ascp_param = ($self->options('ascp_param'), $ascp_param);
    $self->ascp_param($ascp_param);

    return $self;
}


sub run_download {
    my ($self, @args) = @_;

    my ($remote_path, $local_path) = rearrange(
        [ qw(REMOTE_PATH LOCAL_PATH) ],
        @args
    );

    my $log_dir = get_logdir($local_path);
    my $ascp_param = $self->ascp_param;
    $$ascp_param{'L'} = $log_dir;


    # Assemble the Aspera command
    my $cmd = $self->ascp_exe;
    my $ascp_param_str = _get_hash_to_string($ascp_param);
    $cmd .= ' ' . $ascp_param_str if $ascp_param_str;
    my $asp_addr = $self->username . '@' . $self->aspera_url . ':';


    $local_path =~ s{//}{/}g; # Escape any spaces
    make_path(dirname $local_path); ## preserve the directory structure
    $cmd .= ' ' . $asp_addr . $remote_path . ' ' . $local_path;

    # Use the inherited function to execute the command
    $self->execute_command_line($cmd);

    # Cleanup the log file
    _check_log_file($log_dir);
    remove_tree($log_dir);
}


######################################
###        Utilities

sub check_ascp_exe {
    my $self = shift;
    my $ascp_exe = $self->ascp_exe;

    # Check that the Aspera executable is usable
    if ($ascp_exe =~ /\//) {
        die "$ascp_exe doesn't exist", $/ unless -e $ascp_exe;
        die "$ascp_exe not executable", $/ unless -x $ascp_exe;
    }
}

sub get_logdir {
    my ($local_path) = @_;

    my $local_dir = dirname $local_path;
    my $file_name = basename $local_path;
    $file_name =~ s{\\\s}{_}g; # Change spaces to underscores in the file name

    my $log_dir = tempdir($file_name . 'XXXX', DIR => $local_dir, CLEANUP => 0);
    return $log_dir
}


sub _get_hash_to_string {
    my ($ascp_param) = @_;
    my $str;

    foreach my $k (keys %{$ascp_param}) {
        my $v = $$ascp_param{$k} // undef;
        $str .= ' -' . $k;
        $str .= ' ' . $v if $v;
    }
    return $str;
}

sub _check_log_file {
    my ($log_dir) = @_;
    my $log_file = $log_dir . '/aspera-scp-transfer.log';
    open my $fh, '<', $log_file;
    while (<$fh>) {
        if (/LOG - Source file transfers passed\s+:\s+(\d)/) {
            die "file not transferred correctly", $/
                unless $1 > 0;
        }
    }
    close($fh);
}


#################################################
##    Getter / Setter

sub ascp_exe {
    my ($self, $ascp_exe) = @_;
    if (defined($ascp_exe)) {
        $self->{'ascp_exe'} = $ascp_exe;
    }
    return $self->{'ascp_exe'};
}

sub aspera_url {
    my ($self, $aspera_url) = @_;
    if (defined($aspera_url)) {
        $self->{'aspera_url'} = $aspera_url;
    }
    return $self->{'aspera_url'};
}

sub username {
    my ($self, $username) = @_;
    if (defined($username)) {
        $self->{'username'} = $username;
    }
    return $self->{'username'};
}

sub ascp_param {
    my ($self, $ascp_param) = @_;
    if (defined($ascp_param)) {
        $self->{'ascp_param'} = $ascp_param;
    }
    return $self->{'ascp_param'};
}

1;
