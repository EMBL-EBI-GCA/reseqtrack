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
            'Q' => undef,   # Queue the downloads fairly
            'T' => undef,   # Don't use encryption
            'r' => undef,   # ?
            'i' => '~/.aspera/cli/etc/asperaweb_id_dsa.openssh', # Default location for private key file
        },
    };
}


=head2 new

  Arg [-ascp_exe]     :
      string, path to aspera CLI
  Arg [-ascp_param_k] :
      number, the resume level
  Arg [-aspera_url]    :
      string, the root location of the aspera server

  + Arguments for ReseqTrack::Tools::RunProgram parent class

  Function  : Creates a new ReseqTrack::Tools::Aspera object.
  Returntype: ReseqTrack::Tools::Aspera
  Exceptions:
  Example   : my $run_aspera = ReseqTrack::Tools::Aspera->new(
                                -aspera_url=fasp.sra.ebi.ac.uk
                               );

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);

    my ($ascp_exe, $username, $aspera_url,
        $ascp_param, $work_dir,
        $filename, $download_dir, $upload_dir, $trim_path) = rearrange(
        [ qw(ASCP_EXE USERNAME ASPERA_URL
            ASCP_PARAM WORK_DIR
            FILENAME DOWNLOAD_DIR UPLOAD_DIR TRIM_PATH) ],
        @args
    );

    $self->ascp_exe($ascp_exe // $self->options('ascp_exe'));
    $self->username($username);
    $self->aspera_url($aspera_url);

    # merge the parameter hashes (the user's params take priority)
    # This solution doesn't give thr user the ability to disable flags,
    # will need to redesign if this feature is needed in future
    $ascp_param = ($self->options('ascp_param'), $ascp_param);
    $self->ascp_param($ascp_param);

    $self->filename($filename);
    $self->download_dir($download_dir);
    $self->upload_dir($upload_dir);
    $self->trim_path($trim_path);

    # Check that the aspera executable can be accessed
    $self->check_ascp_exe;

    # Check that the user option is either upload or download
    die "require either download_dir or upload_dir", $/ if !$download_dir && !$upload_dir;
    die "mutually exclusive options:  download_dir or upload_dir", $/ if $download_dir && $upload_dir;
    die "trim_path is required for file upload", $/ if $upload_dir && !$trim_path;

    # Set the workdir
    if (!defined($work_dir)){
        if (defined ($download_dir)) {
            $work_dir = dirname $download_dir;
        } else {
            $work_dir = dirname $filename;
        }
    }
    $self->work_dir($work_dir);

    return $self;
}



sub run_program {
    my $self = shift;

    my $filename = $self->filename;
    my $download_dir = $self->download_dir;
    my $upload_dir = $self->upload_dir;

    my $ascp_param = $self->ascp_param;
    $$ascp_param{'L'} = $self->get_logdir;    ## log dir path is required
    $$ascp_param{'d'} = undef if $upload_dir; ## required for creating directory in remote


    # Assemble the Aspera command
    my $cmd = $self->ascp_exe;
    my $ascp_param_str = _get_hash_to_string($ascp_param);
    $cmd .= ' ' . $ascp_param_str if $ascp_param_str;
    my $asp_addr = $self->username . '@' . $self->aspera_url . ':';

    my $trim_path = $self->trim_path;
    my $trim_re = qr/$trim_path/ if $trim_path;

    ### Run Aspera for a Download
    if ($download_dir) {
        my $dest_path = dirname $filename;
        $dest_path =~ s/$trim_re//g if $trim_path; ## trim destination path
        my $download_path = $download_dir . '/' . $dest_path;
        $download_path =~ s{//}{/}g;
        make_path($download_path); ## preserve the directory structure

        $cmd .= ' ' . $asp_addr . $filename . ' ' . $download_path;
    }
    ### Run Aspera for an Upload
    else {
        die "trim_path is required for file upload", $/ unless $trim_path;

        my $upload_path = dirname $filename;         ## remove filename from upload path
        $upload_path =~ s/$trim_re//g if $trim_path; ## trim upload path
        $upload_path = $upload_dir . '/' . $upload_path . '/';
        $upload_path =~ s{//}{/}g;

        $cmd .= ' ' . $filename . ' ' . $asp_addr . $upload_path
    }


    # Use the inherited function to execute the command
    $self->execute_command_line($cmd);

    # Cleanup the log file
    _check_log_file($self->get_logdir);
    remove_tree($self->get_logdir);
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
    my $self = shift;

    my $file_basename = basename($self->filename);
    $file_basename =~ s{\\\s}{_}g; # Change spaces to underscores in the file name

    my $log_dir = tempdir($file_basename . 'XXXX', DIR => $self->work_dir, CLEANUP => 0);
    print "LD: $log_dir , WD: ". $self->work_dir." \n";
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


sub work_dir {
    my ($self, $work_dir) = @_;
    if (defined($work_dir)) {
        $self->{'work_dir'} = $work_dir;
    }
    return $self->{'work_dir'}
}

sub trim_path {
    my ($self, $trim_path) = @_;
    if (defined($trim_path)) {
        $self->{'trim_path'} = $trim_path;
    }
    return $self->{'trim_path'}
}

sub filename {
    my ($self, $filename) = @_;
    if (defined($filename)) {
        $filename =~ s{\s}{\\ }g; # escape spaces in file path
        $self->{'filename'} = $filename;
    }
    return $self->{'filename'}
}

sub download_dir {
    my ($self, $download_dir) = @_;
    if (defined($download_dir)) {
        $self->{'download_dir'} = $download_dir;
    }
    return $self->{'download_dir'}
}

sub upload_dir {
    my ($self, $upload_dir) = @_;
    if (defined($upload_dir)) {
        $self->{'upload_dir'} = $upload_dir;
    }
    return $self->{'upload_dir'}
}

1;
