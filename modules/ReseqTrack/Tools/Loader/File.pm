package ReseqTrack::Tools::Loader::File;
use strict;
use warnings;

use ReseqTrack::Tools::Loader;
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::FileUtils qw(create_objects_from_path_list );
use ReseqTrack::Tools::FileUtils
  qw(create_history assign_type assign_type_by_filename );
use ReseqTrack::Tools::FileSystemUtils
  qw( get_lines_from_file run_md5 get_md5hash);
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Loader;
use Data::Dumper;
use File::Basename;

use vars qw(@ISA);

@ISA = qw(ReseqTrack::Tools::Loader);

sub new {
    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);
    my (
        $md5_file,        $hostname,    $die_for_problems,
        $update_existing, $store_new,   $assign_types,
        $do_md5,          $md5_program, $remote,
        $md5_hash,
      )
      = rearrange(
        [
            qw(
              MD5_FILE        HOSTNAME
              DIE_FOR_PROBLEMS    UPDATE_EXISTING     STORE_NEW
              ASSIGN_TYPES     DO_MD5
              MD5_PROGRAM REMOTE MD5_HASH
              )
        ],
        @args
      );

    #Defaults
    $self->assign_types('1');
    $self->md5_program("md5sum");
#####

    $self->assign_types($assign_types);
    $self->md5_program($md5_program);
    $self->md5_file($md5_file);
    $self->hostname($hostname);
    $self->die_for_problems($die_for_problems);
    $self->update_existing($update_existing);
    $self->store_new($store_new);
    $self->remote($remote);
    $self->do_md5($do_md5);
    $self->md5_hash($md5_hash);

    if ( ( $self->type eq "MUST_FIX" ) && !$self->assign_types ) {
        throw("Must give load_files.pl a file type with -type");
    }

    $self->md5_on_off_check();

    $self->assign_host_object();

    return $self;
}

###############################################################################

sub process_input {
    my ($self) = @_;

    if ( $self->file && scalar( @{ $self->file } ) ) {
        $self->add_files_from_cmd_line();
    }

    if ( $self->list_file ) {
        $self->add_files_from_list_file();
    }

    if ( $self->dir ) {
        $self->add_files_from_dir();
    }

    if ( $self->md5_file ) {
        $self->get_paths_from_md5_file();    # -md5 file
    }

    throw
"Found 0 files specifed using standard options: -file -dir -md5_file -list_file"
      unless ( $self->file_paths && scalar @{ $self->file_paths } );
    return;

}

sub sanity_check_objects {
    my $self  = shift;
    my $files = $self->file_paths;
    my %bad;

    throw "No files specifed" unless $files;

    print "Sanity check\n" if $self->verbose;

    foreach my $file (@$files) {

        if ( !-e $file ) {
            $bad{$file} = "Does not exist    :";
            print STDERR "Does not exist:$file\n";
            next;
        }

        if ( -d $file ) {
            $bad{$file} = "Is directory        :";
            print STDERR "Is directory:$file\n";
            next;
        }

        if ( !( $file =~ /^\// ) ) {
            print STDERR "Not full path $file\n";
            $bad{$file} = "Full path bad      :";
            next;
        }

    }

    #better be unique.
    my %unique;
    my %dups;
    my $tot_dups = 0;

    foreach my $file (@$files) {
        $dups{$file}++;
    }    #just get rid of duplication rather than complain
    my @paths = keys(%dups);
    $files = \@paths;
    if ( keys %bad && $self->hostname eq "1000genomes.ebi.ac.uk" )
    {    # don't throw when it is remote host loading
        warning "Found the following problems:\n";
        foreach my $i ( keys %bad ) {
            print STDERR $bad{$i}, $i, "\n";
        }
        throw "Fix file path problems";
    }

    print "Have " . @$files . " files to load\n" if $self->verbose;
    return;
}

####
sub get_paths_from_md5_file {
    my ( $self, $md5file ) = @_;
    $md5file = $self->md5_file unless ($md5file);

    my $hash = get_md5hash( $self->md5_file );

    my @paths;

    foreach my $key ( keys(%$hash) ) {
        my $md5       = $hash->{$key};
        my $full_path = $key;

        $hash->{$full_path} = $md5;
        push( @paths, $full_path );
    }

    $self->md5_hash($hash);
    $self->file_paths( \@paths );
    return;

}

sub md5_hash {
    my ( $self, $hash ) = @_;

    if ($hash) {
        throw( "Can't pass md5_hash a " . $hash . " it must be a hashref" )
          unless ( ref($hash) eq "HASH" );
        $self->{md5_hash} = $hash;
    }
    elsif ( $self->md5_file ) {
        $self->{md5_hash} = get_md5hash( $self->md5_file );
    }
    return $self->{md5_hash};
}

sub md5_on_off_check {
    my $self = shift;

    if ( !$self->{md5_file} && !$self->{do_md5} && !$self->{md5_hash} ) {

        warning
"do_md5  OFF. Not loading from an md5 list. No md5s for files. You sure??";

        print "Do you want to continue ? y or n\n";
        my $go_ahead = <STDIN>;
        chomp($go_ahead);
        $go_ahead = uc $go_ahead;

        if ( $go_ahead eq 'N' || $go_ahead eq "NO" ) {
            exit;
        }
    }
}

sub create_objects {
    my $self = shift;

    my $objects =
      create_objects_from_path_list( $self->file_paths, $self->type,
        $self->host );

    #my $objs = scalar(@$objects);
    #print "Created $objs file objects\n";

    if ( $self->assign_types && !$self->type ) {

        #print "Assigning types\n";
        $objects = assign_type( $objects, $self->db );
    }
    $self->objects($objects);
    return;
}

sub load_objects {
    my ($self)   = @_;
    my $md5_hash = $self->md5_hash;
    my $files    = $self->objects;

    if ( $self->do_md5 ) {

        if ( @$files >= 50 ) {
            warning("Running md5s for "
                  . @$files
                  . " files this may take a while" );
            print "Or use farm to calculate md5's and load via -md5_file\n";
        }
    }

    my $fa = $self->db->get_FileAdaptor;

    $md5_hash = $self->md5_hash;

    my @storage_problems;
  FILE: foreach my $file (@$files) {
        if ( $self->do_md5 ) {
            my $md5 = run_md5( $file->full_path, $self->md5_program );
            $file->md5($md5);
        }

        if ( $md5_hash && keys %$md5_hash ) {
            my $md5 = $md5_hash->{ $file->full_path };

            if ( $file->md5 ) {
                unless ( $file->md5 eq $md5 ) {
                    print STDERR $file->full_path
                      . " md5 mismatch.Skipping the file\n";
                    print STDERR "\tcalculated = "
                      . $file->md5 . "\n"
                      . "\tmd5hash    = $md5.  \n";
                    throw();
                }
            }
            $file->md5($md5);
        }

        if ( !$file->size && -e $file->full_path ) {
            my $size = -s $file->full_path;
            $file->size($size);
        }

        eval {
            if ( $self->update_existing )
            {
                my $existing = $fa->fetch_by_name( $file->name );
                if ($existing) {
                    $file->dbID( $existing->dbID );
                    my $history = create_history( $file, $existing );
                    $file->history($history) if ($history);
                    unless ($history) {
                        next FILE;
                    }
                }
                else {
                    my $possible_existing =
                      $fa->fetch_by_filename( $file->filename );
                    if ($possible_existing) {
                        if ( @$possible_existing == 1 ) {
                            my $existing = $possible_existing->[0];
                            $file->dbID( $existing->dbID );
                            my $history = create_history( $file, $existing );
                            next FILE unless ($history);
                            $file->history($history) if ($history);
                        }
                        elsif ( @$possible_existing >= 2 ) {
                            my $for_update;
                            foreach my $existing (@$possible_existing) {
                                if ( $existing->type eq $file->type ) {
                                    if ($for_update) {
                                        warning("Can't update "
                                              . $file->filename
                                              . " there are multiple files "
                                              . "which share its name and type"
                                        );
                                    }
                                    $for_update = $existing;
                                }
                            }
                            $file->dbID( $for_update->dbID ) if ($for_update);
                            my $history = create_history( $file, $for_update )
                              if ($for_update);
                            next FILE unless ($history);
                            $file->history($history) if ($history);
                            unless ($for_update) {
                                print STDERR "There are "
                                  . @$possible_existing
                                  . " possible existing files\n";

                                foreach my $file (@$possible_existing) {
                                    print STDERR $file->dbID . " "
                                      . $file->name . "\n";
                                }
                                throw(  "Have multiple files linked to "
                                      . $file->filename
                                      . " not sure how "
                                      . "to update the file" );
                            }
                        }
                    }
                }
            }
            $fa->store( $file, $self->update_existing, $self->store_new );
        };
        if ($@) {
            throw( "Problem storing " . $file . " " . $file->full_path . " $@" )
              if ( $self->die_for_problems );
            push( @storage_problems, $file->full_path . " " . $@ );
        }
    }

    foreach my $problem (@storage_problems) {
        print STDERR $problem . "\n";
    }

    return $files;

}

sub die_for_problems {
    my ( $self, $arg ) = @_;
    $self->{die_for_problems} = $arg if ($arg);
    return $self->{die_for_problems};
}

sub update_existing {
    my ( $self, $arg ) = @_;
    $self->{update_existing} = $arg if ($arg);
    return $self->{update_existing};
}

sub store_new {
    my ( $self, $arg ) = @_;
    $self->{store_new} = $arg if ($arg);
    return $self->{store_new};
}

sub do_md5 {
    my ( $self, $arg ) = @_;
    $self->{do_md5} = $arg if ($arg);
    return $self->{do_md5};
}

sub md5_program {
    my ( $self, $arg ) = @_;
    $self->{md5_program} = $arg if ($arg);
    return $self->{md5_program};
}

sub assign_types {
    my ( $self, $arg ) = @_;
    $self->{assign_types} = $arg if ( defined $arg );
    return $self->{assign_types};
}

sub remote {
    my ( $self, $arg ) = @_;
    $self->{remote} = $arg if ($arg);
    return $self->{remote};
}

sub assign_host_object {
    my $self = shift;

    my $host = get_host_object( $self->hostname, $self->db );

    if ( !$host ) {
        throw "Failed to create host object";
    }

    $self->host($host);
    return;
}

#######################
sub objects {
    my ( $self, $arg ) = @_;
    $self->{objects} = $arg if ($arg);
    return $self->{objects};
}

sub host {
    my ( $self, $arg ) = @_;
    $self->{host} = $arg if ($arg);
    return $self->{host};
}
###
sub hostname {
    my ( $self, $arg ) = @_;
    $self->{hostname} = $arg if ($arg);
    return $self->{hostname};
}
###
sub md5_file {
    my ( $self, $arg ) = @_;
    $self->{md5_file} = $arg if ($arg);
    return $self->{md5_file};
}
###
sub assign_md5_hash {
    my ( $self, $arg ) = @_;
    $self->{md5_hash} = $arg if ($arg);
    return $self->{md5_hash};
}
1;

