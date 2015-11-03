
package ReseqTrack::Hive::PipeSeed::ForeignFiles;

use strict;
use warnings;
use base ('ReseqTrack::Hive::PipeSeed::BasePipeSeed');

use File::Find qw(find);
use File::stat;
use DateTime::Format::MySQL;
use ReseqTrack::Tools::Exception qw(throw warning);
use File::Spec;
=pod

=head1 NAME

ReseqTrack::Hive::PipeSeed::ForeignFiles

=head1 SYNOPSIS

A module for seeding a file release pipeline with foreign files
Child class of ReseqTrack::Hive::PipeSeed::BasePipeSeed

=head1 Options
  
  ignore_foreign_paths: set this option to 1 if you want the module to search the dropbox for the file (i.e. ignore the directory location in the db)
                        set this option to 0 if you want to trust the dropbox directory location is correct

=cut

sub default_options {
  return {
    ignore_foreign_paths => 0,
    sticky_tbi => 1,
    sticky_bai => 1,
    exclude_type => {},
  };
}

sub create_seed_params {
  my ($self) = @_;
  my %options = (%{$self->default_options}, %{$self->options});

  throw('this module will only accept pipelines that work on the file table')
      if $self->table_name ne 'file';

  my $require_attributes = $options{'require_attributes'} || {};
  my $output_attributes  = ref($options{'output_attributes'}) eq 'ARRAY' ? $options{'output_attributes'}
                     : defined $options{'output_attributes'} ? [$options{'output_attributes'}]
                     : [];  

  my $db = $self->db;
  my $pipeline = $self->pipeline;

  my $exclude_type = $options{exclude_type} ? $options{exclude_type} : {};
  
  my $remote_hosts = $db->get_HostAdaptor->fetch_all_remote();
  my $fa = $db->get_FileAdaptor;
  my $psa = $db->get_PipelineSeedAdaptor;
  my @seed_params;
  foreach my $host (@$remote_hosts) {


    FILE:
    foreach my $file (@{$fa->fetch_by_host($host->dbID)}) {
      next FILE if $options{sticky_bai} && $file->name =~ /\.bai$/;
      next FILE if $options{sticky_tbi} && $file->name =~ /\.tbi$/;
      next FILE if exists $$exclude_type{$file->type};


      my $existing_ps = $psa->fetch_by_seed_and_pipeline($file, $pipeline);
      if (@$existing_ps) {
        next FILE if grep {$_->is_running} @$existing_ps;
        next FILE if grep {$_->is_complete} @$existing_ps;
        next FILE if grep {$_->is_futile} @$existing_ps;
      }

      my $dropbox_file;
      if ($options{'ignore_foreign_paths'}) {
        my $filename = $file->filename;
        my @dropbox_files;
        find( sub {push(@dropbox_files, $File::Find::name) if $_ eq $filename}, $host->dropbox_dir);
        throw("multiple files with the same name: @dropbox_files") if @dropbox_files >1;
        next FILE if !@dropbox_files;
        $dropbox_file = $dropbox_files[0];
      }
      else {
        $dropbox_file = File::Spec->rel2abs('./'.$file->name, $host->dropbox_dir);


        next FILE if ! -e $dropbox_file;
      }

      my $index_extension = ($options{sticky_tbi} && $dropbox_file =~ /\.vcf\.gz$/) ? '.tbi'
                    : ($options{sticky_bai} && $dropbox_file =~ /\.bam$/) ? '.bai'
                    : '';
      next FILE if $index_extension && ! -e $dropbox_file.$index_extension;
      my $index_file;
      if ($index_extension) {
        $index_file = $fa->fetch_by_name($file->name . $index_extension);
        next FILE if !$index_file;
      }


      my $st = stat($dropbox_file) or throw("could not stat $dropbox_file: $!");
      my $drop_box_ctime = $st->ctime;
      my $updated = DateTime::Format::MySQL->parse_datetime($file->updated)->set_time_zone('local')->epoch;
      if ($drop_box_ctime > $updated) {
        $updated = $drop_box_ctime;
      }

      my %dropbox_details = ('path' =>$dropbox_file, 'ctime' => $drop_box_ctime);
      my %db_details = (dbID => $file->dbID, updated => $file->updated);

      if ($index_extension) {
        my $index_st = stat($dropbox_file.$index_extension) or throw("could not stat $dropbox_file$index_extension: $!");
        my $index_ctime = $index_st->ctime;
        $dropbox_details{'index_ctime'} = $index_ctime;
        $dropbox_details{'index_ext'} = $index_extension;
        if ($index_ctime > $updated) {
          $updated = $index_ctime;
        }
        my $index_updated = DateTime::Format::MySQL->parse_datetime($index_file->updated)->set_time_zone('local')->epoch;
        $db_details{'index_updated'} = $index_file->updated;
        $db_details{'index_dbID'} = $index_file->dbID;
        if ($index_updated > $updated) {
          $updated = $index_updated;
        }
      }
      next FILE if grep {$_ > $updated} map {DateTime::Format::MySQL->parse_datetime($_->created)->set_time_zone('local')->epoch} @$existing_ps;  
     
      my %output_params = (file => {dropbox => \%dropbox_details, db => \%db_details});

      foreach my $attr_name (keys %$require_attributes) {
        my ($attribute) = grep {$_->attribute_name eq $attr_name} @{$file->attributes};
        next FILE if !$attribute;
        my $values = $require_attributes->{$attr_name};
        $values = ref($values) eq 'ARRAY' ? $values : [$values];
        next FILE if !grep {$attribute->attribute_value eq $_} @$values;
      }
      
      my $file_attributes = $file->attributes;
      ATTRIBUTE:
      foreach my $attribute_name (@$output_attributes) {
         my ($attribute) = grep {$_->attribute_name eq $attribute_name} @$file_attributes;
         next ATTRIBUTE if !$attribute;
         $output_params{$attribute_name} = $attribute->attribute_value;
      }

      push(@seed_params, [$file, \%output_params]);
    }
  }

  $self->seed_params(\@seed_params);
}

1;
