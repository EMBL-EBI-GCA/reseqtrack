package ReseqTrack::Tools::MetaDataUtils;

use strict;
use warnings;
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);

use Exporter;
use vars qw (@ISA  @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(
        create_directory_path );
        
        
=head2 create_directory_path

  Arg [1]   : A metadata object (ReseqTrack::Run,Experiment,Sample or Study)
  Arg [2]   : string, directory layout e.g. 'sample_name/archive_sequence', tokens matching method or attribute names in the meta data will be substituted with that method/attributes value. Search order is method, then attribute, then related objects (run -> (experiment,sample), experiment -> study if no match in sample). 
  Arg [3]   : string, base directory
  Function  : Creates a directory path, combining the base directory path and run_meta_info values
  Returntype: string
  Exceptions: throws if no matches are made 
  Example   : $dir = create_directory_path($run, 'population/sample_name', '/path/to/dir')

=cut

sub create_directory_path{
  my ($meta_data, $directory_layout, $base_directory) = @_;
  my $dir_path = $base_directory;
  my @layout_chunks =  split( /\//, $directory_layout);
  my $method_matches = 0;

  foreach my $layout_chunk (@layout_chunks) {
    my $value = _lookup_property($meta_data,$layout_chunk);
    $dir_path .= '/' if ($dir_path);
		
		if (defined $value){
			$method_matches++; 
			$dir_path .= $value;	
		}
		else{
			$dir_path .= $layout_chunk;
		}
  }

  throw "Directory layout ($directory_layout) did not call any meta data methods or match any attributes"
    if (@layout_chunks && ! $method_matches);

  $dir_path =~ s/\/\//\//g;
  $dir_path =~ s/\s+/_/g;
  return $dir_path;
}

sub _lookup_property {
	my ($meta_data,$property) = @_;
	throw("no property given") unless $property;
	throw("no meta data given") unless $meta_data;
	my $method = $meta_data->can($property);
	return &$method($meta_data) if ($method);
	
	my $attributes = $meta_data->attributes_hash;
	return $attributes->{$property}->attribute_value() if $attributes->{$property};
	
	if ($meta_data->isa('ReseqTrack::Run')){
		my $v = _lookup_property($meta_data->experiment(),$property);
		if (!defined $v){
			$v = _lookup_property($meta_data->sample(),$property);
		}
		return $v;
	}
	if ($meta_data->isa('ReseqTrack::Experiment')){
		return _lookup_property($meta_data->study(),$property);
	}
	if ($meta_data->isa('ReseqTrack::Sample')){
		return undef;
	}
	if ($meta_data->isa('ReseqTrack::Study')){
		return undef;
	}
	
	
}
