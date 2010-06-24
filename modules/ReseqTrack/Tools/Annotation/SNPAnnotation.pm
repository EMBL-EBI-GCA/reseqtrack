=pod

=head1 NAME

ReseqTrack::Tools::Annotation::SNPAnnotation;

=head1 SYNOPSIS

An object to create a file of annotations on a series of given snp calls

This module uses the ensembl-variation code to fetch snp consequences and 
both ensembl and vcftools to discover what sources particular snps can be
found in

=cut

package ReseqTrack::Tools::Annotation::SNPAnnotation;

use strict;
use warnings;
use Exporter;
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::GeneralUtils qw(create_filename);
use ReseqTrack::Tools::FileSystemUtils;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::VariationFeature;
use File::Basename;
use File::Path;
use FileHandle;


=head2 new

  Arg [-input_file] : path to input file
  Arg [-output_file] : path to an output file 
  Arg [-working_dir] : working directory for temp files
  Arg [-vcftools] : path to vcftools program
  Arg [-source_vcf] : arrayref pointing to vcf files to look in for matches
  Arg [-check_ensembl] : check ensembl database for which sources a snp can be
  found int
  Arg [-registry] : path to ensembl registry file
  Arg [-source] : the source of the snps, if not defined the name of the input file
  is used
  Arg [-strict_gvf] : follow gvf guidelines strictly, ie no source info
  Arg [-type] : what type of variation are these, defaults to SNV (single nucleotide
  variant)
  Arg [-species] : what species is the ensembl data base, defauls to human
  Arg [-buffer] : how many variation features to get consequences for at once 
  this defaults to 500
  Function  : create a ReseqTrack::Tools::Annotation::SNPAnnotation object
  Returntype: ReseqTrack::Tools::Annotation::SNPAnnotation
  Exceptions: throws if strict gvf is requested but source data is also requested as
  this two things aren't compatible
  Example   : my $annotate = ReseqTrack::Tools::Annotation::SNPAnnotation->new
  (
   -input_file => $input_file,
   -output_file => $output_file,
   -working_dir => $working_dir,
   -vcftools => $vcftools,
   -source_vcf => \@source_vcfs,
   -check_ensembl => $check_ensembl_sources,
   -registry => $registry_file,
   -source => $source,
   -strict_gvf => $strict_gvf,
   -type => $type,
   -species => $species,
  );

=cut



sub new {
  my ($class, @args) = @_;
  my $self ={};
  bless $self,$class;
  my ($input_file, $output_file, $working_dir, $vcftools, $source_vcf,
      $check_ensembl_sources, $ensembl_registry, $source, $strict_gvf, $type,
      $species, $buffer) =
	rearrange([qw(INPUT_FILE OUTPUT_FILE WORKING_DIR VCFTOOLS SOURCE_VCF
		      CHECK_ENSEMBL REGISTRY SOURCE STRICT_GVF TYPE SPECIES 
		      BUFFER)], @args);
  #setting defaults
  $self->vcftools("vcftools");
  $self->working_dir("/tmp");
  $self->species("human");
  $self->type('SNV');
  $self->buffer(500);
  ###
  $working_dir =~ s/\/$// if($working_dir);
  $self->input_file($input_file);
  $self->output_file($output_file);
  $self->working_dir($working_dir);
  $self->vcftools($vcftools);
  $self->source($source);
  $self->source(basename($self->input_file)) unless($self->source);
  $self->source_vcf($source_vcf);
  $self->type($type);
  $self->check_ensembl_sources($check_ensembl_sources);
  $self->ensembl_registry($ensembl_registry);
  $self->strict_gvf($strict_gvf);
  $self->species($species);
  $self->buffer($buffer);
  if($self->strict_gvf){
    throw("Can't associated given snp locations with sources from ensembl or ".
	  "vcf files as require strict gvf format")
      if($self->check_ensembl_sources || $self->source_vcf);
  }
  return $self;
}

#accessor methods

=head2 accessor methods

  Arg [1]   : ReseqTrack::Tools::Annotation::SNPAnnotation
  Arg [2]   : string/int/arrayref/hashref
  Function  : getter/setter method for a given variable
  Returntype: string/int
  Exceptions: n/a
  Example   : my $vcftools = $self->vcftools

=cut



sub vcftools{
  my ($self, $vcftools) = @_;
  if($vcftools){
    $self->{vcftools} = $vcftools;
  }
  return $self->{vcftools};
}

sub strict_gvf{
  my ($self, $strict_gvf) = @_;
  if($strict_gvf){
    $self->{strict_gvf} = $strict_gvf;
  }
  return $self->{strict_gvf};
}

sub source{
  my ($self, $source) = @_;
  if($source){
    $self->{source} = $source;
  }
  return $self->{source};
}

sub type{
  my ($self, $type) = @_;
  if($type){
    $self->{type} = $type;
  }
  return $self->{type};
}
sub species{
  my ($self, $species) = @_;
  if($species){
    $self->{species} = $species;
  }
  return $self->{species};
}
sub buffer{
  my ($self, $buffer) = @_;
  if($buffer){
    $self->{buffer} = $buffer;
  }
  return $self->{buffer};
}
sub source_vcf{
  my ($self, $source_vcf) = @_;
  if($source_vcf){
    $self->{source_vcf} = $source_vcf;
  }
  return $self->{source_vcf};
}

sub check_ensembl_sources{
  my ($self, $check_ensembl_sources) = @_;
  if($check_ensembl_sources){
    $self->{check_ensembl_sources} = $check_ensembl_sources;
  }
  return $self->{check_ensembl_sources};
}

sub input_lines{
  my ($self, $input_lines) = @_;
  if($input_lines){
    $self->{input_lines} = $input_lines;
  }
  return $self->{input_lines};
}

sub score_hash{
  my ($self, $score_hash) = @_;
  if($score_hash){
    $self->{score_hash} = $score_hash;
  }
  return $self->{score_hash};
}

sub total_reads{
  my ($self, $total_reads) = @_;
  if($total_reads){
    $self->{total_reads} = $total_reads;
  }
  return $self->{total_reads};
}

sub source_hash{
  my ($self, $source_hash) = @_;
  if($source_hash){
    $self->{source_hash} = $source_hash;
  }
  return $self->{source_hash};
}

=head2 input_file

  Arg [1]   : ReseqTrack::Tools::Annotation::SNPAnnotation
  Arg [2]   : string, file path
  Function  : container for input file path
  Returntype: string
  Exceptions: throws if file doesn't exist
  Example   : $self->input_file("/path/to/file");

=cut



sub input_file{
  my ($self, $input_file) = @_;
  if($input_file){
    throw($input_file." doesn't exist") unless(-e $input_file);
    $self->{input_file} = $input_file;
  }
  return $self->{input_file};
}


=head2 output_file

  Arg [1]   : ReseqTrack::Tools::Annotation::SNPAnnotation
  Arg [2]   : string, filepath
  Function  : store output file path and open file if not already openned
  Returntype: string, filepath
  Exceptions: throws if fails to open file
  Example   : my $output_file = $self->output_file;

=cut



sub output_file{
  my ($self, $output_file) = @_;
  if($output_file){
    $self->{output_file} = $output_file;
  }
  unless($self->output_fh){
    if($output_file){
      my $out_fh = new FileHandle;
      $out_fh->open(">".$output_file) or throw("Failed to open ".$output_file." $!");
      $self->output_fh($out_fh)
    }
  }
  return $self->{output_file};
}


=head2 output_fh

  Arg [1]   : ReseqTrack::Tools::Annotation::SNPAnnotation
  Arg [2]   : FileHandle object
  Function  : container for output file handle
  Returntype: FileHandle
  Exceptions: n/a
  Example   : my $fh = $self->output_fh();

=cut



sub output_fh{
  my ($self, $output_fh) = @_;
  if($output_fh){
    $self->{output_fh} = $output_fh;
  }
  return $self->{output_fh};
}


=head2 files_to_cleanup

  Arg [1]   : ReseqTrack::Tools::Annotation::SNPAnnotation
  Arg [2]   : string/arrayref
  Function  : to hold a list of files which need to be deleted
  Returntype: arrayref
  Exceptions: n/a
  Example   : foreach my $file(@{$self->files_to_cleanup){
               unlink $file;
              }

=cut


sub files_to_cleanup{
  my ($self, $file) = @_;
  unless($self->{'files_to_cleanup'}){
    $self->{'files_to_cleanup'} = [];
  }
  if($file){
    if(ref($file) eq 'ARRAY'){
      push(@{$self->{'files_to_cleanup'}}, @$file);
    }else{
      push(@{$self->{'files_to_cleanup'}}, $file);
    }
  }
  return $self->{'files_to_cleanup'};
}



=head2 working_dir

  Arg [1]   : ReseqTrack::Tools::Annotation::SNPAnnotation
  Arg [2]   : string, path to working dir
  Function  : container for working dir
  Returntype: string
  Exceptions: throw if not an directory or doesn't exist
  Example   : my $working_dir = $self->working_dir;

=cut


sub working_dir{
  my ($self, $working_dir) = @_;
  if($working_dir){
    $self->{working_dir} = $working_dir;
  }
  unless(-d $self->{working_dir}){
    mkpath($self->{working_dir});
    throw($self->{working_dir}." does not exist") unless(-d $self->{working_dir});
  }
  return $self->{working_dir};
}



=head2 variation_features

  Arg [1]   : ReseqTrack::Tools::Annotation::SNPAnnotation
  Arg [2]   : arrayref of Bio::EnsEMBL::Variation::VariationFeature objects
  Function  : container for variation feature objects
  Returntype: arrayref
  Exceptions: throws if not passed variation feature objects
  Example   : my $vfs = $self->variation_features;

=cut


sub variation_features{
  my ($self, $vfs) = @_;
  if($vfs){
    unless(ref($vfs) eq 'ARRAY'){
      throw("Must pass variations_features an arrayref of ".
	    "Bio::EnsEMBL::Variation::VariationFeature objects not ".$vfs);
    }
    my @vfs = @$vfs;
    unless($vfs[0]->isa("Bio::EnsEMBL::Variation::VariationFeature")){
      throw("Arrayref must contain Bio::EnsEMBL::Variation::VariationFeature ".
	    "objects not ".$vfs->[0]);
    }
    $self->{variation_features} = $vfs;
  }
  return $self->{variation_features};
}


=head2 ensembl_registry

  Arg [1]   : ReseqTrack::Tools::Annotation::SNPAnnotation
  Arg [2]   : Bio::EnsEMBL::Registry
  Function  : container for a registry object
  Returntype: 
  Exceptions: 
  Example   : 

=cut



sub ensembl_registry{
  my ($self, $registry) = @_;
  if($registry){
    $self->{registry} = $registry;
  }
  return $self->{registry};
}


=head2 get_XXX_adaptor

  Arg [1]   : ReseqTrack::Tools::Annotation::SNPAnnotation
  Function  : get a specific adaptor from the registry object
  Returntype: an ensembl adaptor object
  Exceptions: n/a
  Example   : my $vfa = $self->get_variation_feature_adaptor;

=cut



sub get_variation_feature_adaptor{
  my ($self) = @_;
  return
    $self->ensembl_registry->get_adaptor($self->species, 'variation', 
					 'variationfeature');
}

sub get_transcript_variation_adaptor{
  my ($self) = @_;
  return
    $self->ensembl_registry->get_adaptor($self->species, 'variation', 
					 'transcriptvariation');
}

sub get_slice_adaptor{
  my ($self) = @_;
  return
    $self->ensembl_registry->get_adaptor($self->species, 'core', 'slice');
}



=head2 slice_hash

  Arg [1]   : ReseqTrack::Tools::Annotation::SNPAnnotation
  Arg [2]   : hashref or Bio::EnsEMBL::Slice
  Function  : container for hashref of slice objects
  Returntype: hashref
  Exceptions: throws if it doesn't recognised the given arguement
  Example   : $self->slice_hash($slice);

=cut



sub slice_hash{
  my ($self, $arg) = @_;
  if(!$self->{slice_hash}){
    $self->{slice_hash} = {};
  }
  if($arg){
    if(ref($arg) eq 'HASH'){
      $self->{slice_hash} = $arg;
    }elsif($arg->isa("Bio::EnsEMBL::Slice")){
      $self->{slice_hash}->{$arg->seq_region_name} = $arg;
    }else{
      throw("Don't know how to deal with ".$arg);
    }
  }
  return $self->{slice_hash};
}


=head2 process_input

  Arg [1]   : ReseqTrack::Tools::Annotation::SNPAnnotation
  Function  : parse the input file
  Returntype: arrayref of lines, in the format needed to create VariationFeatures
  and two hashrefs, one mapping name to score and one mapping name to total read 
  depth, both stats are available from vcf files. These are also stored in the
  objects container methods for these variables
  Exceptions: throws if can't open the input file
  Example   : my ($input_lines, $score_hash, $total_reads) = $self->process_input;

=cut



sub process_input{
  my ($self) = @_;
  my $in;
  if($self->input_file =~ /\.gz/){
    open(IN, "gunzip -c ".$self->input_file." | ") or
      throw("Failed to open ".$self->input_file." $!");
    $in = \*IN;
  }else{
    open(IN, "<", $self->input_file) or
      throw("Failed to open ".$self->input_file." $!");
    $in = \*IN;
  }
  my @vfs;
  my %score_hash;
  my %total_reads;
  while(<$in>){
    chomp;
    next if(/^\#/);
    my ($string, $name, $score, $depth) = $self->parse_line($_);
    next unless($string);
    my $vf = $self->create_variation_feature($string);
    push(@vfs, $vf);
    unless($score_hash{$name}){
      $score_hash{$name} = $score;
    }else{
       #throw("Have none unique variant name for ".$self->input_file." ".$name.
#	    " ".$_);
    }
     unless($total_reads{$name}){
      $total_reads{$name} = $depth;
    }else{
      #throw("Have none unique variant name for ".$self->input_file." ".$name);
    }
  }
  close(IN);
  $self->variation_features(\@vfs);
  $self->score_hash(\%score_hash);
  $self->total_reads(\%total_reads);
  return ($self->variation_features, $self->score_hash, $self->total_reads);
}


=head2 calculate_consequences

  Arg [1]   : ReseqTrack::Tools::Annotation::SNPAnnotation
  Function  : uses the transcript variation adaptor to work out consequences for
  the given snps, as the consequences are attached to the variation features which
  go in we can just keep the processed features
  Returntype: n/a
  Exceptions: n/a
  Example   : $self->calculate_consequences

=cut


sub calculate_consequences{
  my ($self) = @_;
  #my @vfs = @{$self->create_variation_features};
  my $tva = $self->get_transcript_variation_adaptor;
  my @temp;
  my @new_vfs;
  my @vfs = @{$self->variation_features};
  while(@vfs){
    my $vf = shift @vfs;
    push(@temp, $vf);
    if(@temp >= $self->buffer){
      my $new = $tva->fetch_all_by_VariationFeatures(\@temp);
      push(@new_vfs, @temp);
      @temp = ();
    }
  }
  if(@temp >= 1){
    my $new = $tva->fetch_all_by_VariationFeatures(\@temp);
    push(@new_vfs, @temp);
  }
  $self->variation_features(\@new_vfs);
}


=head2 generate_unique_vcf_name

  Arg [1]   : ReseqTrack::Tools::Annotation::SNPAnnotation
  Arg [2]   : string, name
  Function  : the function generates a unique name for each line of the gvf file
  it is a recursive method as it calls itself if the name it creates is not unique
  it adds sequential numbers to the end of the name in the form name:001
  Returntype: string, name
  Exceptions: n/a
  Example   : my $unique_name = $self->generate_unique_gvf_name($name);

=cut



sub generate_unique_gvf_name{
  my ($self, $name) = @_;
  if($self->{'gvf_unique_name'}->{$name}){
    my @values = split /\:/, $name;
    unless($values[1]){
      $name .= ":001";
      return $self->generate_unique_gvf_name($name);
    }else{
      my $num =~ s/^0*//;
      $num++;
      $num = sprintf("%02d", $num);
      $name .= ":".$num;
      return $self->generate_unique_gvf_name($name);
    }
  }else{
    $self->{'gvf_unique_name'}->{$name} = 1;
    return $name;
  }
}


=head2 find_in_vcf

  Arg [1]   : ReseqTrack::Tools::Annotation::SNPAnnotation
  Function  : find the given snps in the given vcf files using the vcftools
  --positions filter
  Returntype: hashref of sources, the hashref key is a source and each value is a
  hash keyed on variation name
  Exceptions: throws if a gzipped vcf is provided or if the recode vcf file can't 
  be found or the command line fails
  Example   : my $source_hash = $self->find_in_vcf;

=cut



sub find_in_vcf{
  my ($self) = @_;
  my %hash;
  return undef unless($self->source_vcf && scalar(@{$self->source_vcf}) >= 1);
  #needs to generate appropriate positional file as input to vcf
  #also returns a map of position to name
  my ($tmp_file, $positions_hash) = $self->vcf_positional_input_generator;
  my %source_hash;
  foreach my $in_vcf(@{$self->source_vcf}){
    if($in_vcf =~ /\.gz/){
      throw("Can't give vcftools a gzipped vcf file ".$in_vcf);
    }
    my $basename = basename($in_vcf);
    my $out_stem = create_filename($self->working_dir, $self->source."_".$basename);
    my $vcf_cmd = $self->vcftools;
    $vcf_cmd .= " --vcf ".$in_vcf if($in_vcf =~ /\.vcf$/);
    $vcf_cmd .= " --positions ".$tmp_file." --recode --out ".$out_stem;
    my $exit;
    eval{
      $exit = system($vcf_cmd);
    };
    if($@ || $exit >= 1){
      throw("Failed to run ".$vcf_cmd." $exit $!");
    }
    my ($file_list, $hash) = list_files_in_dir($self->working_dir, 1);
    my $working_files = $hash->{$self->working_dir};
    my @out_files;
    my $vcf_file = $out_stem.".recode.vcf";
    #ensuring all the output files can be cleaned up
    foreach my $file(@$file_list){
      next unless(basename($file) =~ /$out_stem/);
      push(@out_files, $file);
      $self->files_to_cleanup($file);
    }
    unless(-e $vcf_file){
      print STDERR $vcf_file."\n";
      throw("Failed to find a vcf file from ".$self->working_dir." ".$vcf_cmd);
    }
    #uses position to name map to actually generate source hash
    open(VCF, "<", $vcf_file) or throw("Failed to open ".$vcf_file." $!");
    while(<VCF>){
      chomp;
      next if(/^\#/);
      my @values = split /\t/, $_;
      my $position = $values[0]."\t".$values[1];
      my $name = $positions_hash->{$position};
      throw("Failed to find a name for ".$position) unless($name);
      $source_hash{basename($in_vcf)}{$name} = 1;
    }
  }
  return(\%source_hash);
}


=head2 vcf_positional_input_generator

  Arg [1]   : ReseqTrack::Tools::Annotation::SNPAnnotation
  Arg [2]   : arrayref of Bio::EnsEMBL::Variation::VariationFeatures, if not passed
  in it will take what is in the variation_features array from the object
  Function  : to turn the variation features into a file suitable to be used 
  as a positional filter in vcftools and a hashref mapping position string to name
  Returntype: file and hashref of position to name
  Exceptions: throws if it can't open the tmp file
  Example   : my ($position_file, $position_hash) = 
                    $self->vcf_positional_input_generator()

=cut



sub vcf_positional_input_generator{
  my ($self, $vfs) = @_;
  $vfs = $self->variation_features unless($vfs);
  my $tmp_file = create_filename($self->working_dir, "vcf_positional".$self->source, 'out');
  open(FH, ">", $tmp_file) or throw("Failed to open ".$tmp_file." $!");
  my %position_hash;
  foreach my $vf(@$vfs){
    my $string = $vf->seq_region_name."\t".$vf->start;
    print FH $string."\n";
    unless($position_hash{$string}){
      $position_hash{$string} = $vf->variation_name;
    }else{
      throw($string." is not unique it is associated with ".$position_hash{$string}.
	    " as well as ".$vf->variation_name);
    }
  }
  close(FH);
  $self->files_to_cleanup($tmp_file);
  return $tmp_file, \%position_hash;
}


=head2 find_in_ensembl

  Arg [1]   : ReseqTrack::Tools::Annotation::SNPAnnotation
  Arg [2]   : arrayref of Bio::EnsEMBL::Variation::VariationFeatures, if not
  used the object array is used
  Function  : find what sources in the ensembl variation database your
  snps can be found in
  Returntype: hashref key source, value hashref for each name which is found in that
  source
  Exceptions: n/a
  Example   : my $source_hash = $self->find_in_ensembl();

=cut



sub find_in_ensembl{
  my ($self, $vfs) = @_;
  my %source_hash;
  $vfs = $self->variation_features unless($vfs);
  foreach my $vf(@$vfs){
    my $ensembl_sources = $vf->get_all_sources if($vf->variation);
    foreach my $source(@$ensembl_sources){
      $source_hash{$source}{$vf->variation_name} = 1;
    }
  }
  return \%source_hash;
}


=head2 produce_output

  Arg [1]   : ReseqTrack::Tools::Annotation::SNPAnnotation
  Function  : define consequences for the variation feature and create source hash 
  if requested
  Returntype: n/a
  Exceptions: n/a
  Example   : n/a

=cut



sub produce_output{
  my ($self) = @_;
  $self->calculate_consequences();
  my $source_hash;
  if($self->source_vcf){
    $source_hash = $self->find_in_vcf();
  }
  if($self->check_ensembl_sources){
    my $tmp_hash = $self->find_in_ensembl();
    foreach my $key(keys(%$tmp_hash)){
      $source_hash->{$key} = $tmp_hash->{$key};
    }
  }
  $self->source_hash($source_hash);
}


=head2 print_output

  Arg [1]   : ReseqTrack::Tools::Annotation::SNPAnnotation
  Function  : prints the consequence and source information in gvf format
  Returntype: n/a
  Exceptions: n/a
  Example   : $self->print_output

=cut



sub print_output{
  my ($self) = @_;
  my $fh = $self->output_fh;
  my $source_hash = $self->source_hash;
  my $score_hash = $self->score_hash;
 # print STDERR "Printing vf info\n";
  my $count = 1;
  foreach my $vf(@{$self->variation_features}){
 #   print "Looking at ".$count." variation feature\n";
    $count++;
    my $chr = $vf->seq_region_name;
    my $start = $vf->seq_region_start;
    my $end = $vf->seq_region_end;
    my $strand = $vf->seq_region_strand;
    my $name = $vf->variation_name;
    my $score = $score_hash->{$name};
    $score = "." unless($score);
    my $gff3 = join("\t", $chr, $self->source, $self->type, $start, $end,
		$score, $strand, ".");
    my $attrib_strings = $self->create_gvf_attribute_strings($vf);
    my @sources;
    foreach my $source(keys(%$source_hash)){
      if($source_hash->{$source}->{$vf->variation_name}){
	push(@sources, $source);
      }
    }
    foreach my $attrib_string(@$attrib_strings){
      my $output_string = $gff3."\t".$attrib_string."\t";
      $output_string .= join("\t", @sources) unless($self->strict_gvf);
      print $fh $output_string."\n";
    }
  }
}


=head2 create_gvf_attribute_string

  Arg [1]   : ReseqTrack::Tools::Annotation::SNPAnnotation
  Arg [2]   : Bio::EnsEMBL::Variation::VariationFeature
  Function  : the 9th column of gvf format is a series of attributes in key=value
  pairs. This method creates such a string for every consequence a given variant
  has
  Returntype: string
  Exceptions: n/a
  Example   : my $attributes = $self->create_gvf_attribute_strings($vf);

=cut



sub create_gvf_attribute_strings{
  my ($self, $vf) = @_;
  my @strings;

  my $unique_name = $self->generate_unique_gvf_name($vf->variation_name);
  my $total_reads = $self->total_reads->{$vf->variation_name};
  my %feature_type_hash = $self->feature_type;
  #chr1 SOAP SNV 15883 15883 36.5 + . ID=chr1:SOAP:SNV:15883;Variant_seq=G,C;Genotype=heterozygous; Variant_reads=17,16;Total_reads=33;Reference_seq=C;Variant_effect= nonconservative_missense_codon:mRNA:RefSeq:NM_012345,reference;
  my $attribute_string = "ID=".$unique_name.";";
  my $allele_string = $vf->allele_string;
  $allele_string =~ s/\//,/;
  $attribute_string .= "Variant_seq=".$allele_string.";";
  $attribute_string .= "Total_reads=".$total_reads.";"
    if($total_reads && $total_reads =~ /^\d+$/);
  foreach my $con(@{$vf->get_all_TranscriptVariations}){
    foreach my $consequence(@{$con->consequence_type}){
      if(($con->cdna_start && $con->cdna_end) && 
	 ($con->cdna_start > $con->cdna_end)) {
	($con->{'cdna_start'}, $con->{'cdna_end'}) = ($con->{'cdna_end'}, $con->{'cdna_start'});
      }
      if(($con->translation_start && $con->translation_end) &&
	 ($con->translation_start > $con->translation_end)) {
	($con->{'translation_start'}, $con->{'translation_end'}) = ($con->{'translation_end'}, $con->{'translation_start'});
      }
      my $feature_type = $feature_type_hash{$consequence};
      my $effect_string = "Variant_effect=".$consequence.":".$feature_type;
      if($con->transcript){
	my $transcript_name = $con->transcript->stable_id;
	$effect_string .= ":ensembl:".$transcript_name.";";
	
      }else{
	$effect_string .= ";";
      }
      my $string = $attribute_string.$effect_string;
      push(@strings, $string);
    }
  }
  return \@strings;
}



=head2 feature_type

  Function  : This simple returns a hard coded hash mapping consequence type to
  feature type
  Returntype: hashref 
  Exceptions: n/a
  Example   : n/a

=cut



sub feature_type{
  return ('ESSENTIAL_SPLICE_SITE' => 'mRNA',
	  'STOP_GAINED' => 'mRNA',
	  'STOP_LOST' => 'mRNA',
	  'COMPLEX_INDEL' => 'mRNA',
	  'FRAMESHIFT_CODING' => 'mRNA',
	  'NON_SYNONYMOUS_CODING' => 'mRNA',
	  'SPLICE_SITE' => 'mRNA',
	  'PARTIAL_CODON' => 'mRNA',
	  'SYNONYMOUS_CODING' => 'mRNA',
	  'REGULATORY_REGION' => 'genomic',
	  'WITHIN_MATURE_miRNA' => 'miRNA',
	  '5PRIME_UTR' => 'mRNA',
	  '3PRIME_UTR' => 'mRNA',
	  'INTRONIC' => 'mRNA',
	  'NMD_TRANSCRIPT' => 'mRNA',
	  'UPSTREAM' => 'mRNA',
	  'DOWNSTREAM' => 'mRNA',
	  'WITHIN_NON_CODING_GENE' => 'ncRNA',
	  'NO_CONSEQUENCE' => 'genomic',
	  'INTERGENIC' => 'genomic',
	  'HGMD_MUTATION' => 'mRNA',
	 );
}


=head2 find_colocated_variation_features

  Arg [1]   : ReseqTrack::Tools::Annotation::SNPAnnotation
  Arg [2]   : Bio::EnsEMBL::Variation::VariationFeature
  Function  : finds an existing variation features which has exactly the
  same coordinates as the passed in feature if present. This means that
  existing info can be attached to the give feature like sources or names
  Returntype: name, internal dbID, Bio::EnsEMBL::Variation::Variation, string
  Exceptions: n/a
  Example   : n/a

=cut



sub find_colocated_variation_feature{
  my ($self, $vf) = @_;
  my $existing;
  if($vf->adaptor->db){
    my $fs = $vf->feature_Slice;
    if($fs->start > $fs->end) {
      ($fs->{'start'}, $fs->{'end'}) = ($fs->{'end'}, $fs->{'start'});
    }
  EXISTING:foreach my $existing_vf_obj(@{$vf->adaptor->fetch_all_by_Slice($fs)}) {
      if($existing_vf_obj->seq_region_start == $vf->seq_region_start and 
	 $existing_vf_obj->seq_region_end == $vf->seq_region_end){
	$existing = $existing_vf_obj;
	last EXISTING; #this works on the assumption that a properly constructed
        #variation database will only ever have one variation which matches a 
	#specific start and end point
      }
    }
  }
  if($existing){
    return $existing->variation_name, $existing->dbID, $existing->variation, $existing->source ;
  }else{
    return undef;
  }
}


=head2 create_variation_features

  Arg [1]   : ReseqTrack::Tools::Annotation::SNPAnnotation
  Function  : creates variation features given the input lines in the tab delimited 
  format with the order chr start end allele_string strand name
  Returntype: Bio::EnsEMBL::Variation::VariationFeature
  Exceptions: warns if it can't find a suitable slice in the give core database for
  the chromosome name specified, this line is then skipped
  Example   :

=cut


sub create_variation_feature{
  my ($self, $line) = @_;
  my $sa = $self->get_slice_adaptor;
  $sa->dbc->disconnect_when_inactive(1);
  throw("Failed to connect to ensembl database ") unless($sa);
  my $vfa = $self->get_variation_feature_adaptor;
  my ($chr, $start, $end, $allele_string, $strand, $var_name) = split /\t/, $line;
    $chr =~ s/chr//ig;
    $strand = ($strand =~ /\-/ ? "-1" : "1");
    unless($start =~ /^\d+$/ && $end =~ /^\d+$/) {
      warn("WARNING: Start $start or end $end coordinate invalid on line".
	  $line);
      next;
    }

    unless($allele_string =~ /([ACGT-]+\/*)+/) {
      warn("WARNING: Invalid allele string $allele_string on line ".
	  $line);
      next;
    }
  my $slice;
  if($self->slice_hash->{$chr}){
    $slice = $self->slice_hash->{$chr};
    }else{
      $slice = $sa->fetch_by_region('chromosome', $chr);
      unless($slice){
	$slice = $sa->fetch_by_region(undef, $chr);
      }
      unless($slice){
	warning("Failed to find slice for ".$chr." in ".$sa->dbc->dbname);
	print STDERR "Skipping ".$line."\n";
	next LINE;
      }
      $self->slice_hash($slice);
    }
  my $vf = Bio::EnsEMBL::Variation::VariationFeature->new
    (
     -start => $start,
     -end => $end,
     -slice => $slice,           # the variation must be attached to a slice
     -allele_string => $allele_string,
     -strand => $strand,
     -map_weight => 1,
     -adaptor => $vfa,           # we must attach a variation feature adaptor
     -variation_name => $var_name,
    );
  my ($name, $dbID, $variation, $source) = $self->find_colocated_variation_feature($vf);
  $vf->dbID($dbID);
  $vf->variation($variation) if($variation);
  $vf->source($source);
  return $vf;
}


sub create_variation_features{
  my ($self, $input_lines) = @_;
  my $sa = $self->get_slice_adaptor;
  throw("Failed to connect to ensembl database ") unless($sa);
  my %slices;
  my @vfs;
  my $vfa = $self->get_variation_feature_adaptor;
  $input_lines = $self->input_lines unless($input_lines && @$input_lines >= 1);
  LINE:foreach my $input(@$input_lines){
    my ($chr, $start, $end, $allele_string, $strand, $var_name) = split /\t/, $input;
    $chr =~ s/chr//ig;
    $strand = ($strand =~ /\-/ ? "-1" : "1");
    unless($start =~ /^\d+$/ && $end =~ /^\d+$/) {
      warn("WARNING: Start $start or end $end coordinate invalid on line".
	  $input);
      next;
    }

    unless($allele_string =~ /([ACGT-]+\/*)+/) {
      warn("WARNING: Invalid allele string $allele_string on line ".
	  $input);
      next;
    }
    my $slice;
    if($slices{$chr}){
      $slice = $slices{$chr};
    }else{
      $slice = $sa->fetch_by_region('chromosome', $chr);
      unless($slice){
	$slice = $sa->fetch_by_region(undef, $chr);
      }
      unless($slice){
	warning("Failed to find slice for ".$chr." in ".$sa->dbc->dbname);
	print STDERR "Skipping ".$input."\n";
	next LINE;
      }
      $slices{$chr} = $slice;
    }
    my $vf = Bio::EnsEMBL::Variation::VariationFeature->new
      (
       -start => $start,
       -end => $end,
       -slice => $slice,           # the variation must be attached to a slice
       -allele_string => $allele_string,
       -strand => $strand,
       -map_weight => 1,
       -adaptor => $vfa,           # we must attach a variation feature adaptor
       -variation_name => $var_name,
      );
    my ($name, $dbID, $variation, $source) = $self->find_colocated_variation_feature($vf);
    $vf->dbID($dbID);
    $vf->variation($variation) if($variation);
    $vf->source($source);
    push(@vfs, $vf);
  }
  $self->variation_features(\@vfs);
  return \@vfs;
}


=head2 parse_line

  Arg [1]   : ReseqTrack::Tools::Annotation::SNPAnnotation
  Arg [2]   : string, line from input file
  Function  : parse line and produce format needed to create variation feature
  Returntype: string, name, int (score), int (total depth of reads)
  Exceptions: n/a
  Example   : n/a

=cut



sub parse_line{
  my ($self, $line) = @_;
  my @data = (split /\s+/, $_);
  my @output;
  my $name;
  my $score = ".";
  my $total_reads = ".";
  # pileup: chr1 60 T A
  next if(/^\#/);
  if($data[0] =~ /(chr)?\w+/ && $data[1] =~ /\d+/ && $data[2] =~ /[ACGTN-]+/ && $data[3] =~ /[ACGTN-]+/) {
    $name = $data[4];
    my $allele_string = $data[2]."/".$data[3];
    if(!$name || $name eq '.'){
      $name = $self->create_variant_name($data[0], $data[1], $data[1], 
					 $allele_string);
    }
    @output =  ($data[0], $data[1], $data[1], $allele_string, 1, $name);
  }
  
  # VCF: 20      14370   rs6054257 G     A      29    0       NS=58;DP=258;AF=0.786;DB;H2          GT:GQ:DP:HQ
  elsif($data[0] =~ /(chr)?\w+/ && $data[1] =~ /\d+/ && $data[3] =~ /[ACGTN-]+/ && $data[4] =~ /([\.ACGTN-]+\,?)+/) {
    if($data[4] eq '.'){
      return undef;
    }
    my $alt = $data[4];
    $alt =~ s/\,/\//g;
    my $allele_string = $data[3]."/".$alt;
    $name = $data[2];
    if(!$name || $name eq '.'){
      $name = $self->create_variant_name($data[0], $data[1], $data[1], 
					 $allele_string);
    }
    @output = ($data[0], $data[1], $data[1], $allele_string, 1, $name);
    $score = $data[5];
    my @info_elements = split /\;/, $data[7];
  INFO:foreach my $info(@info_elements){
      my ($type, $value) = split /\=/, $info;
      if($type eq 'DP'){
	$total_reads = $value;
	last INFO;
      }
    }
  }
  
  # our format
  else {
    # we allow commas as delimiter so re-split
    @data = (split /\s+|\,/, $_);
    unless($data[5]){
      $data[5] = $self->create_variant_name($data[0], $data[1], $data[2], $data[3]);
    }
    $name = $data[5];
    @output = @data;
  }
  my $string = join("\t", @output);
  return ($string, $name, $score, $total_reads);
}


=head2 create_variant_name

  Arg [1]   : ReseqTrack::Tools::Annotation::SNPAnnotation
  Arg [2]   : string, chromosome name
  Arg [3]   : int, start
  Arg [4]   : int, end
  Arg [5]   : string, allele string N/N
  Function  : create a name from the given values
  Returntype: string, name
  Exceptions: n/a
  Example   : n/a

=cut


sub create_variant_name{
  my ($self, $seq_name, $start, $end, $allele_string) = @_;
  my $name = join("_", ($seq_name, $start, $end, $allele_string));
  return $name
}

sub cleanup_output{
  my ($self) = @_;
  foreach my $file(@{$self->files_to_cleanup}){
    unlink $file;
  }
}

1;


