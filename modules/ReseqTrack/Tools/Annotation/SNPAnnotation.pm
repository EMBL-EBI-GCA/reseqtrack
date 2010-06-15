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
  Function  : 
  Returntype: 
  Exceptions: 
  Example   : 

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
  $self->working_dir("/tmp/");
  $self->species("human");
  $self->type('SNV');
  $self->buffer(500);
  ###

  $self->input_file($input_file);
  $self->output_file($output_file);
  $self->working_dir($working_dir);
  $self->vcftools($vcftools);
  $self->source($source);
  $self->source(basename($self->input_file)) unless($self->source);
  $self->source_vcf($source_vcf);
  $self->type($type);
  $self->check_ensembl_sources($check_ensembl_sources);
  $self->ensembl_registry_file($ensembl_registry);
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

sub input_file{
  my ($self, $input_file) = @_;
  if($input_file){
    throw($input_file." doesn't exist") unless(-e $input_file);
    $self->{input_file} = $input_file;
  }
  return $self->{input_file};
}

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

sub output_fh{
  my ($self, $output_fh) = @_;
  if($output_fh){
    $self->{output_fh} = $output_fh;
  }
  return $self->{output_fh};
}

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

sub ensembl_registry_file{
  my ($self, $ensembl_registry) = @_;
  if($ensembl_registry){
    throw($ensembl_registry." does not exist") unless(-e $ensembl_registry);
    $self->{ensembl_registry} = $ensembl_registry;
  }
  return $self->{ensembl_registry};
}

sub ensembl_registry{
  my ($self, $registry) = @_;
  if($registry){
    $self->{registry} = $registry;
  }
  unless($self->{registry}){
    my $reg = 'Bio::EnsEMBL::Registry';
    $reg->load_all($self->ensembl_registry_file);
    $self->{registry} = $reg;
  }
  return $self->{registry};
}

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

sub get_gene_adaptor{
  my ($self) = @_;
  return
    $self->ensembl_registry->get_adaptor($self->species, 'core', 
					 'gene');
}

sub process_input{
  my ($self) = @_;
  open(IN, "<", $self->input_file) or 
    throw("Failed to open ".$self->input_file." $!");
  my @input_lines;
  my %score_hash;
  my %total_reads;
  while(<IN>){
    chomp;
    next if(/^\#/);
    my ($string, $name, $score, $depth) = $self->parse_line($_);
    push(@input_lines, $string);
    unless($score_hash{$name}){
      $score_hash{$name} = $score;
    }else{
      throw("Have none unique variant name for ".$self->input_file." ".$name.
	    " ".$_);
    }
     unless($total_reads{$name}){
      $total_reads{$name} = $depth;
    }else{
      throw("Have none unique variant name for ".$self->input_file." ".$name);
    }
  }
  close(IN);
  $self->input_lines(\@input_lines);
  $self->score_hash(\%score_hash);
  $self->total_reads(\%total_reads);
  return ($self->input_lines, $self->score_hash, $self->total_reads);
}

sub calculate_consequences{
  my ($self) = @_;
  my @vfs = @{$self->create_variation_features};
  my $tva = $self->get_transcript_variation_adaptor;
  $tva->fetch_all_by_VariationFeatures(\@vfs);
  $self->variation_features(\@vfs);
}


sub generate_unique_gvf_name{
  my ($self, $name) = @_;
  if($self->{'gvf_unique_name'}->{$name}){
    my @values = split /\:/, $name;
    unless($values[1]){
      $name .= ":01";
      return $self->generate_unique_gvf_name($name);
    }else{
      my $num =~ s/^0*//;
      $num++;
      $num = sprintf("%01d", $num);
      $name .= ":".$num;
      return $self->generate_unique_gvf_name($name);
    }
  }else{
    $self->{'gvf_unique_name'}->{$name} = 1;
    return $name;
  }
}

sub find_in_vcf{
  my ($self) = @_;
  my %hash;
  return undef unless($self->source_vcf && scalar(@{$self->source_vcf}) >= 1);
  my $tmp_file = $self->vcf_positional_input_generator;
  foreach my $in_vcf(@{$self->source_vcf}){
    my $basename = basename($in_vcf);
    my $out_stem = create_filename($self->working_dir, $self->source."_".$basename);
    my $vcf_cmd = $self->vcftools." --vcf ".$in_vcf." --positions ".$tmp_file.
      " --recode --out ".$out_stem;
    eval{
      my $exit = system($vcf_cmd);
    };
    if($@ || $exit >= 1){
      throw("Failed to run ".$vcf_cmd." $exit $!");
    }
    my $files = 

 }
}

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
      throw($string." is not unique it is associated with ".$postion_hash{$string}.
	    " as well as ".$vf->variation_name);
    }
  }
  close(FH);
  $self->files_to_cleanup($tmp_file);
  return $tmp_file, \%position_hash;
}

sub print_output{
  my ($self) = @_;
  my $fh = $self->output_fh;
 # print STDERR "Calculating consequences\n";
  $self->calculate_consequences();
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
    if($self->check_ensembl_sources){
      my $ensembl_sources = $vf->get_all_sources if($vf->variation);
      push(@sources, @$ensembl_sources) if($ensembl_sources);
    }
    foreach my $attrib_string(@$attrib_strings){
      my $output_string = $gff3."\t".$attrib_string."\t";
      $output_string .= join("\t", @sources) unless($self->strict_gvf);
      print $fh $output_string."\n";
    }
  }
}

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
      my $transcript_name = $con->transcript->stable_id;
      my $effect_string = "Variant_effect=".$consequence.":".$feature_type.
	":ensembl:".$transcript_name.";";
      my $string = $attribute_string.$effect_string;
      push(@strings, $string);
    }
  }
  return \@strings;
}


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

sub create_variation_features{
  my ($self) = @_;
  my $sa = $self->get_slice_adaptor;
  throw("Failed to connect to ensembl database ") unless($sa);
  my %slices;
  my @vfs;
  my $vfa = $self->get_variation_feature_adaptor;

  LINE:foreach my $input(@{$self->input_lines}){
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
  elsif($data[0] =~ /(chr)?\w+/ && $data[1] =~ /\d+/ && $data[3] =~ /[ACGTN-]+/ && $data[4] =~ /([ACGTN-]+\,?)+/) {
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

sub create_variant_name{
  my ($self, $seq_name, $start, $end, $allele_string) = @_;
  my $name = join("_", ($seq_name, $start, $end, $allele_string));
  return $name
}

sub DESTROY{
  my ($self) = @_;
  foreach my $file(@{$self->files_to_cleanup}){
    unlink $file;
  }
}

1;


