package ReseqTrack::Hive::Process::FileRelease::Move::BlueprintMove;
#package Blueprint::FileRelease::Move;

use strict;
use File::Basename qw(fileparse dirname);
use ReseqTrack::Tools::Exception qw(throw);

use base ('ReseqTrack::Hive::Process::FileRelease::Move');

sub param_defaults {
  my ( $self ) = @_;
  return {
    %{$self->SUPER::param_defaults()},
  };
}

sub derive_path {
  my ( $self, $dropbox_path, $file_object ) = @_;	
  my $derive_path_options = $self->param( 'derive_path_options' );
  
  my $aln_base_dir = $derive_path_options->{aln_base_dir};   ### alignment files direcroty
  throw( "this module needs a aln_base_dir" ) if ! defined $aln_base_dir;
  
  my $results_base_dir = $derive_path_options->{results_base_dir};   ### result files directory
  throw( "this module needs a results_base_dir" ) if ! defined $results_base_dir;
  
  my $run_meta_data_file =  $self->param( 'run_meta_data_file' );
  throw( "this module needs a run_meta_data_file" ) if ! defined $run_meta_data_file;
  
  my $meta_data = _get_meta_data($run_meta_data_file);

  my ( $filename, $incoming_dirname ) = fileparse( $dropbox_path );

  my ( $destination_base_dir, $meta_data_entry, $run_id, $sample_id, $pipeline_name, $suffix, $mark );

  if ( $file_object->type =~ m/BAM/ || $file_object->type =~ m/BAI/ ) {  ### alignment files type
    $destination_base_dir = $aln_base_dir;
  }
  my $destination;
    
  if ( $filename =~ m/(ERR\d{6})/ ) {  ### check for ENA id 
    $run_id = $1;
    $meta_data_entry = $meta_data->{$run_id} if $run_id;
    throw( "No metadata for run id $run_id" ) if ( $run_id && !$meta_data_entry );
  }
  elsif ( $incoming_dirname =~ m{/CNAG/} ) {  ### data in CNAG drop location
  
    $destination = $self->_derive_CNAG_path(  filename => $filename, 
                                              aln_base_dir => $aln_base_dir, 
                                              results_base_dir => $results_base_dir,
                                              meta_data => $meta_data,
                                              file_object => $file_object  
                                            );
   
  }
  elsif ( $incoming_dirname =~ m{/CRG/} ) {  ### data in CRG drop location
   
       $destination = $self->_derive_CRG_path( filename => $filename, 
                                               aln_base_dir => $aln_base_dir, 
                                               results_base_dir => $results_base_dir,
                                               meta_data => $meta_data,
                                               file_object => $file_object  
                                             );
                                         
     
     
     
     
     
     
     
  }
  elsif ( $incoming_dirname =~ m{/NCMLS/} ) {   ### data in NCMLS drop location
    ( $sample_id, $mark ) = split /\.|_/, $filename;
    $meta_data_entry = $meta_data->{$sample_id};
    throw( "No metadata for sample $sample_id" ) if ( $sample_id && !$meta_data_entry );
    
    if ( $file_object->type =~ m/DNASE/ ) {
      $pipeline_name = 'hotspot_ncmls';
    }
  }
  else {
    print STDERR ("Cannot find run id from $filename");
  }
  
  if ( $meta_data_entry->{library_strategy} eq 'RNA-Seq'
    || $meta_data_entry->{library_strategy} eq 'DNase-Hypersensitivity' )
  {
    $meta_data_entry->{experiment_type} = $meta_data_entry->{library_strategy};
  }
  
  my $pipeline_name = undef;
  my $suffix        = undef;

  if ( $expected_file_name =~ m/(\.\w+(\.bai)?(\.gz)?)$/ ) {
    $suffix = $1;
  }
  
  
 
}   

sub _get_meta_data {
  my ( $meta_data_file, ) = @_;
  my %meta_data;
  my @headers;

  open my $mdfh, $meta_data_file or die "Could not open $meta_data_file: $!";

  while (<$mdfh>) {
    chomp;
    my @vals = split "\t", $_;
    if (@headers) {
      my %row;
      @row{@headers}                  = @vals;
     # $meta_data{ $row{run_id} }      = \%row;
      $meta_data{ $row{sample_name} } = \%row;
    }
    else {
      @headers = map { lc($_) } @vals;
    }
  }

  close $mdfh;

  return \%meta_data;
}

sub _derive_CNAG_path {
  my ( $self, %options ) = @_;
  my $destination ;
  
  my $filename = $options{filename} or throw( "missing filename" );
  my $aln_base_dir = $options{aln_base_dir} or throw( "missing aln_base_dir" );
  my $results_base_dir = $options{results_base_dir} or throw( "missing results_base_dir" );
  my $meta_data = $options{meta_data} or throw( "missing meta_data object" );
  my $file_object = $options{file_object} or throw( "missing file object" );
  
  my ( $sample_id, $suffix ) = split /\./,  $filename;
  my $meta_data_entry = $meta_data->{$sample_id};
  throw( "No metadata for sample $sample_id" ) if ( $sample_id && !$meta_data_entry );
  
  if ( $file_object->type =~ m/BAM/ ) {
       $pipeline_name = 'gem_cnag_bs';                                           ### CNAG pipeline name 
       throw("expectinf suffix bam, got $suffix") unless $suffix =~ m/^bam$/i;   ### file suffix check
       $meta_data_entry->{experiment_type} = 'BS';                               ### CNAG experiment type
   
       $destination = $self->_get_new_path( meta_data_entry => $meta_data_entry,
                                            output_base_dir => $aln_base_dir,
                                            filename => $filename
                                           ) or throw("couldn't get new file path");
  }
  return $destination; 
}

sub _derive_CRG_path {
  my ( $self, %options ) = @_;
  
  my $filename = $options{filename} or throw( "missing filename" );
  my $aln_base_dir = $options{aln_base_dir} or throw( "missing aln_base_dir" );
  my $results_base_dir = $options{results_base_dir} or throw( "missing results_base_dir" );
  my $meta_data = $options{meta_data} or throw( "missing meta_data object" );
  my $file_object = $options{file_object} or throw( "missing file object" );
  
  my ( $sample_id, $mark, $pipeline_name, $date, $suffix ) = split /\./,  $filename;
  
  my ($run, $big_wig, $output_dir, $is_summary_file );
  
  my $meta_data_entry = $meta_data->{$sample_id};
  throw( "No metadata for sample $sample_id" ) if ( $sample_id && !$meta_data_entry );
  
  $pipeline_name = 'gem_grape_crg';                                                          ### CRG pipeline name
  
  if ( $file_object->type =~ m/BAM/ || $file_object->type =~ m/BAI/ ) {
      
      $output_dir = $aln_base_dir;
  }
  elsif ( $file_object->type =~ m/SIGNAL/ ) {
      ( $run, $mark, $big_wig ) = split /\./, $filename;
      $pipeline_name = 'gem_grape_crg';                                          
      $suffix = '.bw' if ( $big_wig eq 'bigwig' );                               ### suffix
      if ( $mark eq 'plusRaw' ) {                                                ### RNA-Seq strand info
        $mark = 'plusStrand';
      }
      if ( $mark eq 'minusRaw' ) {
        $mark = 'minusStrand';
      }

      $meta_data_entry->{experiment_type} = $mark;
      
      $output_dir = $results_base_dir;
    } 
    elsif ( $file_object->type =~ m/CONTIGS/ ) {
      $meta_data_entry->{experiment_type} = 'contigs';      
      $output_dir = $results_base_dir;      
    }
    elsif ( $file_object->type =~ m/RNA_JUNCTIONS_CRG/ ) {
      $meta_data_entry->{experiment_type} = 'splice_junctions';
      $output_dir = $results_base_dir;
    }
    elsif ( $file_object->type =~ m/RNA_EXON_QUANT_CRG/ ) {
      $meta_data_entry->{experiment_type} = 'exon_quantification';
      $output_dir = $results_base_dir;
    }
    elsif ( $file_objec->type =~ m/RNA_TRANSCRIPT_QUANT_CRG/ ) {
      $meta_data_entry->{experiment_type} = 'transcript_quantification';
      $output_dir = $results_base_dir;
    }
    elsif ( $file_objec->type =~ m/RNA_GENE_QUANT_CRG/ ) {
      $meta_data_entry->{experiment_type} = 'gene_quantification';
      $output_dir = $results_base_dir;
    }
    elsif ( $file_objec->type =~ m/RNA_SPLICING_RATIOS_CRG/ ) {
      $meta_data_entry->{experiment_type} = 'splice_ratios';
      $output_dir = $results_base_dir;
    }
    else {
      throw( "Unsure how to label file $filename " . $file_objec->type );
    }
}

sub _get_new_path {
 my ( $self, %options ) = @_;
 my $destination;
 
 my $meta_data_entry = $options{meta_data_entry} or throw("No meta_data_entry");
 my $output_base_dir = $options{output_base_dir} or throw("No output_base_dir");
 my $filename = $options{filename} or throw("No filename");
 my $species = $options{species} or 'homo_sapiens';                                 ### species name
 
 my @dir_tokens = (  $output_dir,                 
                     $species,
                     $meta_data_entry->{sample_desc_1},
                     $meta_data_entry->{sample_desc_2},
                     $meta_data_entry->{sample_desc_3},
                     $meta_data_entry->{library_strategy},
                     $meta_data_entry->{center_name} 
                  );
  my $dir = join( '/', @dir_tokens );
  
  $dir =~ s!//!/!g;
  $dir =~ s/ /_/g;
  $dir =~ s/[ ,;()]/_/g;
  $dir =~ s/_\//\//g; ## not  allowing "abc_/def"
  $dir =~ s/_+/_/g;
  
  $destination = $dir . '/' . $file_name;
  return $destination;
}

1;
