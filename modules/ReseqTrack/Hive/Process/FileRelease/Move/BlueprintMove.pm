package ReseqTrack::Hive::Process::FileRelease::Move::BlueprintMove;

use strict;
use File::Basename qw(fileparse dirname);
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::GeneralUtils qw(current_date);

use base ('ReseqTrack::Hive::Process::FileRelease::Move');

sub param_defaults {
  my ( $self ) = @_;
  return {
    %{$self->SUPER::param_defaults()},
    'derive_path_options' => { 'collection_tag' => 'experiment_id' },
  };
}

sub derive_path {
  my ( $self, $dropbox_path, $file_object ) = @_;	
  
  my ( $destination, $collection_name );
  
  my $derive_path_options = $self->param( 'derive_path_options' );
  my $aln_base_dir        = $derive_path_options->{aln_base_dir};                     ## alignment files direcroty
  my $results_base_dir    = $derive_path_options->{results_base_dir};                 ## result files directory
  my $run_meta_data_file  = $derive_path_options->{run_meta_data_file};
  my $species             = $derive_path_options->{species};                          ## species
  my $genome_version      = $derive_path_options->{genome_version};                   ## genome_version
  my $freeze_date         = $derive_path_options->{freeze_date};
  my $collection_tag      = $derive_path_options->{collection_tag};

  throw( "this module needs a aln_base_dir" ) 
       if ! defined $aln_base_dir;
  
  throw( "this module needs a results_base_dir" ) 
       if ! defined $results_base_dir;

  throw( "this module needs a run_meta_data_file" ) 
       if ! defined $run_meta_data_file;

  throw( "this module needs a species" ) 
       if ! defined $species;

  throw( "this module needs a genome_version" ) 
       if ! defined $genome_version;
  
  throw( "this module needs a freeze date" ) 
       if ! defined $freeze_date;
  
  my $meta_data = _get_meta_data( $run_meta_data_file );                         ## get metadata hash from file
  my ( $filename, $incoming_dirname ) = fileparse( $dropbox_path );

  my $db = ReseqTrack::DBSQL::DBAdaptor->new(%{$self->param('reseqtrack_db')});
  
  my $file      = $db->get_FileAdaptor->fetch_by_filename( $filename );
  throw( 'file object not found' ) unless $file;
  throw( "expecting one file got ".scalar @$file ) unless scalar @$file == 1;

  my $host      = $$file[0]->host;
  throw( 'host object not found' ) unless $host;
  my $host_name = $host->name; 

  
  if ( $filename =~ m/^(ERR\d{6})/ ) {                                                 ## check for ENA id 
    throw( 'Expection sample name got ENA run id' );                                   ## ENA is support is disabled
    
  }
  elsif ( $host_name eq'CNAG' ) {                                                      ## data in CNAG drop location
  
   ( $destination, $collection_name ) = $self->_derive_CNAG_path
                                           (  filename         => $filename, 
                                              aln_base_dir     => $aln_base_dir, 
                                              results_base_dir => $results_base_dir,
                                              meta_data        => $meta_data,
                                              file_object      => $file_object,
                                              species          => $species,
                                              genome_version   => $genome_version,
                                              freeze_date      => $freeze_date,
                                              collection_tag   => $collection_tag, 
                                            ); 
  }
  elsif ( $host_name eq 'CRG' ) {                                                      ## data in CRG drop location
   
    ( $destination, $collection_name ) = $self->_derive_CRG_path
                                          ( filename         => $filename, 
                                            aln_base_dir     => $aln_base_dir, 
                                            results_base_dir => $results_base_dir,
                                            meta_data        => $meta_data,
                                            file_object      => $file_object,
                                            species          => $species,
                                            freeze_date      => $freeze_date,
                                            genome_version   => $genome_version,
                                            collection_tag   => $collection_tag,
                                          );                                         
  }
  elsif ( $host_name eq 'NCMLS' ) {                                                    ## data in NCMLS drop location
  
    ( $destination, $collection_name ) = $self->_derive_NCMLS_path
                                            ( filename          => $filename, 
                                              aln_base_dir      => $aln_base_dir, 
                                              results_base_dir  => $results_base_dir,
                                              meta_data         => $meta_data,
                                              file_object       => $file_object,
                                              species           => $species,
                                              genome_version    => $genome_version,
                                              freeze_date       => $freeze_date,
                                              collection_tag    => $collection_tag,  
                                            );  
  }
  elsif ( $host_name eq 'WTSI' ) {                                                    ## data in WTSI drop location
  
    ( $destination, $collection_name ) = $self->_derive_WTSI_path
                                            ( filename         => $filename, 
                                              aln_base_dir     => $aln_base_dir, 
                                              results_base_dir => $results_base_dir,
                                              meta_data        => $meta_data,
                                              file_object      => $file_object,
                                              species          => $species,
                                              genome_version   => $genome_version,
                                              freeze_date      => $freeze_date,
                                              collection_tag   => $collection_tag,
                                            );  
  }
  else {
    throw( "Cannot find host information for $dropbox_path" );
  }
  throw( 'collection name not found' ) if !$collection_name;
  
 return $destination, $collection_name; 
 
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
  
  my $filename         = $options{filename}         or throw( "missing filename" );
  my $aln_base_dir     = $options{aln_base_dir}     or throw( "missing aln_base_dir" );
  my $results_base_dir = $options{results_base_dir} or throw( "missing results_base_dir" );
  my $meta_data        = $options{meta_data}        or throw( "missing meta_data object" );
  my $file_object      = $options{file_object}      or throw( "missing file object" );
  my $species          = $options{species}          or throw( 'missing species name' );  
  my $freeze_date      = $options{freeze_date}      or throw( 'missing freeze date' );       ## species name
  my $genome_version   = $options{genome_version}   or throw( 'missing genome version' );
  my $collection_tag   = $options{collection_tag}   or throw( 'missing collection name tag' );
  my $pipeline_name;
  my $output_dir;
  
  
  my @file_fields = split '\.',  $filename;
  my $sample_id = $file_fields[0];
  my $suffix = $file_fields[-1];
  
  my $meta_data_entry = $meta_data->{$sample_id};
  $meta_data_entry = _get_experiment_names( $meta_data_entry );                              ## reset experiment specific hacks
  my $collection_name = $meta_data_entry->{$collection_tag};
  throw( "no collection name for sample $sample_id" ) unless $collection_name;

  throw( "No metadata for sample $sample_id" ) if ( $sample_id && !$meta_data_entry );
  
  if ( $file_object->type =~ m/BAM/ ) {
       $pipeline_name = 'gem_cnag_bs';                                                      ## CNAG pipeline name 
       throw( "expecting suffix bam, got $suffix" ) unless $suffix =~ m/^bam$/i;            ## file suffix check
       $meta_data_entry->{experiment_type} = 'BS';
       $output_dir = $aln_base_dir;                                                         ## reset CNAG experiment type    
  }
  else {
  #     ( $sample_id, $experiment_type, $pipeline_name, $freeze_date, $suffix ) = ( $filename =~ /(\S+?)\.(\S+?)\.(\S+?)\.(\S+?)\.(\S+)/ );
  #     $meta_data_entry = $meta_data->{$sample_id};
  #     throw( "No metadata for sample $sample_id" ) if ( $sample_id && !$meta_data_entry );
       $output_dir =  $results_base_dir;
  }
  
  my $destination = $self->_get_new_path( meta_data_entry => $meta_data_entry,
                                          output_base_dir => $output_dir,
                                          filename        => $filename,
                                          suffix          => $suffix,
                                          pipeline_name   => $pipeline_name,
                                          species         => $species,
                                          genome_version  => $genome_version,
                                          freeze_date     => $freeze_date 
                                        ) or throw( "couldn't get new file path" );                                           
  return $destination, $collection_name; 
}

sub _derive_CRG_path {
  my ( $self, %options ) = @_;
  
  my $filename         = $options{filename}         or throw( 'missing filename' );
  my $aln_base_dir     = $options{aln_base_dir}     or throw( 'missing aln_base_dir' );
  my $results_base_dir = $options{results_base_dir} or throw( 'missing results_base_dir' );
  my $meta_data        = $options{meta_data}        or throw( 'missing meta_data object' );
  my $file_object      = $options{file_object}      or throw( 'missing file object' );
  my $species          = $options{species}          or throw( 'missing species name' );      
  my $freeze_date      = $options{freeze_date}      or throw( 'missing freeze date' );                 ## species name
  my $genome_version   = $options{genome_version}   or throw( 'missing genome version' );
  my $collection_tag   = $options{collection_tag}   or throw( 'missing collection name tag' );

  my @file_fields = split '\.',  $filename;
  my $sample_id     = $file_fields[0];
  my $mark          = $file_fields[1];
  my $pipeline_name = $file_fields[2];
  my $date          = $file_fields[3];
  my $suffix        = $file_fields[-1];  


  my ($run, $big_wig, $output_dir, $is_summary_file );
  
  my $meta_data_entry = $meta_data->{$sample_id};
  my $collection_name = $meta_data_entry->{$collection_tag};
  throw( "no collection name for sample $sample_id" ) unless $collection_name;

  $meta_data_entry = _get_experiment_names( $meta_data_entry );                                        ## reset experiment specific hacks
  
  throw( "No metadata for sample $sample_id" ) if ( $sample_id && !$meta_data_entry );
  
  $pipeline_name = 'gem_grape_crg';                                                                    ## CRG pipeline name
  
  if ( $file_object->type =~ m/BAM/ || $file_object->type =~ m/BAI/ ) {      
      $output_dir = $aln_base_dir;
  }
  elsif ( $file_object->type =~ m/SIGNAL/ ) {
      ( $run, $mark, $big_wig ) = split '\.', $filename;
      $pipeline_name = 'gem_grape_crg';                                          
      $suffix = '.bw' if ( $big_wig eq 'bigwig' );                                                     ## suffix
      if ( $mark eq 'plusRaw' ) {                                                                      ## RNA-Seq strand info
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
    elsif ( $file_object->type =~ m/RNA_TRANSCRIPT_QUANT_CRG/ ) {
      $meta_data_entry->{experiment_type} = 'transcript_quantification';
      $output_dir = $results_base_dir;
    }
    elsif ( $file_object->type =~ m/RNA_GENE_QUANT_CRG/ ) {
      $meta_data_entry->{experiment_type} = 'gene_quantification';
      $output_dir = $results_base_dir;
    }
    elsif ( $file_object->type =~ m/RNA_SPLICING_RATIOS_CRG/ ) {
      $meta_data_entry->{experiment_type} = 'splice_ratios';
      $output_dir = $results_base_dir;
    }
    else {
      throw( "Unsure how to label file $filename " . $file_object->type );
    }
    
    my $destination = $self->_get_new_path( meta_data_entry => $meta_data_entry,
                                            output_base_dir => $output_dir,
                                            filename        => $filename,
                                            suffix          => $suffix,
                                            pipeline_name   => $pipeline_name,
                                            species         => $species,
                                            genome_version  => $genome_version,
                                            freeze_date     => $freeze_date 
                                           ) or throw( "couldn't get new file path" );
    return $destination, $collection_name;                                     
}

sub _derive_NCMLS_path {
  my ( $self, %options ) = @_;
  
  my $filename         = $options{filename}         or throw( "missing filename" );
  my $aln_base_dir     = $options{aln_base_dir}     or throw( "missing aln_base_dir" );
  my $results_base_dir = $options{results_base_dir} or throw( "missing results_base_dir" );
  my $species          = $options{species}          or throw( 'missing species name' );          ## species name
  my $freeze_date      = $options{freeze_date}      or throw( 'missing freeze date' );
  my $file_object      = $options{file_object}      or throw( 'missing file object' );
  my $genome_version   = $options{genome_version}   or throw( 'missing genome version' );
  my $collection_tag   = $options{collection_tag}   or throw( 'missing collection name tag' );

  my ( $sample_id, $mark, $suffix ) = split /\.|_/, $filename;                                   ## NCMLS file name format: ??? ## FIX_IT
  my $meta_data       = $options{meta_data} or throw( "missing meta_data object" );
  my $meta_data_entry = $meta_data->{$sample_id};
  $meta_data_entry = _get_experiment_names( $meta_data_entry );                                  ## reset experiment specific hacks

  my $pipeline_name;
  my $output_dir; 
  my $collection_name = $meta_data_entry->{$collection_tag};  
  throw( "no collection name for sample $sample_id" ) unless $collection_name;

  if ( $file_object->type =~ m/DNASE/ ) {
      $pipeline_name = 'hotspot_ncmls';
  }
  
  if ( $file_object->type =~ m/BAM/ || $file_object->type =~ m/BAI/ ) {      
      $output_dir = $aln_base_dir;
  }
  
  my $destination = $self->_get_new_path( meta_data_entry => $meta_data_entry,
                                          output_base_dir => $output_dir,
                                          filename        => $filename,
                                          suffix          => $suffix,
                                          pipeline_name   => $pipeline_name,
                                          species         => $species,
                                          genome_version  => $genome_version,
                                          freeze_date     => $freeze_date 
                                        ) or throw( "couldn't get new file path" );                                           
  return $destination, $collection_name;
    
}


sub _derive_WTSI_path {
  my ( $self, %options ) = @_;
  
  my $filename         = $options{filename}         or throw( "missing filename" );
  my $aln_base_dir     = $options{aln_base_dir}     or throw( "missing aln_base_dir" );
  my $results_base_dir = $options{results_base_dir} or throw( "missing results_base_dir" );
  my $species          = $options{species}          or throw( 'missing species name' );           ## species name
  my $freeze_date      = $options{freeze_date}      or throw( 'missing freeze date' );
  my $file_object      = $options{file_object}      or throw( 'missing file object' );
  my $genome_version   = $options{genome_version}   or throw( 'missing genome version' );
  my $collection_tag   = $options{collection_tag}   or throw( 'missing collection name tag' );

  my $pipeline_name;
  my $output_dir; 
  
  my ( $sample_id, $type, $algo, $date, $suffix  ) = split '\.', $filename;                       ## WTSI_proposed file format
  $pipeline_name = $type;
  
  my $meta_data = $options{meta_data} or throw( "missing meta_data object" );
  my $meta_data_entry = $meta_data->{$sample_id};
  my $collection_name = $meta_data_entry->{$collection_tag};
  throw( "no collection name for sample $sample_id" ) unless $collection_name;

  $meta_data_entry = _get_experiment_names( $meta_data_entry );                                   ## reset experiment specific hacks
  
  if ( $file_object->type =~ m/BAM/ || $file_object->type =~ m/BAI/ ) {      
      $output_dir = $aln_base_dir;
  }
  
  
  my $destination = $self->_get_new_path( meta_data_entry => $meta_data_entry,
                                          output_base_dir => $output_dir,
                                          filename        => $filename,
                                          suffix          => $suffix,
                                          pipeline_name   => $pipeline_name,
                                          species         => $species,
                                          genome_version  => $genome_version,
                                          freeze_date     => $freeze_date 
                                        ) or throw("couldn't get new file path");                                           
  return $destination, $collection_name;
}

sub _get_new_path {
  my ( $self, %options ) = @_;
  my $destination;
 
  my $meta_data_entry = $options{meta_data_entry} or throw("No meta_data_entry");
  my $output_base_dir = $options{output_base_dir} or throw("No output_base_dir");
  my $filename        = $options{filename}        or throw("No filename");
  my $species         = $options{species}         or throw( 'missing species name' );                                 ### species name
  my $freeze_date     = $options{freeze_date}     or throw( 'missing freeze date' );
  my $suffix          = $options{suffix}          or throw( 'missing suffix' );
  my $pipeline_name   = $options{pipeline_name}   or throw( 'missing pipeline_name' );
  my $genome_version  = $options{genome_version}  or throw( 'missing genome version' );
 
 $species = lc( $species );
 $species = s{\s+}{_}g;
 
 my @file_tokens = (   $meta_data_entry->{sample_name},
                       $meta_data_entry->{experiment_type},
                       $pipeline_name, $freeze_date,  $suffix
                   );
                                     
 my $new_file_name = join( '.', @file_tokens );
 
 if ( $filename eq $new_file_name ) {
  warning ( "RETAINED file name:",$filename," : ",$new_file_name );
 }
 else {
  warning ( "CHANGED file name:",$filename," : ",$new_file_name );
 }
      
 $$meta_data_entry{sample_desc_1} = "NO_TISSUE" 
              if $meta_data_entry->{sample_desc_1} eq "-";

 $$meta_data_entry{sample_desc_2} = "NO_SOURCE" 
              if $meta_data_entry->{sample_desc_2} eq "-";

 $$meta_data_entry{sample_desc_3} = "NO_CELL_TYPE" 
              if $meta_data_entry->{sample_desc_3} eq "-";

 my @dir_tokens = (  $output_base_dir,                 
                     $species,
                     $genome_version,
                     $meta_data_entry->{sample_desc_1},
                     $meta_data_entry->{sample_desc_2},
                     $meta_data_entry->{sample_desc_3},
                     $meta_data_entry->{library_strategy},
                     $meta_data_entry->{center_name} 
                  );
 
  @dir_tokens = map{ ( my $substr = $_ ) =~ s!//!_!g; $substr; } @dir_tokens;                ## not  allowing "/" in token  

  my $dir = join( '/', @dir_tokens );
  
  $destination = $dir . '/' . $new_file_name;
  
  $destination =~ s!//!/!g;
  $destination =~ s/ /_/g;
  $destination =~ s/[ ,;()]/_/g;
  $destination =~ s/_\//\//g;                                                               ## not  allowing "abc_/def"
  $destination =~ s/_+/_/g;
  
  return $destination;
}

sub _get_experiment_names {
  my ( $meta_data_entry ) = @_;
 
  if ( $meta_data_entry->{library_strategy} eq 'RNA-Seq'
    || $meta_data_entry->{library_strategy} eq 'DNase-Hypersensitivity' )  {
    $meta_data_entry->{experiment_type} = $meta_data_entry->{library_strategy};
  }
  
  $meta_data_entry->{experiment_type} =~ s/\QDNase-Hypersensitivity\E/DNase/;
  $meta_data_entry->{experiment_type} =~ s/\QDNA Methylation\E/BS-Seq/;
  
  return $meta_data_entry;
}

1;
