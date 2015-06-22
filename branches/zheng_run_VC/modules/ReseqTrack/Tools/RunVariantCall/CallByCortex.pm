#package RunVariantCall::CallByCortex;
package ReseqTrack::Tools::RunVariantCall::CallByCortex;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename qw(basename fileparse);
use ReseqTrack::Tools::FileSystemUtils qw( check_file_exists check_executable get_lines_from_file);

use base qw(ReseqTrack::Tools::RunVariantCall);
use File::Path;

=head2 new

  Arg [-executable]    :
          string, optional, path to cortex executable that have been compiled 
  Arg [-parameters]    :
          string, optional, parameters needed for running cortex that are not specified by this module (parameters such as kmer_size are specified here so they do not need 
			to be passed in using -paramters
  Arg [-kmer_size]    :
          integer, required, kmer size 
  Arg [-mem_height]    :
          integer, required, memory height
  Arg [-mem_width]    :
          integer, required, memory width
           
          			
               
  + Arguments for ReseqTrack::Tools::RunProgram parent class

  Function  : Creates a new ReseqTrack::Tools::RunVariantCall::CallByCortex object.
  Returntype: ReseqTrack::Tools::RunVariantCall::CallByCortex
  Exceptions: 
  Example   : my $varCall_byCortex = ReseqTrack::Tools::RunVariantCall::CallByCortex->new(
                -input_files             => ['/path/sam1', '/path/sam2'],
                -program                 => "/path/to/Cortex",
                -working_dir             => '/path/to/dir/',  ## this is working dir and output dir
               
                );

=cut

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my (	$executable,
  		$parameters,
  		$kmer_size,
  		$mem_height,
  		$mem_width,
  		$dump_binary,
  		$dump_covg_distribution,
  		$detect_bubbles1,
  		$sample,
  		$collection,
  		$collection_type,
  		$ref_binary,
  		$population,
		)
    = rearrange( [ qw(	EXECUTABLE
    					PARAMETERS
    					KMER_SIZE
    					MEM_HEIGHT
    					MEM_WIDTH
    					DUMP_BINARY
    					DUMP_COVG_DISTRIBUTION
    					DETECT_BUBBLES1
    					SAMPLE
    					COLLECTION
    					COLLECTION_TYPE
    					REF_BINARY
    					POPULATION	 )], @args);    
  
  ## Set defaults
  $self->program("/nfs/1000g-work/G1K/work/bin/cortex/bin/") if (! $self->program);
  if ( !$executable ) {
      throw("Please provide cortex executable");
  } 
  if ( !$kmer_size ) {
      throw("Please provide cortex kmer_size");
  } 
  if ( !$mem_height ) {
      throw("Please provide cortex mem_height");
  } 
  if ( !$mem_width ) {
      throw("Please provide cortex mem_width");
  }  

  $self->executable($executable);
  $self->options($parameters) if ($parameters); 
  $self->kmer_size($kmer_size);   
  $self->mem_height($mem_height);
  $self->mem_width($mem_width);
  $self->dump_binary($dump_binary);
  $self->dump_covg_distribution($dump_covg_distribution);
  $self->detect_bubbles1($detect_bubbles1);
  $self->sample($sample);
  $self->collection($collection);
  $self->collection_type($collection_type);
  $self->ref_binary($ref_binary);
  $self->population($population);
  return $self;
}


sub _open_fh {

	my ($self, $type) = @_;

	my $outfile_name;

	if (	$type eq "SAMPLE_CTX_TO_POOL" || 
			$type eq "POOLED_CTX" || 
			$type eq "POOL_CTX_CLEAN" || 
			$type eq "SAMPLE_CTX_CLEAN_TO_OL_BUBBLE" ||
			$type eq "FASTQ_LIST" ) {
		$outfile_name = $self->working_dir . "/" . $self->population    . "/" . $self->collection . ".list";
	}  
	elsif ( $type eq "ref_binary" ) {
		$outfile_name = $self->working_dir . "/" . $self->population    . "/ref_binary.list";
	} 
	elsif ($type =~ /colour_list/i ) {
		$outfile_name = $self->working_dir . "/" . $self->population    . "/" . $self->collection . ".colour_list";
	}   		
	elsif ( $type eq "SAMPLE_CTX_TO_CLEAN" ) {
		my $tmp_name = $self->collection;
		$tmp_name =~ s/uncleaned.//;
		$tmp_name =~ s/.ctx//;
		$outfile_name = $self->working_dir . "/" . $self->population    . "/" . $tmp_name . ".list";
	} 		
	elsif ( $type =~ /BUBBLE_FASTA/i ) {
		$outfile_name = $self->working_dir . "/" . $self->population    . "/" . $self->collection . ".se_list";
	}  
	
	my $dir = $self->working_dir . "/" . $self->population;
	unless ( -e $dir ) {
		mkpath $dir;
	}
	
	open (my $fh, ">", $outfile_name) || throw("Cannot open output file $outfile_name");
    return ($outfile_name, $fh);
}	

	
sub _write_to_file {
	my ($self, $fh, $content) = @_;
	print $fh $content ."\n";
    return 1;	
}

### FIXME, remove the temparay files after the pipeline is done

sub _process_input {
	my ($self) = @_;
	
	my $input_lists = $self->input_files; # $input_bams can be an array ref of the fastq list files
    print "Input files are:\n";
    print join ("\n", @$input_lists) . "\n";
    check_file_exists(@$input_lists);
    	
    my	($fh1, 
    	$colour_list_fh,
    	$fh2,
    	$fh_ref, 
    	$out1,
    	$colour_list_filename,
    	$out_ref,
    	@ctx_lists);
    
    ($out1, $fh1) = $self->_open_fh($self->collection_type) if ( $self->collection_type ne "SAMPLE_CTX_CLEAN_BUBBLED_POOL"  &&
    																$self->collection_type ne "MULTI_COLOUR_CTX" );
        	
    foreach my $list ( @$input_lists ) {
        ### If colour list is provided, no need to use collection to create colour list, this is used when run manually
        if ($list =~ /colour_list/i) {  
            $self->_colour_list($list);
        }    
        ### To prepare lists for creating per sample graph
		elsif ( $self->collection_type =~ /FASTQ_LIST/i ) {  
			if ($list =~ /se_list/) {
		    	$self->_se_list($list);
		    }
			elsif ($list =~ /pe1_list/) {
		    	$self->_pe1_list($list);
		    }
			elsif ($list =~ /pe2_list/) {
		    	$self->_pe2_list($list);
		    }  
		}
		### To make ctx_list for creating uncleaned or cleaned pool graph
	    elsif (	$self->collection_type =~ /SAMPLE_CTX_TO_POOL/i || 
	    		$self->collection_type eq "POOLED_CTX"  ) { 
	        $self->_write_to_file($fh1, $list);
	    }      
		### To make ctx list and multicolour_bin for cleaning perl sample graph based on cleaned pool
		elsif (	$self->collection_type =~ /SAMPLE_CTX_TO_CLEAN/i || 
	    		$self->collection_type eq "SAMPLE_CTX_CLEAN_TO_OL_BUBBLE") { 
			if 	($list =~ /.cleaned_pool.ctx/i || 
				($list =~ /bubble/i && $list =~ /.ctx/i ) ) {
				$self->_multi_colour_ctx($list);
			}
			else {	
				$self->_write_to_file($fh1, $list);
			}
		}
		### To make a ctx_list for bubble discovery
		elsif ( $self->collection_type =~ /POOL_CTX_CLEAN/i) { 
			$self->_write_to_file($fh1, $list);	
		}		
	    ### To make a se_list for creating BUBBLE graph
		elsif ( $self->collection_type =~ /BUBBLE_FASTA/i) { 
			$self->_write_to_file($fh1, $list);
			$self->_se_list($out1);
		}
		### To make a ctx_list for EACH sample ctx (including ref ctx), to make multi_colour graph for genotyping
		elsif ( $self->collection_type =~ /SAMPLE_CTX_CLEAN_BUBBLED_POOL/i) {
			my $outlist_name = $self->working_dir . "/" . $self->population    . "/" . basename($list) . ".list"; 
			open ($fh2, ">", $outlist_name) || throw("Cannot open output file $outlist_name");  #### the ref binary has to be the first in the colour list
			$self->_write_to_file($fh2, $list);
			push @ctx_lists, $outlist_name;
		}
		### To make colour list for genotyping
		elsif ( $self->collection_type eq "MULTI_COLOUR_CTX" ) {
			if ($list =~ /ctx/i) {
				$self->_multi_colour_ctx($list);
			}
			elsif ( $list =~ /bubble/i ) {
				$self->_bubble_file($list)
			}
		}
		close($fh2) if ($fh2);	
			
    } # end of foreach input_files  
    
    if ($self->ref_binary) {
  	  ($out_ref, $fh_ref) = $self->_open_fh('ref_binary');
	  $self->_write_to_file($fh_ref, $self->ref_binary);  
    }    	
    
    close($fh1) if ($fh1);
    close($fh_ref) if ($fh_ref);
    
    ($colour_list_filename, $colour_list_fh) = $self->_open_fh('colour_list');
    
 	### To create colour_list for pooling uncleaned per sample graphs - this colour list need a second column of population name
    if ( defined $self->collection_type && $self->collection_type =~ /SAMPLE_CTX_TO_POOL/i  ) { 
    	my $pool_name = $self->collection;
    	$pool_name =~ s/_pool//;
    	$pool_name =~ s/.pool//;
    	$pool_name =~ s/_ctx//;
    	$pool_name =~ s/.ctx//;
    	my $content_to_print = $out1 . "\t" . $pool_name . "_pool";
    	$self->_write_to_file($colour_list_fh, $content_to_print);
    	$self->_colour_list($colour_list_filename);
    }  
    ### To create colour_list for other thing, only one column for this colour list   	
    elsif ( defined $self->collection_type && 
    			( 	$self->collection_type eq "POOLED_CTX" || 
    				$self->collection_type eq "SAMPLE_CTX_TO_CLEAN" ||
    				$self->collection_type eq "SAMPLE_CTX_CLEAN_TO_OL_BUBBLE" ) ) { 
    	$self->_write_to_file($colour_list_fh, $out1);
    	$self->_colour_list($colour_list_filename);
    }
    ### To create colour list for detecting bubbles, need reference binary in it             
    elsif ( defined $self->collection_type && $self->collection_type eq "POOL_CTX_CLEAN" ) {
        $self->_write_to_file($colour_list_fh, $out_ref); # the order of the two lines are important, ref always first line
        $self->_write_to_file($colour_list_fh, $out1);
        $self->_colour_list($colour_list_filename);
    }     
    elsif ( defined $self->collection_type &&  $self->collection_type =~ /SAMPLE_CTX_CLEAN_BUBBLED_POOL/i ) {
        foreach my $l ( @ctx_lists ) {					### the ref binary has to be the first in the colour list
            if ($l =~ /ref/i ) {
                $self->_write_to_file($colour_list_fh, $l);        		
            }
        }
        foreach my $l ( @ctx_lists ) {
            next if ($l =~ /ref/i );       
        	$self->_write_to_file($colour_list_fh, $l);
        }
        $self->_colour_list($colour_list_filename);
    }    	
        
        
    close($colour_list_fh) if ($colour_list_fh);
    
    `rm $out1` if ( defined $out1 && (-s $out1) == 0 );    
    `rm $colour_list_filename` if ( defined $colour_list_filename && (-s $colour_list_filename) == 0 );
    
    return $self;
}    

	
sub run_program {
    my ($self) = @_;  
    
	my $q; # quality score
	my $suffix;
	my $cleaning_threshold;
	    
    my $cortex = $self->program . "/" . $self->executable;
    check_executable($cortex);

	$self->_process_input; ## This sub generates se_list, pe1_list, pe2_list, colour_list and multi_colour_ctx from $self->input_files
     
## FIXME, handle situation when some list do not exist

    my $cmd = $cortex . " \\\n";
    $cmd .= "--kmer_size " . $self->kmer_size . " \\\n";
    $cmd .= "--mem_height " . $self->mem_height . " \\\n";
    $cmd .= "--mem_width " . $self->mem_width . " \\\n";    
    $cmd .= "--se_list " . $self->_se_list . " \\\n" if ($self->_se_list);
	$cmd .= "--pe_list " . $self->_pe1_list . "," . $self->_pe2_list . " \\\n" if ($self->_pe1_list && $self->_pe2_list);  
	$cmd .= "--colour_list " . $self->_colour_list . " \\\n" if ( $self->_colour_list );  

	if ( defined $self->options ) {
		foreach my $p ( keys %{$self->options} ) {
            $cmd .= "--" . $p . " " . $self->options->{$p} . " \\\n";
            if ( $p =~ /quality_score_threshold/ ) {
                $q = $self->options->{$p};
            }        
            elsif ( $p =~ /successively_dump_cleaned_colours/ ) {
                $suffix = $self->options->{$p};
            }   
            elsif ( $p =~ /remove_low_coverage_supernodes/ ) {
                $cleaning_threshold = $self->options->{$p};
            }       
        }
    }       
	
	### Dumping per sample graph to system defined output 
	if ($self->_se_list && $self->dump_binary && $self->_se_list !~ /bubbles/i ) {  
		my $per_sample_colour = $self->working_dir . "/" . $self->population    . "/" . $self->sample . ".uncleaned.q" . $q . ".k" . $self->kmer_size . ".ctx";
        $cmd .= "--dump_binary $per_sample_colour \\\n";
    	$self->output_files($per_sample_colour);
    	$cmd .= "--sample_id " . $self->sample .  " \\\n" if ($self->sample);
    }
    ### Dumping uncleaned or cleaned pooled colour graph to system-defined output 
    elsif ( defined $self->collection_type && 
    		( 	$self->collection_type =~ /SAMPLE_CTX_TO_POOL/i  ||  
    			$self->collection_type eq "POOLED_CTX")  ) {  
    	my $binary_output = $self->working_dir . "/" . $self->population    . "/" . $self->collection;
    	$binary_output =~ s/_to_pool//;
    	$binary_output =~ s/.pool//;
    	$binary_output =~ s/_ctx//;
    	$binary_output =~ s/.ctx//;
    	if ( defined $cleaning_threshold ) { 
    	    $binary_output =~ s/.uncleaned//;   
    	    $binary_output =~ s/.ctx//;      
    		$binary_output .= ".t" . $cleaning_threshold . ".cleaned_pool.ctx";  
    	}
    	elsif ( $self->collection_type eq "POOLED_CTX") {
    		my $tmp;
    		($tmp, $cleaning_threshold) = split (/user_input_T/, $self->collection);
    	    $binary_output =~ s/.uncleaned//;   
    	    $binary_output =~ s/.ctx//;   
    	    $binary_output =~ s/_to//;    
    		$binary_output .=  ".cleaned_pool.ctx";   
    		$cmd .= "--remove_low_coverage_supernodes " . $cleaning_threshold .  " \\\n";	
    	}    
    	else {
    	    $binary_output .= ".k" . $self->kmer_size . ".uncleaned_pool.ctx";
    	}     	
    	$cmd .= "--dump_binary " . $binary_output . " \\\n";
    	if ( $self->dump_covg_distribution ) {
    	    $cmd .= "--dump_covg_distribution " . $binary_output . ".covg.txt \\\n";
    	}    
    	$self->output_files($binary_output);
    }	
    ### Cleaning graph per sample based on the pooled cleaned binary
    ### Overlaying sample ctx with bubble ctx
    elsif ( $self->collection_type =~ /SAMPLE_CTX_TO_CLEAN/i ||
    		$self->collection_type eq "SAMPLE_CTX_CLEAN_TO_OL_BUBBLE" 	 ) { 
    	my $per_sample_colours = get_lines_from_file($self->_colour_list);
    	foreach my $per_sample_colour ( @$per_sample_colours ) {
    	    my $output = $self->working_dir . "/" . $self->population    . "/" . basename($per_sample_colour) . "_" . $suffix . ".ctx";
    	    $self->output_files($output);
    	    print "output of per sample cleaning is $output\n";
    	}    
    	$cmd .= "--multicolour_bin " . $self->_multi_colour_ctx . " \\\n";
    }        
	### Detecting bubbles    
    #elsif ( $self->_colour_list && $self->detect_bubbles1 ) { 
    elsif ( defined $self->collection_type && $self->collection_type eq "POOL_CTX_CLEAN" ) {     
        my $bubbles = $self->working_dir . "/" . $self->population    . "/" . $self->collection;
        $bubbles =~ s/ctx/bubbles/;
        $bubbles =~ s/.cleaned_pool//;
        `rm $bubbles` if (-e $bubbles); ### cortex will throw if bubble output exists
        $self->output_files($bubbles);
        $cmd .= "--detect_bubbles1 " . $self->detect_bubbles1 . " \\\n";
        $cmd .= "--output_bubbles1 $bubbles \\\n";
    }    
    ### Dumping one binary for the one fasta file of all bubbles 
    elsif ( $self->_se_list && $self->dump_binary && $self->_se_list =~ /bubbles/i ) { 
    	my $basename = basename($self->_se_list);
    	$basename =~ s/.se_list//;
    	my $bubble_colour = $self->working_dir . "/" . $self->population    . "/" . $basename . ".k" . $self->kmer_size . ".ctx";
        $cmd .= "--dump_binary $bubble_colour \\\n";
    	$self->output_files($bubble_colour);
    }
    elsif ( $self->collection_type eq "SAMPLE_CTX_CLEAN_BUBBLED_POOL" ) {
        my $mcolour_graph_to_genotype = $self->working_dir . "/" . $self->population    . "/" . $self->collection . ".multiColour_ctx";
        $mcolour_graph_to_genotype =~ s/_to_pool//;
        $cmd .= "--dump_binary $mcolour_graph_to_genotype \\\n";
        $self->output_files($mcolour_graph_to_genotype);
    }  
    ### Genotyping   
    elsif ( $self->collection_type eq "MULTI_COLOUR_CTX" ) {
        my $outfile_name = $self->working_dir . "/" . $self->population    . "/" . $self->collection . ".genotype";
        $outfile_name =~ s/.multiColour_ctx//;
        $outfile_name =~ s/.ctx_pool//;
        $outfile_name =~ s/.bubbles//;
        $outfile_name =~ s/_ctxs_to_pool//;
        
        $cmd .= "--multicolour_bin " . $self->_multi_colour_ctx . " \\\n";
        
        $cmd .= "--gt " . $self->_bubble_file . "," . $outfile_name . ",BC \\\n"; 
        #### [--gt] option requires two filenames, plus either BC or PD, comma separated. The filenames are to be an input (file of cortex bubble calls) and output filename. 
        $cmd .= "--print_colour_coverages ";
        
        $self->output_files($outfile_name);
    }     
      
	print "Running command:\n$cmd\n";      

    $self->execute_command_line($cmd);
    
    return $self;
}


=head2 kmer_size

  Arg [1]   : ReseqTrack::Tools::RunVariantCall::CallByCortex
  Arg [2]   : string, path of kmer_size
  Function  : accessor method for kmer_size
  Returntype: string
  Exceptions: n/a
  Example   : my $kmer_size = $self->kmer_size;

=cut

sub kmer_size {
  my ($self, $kmer_size) = @_;
  if ($kmer_size) {
    $self->{'kmer_size'} = $kmer_size;
  }
  return $self->{'kmer_size'};
}


=head2 mem_height

  Arg [1]   : ReseqTrack::Tools::RunVariantCall::CallByCortex
  Arg [2]   : string, path of mem_height
  Function  : accessor method for mem_height
  Returntype: string
  Exceptions: n/a
  Example   : my $mem_height = $self->mem_height;

=cut

sub mem_height {
  my ($self, $mem_height) = @_;
  if ($mem_height) {
    $self->{'mem_height'} = $mem_height;
  }
  return $self->{'mem_height'};
}



=head2 mem_width

  Arg [1]   : ReseqTrack::Tools::RunVariantCall::CallByCortex
  Arg [2]   : string, path of mem_width
  Function  : accessor method for mem_width
  Returntype: string
  Exceptions: n/a
  Example   : my $mem_width = $self->mem_width;

=cut

sub mem_width {
  my ($self, $mem_width) = @_;
  if ($mem_width) {
    $self->{'mem_width'} = $mem_width;
  }
  return $self->{'mem_width'};
}



=head2 executable

  Arg [1]   : ReseqTrack::Tools::RunVariantCall::CallByCortex
  Arg [2]   : string, path of executable
  Function  : accessor method for executable
  Returntype: string
  Exceptions: n/a
  Example   : my $executable = $self->executable;

=cut

sub executable {
  my ($self, $executable) = @_;
  if ($executable) {
    $self->{'executable'} = $executable;
  }
  return $self->{'executable'};
}



=head2 dump_binary

  Arg [1]   : ReseqTrack::Tools::RunVariantCall::CallByCortex
  Arg [2]   : string, path of dump_binary
  Function  : accessor method for dump_binary
  Returntype: boolean
  Exceptions: n/a
  Example   : my $dump_binary = $self->dump_binary;

=cut

sub dump_binary {
  my ($self, $dump_binary) = @_;
  if ($dump_binary) {
    $self->{'dump_binary'} = $dump_binary;
  }
  return $self->{'dump_binary'};
}

=head2 dump_covg_distribution

  Arg [1]   : ReseqTrack::Tools::RunVariantCall::CallByCortex
  Arg [2]   : string, path of dump_covg_distribution
  Function  : accessor method for dump_covg_distribution
  Returntype: boolean
  Exceptions: n/a
  Example   : my $dump_covg_distribution = $self->dump_covg_distribution;

=cut

sub dump_covg_distribution {
  my ($self, $dump_covg_distribution) = @_;
  if ($dump_covg_distribution) {
    $self->{'dump_covg_distribution'} = $dump_covg_distribution;
  }
  return $self->{'dump_covg_distribution'};
}


=head2 detect_bubbles1

  Arg [1]   : ReseqTrack::Tools::RunVariantCall::CallByCortex
  Arg [2]   : string, to tell which colour to detect bubbles in
  Function  : accessor method for detect_bubbles1
  Returntype: string
  Exceptions: n/a
  Example   : my $detect_bubbles1 = $self->detect_bubbles1;

=cut

sub detect_bubbles1 {
  my ($self, $detect_bubbles1) = @_;
  if ($detect_bubbles1) {
    $self->{'detect_bubbles1'} = $detect_bubbles1;
  }
  return $self->{'detect_bubbles1'};
}


=head2 sample

  Arg [1]   : ReseqTrack::Tools::RunVariantCall::CallByCortex
  Arg [2]   : string, path of sample
  Function  : accessor method for sample
  Returntype: string
  Exceptions: n/a
  Example   : my $sample = $self->sample;

=cut

sub sample {
  my ($self, $sample) = @_;
  if ($sample) {
    $self->{'sample'} = $sample;
  }
  return $self->{'sample'};
}


=head2 collection

  Arg [1]   : ReseqTrack::Tools::RunVariantCall::CallByCortex
  Arg [2]   : string, path of collection
  Function  : accessor method for collection
  Returntype: string
  Exceptions: n/a
  Example   : my $collection = $self->collection;

=cut

sub collection {
  my ($self, $collection) = @_;
  if ($collection) {
    $self->{'collection'} = $collection;
  }
  return $self->{'collection'};
}


=head2 collection_type

  Arg [1]   : ReseqTrack::Tools::RunVariantCall::CallByCortex
  Arg [2]   : string, path of collection_type
  Function  : accessor method for collection_type
  Returntype: string
  Exceptions: n/a
  Example   : my $collection_type = $self->collection_type;

=cut

sub collection_type {
  my ($self, $collection_type) = @_;
  if ($collection_type) {
    $self->{'collection_type'} = $collection_type;
  }
  return $self->{'collection_type'};
}

=head2 ref_binary

  Arg [1]   : ReseqTrack::Tools::RunVariantCall::CallByCortex
  Arg [2]   : string, path of ref_binary
  Function  : accessor method for ref_binary
  Returntype: string
  Exceptions: n/a
  Example   : my $ref_binary = $self->ref_binary;

=cut

sub ref_binary {
  my ($self, $ref_binary) = @_;
  if ($ref_binary) {
    $self->{'ref_binary'} = $ref_binary;
  }
  return $self->{'ref_binary'};
}

=head2 population

  Arg [1]   : ReseqTrack::Tools::RunVariantCall::CallByCortex
  Arg [2]   : string, path of population
  Function  : accessor method for population
  Returntype: string
  Exceptions: n/a
  Example   : my $population = $self->population;

=cut

sub population {
  my ($self, $population) = @_;
  if ($population) {
    $self->{'population'} = $population;
  }
  return $self->{'population'};
}



=head2 _colour_list

  Arg [1]   : ReseqTrack::Tools::RunVariantCall::CallByCortex
  Arg [2]   : string, path of colour_list
  Function  : accessor method for colour_list
  Returntype: string
  Exceptions: n/a
  Example   : my $colour_list = $self->colour_list;

=cut

sub _colour_list {
  my ($self, $colour_list) = @_;
  if ($colour_list) {
    $self->{'colour_list'} = $colour_list;
  }
  return $self->{'colour_list'};
}

=head2 _se_list

  Arg [1]   : ReseqTrack::Tools::RunVariantCall::CallByCortex
  Arg [2]   : string, path of se_list
  Function  : accessor method for se_list
  Returntype: string
  Exceptions: n/a
  Example   : my $se_list = $self->se_list;

=cut

sub _se_list {
  my ($self, $se_list) = @_;
  if ($se_list) {
    $self->{'se_list'} = $se_list;
  }
  return $self->{'se_list'};
}

=head2 _pe1_list

  Arg [1]   : ReseqTrack::Tools::RunVariantCall::CallByCortex
  Arg [2]   : string, path of pe1_list
  Function  : accessor method for pe1_list
  Returntype: string
  Exceptions: n/a
  Example   : my $pe1_list = $self->pe1_list;

=cut

sub _pe1_list {
  my ($self, $pe1_list) = @_;
  if ($pe1_list) {
    $self->{'pe1_list'} = $pe1_list;
  }
  return $self->{'pe1_list'};
}


=head2 _pe2_list

  Arg [1]   : ReseqTrack::Tools::RunVariantCall::CallByCortex
  Arg [2]   : string, path of pe2_list
  Function  : accessor method for pe2_list
  Returntype: string
  Exceptions: n/a
  Example   : my $pe2_list = $self->pe2_list;

=cut

sub _pe2_list {
  my ($self, $pe2_list) = @_;
  if ($pe2_list) {
    $self->{'pe2_list'} = $pe2_list;
  }
  return $self->{'pe2_list'};
}

=head2 _multi_colour_ctx

  Arg [1]   : ReseqTrack::Tools::RunVariantCall::CallByCortex
  Arg [2]   : string, path of multi_colour_ctx
  Function  : accessor method for multi_colour_ctx
  Returntype: string
  Exceptions: n/a
  Example   : my $multi_colour_ctx = $self->_multi_colour_ctx;

=cut

sub _multi_colour_ctx {
  my ($self, $multi_colour_ctx) = @_;
  if ($multi_colour_ctx) {
    $self->{'multi_colour_ctx'} = $multi_colour_ctx;
  }
  return $self->{'multi_colour_ctx'};
}


=head2 _bubble_file

  Arg [1]   : ReseqTrack::Tools::RunVariantCall::CallByCortex
  Arg [2]   : string, path of bubble_file
  Function  : accessor method for bubble_file
  Returntype: string
  Exceptions: n/a
  Example   : my $bubble_file = $self->_bubble_file;

=cut

sub _bubble_file {
  my ($self, $bubble_file) = @_;
  if ($bubble_file) {
    $self->{'bubble_file'} = $bubble_file;
  }
  return $self->{'bubble_file'};
}


1;