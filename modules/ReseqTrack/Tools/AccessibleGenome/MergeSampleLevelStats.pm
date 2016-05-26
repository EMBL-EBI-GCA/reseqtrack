package AccessibleGenome::MergeSampleLevelStats;

use strict;
use ReseqTrack::Tools::Exception;
use ReseqTrack::Tools::CollectionUtils;
use ReseqTrack::DBSQL::DBAdaptor;
use File::Basename;
use Getopt::Long;
use ReseqTrack::Tools::Loader::File;
use ReseqTrack::Tools::Argument qw(rearrange);
use base qw(ReseqTrack::Tools::RunProgram);

=head2 new

  Arg [-file_type]   :
      string, file type
  
  Function  : Creates a new AccessibleGenome::MergeSampleLevelStats object.
  Returntype: AccessibleGenome::MergeSampleLevelStats
  Exceptions: 
  Example   : my $merger = AccessibleGenome::MergeSampleLevelStats->new(
								-db							=> $db,
								-hostname					=> $hostname,
								-dbhost						=> $dbhost,
								-dbname						=> $dbname,
								-dbuser						=> $dbuser,
								-dbpass						=> $dbpass,
								-dbport						=> $dbport,
								-input_files				=> \@file_list,
								-program					=> $merge_program,
								-file_type					=> $file_type,
								-target_count				=> $target_count,
								-goahead					=> $run,
								-verbose					=> $verbose,
								-working_dir				=> $output_dir,
							);

=cut

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my (  $db, $dbhost, $dbname, $dbuser, $dbpass, $dbport, $hostname,
		$goahead, $verbose, $file_type, $target_count, $chr_list) =
    rearrange(
    [
      qw( DB DBHOST DBNAME DBUSER DBPASS DBPORT HOSTNAME
      GOAHEAD VERBOSE FILE_TYPE TARGET_COUNT CHR_LIST)
    ],
    @args
    );
	
	$self->dbhost($dbhost);
	$self->dbname($dbname);
	$self->dbuser($dbuser);
	$self->dbpass($dbpass);
	$self->dbport($dbport);
	$self->goahead($goahead);
	$self->hostname($hostname);
	$self->verbose($verbose);	
	
	if ($db) {
		$self->db($db);
	}	

	$self->target_count($target_count);
	$self->file_type($file_type);
	$self->chr_list($chr_list);

	return $self;
}

=head2 run

  Function  : take the list of input files and organize them into groups to be merged.
  			  For each group, 
  			  1. make a collection, 
  			  2. change the file type to type_ONGOING
  			  3. submit a farm job to merge
  			  	3a. merge
  			  	3b. change file type to type_MERGED upon successful merge
  			  	3c. change collection type to type_MERGED upon successful merge
  			  	
              Output are merged files at different stages, recorded in the db.  The final merge product is the last file that can no longer merge

  Returntype: 
  Exceptions: 
  Example   : $self->run;

=cut

sub run {
	my ( $self ) = @_;
	
	$self->db->dbc->disconnect_when_inactive(1);
	my $fa = $self->db->get_FileAdaptor;
	my $ca = $self->db->get_CollectionAdaptor;

	my $cnt = 1;
	my @merge_list;
	my @merge_fos_list;

	foreach my $file ( @{$self->input_files}) {
		print "input file is $file\n";
		my $fo = $fa->fetch_by_name($file);
		
		if ($cnt < $self->target_count) {
			push @merge_list, $file;
			push @merge_fos_list, $fo;
			$cnt++;
		}
		elsif ($cnt ==  $self->target_count) {
			push @merge_list, $file;
			push @merge_fos_list, $fo;

			$self->load_collection(\@merge_fos_list, $ca);
			$self->update_file_type(\@merge_fos_list, "ONGOING", $fa);
			
			#$self->merge_list(\@merge_list);
			#$self->merge_fos_list(\@merge_fos_list);
			#$self->fa($fa);
			#$self->ca($ca);

			my $merged_file = $self->working_dir . "/merge/" . $self->collection->name . ".stats.gz";			
			my $lsf_log = $self->working_dir . "/lsf_log/" . $self->collection->name . ".stats.gz.log";
			my $lsf_err = $self->working_dir . "/lsf_log/" . $self->collection->name . ".stats.gz.err";
			
			my $tmp_list = $merged_file . ".list";
			
			open(OUT, ">", $tmp_list) || throw("Cannot open tmp file $tmp_list");
			print OUT join ("\n", @merge_list) . "\n";
			
			## Need to add job group 
			## bgadd -L 50 /MyJobGroup
			
			my $script = "/nfs/production/reseq-info/work/zheng/accessible_genome_mask/modules/AccessibleGenome/run_mergeBaseQCSumStats.pl";
			#my $cmd = "bsub -g \/merge_stats -R \"rusage[mem=2000] select[panfs_nobackup_production]\" -q production-rh6 -oo $lsf_log -eo $lsf_err ";
			my $cmd = "bsub -g \/merge_stats -R \"select[panfs_nobackup_production]\" -q production-rh6 -oo $lsf_log -eo $lsf_err ";
			$cmd .= "perl $script ";
			$cmd .=	" -dbhost " . 		$self->dbhost;
			$cmd .=	" -dbname " . 		$self->dbname;
			$cmd .=	" -dbuser " . 		$self->dbuser;
			$cmd .=	" -dbpass " . 		$self->dbpass;
			$cmd .=	" -dbport " . 		$self->dbport;
			$cmd .= " -hostname " . 		$self->hostname;
			$cmd .= " -merge_program " . $self->program;
			$cmd .= " -file_type " . 	$self->file_type;			
			$cmd .=	" -file_list " . 	$tmp_list;
			$cmd .=	" -merged_file_name " . $merged_file;	
			$cmd .= " -chr_list " .		$self->chr_list;		
			$cmd .= " -col_name " . 	$self->collection->name;
			$cmd .= " -col_type " . 	$self->collection->type;
			$cmd .=	" -goahead " 		if ($self->goahead);
			$cmd .= " -verbose	" 		if ($self->verbose);			
			
			#print "$cmd\n";
				
			$self->execute_command_line($cmd) if ($self->goahead);	
			
			my $exit = $?>>8;
	        if ($exit >=1) {
	            throw("merge failed\n"); 
	        }  
			
			$cnt = 1;
			@merge_list = ();
			@merge_fos_list = ();
		}
	}	
}	


### SUBS ###
sub load_file {
	my ($self, $file_name) = @_;

	my $loader = ReseqTrack::Tools::Loader::File->new(						 
						  -file			=> [$file_name],
						  -type			=> $self->file_type,
						  -hostname		=> $self->hostname,
						  -dbhost		=> $self->dbhost,
						  -dbname		=> $self->dbname,
						  -dbuser		=> $self->dbuser,
						  -dbpass		=> $self->dbpass,
						  -dbport		=> $self->dbport,
						  -do_md5			=>1,
						  -update_existing	=>0,
						  -die_for_problems	=>1,
						  -assign_type		=>1,
						  #-store_new		=>1, ## this is to force "store" not "update"
						  -debug			=>1,
						  -verbose			=>$self->verbose,
						 );
	
	$loader->process_input();
	$loader->create_objects();
	$loader->sanity_check_objects();
	$loader->load_objects() if ($self->goahead);

	return $self;
}


sub load_collection {
	my ($self, $file_objs, $col_a) = @_;
	
	## Make up a unique collection_name.
	## It is the sample name of the first sample, with a count;
	## the count is how many times this sample has been the first of a merge list
	my $first_file = basename($file_objs->[0]->name);
	my ($sample_of_first_file, $blah) = split(/\./, $first_file);
	my ($sample_name, $count) = split(/\_/,  $sample_of_first_file);
		 
	if (!$count) {
		$count = 1;
	}
	else {
		$count++;
	}			 
	
	my $collection_name = $sample_name . "_". $count;  
	
	#print "type is " . $self->file_type . "\n";
	print "col name is $collection_name\n";
	
	throw("collection $collection_name exists") if ($col_a->fetch_by_name_and_type($collection_name, $self->file_type));	
	
    my $collection =  ReseqTrack::Collection->new(
                -name => $collection_name,
                -others => $file_objs,
                -type => $self->file_type,
                -table_name => 'file',
	);

    $col_a->store($collection) if ($self->goahead);
    
    $self->collection($collection);

	return $self;
}

sub update_file_type {
	my ($self, $fos, $suffix, $file_a) = @_;
	
	foreach my $fo (@$fos) {

	    my $new_file_type;
	    if ($suffix) {
	    	$new_file_type = $self->file_type . "_" . $suffix;
		}
		else {
		 	$new_file_type = $self->file_type;  ## restore type to original type
		}	 
		$fo->type($new_file_type);

        my $history_ref = $fo->history;
        my $history;

        if (!$history_ref || @$history_ref == 0) {
                $history = ReseqTrack::History->new(
                        -other_id => $fo->dbID,
                        -table_name => 'file',
                        -comment => "This stats file has been changed",
                );
                $fo->history($history);
        }
        else {
                $history = ReseqTrack::History->new(
                        -other_id => $fo->dbID,
                        -table_name => 'file',
                        -comment => "Change file type",
                );
				$fo->history($history);
        }   
                    		
		$file_a->update($fo, 1, 1) if ($self->goahead); # the second 1 is allow change name
	
	}
	return $self;
}	

sub update_collection_type {
	my ($self, $suffix, $col_a, $col_obj) = @_;
	
	my $new_col;
	my $new_type;
	my $dup_col_obj;
	
	## The $col_obj is passed on from a freshly initiated AccessibleGenome::MergeSampleLevelStats object by script run_mergeBaseQCSumStats.pl running as a farm job
	
	if ($col_obj) {
		$new_col = copy_collection_object($col_obj);
		$new_type = $col_obj->type . "_". $suffix;
		$dup_col_obj = $col_a->fetch_by_name_and_type($col_obj->name, $new_type);
		if ( $dup_col_obj ) { # This is to avoid the case when the collection has been updated previously
			     print STDERR "existing collection ". $col_obj->name . " with type $new_type\n";
			     $new_type = $new_type . "_2";
		}    
			 
		$new_col->type($new_type);		 
		my $history = create_collection_history_object($new_col, $col_obj);
		print "Have " . $new_col->name . " " . $history->comment."\n" if ($self->verbose);
		if($history){
		 	$new_col->history($history) if ($self->goahead);
		 	$col_a->update_type($new_col) if ($self->goahead);
		}
	}
	else {
		$new_col = copy_collection_object($self->collection);	
		$new_type = $self->collection->type . "_". $suffix;
			 
		$dup_col_obj = $col_a->fetch_by_name_and_type($self->collection->name, $new_type);
		if ( $dup_col_obj ) { # This is to avoid the case when the collection has been updated previously
			 print STDERR "existing collection ". $self->collection->name . " with type $new_type\n";
			 $new_type = $new_type . "_2";
		}    
			 
		$new_col->type($new_type);		 
		my $history = create_collection_history_object($new_col, $self->collection);
		print "Have " . $new_col->name . " " . $history->comment."\n"  if ($self->verbose);
		if($history){
		 	$new_col->history($history) if ($self->goahead);
		 	$col_a->update_type($new_col) if ($self->goahead);
		}
	}
	return $self;
		 	
}

#### ACCESS SUBS ####

sub db {
 my ( $self, $arg ) = @_;

  $self->{db} = $arg if ( $arg) ;
  
   if(!$self->{db}){
    # print "Have to create DBAdaptor\n";
    $self->{db} = $self->create_DBAdaptor();
  }
 return $self->{db};
}
###
sub dbhost {
 my ( $self, $arg ) = @_;
   $self->{dbhost} = $arg if ( $arg);
 return $self->{dbhost};
}
###
sub dbname {
 my ( $self, $arg ) = @_;
   $self->{dbname} = $arg if ( $arg);
  return $self->{dbname};
}
###
sub dbuser {
 my ( $self, $arg ) = @_;
   $self->{dbuser} = $arg if ( $arg);
  return $self->{dbuser};
}
###
sub dbpass {
 my ( $self, $arg ) = @_;
  $self->{dbpass} = $arg if ( $arg);
 return $self->{dbpass};
}
###
sub dbport {
 my ( $self, $arg ) = @_;
  $self->{dbport} = $arg if ($arg);
 return $self->{dbport};
}
###

###
sub hostname {
 my ( $self, $arg ) = @_;
  $self->{hostname} = $arg if ($arg);
 return $self->{hostname};
}
###
        	
sub target_count {
 my ( $self, $arg ) = @_;
  $self->{target_count} = $arg if ($arg);
 return $self->{target_count};
}

###
sub goahead {
 my ( $self, $arg ) = @_;
  $self->{goahead} = $arg if ($arg);
 return $self->{goahead};
}

###
sub collection {
 my ( $self, $arg ) = @_;
  $self->{collection} = $arg if ($arg);
 return $self->{collection};
}

###
sub file_type {
 my ( $self, $arg ) = @_;
  $self->{file_type} = $arg if ($arg);
 return $self->{file_type};
}

###
sub chr_list {
 my ( $self, $arg ) = @_;
  $self->{chr_list} = $arg if ($arg);
 return $self->{chr_list};
}

###
sub verbose {
 my ( $self, $arg ) = @_;
  $self->{verbose} = $arg if ($arg);
 return $self->{verbose};
}

=head
sub merge_list {
 my ( $self, $arg ) = @_;
   $self->{merge_list} = $arg if ( $arg);
 return $self->{merge_list};
}

sub merge_fos_list {
 my ( $self, $arg ) = @_;
   $self->{merge_fos_list} = $arg if ( $arg);
 return $self->{merge_fos_list};
}

sub fa {
 my ( $self, $arg ) = @_;
   $self->{fa} = $arg if ( $arg);
 return $self->{fa};
}

sub ca {
 my ( $self, $arg ) = @_;
   $self->{ca} = $arg if ( $arg);
 return $self->{ca};
}
=cut
return 1;