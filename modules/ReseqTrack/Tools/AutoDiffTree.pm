package ReseqTrack::Tools::AutoDiffTree;

use strict;
#use warnings;
#use Data::Dumper;
use ReseqTrack::Tools::Exception;
use ReseqTrack::Tools::GeneralUtils;
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename;

sub new {

  my ( $class, @args ) = @_;
  my $self = {};
  bless $self, $class;

  my ( $old_tree_file, $new_tree_file, $changelog, 
       $staging_dir, $verbose,$date, $changelog_header ) = rearrange(
					    [
					     qw(
                          OLD_TREE_FILE
			  NEW_TREE_FILE
			  CHANGELOG
			  STAGING_DIR
			  VERBOSE
                          DATE
                          CHANGELOG_HEADER
			  )
					    ],
					    @args
					   );

  $self->verbose($verbose);
  $self->changelog($changelog);

  if ($changelog_header){ # formatted time stamp for CHANGELOG file
      $self->changelog_header($changelog_header);
    }
  else{
    $self->create_timestamps();
  }
 
  #Should be getting time stamps from run_tree_for_ftp.pl 
  if ($date){
    $self->details_date($date);
  }else{
    my $tmp = current_date;
    $self->details_date($tmp);
  }  

  $self->staging_dir($staging_dir);
	

  $self->new_tree( ( $self->get_tree_hash($new_tree_file) ) );
  $self->old_tree( ( $self->get_tree_hash($old_tree_file) ) );
  $self->assign_g1k_types( $self->new_tree() );
  $self->assign_g1k_types( $self->old_tree() );

  $self->delete_identical( $self->new_tree, $self->old_tree );

  $self->check_no_md5_files ( $self->new_tree, $self->old_tree );
  
	
  $self->flag_moved_replaced_withdrawn();
  $self->flag_new();      
 
  
  throw "staging_dir not set" if ( ! $self->staging_dir);
  
  return ($self);
}

sub verbose{
  my ($self, $arg) = @_;
  if(defined($arg)){
    $self->{'verbose'} = $arg;
  }
  return $self->{'verbose'};
}

sub check_no_md5_files{
  my ( $self, $new, $old ) = @_;

  my $bad_md5 ="NO_MD5_IN_DB";
  my $something_wrong = 0;
  my @bad_old;
  my @bad_new;


  foreach my $f ( keys %$old) {
    if ( $$old{$f}{md5} eq $bad_md5   ) {
      print "$f\n" if($self->verbose);
      $something_wrong++;
      push (@bad_old, $f);
    }
  }

  foreach my $f ( keys %$new) {
    if ( $$new{$f}{md5} eq $bad_md5   ) {
      print "$f\n" if($self->verbose);
      $something_wrong++;
      push (@bad_new, $f);
    }
  }



  return if (!$something_wrong);

  print STDERR "==========================================\n";
  print STDERR "There can be files on ftp site that are not in database\n";
  print STDERR "as long as there is a matching entry in both tree files\n";
  print STDERR "they can be ignored.\n";

  print STDERR "\nList of files in old/new tree files with no md5\n";
  print STDERR "that script cannot account for\n";

  print STDERR "\nFound $something_wrong file(s) on ftp\n";
  print STDERR " with no explainable reason for existence \n";
  
  print STDERR "Will proceed but will ignore these files\n";
  print STDERR "=============================\n";

  #to prevent confusion downstream. 
  #Totally skip bad entries via the "delete em" route.

  print STDERR "File(s) from old tree file\n" if (@bad_old);
  foreach my $of (@bad_old) {
    print STDERR  "Ignoring:$of\n";
    delete $$old{$of};
  }

  print STDERR "\nFile(s) from new tree file\n" if (@bad_new);
  foreach my $nf (@bad_new) {
    print STDERR "Ignoring $nf\n";
    delete $$new{$nf};
  }
  
  print STDERR "==========================================\n";
  return;

}

sub create_log_files {
  my $self = shift;

  print "Creating log files\n\n" if($self->verbose);


  my $new = $self->new_tree();
  my $old = $self->old_tree();

  $self->output_changelog_files( $new, "new" );
  $self->output_changelog_files( $old, "withdrawn" );
  $self->output_changelog_files( $old, "replacement" );
  $self->output_moved_changelog_files();
  if ($self->files_to_archive) {
    $self->amend_CHANGELOG;
  } else {
    #print "No obvious changes detected , not amending CHANGELOG\n";
  }
  
  $self->files_to_archive_hash_to_array();
  $self->change_permissions(); 
  #print "Finished creating log files\n";

  return;
}

sub change_permissions {

  my $self = shift;
  my $log_files         = $self->files_to_archive_array();
  
  foreach my $i ( @$log_files) {
    chmod ( 0775, $i);
  }
  

  return;
}




sub amend_CHANGELOG {
  my $self = shift;
  my @changes;
  
  # print "------------------------\n";
  push( @changes, $self->changelog_header);
  push( @changes,"\n\n");

  my $all_modifications = $self->all_modifications();
  my $log_files         = $self->log_files();



  foreach my $key ( keys %$all_modifications  ) {
    my $line;
    $line .=  "Modification to: ";
    foreach my $i (@{$$all_modifications{$key}}) {    
      $line .= lc ($i). ",";
    }
    $line =~ s/\,$//;
    $line =~ s/\/$//;
    $line .= "\n";
    push( @changes,$line);
    push( @changes,  "\nDetails can be found in :\n");

    my $file = $$log_files{$key};
    my @aa =split /\//,$file;
    $file = $aa[-2] . "\/". $aa[-1];

 
    push( @changes,  $file);	
    push( @changes,"\n\n");
  }
   

  open my $IN, '<', $self->changelog || throw "open CHANGELOG failed";
  my @bot = <$IN>;
  close($IN);
 
  my $amended_changelog = $self->staging_dir . "CHANGELOG";
  $self->files_to_archive($amended_changelog);
  open my $OUT, '>', "$amended_changelog" || throw "no out";
  print $OUT @changes;
  print $OUT @bot;
  close($OUT);
    
  return;
}


sub flag_new {
  my $self = shift;
  my $new  = $self->new_tree;

  foreach my $f ( keys %$new ) {	
    if ( defined ($$new{$f}{change}) ) {
      #print "skipping $f ",$$new{$f}{change},"\n";
      next;
    }
		
    $$new{$f}{change} = "new";
    #print "new:: $$new{$f}{change}\n";
    $self->modified_type_hash( $$new{$f}{g1k_type}, "new" );		
  }

  return;
}


sub output_changelog_files {
  my ( $self, $hash, $action ) = @_;


  my $all_modifications = $self->all_modifications();
  my $altered_types     = $self->altered_types();
  my $timestamp         = $self->details_date();
  my $log_files         = $self->log_files();
  my $out;

 

  foreach my $key ( keys %$altered_types ) {
    if ( ($key =~ /CHANGELOG/) && (defined $$altered_types{$key}{$action}) ) {
      #print "Skipping change for $key $action\n";
      next;
    }

    if ( defined $$altered_types{$key}{$action} ) {
      printf "%-20s", $key if($self->verbose);
      printf "%-15s", "$action: " if($self->verbose);
      printf "%6s",   $$altered_types{$key}{$action} if($self->verbose);

      $out = $self->staging_dir(). "changelog_details/changelog_details_"
	. $timestamp . '_'  . $action;
    

      print "\n$out\n" if($self->verbose);
      
      #    $$log_files{$key}{$action} = $out;
      $$log_files{$action} = $out;
      open my $OUT, '>>', $out || throw "No on open $out";

      $self->files_to_archive($out);

      foreach my $f ( keys %$hash ) {

				#print "$f\n ";
	my $h = $$hash{$f};

	my $change = $$hash{$f}{change};
	if ( !defined($change) ) {
	  #print $f, "\n";
	  #print Dumper($h);
	  throw "ERROR. Residual file in tree hash with no change\n";
	}
				#print $change, "\n";
	next if ( $change ne "$action" );
	next if $$hash{$f}{g1k_type} ne "$key";

	my $i_type = $$hash{$f}{g1k_type};
	my @aa;
	if ( defined( $$all_modifications{$action} ) ) {
	  @aa = grep /$i_type/, @{ $$all_modifications{$action} };
	}
	if ( $#aa == -1 ) {

	  # print "adding $action g1k_type\n";
	  push(
	       @{ $$all_modifications{$action} },
	       $$hash{$f}{g1k_type} );
	}
	$f =~ s/^ftp\///;

	print $OUT $f, "\n";
      }

      close($OUT);

    }
  }

  $self->log_files($log_files);
  $self->all_modifications( $all_modifications);

  return;
}

sub all_modifications {
  my ( $self, $arg, $out ) = @_;
  if ($arg) {
    $self->{all_modifications} = $arg;
  }
  return $self->{all_modifications};

}



sub output_moved_changelog_files {
  my $self              = shift;
  my $new               = $self->new_tree();
  my $old               = $self->old_tree();
  my $verbose           = $self->verbose;
  my $altered_types     = $self->altered_types();
  my $timestamp         = $self->details_date();
  my $all_modifications = $self->all_modifications();
  my $log_files         = $self->log_files();
  my $out;

  foreach my $key ( keys %$altered_types ) {

    if ( defined $$altered_types{$key}{moved} ) {
      printf "%-20s\t", $key if($self->verbose);
      printf "%12s",    "moved: " . $$altered_types{$key}{moved} if($self->verbose);

      $out = $self->staging_dir() .
	"/changelog_details/changelog_details_" . $timestamp . "_moved";

      $out =~ s/\/\//\//;

      print "\t$out\n" if($self->verbose);
      #     $$log_files{$key}{log} = $out;
      $$log_files{moved} = $out;
      open my $OUT, '>>', $out || throw "No on open $out";
      $self->files_to_archive($out);
 
      foreach my $f ( keys %$new ) {

	#print $f, "\n";
	#print $$new{$f}{g1k_type}, "\n";
	next if ( ( $$new{$f}{type} ) eq ("directory") );
	next if ( $$new{$f}{change} ne "moved" );
	next if $$new{$f}{g1k_type} ne "$key";

	my $i_type = $$new{$f}{g1k_type};
	my @aa;
	if ( defined( $$all_modifications{moved} ) ) {
	  @aa = grep /$i_type/, @{ $$all_modifications{moved} };
	}
	if ( $#aa == -1 ) {

	  #  print "adding moved g1k_type\n";
	  push( @{ $$all_modifications{moved} },
		$$new{$f}{g1k_type} );
	}

	#   print Dumper( $new{$f});
	#   print $f,"\t", $new{$f}{change},"\t", $new{$f}{g1k_type},"\n";

	my $old = $$new{$f}{old_location};
	my $new = $$new{$f}{new_location};

	$old =~ s/^ftp\///;
	$new =~ s/^ftp\///;

	print $OUT $old, "\t";
	print $OUT $new, "\n";

      }
      close($OUT);

    }
  }
  $self->all_modifications( $all_modifications);
  $self->log_files($log_files);
  return;
}



sub flag_moved_replaced_withdrawn {
  my ($self) = shift;
  print "Starting Comparison\n\n" if($self->verbose);

  my $new     = $self->new_tree();
  my $old     = $self->old_tree();
  my $verbose = $self->verbose;

  foreach my $of ( keys %$old ) {	
    next if $self->was_replaced( $old, $new, $of, $verbose );
    next if $self->was_moved( $old, $new, $of );
    next if $self->was_withdrawn( $old, $new, $of );
  }

  return;
}

sub was_moved {
  my $self     = shift;
  my $old_hash = shift;
  my $new_hash = shift;
  my $f_name   = shift;

  my $old_md5  = $$old_hash{$f_name}->{md5};
  my $old_base = basename($f_name);
  my $verbose  = $self->verbose();

  return 0 if $$old_hash{$f_name}->{type} eq "directory";

  foreach my $key ( keys %{$new_hash} ) {

    #same md5
    if ( $old_md5 eq $$new_hash{$key}->{md5} ) {
      next if ( $old_md5 eq "d41d8cd98f00b204e9800998ecf8427e" );


      return if ($f_name eq $key);   
      next   if ( ( $f_name =~ /bai$/ ) && ( $key =~ /bai$/ ) );
      next if ( $$new_hash{$key}{size} == 0 );

      my $new_base = basename($key);

      if ( $new_base eq $old_base ) {

	print "++ moved:\n:$f_name:\n:$key:\n\n" if($self->verbose); # if $verbose;

	$$old_hash{$f_name}{change} = "moved";
	$$new_hash{$key}{change}    = "moved";


	$$old_hash{$f_name}{new_location} = $key;
	$$old_hash{$f_name}{old_location} = $f_name;

	$$new_hash{$key}{new_location} = $key;
	$$new_hash{$key}{old_location} = $f_name;


	$self->modified_type_hash( $$new_hash{$key}{g1k_type},
				   $$new_hash{$key}{change} );

	$self->modified_type_hash( $$old_hash{$f_name}{g1k_type},
				   $$old_hash{$f_name}{change} );
	return 1;
      }

      if ( $new_base ne $old_base ) {

	#   print "--was_moved or renamed ???\n";
	#   print "Assuming renamed = moved\n";
	#   print $f_name ,"\n";
	#   print $key ,"\n";
	print
	  "Same md5 for (flagged as moved):\n$old_base\n$new_base\n" if($self->verbose);

				# print "Same md5 for\n";

	$$old_hash{$f_name}{change} = "moved";
	$$new_hash{$key}{change}    = "moved";

	$$new_hash{$key}{new_location} = $key;
	$$new_hash{$key}{old_location} = $f_name;

	$self->modified_type_hash(
				  $$old_hash{$f_name}{g1k_type},
				  $$old_hash{$f_name}{change}
				 );
	$self->modified_type_hash( 
				  $$new_hash{$key}{g1k_type},
				  $$new_hash{$key}{change} );

				#  print "Do not know for this case\n";
	my $x = $$old_hash{$f_name};
	print Dumper ($x) if $self->verbose;
	my $y = $$new_hash{$key};
	print Dumper($y) if $self->verbose;
	return 1;
      }

    }

  }

  return 0;
}

sub altered_types {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{altered_types} = $arg;
  }
  return $self->{altered_types};
}

sub modified_type_hash {
  my ( $self, $type, $action,$debug ) = @_;

  my $altered_types = $self->altered_types();

  print "--$action --$type\n" if $self->verbose;

  if ( !($action) || !($type) ) {
    throw("Failed action $action or type $type");

  }

  if ( !defined $$altered_types{$type} ) {
    $$altered_types{$type}{$action} = 1;

  } else {
    $$altered_types{$type}{$action} += 1;
  }

  print $$altered_types{$type}{$action}, "\n" if $self->verbose;

  # print "-------\n";
  $self->altered_types($altered_types);
  return;
}

sub was_withdrawn {

  #1 old filename not in new hash
  #2 old filename does not have a md5 matching
  #  anything in new hash ( skip bai, bas, 0 size files);
  #3 maybe print warning about 0 size files
  my $self     = shift;
  my $old_hash = shift;
  my $new_hash = shift;
  my $f_name   = shift;
  my $ignore   = 0;

  my $k = keys(%$new_hash);

  if ( $k == 0 ) {
    $$old_hash{$f_name}{change} = "withdrawn";
    $self->modified_type_hash( $$old_hash{$f_name}{g1k_type},
			       $$old_hash{$f_name}{change} );
    return 1;
  }

  if ( exists $new_hash->{$f_name} ) {
    return 0;			# not withdrawn
  }

  my $old_md5 = $$old_hash{$f_name}->{md5};

  foreach my $key ( keys %{$new_hash} ) {

    #    print Dumper $$new_hash{$key};
    #    print $key, "$old_md5 ==== ",$$new_hash{$key}->{md5}."\n";

    if ( $old_md5 eq $$new_hash{$key}->{md5} ) {

      print "Have matching md5 \n" if($self->verbose);
      print $key, "\n$old_md5", $$new_hash{$key}->{md5}, "\n" if($self->verbose);

      $ignore = 1 if ( $$old_hash{$f_name}{g1k_type} =~ /BAI/ );
      $ignore = 1 if ( $$old_hash{$f_name}{g1k_type} =~ /BAS/ );
      $ignore = 1 if ( $$old_hash{$f_name}{size} == 0 );

    }

    if ( $ignore == 0 ) {
      print "withdrawn:\n$f_name   ", $$old_hash{$f_name}{g1k_type},
	"\n\n" if($self->verbose);
      $$old_hash{$f_name}{change} = "withdrawn";

      $self->modified_type_hash( $$old_hash{$f_name}{g1k_type},
				 $$old_hash{$f_name}{change} );

      return 1;
    }

  }

  return 0;
}

###############

sub was_replaced {
  my $self     = shift;
  my $old_hash = shift;
  my $new_hash = shift;
  my $f_name   = shift;
  

  if ( exists $new_hash->{$f_name} ) {
    print "$f_name in new hash\n" if $self->verbose;

    if ( $$new_hash{$f_name}->{md5} ne $$old_hash{$f_name}->{md5} ) {
      print "replaced:$f_name\n"  if $self->verbose;
      print "old md5:"
	. $$new_hash{$f_name}->{md5}
	  . "\nnew md5:"
	    . $$old_hash{$f_name}->{md5}, "\n\n" if $self->verbose;

      $$old_hash{$f_name}{change} = "replacement";
      $$new_hash{$f_name}{change} = "replacement";

      $self->modified_type_hash( $$old_hash{$f_name}{g1k_type},
				 $$old_hash{$f_name}{change} );

      $self->modified_type_hash (
				 $$new_hash{$f_name}{g1k_type},
				 $$new_hash{$f_name}{change}
				);

      return 1;
    }
  } else {

    return 0;
  }

  return 0;


}

sub delete_identical {
  my ( $self, $new, $old ) = @_;

  # Doing it the hard way. Remove identical entries from both
  my $indentical = 0;
  my @clear;
  my $dates_md5_same = 0;
  #print "Deleting identical entries\n";
  foreach my $f ( keys %$old) {
    my $m = 0;
    my $date_match = 0;
    if ( exists $$new{$f} ) {
      $m++ if $$old{$f}{type} eq $$new{$f}{type};
      $m++ if $$old{$f}{size} eq $$new{$f}{size};
      $date_match = 1 if $$old{$f}{date} ne $$new{$f}{date};
      $m++ if $$old{$f}{md5}  eq $$new{$f}{md5};
    }

    if ( $m == 3 ) {
      $indentical++;
      push( @clear, $f );


      if ( $date_match == 1  ) {
	$dates_md5_same++;
      }

      if ( $$old{$f}{md5} =~/NO_MD5_IN_DB/) {
	#	print "$f has no md5 but is same size. Ignoring for changelogs\n";
      }

    }
  }

  foreach (@clear) {
    delete $$old{$_};
    delete $$new{$_};
  }

  print "Total indentical files = $indentical\n" if($self->verbose);
  my $old_keys = keys(%$old);
  print "Have $old_keys keys in old tree hash\n" if($self->verbose);
  my $new_keys = keys(%$new);
  print "Have $new_keys keys in new tree hash\n" if($self->verbose);

  print "Have $dates_md5_same files with same md5 but diff dates\n" if($self->verbose);


  if ( $old_keys == 0 && $new_keys == 0 ) {
    print "No changes to tree files occurred\n" if($self->verbose);
  }



  $old_keys = keys(%$old);
  $new_keys = keys (%$new);
 foreach my $of ( keys %$old) {
	print $of,"\n" if ($old_keys < 10 && $self->verbose);
 }
 
 foreach my $nf ( keys %$new) {
        print $nf,"\n"  if ($new_keys < 10 && $self->verbose);
 }

#  print Dumper ($old) if ( $old_keys < 11) ;
#  print Dumper ($new) if ( $new_keys < 11) ;


  print "\n\n" if($self->verbose);
  return 0;
}

sub assign_g1k_types {

  my ( $self, $tree_hash ) = @_;

  foreach my $key ( keys %$tree_hash ) {

    $$tree_hash{$key}{g1k_type} = assign_type($key);
    if ( $$tree_hash{$key}{g1k_type} =~ /BAI|BAS/ ) {
      $$tree_hash{$key}{g1k_type} = "BAM";
    }

  }

  return 0;
}

sub assign_type {
  my $name = shift;
  my $type = "FRED";

  if ( $name =~ /reference/i ) {
    $type = "REFERENCE";
  } elsif ( $name =~ /\/technical\/working/i ) {
    $type = "INTERNAL";
  } elsif ( $name =~ /changelog/i ) {
    $type = "CHANGELOG";
  } elsif ( $name =~ /\/technical\//i   && $name !~ /\/ncbi_varpipe_data\//i) {
    $type = "TECHNICAL";
  } elsif ( $name =~ /\/release/i ) {
    $type = "RELEASE";
  } elsif ( $name =~ /\.filt\.fastq\.gz$/ ) {
    $type = "FILTERED_FASTQ";
  } elsif ( $name =~ /sequence_staging/ ) {
    $type = "ARCHIVE_FASTQ";
  } elsif ( $name =~ /\.fastq\.gz/i ) {
    $type = "FASTQ";
  } elsif ( $name =~ /\.bam$/ ) {
    $type = "BAM";
  } elsif ( $name =~ /\.bai$/ ) {
    $type = "BAI";
  } elsif ( $name =~ /\.bas$/ && !( $name =~ /README/ ) ) {
    $type = "BAS";
  } elsif ( $name =~ /\.index/ ) {
    $type = "INDEX";
  } else {
    $type = "MISC";
  }
  if ( $name =~ /\/pilot3_/ ) {
    my $tmp = $type;
    $tmp  = "PILOT3_" . $type;
    $type = $tmp;
  } elsif ( $name =~ /\/pilot_data\// ) {
    my $tmp = $type;
    $type = "PILOT_" . $tmp;
    $type = $tmp;
  }
  if ( $name =~ /\/withdrawn\// ) {
    my $tmp = $type;
    $tmp  = "WITHDRAWN_" . $type;
    $type = $tmp;
  }
 if ($name =~ /mosaik/i ){
      my $tmp = $type;
      $tmp = "NCBI_".$type;
      $type =  $tmp;
    }
  exit("WHAT ???") if ( $type =~ /FRED/ );

  return $type;

}

###########################

sub get_tree_hash {
  my ( $self, $file, $verbose ) = @_;

  my %hash        = ();
  my $bad_entries = 0;
  my $line        = 0;
  print "Processing $file\n" if($self->verbose);
  open( my $IN, "<", "$file" ) || throw "\n\nNo open :bad filehandle $file\n\n";

  while (<$IN>) {
    $line++;

    # print;
    chomp;

    #$_ =~ s/\//\_/g;
    my @a = split /\t/;

    next if ( $a[1] eq "directory" );

    my $cols = scalar(@a);
    if ( $cols != 5 ) {
      if ($verbose) {
	print "$a[0]:Wrong values ($cols) expected 5. line $line\n";
      }

      $bad_entries++;
    }

    $a[0] =~ s/\/ftp\///g;

    $hash{ $a[0] }{type} = $a[1];
    $hash{ $a[0] }{size} = $a[2];
    $hash{ $a[0] }{date} = $a[3];

    if ($a[4]) {
      $a[4] =~ s/\s+//g;
      $hash{ $a[0] }{md5} = $a[4];
    } else {
      $hash{ $a[0] }{md5} = "NO_MD5_IN_DB";
      print "missing md5: $a[0] in tree file\n" if($self->verbose);
    }
  }
  my $k = keys(%hash);

  if ( $k == 0 ) {
    throw "ERROR:No keys in hash for $file\n";
  }

  print "Have $k keys in $file hash\n" if($self->verbose);
  print "Also have $bad_entries bad _entries\n\n" if($self->verbose);
  return \%hash;
}


sub old_tree_file {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{old_tree_file} = $arg;
  }
  return $self->{old_tree_file};
}

sub new_tree_file {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{new_tree_file} = $arg;
  }
  return $self->{new_tree_file};
}

sub new_tree {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{new_tree} = $arg;
  }
  return $self->{new_tree};
}

sub old_tree {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{old_tree} = $arg;
  }
  return $self->{old_tree};
}

sub changelog_header {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{changelog_header} = $arg;
  }
  return $self->{changelog_header};
}

sub details_date {
  my ( $self, $arg ) = @_;

  if ($arg) {
    $self->{details_date} = $arg;
  }

  return $self->{details_date};
}

sub create_timestamps {
  my ($self,$arg) = shift;

  my $time = current_time;
 
  my @aa = split /:|\s+/,$time;
 
  $self->changelog_header($aa[0]);
 
  return;
}

sub changelog {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{changelog} = $arg;
  }
  return $self->{changelog};
}



sub log_files {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{log_files} = $arg;
  }
  return $self->{log_files};

}

sub files_to_archive {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{files_to_archive}{$arg} =1;
  }
  return $self->{files_to_archive};

}
sub staging_dir {
  my ( $self, $arg ) = @_;
  if ($arg) {
    if (! ($arg =~ /\/$/)) {
      $arg .= "\/";
      $arg =~ s/\/\//\//;
    }
    $self->{staging_dir} = $arg;
  }
  return $self->{staging_dir};

}

sub files_to_archive_hash_to_array{
  my  $self =shift;

  my $hash = $self->files_to_archive;
  my $files= ();

  foreach my $key (keys %$hash) {
    print $key,"\n" if($self->verbose);
    push (@$files,$key);
  }



  $self->files_to_archive_array($files);
  return;
}

sub files_to_archive_array{
  my ( $self, $arg ) = @_;
  
  if ($arg) {
    $self->{files_to_archive_array}= $arg;
  }
 
  return $self->{files_to_archive_array};
}

1;
