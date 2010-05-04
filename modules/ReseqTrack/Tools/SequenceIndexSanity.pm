package ReseqTrack::Tools::SequenceIndexSanity;

use strict;

use Exporter;
use File::Copy;
use File::Basename;
use Getopt::Long;
use File::Find ();
use ReseqTrack::Tools::Intersection;
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::FileSystemUtils qw(list_files_in_dir);
use ReseqTrack::Tools::SequenceIndexUtils qw(get_index_hash_on_column get_index_hash);
use ReseqTrack::DBSQL::ERADBAdaptor;
use ReseqTrack::Tools::RunMetaInfoUtils qw(index_to_name_hash);

use vars qw (@ISA  @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(
      check_column_count_sanity
      compare_index_to_filesystem
      compare_index_to_ftp
      has_empty_columns
      has_null_strings
      has_trailing_spaces
      check_column_syntax
      check_index_file_existance_to_ftp
      check_paired_files
      compare_index_to_era
      check_analysis_group
      check_population
      check_header
    );


sub check_header{
   my ($file) = @_;
   open(FH, $file) or throw("IndexUtils:check_column_count_sanity failed to open ".
                            $file." $!");
   my $count = 0;
   while(<FH>){
     chomp;
     if($count == 0){
       unless(/SUBMISSION/){
         print STDERR $file." doesn't have an index header it first line is\n";
         print $_."\n";
         $count++;
         next;
       }else{
         $count++;
       }
     }else{
       $count++;
     }
     last unless($count == 0);
   }
   close(FH);
}

sub check_column_count_sanity{
  my ($file) = @_;
  open(FH, $file) or throw("IndexUtils:check_column_count_sanity failed to open ".
                            $file." $!");
  my %sanity_hash;
  while(<FH>){
    next if(/SUBMISSION_ID/i);
    my @values = split /\t/, $_;
    if(@values ne 26){
      $sanity_hash{$values[0]} = "Wrong column number should have 26 entries not ".
          scalar(@values);
    }
  }
  close(FH);
  return \%sanity_hash;
}

sub compare_index_to_filesystem{
  my ($file, $dir) = @_;
  open(FH, $file) or throw("IndexUtils:check_column_count_sanity failed to open ".
                           $file." $!");
  
  my %index_hash;
  my $dir_copy = $dir;
  while(<FH>){
    next if(/SUBMISSION_ID/i);
    chomp;
    my @values = split /\t/, $_;
    next if($values[20] == 1);
    if($dir_copy =~ /data/ && $values[0] =~ /data/){
      $dir_copy =~ s/data//;
    }
    $dir_copy .= "/" unless($dir_copy =~ /\/$/);
    my $full_path = $dir_copy.$values[0];
    $full_path =~ s/\/\//\//;
    unless($index_hash{$full_path}){
      $index_hash{$full_path} = 1;
    }else{
      warning($full_path." already existing in hash will not replace first entry");
    }
  }
  close(FH);
  my ($files, $dir_hash) = list_files_in_dir($dir, 1);
  my @index_list = keys(%index_hash);
  my $index_set = ReseqTrack::Tools::Intersection->new
      (
       -list => \@index_list,
      );
  my $dir_set = ReseqTrack::Tools::Intersection->new
      (
       -list => $files,
      );
  my $not_index = $dir_set->not($index_set);
  my $not_dir = $index_set->not($dir_set);
  return ($not_index, $not_dir);
}

sub check_index_file_existance_to_ftp{
	my ($file, $ftp_root) = @_;
	my $total_count=0;
	my $faulty_count=0;
	my %log;
	open(FH, $file) or throw("IndexUtils:check_column_count_sanity failed to open ".$file." $!");
    
  	while(<FH>){
  		$total_count++;
    	
    	next if(/SUBMISSION_ID/i);
    	chomp;
    	my @values = split /\t/, $_;
    	my $full_path = $ftp_root."/".$values[0];
  		unless(-e $full_path){
    		unless($values[20] == 1){
                  #print STDERR $full_path." doesn't exist\n";
                  $log{"Not in FTP:\t".$full_path} = 1;
                  #print $string."\n";
                  $faulty_count++;
                  
    		}
      		
  		}else{
                  if($values[20] == 1){
                    #print STDERR $full_path." exists when it should be withdrawn\n";
                    $log{"In FTP but withdrawn:\t".$full_path} = 1;
                    #print $string."\n";
                    $faulty_count++;
                    
                  } 
                  
                }
                
                
                #CHECK the PAIRED files, if any...
                if($values[18] eq "PAIRED" && $values[20] == 0)
                {
                  my $full_path = $ftp_root."/".$values[19];
                  unless(-e $full_path){
                    $log{"Not in FTP:\t".$full_path} = 1;
                  }
                }
  	}
	close(FH);
	return ($total_count, $faulty_count, %log);
        
}

sub type_file{
  my ($file) = @_;
  my $name = basename($file);
  my $type = -1;
  if($name =~ /[E|S]RR\d+\.fastq\.gz/i   || $name =~ /[E|S]RR\d+\.recal\.fastq\.gz/i || $name =~ /[E|S]RR\d+\.filt\.fastq\.gz/i ){
    #$frag = $file;
    $type=0;
  }elsif($name =~ /[E|S]RR\d+\_1\.fastq\.gz/i || $name =~ /[E|S]RR\d+\_1\.recal\.fastq\.gz/i || $name =~ /[E|S]RR\d+\_1\.filt\.fastq\.gz/i ){
    #$mate1 = $file;
    $type=1;
  }elsif($name =~ /[E|S]RR\d+\_2\.fastq\.gz/i ||  $name =~ /[E|S]RR\d+\_2\.recal\.fastq\.gz/i || $name =~ /[E|S]RR\d+\_2\.filt\.fastq\.gz/i){
    #$mate2 = $file;
    $type=2;
  }
  return $type;
}

sub compare_index_to_era{
  my( $dbhost,$dbname,$dbuser,$dbpass,$dbport,$index_file) = @_;
  my $db = ReseqTrack::DBSQL::ERADBAdaptor->new(
    -host => $dbhost,
    -user => $dbuser,
    -port => $dbport,
    -dbname => $dbname,
    -pass => $dbpass,
      );
  
  my $run_hash = get_index_hash_on_column($index_file, 2);
  
  my $g1k_index_adaptor = $db->get_ERARunMetaInfoAdaptor;
  my $indexes = $g1k_index_adaptor->fetch_all;
  
  my %era_index_hash;
  my %log;
  my $total_count=0;
  my $faulty_count=0;
  
  foreach my $index(@$indexes){
    $era_index_hash{$index->run_id} = $index;
  }
  
  my $index_to_name_hash = index_to_name_hash();
  
  foreach my $run(keys(%$run_hash)){
    my $lines = $run_hash->{$run};
    my $index = $era_index_hash{$run};
    if(!$index){
      #print STDERR $run." isn't found in the index table\n";
      next;
    }
    foreach my $line(@$lines){
      $total_count++;
      my @values = split /\t/, $line;
      for(my $i=0;$i < @values;$i++){
        
        my $name = $index_to_name_hash->{$i};
        next unless($i == 3);
        next if($i <= 1 || $i >= 19);
        next if($i == 7); #submission data
        next if($i == 12); #insturment platform
        next if($i == 15); #run name
        next if($i == 17); #insert size
        next if($i == 10); #population
        next if($i == 16); #run block name
        my $line_string = $values[$i];
        my $index_string = $index->$name;
        if($index_string && length($index_string) >= 1){
          if($line_string ne $index_string){
            #print $run." has a mismatch for ".$name." ".$i." ftp ".$line_string.
            " is not era ".$index_string."\n";
            $log{$run." has a mismatch for ".$name." ".$i." ftp ".$line_string." is not era ".$index_string} = 1;
            
            
          }
        }
      }
    }
  }
  $faulty_count = scalar(keys %log);  
  
  return ($total_count,$faulty_count,%log);
  

}

sub check_paired_files{
  my ($file, $ftp_root) = @_;
  my $total_count=0;
  my $faulty_count=0;
  my %first_pair =();
  my %second_pair =();
  my %index_hash = ();
  my %log= ();
  open(FH, $file) or throw("IndexUtils:check_column_count_sanity failed to open ".$file." $!");
  
  while(<FH>){
    
    
    next if(/SUBMISSION_ID/i);
    chomp;
    my @values = split /\t/, $_;
    if($values[20] == 1)
    	{
          #print "Withdrawn".$values[20]."\n";
          next;
    	}
    
    if($values[18] eq "PAIRED" && $values[20] == 0)
    {
      next unless($values[0] =~ /\_1\./ || $values[0] =~ /\_2\./);
      $total_count++;
      #FIRST PAIR 
      my $full_path_first = $ftp_root."/".$values[0];
      my $full_path_first_type = type_file($values[0]);
      my $full_path_second = $ftp_root."/".$values[19];
      my $full_path_second_type = type_file($values[19]);    	
      
      
      
      if($full_path_first_type == 0 && $full_path_second_type != -1)
      { 
        #print "Should second pair not be empty?: ".$full_path_first."\n"; 
        $faulty_count++;
        $log{"Second pair not empty:\t".$full_path_first} = 1;
          }elsif($full_path_first_type == 0 && $full_path_second_type== -1)
          {
            #good, no further action    		
          }elsif(($full_path_first_type == 1 && $full_path_second_type == 2) || ($full_path_first_type == 2 && $full_path_second_type == 1) )
          {
            #good
            $first_pair{$full_path_first} = 1;
            $second_pair{$full_path_second} = 1;
            $index_hash{$full_path_first} = $full_path_second;
            
          }
          else
          {
            #bad entry, report. 
            #print STDERR "Mismatch filetype 1 and 2: $full_path_first vs $full_path_second\n"; 
            $faulty_count++;
            $log{"Mismatch filetype PAIRED:\t".$full_path_first." vs ". $full_path_second} = 1;
          }
          
    	}
    	next; 
    	   	
  	
	}
	close(FH);
 
 
	#if file is in one pair it should also be in other pair.
	my @first_missing = grep ! exists $index_hash{$_}, values %index_hash;
	%index_hash = reverse %index_hash;
	my @second_missing = grep ! exists $index_hash{$_}, values %index_hash;
 

	#print STDERR "Files as Pair 2 but not as Pair 1: ".@first_missing."\n";
	foreach (@first_missing) {
		$faulty_count++;
		$log{"Files as Pair 2 but not as Pair 1:\t".$_}=1;
	}
	
	#my @second_missing = grep ! exists $second_pair{$_}, keys %first_pair;
	#print STDERR "Files as Pair 1 but not as Pair 2: ".@second_missing."\n";
	foreach (@second_missing) {
		$faulty_count++;
		$log{"Files as Pair 1 but not as Pair 2:\t".$_}=1;
	}
	
	
	return ($total_count, $faulty_count, %log);

}


sub compare_index_to_ftp{
  my ($file, $dir, $path_must_contain, $dir_must_contain) = @_;
  open(FH, $file) or throw("IndexUtils:check_column_count_sanity failed to open ".
                           $file." $!");
  
  my %index_hash;
  my $dir_copy = $dir;
  my $withdrawn_count = 0;
  $path_must_contain = 'sequence_read' if(!$path_must_contain);
  $dir_must_contain = 'ftp/data' unless($dir_must_contain);
  while(<FH>){
    #print;
    next if(/SUBMISSION_ID/i);
    chomp;
    my @values = split /\t/, $_;
    #print "Skipping ".$values[0]." as withdrawn\n" if($values[20] == 1);
    $withdrawn_count++ if($values[20] == 1);
    next if($values[20] == 1);

    #print "Comparing ".$values[0]." to ".$dir_copy."\n";

    #if($dir_copy =~ /data/ && $values[0] =~ /data/){
     # $dir_copy =~ s/data//;
    #}
    $dir_copy .= "/" unless($dir_copy =~ /\/$/);
    my $full_path = $dir_copy.$values[0];
    $full_path =~ s/\/\//\//;
    unless($index_hash{$full_path}){
      $index_hash{$full_path} = 1;
    }else{
      warning($full_path." already existing in hash will not replace first entry");
    }
  }
  close(FH);
  #print "There are ".$withdrawn_count." withdrawn files\n";
  my ($files, $dir_hash) = list_files_in_dir($dir, 1);
  my @ftp_files;
  foreach my $file(@$files){
    #print $file."\n";
    next unless($file =~ /$path_must_contain/);
    next unless($file =~ /$dir_must_contain/);
    #print "Adding ".$file." to ftp array\n";
    push(@ftp_files, $file);
  }
  my @index_list = keys(%index_hash);
  my $index_set = ReseqTrack::Tools::Intersection->new
      (
       -list => \@index_list,
      );
  my $dir_set = ReseqTrack::Tools::Intersection->new
      (
       -list => \@ftp_files,
      );
  my $not_index = $dir_set->not($index_set);
  my $not_dir = $index_set->not($dir_set);
  return ($not_index, $not_dir);
}

sub has_empty_columns{
  my ($file) = @_;
  open(FH, $file) or throw("IndexUtils:has_empty_columns failed to open ".
                            $file." $!");
  my %sanity_hash;
  while(<FH>){
    next if(/SUBMISSION_ID/i);
    chomp;
    my @values = split /\t/, $_;
    my $count = 0;
    foreach my $value(@values){
      if(!defined $value || $value =~ /^NULL$/i || $value eq '' || 
         $value =~ /^\s+$/){
        $sanity_hash{$values[0]} = [] if(!$sanity_hash{$values[0]});
        push(@{$sanity_hash{$values[0]}}, $count);
      }
      $count++;
    }
  }
  close(FH);
  return \%sanity_hash;
}

sub has_trailing_spaces{
  my ($file) = @_;
  open(FH, $file) or throw("IndexUtils:has_empty_columns failed to open ".
                            $file." $!");
  my %sanity_hash;
  while(<FH>){
    next if(/SUBMISSION_ID/i);
    chomp;
    my @values = split /\t/, $_;
  VALUE:for (my $i = 0; $i < @values; $i++){
    next VALUE if($i == 7);
    next VALUE if($i == 12);
    if($values[$i] =~ /\s+$/){
      unless($i == 12 || $i == 7){
        push(@{$sanity_hash{$values[0]}}, $i);
      }
    }
  }
  }
  close(FH);
  return \%sanity_hash;
}

sub has_null_strings{
  my ($file) = @_;
  open(FH, $file) or throw("IndexUtils:has_empty_columns failed to open ".
                            $file." $!");
  my %sanity_hash;
  while(<FH>){
    next if(/SUBMISSION_ID/i);
    chomp;
    my @values = split /\t/, $_;
  VALUE:for (my $i = 0; $i < @values; $i++){
    if($values[$i] =~ /^NULL$/i){
      push(@{$sanity_hash{$values[0]}}, $i);
    }
  }
  }
  close(FH);
  return \%sanity_hash;
}

sub check_analysis_group{
  my ($file) = @_;
  open(FH, $file) or throw("IndexUtils:has_empty_columns failed to open ".
                           $file." $!");
  my %sanity_hash;
  my %hash;
  $hash{'SRP000032'} = 'high coverage';
  $hash{'SRP000033'} = 'exon targetted';
  while(<FH>){
    chomp;
    next if(/SUBMISSION_ID/);
    my @values = split /\t/, $_;
    my $key = $values[0];
    if($values[25]){
      my $analysis_group = $hash{$values[3]};
      $analysis_group = 'low coverage' unless($analysis_group);
      my $problem = $key." has the incorrect analysis group ".
          $values[3]." should have ".$analysis_group." not ".$values[25];
      unless($analysis_group eq $values[25]){
        $sanity_hash{$key} = $problem;
      }
    }else{
      my $problem = $key." has no analysis group defined";
      $sanity_hash{$key} = $problem;
    }
  }
  close(FH);
  return \%sanity_hash;
}

sub check_population{
  my ($file) = @_;
  open(FH, $file) or throw("IndexUtils:has_empty_columns failed to open ".
                           $file." $!");
  my %sanity_hash;
  my %hash;
  $hash{'YRI'} = 1;
  $hash{'CHB'} = 1;
  $hash{'CHS'} = 1;
  $hash{'CHD'} = 1;
  $hash{'JPT'} = 1;
  $hash{'KHV'} = 1;
  $hash{'CEU'} = 1;
  $hash{'TSI'} = 1;
  $hash{'GBR'} = 1;
  $hash{'FIN'} = 1;
  $hash{'IBS'} = 1;
  $hash{'LWK'} = 1;
  $hash{'GWD'} = 1;
  $hash{'GHN'} = 1;
  $hash{'MAB'} = 1;
  $hash{'ASW'} = 1;
  $hash{'AJM'} = 1;
  $hash{'ACB'} = 1;
  $hash{'MXL'} = 1;
  $hash{'CLM'} = 1;
  $hash{'PEL'} = 1;
  $hash{'PUR'} = 1;
  while(<FH>){
    chomp;
    next if(/SUBMISSION_ID/);
    my @values = split /\t/, $_;
    my $key = $values[0];
    my $pop = $values[10];
    if($pop){
      unless($hash{$pop}){
        my $problem = $key." has a population which isn't recognised";
        $sanity_hash{$key} = $problem;
      }
    }else{
      my $problem = $key." has no population defined";
      $sanity_hash{$key} = $problem;
    }
  }
  close(FH);
  return \%sanity_hash;
}

sub check_column_syntax{
  my ($file) = @_;
  open(FH, $file) or throw("IndexUtils:has_empty_columns failed to open ".
                           $file." $!");
  my %sanity_hash;
  while(<FH>){
    next if(/SUBMISSION_ID/i);
    chomp;
    my @values = split /\t/, $_;
    VALUE:for (my $i = 0; $i < @values; $i++){
      my $key = $values[0];
      my $method_name = 'check_column_'.$i;
      no strict;
      unless(&$method_name($values[$i])){
        if(!$sanity_hash{$key}){
          $sanity_hash{$key} = {};
        }
        $sanity_hash{$key}->{$i} = $values[$i];
      }
      use strict;
    }
  }
  close(FH);
  return \%sanity_hash; 
}

sub check_column_0{
  my ($value) = @_;
  return 0 if(!$value);
  my $name = basename($value);
  unless($name =~ /^[S|E]RR.*\.fastq\.gz$/){
    return 0;
  }
  return 1;
}

sub check_column_1{
  my ($value) = @_;
  return 0 if(!$value);
  return 0 unless(length($value) == 32);
  return 0 unless($value =~ /^\S{32}$/);
  return 1;
}

sub check_column_2{
  my ($value) = @_;
  return 0 if(!$value);
  unless($value =~ /^[S|E]RR\d+$/){
    return 0;
  }
  return 1;
}

sub check_column_3{
  my ($value) = @_;
  return 0 if(!$value);
  return 0 unless($value =~ /^[S|E]RP\d+$/);
  return 1;
}


sub check_column_4{
  my ($value) = @_;
  return 1;
}

sub check_column_5{
  my ($value) = @_;
  return 0 if(!$value);
  return 1;
}

sub check_column_6{
  my ($value) = @_;
  return 0 if(!$value);
  return 0 unless($value =~ /^[S|E]RA\d+$/);
  return 1;
}

sub check_column_7{
  my ($value) = @_;
  if($value){
    return 0 unless($value =~ /^\d{2,4}.\d{1,2}.\d{1,2}/);
  }
  return 1;
}

sub check_column_8{
  my ($value) = @_;
  if($value){
    return 0 unless($value =~ /^[S|E]RS\d+$/);
  }
  return 1;
}

sub check_column_9{
  my ($value) = @_;
  return 0 if(!$value);
  return 0 unless($value =~ /^NA\d+$/ || $value =~ /^HG\d+$/);
  return 1;
}

sub check_column_10{
  my ($value) = @_;
  return 0 unless(length($value) == 3);
  return 0 unless($value =~ /^\w+$/);
  return 1;
}

sub check_column_11{
  my ($value) = @_;
  return 0 if(!$value);
  return 0 unless($value =~ /^[S|E]RX\d+$/);
  return 1;
}

sub check_column_12{
  my ($value) = @_;
  return 0 if(!$value);
  return 1;
}

sub check_column_13{
  my ($value) = @_;
  return 1;
}

sub check_column_14{
  my ($value) = @_;
  return 0 if(!$value);
  return 1;
}


sub check_column_15{
  my ($value) = @_;
  return 0 if(!$value);
  return 1;
}

sub check_column_16{
  my ($value) = @_;
  return 1;
}

sub check_column_17{
  my ($value) = @_;
  if($value){
    return 0 unless($value =~ /^\d+/);
  }
  return 1;
}

sub check_column_18{
  my ($value) = @_;
  return 0 if(!$value);
  return 0 unless($value eq 'PAIRED' || $value eq 'SINGLE');
  return 1;
}


sub check_column_19{
  my ($value) = @_;
  if($value){
    return check_column_0($value);
  }
  return 1;
}

sub check_column_20{
  my ($value) = @_;
  return 0 if(!(defined($value)) || $value eq '');
  return 0 unless($value == 0 || $value == 1);
  return 1;
}

sub check_column_21{
  my ($value) = @_;
  return 1;
}

sub check_column_22{
  return 1;
}

sub check_column_23{
  my ($value) = @_;
  return 0 unless($value);
  return 0 unless($value eq 'not available' || $value =~ /^\d+/);
  return 1;
}

sub check_column_24{
  my ($value) = @_;
  return 0 unless($value);
  return 0 unless($value eq 'not available' || $value =~ /^\d+/);
  return 1;
}

sub check_column_25{
  my ($value) = @_;
  my %analysis_group;
  $analysis_group{'low coverage'} = 1;
  $analysis_group{'high coverage'} = 1;
  $analysis_group{'exon targetted'} = 1;
  return 0 unless($value);
  return 0 unless($analysis_group{$value});
  return 1;
}

1;
