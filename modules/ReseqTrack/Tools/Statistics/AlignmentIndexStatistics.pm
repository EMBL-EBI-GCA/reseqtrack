=pod

=head1 NAME

ReseqTrack::Tools::Statistics::AlignmentIndexStatistics

=head1 SYNOPSIS

This is an object to provide some basic alignment statistics based on bas files

=head1 Example

my $alignment_index_stats = ReseqTrack::Tools::Statistics::AlignmentIndexStatistics
  ->new(
	-db => $db,
	-new_index => $new_index_file,
	-old_index => $old_index_file,
       );

$alignment_index_stats->make_stats;
$alignment_index_stats->print_stats;

=cut

package ReseqTrack::Tools::Statistics::AlignmentIndexStatistics;
#package AlignmentIndexStatistics;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::AlignmentIndexUtils;
use ReseqTrack::Tools::ERAUtils;
use ReseqTrack::Tools::GeneralUtils;
use ReseqTrack::Tools::Statistics;


@ISA = qw(ReseqTrack::Tools::Statistics);


=head2 new

  Arg [1]   : ReseqTrack::Tools::Statistics:AlignmentIndexStatistics
  Arg [NEW_INDEX]   : string, path to alignment.index.bas.gz file
  Arg [OLD_INDEX]   : string, path to alignment.index.bas.gz file
  Function  : create ReseqTrack::Tools::Statistics::AlignmentIndexStatistics object
  Returntype: ReseqTrack::Tools::Statistics:AlignmentIndexStatistics
  Exceptions: n/a
  Example   : my $stats = ReseqTrack::Tools::Statistics::AlignmentIndexStatistics
  ->new(
	-db => $db,
	-new_index => $new_index_bas_file,
	-old_index => $old_index_bas_file,
       );

=cut

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  unless($self->new_index && $self->old_index){
    $self->fetch_index_bas_file();
  }
  return $self;
}

#Accessor methods

=head2 new/old_index

  Arg [1]   : ReseqTrack::Tools::Statistics::AlignmentIndexStatistics
  Arg [2]   : string, filepath
  Function  : accessor method for index.bas.gz file paths
  Returntype: string
  Exceptions: throws if file doesn't exist or if filename doesn't match the
  expected patter YYYYMMDD.alignment.index
  Example   :

=cut


sub new_index{
  my ($self, $new_index) = @_;
  if($new_index){
    throw("AlignmentIndexStatistics:new_index ".$new_index." should exist")
      unless(-e $new_index);
    throw("AlignmentIndexStatistcs:new_index ".$new_index." must match pattern ".
	  "YYYYMMDD.alignment.index.bas.gz") unless($new_index =~ /\d+\.alignment\.index\.bas\.gz/ || $new_index =~ /\d+\.exome\.alignment\.index\.bas\.gz/);
    $self->{new_index} = $new_index;
  }
  return $self->{new_index};
}

sub old_index{
  my ($self, $old_index) = @_;
  if($old_index){
    throw("AlignmentIndexStatistics:old_index ".$old_index." should exist")
      unless(-e $old_index);
    throw("AlignmentIndexStatistcs:new_index ".$old_index." must match pattern ".
	  "YYYYMMDD.alignment.index.bas.gz") unless($old_index =~ /\d+\.alignment\.index\.bas\.gz/ || $old_index =~ /\d+\.exome\.alignment\.index\.bas\.gz/);
    $self->{old_index} = $old_index;
  }
  return $self->{old_index};
}


#Input methods

=head2 fetch_index_bas_file

  Arg [1]   : ReseqTrack::Tools::Statistics::AlignmentIndexStatistics
  Function  : fetch the most recent index file and the previous verion from
  the given database
  Returntype: n/a
  Exceptions: throw if unable to define a new_index file or if index file name
  doesn't match the expected patter
  Example   : $self->fetch_index_bas_file

=cut

sub fetch_index_bas_file{
  my ($self) = @_;
  my $fa = $self->db->get_FileAdaptor;
  #fetch all files of type MISC
  my $files = $fa->fetch_by_type('MISC');
  my %index_files;
  foreach my $file(@$files){

    next unless($file->filename =~ /index/);
    next unless($file->filename =~ /alignment/);
    next if($file->filename eq 'alignment.index');
    next if($file->filename =~ /^\d+\.alignment\.index$/);
    #print "file name is " . $file->filename . "\n";    
    my $new_date;
    if ($file->filename !~ /exome/i) {
    	$file->filename =~ /(\d+)\.alignment\.index\.bas\.gz/;
    	$new_date = $1;
    }
    else {
	$file->filename =~ /(\d+)\.exome\.alignment\.index\.bas\.gz/;
        $new_date = $1;
    }
    if(!$new_date){
      print "Can't parse ".$file->filename."\n";
    }
    #$index_files{$1} = $file->name;
    $index_files{$new_date} = $file->name;
  }
  
  #sort dates numerically, the newest date should always be the largest number
  #provided the YYYYMMDD format is followed

  my @dates = sort {$a <=> $b} keys(%index_files);
  my $new_date;
  if ($self->new_index && $self->new_index !~ /exome/i) {
    $self->new_index =~ /(\d+)\.alignment\.index/ ;
    $new_date = $1;
  }
  elsif ($self->new_index && $self->new_index =~ /exome/i) {
    $self->new_index =~ /(\d+)\.exome\.alignment\.index/ ;
    $new_date = $1;
  }
  my $old_date;
  #The newest date is the new file unless it is already defined
  $new_date = $dates[-1] unless($new_date);
  unless($self->new_index){
    $self->new_index($index_files{$new_date});
  }
  #Iterate through the other dates, once the new date is matched take the previous
  #date to defind the old index
  for(my $i=0; $i<@dates; $i++){
    next unless($new_date eq $dates[$i]);
    unless($i == 0){
      $old_date = $dates[$i-1];
    }
  }
  if($old_date && !$self->old_index){
    $self->old_index($index_files{$old_date});
  }
  throw("SequenceIndexStatistics:fetch_index_bas_file Failed to define a new index ".
	"file") unless($self->new_index);
	
  #print "current bas is $new_date and old bas is on $old_date\n";	
}



=head2 hashref_holders

  Arg [1]   : ReseqTrack::Tools::Statistics::AlignmentIndexStatistics
  Arg [2]   : hashref
  Function  : accessor methods for several hashes used to hold stats
  Returntype: hashref
  Exceptions: throws if not passed a hashref
  Example   : my $new_bp_cnt_hash = $self->new_bp_cnt_hash();

=cut


sub new_bp_cnt_hash{
  my ($self, $hashref) = @_;
  if($hashref){
    throw("Must pass new_per_sample a hashref not ".$hashref)
      unless(ref($hashref) eq 'HASH');
    $self->{new_bp_cnt_hash} = $hashref;
  }
  return $self->{new_bp_cnt_hash};
}

sub old_bp_cnt_hash{
  my ($self, $hashref) = @_;
  if($hashref){
    throw("Must pass new_per_sample a hashref not ".$hashref)
      unless(ref($hashref) eq 'HASH');
    $self->{old_bp_cnt_hash} = $hashref;
  }
  return $self->{old_bp_cnt_hash};
}

sub new_bp_cnt_by_pop_hash{
  my ($self, $hashref) = @_;
  if($hashref){
    throw("Must pass new_per_sample a hashref not ".$hashref)
      unless(ref($hashref) eq 'HASH');
    $self->{new_bp_cnt_by_pop_hash} = $hashref;
  }
  return $self->{new_bp_cnt_by_pop_hash};
}

sub old_bp_cnt_by_pop_hash{
  my ($self, $hashref) = @_;
  if($hashref){
    throw("Must pass new_per_sample a hashref not ".$hashref)
      unless(ref($hashref) eq 'HASH');
    $self->{old_bp_cnt_by_pop_hash} = $hashref;
  }
  return $self->{old_bp_cnt_by_pop_hash};
}

sub new_sample_cnt_hash{
  my ($self, $hashref) = @_;
  if($hashref){
    throw("Must pass new_per_sample a hashref not ".$hashref)
      unless(ref($hashref) eq 'HASH');
    $self->{new_sample_cnt_hash} = $hashref;
  }
  return $self->{new_sample_cnt_hash};
}

sub old_sample_cnt_hash{
  my ($self, $hashref) = @_;
  if($hashref){
    throw("Must pass old_per_sample a hashref not ".$hashref)
      unless(ref($hashref) eq 'HASH');
    $self->{old_sample_cnt_hash} = $hashref;
  }
  return $self->{old_sample_cnt_hash};
}

sub stats_hash{
  my ($self, $hashref) = @_;
  if($hashref){
    throw("Must pass stats a hashref not ".$hashref)
      unless(ref($hashref) eq 'HASH');
    $self->{stats} = $hashref;
    #$self->{stats_hash} = $hashref;
  }
  return $self->{stats_hash};
}

sub pop_hash {
  my ($self, $hashref) = @_;
  if($hashref){
    throw("Must pass stats a hashref not ".$hashref)
      unless(ref($hashref) eq 'HASH');
    $self->{pop_hash} = $hashref;
  }
  return $self->{pop_hash};
}

sub ind_pop_hash{
  my ($self, $hashref) = @_;
  if($hashref){
    throw("Must pass stats a hashref not ".$hashref)
      unless(ref($hashref) eq 'HASH');
    $self->{ind_pop_hash} = $hashref;
  }
  return $self->{ind_pop_hash};
}
  

=head2 make_stats

  Arg [1]   : ReseqTrack::Tools::Statistics::AlignmentIndexStatistics
  Function  : create the statistics given these specified index.bas.gz files
  Returntype: n/a
  Exceptions: n/a
  Example   : $self->make_stats;

=cut

sub make_stats{
  my ($self) = @_;
  my ($pops, $ind2pop) = $self->get_population;
  $self->pop_hash($pops);
  $self->ind_pop_hash($ind2pop);  
  my ($bp_cnt_hash, $bp_cnt_by_pop_hash, $sample_cnt_hash) = $self->parse_bas($self->new_index);
  $self->new_bp_cnt_hash($bp_cnt_hash);
  $self->new_bp_cnt_by_pop_hash($bp_cnt_by_pop_hash);
  $self->new_sample_cnt_hash($sample_cnt_hash);
  if($self->old_index){
    my ($old_bp_cnt_hash, $old_bp_cnt_by_pop_hash, $old_sample_cnt_hash) = 
      $self->parse_bas($self->old_index);
     $self->old_bp_cnt_hash($old_bp_cnt_hash);
 	 $self->old_bp_cnt_by_pop_hash($old_bp_cnt_by_pop_hash);
 	 $self->old_sample_cnt_hash($old_sample_cnt_hash);
  }

=head

 my %h = %$sample_cnt_hash;
 my $c = 1;
 foreach my $s ( keys %{$h{"SOLiD"}{"total"}} ) {
     print "sample $c is $s\n";
     $c++;
 }     

=cut

## FIXME: remove above lines after testing

}

=head2 get_population

  Arg [1]   : ReseqTrack::Tools::Statistics::AlignmentIndexStatistics
  Function  : get population information from run_meta_info table
  Returntype: two hashref
  Exceptions: n/a
  Example   : $self->get_population;

=cut

sub get_population {
	my ($self) = @_;	

	my $rmis = $self->rmis;
	
	my %populations;
	my %ind_to_pop;
	
	foreach my $rmi(@$rmis){
	  next unless($rmi->status eq 'public');
	  next if($rmi->study_id eq 'SRP000032' || $rmi->study_id eq 'SRP000033');
	  $populations{$rmi->population} = 1;
	  $ind_to_pop{$rmi->sample_name} = $rmi->population;
	}
	return (\%populations, \%ind_to_pop);
}	

=head2 parse_bas

  Arg [1]   : ReseqTrack::Tools::Statistics::AlignmentIndexStatistics
  Arg [2]   : string, path to bas file
  Function  : to parse the bas files to produce base pair statistics
  Returntype: array of hashrefs
  Exceptions: throw if a run doesn't have a base count
  Example   : my ($bp_cnt_hash, $bp_cnt_by_pop_hash, $sample_cnt_hash) = $self->parse_bas($self->new_index);

=cut


sub parse_bas {
	my ($self, $bas) = @_;
	my $hash = get_bas_hash_array($bas);
	my %bp_hash;
	my %sample_cnt_hash;
	my %bp_by_pop_hash;
	my $individual_to_pop = $self->ind_pop_hash;
	my $run_id_hash = $self->run_id_hash;
	foreach my $file (keys(%{$hash})){
		foreach my $line ( @{$hash->{$file}} ) {
			next if( $line =~ /bam_filename/i);
			my @values = split /\t/, $line;
			#print STDERR "line is $line\n";
			my $read_group = $values[6];
            #            next unless($run_id_hash->{$values[6]});  ### FIXME, CHECKME, for baylor solid BAMs, lots of read groups are 1, 0.1 etc. but in any event, run id should be unique when used properly
  			my $individual = $values[3];
  			my $population = $individual_to_pop->{$individual};
  			my $study = $values[2];
			
  			#next if ( $line =~ /20101123/ && $line =~ /chrom20/);  ## FIXME: change date!! This line produce the right results 
  														# this way chrom20 data are not double counted in alignment starting from 20101123 
  														# also chrom20 data can be slightly out of date
  			
  			next if ( $line =~ /chrom/ );
  			next if ( $line =~ /exon_targetted|high_coverage/i);
  			my $platform = $values[4];
  			if ($platform eq "illumina") {
  				$platform = "ILLUMINA";
  			}
  			my $mapped_bp = $values[8];
		    $bp_hash{$platform}->{$individual} += $mapped_bp;
		    $bp_hash{$platform}->{'total'} += $mapped_bp;
		    $bp_hash{"total"}->{'total'} += $mapped_bp;
			my $verbose = 0;
		  	if($verbose){
		    	print  "Adding ".$platform." ".$population." ".$read_group." ".$individual." ".$study." ".$mapped_bp.
		      	" to platform total ".convert_to_giga($bp_hash{$platform}->{'total'})." total total " . convert_to_giga($bp_hash{"total"}->{'total'}) ."\n";
		  	}
		  	$bp_by_pop_hash{$platform}->{$population} += $mapped_bp;
		  	$bp_by_pop_hash{total}->{$population} += $mapped_bp;	
		}
	}
 
	foreach my $platform(keys(%bp_hash)){
	  	foreach my $sample(keys(%{$bp_hash{$platform}})){
	  	    next if $sample eq "total"; ## FIXME, check me??
	    	my $gigabase = convert_to_giga($bp_hash{$platform}->{$sample});
	    	$sample_cnt_hash{$platform}{"total"}{$sample} = 1;
	    	$sample_cnt_hash{$platform}{"greater than 10"}{$sample} = 1 if($gigabase >= 10);
	    	$sample_cnt_hash{total}{"total"}{$sample} = 1; 
	    	$sample_cnt_hash{total}{"greater than 10"}{$sample} = 1 if($gigabase >= 10);  
 		 }
	}

	return (\%bp_hash, \%bp_by_pop_hash, \%sample_cnt_hash); 
}		


=head2 calculate_summary_stats

  Arg [1]   : ReseqTrack::Tools::Statistics::AlignmentIndexStatistics
  Function  : create diffs between given stats and populate a hash with the 
  statistics
  Returntype:hashref
  Exceptions: n/a
  Example   : $stats_hash = $self->calculate_summary_stats

=cut


sub calculate_summary_stats{
	my ($self) = @_;

	#get summary hashes
	my $new_hash = $self->new_bp_cnt_hash;
	my $new_hash_by_pop = $self -> new_bp_cnt_by_pop_hash;
	my $new_hash_sample_cnt = $self -> new_sample_cnt_hash;
	my ($old_hash, $old_hash_by_pop, $old_hash_sample_cnt);
	if($self->old_index){
		$old_hash = $self->old_bp_cnt_hash;
		$old_hash_by_pop = $self -> old_bp_cnt_by_pop_hash;
		$old_hash_sample_cnt = $self -> old_sample_cnt_hash;
	}
    
	my %stats_hash;

	#### for new bas  	
	my @platform = sort {$a cmp $b} keys(%$new_hash);
 	$stats_hash{new}{Date} = $self->get_bas_date($self->new_index);
	foreach my $platform(@platform){
    	next if($platform eq 'total');
    	my $mapped_giga_bp = convert_to_giga($$new_hash{$platform}->{'total'});
    	$stats_hash{new}{'mapped_giga_bp'}{$platform} = $mapped_giga_bp;
    	$stats_hash{new}{'num_individuals'}{$platform} = keys(%{$$new_hash_sample_cnt{$platform}->{'total'}});;
    	$stats_hash{new}{'num_ind_gt_10g'}{$platform} = keys(%{$$new_hash_sample_cnt{$platform}->{'greater than 10'}});
	}	
  
  	my $mapped_giga_bp = convert_to_giga($$new_hash{total}->{'total'});
  	$stats_hash{new}{'mapped_giga_bp'}{'Total'} = $mapped_giga_bp;
  	$stats_hash{new}{'num_individuals'}{'Total'} = keys(%{$$new_hash_sample_cnt{total}->{'total'}});
  	$stats_hash{new}{'num_ind_gt_10g'}{'Total'} = keys(%{$$new_hash_sample_cnt{total}->{'greater than 10'}});
    
    #print STDERR "Total ind gt 10g is " . keys(%{$$new_hash_sample_cnt{total}->{'greater than 10'}}) . "\n";
    
  	##### for old bas
  	if ($self->old_index) { 
	  	my @platform = sort {$a cmp $b} keys(%$old_hash);
	  	$stats_hash{old}{Date} = $self->get_bas_date($self->old_index);
		foreach my $platform(@platform){
	    	next if($platform eq 'total');
	    	my $mapped_giga_bp = convert_to_giga($$old_hash{$platform}->{'total'});
	    	$stats_hash{old}{'mapped_giga_bp'}{$platform} = $mapped_giga_bp;
	    	$stats_hash{old}{'num_individuals'}{$platform} = keys(%{$$old_hash_sample_cnt{$platform}->{'total'}});
	    	$stats_hash{old}{'num_ind_gt_10g'}{$platform} = keys(%{$$old_hash_sample_cnt{$platform}->{'greater than 10'}});
		}	
	 
		my $mapped_giga_bp_old = convert_to_giga($$old_hash{total}->{'total'});
		$stats_hash{old}{'mapped_giga_bp'}{'Total'} = $mapped_giga_bp_old;
		$stats_hash{old}{'num_individuals'}{'Total'} = keys(%{$$old_hash_sample_cnt{total}->{'total'}});
		$stats_hash{old}{'num_ind_gt_10g'}{'Total'} = keys(%{$$old_hash_sample_cnt{total}->{'greater than 10'}});
  	}
  	
	###### for diff
	$stats_hash{diff}{Date} = current_date;
	foreach my $pltf(@platform){
    	next if($pltf eq 'total');
		$stats_hash{diff}{'mapped_giga_bp'}{$pltf} = $stats_hash{new}{'mapped_giga_bp'}{$pltf} - $stats_hash{old}{'mapped_giga_bp'}{$pltf};
		$stats_hash{diff}{'num_individuals'}{$pltf} = $stats_hash{new}{'num_individuals'}{$pltf} - $stats_hash{old}{'num_individuals'}{$pltf};
  		$stats_hash{diff}{'num_ind_gt_10g'}{$pltf} = $stats_hash{new}{'num_ind_gt_10g'}{$pltf} - $stats_hash{old}{'num_ind_gt_10g'}{$pltf};
	}
	$stats_hash{diff}{'mapped_giga_bp'}{'Total'} = $stats_hash{new}{'mapped_giga_bp'}{'Total'} - $stats_hash{old}{'mapped_giga_bp'}{'Total'};
	$stats_hash{diff}{'num_individuals'}{'Total'} = $stats_hash{new}{'num_individuals'}{'Total'} - $stats_hash{old}{'num_individuals'}{'Total'};
  	$stats_hash{diff}{'num_ind_gt_10g'}{'Total'} = $stats_hash{new}{'num_ind_gt_10g'}{'Total'} - $stats_hash{old}{'num_ind_gt_10g'}{'Total'};
		
    ##### population breakdowns for each platform 

  	my $population = $self->pop_hash;	
	my @populations = sort {$a cmp $b} keys(%$population);

  	foreach my $platform(@platform){ #platform of the new_hash
    	next if($platform eq 'total');
    	my $count = convert_to_giga($$new_hash{$platform}->{'total'});
    	$stats_hash{new}{$platform}{'total_by_platform'} = $count;
    	foreach my $pop(@populations){
      		# print STDERR "Calling convert for ".$pop." ".$platform." with ".
      		#	$new_hash_by_pop{$platform}->{$pop}."\n";
      		my $bp_count = $$new_hash_by_pop{$platform}->{$pop};
      		my $value = convert_to_giga($bp_count);
      		$value = 0 unless($value);
      		$stats_hash{new}{$platform}{$pop}=$value;
    	}
	}
	
	## handle total count (lump up platforms)
  	my $total_count = $$new_hash{total}->{total};
  	my $total_cnt_in_gb = convert_to_giga($total_count);
  	$stats_hash{new}{'total_by_pop'}{'total_by_platform'}=$total_cnt_in_gb;
  	foreach my $pop(@populations){
    	my $bp_count = $$new_hash_by_pop{total}->{$pop};
    	my $value = convert_to_giga($bp_count);
    	$value = 0 unless($value);
     	$stats_hash{new}{'total_by_pop'}{$pop}=$value;
  	}
  	
  	## for old bas file
  	if ($self->old_index) {
	  	foreach my $platform(@platform){
	    	next if($platform eq 'total');
	    	my $count = convert_to_giga($$old_hash{$platform}->{'total'});
	    	$stats_hash{old}{$platform}{'total_by_platform'} = $count;
	    	foreach my $pop(@populations){
	      		my $bp_count = $$old_hash_by_pop{$platform}->{$pop};
	      		my $value = convert_to_giga($bp_count);
	      		$value = 0 unless($value);
	      		$stats_hash{old}{$platform}{$pop}=$value;
	    	}
		}
		## handle total counts (lump up platforms)
	  	my $total_count = $$old_hash{total}->{total};
	  	my $total_cnt_in_gb = convert_to_giga($total_count);
	  	$stats_hash{old}{'total_by_pop'}{'total_by_platform'}= $total_cnt_in_gb;
	  	foreach my $pop(@populations){
	    	my $bp_count = $$old_hash_by_pop{total}->{$pop};
	    	my $value = convert_to_giga($bp_count);
	    	$value = 0 unless($value);
	    	$stats_hash{old}{'total_by_pop'}{$pop}=$value;
	  	}
  	}
  	 
  $self->stats_hash(\%stats_hash);
  return \%stats_hash;
}


=head2 print_stats

  Arg [1]   : ReseqTrack::Tools::Statistics::AlignmentIndexStatistics
  Arg [2]   : hashref, hash of stats as produced by calculate summary stats
  Function  : to print to STDOUT a csv file with stats
  Returntype: n/a
  Exceptions: throws if doesn't not how to print a particular row of stats
  Example   : $self->print_stats();

=cut

sub print_stats{
  my ($self, $stats_hash) = @_;

  $stats_hash = $self->stats_hash unless($stats_hash);
  $stats_hash = $self->calculate_summary_stats() unless($stats_hash);

  my @headers = ('old', 'new', 'diff');
  my @rows = ('Date','mapped_giga_bp', 'num_individuals','num_ind_gt_10g');    
  my $type_hash = $self->row_type;
  print join(", ", "", @headers)."\n";
  
  my %print_text;
  my %print_text_by_pop;
  my @sorted_keys; #platform
  my @sorted_pops;
  
  foreach my $row(@rows){
    my $type = $type_hash->{$row};
    if($type eq 'SCALAR'){ #print Date
      print ", ";
      foreach my $header(@headers){
		my $value = $stats_hash->{$header}->{$row};
		throw($value. "This isn't a scalar it is a hash") if(ref($value) eq 'HASH');
		print $value.", ";
      }
      print "\n";
    }
	elsif($type eq 'HASH'){ #print stats, platform by platform
      my @keys = keys(%{$stats_hash->{new}->{$row}});
      @sorted_keys = sort {$a cmp $b} @keys;
	  foreach my $key(@sorted_keys){ #platforms
			my $string = "$row" . ", ";
			foreach my $header(@headers){
			  my $value = $stats_hash->{$header}->{$row}->{$key};
			  $value = 0 unless($value && $value >= 1);
			  $string .= $stats_hash->{$header}->{$row}->{$key}.", ";
			}
			$print_text{$key}{$row} = $string;
	    }
	}
	else{
		throw("Not sure how to print ".$row." ".$type);
	}
  }		
  
  foreach my $k (@sorted_keys){
		print "$k";
		foreach my $r (@rows) {
			print $print_text{$k}{$r}."\n";
		}
	}

	print "\n\nPopulation breakdown\n";

	foreach my $header(@headers){	
		next if $header eq "diff";
		foreach my $platform (  keys %{$stats_hash->{$header}} ) {
			my $type = $type_hash->{$platform};
			next if ( $type ); # skip stats_hash part that are NOT related to breakdown by population
			my $string = $platform . ", ";
			@sorted_pops = sort {$a cmp $b} keys%{$stats_hash->{$header}->{$platform}};
			foreach my $pop (@sorted_pops){ 
				if ( $stats_hash->{$header}->{$platform}->{$pop} >=1 ) {
					$string .=  $stats_hash->{$header}->{$platform}->{$pop} . ", ";
				}
				else {
					$string .= "0, ";
				}
			}			
			$print_text_by_pop{$platform} = $string;
		}#end of foreach platform
		print "\n$header - " . $stats_hash->{$header}->{Date} . "\n";
		print join(", ", "", @sorted_pops) . "\n";
		my @sorted_pf = sort {$a cmp $b} keys%print_text_by_pop;
		foreach my $pf (@sorted_pf) {
			print $print_text_by_pop{$pf} . "\n";
		}	
	}
}


=head2 row_type

  Arg [1]   : ReseqTrack::Tools::Statistics::AlignmentIndexStatistics
  Function  : return a hash which defined data types for particular rows of
  statistics, can be either Scalar or HASH
  Returntype: hashref
  Exceptions: 
  Example   : 

=cut

sub row_type{
  my ($self) = @_;
  my %hash;
  $hash{Date} = 'SCALAR';
  $hash{'mapped_giga_bp'} = 'HASH';
  $hash{'num_individuals'} = 'HASH';
  $hash{'num_ind_gt_10g'} = 'HASH';
  return \%hash;
}


=head2 get_bas_date

  Arg [1]   : ReseqTrack::Tools::Statistics::AlignmentIndexStatistics
  Arg [2]   : string, file path
  Function  : parses date from index file name assuming format
  YYYYMMDD.alignment.index
  Returntype: string, date string
  Exceptions: n/a
  Example   : my $date = $self->get_index_date($index_path);

=cut

sub get_bas_date{
  my ($self, $index) = @_;
  #print STDERR "Parsing ".$index."\n";
  $index =~ /(\d+)\.alignment\.index\.bas\.gz/;
  my $date = $1;
  #print STDERR "RETURNING ".$date."\n";
  return $date;
}

1;
