#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use IPC::System::Simple qw(system);
use IPC::Open2;

my ($vcf,
	$sample_panel,
	$user_input_pop_string,
	$region,
	$out_dir,
	$tabix,
	$vcftools_dir,
	$vcftools,
	$vcf_subset,
	$vcf_fill_an_ac,
	$no_tabix,
	$help,
	%result_hash,
  %selected_pop_samples,
);
	
&GetOptions( 
	'vcf=s'					=> \$vcf,
	'sample_panel=s'		=> \$sample_panel,
	'pop=s'					=> \$user_input_pop_string,
	'region=s'				=> \$region,
	'out_dir:s'				=> \$out_dir,
	'tabix=s'				=> \$tabix,
	'vcftools_dir=s'		=> \$vcftools_dir,
	'no_tabix!'				=> \$no_tabix,
	'help!'					=> \$help,
);

if ($help) {
	 exec('perldoc', $0);
}

if (!$region && !$no_tabix) {
	die("Please specify chromosome region. If the input VCF contains small number of sites, please use -no_tabix if you don't want to specify chromosome region\n");
}	 

if($no_tabix){
  unless(-e $vcf){
    die($vcf." must exist if tabix is not being run");
  }
}

	 	
$tabix = "/nfs/1000g-work/G1K/work/bin/tabix/tabix" if (!$tabix);
$vcftools_dir = "/nfs/1000g-work/G1K/work/bin/vcftools" if (!$vcftools_dir);

$vcftools = $vcftools_dir . "/bin/vcftools";
$vcf_subset = $vcftools_dir . "/perl/vcf-subset";
$vcf_fill_an_ac = $vcftools_dir . "/perl/fill-an-ac";

my $panel_hash = parse_sample_panel($sample_panel);
my @selected_pops;

if ( (defined $user_input_pop_string && $user_input_pop_string =~ /ALL/i) || !$user_input_pop_string) {
  push @selected_pops, keys %$panel_hash;
  $user_input_pop_string = "ALL";
	#foreach my $pop (keys %$panel_hash ) {
	#	process(join(",", @{$panel_hash->{$pop}}), $pop);
	#}	
	#print_hash(\%result_hash, "all");
}
else {
	my @user_input_pops = split(/,/, $user_input_pop_string);
	foreach my $pop (keys %$panel_hash ) {
		foreach my $user_pop (@user_input_pops) {
			$user_pop =~ s/^\s+|\s+$//;
			if ($pop =~ /$user_pop/i) {
        push @selected_pops, $user_pop;
				#process(join(",", @{$panel_hash->{$pop}}), $user_pop);
			}
		}	
	}	
	#print_hash(\%result_hash, $user_input_pop_string);  
}	

foreach my $p (@selected_pops){
  #print $p."\n";
}

#my %selected_pop_samples;
foreach my $p (@selected_pops){
  $selected_pop_samples{$p} = join(",", @{$panel_hash->{$p}});
}

process(\%selected_pop_samples);

####################
###### SUBS ########
####################
sub process {
	my (%pop_samples) = %{shift @_};
  my @pops;
  push @pops, keys %pop_samples;


  #print header for output file
  my $file_region = $region =~ s/:/\./r;
  my $outfile = "$out_dir/calculated_fra.process$$" . "." . $file_region . "." . $user_input_pop_string.".txt";
  $outfile =~ s/,/_/g;
  $outfile =~ s/\s+//;
  my $out_h;
  open ($out_h, ">", $outfile) || die("Cannot open output file $outfile");
  print $out_h "CHR\tPOS\tID\tREF\tALT\t";
  print $out_h "ALL_POP_TOTAL_CNT\tALL_POP_ALT_CNT\tALL_POP_FRQ\t";
  foreach my $pop ( @pops ) {
    print $out_h $pop . "_TOTAL_CNT\t";
    print $out_h $pop . "_ALT_CNT\t";
    print $out_h $pop . "_FRQ\t";
  }


	my $vcf_cmd;
	if ($no_tabix) {
		$vcf_cmd = "cat $vcf";
	}
	else {
		$vcf_cmd = "$tabix -f -h $vcf $region";# | $vcf_subset -c $samples_string | $vcf_fill_an_ac | cut -f1-8";
	}

	open(VCF, $vcf_cmd." | ") or die("Failed to open ".$vcf_cmd." $!");

  #Read VCF in batches
  #for each batch
    #process each pop
      #AN,AC,AF calcs
    #print batch output
  #set batch size in n
  my $n = 100000;
  my @header;
  my $head = 1;
  my $line;
  my @body;
  my $body_ln_cnt = 0;


  while(<VCF>){
    chomp;
    $line = $_;
    #print $line."\n";
    #read header
    if($head){
      if($line =~ /^#/){
        my $hline = $line."\n";
        push @header, $hline;
        #print $hline;
      }
      else{
        $head = 0;
      }
    }
    #process body
    if(!$head){
      $body_ln_cnt++;
      my $bline = $line."\n";
      push @body, $bline;

      if(($body_ln_cnt % $n) == 0){
        #print scalar(keys %result_hash)." keys\n";
        print_VCF_block_results($out_h, \@pops, \%result_hash);
        #empty structures before next block
        @body = ();
        %result_hash = ();
        #print scalar(@body)."\n";
      } 

    }
  }
  #process last block of lines from body
  process_VCF_block(\@header, \@body, \%pop_samples, \@pops);
  print_VCF_block_results($out_h, \@pops, \%result_hash);

  return 1;
}

sub process_VCF_block(){
  my @header = @{shift @_};
  my @body = @{shift @_};
  my %pop_samples = %{shift @_};
  my @pops = @{shift @_};

  foreach my $p (@pops){
    #print $p."\n";
    my $samples_string = $pop_samples{$p};
    #print $samples_string."\n";

    $region =~ s/:/\./;
    my $tmpfile = "$out_dir/tmp_calculated_fra.process$$" . "." . $region . "." . $p.".vcf";
    $tmpfile =~ s/,/_/g;
    $tmpfile =~ s/\s+//;
    die "$tmpfile already exists\n" if(-e $tmpfile);
    my $tmp_h;
    open ($tmp_h, ">", $tmpfile) || die("Cannot open output file $tmpfile");

    my @mini_vcf = (@header, @body);
    for(my $v=0;$v<scalar(@mini_vcf);$v++){
      print $tmp_h $mini_vcf[$v];
    }
    close $tmp_h;

    my $pop_vcf_cmd = "cat $tmpfile | $vcf_subset -c $samples_string | $vcf_fill_an_ac | cut -f1-8";
    open(POP_VCF, $pop_vcf_cmd." |") or die("Failed to open ".$pop_vcf_cmd." $!");
    #my $pop_line_cnt = 0;
    POPVCFLINE: while(<POP_VCF>){
      chomp;
      my $ret_line = $_;
      #print $ret_line."\n";
      next POPVCFLINE if($ret_line =~ /^#/);
      #$pop_line_cnt++;
      process_ac_an_line($ret_line, $p);
    }
    #print "pop lines returned from an/ac: $pop_line_cnt\n";
    close(POP_VCF);
    system("rm $tmpfile");
  }
}

sub process_ac_an_line{
  my $line = shift @_;
  my $population = shift @_;
  #print $line."\n";

  
  my @values = split /\t/, $line;
  my @info = split /\;/, $values[7];

  my $ac;
  my $an;
  my $af = 0;
  my $alt_alleles = $values[4];

    foreach my $info(@info){
      if($info =~ /^AC/){
        my ($tag, $num) = split /\=/, $info;
        $ac = $num;
      }
      if($info =~ /^AN/){
        my ($tag, $num) = split /\=/, $info;
        $an = $num;
      }
    }
    
    if(!defined($ac)){
      die("Failed to find ac from ".join("\t", $values[7]));
    }
    if(!$an){
      die("Failed to find an from ".join("\t", $values[7]));
    }

    #multiple alleles
    if($ac =~ m/,/){
      #$cmd can return multiple alleles in varying order T,A or A,T, for example
      #need to account for this for combining populations in output
      #reporting in alphabetical order for consistency across pops
      my @alts = split /,/, $values[4];
      my @sorted_alts = sort @alts;
      $alt_alleles = join ",", @sorted_alts;
      my @acs = split /,/, $ac;
      my @sorted_acs;
      # puts AC values into sorted order
      ALLELE: for (my $i = 0; $i < scalar(@sorted_alts); $i++){
        for(my $j =0; $j < scalar(@alts); $j++){
          if($sorted_alts[$i] eq $alts[$j]){
            $sorted_acs[$i] = $acs[$j];
            next ALLELE;
          } 
        }
      }
      #calc frequencies and convert to string
      my @afs;
      foreach my $sort_ac (@sorted_acs){
        my $allele_freq = sprintf "%.2f", $sort_ac/$an;
        push @afs, $allele_freq;
      }
      $af = join(",",@afs);
      $ac = join(",",@sorted_acs);
    }
    else{
      #not multi-allelic, so just do calc
      $af = sprintf "%.2f", $ac/$an if($ac);
    }

    my $site_info = join ("\t", $values[0], $values[1],$values[2],$values[3], $alt_alleles);
    my $site_stats = join ("\t", $an, $ac, $af);
    #print $population." ".$site_info." ".$site_stats." being added to result hash\n";
    $result_hash{$site_info}{$population} = $site_stats; 
}


sub parse_sample_panel {
	my ($panel_file) = @_;
	my %pop_sample_hash;
	
	my $fh;
	open($fh, "<", $panel_file);
	while (<$fh>) {
		chomp;
		my $sample_line = $_;
		my ($sample, $pop_in_panel, $sup_pop, $platform) = split(/\t/, $sample_line);
    #add unless it looks like the header
		push @{$pop_sample_hash{$pop_in_panel}}, $sample unless $sample eq "sample";
	}
	return \%pop_sample_hash;
}	

sub print_VCF_block_results{
  my $out_h = shift @_;
  my @pops = @{shift @_};
  my $hash = shift @_;
	my %total_cnt_all_pops;
	my %alt_cnt_all_pops;
	
  my @random_sites;	
	my %position_hash;
  #print scalar(keys %{$hash})." keys in hash in print\n";
	foreach my $site (keys %{$hash}) {
    #print $site."\n";
		$total_cnt_all_pops{$site} = 0;
		$alt_cnt_all_pops{$site} = 0;
		foreach my $p (@pops) {
      #print "popoulation p: $p\n";
      #print $hash->{$site}->{$p}."\n";
			my @data = split(/\t/, $hash->{$site}->{$p});
			$total_cnt_all_pops{$site} += $data[0];
			if($data[1] =~ m/,/){#multiple alt alleles
			  if($alt_cnt_all_pops{$site} eq 0){
			    $alt_cnt_all_pops{$site} = $data[1];
			  }
			  else{
				my @curr_totals = split /,/, $alt_cnt_all_pops{$site};
				my @new_cnts = split /,/, $data[1];
				for(my $i=0; $i<scalar(@curr_totals);$i++){
				  $curr_totals[$i] = $curr_totals[$i] + $new_cnts[$i];
				}
				my $new_all_alt_cnt = join(",",@curr_totals);
				$alt_cnt_all_pops{$site} = $new_all_alt_cnt;
			  }
			}
			else{#not multiple alt alleles
			  $alt_cnt_all_pops{$site} += $data[1];
			}
		}
		my @meta_array = split(/\t/, $site);
    #print $site." site\n";
		#my $pos = $meta_array[0] . "." . $meta_array[1];
    #print $pos."\n";
    push @random_sites, {
      'chr' => $meta_array[0],
      'loc' => $meta_array[1],
      'rs' => $meta_array[2],
      'ref' => $meta_array[3],
      'alt' => $meta_array[4]
    };
		#$position_hash{$pos} = $site;
	}	
  #print scalar(keys %position_hash)." keys in position hash\n";	
	my @sorted_sites = sort {
    $a->{'chr'} <=> $b->{'chr'} ||
    $a->{'loc'} <=> $b->{'loc'} ||
    $a->{'rs'} cmp $b->{'rs'} ||
    $a->{'ref'} cmp $b->{'ref'} ||
    $a->{'alt'} cmp $b->{'alt'}
  } @random_sites;
	
	#print scalar(@sorted_sites)." sorted positions in results\n";			
		foreach my $s (@sorted_sites) {
      my $site = $s->{'chr'}."\t".$s->{'loc'}."\t".$s->{'rs'}."\t".$s->{'ref'}."\t".$s->{'alt'};
			print $out_h "$site\t";
			print $out_h $total_cnt_all_pops{$site} . "\t" . $alt_cnt_all_pops{$site} . "\t";
			my $alt_afs_all_pops;
			if($alt_cnt_all_pops{$site} =~ m/,/){
			  #handle multi-allelic site
			  my @alt_cnts = split /,/, $alt_cnt_all_pops{$site};
			  my @alt_afs;
			  foreach my $alt_cnt (@alt_cnts){
				my $alt_af = sprintf "%.2f", $alt_cnt/$total_cnt_all_pops{$site};
				push @alt_afs, $alt_af;
			  }
			  $alt_afs_all_pops = join(',',@alt_afs);
			}
			else{
				$alt_afs_all_pops = sprintf "%.2f", $alt_cnt_all_pops{$site}/$total_cnt_all_pops{$site};
			}
			print $out_h $alt_afs_all_pops;
			print $out_h "\t";
			my @numbers;
			foreach my $pop_key (@pops ) {
				push @numbers, $hash->{$site}->{$pop_key};
			}
			print $out_h join ("\t", @numbers) . "\n";
		
		}
}		


=pod

=head1 NAME

calculate_allele_frq_from_vcf_file.pl

=head1 SYNOPSIS

This script takes a VCF file, a matching sample panel file, a chromosomal region, population names, it then calculates population-wide allele 
frequency for sites within the chromosomal region defined.

When no population is specified, allele fequences will be calcuated for all populations in the VCF files, one at a time.

=head1 Dependency

	To run this script, you need to install the following software
	tabix: http://sourceforge.net/projects/samtools/files/tabix/
	vcftools: http://sourceforge.net/projects/vcftools/files/
	
=head1 Options

	-vcf			Input VCF file that contains genotype data for each samples; this file must be bgzipped and tabix indexed; if '-no_tabix' is used, the 
vcf file can be uncompressed and un-indexed.
	-sample_panel	A tab deliminated file lists mapping between sample and population (see example below) 
	-pop			Populations of interest; separated by ",".  This field can be null
	-region			chromosome region of interest, format is chr_number:start_end (optional). 
	-out_dir		Where temporary file and final output files should be written to. Default is current direcotry.
	-tabix			A path to tabix executable
	-vcftools_dir	A path to the vcftools base directory that contains vcftools executable and perl libraries
	-no_tabix		If the input VCF files is a pre-sliced VCF file containing a small number of sites, this option can be used so the vcf file doesn't have 
					to be tabix indexed or bgzipped.  This is to speed up run time for the web application.
	-help			Print this page when specified

=head1 EXAMPLE lines from a sample panel file. Only the first 2 columns are essential.

	HG00096	GBR	EUR	ILLUMINA
	HG00097	GBR	EUR	ABI_SOLID
	HG00099	GBR	EUR	ABI_SOLID
	HG00100	GBR	EUR	ILLUMINA

=head1 OUTPUT

The allele frequency of an user-specified population for sites within the user-specified chromosomal region is written to a file.  The headers of the 
output file are:

	CHR:		Chromosome
	POS:		Start position of the variant
	ID:		Identification of the variant
	REF:		Reference allele
	ALT:		Alternative allele
	TOTAL_CNT:	Total number of alleles in samples of the chosen population(s) 
	ALT_CNT:	Number of alternative alleles observed in samples of the chosen populations(s)
	FRQ:		Ratio of ALT_CNT to TOTAL_CNT


=head1 EXAMPLE

perl $ZHENG_RP/bin/calculate_allele_frq_from_vcf.pl \
-vcf /nfs/1000g-archive/vol1/ftp/phase1/analysis_results/integrated_call_sets/ALL.chr22.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz \
-out_dir . \
-sample_panel /nfs/1000g-archive/vol1/ftp/phase1/analysis_results/integrated_call_sets/integrated_call_samples.20101123.ALL.panel \
-region 22:17000000-17005000 \
-pop CEU,FIN \

OR (when the input VCF file is pre-sliced small size file)

perl $ZHENG_RP/bin/calculate_allele_frq_from_vcf.pl \
-vcf ALL.chr22_17000000_17005000.test.vcf \
-out_dir ~ \
-sample_panel /nfs/1000g-archive/vol1/ftp/phase1/analysis_results/integrated_call_sets/integrated_call_samples.20101123.ALL.panel \
-pop CEU,FIN \
-no_tabix 

