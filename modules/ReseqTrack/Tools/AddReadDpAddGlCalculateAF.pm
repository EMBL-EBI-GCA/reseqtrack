
=pod

=head1 NAME

ReseqTrack::Tools::AddReadDpAddGlCalculateAF

=head1 SYNOPSIS

This is the core module to annotate a light weighted VCF genotype files with the following information:
1. Per sample per site read depth and total depth (from a pre-calculated depth matrix)
2. Allele frequency by continental super population (and global AF if not already in the light weighted VCF file) per site. AF for multi-allelic variants are 
calculated and reported for each allele independently
3. Optional.  Genotype likehood for each genotype (from two GL files, one for simple variants, one for complex variants including SV and STR etc.; GL from the simple 
variant GL file takes precedent if a site show up in both simple and complex variant GL list).  Do this when the 'gl_flag' is set

when 'gl_flag' is set, GL will be added to the GT columns; when 'per_sample_dp_flag' is set, per sample dp will be added to individual GT columns. However, the total dp will be added to the INFO field 
regardless of this flag.
 
example

	my $object = AddReadDpAddGlCalculateAF->new (
		-program					=> /path/to/tabix/executable,
		-input_files				=> /path/to/input/light_weighted/VCF,
		-depth_matrix				=> /path/to/depth/matrix,
		-supp_dp_mx					=> /path/to/supplement/depth/matrix,
		-region						=> /chrom/region,
		-working_dir				=> /path/to/out_put directory,
		-sample_panel				=> /path/to/sample/panel/file,
		-related_samples			=> /path/to/list/of/samples/that/are/related
		-simple_var_gl_file			=> /path/to/simple_var_gl_file,
		-complex_var_gl_file		=> /path/to/complex_var_gl_file,
		-gl_flag					=> 0,
		-per_sample_dp_flag			=> 0,
		-save_files_from_deletion	=> 0,
	);

=cut

package ReseqTrack::Tools::AddReadDpAddGlCalculateAF;
use strict;
use warnings;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename qw(basename);
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists check_executable get_lines_from_file);
use ReseqTrack::Tools::RunProgram;

use vars qw(@ISA);
@ISA = qw(ReseqTrack::Tools::RunProgram);

sub new {
	my ( $class, @args ) = @_;
	my $self = $class->SUPER::new(@args);

	my (	$depth_matrix, 
			$supp_dp_mx, 
			$region, 
			$sample_panel, 
			$related_sample_list, 
			$simple_var_gl_file, 
			$complex_var_gl_file,
			$gl_flag,
			$per_sample_dp_flag ) = rearrange(
		[
			qw(
			  DEPTH_MATRIX
			  SUPP_DP_MX
			  REGION
			  SAMPLE_PANEL
			  RELATED_SAMPLE_LIST
			  SIMPLE_VAR_GL_FILE
			  COMPLEX_VAR_GL_FILE
			  GL_FLAG
			  PER_SAMPLE_DP_FLAG
			  )
		],
		@args
	);
	
	$self->depth_matrix($depth_matrix);
	$self->supp_dp_mx($supp_dp_mx);
	$self->region($region);
	$self->sample_panel($sample_panel);
	$self->related_sample_list($related_sample_list);	
	$self->simple_var_gl_file($simple_var_gl_file) if ($simple_var_gl_file && $gl_flag);	
	$self->complex_var_gl_file($complex_var_gl_file) if ($complex_var_gl_file  && $gl_flag);	
	$self->gl_flag($gl_flag) if ($gl_flag);
	$self->per_sample_dp_flag($per_sample_dp_flag) if ($per_sample_dp_flag);
	
	return $self;
}

sub run_program {
	my $self = shift;
	
	print "input vcf file is " ;
	print join (" ", @{$self->input_files}) . "\n";
	print "depth matrix file is " . $self->depth_matrix . "\n";
	
	#print "gl_flag is " . $self->gl_flag . "\n";
	
	foreach (@{$self->input_files}) {
		check_file_exists($_);
		check_file_exists($_ . ".tbi");
	}	 
	
	check_file_exists($self->depth_matrix);
	check_file_exists($self->depth_matrix . ".tbi");

	check_file_exists($self->supp_dp_mx);
	check_file_exists($self->supp_dp_mx . ".tbi");
			
	check_file_exists($self->sample_panel);
	
	check_executable($self->program); ### tabix
	
	$self->find_sample_spop;
	$self->get_related_samples;

	if ( $self->gl_flag ) {
		throw("please provide GL files when -gl_flag is set") unless ($self->simple_var_gl_file && $self->complex_var_gl_file);
		check_file_exists($self->simple_var_gl_file);
		check_file_exists($self->simple_var_gl_file . ".tbi");	
		check_file_exists($self->complex_var_gl_file);
		check_file_exists($self->complex_var_gl_file . ".tbi");
		$self->get_complex_var_gl_for_a_region;
	}
	
	$self->get_supp_depth_for_a_region;  ## this is for the supplement depth matrix
	
	$self->run_annotation;
	
	print join ("output file are " , " ", @{$self->output_files})  . "\n";
	return $self;
}

sub run_annotation {
	my ($self) = @_;

	my $out_file = $self->working_dir . "/" . $self->region . ".vcf";
	$out_file =~ s/:/_/g;

	open (my $ofh, ">", $out_file) || throw("Cannot open output file $out_file");

	my $open_vcf = $self->get_tabix_file_handler($self->input_files->[0]);	
	my $open_dp_matrix = $self->get_tabix_file_handler($self->depth_matrix);
	my $open_simple_gl_file = $self->get_tabix_file_handler($self->simple_var_gl_file) if ($self->gl_flag);
				
	my @vcf_header;	
	my @mx_header;
	my @gl_header;
	
	my $num_of_samples_in_mx;
	
	my $matrix_cursor_wait = 0;
	my $gl_cursor_wait = 0;
	
	my $gl_line;
	my $mx_line;
	
	my $seen_vcf_site = "";
	
	my ($user_chrom, $user_region) = split(/:/, $self->region);
	my ($user_start, $user_end) = split(/-/, $user_region);
			
	while (my $vcf_line = <$open_vcf>) {
		chomp $vcf_line;
		
		if ($vcf_line =~ /^\#\#/) {
			if ($vcf_line =~ /^\#\#bcftools_annotate/) {
				next;
			}
			print $ofh $vcf_line . "\n";
			next;
		}	
		elsif ($vcf_line =~ /^\#CHROM/) {
			
			print $ofh "\#\#INFO=<ID=SAS_AF,Number=A,Type=Float,Description=\"Allele frequency in the SAS populations calculated from AC and AN, in the range (0,1)\">\n";
			print $ofh "\#\#INFO=<ID=EUR_AF,Number=A,Type=Float,Description=\"Allele frequency in the EUR populations calculated from AC and AN, in the range (0,1)\">\n";
			print $ofh "\#\#INFO=<ID=AFR_AF,Number=A,Type=Float,Description=\"Allele frequency in the AFR populations calculated from AC and AN, in the range (0,1)\">\n";
			print $ofh "\#\#INFO=<ID=AMR_AF,Number=A,Type=Float,Description=\"Allele frequency in the AMR populations calculated from AC and AN, in the range (0,1)\">\n";
			print $ofh "\#\#INFO=<ID=EAS_AF,Number=A,Type=Float,Description=\"Allele frequency in the EAS populations calculated from AC and AN, in the range (0,1)\">\n";
			print $ofh "\#\#INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth\">\n";
			print $ofh "\#\#FORMAT=<ID=DP,Number=1,Type=String,Description=\"Per sample read depth\">\n" if ($self->per_sample_dp_flag);
			print $ofh "\#\#FORMAT=<ID=GL,Number=1,Type=String,Description=\"Genotype likelihood\">\n" if ($self->gl_flag);
			print $ofh "\#\#FORMAT=<ID=CNL,Number=1,Type=String,Description=\"Genotype likelihood for CNVs\">\n" if ($self->gl_flag);
			print $ofh "\#\#FORMAT=<ID=PL,Number=1,Type=String,Description=\"Genotype likelihood for STRs\">\n" if ($self->gl_flag);
			
			print $ofh $vcf_line . "\n";
			@vcf_header = split(/\t/, $vcf_line);
			next;
		}
		
		my @vcf_data = split(/\t/, $vcf_line);		
		
		my $vcf_chrom = $vcf_data[0];
		my $vcf_start = $vcf_data[1];
		my $vcf_allele = $vcf_data[4];
		print STDERR ("site is $vcf_chrom _ $vcf_start\n");	
		
		if ($vcf_chrom ne $user_chrom ) {
			throw("Wrong chromosome chosen for " . $self->region . "\n");
		}	
		if ($vcf_start < $user_start || $vcf_start > $user_end) {  	### This is to handle tabix SV sites show up in multiple chunks
																	#### If not handled here, the second chunk would have no depth in the supplement matrix. it will get a "." as depth.
																	#### this is a source of duplicated lines.  merge VCF may not be able to merge two identical sites with different information content
																	#### The site in the second chunnk will have gl in the complex gl file as tabix works the same way there. 
																	#### This will introduce duplicated lines		
			print STDERR "site $vcf_chrom _ $vcf_start outside user specified region " . $self->region . "\n";
			next;	
		}
		
		if ($seen_vcf_site && $seen_vcf_site == $vcf_start) {
			print STDERR "Duplicated site in VCF file $vcf_chrom $vcf_start, don't move the depth matrix and GL cursors\n";
			$matrix_cursor_wait=1;
			$gl_cursor_wait = 1;
		}
		else {
			$seen_vcf_site = $vcf_start;
		}
		
		MATRIX:
		if ($matrix_cursor_wait==1) {
			$matrix_cursor_wait = 0;
		}
		else {
			$mx_line = <$open_dp_matrix>;
			chomp $mx_line;
		}

		if ($mx_line =~ /^\#CHROM/) {
			@mx_header = split(/\t/, $mx_line);
			shift @mx_header;
			shift @mx_header;
			$num_of_samples_in_mx = scalar @mx_header;
			$mx_line = <$open_dp_matrix>;  ## this serves as 'next' into the matrix file
			chomp $mx_line;
		}
		
		my @mx_depth_all_samples = split(/\t/, $mx_line);
		my $mx_chrom = shift @mx_depth_all_samples;
		my $mx_start = shift @mx_depth_all_samples;	
		
		if ($vcf_chrom ne $mx_chrom) {
			print STDERR "matrix file on chrom $mx_chrom, input vcf file on $vcf_chrom\n";
		}
		
		
		my %depth_hash_by_line;
		if ($vcf_chrom eq $mx_chrom && $vcf_start == $mx_start) {	
			for (my $i=0; $i < $num_of_samples_in_mx; $i++) {
				$depth_hash_by_line{$mx_header[$i]} = $mx_depth_all_samples[$i];
			}
		}
		elsif ( (	$vcf_chrom eq $mx_chrom && $vcf_start < $mx_start ) || 
					!$mx_line ) {		##if the mx file reach the end in the tatix block and doesn't have a matching site		
			print STDERR "New site $vcf_chrom $vcf_start, get dp from a separate matrix\n";
			for (my $i=9; $i <= $#vcf_header; $i++) {
				my $dp = $self->supp_dp_hash->{$vcf_chrom . "_" . $vcf_start}->{$vcf_header[$i]};
				if ( $dp eq "" ) {	### have to use this as depth "0" would be considered as null, cannot use "defined $dp"
					print STDERR ($vcf_data[0] . "_" . $vcf_data[1] . "_" . $vcf_header[$i] . " do not have depth\n");  
					$dp = ".";   
				}
				my $mx_i = $i - 9;
				$depth_hash_by_line{$mx_header[$mx_i]} = $dp;
				#print STDERR "Found supp dp: $mx_header[$mx_i] $dp\n";			  
			}	
			$matrix_cursor_wait = 1;
		}	
		elsif ($vcf_chrom eq $mx_chrom && $vcf_start > $mx_start)  {
			print STDERR "Skip a site $mx_chrom $mx_start in the depth matrix\n";
			goto MATRIX;
		}	
		else {
			throw("What's wrong with the vcf file $vcf_chrom $vcf_start and matrix file $mx_chrom $mx_start");
		}	
		
		my %gl_hash_by_line;
		if ($self->gl_flag) {
			GL:
			if ($gl_cursor_wait==1) {
				$gl_cursor_wait = 0;
			}
			else {
				$gl_line = <$open_simple_gl_file>;  ### at the end of document, this returns an empty line
				chomp $gl_line;
			}
		
			if ($gl_line =~ /^##/) {
				HEADER:
				$gl_line = <$open_simple_gl_file>;
				if ($gl_line =~ /^##/) {
					goto HEADER;
				}	
				chomp $gl_line;	
			}	
					
			if ($gl_line =~ /^\#CHROM/) {
				@gl_header = split(/\t/, $gl_line);
				$gl_line = <$open_simple_gl_file>;
				chomp $gl_line;	
			}
			
			my @gl_data = split(/\t/, $gl_line);
			my $gl_chrom = $gl_data[0];
			my $gl_start = $gl_data[1];
			
			print STDERR ("simple var GL file on chrom $gl_chrom, input vcf file on $vcf_chrom") if ($vcf_chrom ne $gl_chrom);		
			
			if ( $vcf_chrom eq $gl_chrom && $vcf_start == $gl_start ) {
				$self->get_gl_hash_for_a_line($gl_line, \@gl_header);
				%gl_hash_by_line = %{$self->gl_hash};
			}
			elsif ( (	$vcf_chrom eq $gl_chrom && $vcf_start < $gl_start) || 
						!$gl_line ) { ## if the gl file reach the end in the batix block and doesn't have a matching site
				print STDERR "site $vcf_chrom $vcf_start, see if the site is in the complex var gl file\n";
				for (my $i=9; $i <= $#gl_header; $i++) {				
					if ($self->complex_var_gl_hash->{$vcf_chrom . "_" . $vcf_start}->{$vcf_allele}->{$vcf_header[$i]})  {  ### adding the allele is to make sure the allele is identical in the complex VAR gl file and the input vcf file, like chr11:90692
						my $complex_var_gl = $self->complex_var_gl_hash->{$vcf_chrom . "_" . $vcf_start}->{$vcf_allele}->{$vcf_header[$i]};
						$gl_hash_by_line{$gl_header[$i]} = $complex_var_gl;
						#print STDERR ("Found complex var GL: $gl_header[$i] $complex_var_gl\n");
					}
					else {
						print STDERR ("No simple or complex var GL for site $vcf_chrom  _  $vcf_start\n");
						$gl_hash_by_line{$gl_header[$i]} = "."; 
					}												  
				}	
				$gl_cursor_wait = 1;
			}	
			elsif ( $vcf_chrom eq $gl_chrom && $vcf_start > $gl_start ) {
				print STDERR "Skip a site $gl_chrom $gl_start in the simple var gl file\n";
				goto GL;
			}	
			else {
				throw("What's wrong with the vcf file $vcf_chrom $vcf_start and simple var gl file $gl_chrom $gl_start");
			}	
	
			if ($vcf_data[7] =~ /SVTYPE=CNV/) {
				$vcf_data[8] .= ":DP:CNL";
			}
			elsif ($vcf_data[7] =~ /VT=STR/) {
				$vcf_data[8] .= ":DP:PL";
			}
			else {	
				$vcf_data[8] .= ":DP:GL";  
			}
	
		}
		elsif ($self->per_sample_dp_flag) {
			if ($vcf_data[7] =~ /SVTYPE=CNV/) {
				$vcf_data[8] .= ":DP";
			}
			elsif ($vcf_data[7] =~ /VT=STR/) {
				$vcf_data[8] .= ":DP";
			}
			else {	
				$vcf_data[8] .= ":DP";  
			}
		}
		else {	
			
		}	
			
		my @multi_alt_alleles = split(/,/, $vcf_data[4]);  ### this is to help calculating AF for multi-allelic sites 
		my $multi_allelic_levels = $#multi_alt_alleles + 1;
		
		my $sum_d = 0;
		$self->init_allele_cnts_per_site($multi_allelic_levels);
		
		for (my $j = 9; $j <= $#vcf_header; $j++) {	
			
			### DEPTH ###
			my $d = $depth_hash_by_line{$vcf_header[$j]};
			#print "vcf pos: " . $vcf_data[0] . "_" . $vcf_data[1] . " " . $vcf_header[$j] . " depth " . $d . "\n";
			if ( $d eq "" ) {	
				print STDERR ($vcf_data[0] . "_" . $vcf_data[1] . "_" . $vcf_header[$j] . " does not have depth\n");  
				$d = "."; 
			}
			$sum_d += $d if ($d ne ".");
			$vcf_data[$j] .= ":" . $d if ($self->per_sample_dp_flag); 
		
			### GL ###
			if ($self->gl_flag)	{
				my $gl = $gl_hash_by_line{$vcf_header[$j]};
				if (!$gl) {
					print STDERR ($vcf_data[0] . "_" . $vcf_data[1] . "_" . $vcf_header[$j] . " does not have GL\n");  
					$gl = ".";
				}			
				$vcf_data[$j] .= ":" . $gl ;
			}
			
			#### AF ###
			if ( $self->related_sample_hash->{$vcf_header[$j]} ) { 			
				#print STDERR ("skipping sample $vcf_header[$j] for allele frq calculation as it is one of the related samples\n");
			}
			else {	
				$self->populate_allele_cnt_per_site_hash($vcf_header[$j], $multi_allelic_levels, $vcf_data[$j]);		
			}
		}
		if ($sum_d == 0 && $self->region !~ /Y/i ) {
			$vcf_data[7] .= ";DP=" . ".";	
		}
		elsif ($self->region !~ /Y/i) {
			$vcf_data[7] .= ";DP=" . $sum_d;
		}	 	
		$vcf_data[7] .= ";" . $self->calculate_AF_per_site->AF;  
		print $ofh join ("\t", @vcf_data) . "\n";
	}
	
	close($ofh);
	$self->output_files($out_file);
	return $self;
}

sub get_supp_depth_for_a_region {
	my ($self) = @_;
	my %depth_hash = (); 
	
	my $open_supp_dp_matrix = $self->get_tabix_file_handler($self->supp_dp_mx);
	
	my @h;
	my $number_of_samples;
	while (<$open_supp_dp_matrix>) {
		chomp;
		#print "get depth line is $_\n";
		if ($_ =~ /^\#CHROM/) {
			@h = split(/\t/, $_);
			shift @h;
			shift @h;
			$number_of_samples = scalar @h;
			next;
		}
		
		my @depth_all_samples = split(/\t/, $_);
		my $chrom = shift @depth_all_samples;
		my $start = shift @depth_all_samples;
		
		for (my $i=0; $i < $number_of_samples; $i++) {
			#print "matrix pos is $chrom _ $start, sample is $h[$i], dpeth is $depth_all_samples[$i]\n";
			$depth_hash{$chrom . "_" . $start}{$h[$i]} = $depth_all_samples[$i];
		}
	}
	$self->supp_dp_hash(\%depth_hash);
	return $self;
}		

sub get_gl_hash_for_a_line {
	my ($self, $line, $header) = @_;
	
	my $number_of_samples = scalar @$header - 9;
	my %gl_hash;
			
	my @data = split(/\t/, $line);
	my $chrom = $data[0];
	my $start = $data[1];
		
	my @format_headers = split(/:/, $data[8]);
	my $gl_index = "";	
	my $found = 0;
	if ($data[7] =~ /SVTYPE=CNV/) {
		for (my $k=0; $k<= $#format_headers; $k++) {
			if ($format_headers[$k] eq "CNL") {  
				$gl_index=$k;
				#print "GL index field is $gl_index\n";
				$found = 1;
			}
		}
	}
	elsif ($data[7] =~ /VT=STR/) {
		for (my $k=0; $k<= $#format_headers; $k++) {
			if ($format_headers[$k] eq "PL") {  
				$gl_index=$k;
				#print "GL index field is $gl_index\n";
				$found = 1;
			}
		}
	}
	else {
		for (my $k=0; $k<= $#format_headers; $k++) {
			if ($format_headers[$k] eq "GL" || $format_headers[$k] eq "GL0") {  
				$gl_index=$k;
				#print "GL index field is $gl_index\n";
				$found = 1;
			}
		}
	}		
	throw("Cannot find GL field for site $chrom $start") if ($found == 0);
				
	for (my $n=9; $n < ($number_of_samples+9); $n++) {
		#print "matrix pos is $chrom _ $start, sample is $h[$n], dpeth is $data[$n]\n";
		my @gt_data_per_sam = split(/:/, $data[$n]);
		$gl_hash{$header->[$n]} = $gt_data_per_sam[$gl_index];
	}
	$self->gl_hash(\%gl_hash);
	return $self;
}		
	

#### Only use this hash module for complex var gls
sub get_complex_var_gl_for_a_region {
	my ($self) = @_;  

	my %complex_var_gl_hash = ();
	
	my $open_gl = $self->get_tabix_file_handler($self->complex_var_gl_file);
	
	my @h;
	my $number_of_samples;
	while (<$open_gl>) {
		chomp;
		if ($_ =~ /^\#CHROM/) {
			@h = split(/\t/, $_);;
			$number_of_samples = scalar @h - 9;
			next;
		}
		elsif ($_ =~ /^##/) {
			next;
		}	
		
		my @data = split(/\t/, $_);
		my $chrom = $data[0];
		my $start = $data[1];
		my $allele = $data[4];
		
		my @format_headers = split(/:/, $data[8]);
		my $gl_index = "";	
		my $found = 0;
		if ($data[7] =~ /SVTYPE=CNV/) {
			for (my $k=0; $k<= $#format_headers; $k++) {
				if ($format_headers[$k] eq "CNL") {  
					$gl_index=$k;
					#print "GL index field is $gl_index\n";
					$found = 1;
				}
			}
		}
		elsif ($data[7] =~ /VT=STR/) {
			for (my $k=0; $k<= $#format_headers; $k++) {
				if ($format_headers[$k] eq "PL") {  
					$gl_index=$k;
					#print "GL index field is $gl_index\n";
					$found = 1;
				}
			}
		}
		else {
			for (my $k=0; $k<= $#format_headers; $k++) {
				if ($format_headers[$k] eq 'GL' || $format_headers[$k] eq 'GL0' ) {  
					$gl_index=$k;
					#print "GL index field is $gl_index\n";
					$found = 1;
				}
			}
		}	
		
		throw("Cannot find GL field for site $chrom $start") if ($found == 0);
				
		for (my $n=9; $n < ($number_of_samples+9); $n++) {
			my @gt_data_per_sam = split(/:/, $data[$n]);
			$complex_var_gl_hash{$chrom . "_" . $start}{$allele}{$h[$n]} = $gt_data_per_sam[$gl_index];
			#print "From GL file: $chrom" .  "_" .  "$start $h[$n]\t$gt_data_per_sam[$gl_index]\n";		
		}	
	}

	$self->complex_var_gl_hash(\%complex_var_gl_hash);

	return $self;
}		

sub get_tabix_file_handler {
        my ($self, $input) = @_;
        my $cmd = $self->program . " -f -h ";
        $cmd .= $input . " ";
        $cmd .= $self->region . " ";
        print "command is $cmd\n";

        open(my $open_file_handle, $cmd." | ") or throw("Failed to open ".$cmd." $!");
        return $open_file_handle;
}

sub get_related_samples {
	my ($self) = @_;
	my %related_sample_hash = ();
	my $lines = get_lines_from_file($self->related_sample_list);
	foreach my $line (@$lines ) {
		chomp $line;
		next if ($line =~ /^Sample/);
		my @data = split(/\t/, $line);
		$related_sample_hash{$data[0]} = 1;
	}	
	$self->related_sample_hash(\%related_sample_hash);
	return $self;
}	
	

sub find_sample_spop {
	my ($self) = @_;

	my $sample_lines = get_lines_from_file($self->sample_panel);
		
	my %sample_to_spop;
	foreach my $sample_l ( @$sample_lines ) {
		next if ($sample_l =~ /Sample/i);
		my @data2 = split(/\t/, $sample_l);
		$sample_to_spop{$data2[0]} = $data2[2];
	}	
	$self->sample_to_spop(\%sample_to_spop);
	return $self;
}	
	
sub init_allele_cnts_per_site {
	my ($self, $allelic_levels) = @_;
	my %allele_cnt_pop;
	my @super_pops = qw(ALL AFR AMR SAS EUR EAS);	
	foreach my $p ( @super_pops ) {
		$allele_cnt_pop{$p}{'ref'}{0} = 0;
		for (my $l = 1; $l <= $allelic_levels; $l++) {
			$allele_cnt_pop{$p}{'alt'}{$l} = 0;
		}	
	}		
	$self->allele_cnt_per_site_hash(\%allele_cnt_pop);
	return $self;
}	


sub populate_allele_cnt_per_site_hash { 
	my ($self, $sample, $allelic_level3, $gt) = @_;
	
	my $hash = $self->allele_cnt_per_site_hash;
	my $spop = $self->sample_to_spop->{$sample};
	
	my @tmp = split(/:/, $gt);
	my ($a1, $a2) = split (/\|/, $tmp[0]);	
	
	if ($a1 == 0) {
		$hash->{'ALL'}->{'ref'}->{0}++;
		$hash->{$spop}->{'ref'}->{0}++;
	}
	
	if ($a2 == 0) {
		$hash->{'ALL'}->{'ref'}->{0}++;
		$hash->{$spop}->{'ref'}->{0}++;
	}
			
	for (my $m = 1; $m <= $allelic_level3; $m++) {
		if ($a1 == $m) {
			$hash->{'ALL'}->{'alt'}->{$m}++;
			$hash->{$spop}->{'alt'}->{$m}++;	
		}
		if ($a2 == $m) {
			$hash->{'ALL'}->{'alt'}->{$m}++;
			$hash->{$spop}->{'alt'}->{$m}++;	
		}
	}		
	
	$self->allele_cnt_per_site_hash($hash);
	return $self;	
}			

sub calculate_AF_per_site {
	my ($self) = @_;

	my $af_hash = $self->allele_cnt_per_site_hash;

	my @af_fields = ();
	
	foreach my $spop (keys %$af_hash) {	
		
		my $ref_cnt = $af_hash->{$spop}->{'ref'}->{0};
		
		my %alt_cnt_by_allelic_level;
		my $all_level_alt_cnt = 0;
		my $all_allele_cnt = 0;
		my @af_array = ();
		foreach my $allelic_level (keys %{$af_hash->{$spop}->{'alt'}}) {
			$alt_cnt_by_allelic_level{$allelic_level} = $af_hash->{$spop}->{'alt'}->{$allelic_level};	
			$all_level_alt_cnt += $alt_cnt_by_allelic_level{$allelic_level};
			#print STDERR "$spop Alt_allele_cnt " . $af_hash->{$spop}->{'alt'}->{$allelic_level} . " level $allelic_level " . "ref_allele cnt $ref_cnt\n";
		}
		
		$all_allele_cnt = $ref_cnt + $all_level_alt_cnt;		
		
		foreach my $allelic_level2 (sort {$a<=>$b} keys %alt_cnt_by_allelic_level) {
			my $af_by_allelic_level;
			if ($all_allele_cnt == 0 ) {
				$af_by_allelic_level = 0;
			}
			else {	
				$af_by_allelic_level = sprintf("%.4f", $alt_cnt_by_allelic_level{$allelic_level2}/$all_allele_cnt);
			}
			
			#print STDERR "ordered allelic level is $allelic_level2\n";
			push @af_array, $af_by_allelic_level;
		}	
	
		my $af_string_by_spop = join(',', @af_array);
		if ($spop ne 'ALL') {  ### FIXME: comment this out if the GT file doesn't have global AF
			push @af_fields, "$spop" . "_AF=" . $af_string_by_spop;
		}		
	}	

	my $af_string = join(";", @af_fields);	
	$self->AF($af_string);
	return $self;
}	

###################
#### Accessors ####
###################

sub depth_matrix {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{depth_matrix} = $arg;
  }
  return $self->{depth_matrix};
}

sub supp_dp_mx {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{supp_dp_mx} = $arg;
  }
  return $self->{supp_dp_mx};
}

sub region {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{region} = $arg;
  }
  return $self->{region};
}

sub related_sample_list {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{related_sample_list} = $arg;
  }
  return $self->{related_sample_list};
}

sub related_sample_hash {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{related_sample_hash} = $arg;
  }
  return $self->{related_sample_hash};
}

sub simple_var_gl_file {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{simple_var_gl_file} = $arg;
  }
  return $self->{simple_var_gl_file};
}

sub complex_var_gl_file {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{complex_var_gl_file} = $arg;
  }
  return $self->{complex_var_gl_file};
}

sub sample_panel {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{sample_panel} = $arg;
  }
  return $self->{sample_panel};
}

sub allele_cnt_per_site_hash {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{allele_cnt_per_site_hash} = $arg;
  }
  return $self->{allele_cnt_per_site_hash};
}

sub AF {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{AF} = $arg;
  }
  return $self->{AF};
}

sub sample_to_spop {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{sample_to_spop} = $arg;
  }
  return $self->{sample_to_spop};
}

sub gl_hash {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{gl_hash} = $arg;
  }
  return $self->{gl_hash};
}

sub gl_flag {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{gl_flag} = $arg;
  }
  return $self->{gl_flag};
}

sub per_sample_dp_flag {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{per_sample_dp_flag} = $arg;
  }
  return $self->{per_sample_dp_flag};
}

sub supp_dp_hash {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{supp_dp_hash} = $arg;
  }
  return $self->{supp_dp_hash};
}

sub complex_var_gl_hash {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{complex_var_gl_hash} = $arg;
  }
  return $self->{complex_var_gl_hash};
}

