package ReseqTrack::Tools::QAandCntFastq;
  
# This package is to QA fastq files and calculate the base and reads counts
# Holly Zheng, EBI 1000 Genomes Project, August 2009

use Exporter;
@ISA 	= qw(Exporter);
@EXPORT	= qw(qa_and_cnt getReadLen);

use strict;
use warnings;
use IO::Handle;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);
use File::Basename;
use ReseqTrack::Tools::Exception;
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::Tools::Counts;

sub qa_and_cnt {

	my ($files_ref, $instrument, $output_dir, $len_limit, $run_id) = @_;

	throw("No instrument information is available for files\n") if (!$instrument);
	
	my %files = %$files_ref;
   
	my $file_num = keys(%files);
	my (
		$mate1, 
		$mate2, 
		$frag, 
		$mate1_fh, 
		$mate2_fh, 
		$frag_fh, 
		$mate1_file, 
		$mate2_file,
		$frag_file, 
		$mate1_stats, 
		$mate2_stats, 
		$frag_stats,
		$ori_mate1_obj, 
		$ori_mate2_obj, 
		$ori_frag_obj, 
		);
	
	my $log = $$output_dir . "/" . $run_id . ".log";
        throw($$output_dir." does not exist") unless (-e $$output_dir);
	open (my $log_fh,  ">", $log) || throw("cannot open log file $log\n");
	
	if ( $file_num == 3 ) {

		my ($tmp_frag_file, $frag_fh_tmp, $frag_stats_tmp);
		
		foreach my $f (keys %files) {
			if ($f =~ /\_1.fastq/) {
				$mate1 = $f;
				$ori_mate1_obj = $files{$f};
				($mate1_fh, $mate1_file) = open_fh($output_dir, \$mate1); 
				
			}
			elsif ( $f =~ /\_2.fastq/) {
				$mate2 = $f;
				$ori_mate2_obj = $files{$f};
				($mate2_fh, $mate2_file) = open_fh($output_dir, \$mate2); 
			}
			elsif ( $f =~ /.fastq/) {
				$frag = $f;
				$ori_frag_obj = $files{$f};
				($frag_fh, $frag_file) = open_fh($output_dir, \$frag); 
				
				$tmp_frag_file = $$frag_file . ".tmp"; #$tmp_frag_file is a full file path
				open ($frag_fh_tmp, ">", $tmp_frag_file) || throw("Cannot open temperary output fastq file $tmp_frag_file for fragments\n");
				
			}
		}					
		
		throw("This run supposes to have three files. Either that or file names do not fit convention\n") if(!$mate1 || !$mate2 || !$frag);
		
		($mate1_stats, $mate2_stats, $frag_stats_tmp) = check_paired_fastq(\$mate1, \$mate2, \$instrument, $frag_fh_tmp, $mate1_fh, $mate2_fh, 
$len_limit, $run_id, $log_fh);	
	    
	    #### cat together the original fragment file with the fragment file generated from singled out mate reads
	    `cat $frag $tmp_frag_file > ./tmp.concatenated.fastq`;
	    
	    #print "original frag file: " . `ls -l $frag` . "\n";
	    #print "temparay frag file: " . `ls -l $tmp_frag_file` . "\n";
	    #print "concatenated frag file: " . `ls -l ./tmp.concatenated.fastq` . "\n";
	    
	    `mv ./tmp.concatenated.fastq $$frag_file`;

	    #print "new frag file: " . `ls -l $$frag_file` . "\n";
	    
	    open (my $frag_fh_new, ">", "./tmp.filtered.concatenated.fastq") || throw("Cannot open output filtered fastq file\n");
		
		#### QA the concatenated fragment file, output fragment file will be in ./tmp.filtered.concatenated.fastq
		$frag_stats = check_single_fastq($frag_file, \$instrument, $frag_fh_new, $len_limit, $run_id, $log_fh);
		`mv ./tmp.filtered.concatenated.fastq $$frag_file`;
		`rm $tmp_frag_file`;
		close $frag_fh_new;
	}
	elsif ($file_num == 2 ) {
		foreach my $f (keys %files) {
			if ($f =~ /\_1.fastq/) {
				$mate1 = $f;
				$ori_mate1_obj = $files{$f};
				($mate1_fh, $mate1_file) = open_fh($output_dir, \$mate1); 
			}
			elsif ( $f =~ /\_2.fastq/) {
				$mate2 = $f;
				$ori_mate2_obj = $files{$f};
				($mate2_fh, $mate2_file) = open_fh($output_dir, \$mate2); 
			}
		}					
		
		my ($tmp1, $tmp2) = split (/\_1.fastq/, $mate1);
		my $frag = $tmp1 . ".fastq";
		($frag_fh, $frag_file) = open_fh($output_dir, \$frag); 	
	       
		throw("This run supposes to have two files.  File names do not fit convention\n") if(!$mate1 || !$mate2);
		($mate1_stats, $mate2_stats, $frag_stats) = check_paired_fastq(\$mate1, \$mate2, \$instrument, $frag_fh, $mate1_fh, $mate2_fh, $len_limit, $run_id, $log_fh);
	}			
	elsif ( $file_num == 1 ) {
		my @array = keys %files;
		my $f = $array[0];
		$ori_frag_obj = $files{$f};
		($frag_fh, $frag_file) = open_fh($output_dir, \$f);
		$frag_stats = check_single_fastq(\$f, \$instrument, $frag_fh, $len_limit, $run_id, $log_fh);
	}	
	else {
		throw("The # of fastq files associated with the specified run_id is not 1, 2 or 3\n");
	}	
	
	#print "For mate1 $$mate1_file, stats is " . $mate1_stats->filt_read_cnt . "\n";
	
	my %results; #key is outfile path, value is an anonymous array, first element is the Counts object, second element is the original file object
	if ($mate1_file && $mate1_stats) { $results{$$mate1_file}=[$mate1_stats, $ori_mate1_obj] }; #anonymous array, need to de-reference later
	if ($mate2_file && $mate2_stats) { $results{$$mate2_file}=[$mate2_stats, $ori_mate2_obj] };
	if ($frag_file && $frag_stats)  { $results{$$frag_file}=[$frag_stats, $ori_frag_obj] };
	
	close $frag_fh if ($frag_fh);
	close $mate1_fh if ($mate1_fh);
	close $mate2_fh if ($mate2_fh);
	
	return (\%results);	 #### key is outfile path, value is an array, first element is the Counts object, second element is the original file object
}

sub open_fh { #### open output handle 
	my ($output_dir, $infile_path) = @_;
	
	my $file = basename($$infile_path);
	if ($file =~ /.gz/) {
		$file =~ s/.fastq.gz/.filt.fastq/i;
	}
	elsif ( $file =~ /.fastq$/ ) {	
		$file =~ s/.fastq/.filt.fastq/i;
	} 
	my $output_path = $$output_dir . "/" . $file;
	open (my $fh, ">", $output_path) || throw("cannot open output file $output_path\n");
	return ($fh, \$output_path);
}	

sub open_infile_handle {
	my ($file) = @_;
	my ($FH, $line_cnt);
	
	if ($$file =~ /.gz/ ) {
		open($FH, "gunzip -c $$file |") or throw( "could not open $$file" );	
		$line_cnt = `gunzip -c $$file | wc -l`;		
	}
	else {
		open($FH, "<", $$file) or throw( "could not open $$file" );
		$line_cnt = `cat $$file|wc -l`;
	}
	chomp $line_cnt;
	return ($FH, $line_cnt);	
}

sub check_paired_fastq {
	my ($mate1, $mate2, $instrument, $FH_frag, $OUT_mate1, $OUT_mate2, $len_limit, $run_id, $log_fh) = @_;
	
	###### open file handles and check if the two files have same number of reads ##########	
	my ($FILE1, $mate1_line_cnt) = open_infile_handle($mate1);
	my ($FILE2, $mate2_line_cnt) = open_infile_handle($mate2);
	
	throw("mate1 and mate2 do not have identical number of lines\n") if ($mate1_line_cnt != $mate2_line_cnt);
	#print "line cnts for two files are: $mate1_line_cnt and $mate2_line_cnt\n";
		
	my $stats1 = ReseqTrack::Tools::Counts->new(); 
	my $stats2 = ReseqTrack::Tools::Counts->new(); 
	my $stats_frag = ReseqTrack::Tools::Counts->new();
			
	my $flt_read_cnt1 = 0;
	my $flt_base_cnt1 = 0;
	my $unflt_read_cnt1 = 0;
	my $unflt_base_cnt1 = 0;
	my $flt_read_cnt2 = 0;
	my $flt_base_cnt2 = 0;
	my $unflt_read_cnt2 = 0;
	my $unflt_base_cnt2 = 0;
	my $frag_read_cnt=0;
	my $frag_base_cnt=0;
		
	while ( my $line1 = <$FILE1> ) {
		my $line2 = <$FILE1>;
		my $line3 = <$FILE1>;
		my $line4 = <$FILE1>;
		
		my $line_a = <$FILE2>;
		my $line_b = <$FILE2>;
		my $line_c = <$FILE2>;
		my $line_d = <$FILE2>;
			
		chomp ($line1, $line2, $line3, $line4, $line_a, $line_b, $line_c, $line_d);
		$line1 =~ s/^\s+|\s+$//g;
		$line2 =~ s/^\s+|\s+$//g;
		$line3 =~ s/^\s+|\s+$//g;
		$line4 =~ s/^\s+|\s+$//g;
	
		$line_a =~ s/^\s+|\s+$//g;
		$line_b =~ s/^\s+|\s+$//g;
		$line_c =~ s/^\s+|\s+$//g;
		$line_d =~ s/^\s+|\s+$//g;

=head		
		print "processed mate1 line 1: $line1\n";
		print "processed mate2 line a: $line_a\n";
=cut	
		my $mate1_read_len = getReadLen(\$line2, $instrument);
		my $mate2_read_len = getReadLen(\$line_b, $instrument);	
			
		$unflt_read_cnt1++;
		$unflt_base_cnt1 = $unflt_base_cnt1 + $$mate1_read_len;
		
		$unflt_read_cnt2++;
		$unflt_base_cnt2 = $unflt_base_cnt2 + $$mate2_read_len;

=head						
		print "mate1 read len: $$mate1_read_len\n";
		print "mate1 unfilter read cnt: $unflt_read_cnt1\n";
		print "mate1 unfilter base cnt: $unflt_base_cnt1\n";
		
		print "mate2 read len: $$mate2_read_len\n";
		print "mate2 unfilter read cnt: $unflt_read_cnt2\n";
		print "mate2 unfilter base cnt: $unflt_base_cnt2\n";	
=cut	
	
		if ( QA(\$line1, \$line2, \$line3, \$line4, \$unflt_read_cnt1, $instrument, $len_limit, $run_id, $log_fh) eq "Fail" ) {
		    if ( QA(\$line_a, \$line_b, \$line_c, \$line_d, \$unflt_read_cnt2, $instrument, $len_limit, $run_id, $log_fh) eq "Fail") {
		        print $log_fh "ERROR: both mates failed syntax and sequence check\n\n";
			    next;
		    }
		    else {
		     	print $log_fh "ERROR: mate1 failed mate2 passed, move mate2 to fragment file\n\n";
		        print  $FH_frag "$line_a\n$line_b\n$line_c\n$line_d\n";
		        $frag_read_cnt++;
		        $frag_base_cnt = $frag_base_cnt + $$mate2_read_len;
		    	next;
		    }       
		}
		else {
			# case where mate1 pass and mate 2 failed, move mate1 to fragment file
		  	if ( QA(\$line_a, \$line_b, \$line_c, \$line_d, \$unflt_read_cnt2,  $instrument, $len_limit, $run_id, $log_fh) eq "Fail") {
		  	    print  $FH_frag "$line1\n$line2\n$line3\n$line4\n";
			    print $log_fh "ERROR: mate1 passed mate2 failed, moved mate1 to fragment file\n\n";
			    $frag_read_cnt++;
		        $frag_base_cnt = $frag_base_cnt + $$mate1_read_len;
			    next;
			}
		}
		
		# anything that passes above checks will come here:	 	
		$flt_read_cnt1++;
		$flt_base_cnt1 = $flt_base_cnt1 + $$mate1_read_len;
		$flt_read_cnt2++;
		$flt_base_cnt2 = $flt_base_cnt2 + $$mate2_read_len;
=head		
		print "mate1 filtered read cnt: $flt_read_cnt1\n";
		print "mate1 filter base cnt: $flt_base_cnt1\n";
		print "mate2 filtered read cnt: $flt_read_cnt2\n";
		print "mate2 filter base cnt: $flt_base_cnt2\n";	
=cut		
		print $OUT_mate1 "$line1\n$line2\n$line3\n$line4\n";	
		print $OUT_mate2 "$line_a\n$line_b\n$line_c\n$line_d\n";	
	}
	
	$stats1->read_cnt($unflt_read_cnt1);
	$stats1->base_cnt($unflt_base_cnt1);
	
	$stats2->read_cnt($unflt_read_cnt2);
	$stats2->base_cnt($unflt_base_cnt2);
		
	$stats1->filt_read_cnt($flt_read_cnt1);
	$stats1->filt_base_cnt($flt_base_cnt1);

	$stats2->filt_read_cnt($flt_read_cnt2);
	$stats2->filt_base_cnt($flt_base_cnt2);
				
	$stats_frag->filt_read_cnt($frag_read_cnt);
	$stats_frag->filt_base_cnt($frag_base_cnt);

	close $FILE1;
	close $FILE2;
	return ($stats1, $stats2, $stats_frag);	
}	

sub check_single_fastq {
	my ($file, $instrument,  $OUT_fastq, $len_limit, $run_id, $log_fh) = @_;
	
	my ($FILE1, $frag_file_line_cnt) = open_infile_handle($file);
	
	my $flt_read_cnt = 0;
	my $flt_base_cnt = 0;
	my $unflt_read_cnt = 0;
	my $unflt_base_cnt = 0;

	my $stats = ReseqTrack::Tools::Counts->new();

	while ( my $line1 = <$FILE1> ) {
		my $line2 = <$FILE1>;
		my $line3 = <$FILE1>;
		my $line4 = <$FILE1>;
			
		chomp ($line1, $line2, $line3, $line4);
		$line1 =~ s/^\s+|\s+$//g;
		$line2 =~ s/^\s+|\s+$//g;
		$line3 =~ s/^\s+|\s+$//g;
		$line4 =~ s/^\s+|\s+$//g;
		
		#print "processed line 1: $line1\n";
			
		my $read_len = getReadLen(\$line2, $instrument);
			
		$unflt_read_cnt++;
		$unflt_base_cnt = $unflt_base_cnt + $$read_len;
=head				
		print "read len: $$read_len\n";
		print "unfilter read cnt: $unflt_read_cnt\n";
		print "unfilter base cnt: $unflt_base_cnt\n";
=cut				
		if ( QA(\$line1, \$line2, \$line3, \$line4, \$unflt_read_cnt,  $instrument, $len_limit, $run_id, $log_fh) eq "Fail" ) {		
		    print $log_fh "ERROR: fragment $unflt_read_cnt failed QA\n";
		    next;
		}    
		 	
		$flt_read_cnt++;
		$flt_base_cnt = $flt_base_cnt + $$read_len;		
=head	
		print "read len: $$read_len\n";
		print "filtered read cnt: $flt_read_cnt\n";
		print "filter base cnt: $flt_base_cnt\n\n";					
=cut		
		print $OUT_fastq "$line1\n$line2\n$line3\n$line4\n";
	}

	$stats->read_cnt($unflt_read_cnt);
	$stats->base_cnt($unflt_base_cnt);
				
	$stats->filt_read_cnt($flt_read_cnt);
	$stats->filt_base_cnt($flt_base_cnt);
	
	close $FILE1;
	return $stats;
}		
	
sub QA {
	my ($line_1_r, $line_2_r, $line_3_r, $line_4_r, $read_num_r,  $instrument, $len_limit, $run_id, $log_fh) = @_;
	my $line_1 = $$line_1_r;
	my $line_2 = $$line_2_r;
	my $line_3 = $$line_3_r;
	my $line_4 = $$line_4_r;
	my $read_num = $$read_num_r;
	
	my $flag = "Pass";
	
	###### check if a read block start with @ and if the 3rd line start with a + #####
	unless ( $line_1 =~ /^\@/ && $line_3 =~ /^\+/ ) {
		throw("FATAL ERROR: block $read_num, the first line does not start with @ or the third line does not start with a +\n");
	}	
	
	unless ( $line_1 =~ /$run_id/ ) {
		warning("read name does not contain run_id\n");
		print $log_fh "ERROR: in Block $read_num, the read name $line_1 does not contain run id $run_id\n";
		$flag = "Fail";
		goto END_OF_QA;
	}	
		
	##### check if the read id on line1 and line_3 agree; if line_3 has a read id. ######
	$line_1 =~ s/^@//;
	$line_1 =~ s/^\s+|\s+$//;

	my $read_id1 = $line_1;
	
	my $read_id2;
	if (length($line_3) >1) {
		$line_3 =~ s/^\+//g;
		$line_3 =~ s/^\s+|\s+$//;
		$read_id2 = $line_3;
		
		if ($read_id1 ne $read_id2) {
			print $log_fh "ERROR: in Block $read_num, line_1 and line_3 do not have identical read id\n";
			$flag = "Fail";
		}
	}	
	
	##### check to see if the sequence line and the quality line has the same number of bits #####
	if (length($line_2) != length($line_4)) {
		$flag = "Fail";
		print $log_fh "ERROR: in Block $read_num, the sequence length is not equal to the quality length\n";
	}	
		
	##### If a read pass all above check, go on checking more below ######
	##### check if the seq length is > 35bp for Solexa, 25bp for solid, 30bp for 454 #######
	##### check to see if the first 25,30 or 35 bp contain no Ns and are not single type of base #####	

	if ($flag eq "Pass") {	
		if ($$instrument eq "ABI_SOLID") {

                  if (!$len_limit) { # default len limit is 25 bp
				$flag = runSolidSeqCheck($line_2, $line_4, $read_num,  $instrument, 25, $log_fh);
			}
			else {
				$flag = runSolidSeqCheck($line_2, $line_4, $read_num,  $instrument, $len_limit, $log_fh);
			}
			
		}	
		elsif ($$instrument eq "ILLUMINA") {
			if (!$len_limit) {
				$flag = runSeqCheck($line_2, $line_4, $read_num,  $instrument, 35, $log_fh);
			}
			else {
				$flag = runSeqCheck($line_2, $line_4, $read_num,  $instrument, $len_limit, $log_fh);
			}
		}
		elsif ($$instrument eq "LS454") {
			if (!$len_limit) {
				$flag = runSeqCheck($line_2, $line_4, $read_num,  $instrument, 30, $log_fh);
			}
			else {
				$flag = runSeqCheck($line_2, $line_4, $read_num,  $instrument, $len_limit, $log_fh);
			}
		}	
	}	
	
	##### If the seq passes QA so far, check % of Ns in the whole length of the read ######
	if ($flag eq "Pass") {	
		if ($$instrument eq "ABI_SOLID") {
			$flag = tooManyNsFullLenSolid($line_2, $read_num, $log_fh);
		}
		else {
			$flag = tooManyNsFullLen_OddChar($line_2, $read_num, $log_fh);
		}
	}
				
	END_OF_QA:
	return $flag;
}	

sub tooManyNsFullLenSolid {
	
	my ($seq, $read_num, $log_fh) = @_;
	my $flag;
	my $N_cnt = 0;
	my $read_length = length($seq) -1;
			
	my @characters = split (//, $seq); #split string into characters
	foreach my $char (@characters) {
		$N_cnt++ if ($char eq "\." || $char eq "N");
	}

	if ($N_cnt/$read_length > 0.5) {
		$flag = "Fail";
		print $log_fh "ERROR: The colorspace string contains more than 50% Ns in the whole length of the read $read_num\n";
	}
	return $flag;
}

sub tooManyNsFullLen_OddChar {
	my ($seq,  $read_num, $log_fh) = @_;
	my $flag;
	my $N_cnt = 0;
	my $read_length = length($seq);
	
	my @characters = split (//, $seq); #split string into characters
	foreach my $char (@characters) {
		$N_cnt++ if ($char eq "N");
	}
	
	if ($N_cnt/$read_length > 0.5) {
		$flag = "Fail";
		print $log_fh "ERROR: The seq contains more than 50% Ns in the whole length of the read $read_num\n";
	}
	
	if ( $seq =~ /[^ATGCN]/i) { #match anything other than A, T, G, C, N
		print $log_fh "ERROR: The seq contains characters other than ATGCN in its full length\n";
		$flag = "Fail";
	}
		
	return $flag;
}

sub runSolidSeqCheck {
	my ($line_2, $line_4, $read_num,  $instrument, $len_limit, $log_fh) = @_;
	my $flag = "Pass";
	my $read_len = length($line_2) -1;
	$flag = "Fail" if (checkSeqLen($read_len, $len_limit, $read_num, $log_fh) eq "Too short");
	my $seq = substr($line_2, 0, $len_limit+1);
	my $qual = substr($line_4, 1, $len_limit);
	#print "len limit is $len_limit\n";
	$flag = "Fail" if (checkSolidSeq($seq, $len_limit, $read_num, $log_fh) eq "Bad");
	$flag = "Fail" if ( checkQualityScore(\$qual, $instrument, $read_num, $log_fh) );
	return $flag;
}

sub runSeqCheck {
	my ($line_2, $line_4, $read_num,  $instrument, $len_limit, $log_fh) = @_;
	my $flag = "Pass";
	my $read_len = length($line_2);
	$flag = "Fail" if (checkSeqLen($read_len, $len_limit, $read_num, $log_fh) eq "Too short");
	my $seq = substr($line_2, 0, $len_limit);
	my $qual = substr($line_4, 0, $len_limit);
	#print "len limit is $len_limit\n";
	$flag = "Fail" if (checkSeq($seq, $len_limit, $read_num, $log_fh) eq "Bad");
	$flag = "Fail" if ( checkQualityScore(\$qual, $instrument, $read_num, $log_fh) );
	return $flag;
}	

sub checkQualityScore {
	my ($q_string, $instrument, $block, $log_fh) = @_;
	my $flag = 0;
	my $len = length($$q_string);
	for (my $i=0; $i < $len; $i++) {
		my $q_num;
		my $q_char = substr($$q_string, $i, 1);
		#if ($$instrument eq "ILLUMINA") {  #the short read archives have changed the quality score when fastq files were genreated
		#	$q_num = 10 * log(1 + 10 ** (ord($q_char) - 64) / 10.0) / log(10); # for illumina
			#print "quality score for illumina: $q_num\n";
		#}
		#else {
                #print "QUAL character $q_char\n";
		$q_num = ord($q_char) - 33; # for others	
		#print "quality char: $q_char\tquality score: $q_num\n";
		#}		
                
		$flag = 1 if ($q_num < 2);		
	}		
	print $log_fh "ERROR: in block $block, the quality value for one or more base is lower than 2 in the first $len bp region\n" if ($flag ==1);
	return $flag;		
}

sub checkSolidSeq {
	my ($seq, $len, $block, $log_fh) = @_;
	my $tag = "Good";
	
	########## test to see if the sequence line start with a letter and end with a string of numbers ####
	my $color_space = substr($seq, 1);
	unless ($seq =~ /^[A-Z]/i && looks_like_number($color_space) ) {
		$tag = "Bad"; 
		print $log_fh "ERROR: in Block $block, the sequence is not represented by proper color space code\n";
		#print $log_fh "ERROR: Here is the first 26 positions of the line:$seq\n"; 
	}
	
	######### test if the sequence in the first 25bp contains any N or are a single type of base such as AAAAAAAAAAA #####
	if ($seq =~ /\./ || $seq =~ /N/i ) { #BGI solid fastq have N as ambiguous base
		$tag = "Bad";
		print $log_fh "ERROR: The colorspace string contains Ns in the first $len bp\n";
	}
	
	if ($color_space =~ /0{$len}/ ) {
		$tag = "Bad";	
		print $log_fh "ERROR: The colorspace string contains only one type of nucleotide in the first $len bp\n";		
	}
	return $tag;	
}
	
sub checkSeq {
	my ($seq, $len, $block, $log_fh) = @_;
	my $tag = "Good";
	#print "IN checkSeq block: $seq\n";
	
	if ($seq =~ /N/i) {
		print $log_fh "ERROR: in Block $block, the sequence contains Ns in the first $len bp region\n";
		$tag = "Bad";
	}
	elsif ($seq =~ /A{$len}/i ){ 
		print $log_fh "ERROR: in Block $block, the sequence contains only AAAAAAAAAAAAA in the first $len bp\n";
		$tag = "Bad";
	}
	elsif ($seq =~ /T{$len}/i ) {
		print $log_fh "ERROR: in Block $block, the sequence contains only TTTTTTTTTTTT in the first $len bp\n";
		$tag = "Bad";
	}
	elsif ($seq =~ /G{$len}/i ) {
		print $log_fh "ERROR: in Block $block, the sequence contains only GGGGGGGGGGGGGG in the first $len bp\n";
		$tag = "Bad";
	}
	elsif ($seq =~ /C{$len}/i ) {
		print $log_fh "ERROR: in Block $block, the sequence contains only CCCCCCCCCCCCCCC in the first $len bp\n";
		$tag = "Bad";
	}
	return $tag;		 
}

sub checkSeqLen {
	my ($len, $min_len, $block, $log_fh) = @_;
	my $tag = "good";
	if ($len < $min_len) {
			print $log_fh "ERROR: in Block $block, the sequence is not longer than the minimal length $min_len bp\n";
			$tag = "Too short";
	}
	return $tag;
}		
	

sub getReadLen {
	my ($seq, $instrument) = @_;
	my $len;

	#color space
	if ($$instrument eq "ABI_SOLID") {
		$len = length($$seq) -1;		
	}
	#base space
	else {
		$len = length($$seq);
	}
	return \$len;
}	
	
1;
