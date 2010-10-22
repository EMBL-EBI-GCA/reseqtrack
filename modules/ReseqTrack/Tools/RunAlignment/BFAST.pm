package ReseqTrack::Tools::RunAlignment::BFAST;

use strict;
use warnings;
use vars qw(@ISA);
use Data::Dumper;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::RunAlignment;
use ReseqTrack::Tools::FileSystemUtils qw(delete_directory check_files_exists);

#use FastQ;
use File::Basename;
use File::stat;

@ISA = qw(ReseqTrack::Tools::RunAlignment);

sub new {

	my ( $class, @args ) = @_;
	my $self = $class->SUPER::new(@args);

	my ( $read_length, $program,
		$preprocess_exe, $verbose, $sub_sample, $max_bases, )
	  = rearrange(
		[
			qw(
			  READ_LENGTH
			  PROGRAM
			  PREPROCESS_EXE
			  VERBOSE
			  SUB_SAMPLE
			  MAX_BASES
			  )
		],
		@args
	  );

	my $def_max_bases = 1.2 * 100000000;
	$self->max_bases($def_max_bases);

	$self->max_bases($max_bases) if ($max_bases);
	$self->verbose($verbose);

	$self->retrieve_base_counts($sub_sample);

	$self->program($program) if ( !$self->program );
	$self->preprocess_exe($preprocess_exe);

	$self->set_read_length();
   

	print "Created BFAST object\n";
	$self->check_reference();
	
	$self->create_tmp_process_dir();
	$self->create_cmds();
	return $self;

}

sub create_cmds {
	my $self = shift;
	$self->set_options();
	$self->set_read_length() if ( !$self->read_length );
	#$self->set_reference();
	$self->create_fastq_cmd_line();
	$self->create_match_cmd_line();
	$self->create_localalign_cmd_line();
	$self->create_postprocess_cmd_line();
}

sub check_reference {
	my $self = shift;
	
	if ($self->{read_length} < 40 ){
		if ( ! ($self->{reference} =~ /short/i)){
		throw "Reference name looks wrong.No 'short' but read length < 40";
		}
	}
	
	   if ($self->{read_length} >= 40 ){
        if ( ! ($self->{reference} =~ /long/i)){
        throw "Reference name looks wrong.No 'long' but read length >=40";
        }
    }
    
	return;
}



sub subsample_fastq {
	my $self = shift;
	my $file = shift;

	my $seq_list;
	my $size = $self->max_bases;

	my $copy_fastq;
	my $outfile = $self->working_dir . '/' . $$ . basename($file);

	my $check = $self->files_to_subsample;

	if ( !defined $$check{$file} ) {
		print "$file does not need subsampling\n";

		# $seq_list = FastQ::sample( $file, $outfile,  $size);
		return $outfile;
	}

	print "sampling: $file to $size bases\n";    # if $verbose;

	my $tmp_size1 = stat($file)->size;
	$tmp_size1 /= ( 1024 * 1000 );
	printf "Size of org  ~ %6.3f Mb (bytes)\n", $tmp_size1;    # if $verbose;

	#$seq_list = FastQ::sample( $file, $outfile,  $size);#  if $run;
	#FastQ::sample( $file , $outfile, $size,  $seq_list);#  if $run;

	$tmp_size1 = stat($outfile)->size;
	$tmp_size1 /= ( 1024 * 1000 );
	printf "Size now     ~ %6.3f Mb (bytes)\n", $tmp_size1;
	print "Processing $outfile\n";

	return $outfile;

}

sub perform_subsample {
	my $self = shift;

	if ( $self->fragment_file ) {
		$self->fragment_file( $self->subsample_fastq( $self->fragment_file ) );
	}

	if ( $self->mate1_file ) {
		$self->mate1_file( $self->subsample_fastq( $self->mate1_file ) );
	}

	if ( $self->mate2_file ) {
		$self->mate2_file( $self->subsample_fastq( $self->mate2_file ) );
	}

	return;

}

sub run {
	my ($self) = @_;
	my $exit;

	if ( $self->should_subsample ) {
		$self->perform_subsample;
	}

	$self->create_cmds();
	$self->change_dir();

	print "Working in ", $self->working_dir(), "\n";

	exit;

	# 'match procedure, incorporates preprocess pipe
	eval {
		$exit = system( $self->match_cmd );
		if ( $exit && $exit >= 1 ) {
			throw( "Failed to run " . $self->match_cmd );
		}
	};
	if ($@) {
		throw("Failed to run bfast match: $@");
	}

	#localalign
	eval {
		$exit = system( $self->localalign_cmd );
		if ( $exit && $exit >= 1 ) {
			print STDERR ( "Failed to run " . $self->localalign_cmd );
		}
	};
	if ( $@ || $exit >= 1 ) {
		throw("Failed to run bfast localalign: $@");
	}

	#postprocess
	eval {
		$exit = system( $self->postprocess_cmd );
		if ( $exit && $exit >= 1 ) {
			print STDERR( "Failed to run " . $self->postprocess_cmd );
		}
	};
	if ( $@ || $exit > 1 ) {
		throw("Failed to bfast postprocess:  $@");
	}

	$self->create_bam();

	print "REALLY ??????? I made it to here . Wow";

	print "Delete tmp dir\n";

	#delete_directory( $self->working_dir() );

}

=head2 base_counts

=cut

sub retrieve_base_counts {

	my ( $self, $do_subsample ) = @_;

	my %base_counts;
	my %read_counts;
	my %files_to_subsample;

	if ( !( $self->input->isa('ReseqTrack::Collection') ) ) {
		throw "Only configged for collections at the moment";
	}

	my $others = $self->input->others;

    if ( ! $others){
        print "No 'others' associated.Cannot calc read length this way\n";
        return;
        
    }

	my @names;

	foreach my $other ($others) {

		foreach my $x (@$others) {

			#get file name
			#print $x->name,"\n";

			my $stats = $x->statistics;

			#get base counts for this file
			foreach my $y (@$stats) {
				if ( $y->{attribute_name} eq "base_count" ) {
			          printf "%-12d %s\n", $y->{attribute_value}, $x->name;
					$base_counts{ $x->name } = $y->{attribute_value};
				}
				
			 if ( $y->{attribute_name} eq "read_count" ) {
                      printf "%-12d %s\n", $y->{attribute_value}, $x->name;
                    $read_counts{ $x->name } = $y->{attribute_value};
                }	
				
			}
		}
	}

	foreach my $key ( keys %base_counts ) {
		my $b= $base_counts{$key};
            my $r= $read_counts{$key};
            my $length = $b/$r;
           
            
            if (! ($self->read_length) ){               
                $self->read_length ($length);
            }
		
		
		
		if ( $base_counts{$key} > $self->max_bases ) {
			my $count   = $base_counts{$key} - $self->max_bases;
			my $over_by = ( $count / $self->max_bases ) * 100;

			# arbitrary cut of 10x100MB
			# covers 85% of Illumin runs from Aug 2010 sequence.index
			# covers 100% of 454
			# covers > 97% of SOLID runs
			if ( $over_by > 33 ) {
				print $key , " exceeds max bases by ";
				printf "%6.2f", $over_by;
				print "%. sub sample\n";
				$self->{should_subsample} = 1;
				$files_to_subsample{$key} = 1;
			}
			
			
			
		}

	}
    print "Read length= ", $self->read_length if ($self->read_length);
	print "\n\n";
	$self->base_counts( \%base_counts );
	$self->files_to_subsample( \%files_to_subsample );
	return;
}

############
sub create_bam {
	my $self = shift;
	my $exit;

	my $sam = $self->sam();

	my $make_bam =
	  $self->samtools . "  import " . $self->reference . " $sam  $$.bam ";

	print "$make_bam\n";

	eval {
		$exit = system($make_bam );
		if ( $exit && $exit >= 1 ) {
			print STDERR( "Failed to run " . $make_bam );
		}
	};
	if ( $@ || $exit > 1 ) {
		throw("Failed to run $make_bam:  $@");
	}

	$self->bam("$$.bam");
	print "Made " . $self->bam, "\n";
}

sub set_read_length {
	my $self = shift;
	my $r_length;

	if ( $self->read_length ) {
		print "Read Length already set. Skipping\n";
		return;
	}

	print "Determining read length\n";

	# Assuming all reads are about the same length.
	if ( -e $self->fragment_file ) {

		my $frag = $self->fragment_file();

		open( my $FQ, "zcat $frag |" ) || die "Could not open: $frag";
		my $line_id   = <$FQ>;
		my $line_seq  = <$FQ>;
		my $line_sep  = <$FQ>;
		my $line_qual = <$FQ>;
		close($FQ);
		$r_length = length($line_seq);
		print "read length = $r_length \n" if $self->verbose;
		$self->read_length($r_length);
		return;
	}

	if ( -e $self->mate1_file() ) {
		my $frag = $self->mate1_file();
		open( my $FQ, "zcat $frag |" ) || die "Could not open: $frag";
		my $line_id   = <$FQ>;
		my $line_seq  = <$FQ>;
		my $line_sep  = <$FQ>;
		my $line_qual = <$FQ>;
		close($FQ);
		$r_length = length($line_seq);
		print "read length = $r_length \n";
		$self->read_length($r_length);
		return;
	}

	if ( -e $self->mate2_file() ) {
		my $frag = $self->mate2_file();
		open( my $FQ, "zcat $frag |" ) || die "Could not open: $frag";
		my $line_id   = <$FQ>;
		my $line_seq  = <$FQ>;
		my $line_sep  = <$FQ>;
		my $line_qual = <$FQ>;
		close($FQ);
		$r_length = length($line_seq);
		print "read length = $r_length \n";
		$self->read_length($r_length);
		return;
	}

	throw "Could not set read_length" if ( !$self->read_length() );

}

sub create_postprocess_cmd_line {
	my $self = shift;
	my $cmd_line;
	my $tmp_dir = $self->working_dir() . '/';

	$cmd_line = $self->program() . " postprocess ";
	$cmd_line .= " -f " . $self->reference() . " ";
	$cmd_line .= " -i  ${tmp_dir}$$.baf  ";
	$cmd_line .= " -n 4 -Q 1000 -t  > ${tmp_dir}$$.sam ";
	print "$cmd_line\n\n" if $self->verbose;
	$self->postprocess_cmd($cmd_line);
	$self->sam("$$.sam");
}

sub create_localalign_cmd_line {
	my $self = shift;
	my $cmd_line;

	my $tmp_dir = $self->working_dir() . '/';

	$cmd_line = $self->program() . " localalign ";
	$cmd_line .= " -f " . $self->reference() . " ";

	$cmd_line .= " -m  ${tmp_dir}$$.bmf ";

	$cmd_line .= " -A 1 -o 20 -n 4 -t > $tmp_dir$$.baf ";
	print "$cmd_line\n\n" if $self->verbose;
	$self->localalign_cmd($cmd_line);

}

sub create_match_cmd_line {
	my $self = shift;
	my $cmd_line;

	my $tmp_dir = $self->working_dir() . '/';

	$cmd_line = $self->fastq_cmd();
	$cmd_line .= $self->program() . " match ";
	$cmd_line .= " -f " . $self->reference() . " ";
	$cmd_line .= $self->ali_options() . " ";
	$cmd_line .= ">  ${tmp_dir}$$.bmf";

	#print "\n\n", $cmd_line, "\n";
	print "$cmd_line\n\n" if $self->verbose;
	$self->match_cmd($cmd_line);

}

sub create_fastq_cmd_line {
	my $self = shift;
	my $cmd_line;

	my $tmp_dir = $self->working_dir() . '/';

    throw ("No preprocess exe set") if (! (defined $self->{preprocess_exe}) );


	$cmd_line = $self->preprocess_exe() . " ";
	$cmd_line .= $self->mate1_file() . " " if ( defined $self->mate1_file );
	$cmd_line .= $self->mate2_file() . " " if ( defined $self->mate2_file );
	$cmd_line .= $self->fragment_file() . " "
	  if ( defined $self->fragment_file );

	$cmd_line .= " \| ";

	print "$cmd_line\n\n" if $self->verbose;

	$self->fastq_cmd($cmd_line);
}

sub set_options {
	my $self = shift;
	
	$self->ali_options('-A 1 -K 8 -M 384 -n 4 -Q 25000 ');
}
######################
sub base_counts {
	my ( $self, $arg ) = @_;

	if ($arg) {
		$self->{base_counts} = $arg;
	}
	return $self->{base_counts};
}

sub files_to_subsample {
	my ( $self, $arg ) = @_;

	if ($arg) {
		$self->{files_to_subsample} = $arg;
	}
	return $self->{files_to_subsample};
}

sub max_bases {
	my ( $self, $arg ) = @_;

	if ($arg) {
		$self->{max_bases} = $arg;
	}
	return $self->{max_bases};

}

sub should_subsample {
	my ( $self, $arg ) = @_;

	if ($arg) {
		$self->{should_subsample} = $arg;
	}
	return $self->{should_subsample};

}

#sub set_reference {
#	my $self = shift;
#
#	throw "No read length given" if ( !defined( $self->read_length() ) );
#
#	if ( $self->read_length() < 40 ) {
#		my $ref = $self->shreference;
#		$self->reference($ref);
#	}
#	else {
#		my $ref = $self->long_reference;
#		$self->reference($ref);
#	}
#}

sub ali_options {
	my ( $self, $arg ) = @_;

	if ($arg) {
		$self->{ali_options} = $arg;
	}
	return $self->{ali_options};
}

sub preprocess_exe {
	my ( $self, $arg ) = @_;

	if ($arg) {
		$self->{preprocess_exe} = $arg;
	}
	return $self->{preprocess_exe};
}

sub read_length {
	my ( $self, $arg ) = @_;

	if ($arg) {
		$self->{read_length} = $arg;
	}
	return $self->{read_length};
}

sub fastq_cmd {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->{fastq_cmd} = $arg;
	}
	return $self->{fastq_cmd};
}

sub match_cmd {
	my ( $self, $arg ) = @_;

	if ($arg) {
		$self->{match_cmd} = $arg;
	}
	return $self->{match_cmd};
}

sub localalign_cmd {
	my ( $self, $arg ) = @_;

	if ($arg) {
		$self->{localalign_cmd} = $arg;
	}
	return $self->{localalign_cmd};
}

sub postprocess_cmd {
	my ( $self, $arg ) = @_;

	if ($arg) {
		$self->{postprocess_cmd} = $arg;
	}
	return $self->{postprocess_cmd};
}

sub space {
	my ( $self, $arg ) = @_;

	if ($arg) {
		$self->{space} = $arg;
	}
	return $self->{space};
}

sub sam {
	my ( $self, $arg ) = @_;

	if ($arg) {
		$self->{sam} = $arg;
	}
	return $self->{sam};
}

sub bam {
	my ( $self, $arg ) = @_;

	if ($arg) {
		$self->{bam} = $arg;
	}
	return $self->{bam};
}

sub reference {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->{reference} = $arg;
	}
	return $self->{reference};
}

#sub short_reference {
#	my ( $self, $arg ) = @_;
#	if ($arg) {
#		$self->{short_reference} = $arg;
#	}
#	return $self->{short_reference};
#}

#sub long_reference {
#	my ( $self, $arg ) = @_;
#	if ($arg) {
#		$self->{long_reference} = $arg;
#	}
#	return $self->{long_reference};
#}

sub samtools {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->{samtools} = $arg;
	}
	return $self->{samtools};
}

sub snpbin {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->{snpbin} = $arg;
	}
	return $self->{snpbin};
}

sub program {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->{program} = $arg;
	}
	return $self->{program};
}

sub verbose {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->{verbose} = $arg;
	}
	return $self->{verbose};
}

1;

