package ReseqTrack::Tools::RunAlignment::BFAST;

use strict;
use warnings;
use vars qw(@ISA);
use Data::Dumper;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::RunAlignment;
use ReseqTrack::Tools::FileSystemUtils qw(delete_directory check_files_exists create_tmp_process_dir );

use File::Basename;
use File::stat;

@ISA = qw(ReseqTrack::Tools::RunAlignment);

sub new {

	my ( $class, @args ) = @_;
	my $self = $class->SUPER::new(@args);

	my (
		$read_length, $program,    $preprocess_exe,
		$verbose,     $sub_sample_size , $max_bases,
	  )
	  = rearrange(
		[
			qw(
			  READ_LENGTH
			  PROGRAM
			  PREPROCESS_EXE
			  VERBOSE			
			  MAX_BASES
			  )
		],
		@args
	  );

	$self->ali_options('-A 1 -K 8 -M 384 -n 4 -Q 25000 ');
	$self->space ("A -1"); # color space by default

	$self->read_length($read_length);

	$self->verbose($verbose);

	$self->preprocess_exe($preprocess_exe);



	return $self;

}
#################################################################


sub run {
	my ($self) = @_;
	my $exit;

	$self->change_dir();

       	$self->check_reference();

	$self->create_cmds();

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

	$self->convert_sam_to_bam();
	$self->get_bfast_tmp_files;


	return 0 ;
}

sub create_cmds {

    my $self = shift;
    $self->create_fastq_cmd_line();
    $self->create_match_cmd_line();
    $self->create_localalign_cmd_line();
    $self->create_postprocess_cmd_line();
    return;
}

sub check_reference {
    my $self = shift;
    print "Check:",$self->{reference},"\n";
    print "Check:",$self->{max_read_length},"\n";
 

    if ( $self->{max_read_length} < 40 ) {
        if ( !( $self->{reference} =~ /short/i ) ) {
            warning "Reference name looks wrong.No 'short' but read length < 40";
	    print "Automatically changing\n";
	    my $ref =  $self->{reference};
	    print $ref,"\n";
	    $ref =~ s/long/short/;
	    $self->reference($ref);
        }
    }

    if ( $self->{max_read_length} >= 40 ) {
        if ( !( $self->{reference} =~ /long/i ) ) {
	  warning "Reference name looks wrong.No 'long' but read length < 40";
	  print "Automatically changing\n";
	   my $ref =  $self->{reference};
	   $ref =~ s/short/long/;
	   print $ref,"\n";
	   $self->reference($ref);
        }
    }
    print "reference looks OK\n";
    return;
}



############

sub convert_sam_to_bam {
	my $self = shift;
	my $exit;

	my $sam = $self->sam();

	print "Converting sam file to bam format\n";

	my $out_bam = $self->working_dir() . '/'. $$ . '.bam';


	$out_bam =~ s/\/\//\//;
	#print $out_bam,"\n";

	my $make_bam =
	  $self->samtools . " import " . $self->reference . " $sam  $out_bam";

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

	$self->bam($out_bam);
        $self->output_files ($out_bam);
	if (-e	$self->bam) {
	  print "Made " . $self->bam, "\n";
	}
	else{
	  throw "Failed to create $out_bam";
	}

	return;
}

sub create_postprocess_cmd_line {
	my $self = shift;
	my $cmd_line;
	my $tmp_dir = $self->working_dir() . '/';

	$cmd_line = $self->program() . " postprocess ";
	$cmd_line .= " -f " . $self->reference() . " ";
	my $tmp_baf =  $tmp_dir . $$ . ".baf"; 
	$cmd_line .= " -i  $tmp_baf ";
	my $tmp_sam =  $tmp_dir . $$ . ".sam"; 
	$cmd_line .= " -n 4 -Q 1000 -t -U  > $tmp_sam";
	print "$cmd_line\n\n";# if $self->verbose;
	$self->postprocess_cmd($cmd_line);
	$self->sam("$tmp_sam");
	return;
}

sub create_localalign_cmd_line {
	my $self = shift;
	my $cmd_line;

	my $tmp_dir = $self->working_dir() . '/';

	$cmd_line = $self->program() . " localalign ";
	$cmd_line .= " -f " . $self->reference() . " ";
	
	my $tmp_bmf =  $tmp_dir . $$ . ".bmf";
	$cmd_line .= " -m  $tmp_bmf ";

	my $tmp_baf =  $tmp_dir . $$ . ".baf"; 
	$cmd_line .= " -A 1 -o 20 -n 4 -t > $tmp_baf ";
	
	$self->files_to_delete ("$tmp_baf");
	
	print "$cmd_line\n\n";#  if $self->verbose;
	
	$self->localalign_cmd($cmd_line);
    return;
}

sub create_match_cmd_line {
	my $self = shift;
	my $cmd_line;

	my $tmp_dir = $self->working_dir() . '/';
	$tmp_dir =~ s/\/\//\//;
	$cmd_line = $self->fastq_cmd();
	$cmd_line .= $self->program() . " match ";
	$cmd_line .= " -f " . $self->reference() . " ";
	$cmd_line .= $self->ali_options() . " ";
	my $tmp_bmf =  $tmp_dir . $$. ".bmf";
	$cmd_line .= ">  $tmp_bmf";
	$self->localalign_cmd($cmd_line);
    

	$self->files_to_delete ("$tmp_bmf");

	print "$cmd_line\n\n";# if $self->verbose;
	$self->match_cmd($cmd_line);
    return;
}

sub create_fastq_cmd_line {
	my $self = shift;
	my $cmd_line;

	my $tmp_dir = $self->working_dir() . '/';

	throw("No preprocess exe set") if ( !( defined $self->{preprocess_exe} ) );

	$cmd_line = $self->preprocess_exe() . " ";

	if (! $self->skip_mate_files){

	  $cmd_line .= $self->mate1_file() . " " if ( $self->mate1_file );
	  $cmd_line .= $self->mate2_file() . " " if ( $self->mate2_file );
	}

	if (! $self->skip_fragment){
	  $cmd_line .= $self->fragment_file() . " "
	    if ( $self->fragment_file );
	}

	$cmd_line .= " \| ";

	print "$cmd_line\n\n"  if $self->verbose;

	$self->fastq_cmd($cmd_line);
	return;
}

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

sub verbose {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->{verbose} = $arg;
	}
	return $self->{verbose};
}


sub get_bfast_tmp_files {
  my ( $self ) = @_;

  my $dir = $self->working_dir;
  my $tmp_files_found = 0; 
  print $dir,"\n";

  if ( !(-d $dir) ||  !( -e $dir) ||  !( $dir =~ /^\// ) ){
	warning "Bad directory info for bfast tmp file search.Skipping\n";
	return;
      }


  opendir ( DIR, $dir ) || die "Error in opening dir $dir\n";

  my @files = readdir(DIR);

  foreach my $file (@files) {
    if ( $file =~  /\.bfast\.tmp/ ){
      print $file,"\n";
      my $del_file = $dir . '/'. $file;
      $del_file =~ s/\/\//\//;
      $self->file_to_delete ($del_file);
      $tmp_files_found ++;
    }
  }
 
  print "Number of temporary bfast files found = $tmp_files_found\n";
  
  return;
}


1;

