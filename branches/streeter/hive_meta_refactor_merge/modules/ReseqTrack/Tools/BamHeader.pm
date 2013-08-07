package ReseqTrack::Tools::BamHeader;
=pod

=head1 NAME

ReseqTrack::Tools::BamUtils::BamHeader

=head1 SYNOPSIS

	This is a module should provide facility to add headers to bam files
	and append headers with appropriate @RG, @PG and @CO header lines.
	Header information can be extracted from a text file 
	(eg @HD and @SQ for a given reference genome) or an existing bam file. 
	Any information appended should be supplied via a hash ref in the 
	standard samtools format.
	The header lines should be output following the tag order
	my @tag_order = ( '@HD', '@SQ', '@RG', '@PG', '@CO' ); 
	
=cut

use strict;
use warnings;

use Data::Dumper;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename;

use ReseqTrack::Tools::RunProgram;
use ReseqTrack::Tools::RunAlignment;
use vars qw(@ISA); 

@ISA = qw(ReseqTrack::Tools::RunAlignment);

sub new {
	my ( $class, @args ) = @_;
	my $self = $class->SUPER::new(@args);

	my ( $header_file, $info_hash, $add_header, $header,
		 $delete_old_bams ,$bam) =
	  rearrange(
		[
			qw(
			  HEADER_FILE
			  INFO_HASH
			  ADD_HEADER
			  HEADER
			  DELETE_OLD_BAMS
			  BAM
			  )
		],
		@args
	  );

	$self->header_file($header_file);
	$self->info_hash($info_hash);
	$self->header($header);
	$self->add_header($add_header);
	$self->delete_old_bams($delete_old_bams);
	$self->bam($bam);
	
	return $self;
}




=head2 reheader_bam_file 

  Arg [1]   : ReseqTrack::Tools::BamUtils::BamHeader
  Arg [2]   : Path to bam header text file
  Function  : replace header of bam file with contents of text file
  Returntype: Path to reheaded bam file
  Exceptions: 
  Example   : $BH->reheader_bam_file($BH->tmp_header_file());
  Example   : $BH->reheader_bam_file($new_header_file); 
=cut

sub reheader_bam_file {
	my ( $self, $header_file ) = @_;

	my $aln_file = $self->bam;
	my $samtools = $self->samtools;

	throw("Could not find:$header_file") if ( !-e $header_file );

	print "Reheadering bam using $header_file\n";

	my $new_bam = $aln_file . ".$$.bam";
	my $cmd     = "$samtools reheader $header_file $aln_file > $new_bam";

	print $cmd, "\n";
	eval { `$cmd`; };

	if ($@) {
		throw("Failed: $cmd");
	}
	print "See $new_bam\n";

	if ( $self->delete_old_bams ) {
		$self->files_to_delete($aln_file);
		$self->delete_files() if $self->delete_old_bams ;
	}

	$self->bam($new_bam);

	return $new_bam;
}


=head2 get_bam_header_from_text_file

  Arg [1]   : ReseqTrack::Tools::BamUtils::BamHeader
  Arg [2]   : Optional path to header file, if not already specified
  Function  : Read text file into array to be processed for adding to bam 
  Returntype: None
  Exceptions: 
  Example   : $BH->get_bam_header_from_text_file($my_generic_header_file); 
=cut

sub get_bam_header_from_text_file {
	my ( $self, $header_file ) = @_;

	if ( !$header_file ) {
		$header_file = $self->header_file;
	}

	if ( !-e $header_file ) {
		throw "Reference sequence header file does not exist:$header_file";
	}

	open my $IN, '<', $header_file || throw "Failed to open bam header file";
	my @header_lines = <$IN>;
	close($IN);

	print "Got ", scalar(@header_lines), " lines from header file\n";

	$self->header_lines( \@header_lines );

	return;
}


=head2 get_header_from_bam

  Arg [1]   : ReseqTrack::Tools::BamUtils::BamHeader
  Arg [2]   : Optional path to bam, if not already specified
  Function  : Use samtools to extract header from bam into array to processing 
  Returntype: None
  Exceptions: Bam does not exist
  Example   : $BH->get_header_from_bam( ""my.bam");   
=cut

sub get_header_from_bam {
	my ( $self, $aln_file ) = @_;

	if ( !$aln_file ) {
		$aln_file = $self->bam;
	}

	if ( !-e $aln_file ) {
		throw "Alignment file not exist:$aln_file";
	}

	my $samtools = $self->samtools;

	my @header_lines;

	@header_lines = `$samtools view -H $aln_file`;
	print "Got ", scalar(@header_lines), " lines from $aln_file\n";
	warning "No header line found in $aln_file" if ( !scalar(@header_lines) );

	$self->header_lines( \@header_lines );

	return \@header_lines;
}

=head2 convert_bam_header_array_to_hash

  Arg [1]   : ReseqTrack::Tools::BamUtils::BamHeader
  Function  : Convert @header lines into hash for sequential output
  Returntype: hash ref to header contents 
  Exceptions: 
  Example   : $BH->convert_bam_header_array_to_hash();  
=cut

sub convert_bam_header_array_to_hash {
	my ($self) = @_;
	my %header_hash;

	my $header_lines = $self->header_lines;

	foreach my $line (@$header_lines) {
		my @aa = split /\t/, $line;
		push( @{ $header_hash{"$aa[0]"} }, $line );
	}

	$self->header( \%header_hash );

	return \%header_hash;
}


=head2 convert_bam_header_array_to_hash

  Arg [1]   : ReseqTrack::Tools::BamUtils::BamHeader
  Function  : Prints tmp file for use in reheading bam file
  Returntype: 
  Exceptions: No header line available
  Example   : $BH->print_header_to_file(); 
=cut

sub print_header_to_file {
	my ( $self, $outfile_name ) = @_;

	my @tag_order = ( '@HD', '@SQ', '@RG', '@PG', '@CO' );

	my $header = $self->header;

	if ( !defined $header ) {
		throw "No bam file header to ouput\n";
	}

	my $tmp_file =
	  !( defined($outfile_name) ) ? "$$\_tmp_header" : $outfile_name;

	open my $OUT, '>',
	  $tmp_file || throw "Failed to open $tmp_file for write\n";

	foreach my $tag (@tag_order) {

		next if ( !defined $$header{$tag} );

		if ( ref $$header{$tag} eq "ARRAY" ) {
			foreach my $line ( @{ $$header{$tag} } ) {
				print $OUT $line;
			}
		}
		else {
			print "Doing this\n";
			print $OUT $$header{$tag};
		}
	}

	close $OUT;

	$self->tmp_header_file($tmp_file);

	print "See $tmp_file for header information\n";

	return;
}


=head2 add_info_line_to_header

  Arg [1]   : ReseqTrack::Tools::BamUtils::BamHeader
  Arg [2]   : hash ref. Contains formatted lines to add to header
  Function  : Add header information to header hash
  Returntype: ref to hash containing header lines 
  Exceptions: 
  Example   : $BH->add_info_line_to_header(\%info);
=cut

sub add_info_line_to_header {
	
	my ( $self, $info ) = @_;
	my @RG_TAGS = qw(ID SM LB CN PI PL PU DS DT FO KS PG  );
	my @PG_TAGS = qw(ID PN CL PP VN  );

	my $header;

	if ( defined $self->header ) {
		$header = $self->header;
	}

	if ( ref($header) ne "HASH" || ref($info) ne "HASH" ) {
		throw("Need hash refs to function properly");
	}

	if ( defined $$info{'@RG'} ) {
		my @kv_pairs;
		push( @kv_pairs, '@RG' );
		foreach my $tag (@RG_TAGS) {
			if ( defined $$info{$tag} ) {
				my $pair = $tag . ":" . $$info{$tag};
				push( @kv_pairs, $pair );
			}
		}
		my $new_line = join( "\t", @kv_pairs );
		$new_line .= "\n";

		push( @{ $$header{'@RG'} }, $new_line );
	}

	if ( defined $$info{'@PG'} ) {
		my @kv_pairs;
		push( @kv_pairs, '@PG' );
		foreach my $tag (@PG_TAGS) {
			if ( defined $$info{$tag} ) {
				my $pair = $tag . ":" . $$info{$tag};
				push( @kv_pairs, $pair );
			}
		}
		my $line = join( "\t", @kv_pairs );
		$line .= "\n";

		push( @{ $$header{'@PG'} }, $line );
	}

	if ( defined $$info{'@CO'} ) {
		if ( ref( $$info{'@CO'} ) eq "ARRAY" ) {
			foreach my $line ( @{ $$info{'@CO'} } ) {
				my $new_line = "\@CO\t" . $line . "\n";
				push( @{ $$header{'@CO'} }, $new_line );
			}
		}
		else {

			my $new_line = "\@CO\t" . $$info{'@CO'};
			push( @{ $$header{'@CO'} }, $new_line );
		}
	}

	$self->header($header);

	return %$header;
}





sub text_file {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->{text_file} = $arg;
	}
	return $self->{text_file};
}

sub info_hash {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->{info_hash} = $arg;
	}
	return $self->{info_hash};
}

sub header_file {
	my ( $self, $arg ) = @_;
	if ($arg) {
		throw "Header file:$arg . Does not exist" if ( !-e $arg );
		$self->{header_file} = $arg;
	}

	return $self->{header_file};
}

sub header_lines {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->{header_lines} = $arg;
	}
	return $self->{header_lines};
}

sub header {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->{header} = $arg;
	}
	return $self->{header};

}

sub add_header {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{add_header} = $arg;
	}
	return $self->{add_header};
}

sub tmp_header_file {
	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		throw "Imaginary header file name: -$arg-" if ( !-e $arg );
		$self->{tmp_header_file} = $arg;
	}

	return $self->{tmp_header_file};
}

sub delete_old_bams {

	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{delete_old_bams} = $arg;
	}
	return $self->{delete_old_bams};
}
sub bam {

	my ( $self, $arg ) = @_;
	if ( defined $arg ) {
		$self->{bam} = $arg;
	}
	return $self->{bam};
}
1;
