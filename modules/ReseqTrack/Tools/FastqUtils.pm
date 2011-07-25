package ReseqTrack::Tools::FastqUtils;
use strict;
use warnings;
use IO::Compress::Gzip qw(gzip $GzipError) ;
use File::Basename;
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);

use vars qw (@ISA  @EXPORT);

@ISA    = qw(Exporter);
@EXPORT = qw(sample  estimate_average_read_length);

=pod

=head1 NAME

FastqUtils

=head1 SYNOPSIS

A module for fastq files manipulation.
Extracted from Sanger code and modified for purposes of 1000 Genomes
Project DCC

=cut

=head1 AUTHORS

Petr Danecek, I<pd3@sanger.ac.uk>
Richard Smith <smithre.ebi.ac.uk>

=cut

=head1 METHODS



=head2 sample_files

=cut


=head2 sample

        Arg [1]    : input fastq file (may or may not be gzipped)
        Arg [2]    : output fastq file (may or may not be gzipped)
        Arg [3]    : sample size (bases)
        Arg [4]    : Hash reference with the list of sequence indexes to extract (indexes from 0). [optional]
                     If the hash is empty, all reads will be selected (the fastq file will be symlinked). 
                     If the hash is not present, the routine will select the reads and return the hash
                     of indexes, which can be passed to subsequent calls of sample.
        Description: Selects sequences evenly from the fastq file 
        Returntype : hash reference with the list of sequence indexes to extract (indexes from 0)

=cut

sub sample
  {
    my ($in_file, $out_file, $sample_size, $idxs, $do_sample) = @_;

    my $fh_in;

    # With the list of indexes, sampling is easy - just repeat what has been done for the other lane.
    #   Otherwise, see how big is the file and find out how dense should be the sampling.
    if ( !defined $idxs ) {


      print "Indices not defined.Collecting\n";

      if ( $in_file =~ /\.gz$/i ) {
	open($fh_in, "zcat $in_file | ") or die("failed zcat $in_file: $!");
      } else {
	open($fh_in, '<' ,$in_file ) or die ("failed open $in_file: $!");
      }
      # Read only the first four lines:
      #   @IL17_2470:1:1:11:194/2
      #   GNAAGGCGGNNGNANNNNAGCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
      #   +
      #   +%+99666+%%+%+%%%%+0+%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      #   
      my $line_id   = <$fh_in>;
      my $line_seq  = <$fh_in>;
      my $line_sep  = <$fh_in>;
      my $line_qual = <$fh_in>;
      close $fh_in;

      # How many sequences must be read for the sample
      my $nseqs_to_read = $sample_size / (length($line_seq)-1);

      #  Estimate (assume that all fields have the same length).
      #
      my  $nseqs_total;
        
      my $file_size ;

      if ( $in_file =~ /\.gz$/i ) {
	$file_size  = uncompressed_gz_size($in_file)
      } else {
	$file_size = (-s $in_file);
      }

      $nseqs_total = $file_size / ( length($line_id) + length($line_seq) + length($line_sep) + length($line_qual)  );
      print "Total reads $nseqs_total\n";
 

      my $nseqs_to_skip = int($nseqs_total / $nseqs_to_read);

      # If equal to 1, all sequences will be included. 
      # For the other pair, if the returned $idxs is empty, do the same.
      if ( $nseqs_to_skip <= 0 ) {
	$nseqs_to_skip = 1;
      }
      if ( $nseqs_to_skip == 1 ) {
	$idxs = {};
      } else {
	for (my $i=$nseqs_to_skip; $i<$nseqs_total; $i+=$nseqs_to_skip ) {
	  $$idxs{$i} = 1;
	}
      }
    

      if ( ! (scalar keys %$idxs) ) {
	print "Base in fastq < subsample limit. Not sampling\n";
	return;
      } else {
	print "Have indices " , scalar (keys %$idxs) , "\n";
	return ($idxs);
      }
    }

    my $fastq_compressed;
    my $z;
    my $fh_out;

    if ( $in_file =~ /\.gz$/i ) {
      open($fh_in, "zcat $in_file | ") or error("zcat $in_file: $!");
      $fastq_compressed = 1;
      $z = new IO::Compress::Gzip  $out_file 
	or die "IO::Compress::Gzip failed: $GzipError\n";
    

    } else {
      open($fh_in, '<', $in_file) or error("$in_file: $!");
      open $fh_out,'>', $out_file || die "Failed to open $out_file\n";
      $fastq_compressed = 0;
    }


    print "Creating $out_file\n";
 
    my $iblock = 0;
    while (my $line_id=<$fh_in>) {
      my $line_seq  = <$fh_in>;
      my $line_sep  = <$fh_in>;
      my $line_qual = <$fh_in>;

      if ( exists($$idxs{$iblock}) ) {
	if ( !$line_seq || !$line_sep || !$line_qual ) {
	  error("Sanity check failed.\n");
	}

	$z->print ($line_id, $line_seq, $line_sep, $line_qual)     if ( $fastq_compressed) ;
	print $fh_out ($line_id, $line_seq, $line_sep, $line_qual) if (!$fastq_compressed );
      }

      $iblock++;
    }
    close $fh_in;
  
    return;
  }


=head2 uncompressed_gz_size

    Arg[1]          : The gzipped file name
    Description     : Determines what is the size of the uncompressed file.
    Returntype      : The size of the uncompressed file in bytes.

=cut

sub uncompressed_gz_size
  {
    my ($file) = @_;

    # From some reason, this does not work ("Use of uninitialized value in <HANDLE>")
    #
    #   my ($writer, $reader, $err);
    #   open3($writer, $reader, $err, "gzip -l $file");
    #   my @output = <$reader>;
    #   my @errors = <$err>;

    # Uh, gzip -l is not reliable. Occasionally it gives wrong uncompressed size, as it
    #   happened e.g. for 3290_5_2.fastq.gz. In that case, the ratio was -2208.2%
    #
    #    my ($writer, $reader);
    #    open3($writer, $reader, \*ERR, "gzip -l $file");
    #    my @output = <$reader>;
    #    my @errors = <ERR>;
    #
    #    if ( scalar @errors ) { error("gzip -l $file:\n", @errors) }
    #
    #    # Expected output:
    #    #   compressed        uncompressed  ratio uncompressed_name
    #    #   1772198718          4147914907  57.3% mouse-2470_1_2.fastq
    #    if ( scalar @output != 2 || !($output[1]=~/^\s+\d+\s+(\d+)/ ) ) { error("Uh, expected something else:\n", @output) }

    my @output = `zcat $file | wc -c`;
    if ( scalar @output != 1 ) {
      die("Expected something else from zcat $file | wc -c, got: ", join('',@output), "\n");
    }
    if ( !($output[0] =~ /^(\d+)$/) ) {
      die("Uh, expeceted something else from zcat $file | wc -c, not \"$output[0]\"\n");
    }
    return $1;
  }


sub estimate_average_read_length{
  my ($in_file, $platform) = @_;
 
  my $fh_in;
  my $ctr;
 
  if ( $in_file =~ /\.gz$/i ) {
    open($fh_in, "zcat $in_file | ") or die("failed zcat $in_file: $!");
  } else {
    open($fh_in, '<' ,$in_file ) or die ("failed open $in_file: $!");
  }

  my ($read_count,$base_count,$ave_read_length);
 
  my ($line1, $line2, $line3, $line4);
  my $read_length = 0;

  while (my $line1 = <FH>) {
  
    my $line2 = <FH>;
    my $line3 = <FH>;
    my $line4 = <FH>;
 

    $line2 =~ s/\s+$//;
    chomp $line2;
 
    $ctr++;

    $read_count++;
    my $length += length($line2);

    #  print $length,"\n";

    if ( $platform eq 'ABI_SOLID') {
      $length -= 1;
    }

    $base_count += $length;

    if ( $ctr == 4000){
      close ($fh_in);
      $ave_read_length = int ($base_count/$read_count);
      print "Average read length of 1000 reads = $ave_read_length\n";
      return $ave_read_length;
    }
  }

  throw "Failed to get average read length\n";
}



1;

