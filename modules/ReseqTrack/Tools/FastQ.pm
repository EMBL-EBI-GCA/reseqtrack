package ReseqTrack::Tools::FastQ;
use strict;
use warnings;
use IO::Compress::Gzip qw(gzip $GzipError) ;

use ReseqTrack::Tools::Utils qw(CMD error cleanup exit_cleanup uncompressed_gz_size shadow_hierarchy relative_symlink basename log_msg );

use vars qw (@ISA  @EXPORT);

@ISA    = qw(Exporter);
@EXPORT = qw(sample);

=pod

=head1 NAME

FastQ

=head1 SYNOPSIS

A module for fastq files manipulation.

=cut

=head1 AUTHORS

Petr Danecek, I<pd3@sanger.ac.uk>

=cut

=head1 METHODS

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
    my ($in_file, $out_file, $sample_size, $idxs) = @_;

    my $fh_in;

    # With the list of indexes, sampling is easy - just repeat what has been done for the other lane.
    #   Otherwise, see how big is the file and find out how dense should be the sampling.
    if ( !defined $idxs )
    {
      print "indices not defined\n";
        open($fh_in, "zcat $in_file | ") or error("zcat $in_file: $!");

        # Read only the first four lines:
        #   @IL17_2470:1:1:11:194/2
        #   GNAAGGCGGNNGNANNNNAGCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
        #   +
        #   +%+99666+%%+%+%%%%+0+%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        #   
        my $line_id   = <$fh_in>;
        my $line_seq  = <$fh_in>;
        my $line_sep  = <$fh_in>;
        my $line_qual = <$fh_in>;
        close $fh_in;

        # How many sequences must be read for the sample
        my $nseqs_to_read = $sample_size / (length($line_seq)-1);

        # How many blocks (sequences) are in the file. First look in the fastqcheck file.
        #   If it does not exist, try to estimate (assumes that all fields have the same length).
        #
        my ($fqcheck_data, $nseqs_total);
        eval { $fqcheck_data = parse_fastqcheck($in_file . '.fastqcheck'); };
        if ( $@ )
        {
            my $file_size = ( $in_file =~ /\.gz$/i ) ? uncompressed_gz_size($in_file) : -s $in_file;
            $nseqs_total = $file_size / ( length($line_id) + length($line_seq) + length($line_sep) + length($line_qual)  );
	    print "1. Total reads $nseqs_total\n";
        }
        else
        {
            $nseqs_total = $$fqcheck_data{nseqs};
	    print "2. Total reads $nseqs_total\n";
        }

        my $nseqs_to_skip = int($nseqs_total / $nseqs_to_read);

        # If equal to 1, all sequences will be included. Just symlink the original
        #   file. For the other pair, if the returned $idxs is empty, do the same.
        if ( $nseqs_to_skip <= 0 ) { $nseqs_to_skip = 1; }
        if ( $nseqs_to_skip == 1 )
        {
            $idxs = {};
        }
        else
        {
            for (my $i=$nseqs_to_skip; $i<$nseqs_total; $i+=$nseqs_to_skip )
            {
                $$idxs{$i} = 1;
            }
        }
    

      if ( !scalar keys %$idxs )
	{
	  # Utils::relative_symlink($in_file,$out_file) unless -e $out_file;
	  return;
	}
      else{
	print "Have indices returning\n";
	return ($idxs);
      }
    }


    if ( $in_file =~ /\.gz$/i ) { open($fh_in, "zcat $in_file | ") or error("zcat $in_file: $!"); }
    else { open($fh_in, '<', $in_file) or error("$in_file: $!"); }

    my $tmp_out_file = "$out_file";
    print "Creating $tmp_out_file\n";
    $Utils::remove_on_exit->{$tmp_out_file} = 1;    # Clean this file if we get killed
    
     my $z = new IO::Compress::Gzip  $tmp_out_file 
      or die "IO::Compress::Gzip failed: $GzipError\n";
       
    my $iblock = 0;
    while (my $line_id=<$fh_in>)
    {
        my $line_seq  = <$fh_in>;
        my $line_sep  = <$fh_in>;
        my $line_qual = <$fh_in>;

        if ( exists($$idxs{$iblock}) )
        {
            if ( !$line_seq || !$line_sep || !$line_qual ) { error("Sanity check failed.\n") }
            $z->print ($line_id, $line_seq, $line_sep, $line_qual);
        }

        $iblock++;
    }
    close $fh_in;
  
    print "done\n";
    return;
}



=head2 parse_fastqcheck

        Arg [1]    : input fastqcheck file
        Returntype : hash reference with the data, for use with Graphs::plot_stats

=cut

sub parse_fastqcheck
{
    my ($in_file) = @_;

    my @out = ();

    open(my $fh,'<',$in_file) or error("$in_file: $!");

    # 23891022 sequences, 1815717672 total length, 76.00 average, 76 max
    my $line  = <$fh>; 
    if ( !($line=~/^(\d+)\s+sequences/) ) { error("Could not parse $in_file: $line"); }
    my $nseqs = $1;

    <$fh>; 

    # A    C    G    T    N    0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  2
    $line = <$fh>;
    my @header = split(/\s+/, $line);
    my @xvals  = splice(@header,6,-1);

    while ($line=<$fh>)
    {
        # Note: Only first 50 positions are listed in the fastqcheck file. The line Total contains
        #   histogram for the whole read. If the reads have 50 bases, the line Total would be the
        #   averages of each column.
        #
        # Total    28.4 21.4 21.5 28.6  0.1    0   0   ..   2   1   2   2  -72-136 0 25.0
        # base  1  29.6 19.6 21.7 28.9  0.3    0   0   ..   0   3   1   4  175 639 0 27.3
        #
        my @items = split(/\s+/, $line);

        my $title = $items[0];
        if ( $title=~/^base/ ) 
        { 
            # If we are here, we are parsing one of the 'base X' lines.
            $title .= ' ' . shift(@items); 
        } 

        my @data  = splice(@items,6,-1);
        if ( scalar @xvals != scalar @data ) 
        { 
            error("Could not parse the file, different data size?? ", scalar @xvals," ",scalar @data,"\nThe file .. $in_file\nThe line .. $line");
        }
        push @out, {'title'=>$title, 'xvals'=>\@xvals, 'yvals'=>\@data };
    }
    close($fh);

    return {'data'=>\@out, 'nseqs'=>$nseqs };
}
1;

