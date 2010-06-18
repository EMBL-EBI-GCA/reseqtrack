package ReseqTrack::Tools::AlignmentIndexUtils;

#package AlignmentIndexUtils;

use strict;
use warnings;
use Exporter;
use File::Copy;
use File::Basename;
use File::Find ();
use ReseqTrack::Tools::Exception;
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::Tools::GeneralUtils;

use vars qw (@ISA  @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(
             get_bas_hash_array
	    );

sub get_bas_hash_array{
  my ($file, $trim) = @_;
  my $line_count = 0;
  #print "bas file is $file\n"; 
  open(FH, "gunzip -c $file |") or throw("IndexUtils:get_index_hash failed to open ".$file." $!");
  my %hash;
  while(<FH>){
    chomp;
    my $line = $_;
    next if(/filename/i);
    my @values = split /\t/, $line;
    my $name = $values[0];
    $name = basename($values[0]) if($trim);
    $name =~ s/\s+$//;
    $name =~ s/^\s+//;

    push @{$hash{$name}}, $line;

  }
  close(FH);
  return (\%hash);
}

1;
