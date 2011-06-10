package ReseqTrack::Tools::QC::QCUtils;

use strict;
use warnings;
use Exporter;
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use File::Basename;

use vars qw (@ISA  @EXPORT_OK);

@ISA       = qw(Exporter);
@EXPORT_OK = qw(
  get_params
);



sub get_params {

    my $file  = shift;
    my $input = shift;

   
    throw ("Could not open $file") if ( !-e $file);

    open my $IN, '<', $file || die "No config file found";

    while (<$IN>) {
        chomp $_;
        next if ( !$_ );
        my @aa = split /=/;
	$aa[1] =~ s/\s+//g;
        if (defined  $$input{ $aa[0] }){
          print "$aa[0] already set. Skipping cfg entry\n";
          next;
        }
        $$input{ $aa[0] } = $aa[1] ;
    }
    close $IN;
    
    return ( \%$input );
}





1;
