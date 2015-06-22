package ReseqTrack::Hive::Process::SampleFactory;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::FileSystemUtils qw(check_executable check_file_exists);

sub param_defaults {
  return {  
    vcf => undef,  
    tabix => undef,
  };
}

sub run {
    my $self = shift @_;
    my $vcf = $self->param_required('vcf');

    my $tabix = $self->param('tabix');

	throw("Please bgzip and index the genotype VCF file") unless ($vcf =~ /.vcf.gz$/);
    
    check_file_exists($vcf);
    #print "input file to sample factory is $vcf\n";
     
    my $index = $vcf . ".tbi";
    check_file_exists($index);
    
    check_executable($tabix);
    
    open (my $IN, "$tabix $vcf -H |") || throw("Cannot open VCF file $vcf by tabix $!");
    
    while ( <$IN> ){
        chomp;
        next unless $_ =~ /^\#CHROM/;
        my @column_headers = split(/\t/, $_);
        foreach my $i ( 9..$#column_headers ) {  ### fan_index doesn't have to start from 0
             $self->prepare_factory_output_id({'sample' 	=> $column_headers[$i],
             								   'fan_index'	=> $i-9 });
        }
    }        
}


1;

