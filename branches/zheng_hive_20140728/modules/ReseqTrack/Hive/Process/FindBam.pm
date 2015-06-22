package ReseqTrack::Hive::Process::FindBam;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists);

sub param_defaults {
  return {
    sample => undef,
    bam_type => undef,
  };
}

sub run {
    my $self = shift @_;
    
	my $db_params = $self->param_required('reseqtrack_db');
	my $sample = $self->param_required('sample');
	
	my $bam_type = $self->param_required('bam_type');
	throw("BAM type needs to be either BAM or EXOME_BAM") if ($bam_type ne "BAM" && $bam_type ne "EXOME_BAM");
	
	my $db = ReseqTrack::DBSQL::DBAdaptor->new(%{$db_params});
	
	my $f_objs = $db->get_FileAdaptor->fetch_all_like_name($sample);
	
	$db->dbc->disconnect_when_inactive(1);
	
	my $bam;
	my %bams;
	foreach my $fo ( @$f_objs ) {
		#print "fetched file: " . $fo->name . "\n";
		next if ($fo->type ne $bam_type);
		next unless ($fo->name =~ /\.mapped/);
		next if ($fo->name =~ /high_coverage/);
		$bam = $fo->name;
		$bams{$bam} = 1;
	} 

	check_file_exists($bam);
	
	throw("More than one BAM of type $bam_type found for $sample") if ( keys %bams > 1);
	$self->output_param('bam', $bam);  
}

1;