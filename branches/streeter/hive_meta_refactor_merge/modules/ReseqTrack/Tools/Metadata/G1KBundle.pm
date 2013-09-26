package ReseqTrack::Tools::Metadata::G1KBundle;

use strict;
use warnings;

use base qw(ReseqTrack::Tools::Metadata::AddInBundle);

sub addins_to_load {
	return qw(
		ReseqTrack::Tools::Metadata::PopulationRulesManipulator
		ReseqTrack::Tools::Metadata::G1KManipulator
		ReseqTrack::Tools::Metadata::G1KCheckRunStatus
		ReseqTrack::Tools::Metadata::G1KCheckSampleNameInFile
		ReseqTrack::Tools::Metadata::G1KMetaDataClashCheck
	);
}
1;