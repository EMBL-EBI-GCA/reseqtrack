package ReseqTrack::Tools::ERARunMetaInfoWhere;

use strict;
use warnings;
use vars qw(@ISA);

use Exporter;

use vars qw (@ISA  @EXPORT);

@ISA = qw(Exporter);

@EXPORT = qw(era_run_meta_info_where);


sub era_run_meta_info_where{
  my $study_ids = shift;
  if (!$study_ids) {
      $study_ids = [ qw(SRP001522 SRP004069 SRP004375 SRP004376 SRP001521
                SRP004068 SRP000805 SRP004055 SRP001515 SRP004062 SRP000547
                SRP004078 SRP000546 SRP004364 SRP001517 SRP004064 SRP001293
                SRP004365 SRP001523 SRP004070 SRP000808 SRP004058 SRP001294
                SRP004060 SRP001519 SRP004066 SRP000806 SRP004056 SRP001518
                SRP004065 SRP001514 SRP004061 SRP000544 SRP004076 SRP004374
                SRP004373 SRP001516 SRP004063 SRP000543 SRP004075 SRP001520
                SRP004067 SRP000807 SRP004057 SRP004369 SRP004370 SRP000803
                SRP004054 SRP001524 SRP004071 SRP004077 SRP004059 SRP001525
                SRP004072 SRP004372 SRP004371 SRP000540 SRP004073 SRP000542
                SRP004074 SRP000031 SRP000032 SRP000033) ];
  }

  my $string = "era.g1k_sequence_index_all.study_id in ('"
            . join("', '", @$study_ids) . "')";
  return $string;
}

1;
