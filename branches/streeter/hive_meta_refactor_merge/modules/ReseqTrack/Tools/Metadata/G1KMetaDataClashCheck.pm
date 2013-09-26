package ReseqTrack::Tools::Metadata::G1KMetaDataClashCheck;

use strict;
use warnings;
use base qw(ReseqTrack::Tools::Metadata::BaseMetaDataClashCheck);

1;

sub attribute_2_attribute_whitelist {
  return {
    'GSSR_ID'        => { 'experiment' => { 'run' => 1 } },
    'LSID'           => { 'experiment' => { 'run' => 1 } },
    'PROJECT'        => { 'experiment' => { 'run' => 1 } },
    'ROOT_SAMPLE_ID' => { 'experiment' => { 'run' => 1 } },
    'WORK_REQUEST'   => { 'experiment' => { 'run' => 1 } },
    'SAMPLE_ID'      => { 'experiment' => { 'run' => 1 } },
  };
}

sub attribute_2_column_whitelist {
  return {
    'study' => { 'type' => { 'sample' => { 'TYPE' => 1 } } },
    'experiment' =>
      { 'center_name' => { 'experiment' => { 'CENTER_NAME' => 1 } } },
    'run' => {
      'center_name' => { 'experiment' => { 'CENTER_NAME' => 1 } },
      'sample_id' => { 
        'experiment' => { 'SAMPLE_ID' => 1 }, 
        'run' => { 'SAMPLE_ID' => 1 } 
        }
    },
    'sample' => {
      'center_name' => { 'experiment' => { 'CENTER_NAME' => 1 } },
      'sample_id' =>
        { 'experiment' => { 'SAMPLE_ID' => 1 }, 'run' => { 'SAMPLE_ID' => 1 } }
    },
  };
}

