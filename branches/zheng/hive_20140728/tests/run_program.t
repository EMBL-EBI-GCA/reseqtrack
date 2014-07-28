use strict;

use Test::More;
use ReseqTrack::Tools::RunProgram;

# Monkey patch run program to implement a method normally left to the implementing class
local *ReseqTrack::Tools::RunProgram::DEFAULT_OPTIONS = sub {
	return {
		'thingy' => 15,
    'whatsit' => 'what',
    'gubbins' => 1, 
    'monkey' => 1,
	}	
};

my $opts = {
        'thingy' => 10,
        'whatsit' => 'foo',
        'monkey' => 1,
};

my $rp = ReseqTrack::Tools::RunProgram->new(
	-options => $opts,
);
$opts->{gubbins} = 1;

is('HASH',ref($rp->options),'return all type');
# RP now has a merge of the defaults and those given in the constructor. The constructor values override the defaults
is_deeply($rp->options, $opts, "construct+default"); 
is('foo',$rp->options('whatsit'),'retrieve_one');

my $extra_opts = {
	extra_opt => 'bar',
};
$rp->options($extra_opts);
$opts->{extra_opt} = 'bar';
is_deeply($rp->options, $opts, "merge"); # additional values will be merged in without losing the old values

$rp->options('mismatch_penalty','bar'); 
is('bar',$rp->options('mismatch_penalty'),'set_one');

done_testing();