use strict;

use Test::More;
use ReseqTrack::Tools::RunProgram;

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
        # 'gubbins' => 1, this is set by the defaults
        'monkey' => 1,
};

my $rp = ReseqTrack::Tools::RunProgram->new(
	-options => $opts,
);
$opts->{gubbins} = 1;

is('HASH',ref($rp->options),'return all type');
is_deeply($rp->options, $opts, "construct+default");
is('foo',$rp->options('whatsit'),'retrieve_one');

my $extra_opts = {
	extra_opt => 'bar',
};
$rp->options($extra_opts);
$opts->{extra_opt} = 'bar';
is_deeply($rp->options, $opts, "merge");

$rp->options('mismatch_penalty','bar');
is('bar',$rp->options('mismatch_penalty'),'set_one');



done_testing();