package ReseqTrack::Tools::myTIME;

use strict;
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::Argument qw(rearrange);
use Time::Local;

sub get_time {
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime();
	$mon++; 
	$year = $year + 1900;
	my $time = $year . "-" . $mon . "-" . $mday . "." . $hour . "." . $min;
	my $month = $year . "-" . $mon ;
	my $day = $year . "-" . $mon . "-" . $mday;
	return ($time, $month, $day);
}

1;
