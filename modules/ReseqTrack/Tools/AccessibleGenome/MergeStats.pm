
package AccessibleGenome::MergeStats;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::FileSystemUtils qw(check_directory_exists check_executable check_file_exists);
use File::Basename qw(basename);

sub param_defaults {
  return {
    bgzip => 'bgzip',
    tabix => 'tabix',
    run_tabix => 0,
  };
}


sub run {
    my $self = shift @_;
    $self->param_required('stats');
    my $bgzip = $self->param('bgzip') || $self->param_defaults->bgzip;
    my $tabix = $self->param('tabix') || $self->param_defaults->tabix;
    my $run_tabix = $self->param('run_tabix');

    my $stats = $self->param_as_array('stats');

    $self->dbc->disconnect_when_inactive(1);

    foreach my $stats_path (grep {defined $_} @$stats) {
    	print "file path is $stats_path\n"; 
    	#throw();
		next if ((basename $stats_path) =~ /\.HLA|decoy/); ##CHECKME
		check_file_exists($stats_path);
	}

	my $output_dir = $self->output_dir;
	my $job_name = $self->job_name;
    my $output_file = "$output_dir/$job_name.stats.gz";
    check_directory_exists($output_dir);
    check_executable($bgzip);
    if ($run_tabix) {
      check_executable($tabix);
    }
    open my $OUT, "| $bgzip -c > $output_file" or throw ("cannot bgzip to $output_file: $!");
    my $first_stats = 1;

    STATS:
    foreach my $i (0..$#{$stats}) {
      my $stats_path = $stats->[$i];
      next STATS if ! defined $stats_path;
      my $IN;
      if ($stats_path =~ /\.b?gz(?:ip)?$/){
        if ($i == 0) {
	        open $IN, "gzip -cd $stats_path |" or throw("cannot open $stats_path: $!");
        }
        else {
        	open $IN, "gzip -cd $stats_path | grep -v ZeroMapQual |" or throw("cannot open $stats_path: $!");
        }
      }
      else {
        if ($i == 0 ) {
	        open $IN, '<', $stats_path or throw("cannot open $stats_path: $!");
        }
        else {
	        open $IN, "grep -v ZeroMapQual $stats_path |" or throw("cannot open $stats_path: $!");
        }
      }
	
      if ($first_stats) {
        $first_stats = 0;
      }
      else {
        <$IN>;
      }
	
      LINE:
      while (my $line = <$IN>) {
        print $OUT $line;
      }
      close $IN;
    }
    close $OUT;

    if ($run_tabix) {
      my $cmd = "$tabix -s1 -b2 -e2 $output_file";
      print $cmd, "\n";
      system($cmd) ==0 or throw("tabix failed $!");
    }

    $self->dbc->disconnect_when_inactive(0);

    $self->output_param('stats' => $output_file);
    if ($run_tabix) {
      $self->output_param('tbi' => "$output_file.tbi");
    }
}


1;

