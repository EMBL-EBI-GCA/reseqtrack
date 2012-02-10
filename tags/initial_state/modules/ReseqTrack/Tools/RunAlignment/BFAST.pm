package ReseqTrack::Tools::RunAlignment::BFAST;

use strict;
use warnings;

use Data::Dumper;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);

use base qw(ReseqTrack::Tools::RunAlignment);

sub new {

	my ( $class, @args ) = @_;
	my $self = $class->SUPER::new(@args);

	my ( $preprocess_exe,)
	  = rearrange( [
			qw( PREPROCESS_EXE )
		], @args);

	$self->ali_options('-A 1 -K 8 -M 384 -n 4 -Q 25000 ');

	$self->preprocess_exe($preprocess_exe);



	return $self;

}
#################################################################


sub run {
    my ($self) = @_;
    my $exit;

    $self->change_dir();

    my $tmp_bmf = $self->run_match();
    my $tmp_baf = $self->run_localalign($tmp_bmf);

    my $sam = $self->run_postprocess($tmp_baf);
    $self->sam_files($sam);

    $self->run_samtools(1);

    $self->files_to_delete( [$tmp_bmf, $tmp_baf] );

    $self->get_bfast_tmp_files;

    return;
}




############

sub run_postprocess {
	my ($self, $tmp_baf) = @_;
	my $cmd_line;

        throw( "no baf file\n" )
            if (! $tmp_baf );

	my $tmp_sam =  $self->working_dir . '/' . $$. ".sam";
        $tmp_sam =~ s{//}{/};

	$cmd_line = $self->program() . " postprocess ";
	$cmd_line .= " -f " . $self->reference() . " ";
	$cmd_line .= " -i  $tmp_baf ";
	$cmd_line .= " -n 4 -Q 1000 -t -U  > $tmp_sam";

        $self->execute_command_line($cmd_line);

	return $tmp_sam;
}

sub run_localalign {
	my ($self, $tmp_bmf) = @_;

        throw( "no bmf file\n" )
            if (! $tmp_bmf );

	my $cmd_line;

	my $tmp_baf =  $self->working_dir . '/' . $$. ".baf";
        $tmp_baf =~ s{//}{/};

	$cmd_line = $self->program() . " localalign ";
	$cmd_line .= " -f " . $self->reference() . " ";
	$cmd_line .= " -m  $tmp_bmf ";
	$cmd_line .= " -A 1 -o 20 -n 4 -t > $tmp_baf ";
        
        $self->execute_command_line($cmd_line);
	
    return $tmp_baf;
}

sub run_match {
	my $self = shift;
	my $cmd_line;

	my $tmp_bmf =  $self->working_dir . '/' . $$. ".bmf";
        $tmp_bmf =~ s{//}{/};

	$cmd_line = $self->preprocess_exe() . " ";

        foreach my $file (@{$self->input_files}) {
            $cmd_line .= $file . " ";
        }

	$cmd_line .= " \| ";

	$cmd_line .= $self->program() . " match ";
	$cmd_line .= " -f " . $self->reference() . " ";
	$cmd_line .= $self->ali_options() . " ";
	$cmd_line .= ">  $tmp_bmf";

        $self->execute_command_line($cmd_line);

    return $tmp_bmf;
}

sub ali_options {
	my ( $self, $arg ) = @_;

	if ($arg) {
		$self->{ali_options} = $arg;
	}
	return $self->{ali_options};
}

sub preprocess_exe {
	my ( $self, $arg ) = @_;

	if ($arg) {
		$self->{preprocess_exe} = $arg;
	}
	return $self->{preprocess_exe};
}


sub get_bfast_tmp_files {
  my ( $self ) = @_;

  my $dir = $self->working_dir;
  my $tmp_files_found = 0; 
  print $dir,"\n";

  if ( !(-d $dir) ||  !( -e $dir) ||  !( $dir =~ /^\// ) ){
	warning "Bad directory info for bfast tmp file search.Skipping\n";
	return;
      }


  opendir ( DIR, $dir ) || die "Error in opening dir $dir\n";

  my @files = readdir(DIR);

  foreach my $file (@files) {
    if ( $file =~  /\.bfast\.tmp/ ){
      print $file,"\n";
      my $del_file = $dir . '/'. $file;
      $del_file =~ s/\/\//\//;
      $self->files_to_delete ($del_file);
      $tmp_files_found ++;
    }
  }
 
  print "Number of temporary bfast files found = $tmp_files_found\n";
  
  return;
}


1;

