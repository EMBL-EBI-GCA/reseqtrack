package ReseqTrack::Tools::RunAlignment::BFAST;

use strict;
use warnings;

use Data::Dumper;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use Env qw( @PATH );
use List::Util qw (first);

use base qw(ReseqTrack::Tools::RunAlignment);

sub new {

	my ( $class, @args ) = @_;
	my $self = $class->SUPER::new(@args);

	my ( $preprocess_exe,)
	  = rearrange( [
			qw( PREPROCESS_EXE )
		], @args);

	#setting defaults
        if (!$self->program) {
          if ($ENV{BFAST}) {
            $self->program($ENV{BFAST} . '/bfast');
          }
          else {
            $self->program(first {-x $_} map {"$_/bfast"} @PATH);
          }
        }

	$self->ali_options('-A 1 -K 8 -M 384 -n 4 -Q 25000 ');

	$self->preprocess_exe($preprocess_exe);



	return $self;

}
#################################################################


sub run_alignment {
    my ($self) = @_;

    my $tmp_dir = $self->working_dir()
          .'/'.$self->job_name.'.'.$$.'.tmp/';
    check_file_does_not_exist($tmp_dir);
    make_directory($tmp_dir);
    $self->created_files($tmp_dir);
    $self->tmp_dir($tmp_dir);

    my $tmp_bmf = $self->run_match();
    my $tmp_baf = $self->run_localalign($tmp_bmf);

    $self->run_postprocess($tmp_baf);

    return;
}




############

sub run_postprocess {
	my ($self, $tmp_baf) = @_;
	my $cmd_line;

        throw( "no baf file\n" )
            if (! $tmp_baf );

	my $output_sam =  $self->working_dir . '/' . $self->job_name. ".sam";
        $output_sam =~ s{//}{/};

	$cmd_line = $self->program() . " postprocess ";
	$cmd_line .= " -f " . $self->reference() . " ";
	$cmd_line .= " -i  $tmp_baf ";
        $cmd_line .= " -T " . $self->tmp_dir;
	$cmd_line .= " -n 4 -Q 1000 -t -U  > $output_sam";

        $self->sam_files($output_sam);
        $self->execute_command_line($cmd_line);

}

sub run_localalign {
	my ($self, $tmp_bmf) = @_;

        throw( "no bmf file\n" )
            if (! $tmp_bmf );

	my $cmd_line;

	my $tmp_baf =  $self->working_dir . '/' . $self->job_name . ".baf";
        $tmp_baf =~ s{//}{/};

	$cmd_line = $self->program() . " localalign ";
	$cmd_line .= " -f " . $self->reference() . " ";
	$cmd_line .= " -m  $tmp_bmf ";
	$cmd_line .= " -A 1 -o 20 -n 4 -t > $tmp_baf ";
        
        $self->created_files($tmp_baf);
        $self->execute_command_line($cmd_line);
	
    return $tmp_baf;
}

sub run_match {
	my $self = shift;
	my $cmd_line;

	my $tmp_bmf =  $self->working_dir . '/' . $self->job_name. ".bmf";
        $tmp_bmf =~ s{//}{/};

	$cmd_line = $self->preprocess_exe() . " ";

        if ( $self->mate1 ) {
            $cmd_line .= $self->get_fastq_cmd_string('mate1') . " ";
        }
        if ( $self->mate2 ) {
            $cmd_line .= $self->get_fastq_cmd_string('mate2') . " ";
        }
        if ( $self->fragment_file ) {
            $cmd_line .= $self->get_fastq_cmd_string('frag') . " ";
        }

	$cmd_line .= " \| ";

	$cmd_line .= $self->program() . " match ";
	$cmd_line .= " -f " . $self->reference() . " ";
	$cmd_line .= $self->ali_options() . " ";
        $cmd_line .= " -T " . $self->tmp_dir;
	$cmd_line .= ">  $tmp_bmf";

        $self->created_files($tmp_bmf);
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

sub tmp_dir {
	my ( $self, $arg ) = @_;

	if ($arg) {
		$self->{tmp_dir} = $arg;
	}
	return $self->{tmp_dir};
}

1;

