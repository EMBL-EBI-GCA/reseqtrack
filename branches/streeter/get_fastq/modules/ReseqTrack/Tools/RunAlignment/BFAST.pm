package ReseqTrack::Tools::RunAlignment::BFAST;

use strict;
use warnings;

use Data::Dumper;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);

use base qw(ReseqTrack::Tools::RunAlignment);

=pod

=head1 NAME

ReseqTrack::Tools::RunAlignment::BFAST

=head1 SYNOPSIS

class for running BFAST. Child class of ReseqTrack::Tools::RunAlignment

=head1 Example

my $bfast = ReseqTrack::Tools::RunAlignment::BFAST(
                      -input => '/path/to/file'
                      -reference => '/path/to/reference',
                      -options => {'threads' => 4},
                      -paired_length => 3000,
                      -read_group_fields => {'ID' => 1, 'LB' => 'my_lib'},
                      -first_read => 1000,
                      -last_read => 2000,
                      -preprocess_exe => '/path/to/1kg_2_bfast_fastq',
                      );
$bfast->run;
my $output_file_list = $bfast->output_files;

=cut


sub DEFAULT_OPTIONS { return {
        'threads' => 1,
        'colour_space' => 1,
        'offset' => 20,
        'postprocess_algorithm' => 3,
        'pairing' => undef,
        };
}

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
            $self->program('bfast');
          }
        }

	$self->preprocess_exe($preprocess_exe);

	return $self;

}
#################################################################


sub run_alignment {
    my ($self) = @_;

    my $tmp_bmf = $self->run_match();
    my $tmp_baf = $self->run_localalign($tmp_bmf);

    $self->run_postprocess($tmp_baf);

    return;
}




############

sub run_postprocess {
    my ($self, $tmp_baf) = @_;

    throw( "no baf file\n" )
        if (! $tmp_baf );

    my $output_file =  $self->working_dir . '/' . $self->job_name;
    $output_file .= ($self->output_format eq 'BAM' ? '.bam' : '.sam');
    $output_file =~ s{//}{/};

    my @cmd_words;
    push(@cmd_words, $self->program, 'postprocess');
    push(@cmd_words, '-f', $self->reference);
    push(@cmd_words, '-i', $tmp_baf);
    push(@cmd_words, '-a', $self->options('postprocess_algorithm')) if (defined $self->options('postprocess_algorithm'));
    push(@cmd_words, '-Y', $self->options('pairing')) if (defined $self->options('pairing'));
    push(@cmd_words, '-v', $self->paired_length) if ($self->paired_length);

    if ($self->read_group_fields->{'ID'}) {
      my $rg_string = "\@RG\tID:" . $self->read_group_fields->{'ID'};
      RG:
      while (my ($tag, $value) = each %{$self->read_group_fields}) {
        next RG if ($tag eq 'ID');
        next RG if (!$value);
        $rg_string .= "\t$tag:$value";
      }
      my $RG_file = $self->working_dir . '/' . $self->job_name. ".RG_string";
      $self->created_files($RG_file);
      open (my $fh, '>', $RG_file) or throw ("could not open file $RG_file");
      print $fh $rg_string;
      close $fh;
      push(@cmd_words, '-r', $RG_file);
    }

    if ($self->output_format eq 'BAM') {
      push(@cmd_words, '|', $self->samtools, 'view -bS -');
    }
    push(@cmd_words, '>', $output_file);

    my $cmd_line = join(' ', @cmd_words);

    $self->output_files($output_file);
    $self->execute_command_line($cmd_line);

}

sub run_localalign {
    my ($self, $tmp_bmf) = @_;

    throw( "no bmf file\n" )
        if (! $tmp_bmf );

    my $tmp_baf =  $self->working_dir . '/' . $self->job_name . ".baf";
    $tmp_baf =~ s{//}{/};

    my @cmd_words;
    push(@cmd_words, $self->program, 'localalign');
    push(@cmd_words, '-f', $self->reference);
    push(@cmd_words, '-n', $self->options('threads') || 1);
    push(@cmd_words, '-A', $self->options('colour_space') ? 1 : 0);
    push(@cmd_words, '-o', $self->options('offset')) if ($self->options('offset'));
    push(@cmd_words, '-m', $tmp_bmf);

    push(@cmd_words, '>', $tmp_baf);

    my $cmd_line = join(' ', @cmd_words);
    
    $self->created_files($tmp_baf);
    $self->execute_command_line($cmd_line);
	
    return $tmp_baf;
}

sub run_match {
    my $self = shift;

    my $tmp_bmf =  $self->working_dir . '/' . $self->job_name. ".bmf";
    $tmp_bmf =~ s{//}{/};

    my @cmd_words = ("bash -c '");
    push(@cmd_words, $self->preprocess_exe());
    push(@cmd_words, $self->get_fastq_cmd_string('mate1')) if ($self->mate1_file);
    push(@cmd_words, $self->get_fastq_cmd_string('mate2')) if ($self->mate2_file);
    push(@cmd_words, $self->get_fastq_cmd_string('frag')) if ($self->fragment_file);

    push(@cmd_words, '|', $self->program, 'match');
    push(@cmd_words, '-f', $self->reference);
    push(@cmd_words, '-n', $self->options('threads') || 1);
    push(@cmd_words, '-A', $self->options('colour_space') ? 1 : 0);
    push(@cmd_words, '-T', $self->get_temp_dir);

    push(@cmd_words, '>', $tmp_bmf);
    push(@cmd_words, "'");

    my $cmd_line = join(' ', @cmd_words);

    $self->created_files($tmp_bmf);
    $self->execute_command_line($cmd_line);

    return $tmp_bmf;
}

sub preprocess_exe {
	my ( $self, $arg ) = @_;

	if ($arg) {
		$self->{preprocess_exe} = $arg;
	}
	return $self->{preprocess_exe};
}

1;

