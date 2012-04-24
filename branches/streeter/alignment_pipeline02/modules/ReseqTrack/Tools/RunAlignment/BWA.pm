package ReseqTrack::Tools::RunAlignment::BWA;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::RunAlignment;
use Env qw( @PATH );
use List::Util qw (first);

@ISA = qw(ReseqTrack::Tools::RunAlignment);

=pod

=head1 NAME

=head1 SYNOPSIS

=head1 Example

=head2 new

  Arg [1]   : ReseqTrack:Tools::RunAlignment::BWA
  Arg [2]   : string, options for the single ended bwa mod
  Arg [3]   : string, options for paired end bwa mod
  Function  : create a BWA object, defaults program name to bwa and
  defaults to using ReseqTrack::Tools::RunAlignment::options for sampe options 
  if sampe options aren't defined and options are
  Returntype: ReseqTrack::Tools::RunAlignment::BWA
  Exceptions: n/a
  Example   : my $bwa = ReseqTrack::Tools::RunAlignment::BWA(
                      -program => "bwa",
                      -input => '/path/to/file'
                      -reference => '/path/to/reference',
                      -samtools => '/path/to/samtools',
                      );

=cut

sub DEFAULT_OPTIONS { return [
        'aln_q' => 15,
        'threads' => '',
        ];
}

sub new {
	my ( $class, @args ) = @_;
	my $self = $class->SUPER::new(@args);

	#setting defaults
        if (!$self->program) {
          if ($ENV{BWA}) {
            $self->program($ENV{BWA} . '/bwa');
          }
          else {
            $self->program(first {-x $_} map {"$_/bwa"} @PATH);
          }
        }

	return $self;
}

=head2 accessor methods

  Arg [1]   : ReseqTrack::Tools::RunAlignment::BWA
  Arg [2]   : string, options string
  Function  : samse_options and sampe_options are specific bwa commandline
  options and these are the accessor methods for those variables
  Returntype: string
  Exceptions: n/a
  Example   : my $options = $self->samse_options;

=cut


sub run_alignment {
    my ($self) = @_;
    
    if ( $self->fragment_file ) {
        $self->run_samse_alignment();
    }
    
    if ( $self->mate1_file && $self->mate2_file ) {
        $self->run_sampe_alignment();
    }

    return;
}

sub run_sampe_alignment {
	my ($self) = @_;
	my $mate1_sai =
	  $self->run_aln_mode( "mate1" );
	my $mate2_sai =
	  $self->run_aln_mode( "mate2" );
        my $output_sam = $self->working_dir() . '/'
            . $self->job_name
            . '_pe.sam';
        $output_sam =~ s{//}{/};

	if ( $self->paired_length() ) {
		my $paired_length = " -a " . $self->paired_length;
		$self->sampe_options($paired_length);
	}

	my $sampe_cmd = $self->program . " sampe ";
	$sampe_cmd .= $self->sampe_options . " " if ($self->sampe_options);

        if ($self->read_group_fields->{'ID'}) {
          my $rg_string = q('@RG\tID:) . $self->read_group_fields->{'ID'};
          RG:
          while (my ($tag, $value) = each %{$self->read_group_fields}) {
            next RG if ($tag eq 'ID');
            next RG if (!$value);
            $rg_string .= '\t' . $tag . ':' . $value;
          }
          $rg_string .= q(');
          $sampe_cmd .= " -r $rg_string ";
        }

	$sampe_cmd .= $self->reference . " "
	  . $mate1_sai . " "
	  . $mate2_sai . " "
	  . $self->get_fastq_cmd_string('mate1') . " "
	  . $self->get_fastq_cmd_string('mate2') . " > "
          . $output_sam;

        $self->output_files($output_sam);
        $self->execute_command_line($sampe_cmd);

	return $output_sam;
}

sub run_samse_alignment {
	my ($self) = @_;

	my $sai_file =
	  $self->run_aln_mode( "frag" );

        my $output_sam = $self->working_dir() . '/'
            . $self->job_name
            . '_se.sam';
        $output_sam =~ s{//}{/};

	my $samse_cmd = $self->program . " samse ";
	$samse_cmd .= $self->samse_options . " " if ($self->samse_options);

        if ($self->read_group_fields->{'ID'}) {
          my $rg_string = q('@RG\tID:) . $self->read_group_fields->{'ID'};
          RG:
          while (my ($tag, $value) = each %{$self->read_group_fields}) {
            next RG if ($tag eq 'ID');
            next RG if (!$value);
            $rg_string .= '\t' . $tag . ':' . $value;
          }
          $rg_string .= q(');
          $samse_cmd .= " -r $rg_string ";
        }

	$samse_cmd .= $self->reference . " "
	  . $sai_file . " "
	  . $self->get_fastq_cmd_string('frag') . " > "
	  . $output_sam;
	  
        $self->output_files($output_sam);
        $self->execute_command_line($samse_cmd);

	return $output_sam;
}

sub run_aln_mode {
	my ( $self, $fastq_type ) = @_;

	my $bwa_log_file =  $self->working_dir. '/'. $self->job_name .".$$.log";
	$bwa_log_file =~ s/\/\//\//;
	my $do_bwa_log_file = "2>> ". $bwa_log_file;

	my $options = $self->aln_options;

	my $output_file = $self->working_dir . "/" . $self->job_name;
	$output_file .= ".$fastq_type.sai";

	my $aln_command;
	$aln_command .= $self->program;
	$aln_command .= " aln ";
	$aln_command .= $options . "  " if $options;
        $aln_command .= '-t ' . $self->options('threads') . ' ' if ($self->options('threads');
	$aln_command .= $self->reference . " ";
	$aln_command .= get_fastq_cmd_string($fastq_type);
	$aln_command .= " $do_bwa_log_file > ";
	$aln_command .= $output_file;

        $self->created_files($output_file);
        $self->created_files($bwa_log_file);
        $self->execute_command_line($aln_command);

	return $output_file;
}

sub samse_options {
    my ( $self, $arg ) = @_;
    unless ( $self->{samse_options} ) {
        $self->{samse_options} = '';
    }
    if ($arg) {
        $self->{samse_options} = $arg;
    }
    return $self->{samse_options};
}

sub sampe_options {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{sampe_options} = $arg;
    }
    return $self->{sampe_options};
}

sub aln_options {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{aln_options} = $arg;
    }
    return $self->{aln_options};
}

sub threads {
    my $self = shift;
    if (@_) {
        $self->{threads} = shift;
    }
    return $self->{threads};
}


1;

