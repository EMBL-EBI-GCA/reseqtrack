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
        'read_trimming' => 15,
        'mismatch_penalty' => '',
        'gap_open_penalty' => '',
        'gap_extension_penalty' => '',
        'max_gap_opens' => '',
        'max_gap_extensions' => '',
        'threads' => 1,
        'load_fm_index' => 1,
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

    my @cmd_words = ($self->program, 'sampe');
    push(@cmd_words, '-P') if $self->options('load_fm_index');

    if ($self->read_group_fields->{'ID'}) {
      my $rg_string = q('@RG\tID:) . $self->read_group_fields->{'ID'};
      RG:
      while (my ($tag, $value) = each %{$self->read_group_fields}) {
        next RG if ($tag eq 'ID');
        next RG if (!$value);
        $rg_string .= '\t' . $tag . ':' . $value;
      }
      $rg_string .= q(');
      push(@cmd_words, '-r', $rg_string);
    }
    push(@cmd_words, $self->reference, $mate1_sai, $mate2_sai);
    push(@cmd_words, $self->get_fastq_cmd_string('mate1'));
    push(@cmd_words, $self->get_fastq_cmd_string('mate2'));
    push(@cmd_words, $output_sam);

    my $sampe_cmd = join(' ', @cmd_words);

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

    my @cmd_words = ($self->program, 'samse');

    if ($self->read_group_fields->{'ID'}) {
      my $rg_string = q('@RG\tID:) . $self->read_group_fields->{'ID'};
      RG:
      while (my ($tag, $value) = each %{$self->read_group_fields}) {
        next RG if ($tag eq 'ID');
        next RG if (!$value);
        $rg_string .= '\t' . $tag . ':' . $value;
      }
      $rg_string .= q(');
      push(@cmd_words, '-r', $rg_string);
    }
    push(@cmd_words, $self->reference, $sai_file);
    push(@cmd_words, $self->get_fastq_cmd_string('frag'));
    push(@cmd_words, $output_sam);

    my $samse_cmd = join(' ', @cmd_words);
      
    $self->output_files($output_sam);
    $self->execute_command_line($samse_cmd);

    return $output_sam;
}

sub run_aln_mode {
    my ( $self, $fastq_type ) = @_;

    my $options = $self->aln_options;

    my $output_file = $self->working_dir . "/" . $self->job_name;
    $output_file .= ".$fastq_type.sai";

    my @cmd_words = ($self->program, 'aln');

    push(@cmd_words, '-q', $self->options('read_trimming'))
            if ($self->options('read_trimming'));
    push(@cmd_words, '-M', $self->options('mismatch_penalty'))
            if ($self->options('mismatch_penalty'));
    push(@cmd_words, '-O', $self->options('gap_open_penalty'))
            if ($self->options('gap_open_penalty'));
    push(@cmd_words, '-E', $self->options('gap_extension_penalty'))
            if ($self->options('gap_extension_penalty'));
    push(@cmd_words, '-o', $self->options('max_gap_opens'))
            if ($self->options('max_gap_opens'));
    push(@cmd_words, '-e', $self->options('max_gap_extensions'))
            if ($self->options('max_gap_extensions'));

    push(@cmd_words, '-t', $self->options('threads') || 1);
    push(@cmd_words, $self->reference);
    push(@cmd_words, $self->get_fastq_cmd_string($fastq_type));
    push(@cmd_words, $output_file);

    my $aln_command = join(' ', @cmd_words);

    $self->created_files($output_file);
    $self->execute_command_line($aln_command);

    return $output_file;
}

1;

