=pod

=head1 NAME

ReseqTrack::Tools::RunVcfTools

=head1 SYNOPSIS

This is a class for running vcftools
It is a sub class of a ReseqTrack::Tools::RunProgram.

=head1 Example


=cut

package ReseqTrack::Tools::RunVcfTools;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::FileSystemUtils qw (check_executable);
use List::Util qw (first);
use Env qw( @PATH );

use base qw(ReseqTrack::Tools::RunProgram);


=head2 DEFAULT_OPTIONS

  Function  : Called by the RunProgram parent object in constructor
  Returntype: hashref
  Example   : my %options = %{&ReseqTrack::Tools:RunSamtools::DEFAULT_OPTIONS};

=cut

sub DEFAULT_OPTIONS { return {
        'remove_duplicates' => 0, # use the -d option when merging vcf
        };
}

sub CMD_MAPPINGS { return {
    'merge'    => \&run_merge,
    'concat'   => \&run_concat,
    };    
}

=head2 new

  Arg [-vcftools_dir]   :
      string, directory containing vcftools executables (not needed if vcftools is in $PATH)
  Arg [-bgzip]   :
      string, path of the bgzip executable (not needed if bgzip is in $PATH)
  Arg [-tabix]   :
      string, path of the tabix executable (not needed if tabix is in $PATH)
  Arg [-create_index]   :
      boolean, flag to create index files for all outputs
  Arg [-reference_index]   :
      string, path of a .fai file.  Needed when concatenating / sorting over multiple chromosomes
  + Arguments for ReseqTrack::Tools::RunProgram parent class

  Function  : Creates a new ReseqTrack::Tools::RunVcfTools object.
  Returntype: ReseqTrack::Tools::RunVcfTools
  Exceptions: 
  Example   : my $run_alignment = ReseqTrack::Tools::RunVcfTools->new(
                -input_files => ['/path/vcf1', '/path/vcf2'],
                -create_index => 1,
                -working_dir => '/path/to/dir/');

=cut

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my ( $vcftools_dir, $create_index, $reference_index, $bgzip, $tabix,
        )
    = rearrange( [
         qw( VCFTOOLS_DIR CREATE_INDEX REFERENCE_INDEX BGZIP TABIX
                ) ], @args);

  #setting defaults
  $self->vcftools_dir($vcftools_dir
              || $self->program
              || $ENV{VCFTOOLS}
              || first {/vcftools/} @PATH);

  $self->create_index( $create_index );
  $self->reference_index( $reference_index );
  $self->bgzip($bgzip || 'bgzip');
  $self->tabix($tabix || 'tabix');

  return $self;
}

sub run_merge {
    my $self = shift;

    throw("need two or more files to merge") if (@{$self->input_files} <2);

    my $output_vcf = $self->working_dir . '/' . $self->job_name . '.merged.vcf';
    $output_vcf =~ s{//}{/}g;

    my $executable = $self->vcftools_dir . '/vcf-merge';
    check_executable($executable);

    my @cmd_words = ($executable);
    push(@cmd_words, '-d') if ($self->options('remove_duplicates'));
    push(@cmd_words, @{$self->input_files});

    $output_vcf .= '.gz';
    push(@cmd_words, '|', $self->bgzip, '-c');
    push(@cmd_words, '>', $output_vcf);

    my $cmd = join(' ', @cmd_words);

    $self->output_files($output_vcf);
    $self->execute_command_line($cmd);

}

sub run_concat {
    my $self = shift;

    throw("need two or more files to concat") if (@{$self->input_files} <2);

    my $output_vcf = $self->working_dir . '/' . $self->job_name . '.concat.vcf';
    $output_vcf =~ s{//}{/}g;

    my $executable = $self->vcftools_dir . '/vcf-concat';
    check_executable($executable);

    my $ordered_inputs = $self->_order_input_files();

    my @cmd_words = ($executable);
    push(@cmd_words, @$ordered_inputs);

    $output_vcf .= '.gz';
    push(@cmd_words, '|', $self->bgzip, '-c');
    push(@cmd_words, '>', $output_vcf);

    my $cmd = join(' ', @cmd_words);

    $self->output_files($output_vcf);
    $self->execute_command_line($cmd);

}

sub run_tabix {
  my $self = shift;
  foreach my $vcf (@{$self->output_files}) {
    my $tbi_file = $vcf . '.tbi';
    my @cmd_words = ($self->tabix, '-p vcf', $vcf);
    my $cmd = join(' ', @cmd_words);
    $self->output_files($tbi_file);
    $self->execute_command_line($cmd);
  }
}




=head2 run_program

  Arg [1]   : ReseqTrack::Tools::RunVcfTools
  Arg [2]   : string, command, must be one of the following:
              merge, concat,
  Function  : uses vcftools to process the files in $self->input_files.
              Output is files are stored in $self->output_files
              This method is called by the RunProgram run method
  Returntype: 
  Exceptions: Throws if the command is not recognised.
  Example   : $self->run();

=cut

sub run_program {
    my ($self, $command) = @_;

    my $subs = CMD_MAPPINGS();

    throw("Did not recognise command $command") if (!defined $subs->{$command});

    print("Checking that tabix and bgzip are in your PATH (requirement of vcftools)\n");
    if (check_executable('tabix') && check_executable('bgzip')) {
      print "PATH verified\n";
    }

    &{$subs->{$command}}($self);

    if ($self->create_index) {
      $self->run_tabix;
    }

    return;
}

sub _order_input_files {
  my $self = shift;

  # Open each vcf file to find position of first record
  my %file_regions;
  foreach my $file(@{$self->input_files}) {
    my $fh;
    if ($file =~ /\.b?gz$/) {
      open $fh, "gunzip -c $file |" or throw("cannot open $file $!");
    }
    else {
      open $fh, '<', $file or throw("cannot open $file $!");
    }
    my ($chr, $pos);
    LINE:
    while(my $line = <$fh>) {
      next LINE if $line =~ /^#/;
      ($chr, $pos) = split(/\t/, $line);
      last LINE;
    }
    throw("did not find a valid record in $file") if (!$chr || !$pos);
    print STDERR "About to close a pipe.  This might give a 'broken pipe' warning but it can be ignored:\n";
    close $fh;
    $file_regions{$chr}{$file} = $pos;
  }

  # Open reference index file to find correct order of chromosomes
  my @ordered_chrs;
  if (keys %file_regions > 1) {
    throw("do not have a reference index file") if !$self->reference_index;
    open my $fh, '<', $self->reference_index or throw("cannot open ".$self->reference_index." $!");
    LINE:
    while (my $line = <$fh>) {
      next LINE if $line =~ /^#/;
      my ($chr) = split(/\t/, $line);
      next LINE if !$chr;
      push(@ordered_chrs, $chr);
    }
    close $fh;
    foreach my $chr (keys %file_regions) {
      throw("sequence $chr is not in the reference index file") if !grep {$chr eq $_} @ordered_chrs;
    }
  }
  else {
    @ordered_chrs = keys %file_regions;
  }

  my @ordered_inputs;
  foreach my $chr (@ordered_chrs) {
    push(@ordered_inputs, sort {$file_regions{$chr}{$a} <=> $file_regions{$chr}{$b}} keys %{$file_regions{$chr}});
  }
  return \@ordered_inputs;
}

sub get_valid_commands{
	my $subs = CMD_MAPPINGS();
	return keys %$subs;
}


sub vcftools_dir {
  my $self = shift;
  return $self->program(@_);
}

sub bgzip {
    my ($self, $arg) = @_;
    if ($arg) {
        $self->{'bgzip'} = $arg;
    }
    return $self->{'bgzip'};
}

sub tabix {
    my ($self, $arg) = @_;
    if ($arg) {
        $self->{'tabix'} = $arg;
    }
    return $self->{'tabix'};
}

sub output_vcf_files {
  my $self = shift;
  my @files = grep {/\.vcf(?:\.gz)$/} @{$self->output_files};
  return \@files;
}

sub output_tbi_files {
  my $self = shift;
  my @files = grep {/\.tbi$/} @{$self->output_files};
  return \@files;
}

sub create_index {
  my $self = shift;
  if (@_) {
    $self->{'create_index'} = (shift) ? 1 : 0;
  }
  return $self->{'create_index'};
}

sub reference_index {
  my ($self, $reference_index) = @_;
  if ($reference_index) {
    $self->{'reference_index'} = $reference_index;
  }
  return $self->{'reference_index'};
}


1;

