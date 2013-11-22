=pod

=head1 NAME

ReseqTrack::Tools::RunMergePhasedVcf

=head1 SYNOPSIS

This is a class for running mergevcf.jar from Beagle utilities
It is a sub class of a ReseqTrack::Tools::RunProgram.

=head1 Example


=cut

package ReseqTrack::Tools::RunMergePhasedVcf;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename qw(fileparse);
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists);

use base qw(ReseqTrack::Tools::RunProgram);

=head2 new

  Arg [-chrom]   :
      string, chrom name
  Arg [-java_exe]   :
      string, Java exe path
  + Arguments for ReseqTrack::Tools::RunProgram parent class

  Function  : Creates a new ReseqTrack::Tools::RunMergePhasedVcf object.
  Returntype: ReseqTrack::Tools::RunMergePhasedVcf
  Exceptions: 
  Example   : my $run_alignment = ReseqTrack::Tools::RunMergePhasedVcf->new(
                -input_files => ['/path/vcf1', '/path/vcf2'],
                -java_exe => "/usr/bin/java",
                -program => "mergevcf.jar",
                -chrom => 11, ## CHROM name
                -create_index => 1,
                -bgzip => “/path/bgzip“,
                -tabix => “/path/tabix“,
                -working_dir => '/path/to/dir/');

=cut

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my ( $java_exe, $chrom, $bgzip, $tabix, $create_index
        )
    = rearrange( [
         qw( JAVA_EXE CHROM BGZIP TABIX CREATE_INDEX
                ) ], @args);


  $self->java_exe($java_exe  || '/usr/bin/java' );
  $self->bgzip($bgzip  || '/nfs/1000g-work/G1K/work/bin/tabix/bgzip' );
  $self->tabix($tabix  || '/nfs/1000g-work/G1K/work/bin/tabix/tabix' );
  $self->chrom($chrom);
  $self->create_index( $create_index );

  return $self;
}

=head2 run_program

  Arg [1]   : ReseqTrack::Tools::RunPhase::RunMergePhasedVcf
  Function  : uses Beagle mergevcf.jar to combine overlaped phased vcf files in $self->input_files.
              Output file is stored in $self->output_files
              This method is called by the RunProgram run method
  Returntype: 
  Exceptions: 
  Example   : $self->run();

=cut

sub run_program {
    my ($self) = @_;
    
    check_executable($self->java_exe);
    check_executable($self->bgzip);
    
    throw("chrom name required") if(!$self->chrom);
      
    my @cmd_words = ($self->java_exe, $self->jvm_args, '-jar');
    push(@cmd_words,$self->program);
    
    push(@cmd_words,$self->chrom);
    
    foreach my $input (@{$self->input_files}) {
      push(@cmd_words, $input);
    }
    
    my $output_vcf = $self->working_dir .'/'. $self->job_name;
    $output_vcf =~ s{//}{/};
    $output_vcf .= '.vcf.gz';
    $self->output_files($output_vcf);
    
    push(@cmd_words, '|', $self->bgzip ,'-c','>',$output_vcf);
    my $cmd = join(' ', @cmd_words);
    
    $self->execute_command_line($cmd);
    
    if ($self->create_index) {
      $self->run_tabix;
    }
    
    return;
}

sub run_tabix {
  my $self = shift;
  foreach my $vcf (@{$self->output_files}) {
    my $tbi_file = $vcf . '.tbi';
    my @cmd_words = ($self->tabix, '-p vcf', $vcf) unless(-e $tbi_file);
    my $cmd = join(' ', @cmd_words);
    $self->output_files($tbi_file);
    $self->execute_command_line($cmd);
  }
}

sub output_vcf_files {
    my $self = shift;
    my @files = grep { /\.vcf$|\.vcf\.gz$/ } @{ $self->output_files };
    return \@files;
}

sub java_exe {
    my ($self, $arg) = @_;
    if ($arg) {
        $self->{'java_exe'} = $arg;
    }
    return $self->{'java_exe'};
}

sub chrom {
    my ($self, $arg) = @_;
    if ($arg) {
        $self->{'chrom'} = $arg;
    }
    return $self->{'chrom'};
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

sub create_index {
  my $self = shift;
  if (@_) {
    $self->{'create_index'} = (shift) ? 1 : 0;
  }
  return $self->{'create_index'};
}

1;