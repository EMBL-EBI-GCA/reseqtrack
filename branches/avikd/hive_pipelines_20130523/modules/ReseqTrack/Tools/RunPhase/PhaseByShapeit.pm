package ReseqTrack::Tools::RunPhase::PhaseByShapeit;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists check_executable);

use base qw(ReseqTrack::Tools::RunPhase);
=pod

=head1 NAME

ReseqTrack::Tools::RunPhase::PhaseByShapeit

=head1 SYNOPSIS

This is a class for running Shapeit Phasing 
It is a sub class of a ReseqTrack::Tools::RunPhase

It generates and run command to phase variants from multiple diploid individuals.
If running with a reference panel, it can check for alignment between panel and 
variant call sets. It can exclude SNPs which are missing in ref panel and having
incompatible allele types (strand issues) OR only SNPs with strand issue.
Phasing can be done with or without alignment checking, and with or without 
reference panel


Case 1: CHECK with panel, REMOVE all incompatible SNPs, PHASE with reference panel
    
    Hints: run with -reference_config
        
Case 2: CHECK with panel, REMOVE SNPs with strand issue, PHASE without reference panel

    Hints: use -exclude_only_strand and -phase_without_ref and -reference_config

Case 3: No check, PHASE without reference panel  

    Hints: use -no_check and -phase_without_ref or do not provide a -reference_config

Note:
	Reference file format: chrom <tab> map <tab> haps <tab> legend <tab> samples

=cut

=head2 new

  Arg [-parameters]   :
      hashref, command line options to use with Shapeit2
            
  + Arguments for ReseqTrack::Tools::RunProgram parent class

  Function  : Creates a new ReseqTrack::Tools::RunPhase::PhaseByShapeit object.
  Returntype: ReseqTrack::Tools::RunPhase::PhaseByShapeit
  Exceptions: 
  Example   : my $varPhase_byShapeit = ReseqTrack::Tools::RunPhase::PhaseByShapeit->new(
                -program        =>"/path/shapeit" ,
                -input_files            => ‘/path/to/gen‘,
                -input_samples             => ‘/path/to/samples‘,
                -working_dir             => '/path/to/dir/',  ## this is working dir and output dir
                -reference_config               => '/path/to/ref_config', ###format:chrom\tmap\thaps\tlegend\tsamples
                -options              => {	phase => "--states 100 --window 0.5", }, 
                -chrom                   => '1',
		-phase_without_ref	 =>1, ## no ref will be used for phasing
		-exclude_only_strand	 =>1, ## only strand flip snps will be removed
                -no_check                => 1  ## only perform phase, not strand check with ref
                );

=cut

# This is to set default parameters
sub DEFAULT_OPTIONS { return {
       phase=>'--window 0.5 --states 100',
        };
}

sub new {
	my ( $class, @args ) = @_;
	my $self = $class->SUPER::new(@args);
	
    my ( $input_samples, $no_check, $phase_without_ref, $exclude_only_strand )= rearrange([ qw( INPUT_SAMPLES NO_CHECK PHASE_WITHOUT_REF EXCLUDE_ONLY_STRAND )], @args); ### ?
    
    
    $self->input_samples($input_samples);
    $self->no_check($no_check);
    $self->phase_without_ref($phase_without_ref);
    $self->exclude_only_strand($exclude_only_strand);
    
    unless($self->reference_config)
    {
        $self->no_check(1);
        $self->phase_without_ref(1);
    }
    
    if($self->exclude_only_strand)
    {
        $self->phase_without_ref(1);
    }

    if($self->no_check)
    {
        $self->phase_without_ref(1);
    }

    return $self;
}

=head2 run_program

  Arg [1]   : ReseqTrack::Tools::RunPhase::PhaseByShapeit
  Function  : uses Shapeit2 to phase variants in $self->input_files.
              Output is files are stored in $self->output_files
              This method is called by the RunProgram run method
  Returntype: 
  Exceptions: 
  Example   : $self->run();

=cut

sub run_program {
	my $self = shift;
    
    check_executable($self->program);
    check_file_exists($self->reference_config) if $self->reference_config; # ref is optional, only required for check
  
    check_file_exists($self->input_samples);

    throw "Can't accept multiple input files" if(scalar @{$self->input_files} >1);   ### ?
    
    if($self->no_check)
    {
        $self->run_shapeit();
    }
    else
    {

        if($self->reference_config)
        {
            $self->run_check();

	   if($self->exclude_only_strand){
            	$self->make_exclude_list();
	   }

            $self->run_shapeit();
        }
        else
        {
            throw "Reference required to preform alignment check";
        }
    }        
    return;
}

sub run_check {
    my $self = shift;


    my @cmd_words = ($self->program, "-check");
    
    push(@cmd_words, "-G", ${$self->input_files}[0], $self->input_samples);
    
    my $output_check = $self->working_dir .'/'. $self->job_name.".alignments";
    $output_check =~ s{//}{/};
    push(@cmd_words, "--output-log", $output_check);
    
    my $exclude_snp_list=$output_check.".snp.strand.exclude";
    
    $self->exclude_snp_list($exclude_snp_list);
        
    push(@cmd_words, $self->options->{'check'});
    
    push(@cmd_words,"--chrX") if($self->chrom eq 'X' || $self->chrom eq "chrX");

    $self->get_ref;
	   
 
    push(@cmd_words, "-M", $self->ref_map);
   
    my $ref_files_line=$self->ref_hap." ".$self->ref_legend." ".$self->ref_samples;
    
    push(@cmd_words, "-R", $ref_files_line);
   
    my $cmd = join(' ', @cmd_words);
    
    $self->execute_command_line ($cmd,\&allow_exit_code_1);
    
    return;
}

sub allow_exit_code_1 {
  my ($self, $exit, $command_line) = @_;
  throw("unexpected exit code $exit $command_line") if $exit != 1;
}

sub make_exclude_list {
    my $self = shift;
    
    my $check_snp_list=$self->working_dir .'/'. $self->job_name.".alignments.snp.strand";
    
    my $exclude_snp_list=$self->working_dir .'/'. $self->job_name.".alignments.snp.strand_only.exclude";
    
    my @cmd_words = ("cat",$check_snp_list,"|grep strand|cut -f2 >",$exclude_snp_list);
    
    $self->exclude_snp_list($exclude_snp_list);
    
    my $cmd = join(' ', @cmd_words);
    
    $self->execute_command_line ($cmd);
    
    return;
}

sub run_shapeit {
    my $self = shift;
    
    my @cmd_words = ($self->program);
    push(@cmd_words, "-G", ${$self->input_files}[0]." ".$self->input_samples);
    
    my $output_phase = $self->working_dir .'/'. $self->job_name;
    $output_phase =~ s{//}{/};
    push(@cmd_words, "--output-max", $output_phase);
    
    my $output_haps=$output_phase.".haps";
    my $output_samples=$output_phase.".samples";
    $self->output_files($output_haps);
    
    push(@cmd_words,"--chrX") if($self->chrom eq 'X' || $self->chrom eq "chrX");

    push(@cmd_words, $self->options->{'phase'});

    unless($self->phase_without_ref){

    	if($self->reference_config)
    	{
        	$self->get_ref;
        	push(@cmd_words, "-M", $self->ref_map);
	
		my $ref_files_line=$self->ref_hap." ".$self->ref_legend." ".$self->ref_samples;
		push(@cmd_words, "-R", $ref_files_line);
    	}
    }

    if($self->exclude_snp_list)
    {
        push(@cmd_words,"--exclude-snp",$self->exclude_snp_list);
    }
    
       
    my $cmd = join(' ', @cmd_words);
    
    $self->execute_command_line ($cmd);
    
    return;
}

sub get_ref {
    my ( $self, $arg ) = @_;
    my $reference_config = $self->reference_config; ###format:chrom\tmap\thaps\tlegend\tsamples
    my $chrom= $self->chrom;
   
    if($chrom eq 'X' || $chrom eq 'chrX')
    {
     	throw("No support for ChrX yet");
    }
    else
    {
    	open(FH,$reference_config)||throw("can't read $reference_config\n");
    	while(my $ref_line=<FH>)
    	{
        	chomp($ref_line);
        	next if($ref_line=~ /^#/);
        
        	my @field=split(/\t/,$ref_line);
        
        	if($field[0]==$chrom)
        	{
           		$self->ref_map($field[1]);
           		$self->ref_hap($field[2]);
           		$self->ref_legend($field[3]);
           		$self->ref_samples($field[4]);
        	}
	}
    }
    return $self->{get_ref};
}

sub ref_map {
   my ( $self, $arg ) = @_;

  if ($arg) {
    $self->{ref_map} = $arg;
  }
  return $self->{ref_map};
}

sub ref_hap {
   my ( $self, $arg ) = @_;

  if ($arg) {
    $self->{ref_hap} = $arg;
  }
  return $self->{ref_hap};
}

sub ref_legend {
   my ( $self, $arg ) = @_;

  if ($arg) {
    $self->{ref_legend} = $arg;
  }
  return $self->{ref_legend};
}

sub ref_samples {
   my ( $self, $arg ) = @_;

  if ($arg) {
    $self->{ref_samples} = $arg;
  }
  return $self->{ref_samples};
}

sub exclude_snp_list {
  my ( $self, $arg ) = @_;

  if ($arg) {
    $self->{exclude_snp_list} = $arg;
  }
  return $self->{exclude_snp_list};
}

sub no_check {
    my $self = shift;
    if (@_) {
        $self->{'no_check'} = (shift) ? 1 : 0;
    }
    return $self->{'no_check'};
}

sub phase_without_ref {
    my $self = shift;
    if (@_) {
        $self->{'phase_without_ref'} = (shift) ? 1 : 0;
    }
    return $self->{'phase_without_ref'};
}

sub exclude_only_strand {
	my $self = shift;
    if (@_) {
        $self->{'exclude_only_strand'} = (shift) ? 1 : 0;
    }
    return $self->{'exclude_only_strand'};
}

sub input_samples {
  my ( $self, $arg ) = @_;

  if ($arg) {
    $self->{input_samples} = $arg;
  }
  return $self->{input_samples};
}

1;
