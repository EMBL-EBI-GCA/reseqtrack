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

Case 4: For ChrX, use -chrX ,-run_from and -run_to,  to phase PAR1, PAR2 or nonPAR regions separately, or phase without ref

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
                -chrom                   => '1', # or ‘X_PAR1’, required for fetching reference files
		-phase_without_ref	 =>1, ## no ref will be used for phasing
		-exclude_strand_flip	 =>1, ## only strand flip snps will be removed
                -no_check                => 1  ## only perform phase, not strand check with ref
                -chrX =>1, # required for ChrX
                -region_start => 150118, # required for ChrX
                -region_end => 2695340, # required for ChrX
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
	
  my ( $input_samples, $no_check, $phase_without_ref, $exclude_strand_flip, $chrX )= rearrange([ qw( INPUT_SAMPLES NO_CHECK PHASE_WITHOUT_REF EXCLUDE_STRAND_FLIP CHRX)], @args); ### ?
    
  $self->input_samples($input_samples);
  $self->no_check($no_check);
  $self->phase_without_ref($phase_without_ref);
  $self->exclude_strand_flip($exclude_strand_flip);
  $self->chrX($chrX);
  
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
    
    check_file_exists($self->reference_config) if $self->reference_config; # ref is optional, only required for check
  
    throw "Can't accept multiple input files" if(scalar @{$self->input_files} >1);   
    throw "expecting single samples file" if(scalar @{$self->input_samples} >1);
    throw("samples input file required") if(!$self->input_samples);
    
    foreach my $sample(@{$self->input_samples}){
      check_file_exists($sample);
    }

    
    
    if($self->chrX)
    {
      throw("Coordinate range required for ChrX") unless($self->region_start && $self->region_end);
    }
    
   ###Check options
    unless($self->reference_config)
    {
      $self->no_check(1);
      $self->phase_without_ref(1);
    }
    
    if($self->exclude_strand_flip)
    {
      $self->phase_without_ref(1);
    }

    if($self->no_check)
    {
      $self->phase_without_ref(1);
    }
    
    ###run
    if(!$self->reference_config && $self->no_check)
    {
      $self->run_shapeit();
    }
    elsif($self->reference_config && !$self->no_check)
    {
      $self->run_check();
      
      $self->make_exclude_list() if($self->exclude_strand_flip);
      
      $self->run_shapeit();
      
    } 
    elsif(!$self->reference_config && !$self->no_check)
    {
      throw "Reference required to preform alignment check";
    }
           
    return;
}

=head2 run_check

  Arg [1]   : ReseqTrack::Tools::RunPhase::PhaseByShapeit
  Function  :uses Shapeit2 to check variants in $self->input_files.
             Output is files are stored in $self->created_files and 
             $self->exclude_snp_list
            
              
  Returntype: 
  Exceptions: 
  Example   : $self->run_check()

=cut

sub run_check {
    my $self = shift;


    my @cmd_words = ($self->program, "-check");
    
    push(@cmd_words, "-G", ${$self->input_files}[0], ${$self->input_samples}[0]);
    
    my $output_check = $self->working_dir .'/'. $self->job_name.".alignments";
    $output_check =~ s{//}{/};
    push(@cmd_words, "--output-log", $output_check);
    
    my $check_snp_list=$output_check.".snp.strand";
    my $exclude_snp_list=$output_check.".snp.strand.exclude";
    
    $self->exclude_snp_list($exclude_snp_list);
        
    push(@cmd_words, $self->options->{'check'}) if($self->options->{'check'});
    
    push(@cmd_words,"--chrX") if($self->chrX eq 1);  ### check cannot be performed without ref
    
    $self->get_ref;
	   
    push(@cmd_words, "-M", $self->ref_map);
   
    my $ref_files_line=$self->ref_hap." ".$self->ref_legend." ".$self->ref_samples;
    
    push(@cmd_words, "-R", $ref_files_line);
    
    push(@cmd_words, "--input-from", $self->region_start ) if($self->region_start);
    push(@cmd_words, "--input-to", $self->region_end ) if($self->region_end);
    push(@cmd_words, "--output-from", $self->region_start ) if($self->region_start);
    push(@cmd_words, "--output-to", $self->region_end ) if($self->region_end);
   
    my $cmd = join(' ', @cmd_words);
    
    $self->created_files($exclude_snp_list);
    $self->created_files($check_snp_list);
    $self->created_files($output_check.".log");
    
    $self->execute_command_line ($cmd,\&allow_exit_code_1);
    
    return;
}

sub allow_exit_code_1 {
  my ($self, $exit, $command_line) = @_;
  throw("unexpected exit code $exit $command_line") if $exit != 1;
}

=head2 make_exclude_list

  Arg [1]   : ReseqTrack::Tools::RunPhase::PhaseByShapeit
  Function  : Create a list of variants with strand flip error 
              Output is files are stored in $self->exclude_snp_list
              
  Returntype: 
  Exceptions: 
  Example   : $self->make_exclude_list()

=cut

sub make_exclude_list {
    my $self = shift;
    
    my $check_snp_list=$self->working_dir .'/'. $self->job_name.".alignments.snp.strand";
    
    #if( -e $check_snp_list ) {
    
        my $exclude_snp_list=$self->exclude_snp_list.".strand_only";
    
        my @cmd_words = ("cat",$check_snp_list,"|grep strand|cut -f2 >", $exclude_snp_list);
    
        $self->exclude_snp_list($exclude_snp_list) ;
    
        my $cmd = join(' ', @cmd_words);
    
        $self->created_files($exclude_snp_list);
    
        $self->execute_command_line ($cmd);
   # }
    
    return;
}


=head2 run_shapeit

  Arg [1]   : ReseqTrack::Tools::RunPhase::PhaseByShapeit
  Function  : uses Shapeit2 to check variants in $self->input_files.
              Output is files are stored $self->output_files
            
              
  Returntype: 
  Exceptions: 
  Example   :  $self->run_shapeit();

=cut

sub run_shapeit {
    my $self = shift;
    
    my @cmd_words = ($self->program);
    push(@cmd_words, "-G", ${$self->input_files}[0]." ".${$self->input_samples}[0]);
    
    my $output_phase = $self->working_dir .'/'. $self->job_name;
    $output_phase =~ s{//}{/};
    push(@cmd_words, "--output-max", $output_phase);
    
    my $output_haps=$output_phase.".haps";
    my $output_samples=$output_phase.".sample";
    $self->output_files($output_haps);
    $self->output_files($output_samples);
    
    push(@cmd_words,"--chrX") if($self->chrX eq 1);

    push(@cmd_words, $self->options->{'phase'});

    unless($self->phase_without_ref)
    {
      if($self->reference_config)
      {   	   
        $self->get_ref;
        push(@cmd_words, "-M", $self->ref_map);
	  
	    my $ref_files_line=$self->ref_hap." ".$self->ref_legend." ".$self->ref_samples;
	    push(@cmd_words, "-R", $ref_files_line);
	  }
    }

    if($self->exclude_snp_list) {
        push(@cmd_words,"--exclude-snp",$self->exclude_snp_list) if( -e $self->exclude_snp_list);
    }
        
    push(@cmd_words, "--input-from",$self->region_start ) if($self->region_start);
    push(@cmd_words, "--input-to",$self->region_end ) if($self->region_end);
    push(@cmd_words, "--output-from",$self->region_start ) if($self->region_start);
    push(@cmd_words, "--output-to",$self->region_end ) if($self->region_end);
       
    my $cmd = join(' ', @cmd_words);
    
    $self->execute_command_line ($cmd);
    
    return;
}

=head2 get_ref

  Arg [1]   : ReseqTrack::Tools::RunPhase::PhaseByShapeit
  Function  : Get reference files from $self->reference_config for $self->chrom
              Reference config file format:chrom\tmap\thaps\tlegend\tsamples
              Output stored in $self->ref_map, $self->ref_hap, $self->ref_legend and ref_samples
                          
  Returntype: 
  Exceptions: 
  Example   : $self->get_ref

=cut

sub get_ref {
  my ( $self, $arg ) = @_;
  my $reference_config = $self->reference_config; ###format:chrom\tmap\thaps\tlegend\tsamples
  my $chrom= $self->chrom;
  my $chr_flag=0;

  open(FH,$reference_config)||throw("can't read $reference_config\n");
  while(my $ref_line=<FH>)
  {
    chomp($ref_line);
    next if($ref_line=~ /^#/);
    
    
     
    my @field=split(/\s+/,$ref_line);
        
    if($field[0] eq $chrom)
    {
      $self->ref_map($field[1]);
      $self->ref_hap($field[2]);
      $self->ref_legend($field[3]);
      $self->ref_samples($field[4]);
      $chr_flag++;
    }
       
  }
  
  throw("Chromosome name not present in Reference Config file") unless($chr_flag >0);
  
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

sub exclude_strand_flip {
  my $self = shift;
  if (@_) {
    $self->{'exclude_strand_flip'} = (shift) ? 1 : 0;
  }
  return $self->{'exclude_strand_flip'};
}

=head2 chrX

  Arg [1]   : ReseqTrack::Tools::RunPhase
  Arg [2]   : Boolean (0 or 1)
  Function  : 
  Returntype: Boolean
  Exceptions: n/a
  Example   : 

=cut

sub chrX {
  my $self = shift;
  if (@_) {
    $self->{'chrX'} = (shift) ? 1 : 0;
  }
  return $self->{'chrX'};
}

sub input_samples {
  my ( $self, $arg ) = @_;

  if ($arg) {
    foreach my $file (@{ref($arg) eq 'ARRAY' ? $arg : [$arg]}) {
      $file =~ s{//}{/}g;
      $self->{'input_samples'}->{$file} = 1;
    }
  }

  my @files = keys %{$self->{'input_samples'}};
  return \@files;
}

sub output_haps_files {
    my $self = shift;
    my @files = grep { /\.haps$/ } @{ $self->output_files };
    return \@files;
}

sub output_samples_files {
    my $self = shift;
    my @files = grep { /\.sample$/ } @{ $self->output_files };
    return \@files;
}

1;
