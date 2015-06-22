=pod

=head1 NAME

ReseqTrack::Tools::GATKTools::LeftAlignVariants

=head1 SYNOPSIS

Object to create a bam file that in realigned around
known indel sites using GATK RealignerTargetCreator and
IndelRealigner

example

my $REALIGN_AROUND_INDELS = ReseqTrack::Tools::GATKTools::LeftAlignVariants->new(
     -variants_file     => '/path/to/input.vcf'
     -java_exe        =>"/usr/bin/java" ,
     -jvm_args        =>"-Xmx4g",
     -GATK_PATH       =>$GATK/GenomeAnalysisTK/",
     -working_dir     => '/path/to/dir/',
     -reference       => '/path/to/reference',
);


=cut

package ReseqTrack::Tools::GATKTools::LeftAlignVariants;

use strict;
use warnings;

use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists check_executable);

use base qw(ReseqTrack::Tools::GATKTools);

sub DEFAULT_OPTIONS { return {
        };
}

sub new {
	my ( $class, @args ) = @_;
	my $self = $class->SUPER::new(@args);

	my ( $variants_file,)
        	= rearrange( [ qw( VARIANTS_FILE )], @args);

        $self->variants_file($variants_file);
        
	return $self;
}

sub run_program {
	my $self = shift;

        throw "no input vcf to left align" if (!$self->variants_file);
        
        check_file_exists($self->variants_file);
        $self->check_jar_file_exists;
        check_file_exists($self->reference);
        check_executable($self->java_exe);

        my $leftaligned_vcf = $self->working_dir."/".$self->job_name.".leftaligned.vcf.gz";
        $leftaligned_vcf =~ s{//}{/}g;

        my @cmd_words = ($self->java_exe, $self->jvm_args, '-jar');
        push(@cmd_words, $self->gatk_path . '/' . $self->jar_file);
        push(@cmd_words, '-T', 'LeftAlignVariants');
        push(@cmd_words, '--variant', $self->variants_file);
        push(@cmd_words, '-R', $self->reference);
        push(@cmd_words, '-o', $leftaligned_vcf);

        my $cmd = join(' ', @cmd_words);

        $self->output_files($leftaligned_vcf);
        $self->execute_command_line ($cmd);
	    return;
}

sub variants_file {

  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{variants_file} = $arg;
  }
  return $self->{variants_file};
}

1;
