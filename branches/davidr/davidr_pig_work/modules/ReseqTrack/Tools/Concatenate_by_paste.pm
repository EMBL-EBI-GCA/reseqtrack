
=pod

=head1 NAME

ReseqTrack::Tools::Concatenate_by_paste

=head1 SYNOPSIS

This is a wrapper to run unix command 'paste'.

It takes a list of files and concatenate them by columns

example

my $paste = Concatenate_by_paste->new(
     -input_filess        =>"/path/to/bam1 /path/to/bam2" ,
     -program       =>"/path/to/unix/paste",
     -output_dir     => '/path/to/dir/',
);
=cut

package ReseqTrack::Tools::Concatenate_by_paste;

use strict;
use warnings;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists check_executable);
use ReseqTrack::Tools::RunProgram;
use vars qw(@ISA);
@ISA = qw(ReseqTrack::Tools::RunProgram);

sub new {
	my ( $class, @args ) = @_;
	my $self = $class->SUPER::new(@args);
	return $self;
}

sub run_program {
	my $self = shift;

	#print "input files are " ;
	#print join (" ", @{$self->input_files}) . "\n";
	throw "No input files" if (! @{$self->input_files} || @{$self->input_files} == 0 );

	check_file_exists($_) foreach (@{$self->input_files});
	check_executable($self->program);
	
	my @sorted_input_files = sort {$a cmp $b} @{$self->input_files};
	my @sh_input_names = ();
	my $short_name_hash = $self->get_short_input_names(2);
	foreach my $long_input_name ( @sorted_input_files ) {
		push @sh_input_names, $short_name_hash->{$long_input_name};
	}
		
	#print "shortened file names are \n";
	#print join("\n", @sh_input_names) . "\n";
	
	my $output_file = $self->working_dir . "/" . $self->job_name . ".$$.matrix";
	my $cmd = join( ' ', $self->program, @sh_input_names, ' > ', $output_file );

    $self->execute_command_line ($cmd);
    $self->output_files($output_file);
}    