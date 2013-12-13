package ReseqTrack::Tools::QC::ChkIndelsBAM;

use strict;
use warnings;
use vars qw(@ISA);
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::FileSystemUtils qw(get_lines_from_file);
use File::Path;
use File::Basename;
use base qw(ReseqTrack::Tools::RunProgram);
use ReseqTrack::DBSQL::StatisticsAdaptor;

=pod

=head1 NAME

ReseqTrack::Tools::QC::ChkIndelsBAM

=head1 SYNOPSIS

Object for running Heng's chk_indel_rg code, which checks BAM files for excessive amount of short indels, which indicates 
low data quality.

The output file of the module is input_file.out in the outdir specified. The header of the out put file is:
file name, RG-ID, #<=6bp-ins, #<=6bp-del, #>6bp-ins and #>6bp-del, 

All results are also stored in the statistic table. 

For 1KG BAMs, anything with ratio of col3/col4 or col4/col3 > 5 is classified as bad. The results can be queried using script query_chk_indel_results.pl

=head1 Example

my $chk_indels = ReseqTrack::Tools::QC::ChkIndelsBAM (
					  -db				=> $db,
                      -input_files 		=> '/path/to/file',
                      -working_dir		=> '/path/to/output_dir
                      -store			=> $store,
                      );
$chk_indels->run;
my $output_file = $chk_indels->output_files->[0];

=cut


sub new {
	my ( $class, @args ) = @_;
	my $self = $class->SUPER::new(@args);
	my ($db, $store) = rearrange([qw(DB STORE)], @args);
	$self->db($db);
	$self->store($store);
	return $self;
}

sub run_program {
    my ($self) = @_;
       
    mkpath $self->working_dir unless (-e $self->working_dir);

	print "input file is " . $self->input_files->[0] . "\n";
	my $outfile = $self->working_dir . "/" . basename($self->input_files->[0]) . ".out";
	my $command = $self->program . " " . $self->input_files->[0] . " > " . $outfile;
	
	ReseqTrack::Tools::RunProgram->execute_command_line($command);
	
	$self->output_files($outfile);
	
	$self->_insert_results_into_db;
    return $self;
}

sub _insert_results_into_db {
	my ($self) = @_;
    
	my $rmi_a = $self->db->get_RunMetaInfoAdaptor;
	my $stats_a = $self->db->get_StatisticsAdaptor;
	
	my $result_file = $self->output_files->[0];
	my $lines = get_lines_from_file($result_file);

	foreach my $line ( @$lines ) {
		
		my @data = split (/\t/, $line);
		my $rmi_obj = $rmi_a->fetch_by_run_id($data[1]);	
				
		my $stats_obj_1 = ReseqTrack::Statistic->new(
			-table_name => 'run_meta_info',
			-other_id	=> $rmi_obj->dbID,
			-attribute_name	=> 'sh_insert',
			-attribute_value	=> $data[2],
		);	
		my $stats_obj_2 = ReseqTrack::Statistic->new(
			-table_name => 'run_meta_info',
			-other_id	=> $rmi_obj->dbID,
			-attribute_name	=> 'sh_del',
			-attribute_value	=> $data[3],
		);	
		my $stats_obj_3 = ReseqTrack::Statistic->new(
			-table_name => 'run_meta_info',
			-other_id	=> $rmi_obj->dbID,
			-attribute_name	=> 'reg_insert',
			-attribute_value	=> $data[4],
		);	
		my $stats_obj_4 = ReseqTrack::Statistic->new(
			-table_name => 'run_meta_info',
			-other_id	=> $rmi_obj->dbID,
			-attribute_name	=> 'reg_del',
			-attribute_value	=> $data[5],
		);	
			
		$stats_a->store($stats_obj_1) if $self->store;	
		$stats_a->store($stats_obj_2) if $self->store;	
		$stats_a->store($stats_obj_3) if $self->store;	
		$stats_a->store($stats_obj_4) if $self->store;	
		
	}	
		return $self;
}

sub db{
  my ($self, $db) = @_;
  if($db){
    $self->{db} = $db;
  }
  return $self->{db};
}

sub store{
  my ($self, $store) = @_;
  if($store){
    $self->{store} = $store;
  }
  return $self->{store};
}

1;