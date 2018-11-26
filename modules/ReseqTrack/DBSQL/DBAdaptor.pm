package ReseqTrack::DBSQL::DBAdaptor;
 
use ReseqTrack::DBSQL::DBConnection;

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);

use ReseqTrack::DBSQL::FileAdaptor;
use ReseqTrack::DBSQL::HostAdaptor;
use ReseqTrack::DBSQL::HistoryAdaptor;
use ReseqTrack::DBSQL::EventAdaptor;
use ReseqTrack::DBSQL::EventCompleteAdaptor;
use ReseqTrack::DBSQL::JobAdaptor;
use ReseqTrack::DBSQL::WorkflowAdaptor;
use ReseqTrack::DBSQL::CollectionAdaptor;
use ReseqTrack::DBSQL::RunMetaInfoAdaptor;
use ReseqTrack::DBSQL::HistoryAdaptor;
use ReseqTrack::DBSQL::AttributeAdaptor;
use ReseqTrack::DBSQL::ArchiveAdaptor;
use ReseqTrack::DBSQL::ArchiveActionAdaptor;
use ReseqTrack::DBSQL::ArchiveLocationAdaptor;
use ReseqTrack::DBSQL::InputStringAdaptor;
use ReseqTrack::DBSQL::MetaAdaptor;
use ReseqTrack::DBSQL::GenotypeResultsAdaptor;
use ReseqTrack::DBSQL::VerifyBamIDReadGroupAdaptor;
use ReseqTrack::DBSQL::VerifyBamIDSampleAdaptor;
use ReseqTrack::DBSQL::FileTypeRuleAdaptor;
use ReseqTrack::DBSQL::PopulationRuleAdaptor;
use ReseqTrack::DBSQL::StudyIDAdaptor;
use ReseqTrack::DBSQL::VerifyBamIDAdaptor;
use ReseqTrack::DBSQL::RejectLogAdaptor;
use ReseqTrack::DBSQL::StudyAdaptor;
use ReseqTrack::DBSQL::ExperimentAdaptor;
use ReseqTrack::DBSQL::SampleAdaptor;
use ReseqTrack::DBSQL::RunAdaptor;
use ReseqTrack::DBSQL::PipelineAdaptor;
use ReseqTrack::DBSQL::HiveDBAdaptor;
use ReseqTrack::DBSQL::PipelineSeedAdaptor;
use ReseqTrack::DBSQL::PipelineOutputAdaptor;

sub new {
  my($class, @args) = @_;
  my $self ={};
  bless $self,$class;

  my ($con) = rearrange([qw(DBCONN)], @args);

  if(defined($con)){
    $self->dbc($con);
  }else{
    $self->dbc(new ReseqTrack::DBSQL::DBConnection(@args));
  }
 
  return $self;
}

=head2 dbc

  Arg[1]    : (optional) Bio::EnsEMBL::DBSQL::DBConnection

  Exmaple    : $dbc = $dba->dbc();
  Description: Getter/Setter for DBConnection.
  Returntype : Bio::EnsEMBL::DBSQL::DBConnection
  Exceptions : throws if argument not a Bio::EnsEMBL::DBSQL::DBConnection
  Caller     : general
  Status     : Stable

=cut

sub dbc{
  my $self  = shift;
  
  if(@_){
    my $arg = shift;
    if(defined($arg)){
      if(!$arg->isa('ReseqTrack::DBSQL::DBConnection')){
	throw("$arg is no a DBConnection\n");
      }
    }
    $self->{_dbc} = $arg;
  }
  print $self->{_dbc};
  return $self->{_dbc};
}


=head2 get_XXXAdaptor

  Arg [1]   : ReseqTrack::DBSQL::DBAdaptor
  Function  : return appropriate adaptor object
  Returntype: ReseqTrack::DBSQL::XXXXAdaptor
  Exceptions: none
  Example   : my $fa = $db->get_FileAdaptor;

=cut


sub get_FileAdaptor{
		my ($self) = @_;
	return $self->_get_adaptor('ReseqTrack::DBSQL::FileAdaptor');
}

sub get_AttributeAdaptor{
		my ($self) = @_;
	return $self->_get_adaptor('ReseqTrack::DBSQL::AttributeAdaptor');
  
}

sub get_HistoryAdaptor{
		my ($self) = @_;
	return $self->_get_adaptor('ReseqTrack::DBSQL::HistoryAdaptor');
 
}

sub get_HostAdaptor{
		my ($self) = @_;
	return $self->_get_adaptor('ReseqTrack::DBSQL::HostAdaptor');

}

sub get_EventAdaptor{
		my ($self) = @_;
	return $self->_get_adaptor('ReseqTrack::DBSQL::EventAdaptor');
}

sub get_EventCompleteAdaptor{
		my ($self) = @_;
	return $self->_get_adaptor('ReseqTrack::DBSQL::EventCompleteAdaptor');
}

sub get_JobAdaptor{
		my ($self) = @_;
	return $self->_get_adaptor('ReseqTrack::DBSQL::JobAdaptor');
}

sub get_WorkflowAdaptor{
		my ($self) = @_;
	return $self->_get_adaptor('ReseqTrack::DBSQL::WorkflowAdaptor');
}

sub get_CollectionAdaptor{
		my ($self) = @_;
	return $self->_get_adaptor('ReseqTrack::DBSQL::CollectionAdaptor');
}

sub get_RunMetaInfoAdaptor{
		my ($self) = @_;
	return $self->_get_adaptor('ReseqTrack::DBSQL::RunMetaInfoAdaptor');
}

sub get_ArchiveAdaptor{
  my ($self, $nolock) = @_;
  if(!$self->{archive_adaptor}){
    $self->{archive_adaptor} = ReseqTrack::DBSQL::ArchiveAdaptor->
        new($self, $nolock);
  }
  return $self->{archive_adaptor};
}


sub get_ArchiveLocationAdaptor{
		my ($self) = @_;
	return $self->_get_adaptor('ReseqTrack::DBSQL::ArchiveLocationAdaptor');
}

sub get_ArchiveActionAdaptor{
		my ($self) = @_;
	return $self->_get_adaptor('ReseqTrack::DBSQL::ArchiveActionAdaptor');
}

sub get_InputStringAdaptor{
		my ($self) = @_;
	return $self->_get_adaptor('ReseqTrack::DBSQL::InputStringAdaptor');
}

sub get_MetaAdaptor{
		my ($self) = @_;
	return $self->_get_adaptor('ReseqTrack::DBSQL::MetaAdaptor');
}

sub get_GenotypeResultsAdaptor{
		my ($self) = @_;
	return $self->_get_adaptor('ReseqTrack::DBSQL::GenotypeResultsAdaptor');
}

sub get_RejectLogAdaptor{
		my ($self) = @_;
	return $self->_get_adaptor('ReseqTrack::DBSQL::RejectLogAdaptor');
}

sub get_VerifyBamIDSampleAdaptor{
		my ($self) = @_;
	return $self->_get_adaptor('ReseqTrack::DBSQL::VerifyBamIDSampleAdaptor');
}

sub get_VerifyBamIDReadGroupAdaptor{
		my ($self) = @_;
	return $self->_get_adaptor('ReseqTrack::DBSQL::VerifyBamIDReadGroupAdaptor');
}

sub get_FileTypeRuleAdaptor{
		my ($self) = @_;
	return $self->_get_adaptor('ReseqTrack::DBSQL::FileTypeRuleAdaptor');
}

sub get_PopulationRuleAdaptor{
		my ($self) = @_;
	return $self->_get_adaptor('ReseqTrack::DBSQL::PopulationRuleAdaptor');
}

sub get_StudyIDAdaptor{
		my ($self) = @_;
	return $self->_get_adaptor('ReseqTrack::DBSQL::StudyIDAdaptor');
}

sub get_VerifyBamIDAdaptor{
		my ($self) = @_;
	return $self->_get_adaptor('ReseqTrack::DBSQL::VerifyBamIDAdaptor');
}

sub get_StudyAdaptor{
	my ($self) = @_;
	return $self->_get_adaptor('ReseqTrack::DBSQL::StudyAdaptor');
}
sub get_ExperimentAdaptor{
	my ($self) = @_;
	return $self->_get_adaptor('ReseqTrack::DBSQL::ExperimentAdaptor');
}
sub get_RunAdaptor{
	my ($self) = @_;
	return $self->_get_adaptor('ReseqTrack::DBSQL::RunAdaptor');
}
sub get_SampleAdaptor{
	my ($self) = @_;
	return $self->_get_adaptor('ReseqTrack::DBSQL::SampleAdaptor');
}
sub get_PipelineAdaptor{
	my ($self) = @_;
	return $self->_get_adaptor('ReseqTrack::DBSQL::PipelineAdaptor');
}
sub get_HiveDBAdaptor{
	my ($self) = @_;
	return $self->_get_adaptor('ReseqTrack::DBSQL::HiveDBAdaptor');
}
sub get_PipelineSeedAdaptor{
	my ($self) = @_;
	return $self->_get_adaptor('ReseqTrack::DBSQL::PipelineSeedAdaptor');
}
sub get_PipelineOutputAdaptor{
	my ($self) = @_;
	return $self->_get_adaptor('ReseqTrack::DBSQL::PipelineOutputAdaptor');
}
sub _get_adaptor{
	my ($self,$class,@args) = @_;
	if (!$self->{$class}){
		$self->{$class} = $class->new($self,@args);
	}
	return $self->{$class};
}

sub get_adaptor_for_table{
  my ($self, $table_name) = @_;
  throw("undefined table name") if ! defined $table_name;
  return $self->get_FileAdaptor if ($table_name eq 'file');
  return $self->get_CollectionAdaptor if ($table_name eq 'collection');
  return $self->get_EventAdaptor if ($table_name eq 'event');
  return $self->get_RunMetaInfoAdaptor if ($table_name eq 'run_meta_info');
  return $self->get_InputStringAdaptor if ($table_name eq 'input_string');
  return $self->get_SampleAdaptor if ($table_name eq 'sample');
  return $self->get_RunAdaptor if ($table_name eq 'run');
  return $self->get_ExperimentAdaptor if ($table_name eq 'experiment');
  return $self->get_StudyAdaptor if ($table_name eq 'study');
  throw("do not know what to return for $table_name (yet)");
}

1;
