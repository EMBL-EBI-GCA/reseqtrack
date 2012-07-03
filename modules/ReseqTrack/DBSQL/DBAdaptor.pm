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
use ReseqTrack::DBSQL::AlignmentMetaInfoAdaptor;
use ReseqTrack::DBSQL::HistoryAdaptor;
use ReseqTrack::DBSQL::StatisticsAdaptor;
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
  if(!$self->{file_adaptor}){
    $self->{file_adaptor} = ReseqTrack::DBSQL::FileAdaptor->
        new($self);
  }
  return $self->{file_adaptor};
}

sub get_StatisticsAdaptor{
  my ($self) = @_;
  if(!$self->{statistics_adaptor}){
    $self->{statistics_adaptor} = ReseqTrack::DBSQL::StatisticsAdaptor->
        new($self);
  }
  return $self->{statistics_adaptor};
}

sub get_HistoryAdaptor{
  my ($self) = @_;
  if(!$self->{history_adaptor}){
    $self->{history_adaptor} = ReseqTrack::DBSQL::HistoryAdaptor->
        new($self);
  }
  return $self->{history_adaptor};
}

sub get_HostAdaptor{
  my ($self) = @_;
  if(!$self->{host_adaptor}){
    $self->{host_adaptor} = ReseqTrack::DBSQL::HostAdaptor->
        new($self);
  }
  return $self->{host_adaptor};
}

sub get_EventAdaptor{
  my ($self) = @_;
  if(!$self->{event_adaptor}){
    $self->{event_adaptor} = ReseqTrack::DBSQL::EventAdaptor->
        new($self);
  }
  return $self->{event_adaptor};
}

sub get_EventCompleteAdaptor{
  my ($self) = @_;
  if(!$self->{event_complete_adaptor}){
    $self->{event_complete_adaptor} = ReseqTrack::DBSQL::EventCompleteAdaptor->
        new($self);
  }
  return $self->{event_complete_adaptor};
}

sub get_JobAdaptor{
  my ($self) = @_;
  if(!$self->{job_adaptor}){
    $self->{job_adaptor} = ReseqTrack::DBSQL::JobAdaptor->
        new($self);
  }
  return $self->{job_adaptor};
}

sub get_WorkflowAdaptor{
  my ($self) = @_;
  if(!$self->{workflow_adaptor}){
    $self->{workflow_adaptor} = ReseqTrack::DBSQL::WorkflowAdaptor->
        new($self);
  }
  return $self->{workflow_adaptor};
}

sub get_CollectionAdaptor{
  my ($self) = @_;
  if(!$self->{collection_adaptor}){
    $self->{collection_adaptor} = ReseqTrack::DBSQL::CollectionAdaptor->
        new($self);
  }
  return $self->{collection_adaptor};
}

sub get_RunMetaInfoAdaptor{
  my ($self) = @_;
  if(!$self->{run_meta_info_adaptor}){
    $self->{run_meta_info_adaptor} = ReseqTrack::DBSQL::RunMetaInfoAdaptor->
        new($self);
  }
  return $self->{run_meta_info_adaptor};
}

sub get_AlignmentMetaInfoAdaptor{
  my ($self) = @_;
  if(!$self->{alignment_meta_info_adaptor}){
    $self->{alignment_meta_info_adaptor} = ReseqTrack::DBSQL::AlignmentMetaInfoAdaptor->
        new($self);
  }
  return $self->{alignment_meta_info_adaptor};
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
  if(!$self->{archive_location_adaptor}){
    $self->{archive_location_adaptor} = ReseqTrack::DBSQL::ArchiveLocationAdaptor->
        new($self);
  }
  return $self->{archive_location_adaptor};
}


sub get_ArchiveActionAdaptor{
  my ($self) = @_;
  if(!$self->{archive_action_adaptor}){
    $self->{archive_action_adaptor} = ReseqTrack::DBSQL::ArchiveActionAdaptor->
        new($self);
  }
  return $self->{archive_action_adaptor};
}

sub get_InputStringAdaptor{
  my ($self) = @_;
  if(!$self->{input_string_adaptor}){
    $self->{input_string_adaptor} = ReseqTrack::DBSQL::InputStringAdaptor->
        new($self);
  }
  return $self->{input_string_adaptor};
}

sub get_MetaAdaptor{
  my ($self) = @_;
  if(!$self->{meta_adaptor}){
    $self->{meta_adaptor} = ReseqTrack::DBSQL::MetaAdaptor->
        new($self);
  }
  return $self->{meta_adaptor};
}

sub get_GenotypeResultsAdaptor{
  my ($self) = @_;
  if(!$self->{genotype_results_adaptor}){
    $self->{genotype_results_adaptor} = ReseqTrack::DBSQL::GenotypeResultsAdaptor->
        new($self);
  }
  return $self->{genotype_results_adaptor};
}

sub get_RejectLogAdaptor{
  my ($self) = @_;
  if(!$self->{rejectlog_adaptor}){
    $self->{rejectlog_adaptor} = ReseqTrack::DBSQL::RejectLogAdaptor->
        new($self);
  }
  return $self->{rejectlog_adaptor};
}

sub get_VerifyBamIDSampleAdaptor{
  my ($self) = @_;
  if(!$self->{VerifyBamIDSampleAdaptor}){
    $self->{VerifyBamIDSampleAdaptor} = ReseqTrack::DBSQL::VerifyBamIDSampleAdaptor->
        new($self);
  }
  return $self->{VerifyBamIDSampleAdaptor};
}

sub get_VerifyBamIDReadGroupAdaptor{
  my ($self) = @_;
  if(!$self->{VerifyBamIDReadGroupAdaptor}){
    $self->{VerifyBamIDReadGroupAdaptor} = ReseqTrack::DBSQL::VerifyBamIDReadGroupAdaptor->
        new($self);
  }

  return $self->{VerifyBamIDReadGroupAdaptor};
}

sub get_FileTypeRuleAdaptor{
  my ($self) = @_;
  if(!$self->{file_type_rule_adaptor}){
    $self->{file_type_rule_adaptor} = ReseqTrack::DBSQL::FileTypeRuleAdaptor->
        new($self);
  }
  return $self->{file_type_rule_adaptor};
}

sub get_PopulationRuleAdaptor{
  my ($self) = @_;
  if(!$self->{population_rule_adaptor}){
    $self->{population_rule_adaptor} = ReseqTrack::DBSQL::PopulationRuleAdaptor->
        new($self);
  }
  return $self->{population_rule_adaptor};
}

sub get_StudyIDAdaptor{
  my ($self) = @_;
  if(!$self->{study_id_adaptor}){
    $self->{study_id_adaptor} = ReseqTrack::DBSQL::StudyIDAdaptor->
        new($self);
  }
  return $self->{study_id_adaptor};
}

sub get_VerifyBamIDAdaptor{
  my ($self) = @_;
  if(!$self->{VerifyBamIDAdaptor}){
    $self->{VerifyBamIDAdaptor} = ReseqTrack::DBSQL::VerifyBamIDAdaptor->
        new($self);
  }

  return $self->{VerifyBamIDAdaptor};
}

1;
