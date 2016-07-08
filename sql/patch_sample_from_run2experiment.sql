#safety check - problem if this returns anything, do not proceed
select
count(distinct r.sample_id),
group_concat(distinct r.sample_id),
r.experiment_id
from run r
group by r.experiment_id
having count(distinct r.sample_id) > 1;

#add to experiment table
alter table experiment add column sample_id  int(10) unsigned ;
create index experiment_fk2 on experiment(sample_id);
#update experiment with existing data
update experiment e
set sample_id = (select min(r.sample_id) from run r where r.experiment_id = e.experiment_id);
#add constraint now we have data
alter table experiment modify sample_id int(10) unsigned not null;
#drop old column on run
drop index run_fk1 on run;
alter table run drop column sample_id;

#extend sample to have ref to BioSamples ID
alter table sample add column biosample_id varchar(15);
alter table sample add column biosample_authority varchar(2);
create index sample_bs_idx on sample(biosample_id);
  
#alter attribute to have units
alter table attribute add column attribute_units varchar(4000);

create or replace view run_meta_info_vw as
select
r.run_id as run_meta_info_id,
r.run_source_id as run_id,
st.study_source_id as study_id,
st.study_alias as study_name,
r.center_name,
r.submission_id,
r.submission_date,
s.sample_source_id as sample_id,
s.sample_alias as sample_name,
pop.attribute_value as population,
e.experiment_source_id as experiment_id,
e.instrument_platform,
e.instrument_model,
e.library_name,
r.run_alias as run_name,
null as run_block_name,
e.paired_nominal_length as paired_length,
e.library_layout,
r.status, 
bc.attribute_value as archive_base_count,
rc.attribute_value as archive_read_count,
e.library_strategy
from study st
join experiment e on st.study_id = e.study_id
join run r on e.experiment_id = r.experiment_id
join sample s on e.sample_id = s.sample_id
left outer join attribute pop on pop.table_name = 'sample' and pop.attribute_name = 'POPULATION' and pop.other_id = s.sample_id
left outer join attribute rc on rc.table_name = 'run' and rc.attribute_name = 'READ_COUNT' and rc.other_id = r.run_id
left outer join attribute bc on bc.table_name = 'run' and bc.attribute_name = 'BASE_COUNT' and bc.other_id = r.run_id
;

#finally, run this command to add units to your attributes
#perl $RESEQTRACK/scripts/metadata/load_from_ena.pl $RESEQTRACK_DB_ARGS -era_dbuser $ERADBUSER -era_dbpass $ERADBPASS -era_dbname $ERADBNAME -update_existing -force_update