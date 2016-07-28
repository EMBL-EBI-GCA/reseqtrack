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

create or replace view running_jobs_vw as 
select 
e.event_id,
e.name,
count(0) AS total,
(count(0) - count(j.exec_host)) AS queued,
count(j.exec_host) AS running,
min(j.output_file) AS example_output,
max(j.retry_count) AS max_retry,
group_concat(j.input_string separator ',') AS input_strings 
from job j
join event e on j.event_id = e.event_id
group by e.event_id,e.name order by e.name;

CREATE OR REPLACE VIEW complete_jobs AS
select
e.event_id,
e.name,
e.farm_options,
substr(
	e.farm_options,(locate('-M',e.farm_options) + 2),
	(case locate(' ',e.farm_options,locate('-M',e.farm_options))
	when 0 then length(e.farm_options)
	else ((locate(' ',e.farm_options,locate('-M',e.farm_options)) - locate('-M',e.farm_options)) - 2) end)
) mem_limit,
max(ec.memory_usage) max_mem_used,
round(((sum(ec.time_elapsed) / count(0)) / 60),1) avg_time_elapsed,
count(*) jobs_complete
from event e
join event_complete ec on e.event_id = ec.event_id
group by e.name,e.event_id,e.farm_options order by e.name;