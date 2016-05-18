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