drop table if exists run;
drop table if exists experiment_attribute;
drop table if exists sample_attribute;
drop table if exists sample;
drop table if exists experiment;
drop table if exists study;
drop table if exists attribute;

CREATE TABLE attribute(
   attribute_id int(10) unsigned NOT  NULL AUTO_INCREMENT,
   table_name  enum('file', 'event', 'run_meta_info', 'alignment_meta_info', 'collection','run','sample','experiment','study'),
   other_id int(10) unsigned NOT  NULL,      
   attribute_name VARCHAR(100) NOT NULL,
   attribute_value VARCHAR(4000) NOT NULL,
   PRIMARY KEY (attribute_id),
   key(attribute_name),
   key(other_id, table_name),
   unique(other_id, table_name, attribute_name)
) ENGINE=MYISAM; 

insert into attribute (attribute_id,table_name,other_id,attribute_name,attribute_value)
select statistics_id,table_name,other_id,attribute_name,attribute_value from statistics;


create table study(
    study_id int(10) unsigned primary key auto_increment, 
    source_id  varchar(15) not null,
    status varchar(50) not null ,
    md5     varchar(32),
    type      varchar(100) not null ,
    submission_id varchar(15),
    submission_date datetime,
    title varchar(2000),
    study_alias varchar(500)
);

create unique index study_src_idx on study(source_id);

create table experiment (
    experiment_id int(10) unsigned primary key auto_increment,
    source_id varchar(15) not null,
    study_id  int(10) unsigned not null ,
    status varchar(50) not null ,
    md5     varchar(32),
    center_name           varchar(100) ,
    experiment_alias      varchar(500) ,
    instrument_platform   varchar(50) not null ,
    instrument_model      varchar(50) ,
    library_layout        varchar(50) not null ,
    library_name          varchar(500) ,
    library_strategy      varchar(50) ,
    library_source        varchar(50) not null ,
    library_selection     varchar(50) not null ,
    paired_nominal_length int(10) ,
    paired_nominal_sdev   int(10),
    submission_id varchar(15),
    submission_date datetime,
    constraint foreign key (study_id) references study(study_id)
);

create unique index experiment_src_idx on experiment(source_id);
create index experiment_fk1 on experiment(study_id);

create table sample
(
    sample_id int(10) unsigned primary key auto_increment,
    source_id varchar(15) not null,
    status varchar(50),
    md5     varchar(32),
    center_name     varchar(100) ,
    sample_alias    varchar(500) ,
    tax_id          varchar(15) ,
    scientific_name varchar(500) ,
    common_name     varchar(4000) ,
    anonymized_name varchar(4000) ,
    individual_name varchar(4000) ,
    submission_id varchar(15),
    submission_date datetime,
    sample_title varchar(4000)
);


create unique index sample_src_idx on sample(source_id);

create table run
(
    run_id        int(10) unsigned primary key auto_increment,
    source_id	varchar(15) not null,
    experiment_id  int(10) unsigned not null,
    sample_id  int(10) unsigned not null,
    run_alias varchar(500) not null ,
    status varchar(50) not null ,
    md5 varchar(32),
    center_name         varchar(100),
    run_center_name     varchar(100),
    instrument_platform varchar(50),
    instrument_model    varchar(50),
    submission_id varchar(15),
    submission_date datetime,
    constraint foreign key (sample_id) references sample(sample_id),
    constraint foreign key (experiment_id) references experiment(experiment_id)
  );

create index run_fk1 on run(sample_id);
create index run_fk2 on run(experiment_id);
create unique index run_src_idx on run(source_id);

create or replace view run_meta_info_vw as
select
r.run_id as run_meta_info_id,
r.source_id as run_id,
st.source_id as study_id,
st.study_alias as study_name,
r.center_name,
r.submission_id,
r.submission_date,
s.source_id as sample_id,
s.sample_alias as sample_name,
pop.attribute_value as population,
e.source_id as experiment_id,
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
join sample s on r.sample_id = s.sample_id
left outer join attribute pop on pop.table_name = 'sample' and pop.attribute_name = 'POPULATION' and pop.other_id = s.sample_id
left outer join attribute rc on rc.table_name = 'run' and rc.attribute_name = 'READ_COUNT' and rc.other_id = r.run_id
left outer join attribute bc on bc.table_name = 'run' and bc.attribute_name = 'BASE_COUNT' and bc.other_id = r.run_id
;
