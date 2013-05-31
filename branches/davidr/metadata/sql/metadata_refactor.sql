drop table if exists run;
drop table if exists experiment_attribute;
drop table if exists sample_attribute;
drop table if exists sample;
drop table if exists experiment;
drop table if exists study;

create table study(
    study_id      varchar(15) primary key,
    status varchar(50) not null ,
    md5     varchar(32),
    type      varchar(100) not null ,
    title varchar(4000)
);

create table experiment (
    experiment_id varchar(15) primary key,
    study_id  varchar(15) not null ,
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
    constraint foreign key (study_id) references study(study_id)
);

create index experiment_fk1 on experiment(study_id);
        
create table experiment_attribute (
  experiment_attribute_id int(10) unsigned primary key auto_increment,
  experiment_id varchar(15) not null,
  name varchar(512) not null,
  value varchar(4000),
  constraint foreign key (experiment_id) references experiment(experiment_id)
);

create index experiment_attribute_fk1 on experiment_attribute(experiment_id);
create unique index experiment_attribute_idx1 on experiment_attribute(experiment_id,name);

create table sample
(
    sample_id     varchar(15) primary key ,
    status varchar(50),
    md5     varchar(32),
    center_name     varchar(100) ,
    sample_alias    varchar(500) ,
    tax_id          varchar(15) ,
    scientific_name varchar(500) ,
    common_name     varchar(4000) ,
    anonymized_name varchar(4000) ,
    individual_name varchar(4000) ,
    sample_title varchar(4000)
);

create table sample_attribute (
  sample_attribute_id int(10) unsigned primary key auto_increment,
  sample_id     varchar(15) not null,
  name varchar(512) not null,
  value varchar(4000),
  constraint foreign key (sample_id) references sample(sample_id)
);

create index sample_attribute_fk1 on sample_attribute(sample_id);
create unique index sample_attribute_idx1 on sample_attribute(sample_id,name);

create table run
(
    run_id        varchar(15) primary key ,
    experiment_id  varchar(15) not null,
    sample_id  varchar(15) not null,
    run_alias varchar(500) not null ,
    status varchar(50) not null ,
    md5 varchar(32),
    center_name         varchar(100),
    run_center_name     varchar(100),
    instrument_platform varchar(50),
    instrument_model    varchar(50),
    constraint foreign key (sample_id) references sample(sample_id),
    constraint foreign key (experiment_id) references experiment(experiment_id)
  );

create index run_fk1 on run(sample_id);
create index run_fk2 on run(experiment_id);