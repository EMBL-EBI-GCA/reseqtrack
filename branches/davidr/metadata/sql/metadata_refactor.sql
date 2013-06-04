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
   attribute_name VARCHAR(50) NOT NULL,
   attribute_value VARCHAR(255) NOT NULL,
   PRIMARY KEY (attribute_id),
   key(attribute_name),
   key(other_id, table_name),
   unique(other_id, table_name, attribute_name, attribute_value)
) ENGINE=MYISAM; 




create table study(
    study_id int(10) unsigned primary key auto_increment, 
    source_id  varchar(15) not null,
    status varchar(50) not null ,
    md5     varchar(32),
    type      varchar(100) not null ,
    title varchar(4000)
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
    constraint foreign key (sample_id) references sample(sample_id),
    constraint foreign key (experiment_id) references experiment(experiment_id)
  );

create index run_fk1 on run(sample_id);
create index run_fk2 on run(experiment_id);
create unique index run_src_idx on run(source_id);