CREATE TABLE file(
       file_id int(10) unsigned NOT NULL AUTO_INCREMENT,
       name VARCHAR(968) NOT NULL,
       md5   VARCHAR(32),       
       type  VARCHAR(50),
       size  BIGINT unsigned NOT NULL default 0, 
       host_id int(10) unsigned NOT NULL,
       withdrawn tinyint(1) NOT NULL DEFAULT '0',
       created   datetime NOT NULL,
       updated   datetime,    

       PRIMARY KEY (file_id),
       KEY file_path_idx (name),
       UNIQUE (name, md5)
  
);

CREATE TABLE host(
       host_id int(10) unsigned NOT NULL AUTO_INCREMENT,
       name  VARCHAR(255),
       remote tinyint(1) NOT NULL DEFAULT '0',
       dropbox_dir VARCHAR(200),
       PRIMARY KEY(host_id),
       UNIQUE(name)
);

create table history(
   history_id int(10) unsigned NOT NULL AUTO_INCREMENT,
   other_id int(10) unsigned NOT NULL,
   table_name enum('file', 'collection', 'event', 'run_meta_info', 'alignment_meta_info'),
   comment VARCHAR(255) NOT NULL,  
   time   datetime NOT NULL, 
   PRIMARY KEY(history_id), 
   key(other_id, table_name)
);



CREATE TABLE event(
       event_id int(10) unsigned NOT NULL AUTO_INCREMENT,
       name VARCHAR(100) NOT NULL,
       program VARCHAR(255),
       program_version VARCHAR(50),
       options VARCHAR(1000),
       input_flag  VARCHAR(40),
       farm_options VARCHAR(255),
       batch_size   MEDIUMINT(10) unsigned NOT NULL default 0,       
       output_path  VARCHAR(255),
       type    VARCHAR(50),
       table_name enum('file', 'collection', 'run_meta_info', 'input_string'),
       created   datetime NOT NULL,
       updated   datetime, 
       
       PRIMARY KEY(event_id),
       UNIQUE(name)   
);


CREATE TABLE event_complete(
   event_id  int(10) unsigned NOT NULL,
   other_id  int(10) unsigned NOT NULL,
   table_name enum('file', 'collection', 'run_meta_info', 'input_string'),
   type VARCHAR(50),       
   success int(1) unsigned NOT NULL,
   time datetime NOT NULL,

   unique(event_id, other_id, table_name)
);

CREATE TABLE input_string(
   input_string_id  int(10) unsigned NOT NULL AUTO_INCREMENT,
   name VARCHAR(949) NOT NULL,
   type VARCHAR(50),       
   
   primary key(input_string_id),
   unique(name, type)
);
CREATE TABLE workflow_goal(
    workflow_id int(10) unsigned NOT NULL AUTO_INCREMENT,
    goal_event_id int(10) unsigned NOT NULL,
    PRIMARY KEY(workflow_id),
    KEY(goal_event_id),
    UNIQUE(goal_event_id)
);

CREATE TABLE workflow_conditions(
    workflow_id int(10) unsigned NOT  NULL,
    conditional_event_id int( 10) unsigned NOT NULL,
    PRIMARY KEY(workflow_id),
    KEY(conditional_event_id),
    KEY workflow_condition(workflow_id, conditional_event_id)
);

CREATE TABLE collection(
   collection_id int(10) unsigned NOT  NULL AUTO_INCREMENT,
   name VARCHAR(255) NOT NULL,
   type  VARCHAR(50) NOT NULL,
   table_name  enum('file', 'run_meta_info', 'alignment_meta_info'),
   PRIMARY KEY (collection_id),
   key(name),
   unique(name, table_name, type)
);

CREATE TABLE collection_group(
  collection_id int(10) unsigned NOT  NULL,
  other_id int(10) unsigned NOT  NULL,
  unique(collection_id, other_id)
);


CREATE TABLE  run_meta_info(
       run_meta_info_id int(10) unsigned NOT NULL AUTO_INCREMENT,
       run_id VARCHAR(15) NOT NULL,
       study_id VARCHAR(20) NOT NULL,
       study_name VARCHAR(500),
       center_name VARCHAR(15),
       submission_id VARCHAR(20) NOT NULL,
       submission_date datetime,
       sample_id VARCHAR(20) NOT NULL,
       sample_name VARCHAR(20) NOT NULL,
       population VARCHAR(50),
       experiment_id VARCHAR(20) NOT NULL,
       instrument_platform VARCHAR(50) NOT NULL,
       instrument_model VARCHAR(100),
       library_name VARCHAR(255) NOT NULL,
       run_name VARCHAR(255),
       run_block_name VARCHAR(255),
       paired_length int(10),
       library_layout VARCHAR(10),
       status VARCHAR(50),     
       archive_base_count bigint,
       archive_read_count bigint,
       PRIMARY KEY(run_meta_info_id),
       KEY (run_id),
       KEY sample_run_idx(run_id, sample_name),
       UNIQUE(run_id)           
);


create table alignment_meta_info(
      alignment_meta_info_id int(10) unsigned NOT NULL AUTO_INCREMENT,
      file_id int(10) unsigned NOT NULL,
      index_file_id int(10) unsigned NOT NULL,
      sample_name  VARCHAR(20) NOT NULL,
      region VARCHAR(20) NOT NULL,
      assembly VARCHAR(20) NOT NULL,
      program VARCHAR(50) NOT NULL,
      technology VARCHAR(30) NOT NULL,     
      mapped_basecount int(10),
      PRIMARY KEY(alignment_meta_info_id),
      KEY(file_id),
      KEY sample_file_idx(sample_name, file_id)
);


CREATE TABLE job(
   job_id            int(10) unsigned NOT NULL AUTO_INCREMENT,
   submission_id     mediumint(10) unsigned,   
   event_id        SMALLINT UNSIGNED NOT NULL,
   input_string      VARCHAR(255) NOT NULL, 
   stdout_file       varchar(200) NOT NULL,
   stderr_file       varchar(200) NOT NULL,
   exec_host           varchar(20),
   retry_count       tinyint(2) unsigned default 0,

   PRIMARY KEY(job_id),
   KEY(input_string),
   KEY(event_id)

);

CREATE TABLE job_status (
  job_id            int(10) unsigned NOT NULL,
  status            varchar(40) DEFAULT 'CREATED' NOT NULL,
  time              datetime NOT NULL,
  is_current        enum('n', 'y') DEFAULT 'n',

  KEY (job_id),
  KEY (status),
  KEY (is_current)
);

CREATE TABLE statistics(
   statistics_id int(10) unsigned NOT  NULL AUTO_INCREMENT,
   table_name  enum('file', 'event', 'run_meta_info', 'alignment_meta_info', 'collection'),
   other_id int(10) unsigned NOT  NULL,      
   attribute_name VARCHAR(50) NOT NULL,
   attribute_value VARCHAR(255) NOT NULL,
   PRIMARY KEY (statistics_id),
   key(attribute_name),
   key(other_id, table_name),
   unique(other_id, table_name, attribute_name, attribute_value)
); 

CREATE TABLE archive(
       archive_id int(10) unsigned NOT NULL AUTO_INCREMENT,
       name  VARCHAR(255) NOT NULL,
       file_id  int(10) unsigned NOT NULL,  
       md5   CHAR(32) NOT NULL,
       size  BIGINT unsigned NOT NULL default 0,       
       relative_path  VARCHAR(700) NOT NULL,
       volume_name VARCHAR(10),
       created   datetime NOT NULL,
       updated   datetime,      
       priority  smallint unsigned default null,
       new_name varchar(255),      
       new_relative_path VARCHAR(700),
       archive_action_id TINYINT(2) unsigned,
       archive_location_id TINYINT(2) unsigned NOT NULL,
       fire_action_id int(15) unsigned unique default null,
       fire_exit_code int(5) unsigned default null,
       fire_exit_reason varchar(4000) default null,
       PRIMARY KEY (archive_id),
       KEY file_name_idx (name),
       UNIQUE (name, relative_path)
);

CREATE TABLE archive_action(
      archive_action_id int(10) unsigned NOT NULL,
      action VARCHAR(50)
);

CREATE TABLE archive_location(
      archive_location_id int(10) unsigned NOT NULL,
      location VARCHAR(50),
      location_name VARCHAR(10)
);


CREATE TABLE era_meta_info( 
       id_string VARCHAR(50),  
       column_name VARCHAR(50)
       );

CREATE TABLE meta (

  meta_id                     INT NOT NULL AUTO_INCREMENT,
  meta_key                    VARCHAR(40) NOT NULL,
  meta_value                  VARCHAR(255) BINARY NOT NULL,

  PRIMARY  KEY (meta_id),    
  unique (meta_key)

) ;

CREATE TABLE reject_log (
    id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
    file_id INT NOT NULL,
    is_reject enum("y", "n") NOT NULL DEFAULT "n",
    reject_reason VARCHAR(500),
    created TIMESTAMP(8)
    );


#Now to add entries to the two standard tables

INSERT INTO archive_action (archive_action_id, action) values(1, 'archive');
INSERT INTO archive_action (archive_action_id, action) values(2, 'dearchive');
INSERT INTO archive_action (archive_action_id, action) values(3, 'replace');
INSERT INTO archive_action (archive_action_id, action) values(4, 'move_within_volume');
INSERT INTO archive_location (archive_location_id, location, location_name) values(1, '/nfs/1000g-work/G1K/archive_staging', 'staging');
INSERT INTO archive_location (archive_location_id, location, location_name) values(2, '/nfs/1000g-archive', 'archive');

