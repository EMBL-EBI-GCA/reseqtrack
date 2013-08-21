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
  
) ENGINE=MYISAM;

CREATE TABLE host(
       host_id int(10) unsigned NOT NULL AUTO_INCREMENT,
       name  VARCHAR(255),
       remote tinyint(1) NOT NULL DEFAULT '0',
       dropbox_dir VARCHAR(200),
       PRIMARY KEY(host_id),
       UNIQUE(name)
) ENGINE=MYISAM;

create table history(
   history_id int(10) unsigned NOT NULL AUTO_INCREMENT,
   other_id int(10) unsigned NOT NULL,
   table_name enum('file', 'collection', 'event', 'run_meta_info', 'alignment_meta_info', 'pipeline'),
   comment VARCHAR(65000) NOT NULL,  
   time   datetime NOT NULL, 
   PRIMARY KEY(history_id), 
   key(other_id, table_name)
) ENGINE=MYISAM;



CREATE TABLE event(
       event_id int(10) unsigned NOT NULL AUTO_INCREMENT,
       name VARCHAR(100) NOT NULL,
       program VARCHAR(255),
       options VARCHAR(1000),
       input_flag  VARCHAR(40),
       farm_options VARCHAR(255),
       runner_options VARCHAR(255),
       max_array_size MEDIUMINT(10) unsigned NOT NULL default 0,
       job_slot_limit MEDIUMINT(10) unsigned,
       output_path  VARCHAR(255),
       type    VARCHAR(50),
       table_name enum('file', 'collection', 'run_meta_info', 'input_string'),
       created   datetime NOT NULL,
       updated   datetime, 
       
       PRIMARY KEY(event_id),
       UNIQUE(name)   
) ENGINE=MYISAM;


CREATE TABLE event_complete(
   event_id  int(10) unsigned NOT NULL,
   other_id  int(10) unsigned NOT NULL,
   table_name enum('file', 'collection', 'run_meta_info', 'input_string'),
   success int(1) unsigned NOT NULL,
   time datetime NOT NULL,
   time_elapsed int(10) unsigned,
   memory_usage MEDIUMINT(10) unsigned,
   swap_usage MEDIUMINT(10) unsigned,
   exec_host varchar(20),

   unique(event_id, other_id, table_name)
) ENGINE=MYISAM;

CREATE TABLE input_string(
   input_string_id  int(10) unsigned NOT NULL AUTO_INCREMENT,
   name VARCHAR(949) NOT NULL,
   type VARCHAR(50),       
   
   primary key(input_string_id),
   unique(name, type)
) ENGINE=MYISAM;
CREATE TABLE workflow_goal(
    workflow_id int(10) unsigned NOT NULL AUTO_INCREMENT,
    goal_event_id int(10) unsigned NOT NULL,
    PRIMARY KEY(workflow_id),
    KEY(goal_event_id),
    UNIQUE(goal_event_id)
) ENGINE=MYISAM;

CREATE TABLE workflow_conditions(
    workflow_id int(10) unsigned NOT  NULL,
    conditional_event_id int( 10) unsigned NOT NULL,
    PRIMARY KEY(workflow_id,conditional_event_id)
) ENGINE=MYISAM;

CREATE TABLE collection(
   collection_id int(10) unsigned NOT  NULL AUTO_INCREMENT,
   name VARCHAR(255) NOT NULL,
   type  VARCHAR(50) NOT NULL,
   table_name  enum('collection', 'file', 'run_meta_info', 'alignment_meta_info'),
   PRIMARY KEY (collection_id),
   key(name),
   unique(name, table_name, type)
) ENGINE=MYISAM;

CREATE TABLE collection_group(
  collection_id int(10) unsigned NOT  NULL,
  other_id int(10) unsigned NOT  NULL,
  unique(collection_id, other_id)
) ENGINE=MYISAM;


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
	   library_strategy varchar(32),
       PRIMARY KEY(run_meta_info_id),
       KEY (run_id),
       KEY sample_run_idx(run_id, sample_name),
       UNIQUE(run_id)           
) ENGINE=MYISAM;


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
) ENGINE=MYISAM;


CREATE TABLE job(
   job_id            int(10) unsigned NOT NULL AUTO_INCREMENT,
   submission_id     mediumint(10) unsigned,   
   submission_index  mediumint(10) unsigned,   
   event_id        SMALLINT UNSIGNED NOT NULL,
   input_string      VARCHAR(255) NOT NULL, 
   output_file   varchar(200) NOT NULL,
   farm_log_file varchar(200) NOT NULL,
   exec_host           varchar(20),
   retry_count       tinyint(2) unsigned default 0,

   PRIMARY KEY(job_id),
   KEY(input_string),
   KEY(event_id)

) ENGINE=MYISAM;

CREATE TABLE job_status (
  job_id            int(10) unsigned NOT NULL,
  status            varchar(40) DEFAULT 'CREATED' NOT NULL,
  time              datetime NOT NULL,
  is_current        enum('n', 'y') DEFAULT 'n',

  KEY (job_id),
  KEY (status),
  KEY (is_current)
) ENGINE=MYISAM;

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
) ENGINE=MYISAM; 

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
) ENGINE=MYISAM;

CREATE TABLE archive_action(
      archive_action_id int(10) unsigned NOT NULL,
      action VARCHAR(50)
) ENGINE=MYISAM;

CREATE TABLE archive_location(
      archive_location_id int(10) unsigned NOT NULL,
      location VARCHAR(50),
      location_name VARCHAR(10)
) ENGINE=MYISAM;


CREATE TABLE era_meta_info( 
       id_string VARCHAR(50),  
       column_name VARCHAR(50)
       ) ENGINE=MYISAM;

CREATE TABLE meta (

  meta_id                     INT NOT NULL AUTO_INCREMENT,
  meta_key                    VARCHAR(40) NOT NULL,
  meta_value                  VARCHAR(255) BINARY NOT NULL,

  PRIMARY  KEY (meta_id),    
  unique (meta_key)

)  ENGINE=MYISAM;

#CREATE TABLE reject_log (
#    id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
#    file_id INT NOT NULL,
#    is_reject enum("y", "n") NOT NULL DEFAULT "n",
#    reject_reason VARCHAR(500),
#    created TIMESTAMP
#    ) ENGINE=MYISAM;


CREATE TABLE `genotype_results` (
  `genotype_results_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `table_name` enum('file','collection','run_meta_info','input_string') DEFAULT NULL,
  `other_id` int(10) unsigned NOT NULL,
  `name` varchar(15) NOT NULL,
  `claimed` varchar(15) NOT NULL,
  `top_hit` varchar(15) NOT NULL,
  `second_hit` varchar(15) NOT NULL,
  `ratio_2_to_1` decimal(6,2) NOT NULL,
  `ratio_claimed` decimal(6,2) NOT NULL,
  `reference` varchar(1000) NOT NULL,
  `snps_bin` varchar(1000) NOT NULL,
  `aligner` varchar(1000) NOT NULL,
  `validation_method` varchar(100) NOT NULL,
  `max_bases` varchar(40) NOT NULL,
  `percent_mapped` int(4) unsigned NOT NULL,
  `verdict` varchar(30) NOT NULL,
  `cfg_file` varchar(350) NOT NULL,
  `performed` datetime NOT NULL,
  PRIMARY KEY (`genotype_results_id`),
  UNIQUE KEY `name` (`name`)
)  ENGINE=MYISAM;

CREATE TABLE file_type_rule (
       rule_block_order int(10) unsigned NOT NULL,
       rule_order int(10) unsigned NOT NULL,
       file_type  VARCHAR(50),
       match_regex  VARCHAR(1000),

       PRIMARY KEY (rule_block_order, rule_order)
) ENGINE=MYISAM;

CREATE TABLE population_rule (
       rule_order int(10) unsigned NOT NULL,
       population  VARCHAR(50),
       match_regex  VARCHAR(1000),
       PRIMARY KEY (rule_order)
) ENGINE=MYISAM;

CREATE TABLE study_id (
       study_id  VARCHAR(50),
       PRIMARY KEY (study_id)
) ENGINE=MYISAM;


CREATE TABLE `verifybamid_readgroup` (
  `verifybamid_readgroup_id` int(10) NOT NULL AUTO_INCREMENT,
  `other_id` int(10) unsigned NOT NULL,
  `run_id` varchar(12) NOT NULL,
  `selfibd` decimal(4,2) DEFAULT NULL,
  `selfmix` decimal(4,2) DEFAULT NULL,
  `best_sample` varchar(12) DEFAULT NULL,
  `bestibd` decimal(4,2) DEFAULT NULL,
  `bestmix` decimal(4,2) DEFAULT NULL,
  `status` tinyint(1) DEFAULT '0',
  PRIMARY KEY (`verifybamid_readgroup_id`),
  UNIQUE KEY `other_id` (`other_id`,`run_id`)
) ENGINE=MYISAM;


CREATE TABLE `verifybamid_sample` (
  `verifybamid_sample_id` int(10) NOT NULL AUTO_INCREMENT,
  `other_id` int(10) unsigned NOT NULL,
  `table_name` enum('file','collection','run_meta_info','input_string') DEFAULT NULL,
  `sample_name` varchar(10) DEFAULT NULL,
  `selfibd` float(4,2) DEFAULT NULL,
  `selfibdllk` decimal(10,0) DEFAULT NULL,
  `selfibdllkdiff` float DEFAULT NULL,
  `het_a1` float DEFAULT NULL,
  `alt_a1` float DEFAULT NULL,
  `dp` float unsigned DEFAULT NULL,
  `mix` float(4,2) DEFAULT NULL,
  `hom` float DEFAULT NULL,
  `besthommixllk` float DEFAULT NULL,
  `besthommixllkdiff` float DEFAULT NULL,
  `num_run_ids` int(4) NOT NULL,
  `num_low_selfibd_run_ids` int(4) DEFAULT NULL,
  `sequence_index` varchar(10) DEFAULT NULL,
  `analysis_group` varchar(25) DEFAULT NULL,
  `chr20` tinyint(4) DEFAULT NULL,
  `failed` varchar(15) DEFAULT NULL,
  `status` tinyint(1) DEFAULT '0',
  `performed` datetime NOT NULL,
  PRIMARY KEY (`verifybamid_sample_id`),
  UNIQUE KEY `other_id` (`other_id`)
) ENGINE=MYISAM;


CREATE TABLE  verifybamid(
        verifybamid_id          INT  unsigned    NOT NULL AUTO_INCREMENT,
        file_id                 INT  unsigned    NOT NULL,
        sample                  VARCHAR (15)     NOT NULL,
        read_group              VARCHAR (15)     NOT NULL,
        chip_id                 VARCHAR (15),
        snps                    INT ,
        num_reads               INT ,
        avg_depth               float,
        free_contam             float,
        free_mlogl_est_contam   float ,
        free_mlogl_zero_contam  float,
        free_ref_bias_ref_het   float,
        free_ref_bias_refhomalt float,
        chip_contam             float,
        chip_mlogl_est_contam   float,
        chip_mlogl_zero_contam  float,
        chip_ref_bias_ref_het   float,
        chip_ref_bias_refhomalt float,
        depth_homref_site       float,
        rel_depth_het_site      float ,
        rel_depth_homalt_site   float,
        run_mode                VARCHAR(10),
        used_genotypes          INT,
        target_region           VARCHAR(15),
        vcf                     VARCHAR(50) NOT NULL ,
        verdict                 VARCHAR(20) ,
        performed               datetime    NOT NULL,

     PRIMARY KEY (verifybamid_id),
     UNIQUE(file_id, read_group)
) ENGINE=MYISAM;

CREATE TABLE pipeline(
       pipeline_id int(10) unsigned NOT NULL AUTO_INCREMENT,
       name VARCHAR(100) NOT NULL,
       table_name VARCHAR(50) NOT NULL,
       config_module VARCHAR(255) NOT NULL,
       config_options VARCHAR(30000),
       seeding_module VARCHAR(255) NOT NULL,
       seeding_options VARCHAR(30000),
       created   datetime NOT NULL,
       PRIMARY KEY(pipeline_id),
       UNIQUE(name)   
) ENGINE=MYISAM;

CREATE TABLE hive_db(
       hive_db_id int(10) unsigned NOT NULL AUTO_INCREMENT,
       pipeline_id int(10) unsigned NOT NULL,
       name VARCHAR(255) NOT NULL,
       host VARCHAR(255) NOT NULL,
       port smallint unsigned NOT NULL,
       created   datetime NOT NULL,
       retired   datetime,
       hive_version VARCHAR(255) NOT NULL,
       is_seeded tinyint NOT NULL DEFAULT 0,
       PRIMARY KEY(hive_db_id),
       UNIQUE(name,host,port,created)   
) ENGINE=MYISAM;

CREATE TABLE pipeline_seed(
       pipeline_seed_id int(10) unsigned NOT NULL AUTO_INCREMENT,
       seed_id int(10) unsigned NOT NULL,
       hive_db_id int(10) unsigned NOT NULL,
       is_running tinyint NOT NULL,
       is_complete tinyint NOT NULL default 0,
       is_failed tinyint NOT NULL default 0,
       is_futile tinyint NOT NULL default 0,
       created datetime NOT NULL,
       completed datetime,
       PRIMARY KEY(pipeline_seed_id)
) ENGINE=MYISAM;

CREATE TABLE pipeline_output(
       pipeline_output_id int(10) unsigned NOT NULL AUTO_INCREMENT,
       pipeline_seed_id int(10) unsigned NOT NULL,
       table_name VARCHAR(50) NOT NULL,
       output_id int(10) unsigned NOT NULL,
       action VARCHAR(50) NOT NULL,     
       PRIMARY KEY(pipeline_output_id)
) ENGINE=MYISAM;




#Now to add entries to the two standard tables

INSERT INTO archive_action (archive_action_id, action) values(1, 'archive');
INSERT INTO archive_action (archive_action_id, action) values(2, 'dearchive');
INSERT INTO archive_action (archive_action_id, action) values(3, 'replace');
INSERT INTO archive_action (archive_action_id, action) values(4, 'move_within_volume');
INSERT INTO archive_location (archive_location_id, location, location_name) values(1, '/nfs/1000g-work/G1K/archive_staging', 'staging');
INSERT INTO archive_location (archive_location_id, location, location_name) values(2, '/nfs/1000g-archive', 'archive');

