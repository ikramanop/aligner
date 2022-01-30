use aligner::align::heuristic_alignment::HeuristicPairwiseAlignmentTool;
use aligner::database::get_connection;
use aligner::{align::enums::Protein, web::models::AlignJob};

use rdkafka::{consumer::CommitMode, message::Message};
use rdkafka::{
    consumer::{Consumer, StreamConsumer},
    ClientConfig,
};
use serde_json::from_slice;
use std::env;

extern crate pretty_env_logger;
#[macro_use]
extern crate log;

fn load_config(path: &'_ str) {
    debug!("Loading config from {:?}.", path);
    dotenv::from_filename(path).unwrap();
    debug!("Config successfully loaded.")
}

#[tokio::main]
async fn main() {
    pretty_env_logger::init();

    match env::var("CONFIG_PATH") {
        Ok(value) => load_config(&value),
        Err(err) => panic!("Unable to load config: {}", err),
    }

    let consumer: &StreamConsumer = &ClientConfig::new()
        .set("group.id", env::var("KAFKA_CONSUMER_GROUP").unwrap())
        .set("bootstrap.servers", env::var("KAFKA_HOST").unwrap())
        .set("enable.partition.eof", "false")
        .set("session.timeout.ms", "10000")
        .set("enable.auto.commit", "true")
        .set("auto.offset.reset", "latest")
        .create()
        .unwrap();

    consumer
        .subscribe(&[&env::var("KAFKA_TOPIC_JOBS").unwrap()])
        .unwrap();

    let mut conn = get_connection(false);

    loop {
        match consumer.recv().await {
            Ok(m) => {
                m.payload();

                debug!(
                    "Consuming {}:{}@{}: Subtask with uuid {}.",
                    m.topic(),
                    m.partition(),
                    m.offset(),
                    String::from_utf8(m.key().unwrap().to_vec()).unwrap(),
                );

                let job: AlignJob = from_slice(m.payload().unwrap()).unwrap();
                let mut tool = HeuristicPairwiseAlignmentTool::from_sequences_pair((
                    Protein::u8_vec_to_protein_vec(&job.sequence_1),
                    Protein::u8_vec_to_protein_vec(&job.sequence_2),
                ));

                let (current_f, optimal, sequences) = tool.local_alignment(
                    job.del_value,
                    job.del_value,
                    job.kd_value,
                    job.r_squared_value,
                    &job.matrix.unwrap(),
                    &job.frequences,
                );

                consumer.commit_message(&m, CommitMode::Async).unwrap();

                if conn
                    .insert_align_subtask(job.hash.clone(), current_f, &optimal, sequences.clone())
                    .is_err()
                {
                    conn = get_connection(false);
                    if let Err(err) =
                        conn.insert_align_subtask(job.hash, current_f, &optimal, sequences)
                    {
                        debug!("Failed with {}. Error with inserting subtask.", err)
                    }
                };
            }
            Err(err) => panic!("Kafka error: {}", err),
        };
    }
}
