use aligner::align::heuristic_alignment::HeuristicPairwiseAlignmentTool;
use aligner::{
    align::enums::Protein,
    web::models::{AlignJob, AlignJobResult},
};
use dotenv;
use rdkafka::{consumer::CommitMode, message::Message, producer::FutureRecord};
use rdkafka::{
    consumer::{Consumer, StreamConsumer},
    producer::FutureProducer,
    ClientConfig,
};
use serde_json::from_slice;
use std::{env, time::Duration};

fn load_config<'a>(path: &'a str) {
    println!("{:?}", path);
    dotenv::from_filename(path).unwrap();
}

#[tokio::main]
async fn main() {
    match env::var("CONFIG_PATH") {
        Ok(value) => load_config(&value),
        Err(err) => panic!("Unable to load config: {}", err),
    }

    println!("{:?}", env::var("SERVER_ENV").unwrap());

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

    let producer: &FutureProducer = &ClientConfig::new()
        .set("bootstrap.servers", env::var("KAFKA_HOST").unwrap())
        .set("message.timeout.ms", "1000")
        .create()
        .unwrap();

    loop {
        match consumer.recv().await {
            Ok(m) => {
                m.payload();

                println!(
                    "{}:{}@{}: Matrix {}",
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
                    job.kd_value,
                    job.r_squared_value,
                    &job.matrix.unwrap(),
                    &job.frequences,
                );

                consumer.commit_message(&m, CommitMode::Async).unwrap();

                println!("{:?}", sequences);

                producer
                    .send(
                        FutureRecord::to(&env::var("KAFKA_TOPIC_RESULTS").unwrap())
                            .key(&String::from_utf8(m.key().unwrap().to_vec()).unwrap())
                            .payload(
                                &serde_json::to_string(&AlignJobResult {
                                    matrix: optimal,
                                    max_f: current_f,
                                    matrices_volume_value: job.matrices_volume_value,
                                    uuid: job.uuid,
                                })
                                .unwrap(),
                            ),
                        Duration::from_secs(0),
                    )
                    .await
                    .ok();
            }
            Err(err) => panic!("Kafka error: {}", err),
        };
    }
}
