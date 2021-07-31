use rdkafka::config::ClientConfig;
use rdkafka::producer::{FutureProducer, FutureRecord};
use std::time::Duration;

#[tokio::main]
async fn main() {
    spawn_jobs().await
}

// A function to spawn jobs for alignment in Kafka
async fn spawn_jobs() {
    let producer: &FutureProducer = &ClientConfig::new()
        .set("bootstrap.servers", "localhost:9092")
        .set("message.timeout.ms", "1000")
        .create()
        .unwrap();

    producer
        .send(
            FutureRecord::to("my.topic").key("Bella").payload("Bobus"),
            Duration::from_secs(0),
        )
        .await
        .ok();
}
