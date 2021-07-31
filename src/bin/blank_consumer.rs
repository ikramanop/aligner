use rdkafka::{consumer::Consumer, ClientConfig};
use rdkafka::{
    consumer::{BaseConsumer, CommitMode},
    message::Message,
};

#[tokio::main]
async fn main() {
    let consumer: &BaseConsumer = &ClientConfig::new()
        .set("group.id", "aligner.results.group")
        .set("bootstrap.servers", "localhost:9092")
        .set("enable.partition.eof", "false")
        .set("session.timeout.ms", "10000")
        .set("enable.auto.commit", "true")
        .set("auto.offset.reset", "earliest")
        .create()
        .unwrap();

    consumer.subscribe(&["aligner.results"]).unwrap();

    loop {
        match consumer.poll(None) {
            Some(m) => match m {
                Ok(message) => {
                    let value = match message.payload_view::<str>() {
                        Some(Ok(value)) => value,
                        Some(Err(e)) => {
                            panic!("Error while deserializing message payload: {:?}", e);
                        }
                        None => "",
                    };

                    println!(
                        "key: '{}', payload: '{}', topic: {}, partition: {}",
                        String::from_utf8(message.key().unwrap().to_vec()).unwrap(),
                        value,
                        message.topic(),
                        message.partition()
                    );

                    consumer.commit_message(&message, CommitMode::Sync).unwrap();
                }
                Err(err) => panic!("Kafka error: {}", err),
            },
            None => println!("No messages fam"),
        }
    }
}
