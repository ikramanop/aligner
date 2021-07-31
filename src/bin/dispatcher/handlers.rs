use crate::models::{
    AlignJob, AlignJobRequest, AlignJobResult, AlignmentResultEventResponse,
    EmptySuccessfulResponse, ErroneousResponse, HealthCheck, ProgressEventResponse,
};
use aligner::{files::get_db, matrices::get_population};
use futures::{Stream, StreamExt};
use ndarray::{arr1, arr2, Array2};
use rdkafka::{
    consumer::{BaseConsumer, CommitMode, Consumer},
    message::Message,
    producer::{BaseProducer, BaseRecord},
    ClientConfig,
};
use scopeguard::defer;
use seq_io::fasta::Reader;
use std::env;
use std::{collections::HashMap, thread};
use std::{convert::Infallible, str::FromStr};
use tokio::sync::mpsc;
use tokio_stream::wrappers::UnboundedReceiverStream;
use uuid::Uuid;
use warp::{
    http::StatusCode,
    sse::{reply, Event},
};

// GET health/check handler
pub async fn check_health() -> Result<impl warp::Reply, Infallible> {
    let status = HealthCheck { status: true };

    Ok(warp::reply::json(&status))
}

// POST /validate handler
pub async fn validate(job: AlignJobRequest) -> Result<impl warp::Reply, Infallible> {
    println!("received body\n{:?}", job);

    let mut reader = Reader::new(job.sequences.as_bytes());

    let seqs: Vec<_> = match reader.records().collect() {
        Ok(seqs) => seqs,
        Err(_) => {
            return Ok(warp::reply::with_status(
                warp::reply::json(&ErroneousResponse {
                    message: "Ошибка при обработке данных".to_owned(),
                }),
                StatusCode::BAD_REQUEST,
            ))
        }
    };

    if seqs.len() != 2 {
        log::debug!("There's should be 2 sequences, not {}", seqs.len());
        return Ok(warp::reply::with_status(
            warp::reply::json(&ErroneousResponse {
                message: "Передано отличное от двух число последовательностей".to_owned(),
            }),
            StatusCode::BAD_REQUEST,
        ));
    }

    defer!(spawn_jobs((seqs[0].seq.clone(), seqs[1].seq.clone()), &job));

    Ok(warp::reply::with_status(
        warp::reply::json(&EmptySuccessfulResponse {}),
        StatusCode::OK,
    ))
}

// GET /progress handler
pub async fn progress(uuid: Uuid) -> Result<impl warp::Reply, Infallible> {
    let stream = consume_results(uuid);

    Ok(reply(warp::sse::keep_alive().stream(stream)))
}

fn consume_results(uuid: Uuid) -> impl Stream<Item = Result<Event, warp::Error>> + Send + 'static {
    let (tx, rx) = mpsc::unbounded_channel();
    let rx = UnboundedReceiverStream::new(rx);

    thread::spawn(move || {
        let consumer: BaseConsumer = ClientConfig::new()
            .set("group.id", env::var("KAFKA_CONSUMER_GROUP").unwrap())
            .set("bootstrap.servers", env::var("KAFKA_HOST").unwrap())
            .set("enable.partition.eof", "false")
            .set("session.timeout.ms", "10000")
            .set("enable.auto.commit", "true")
            .set("auto.offset.reset", "earliest")
            .create()
            .unwrap();

        consumer
            .subscribe(&[&env::var("KAFKA_TOPIC_RESULTS").unwrap()])
            .unwrap();

        let mut counter = 0;
        let mut progress = 0f64;

        let mut best_matrix: Array2<f64> = arr2(&vec![[0f64]]);
        let mut max_f = 0f64;

        loop {
            match consumer.poll(None) {
                Some(m) => {
                    let message = m.unwrap();

                    if Uuid::from_str(&String::from_utf8(message.key().unwrap().to_vec()).unwrap())
                        .unwrap()
                        == uuid
                    {
                        println!(
                            "{}@{}: Receiving data for job {}",
                            message.partition(),
                            message.offset(),
                            String::from_utf8(message.key().unwrap().to_vec()).unwrap(),
                        );

                        let result: AlignJobResult =
                            serde_json::from_slice(message.payload().unwrap()).unwrap();

                        if result.max_f > max_f {
                            max_f = result.max_f;
                            best_matrix = result.matrix;
                        }
                        counter += 1;
                        progress = counter as f64 / result.matrices_volume_value as f64;

                        consumer.commit_message(&message, CommitMode::Sync).unwrap();

                        if progress == 1f64 {
                            &tx.send(
                                serde_json::to_string(&AlignmentResultEventResponse {
                                    progress: progress * 100f64,
                                    matrix: best_matrix,
                                    max_f,
                                })
                                .unwrap(),
                            )
                            .unwrap();

                            break;
                        }

                        &tx.send(
                            serde_json::to_string(&ProgressEventResponse {
                                progress: progress * 100f64,
                                message: "Выравнивание вычисляется".to_owned(),
                            })
                            .unwrap(),
                        )
                        .unwrap();
                    }
                }
                None => tx
                    .send(
                        serde_json::to_string(&ProgressEventResponse {
                            progress: progress * 100f64,
                            message: "Выравнивание вычисляется".to_owned(),
                        })
                        .unwrap(),
                    )
                    .unwrap(),
            }
        }
    });

    rx.map(|msg| Ok(Event::default().data(msg)))
}

// A function to spawn jobs for alignment in Kafka
fn spawn_jobs(sequences: (Vec<u8>, Vec<u8>), job: &AlignJobRequest) {
    let producer: &BaseProducer = &ClientConfig::new()
        .set("bootstrap.servers", env::var("KAFKA_HOST").unwrap())
        .set("message.timeout.ms", "1000")
        .create()
        .unwrap();

    let db = &get_db(&"database/matrices").unwrap();

    let matrices = get_population(
        &job.matrices_volume_value,
        (job.dim_value as usize, job.dim_value as usize),
        db,
    );

    let mut counters = HashMap::<u8, i32>::new();
    let mut freqs = Vec::<f64>::new();

    for symbol in sequences.0.iter() {
        let counter = counters.entry(*symbol).or_insert(0);
        *counter += 1;
    }
    for symbol in sequences.1.iter() {
        let counter = counters.entry(*symbol).or_insert(0);
        *counter += 1;
    }
    for (_, counter) in counters.iter() {
        freqs.push(*counter as f64 / (sequences.0.len() + sequences.1.len()) as f64);
    }

    let mut data = AlignJob {
        sequence_1: sequences.0,
        sequence_2: sequences.1,
        matrix: None,
        frequences: arr1(&freqs[..]),
        kd_value: job.kd_value,
        r_squared_value: job.r_squared_value,
        del_value: job.del_value,
        matrices_volume_value: job.matrices_volume_value,
    };

    for matrix in matrices {
        data.matrix = Some(matrix);

        producer
            .send(
                BaseRecord::to(&env::var("KAFKA_TOPIC_JOBS").unwrap())
                    .key(&job.uuid.to_string())
                    .payload(&serde_json::to_string(&data).unwrap()),
            )
            .unwrap();

        producer.poll(None);
    }
}
