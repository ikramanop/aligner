use aligner::database::get_connection;
use aligner::matrices::get_population;
use aligner::web::models::{
    AlignJob, AlignJobRequest, EmptySuccessfulResponseWithHashes, ErroneousResponse, HealthCheck,
    HealthCheckUnit, ProgressEventResponse,
};
use futures::StreamExt;
use ndarray::arr1;
use rdkafka::{
    client::Client,
    groups::GroupList,
    producer::{BaseProducer, BaseRecord},
    types::RDKafkaType,
    ClientConfig, ClientContext,
};
use seq_io::fasta::{OwnedRecord, Reader};
use std::collections::HashMap;
use std::convert::Infallible;
use std::env;
use std::time::Duration;
use tokio::time::interval;
use tokio_stream::wrappers::IntervalStream;
use uuid::Uuid;
use warp::http::StatusCode;
use warp::sse::Event;

#[derive(Clone, Default)]
pub struct DefaultClientContext;
impl ClientContext for DefaultClientContext {}

// GET health/check handler
pub async fn check_health() -> Result<impl warp::Reply, Infallible> {
    let mut nodes = Vec::<HealthCheckUnit>::new();

    let mut config = ClientConfig::new();
    config.set("bootstrap.servers", env::var("KAFKA_HOST").unwrap());

    let native_config = config.create_native_config().unwrap();

    let client = match Client::new(
        &config,
        native_config,
        RDKafkaType::RD_KAFKA_CONSUMER,
        DefaultClientContext,
    ) {
        Ok(value) => value,
        Err(err) => panic!("{}", err),
    };

    let group_list: GroupList;

    loop {
        match client.fetch_group_list(Some("aligner.jobs.group"), None) {
            Ok(value) => {
                group_list = value;
                break;
            }
            Err(_) => continue,
        }
    }

    for member in group_list.groups()[0].members() {
        nodes.push(HealthCheckUnit {
            consumer_name: member.id().to_owned(),
            status: true,
        })
    }

    Ok(warp::reply::json(&HealthCheck { nodes }))
}

// POST /validate handler
pub async fn validate(job: AlignJobRequest) -> Result<impl warp::Reply, Infallible> {
    debug!("Received body {:?}.", job);

    let mut reader = Reader::new(job.sequences.as_bytes());

    let seqs: Vec<_> = match reader.records().collect() {
        Ok(seqs) => seqs,
        Err(err) => {
            debug!("Failed with {}. Some error with FASTA parsing.", err);

            return Ok(warp::reply::with_status(
                warp::reply::json(&ErroneousResponse {
                    message: "Ошибка при обработке данных. Проверьте входные последовательности."
                        .to_owned(),
                }),
                StatusCode::BAD_REQUEST,
            ));
        }
    };

    if seqs.len() < 2 {
        debug!("There should be 2+ sequences, not {}.", seqs.len());

        return Ok(warp::reply::with_status(
            warp::reply::json(&ErroneousResponse {
                message: "Передано отличное от число последовательностей меньше двух.".to_owned(),
            }),
            StatusCode::BAD_REQUEST,
        ));
    }

    let sequence_pairs = generate_pairs(&seqs);
    let mut hashes = Vec::<String>::new();

    let mut conn = get_connection(false);

    for (query, target) in sequence_pairs.iter() {
        let hash = match conn.insert_align_task(
            &job,
            String::from_utf8(query.head.clone()).unwrap(),
            String::from_utf8(query.seq.clone()).unwrap(),
            String::from_utf8(target.head.clone()).unwrap(),
            String::from_utf8(target.seq.clone()).unwrap(),
        ) {
            Ok(hash) => hash,
            Err(err) => {
                debug!("Failed with {}. Skipping for now.", err);
                continue;
            }
        };

        hashes.push(hash.clone());

        spawn_jobs(
            (seqs[0].seq.clone(), seqs[1].seq.clone()),
            &job,
            hash.clone(),
        );
    }

    if hashes.is_empty() {
        return Ok(warp::reply::with_status(
            warp::reply::json(&ErroneousResponse {
                message:
                    "Ошибка при создании запроса. Не было передано новых задач на выравнивание."
                        .to_owned(),
            }),
            StatusCode::INTERNAL_SERVER_ERROR,
        ));
    }

    Ok(warp::reply::with_status(
        warp::reply::json(&EmptySuccessfulResponseWithHashes { hashes }),
        StatusCode::OK,
    ))
}

// GET /progress?hashes=[] handler
pub async fn progress(hashes: Vec<String>) -> Result<impl warp::Reply, Infallible> {
    let interval = interval(Duration::from_secs(1));
    let stream = IntervalStream::new(interval);

    let mut conn = get_connection(false);

    let event_stream = stream.map(move |_| {
        let mut percentages = HashMap::<String, f64>::new();

        for hash in hashes.iter() {
            let percentage = conn.get_percentage_by_hash(hash.clone()).unwrap();

            debug!("Calculated {}% for task with hash {}.", percentage, hash);

            if (percentage - 100f64).abs() < f64::EPSILON {
                let subtask = conn
                    .get_subtask_with_max_f_value_by_hash(hash.clone())
                    .unwrap();
                conn.insert_result_matrix_by_hash(subtask, hash.clone())
                    .unwrap();
                conn.delete_subtasks_by_hash(hash.clone()).unwrap();
            }

            percentages.insert(hash.clone(), percentage).unwrap();
        }

        sse_event(percentages)
    });

    Ok(warp::sse::reply(event_stream))
}

fn sse_event(percentages: HashMap<String, f64>) -> Result<Event, Infallible> {
    Ok(Event::default()
        .json_data(&ProgressEventResponse {
            progress: percentages,
            message: "Выравнивание вычисляется".to_owned(),
        })
        .unwrap())
}

// A function to spawn jobs for alignment in Kafka
fn spawn_jobs(sequences: (Vec<u8>, Vec<u8>), job: &AlignJobRequest, hash: String) {
    let producer: &BaseProducer = &ClientConfig::new()
        .set("bootstrap.servers", env::var("KAFKA_HOST").unwrap())
        .set("message.timeout.ms", "1000")
        .create()
        .unwrap();

    let mut conn = get_connection(false);

    let matrices = get_population(
        &mut conn,
        job.dim_value as usize,
        job.matrices_volume_value as usize,
    )
    .unwrap();

    let mut counters = HashMap::<u8, i32>::new();
    let mut freqs = vec![0f64; job.dim_value as usize];

    for symbol in sequences.0.iter() {
        let counter = counters.entry(*symbol).or_insert(0);
        *counter += 1;
    }
    for symbol in sequences.1.iter() {
        let counter = counters.entry(*symbol).or_insert(0);
        *counter += 1;
    }
    for (i, (_, counter)) in counters.iter().enumerate() {
        freqs[i] = *counter as f64 / (sequences.0.len() + sequences.1.len()) as f64;
    }

    debug!("{:?}", freqs.len());

    let mut data = AlignJob {
        sequence_1: sequences.0,
        sequence_2: sequences.1,
        matrix: None,
        frequences: arr1(&freqs[..]),
        kd_value: job.kd_value,
        r_squared_value: job.r_squared_value,
        del_value: job.del_value,
        matrices_volume_value: job.matrices_volume_value,
        hash,
    };

    for matrix in matrices {
        data.matrix = Some(matrix);

        producer
            .send(
                BaseRecord::to(&env::var("KAFKA_TOPIC_JOBS").unwrap())
                    .key(&Uuid::new_v4().to_string())
                    .payload(&serde_json::to_string(&data).unwrap()),
            )
            .unwrap();

        producer.poll(None);
    }
}

fn generate_pairs(sequences: &[OwnedRecord]) -> Vec<(OwnedRecord, OwnedRecord)> {
    let mut result = Vec::<(OwnedRecord, OwnedRecord)>::new();

    for (i, rec) in sequences.iter().enumerate() {
        for (j, ord) in sequences[i + 1..].iter().enumerate() {
            result.push((rec.clone(), ord.clone()));
            debug!("Got {} {}", i, j);
        }
    }

    result
}
