use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::path::Path;

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Record {
    pub name: String,
    pub z_value: f64,
    pub left_coord: usize,
    pub right_coord: usize,
}

pub fn read_csv(path: &Path) -> Result<HashMap<String, Vec<Record>>, impl Error> {
    let mut rdr = match csv::Reader::from_path(path) {
        Ok(rdr) => rdr,
        Err(err) => return Err(err),
    };

    let mut map: HashMap<String, Vec<Record>> = HashMap::new();
    for result in rdr.deserialize() {
        let record: Record = match result {
            Ok(record) => record,
            Err(err) => return Err(err),
        };

        if !map.contains_key(&record.name) {
            map.insert(record.name.clone(), vec![]);
        }

        let vec = map.get_mut(&record.name).unwrap();
        vec.push(record);
    }

    Ok(map)
}

pub struct CsvInput {
    wtr: csv::Writer<File>,
}

impl CsvInput {
    pub fn new(path: &Path) -> Result<CsvInput, impl Error> {
        let wtr = match csv::Writer::from_path(path) {
            Ok(wtr) => wtr,
            Err(err) => return Err(err),
        };

        Ok(CsvInput { wtr })
    }

    pub fn write(&mut self, record: &Record) -> Result<(), impl Error> {
        self.wtr.serialize(record)
    }
}
