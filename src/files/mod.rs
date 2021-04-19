use csv::ReaderBuilder;
use ndarray::Array2;
use std::fs::File;
use std::fs::OpenOptions;
use std::io::prelude::*;
use std::path::{Display, Path};

struct FileStruct<'a>(&'a Path, Display<'a>, File);

fn read_file(path: &str, write: bool) -> FileStruct {
    let path = Path::new(path);
    let display = path.display();

    let mut file_options = OpenOptions::new();
    file_options.write(write).read(!write).truncate(write);

    if !Path::new(path).exists() && write {
        file_options.create_new(true);
    }

    let file = match file_options.open(&path) {
        Err(why) => panic!("couldn't open {}: {}", display, why),
        Ok(file) => file,
    };

    FileStruct(path, display, file)
}

pub fn write_to_file(path: &str, contents: &str) {
    let FileStruct(_, display, mut file) = read_file(path, true);

    match file.write_all(contents.as_bytes()) {
        Err(why) => panic!("error writing to {}: {}", display, why),
        _ => {}
    }
}

pub fn load_file_contents(path: &str) -> String {
    let FileStruct(_, display, mut file) = read_file(path, false);

    let mut contents = String::new();
    match file.read_to_string(&mut contents) {
        Err(why) => panic!("couldn't read {}: {}", display, why),
        Ok(_) => {}
    };

    contents
}

pub fn load_pam250() -> Array2<i32> {
    let raw_data = load_file_contents("static/PAM250.csv");

    let mut reader = ReaderBuilder::new()
        .delimiter(b' ')
        .has_headers(false)
        .from_reader(raw_data.as_bytes());

    let mut result_matrix = Array2::<i32>::zeros((25, 25));
    for (i, record) in reader.records().enumerate() {
        let matrix_row = match record {
            Ok(row) => row,
            Err(_) => panic!("error reading PAM250 matrix"),
        };
        for (j, elem) in matrix_row.iter().enumerate() {
            result_matrix[[i, j]] = elem.parse::<i32>().unwrap();
        }
    }

    result_matrix
}

pub fn load_blosum50() -> Array2<i32> {
    let raw_data = load_file_contents("static/BLOSUM50.csv");

    let mut reader = ReaderBuilder::new()
        .delimiter(b' ')
        .has_headers(false)
        .from_reader(raw_data.as_bytes());

    let mut result_matrix = Array2::<i32>::zeros((25, 25));
    for (i, record) in reader.records().enumerate() {
        let matrix_row = match record {
            Ok(row) => row,
            Err(_) => panic!("error reading BLOSUM50 matrix"),
        };
        for (j, elem) in matrix_row.iter().enumerate() {
            result_matrix[[i, j]] = elem.parse::<i32>().unwrap();
        }
    }

    result_matrix
}
