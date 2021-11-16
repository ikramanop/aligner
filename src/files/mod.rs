use csv::ReaderBuilder;
use ndarray::Array2;
use ndarray::Axis;
use std::fs::File;
use std::fs::OpenOptions;
use std::io::prelude::*;
use std::path::{Display, Path};

struct FileStruct<'a>(&'a Path, Display<'a>, File);

fn read_file(path: &Path, write: bool) -> FileStruct {
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

pub fn write_to_file(path: &Path, contents: &str) {
    let FileStruct(_, _, mut file) = read_file(path, true);

    file.write_all(contents.as_bytes()).unwrap();
}

pub fn load_file_contents(path: &Path) -> String {
    let FileStruct(_, _, mut file) = read_file(path, false);

    let mut contents = String::new();
    file.read_to_string(&mut contents).unwrap();

    contents
}

pub fn convert_csv_to_matrix(csv_matrix: &[u8], dim: (usize, usize)) -> Array2<f64> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b' ')
        .has_headers(false)
        .from_reader(csv_matrix);

    let mut result_matrix = Array2::<f64>::zeros(dim);
    for (i, record) in reader.records().enumerate() {
        let matrix_row = match record {
            Ok(row) => row,
            Err(_) => panic!("error reading matrix"),
        };
        for (j, elem) in matrix_row.iter().enumerate() {
            result_matrix[[i, j]] = elem.parse::<f64>().unwrap();
        }
    }

    result_matrix
}

pub fn convert_matrix_to_csv(matrix: &Array2<f64>) -> Vec<u8> {
    let mut result: String = String::from("");

    for (i, elem) in matrix.iter().enumerate() {
        result += &format!("{}", elem);

        if (i + 1) % matrix.len_of(Axis(1)) == 0 {
            result += "\n";
        } else {
            result += " ";
        }
    }

    result.as_bytes().to_vec()
}

pub fn get_db(path: &str) -> Result<sled::Db, &'static str> {
    let tree = match sled::open(path) {
        Ok(tree) => tree,
        Err(_) => return Err("Something wrong happened with db"),
    };

    match tree.get("BLOSUM50") {
        Ok(None) => {
            let raw_data = load_file_contents(Path::new("static/BLOSUM50.csv"));

            match tree.insert("BLOSUM50", raw_data.as_bytes()) {
                Ok(_) => {}
                Err(err) => {
                    panic!("{}", err)
                }
            };
        }
        Err(_) => return Err("Something wrong happened with db"),
        _ => {}
    };

    match tree.get("PAM250") {
        Ok(None) => {
            let raw_data = load_file_contents(Path::new("static/PAM250.csv"));

            match tree.insert("PAM250", raw_data.as_bytes()) {
                Ok(_) => {}
                Err(err) => {
                    panic!("{}", err)
                }
            };
        }
        Err(_) => return Err("Something wrong happend with db"),
        _ => {}
    };

    Ok(tree)
}
