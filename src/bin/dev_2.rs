use aligner::files::get_db;
use aligner::matrices::get_population;

fn main() {
    let db = match get_db(&"database/matrices") {
        Ok(db) => db,
        Err(err) => panic!("{}", err),
    };

    let shape = (20, 20);
    let _matrices = get_population(&300, shape, &22, &db);

    println!("{:?}", _matrices[17]);

    match db.flush() {
        Ok(_) => {}
        Err(err) => {
            panic!("{}", err)
        }
    };
}
