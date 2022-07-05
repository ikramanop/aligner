use aligner_core::enums::Protein;
use aligner_core::get_blosum62;
use aligner_core::simple::SimpleLocalAligner;
use aligner_core::statistics::calculate_p_value;
use aligner_web::repository::get_connection;

extern crate pretty_env_logger;
#[macro_use]
extern crate log;

fn main() {
    pretty_env_logger::init();

    let mut conn = get_connection(false);

    let seqs_vec = conn.get_cmp_sequence_ids_with_null_p_value().unwrap();

    let matrix = get_blosum62();

    for seqs in seqs_vec.iter() {
        debug!("Calculating p-value for task with id {}", seqs.0);

        let query = Protein::string_to_protein_vec(&match conn
            .get_sequence_by_identifier(seqs.1.clone())
        {
            Ok(query) => query,
            Err(_) => {
                debug!("Query sequence not found. Skipping...");
                continue;
            }
        })
        .unwrap();
        let target = Protein::string_to_protein_vec(&match conn
            .get_sequence_by_identifier(seqs.2.clone())
        {
            Ok(target) => target,
            Err(_) => {
                debug!("Target sequence not found. Skipping...");
                continue;
            }
        })
        .unwrap();

        let mut aligner = SimpleLocalAligner::from_protein_seqs(&query, &target).unwrap();

        let result = aligner
            .perform_alignment(11f64, 1f64, &matrix, None)
            .unwrap();

        let p_value =
            calculate_p_value(query, target, result.alignment.f, 11f64, 1f64, &matrix).unwrap();

        debug!("Calculated p-value {} for task with id {}", p_value, seqs.0);

        match conn.add_cmp_p_value_by_id(p_value, seqs.0) {
            Ok(_) => debug!("Successfully written p_value to DB"),
            Err(err) => debug!("{} happened, skipping", err),
        };
    }
}
