use clap::Parser;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
pub(crate) struct Args {
    #[clap(short, long)]
    pub(crate) input: String,

    #[clap(short, long, default_value_t = 11f64)]
    pub(crate) deletions: f64,

    #[clap(short, long, default_value_t = 2f64)]
    pub(crate) extension: f64,

    #[clap(short, long)]
    pub(crate) global: bool,

    #[clap(short, long, default_value_t = String::from("out/result.txt"))]
    pub(crate) output: String,
}
