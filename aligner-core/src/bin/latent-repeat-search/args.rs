use clap::Parser;

#[derive(Parser, Debug, Clone)]
#[clap(author, version, about, long_about = None)]
pub(crate) struct Args {
    #[clap(short, long)]
    pub(crate) input: Option<String>,

    #[clap(short, long)]
    pub(crate) output: Option<String>,

    #[clap(long)]
    pub(crate) csv: Option<String>,

    #[clap(short, long, default_value_t = 30f64)]
    pub(crate) deletions: f64,

    #[clap(short, long, default_value_t = 7f64)]
    pub(crate) extension: f64,

    #[clap(long, default_value_t = 100000f64)]
    pub(crate) rsquared: f64,

    #[clap(long, default_value_t = 0f64)]
    pub(crate) kd: f64,

    #[clap(short, long, default_value_t = 30)]
    pub(crate) query_offset: usize,

    #[clap(short, long, default_value_t = 300)]
    pub(crate) repeat_length: usize,

    #[clap(long, default_value_t = 1)]
    pub(crate) threads: usize,

    #[clap(long)]
    pub(crate) simple_init: bool,

    #[clap(long, default_value_t = 10)]
    pub(crate) repeats: usize,

    #[clap(long)]
    pub(crate) reverse: bool,
}
