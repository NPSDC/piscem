use anyhow::{anyhow, bail};
use clap::{arg, builder::ArgGroup, crate_authors, crate_version, value_parser, Command};
use itertools::Itertools;
use slog::{crit, o, warn, Drain};
use std::path::PathBuf;

use piscem_atac::cellfilter::{generate_permit_list, CellFilterMethod};
use piscem_atac::cmd_parse_utils::{
    pathbuf_directory_exists_validator, pathbuf_file_exists_validator,
};
use piscem_atac::prog_opts::GenPermitListOpts;
const VERSION: &str = env!("CARGO_PKG_VERSION");

fn main() -> anyhow::Result<()> {
    let num_hardware_threads = num_cpus::get() as u32;
    let max_num_threads: String = (num_cpus::get() as u32).to_string();

    let crate_authors = crate_authors!("\n");
    let version = crate_version!();
    let cmdline = std::env::args().join(" ");

    let gen_app = Command::new("generate-permit-list")
        .about("Generate a permit list of barcodes from a whitelist file")
        .version(version)
        .author(crate_authors)
        .arg(arg!(-i --input <INPUT>  "input directory containing the map.bed BED file")
        .required(true)
        .value_parser(pathbuf_directory_exists_validator))
    .arg(arg!(-o --"output-dir" <OUTPUTDIR>  "output directory").required(true).value_parser(value_parser!(PathBuf)))
    .arg(
        arg!(-u --"unfiltered-pl" <UNFILTEREDPL> "uses an unfiltered external permit list")
        .value_parser(pathbuf_file_exists_validator)
    )
    .group(ArgGroup::new("filter-method")
           .args(["knee-distance", "expect-cells", "force-cells", "valid-bc", "unfiltered-pl"])
           .required(true)
           )
    .arg(
        arg!(-m --"min-reads" <MINREADS> "minimum read count threshold; only used with --unfiltered-pl")
            .value_parser(value_parser!(usize))
            .default_value("10"));

    let opts = Command::new("alevin-fry")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .version(version)
        .author(crate_authors)
        .about("Process RAD files from the command line")
        .subcommand(gen_app)
        .get_matches();

    let decorator = slog_term::TermDecorator::new().build();
    let drain = slog_term::CompactFormat::new(decorator)
        .use_custom_timestamp(|out: &mut dyn std::io::Write| {
            write!(out, "{}", chrono::Local::now().format("%Y-%m-%d %H:%M:%S")).unwrap();
            Ok(())
        })
        .build()
        .fuse();
    let drain = slog_async::Async::new(drain).build().fuse();
    let log = slog::Logger::root(drain, o!());

    if let Some(t) = opts.subcommand_matches("generate-permit-list") {
        let input_dir: &PathBuf = t.get_one("input").expect("no input directory specified");
        let output_dir: &PathBuf = t
            .get_one("output-dir")
            .expect("no output directory specified");
        let mut fmeth = CellFilterMethod::KneeFinding;
        if let Some(v) = t.get_one::<PathBuf>("unfiltered-pl") {
            let min_reads: usize = *t
                .get_one("min-reads")
                .expect("min-reads must be a valid integer");
            if min_reads < 1 {
                crit!(
                    log,
                    "min-reads < 1 is not supported, the value {} was provided",
                    min_reads
                );
                std::process::exit(1);
            }
            fmeth = CellFilterMethod::UnfilteredExternalList(v.clone(), min_reads);
        };
        let gpl_opts = GenPermitListOpts::builder()
            .input_dir(input_dir)
            .output_dir(output_dir)
            .fmeth(fmeth)
            .version(VERSION)
            .cmdline(&cmdline)
            .log(&log)
            .build();

        match generate_permit_list(gpl_opts) {
            Ok(nc) if nc == 0 => {
                warn!(log, "found 0 corrected barcodes; please check the input.");
            }
            Err(e) => return Err(e),
            _ => (),
        };
    }

    Ok(())
}
