use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::Parser;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::tpool::ThreadPool;

use v8bam::JsBamFilterEngine;

#[derive(Parser, Debug)]
struct Args {
    /// Input BAM ("-" for stdin)
    input: PathBuf,

    /// Output BAM ("-" for stdout)
    #[arg(short = 'o', long)]
    output: PathBuf,

    /// JS expression/body, e.g.:
    ///   'aln.mapq > 10 && aln.qname.startsWith("q23")'
    /// or:
    ///   'return aln.mapq > 10 && hasFlag(aln.flag, 0x2);'
    #[arg(short = 'e', long)]
    expr: String,

    /// Number of threads for BAM I/O
    #[arg(short = 't', long, default_value = "3")]
    threads: u32,
}

fn main() -> Result<()> {
    let args = Args::parse();

    // Create shared threadpool for BAM I/O
    let tpool = ThreadPool::new(args.threads)?;

    // Open BAM reader & writer
    let mut reader = if args.input.to_string_lossy() == "-" {
        bam::Reader::from_stdin().context("failed to open stdin as BAM")?
    } else {
        bam::Reader::from_path(&args.input)
            .with_context(|| format!("failed to open BAM {}", args.input.display()))?
    };
    reader.set_thread_pool(&tpool)?;

    let header = bam::Header::from_template(reader.header());
    let mut writer = if args.output.to_string_lossy() == "-" {
        bam::Writer::from_stdout(&header, bam::Format::Bam).context("failed to open BAM writer")?
    } else {
        bam::Writer::from_path(&args.output, &header, bam::Format::Bam)
            .context("failed to open BAM writer")?
    };
    writer.set_thread_pool(&tpool)?;
    let header_view = reader.header().clone();

    // Create JS filter engine
    let mut engine = JsBamFilterEngine::new(&args.expr)?;

    // Reuse record buffer
    let mut record = bam::Record::new();

    loop {
        match reader.read(&mut record) {
            Some(Ok(())) => {
                if engine.record_passes(&record, &header_view)? {
                    writer.write(&record)?;
                }
            }
            Some(Err(e)) => return Err(e.into()),
            None => break,
        }
    }

    Ok(())
}
