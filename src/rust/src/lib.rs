use extendr_api::prelude::*;
use bio::io::fastq;
use std::path::Path;


/// Return string `"Hello world!"` to R.
/// @export
#[extendr]
fn hello_world() -> &'static str {
    "Hello world!"
}
#[extendr]
fn add(x: i32, y: i32) -> i32 {
    x + y
}
#[extendr]
fn analyze_fastq(file_path: &str, min_length: Option<usize>, min_avg_qual: Option<f32>) -> Robj {
    let reader = fastq::Reader::from_file(Path::new(file_path));
    if reader.is_err() {
        return list!(total_reads = 0, filtered_reads = 0).into_robj();
    }
    let reader = reader.unwrap();

    let mut total_reads = 0;
    let mut filtered_reads = 0;

    for result in reader.records() {
        if let Ok(record) = result {
            total_reads += 1;

            let seq_len = record.seq().len();
            let avg_qual = if !record.qual().is_empty() {
                let sum: u32 = record.qual().iter().map(|&q| q as u32).sum();
                sum as f32 / record.qual().len() as f32
            } else {
                0.0
            };

            let length_ok = min_length.map_or(true, |min| seq_len >= min);
            let qual_ok = min_avg_qual.map_or(true, |min| avg_qual >= min);

            if length_ok && qual_ok {
                filtered_reads += 1;
            }
        }
    }

    list!(total_reads, filtered_reads).into_robj()
}


#[extendr]
fn analyze_fastq_r(
    file_path: &str,
    min_length: Option<usize>,
    min_avg_qual: Option<f64>,
    min_gc_content: Option<f64>
) -> Robj {
    let reader = fastq::Reader::from_file(Path::new(file_path));
    if reader.is_err() {
        return list!(
            id = Vec::<String>::new(),
            length = Vec::<i32>::new(),
            avg_quality = Vec::<f64>::new(),
            gc_content = Vec::<f64>::new()
        ).into_robj();
    }
    let reader = reader.unwrap();

    let mut ids = Vec::new();
    let mut lengths = Vec::new();
    let mut avg_quals = Vec::new();
    let mut gc_contents = Vec::new();

    for result in reader.records() {
        if let Ok(record) = result {
            let seq = record.seq();
            let seq_len = seq.len();
            if seq_len == 0 { continue; }

            let avg_qual = if !record.qual().is_empty() {
                let sum: u32 = record.qual().iter().map(|&q| q as u32).sum();
                sum as f64 / record.qual().len() as f64
            } else { 0.0 };

            let gc_count = seq.iter().filter(|&&b| b == b'G' || b == b'g' || b == b'C' || b == b'c').count();
            let gc_content = (gc_count as f64 / seq_len as f64) * 100.0;

            let length_ok = min_length.map_or(true, |min| seq_len >= min);
            let qual_ok = min_avg_qual.map_or(true, |min| avg_qual >= min);
            let gc_ok = min_gc_content.map_or(true, |min| gc_content >= min);

            if length_ok && qual_ok && gc_ok {
                ids.push(record.id().to_string());
                lengths.push(seq_len as i32);
                avg_quals.push(avg_qual);
                gc_contents.push(gc_content);
            }
        }
    }

    list!(
        id = ids,
        length = lengths,
        avg_quality = avg_quals,
        gc_content = gc_contents
    ).into_robj()
}



// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod longReads;
    fn hello_world;
    fn add;
    fn analyze_fastq;
    fn analyze_fastq_r;
}
