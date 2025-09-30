use extendr_api::prelude::*;
use bio::io::fastq;
use bio::io::fastq::Record;
use bio::io;
use std::path::Path;


fn mean_qscore(record: &Record) -> f64 {
    if record.qual().is_empty() {
        return 0.0;
    }

    let errors: Vec<f64> = record
        .qual()
        .iter()
        .map(|&q| {
            let phred = (q as i32 - 33) as f64;
            10f64.powf(-phred / 10.0)
        })
        .collect();

    let mean_error = errors.iter().sum::<f64>() / errors.len() as f64;
    -10.0 * mean_error.log10()
}

/// Analyze a FASTQ file and return per-read statistics
///
/// # Arguments
/// * `file_path` - Path to the FASTQ file
/// * `min_length` - Optional minimum read length filter
/// * `min_avg_qual` - Optional minimum average Q-score filter
/// * `min_gc_content` - Optional minimum GC content filter
#[extendr]
fn analyze_fastq_r(
    file_path: &str,
    min_length: Option<usize>,
    min_avg_qual: Option<f64>,
    min_gc_content: Option<f64>,
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

            let avg_qual = mean_qscore(&record);

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
    fn analyze_fastq_r;
}
