# longReads

`analyze_fastq_r` is an R function from the `longReads` package that provides a fast read-level summary of FASTQ files using Rust for performance.

## Features

Processes FASTQ files and outputs a dataframe with the following columns:

| Column        | Description                      |
|---------------|----------------------------------|
| `read_id`     | Identifier of the read           |
| `length`      | Length of the read in base pairs |
| `avg_quality` | Average Phred quality score      |
| `gc_content`  | GC content percentage            |

Built in Rust using [`bio::io::fastq`](https://docs.rs/bio/latest/bio/io/fastq/) for efficient parsing and exposed to R via [`extendr`](https://extendr.github.io/).

## Installation

1.  Make sure **Rust** is installed and in your `$PATH`:

``` bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

2.  Install the R package from local source or GitHub:

``` r
# From local source
devtools::install_local("path/to/analyze_fastq_r")

# Or from GitHub
devtools::install_github("username/analyze_fastq_r")
```

3.  Usage

``` r
library(analyze_fastq_r)

# Path to your FASTQ file
fastq_file <- "example.fastq"

# Generate read-level summary
summary_df <- analyze_fastq_r(fastq_file)

# View results
head(summary_df)
```