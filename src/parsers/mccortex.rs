use nom::{IResult};
use nom::number::complete::{le_u8, le_u32, le_u64};
use nom::multi::{count, length_count};
use nom::bytes::complete::{tag, take_until};
use nom::sequence::tuple;

// Spec: https://github.com/mcveanlab/mccortex/blob/master/docs/file_formats/graph_file_format.txt
const SIGNATURE: &str = "CORTEX";

#[derive(Debug, PartialEq)]
pub struct McCortex {
    pub header: Header,
    pub sample_names: Vec<SampleName>,
    pub cleaning: Vec<Cleaning>
}

#[derive(Debug, PartialEq)]
pub struct Header {
    pub version: u32,
    pub kmer_size: u32,
    pub W: u32,
    pub cols: u32,
    pub mean_read_lens: Vec<u32>,
    pub total_seq_loaded: Vec<u64>,
}

#[derive(Debug, PartialEq)]
pub struct SampleName {
    pub len: u32,
    pub value: String,
}

#[derive(Debug, PartialEq)]
pub struct Cleaning {
    pub top_clip: bool,
    pub remove_low_covg_supernodes: bool,
    pub remove_low_covg_kmers: bool,
    pub cleaned_against_graph: bool,
    pub remove_low_coverage_supernodes_threshold: u32,
    pub remove_low_coverage_kmer_threshold: u32,
    pub graph_name: String,
}

pub fn header(input: &[u8]) -> IResult<&[u8], Header> {
    let (remain, _) = tag(SIGNATURE)(input)?;
    let (remain, (version, kmer_size, W, cols)) = tuple((le_u32, le_u32, le_u32, le_u32))(remain)?;
    let (remain, mean_read_lens) = count(le_u32, cols as usize)(remain)?;
    let (remain, total_seq_loaded) = count(le_u64, cols as usize)(remain)?;

    Ok((remain, Header {
        version,
        kmer_size,
        W,
        cols,
        mean_read_lens,
        total_seq_loaded,
    }))
}

pub fn sample_name(input: &[u8]) -> IResult<&[u8], SampleName> {
    let (remain, len) = le_u32(input)?;
    let (remain, value) = count(le_u8, len as usize)(remain)?;

    Ok((remain, SampleName {
        len,
        value: String::from_utf8(value).unwrap(),
    }))
}

pub fn cleaning(input: &[u8]) -> IResult<&[u8], Cleaning> {
    let (remain, top_clip) = le_u8(input)?;
    let (remain, remove_low_covg_supernodes) = le_u8(input)?;
    let (remain, remove_low_covg_kmers) = le_u8(input)?;
    let (remain, cleaned_against_graph) = le_u8(input)?;
    let (remain, remove_low_coverage_supernodes_threshold) = le_u32(input)?;
    let (remain, remove_low_coverage_kmer_threshold) = le_u32(input)?;
    let (remain, graph_name) = length_count(le_u32, le_u8)(remain)?;

    Ok((remain, Cleaning {
        top_clip: top_clip != 0,
        remove_low_covg_supernodes: remove_low_covg_supernodes != 0,
        remove_low_covg_kmers: remove_low_covg_kmers != 0,
        cleaned_against_graph: cleaned_against_graph != 0,
        remove_low_coverage_supernodes_threshold,
        remove_low_coverage_kmer_threshold,
        graph_name: String::from_utf8(graph_name).unwrap(),
    }))
}

pub fn mccortex(input: &[u8]) -> IResult<&[u8], McCortex> {
    let (remain, header) = header(input)?;
    let (remain, sample_names) = count(sample_name, header.cols as usize)(remain)?;

    // Skip error rate because Rust doesn't have long double (either 80 or 128 bits) yet
    let to_skip = 16 * header.cols as usize;
    let remain = &remain[to_skip..];

    // Skip cleaning information temporary
    // let (remain, cleaning) = count(cleaning, header.cols as usize)(remain)?;
    let (remain, _) = take_until(SIGNATURE)(remain)?;
    let cleaning = vec![Cleaning{
        top_clip: false,
        remove_low_covg_supernodes: false,
        remove_low_covg_kmers: false,
        cleaned_against_graph: false,
        remove_low_coverage_supernodes_threshold: 0,
        remove_low_coverage_kmer_threshold: 0,
        graph_name: "".to_string()
    }];

    let (remain, _) = tag(SIGNATURE)(remain)?;

    Ok((remain, McCortex {
        header,
        sample_names,
        cleaning,
    }))
}

const TEST_FILE: &'static [u8] = include_bytes!("../../tests/data/sample.ctx");

#[test]
fn parse_header() {
    let parsed = match header(&TEST_FILE[..]) {
        Ok((_, parsed)) => parsed,
        Err(e) => panic!(e)
    };
    assert_eq!(
        parsed,
        Header { version: 6, kmer_size: 31, W: 1, cols: 1, mean_read_lens: vec![177], total_seq_loaded: vec![3918788033] }
    );
}

#[test]
fn parse_sample_name() {
    let parsed = match sample_name(&TEST_FILE[34..]) {
        Ok((_, parsed)) => parsed,
        Err(e) => panic!(e)
    };
    assert_eq!(
        parsed,
        SampleName { len: 6, value: "sample".parse().unwrap() }
    );
}