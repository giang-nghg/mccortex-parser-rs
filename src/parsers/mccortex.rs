use nom::{IResult};
use nom::number::complete::{le_u8, le_u32, le_u64};
use nom::multi::{count, length_count, many0};
use nom::bytes::complete::{tag};
use nom::sequence::tuple;
use nom::combinator::map;

// Spec: https://github.com/mcveanlab/mccortex/blob/master/docs/file_formats/graph_file_format.txt
const SIGNATURE: &str = "CORTEX";

#[derive(Debug, PartialEq)]
pub struct McCortex {
    pub header: Header,
    pub kmers: Vec<Kmer>,
}

#[derive(Debug, PartialEq)]
pub struct Header {
    pub version: u32,
    pub kmer_size: u32,
    pub W: u32,
    pub cols: u32,
    pub mean_read_lens: Vec<u32>,
    pub total_seq_loaded: Vec<u64>,
    pub sample_names: Vec<String>,
    pub cleaning: Vec<Cleaning>
}

#[derive(Debug, PartialEq)]
pub struct Kmer {
    pub raw: Vec<u64>,
    pub colour_covs: Vec<u32>,
    pub colour_edges: Vec<u8>,
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

    let (remain, sample_names) = count(string_parser, cols as usize)(remain)?;

    // Skip error rate because Rust doesn't have long double (128 bits) yet
    let to_skip = 16 * cols as usize;
    let remain = &remain[to_skip..];

    let (remain, cleaning) = count(cleaning, cols as usize)(remain)?;

    let (remain, _) = tag(SIGNATURE)(remain)?;

    Ok((remain, Header {
        version,
        kmer_size,
        W,
        cols,
        mean_read_lens,
        total_seq_loaded,
        sample_names,
        cleaning,
    }))
}

pub fn cleaning(input: &[u8]) -> IResult<&[u8], Cleaning> {
    let (remain, top_clip) = le_u8(input)?;
    let (remain, remove_low_covg_supernodes) = le_u8(remain)?;
    let (remain, remove_low_covg_kmers) = le_u8(remain)?;
    let (remain, cleaned_against_graph) = le_u8(remain)?;
    let (remain, remove_low_coverage_supernodes_threshold) = le_u32(remain)?;
    let (remain, remove_low_coverage_kmer_threshold) = le_u32(remain)?;
    let (remain, graph_name) = string_parser(remain)?;

    Ok((remain, Cleaning {
        top_clip: top_clip != 0,
        remove_low_covg_supernodes: remove_low_covg_supernodes != 0,
        remove_low_covg_kmers: remove_low_covg_kmers != 0,
        cleaned_against_graph: cleaned_against_graph != 0,
        remove_low_coverage_supernodes_threshold,
        remove_low_coverage_kmer_threshold,
        graph_name,
    }))
}

pub fn kmer(input: &[u8], W: u32, cols: u32) -> IResult<&[u8], Kmer> {
    let (remain, raw) = count(le_u64, W as usize)(input)?;
    let (remain, colour_covs) = count(le_u32, cols as usize)(remain)?;
    let (remain, colour_edges) = count(le_u8, cols as usize)(remain)?;

    Ok((remain, Kmer{
        raw, colour_covs, colour_edges
    }))
}

pub fn mccortex(input: &[u8]) -> IResult<&[u8], McCortex> {
    let (remain, header) = header(input)?;
    let (remain, kmers) = many0(|b| kmer(b, header.W, header.cols))(remain)?;

    Ok((remain, McCortex {
        header,
        kmers
    }))
}

fn string_parser(input: &[u8]) -> IResult<&[u8], String> {
    let (remain, parsed) = map(length_count(le_u32, le_u8), |s| String::from_utf8(s).unwrap())(input)?;
    Ok((remain, parsed))
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
        Header {
            version: 6,
            kmer_size: 31,
            W: 1,
            cols: 1,
            mean_read_lens: vec![177],
            total_seq_loaded: vec![3918788033],
            sample_names: vec!["sample".to_string()],
            cleaning: vec![Cleaning { top_clip: true, remove_low_covg_supernodes: true, remove_low_covg_kmers: false, cleaned_against_graph: false, remove_low_coverage_supernodes_threshold: 26, remove_low_coverage_kmer_threshold: 4294967295, graph_name: "undefined".to_string() }]
        }
    );
}

#[test]
fn parse_kmer() {
    let (remain, header) = match header(&TEST_FILE[..]) {
        Ok((remain, header)) => (remain, header),
        Err(e) => panic!(e)
    };
    let parsed = match kmer(&remain[..13], header.W, header.cols) {
        Ok((_, parsed)) => parsed,
        Err(e) => panic!(e)
    };
    assert_eq!(
        parsed,
        Kmer { raw: vec![2136700052237066417], colour_covs: vec![759], colour_edges: vec![68] }
    );
}
