use nom::{IResult};
use nom::number::streaming::{le_u8, le_u32, le_u64};
use nom::multi::{count, length_count, many0};
use nom::bytes::streaming::{tag, take};
use nom::sequence::tuple;
use nom::combinator::map;
use std::fs::File;
use std::io::Read;

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
    pub n_words_per_kmer: u32,
    pub cols: u32,
    pub mean_read_lens: Vec<u32>,
    pub total_seq_loaded: Vec<u64>,
    pub sample_names: Vec<String>,
    pub cleaning: Vec<Cleaning>
}

#[derive(Debug, PartialEq)]
pub struct Kmer {
    pub raw: String,
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

pub fn header(remain: &[u8]) -> IResult<&[u8], Header> {
    let (remain, _) = tag(SIGNATURE)(remain)?;

    let (remain, (version, kmer_size, n_words_per_kmer, cols)) = tuple((le_u32, le_u32, le_u32, le_u32))(remain)?;
    let (remain, mean_read_lens) = count(le_u32, cols as usize)(remain)?;
    let (remain, total_seq_loaded) = count(le_u64, cols as usize)(remain)?;

    let (remain, sample_names) = count(length_string, cols as usize)(remain)?;

    // Skip error rate because Rust doesn't have long double (128 bits) yet
    let (remain, _) = take(16 * cols)(remain)?;

    let (remain, cleaning) = count(cleaning, cols as usize)(remain)?;

    let (remain, _) = tag(SIGNATURE)(remain)?;

    Ok((remain, Header {
        version,
        kmer_size,
        n_words_per_kmer,
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
    let (remain, graph_name) = length_string(remain)?;

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

fn nu(b: u64) -> char {
    if b < 2u64.pow(62) { 'A' }
    else if b >= 2u64.pow(62) && b < 2u64.pow(63) { 'C' }
    else if b >= 2u64.pow(63) && b < 2u64.pow(63) + 2u64.pow(62) { 'G' }
    else { 'T' }
}

pub fn kmer(n_words_per_kmer: u32, cols: u32) -> Box<dyn Fn(&[u8]) -> IResult<&[u8], Kmer>> {
    Box::new(move |input| {
        let (remain, raw) = count(le_u64, n_words_per_kmer as usize)(input)?;
        let (remain, colour_covs) = count(le_u32, cols as usize)(remain)?;
        let (remain, colour_edges) = count(le_u8, cols as usize)(remain)?;

        let mut bases = Vec::new();
        for word in &raw {
            let mut word_copy = *word;
            for _ in 0..32 {
                word_copy = word_copy << 2;
                bases.push(nu(word_copy));
            }
        }

        Ok((remain, Kmer{
            raw: bases.into_iter().collect::<String>(), colour_covs, colour_edges
        }))
    })
}

pub fn mccortex(input: &[u8]) -> IResult<&[u8], McCortex> {
    let (remain, header) = header(input)?;
    let (remain, kmers) = many0(kmer(header.n_words_per_kmer, header.cols))(remain)?;

    Ok((remain, McCortex {
        header,
        kmers
    }))
}

pub fn mccortex_stream(on_header: fn(Header), on_kmer: fn(Kmer)) -> Box<dyn Fn(&File) -> IResult<Vec<u8>, ()>> {
    Box::new(move |input| {
        let mut buf = Vec::new();
        let (mut n_words_per_kmer, mut cols) = (0, 0);
        let mut is_header_parsed = false;

        for byte in input.bytes() {
            buf.push(byte.unwrap());

            if !is_header_parsed {
                match header(&*buf) {
                    Err(nom::Err::Incomplete(..)) => continue,
                    Ok((remain, header)) => {
                        is_header_parsed = true;
                        n_words_per_kmer = header.n_words_per_kmer;
                        cols = header.cols;
                        buf = Vec::from(remain);

                        on_header(header);
                    },
                    Err(nom::Err::Error(e)) => {
                        return Err(nom::Err::Error(nom::error::Error { input: Vec::from(e.input), code: e.code }));
                    }
                    Err(nom::Err::Failure(e)) => {
                        return Err(nom::Err::Failure(nom::error::Error { input: Vec::from(e.input), code: e.code }));
                    }
                }
            }

            match kmer(n_words_per_kmer, cols)(&*buf) {
                Err(nom::Err::Incomplete(..)) => continue,
                Ok((remain, kmer)) => {
                    buf = Vec::from(remain);
                    on_kmer(kmer);
                },
                Err(nom::Err::Error(e)) => {
                    return Err(nom::Err::Error(nom::error::Error { input: Vec::from(e.input), code: e.code }));
                }
                Err(nom::Err::Failure(e)) => {
                    return Err(nom::Err::Failure(nom::error::Error { input: Vec::from(e.input), code: e.code }));
                }
            }
        }

        Ok((buf, ()))
    })
}

fn length_string(input: &[u8]) -> IResult<&[u8], String> {
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
            n_words_per_kmer: 1,
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
    let parsed = match kmer(header.n_words_per_kmer, header.cols)(&remain[..13]) {
        Ok((_, parsed)) => parsed,
        Err(e) => panic!(e)
    };
    assert_eq!(
        parsed,
        Kmer { raw: String::from("CTCGGCTACCCCGAACTCCAGCGAGAAGTACA"), colour_covs: vec![759], colour_edges: vec![68] }
    );
}
