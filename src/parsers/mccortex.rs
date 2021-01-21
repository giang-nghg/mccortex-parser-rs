use nom::{IResult};
use nom::number::complete::{le_u8, le_u32, le_u64};
use nom::multi::count;
use nom::bytes::complete::tag;
use nom::sequence::tuple;

// Spec: https://github.com/mcveanlab/mccortex/blob/master/docs/file_formats/graph_file_format.txt
const SIGNATURE: &str = "CORTEX";

#[derive(Debug, PartialEq)]
pub struct McCortex {
    pub header: Header,
    pub sample_names: Vec<SampleName>,
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

pub fn mccortex(input: &[u8]) -> IResult<&[u8], McCortex> {
    let (remain, header) = header(input)?;
    let (remain, sample_names) = count(sample_name, header.cols as usize)(remain)?;

    Ok((remain, McCortex {
        header,
        sample_names,
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