use nom::{IResult};
use nom::number::complete::{le_u32, le_u64};
use nom::multi::count;
use nom::bytes::complete::tag;
use nom::sequence::tuple;
use std::error::Error;

// Spec: https://github.com/mcveanlab/mccortex/blob/master/docs/file_formats/graph_file_format.txt
const SIGNATURE: &str = "CORTEX";

#[derive(Debug, PartialEq)]
pub struct Header {
    pub version: u32,
    pub kmer_size: u32,
    pub W: u32,
    pub cols: u32,
    pub mean_read_lens: Vec<u32>,
    pub total_seq_loaded: Vec<u64>,
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

const TEST_FILE: &'static [u8] = include_bytes!("../../tests/data/sample.ctx");

#[test]
fn parse_header() -> Result<(), Box<dyn Error>> {
    let (_, parsed) = header(&TEST_FILE[..])?;
    assert_eq!(
        parsed,
        Header { version: 6, kmer_size: 31, W: 1, cols: 1, mean_read_lens: vec![177], total_seq_loaded: vec![3918788033] }
    );
    Ok(())
}