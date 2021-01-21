use std::fs::{File};
use crate::parsers::mccortex::{mccortex_stream};

pub mod parsers;

fn main() {
  let file = File::open("tests/data/sample.ctx").unwrap();
  match mccortex_stream(&file, |s| print!("{:?}\n", s), |s| ()) {
    Err(e) => panic!(e),
    _ => {}
  };
}