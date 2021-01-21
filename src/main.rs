use std::fs::{File};
use crate::parsers::mccortex::{mccortex_stream};

pub mod parsers;

fn main() {
  let file = File::open("tests/data/sample.ctx").unwrap();
  match mccortex_stream(|s| print!("{:?}\n", s), |s| print!("{:?}\n", s))(&file) {
    Err(e) => panic!(e),
    _ => {}
  };
}