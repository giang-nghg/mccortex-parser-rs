use std::fs::{File};
use crate::parsers::mccortex::{mccortex};
use std::env;

pub mod parsers;

fn main() {
  let args: Vec<String> = env::args().collect();

  let file = File::open(&args[1]).unwrap();
  match mccortex(|s| print!("{:?}\n", s), |s| print!("{:?}\n", s))(&file) {
    Err(e) => panic!(e),
    _ => {}
  };
}