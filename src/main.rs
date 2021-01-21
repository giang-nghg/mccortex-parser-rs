use std::fs::read;
use std::error;

pub mod parsers;

fn main() -> Result<(), Box<dyn error::Error>>{
  let file = read("tests/data/sample.ctx")?;
  let result = parsers::mccortex::header(&file[..]);
  match result {
    Ok((remain, parsed)) => print!("DEBUG: {:?}", parsed),
    Err(_) => {}
  }

  Ok(())
}