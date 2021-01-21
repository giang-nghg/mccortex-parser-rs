use std::fs::read;
use std::error;
use std::str::from_utf8;

pub mod parsers;

fn main() -> Result<(), Box<dyn error::Error>>{
  let file = read("tests/data/sample.ctx")?;
  let result = parsers::mccortex::mccortex(&file[..]);
  match result {
    Ok((_remain, parsed)) => {
      print!("DEBUG: {:?}", parsed);
      print!("DEBUG: {:?}", from_utf8(&*parsed.sample_names[0].value));
    },
    Err(_) => {}
  }

  Ok(())
}