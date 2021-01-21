# mccortex-parser-rs
An eager McCortex parser. It streams a kmer out as soon as it's parsed. Currently it relies on callbacks (Rust closures) to do so.
The drawback of this approach is subsequent code will be forced to built upon the callback structure as well. However, once generators make it to Rust stable,
this should be rewritten to allow more freedom for users.

For example:
```rust
let file = File::open("sample.ctx").unwrap();

// 1st closure receives the header, 2nd closure receives the kmers
match mccortex(|s| print!("{:?}\n", s), |s| print!("{:?}\n", s))(&file) {
  Err(e) => panic!(e),
  _ => {}
};
```
