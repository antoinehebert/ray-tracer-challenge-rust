# Introduction
Learning Rust by implementing a basic ray tracer, following the book [The Ray Tracer Challenge](http://www.raytracerchallenge.com/) by [Jamis Buck](https://github.com/jamis).

# Prerequisite

-   [Rust](https://www.rust-lang.org/)

# Running tests

```bash
cargo test
```

# Running app

The app expects a `filename` parameter.

```bash
cargo run --release your_file_name.ppm
```

# TODOs
-   Optimize! This is slow.
    -   pre-compute transforms
    -   pre-compute bounding boxes
    -   parallelize using [rayon](https://github.com/rayon-rs/rayon)?
