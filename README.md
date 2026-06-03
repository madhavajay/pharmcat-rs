# PharmCAT Rust Port

This repository is a pure Rust port of PharmCAT. The Java sources under
`repos/PharmCAT` and `repos/pgkb-common` are behavioral references for tests and
parity work, not runtime dependencies for the Rust implementation.

The Rust port should preserve Java-observable behavior first. Any intentional
behavior change needs a fixture-backed test and an explicit note in `TODO.md`.

## Local Gates

```sh
cargo fmt --all --check
cargo clippy --workspace --all-targets --all-features -- -D warnings
cargo test -p pharmcat --all-features
```

The Java reference suite is run from the submodule:

```sh
cd repos/PharmCAT
./gradlew test
```

