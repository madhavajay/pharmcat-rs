# PharmCAT Rust Port TODO

Goal: port the Java PharmCAT library in `repos/PharmCAT/` to Rust, preserving behavior until the Rust implementation is byte-for-byte identical for covered inputs and outputs.

## Guiding Principles

- [ ] Treat tests as the spec. Every Rust behavior area should trace back to a Java test, Java fixture, or explicit parity fixture.
- [ ] Keep the Java source repos read-only references under `repos/`.
- [ ] Keep the production implementation pure Rust. Do not bind to Java or call the Java jar at runtime.
- [ ] Use Rust bioinformatics libraries where they cover file-format semantics, especially `noodles-vcf` and `noodles-bgzf`.
- [ ] Reimplement PharmCAT and the needed `pgkb-common` behavior inside this crate first; split helper crates only after real reuse appears.
- [ ] Preserve Java-observable behavior first, including quirks or bugs, then make explicit later decisions about whether to fix them.
- [ ] Document every intentional non-parity decision with the fixture or test that proves the difference is acceptable.
- [ ] Keep parity work in small slices. Do not mix broad refactors, dependency upgrades, and behavior changes in the same slice.

## 0. Repository Baseline

- [x] Add Java PharmCAT source as a git submodule at `repos/PharmCAT/`.
- [x] Add `noodles` source as a git submodule at `repos/noodles/`, branch `madhava/bioscript`.
- [x] Add Java `pgkb-common` source as a git submodule at `repos/pgkb-common/`.
- [x] Record the exact baseline commits in this file before the first Rust parity test:
  - `repos/PharmCAT`: `55e3cb30a078537b4bec63b8d2b5035a20bc2fc0` on `development`.
  - `repos/noodles`: `5868f00a1f4fa6a0e0c32b685819fd3ef67b6473` on `madhava/bioscript`.
  - `repos/pgkb-common`: `70756e19aef8b46c63f4757f874c9c5ed31e3908` on `main`.
- [ ] Decide whether CI should use submodules recursively or explicit `git clone` commands. Current scaffold clones into `repos/` like the Kestrel port.
- [ ] Add a contributor note explaining that `repos/PharmCAT` and `repos/pgkb-common` are behavioral references, not runtime dependencies.

## 1. CI And Test Gates

- [x] Add `.github/workflows/ci.yml` with independent jobs for source repo checkout, Java tests, and Rust tests.
- [x] Make CI clone reference repos into `repos/`:
  - `https://github.com/madhavajay/PharmCAT.git`, branch `development`.
  - `https://github.com/madhavajay/noodles.git`, branch `madhava/bioscript`.
  - `https://github.com/madhavajay/pgkb-common.git`, default branch.
- [x] Run Java and Rust jobs in parallel, following the Kestrel `rust.yml` / `java.yml` pattern.
- [x] Add Rust CI commands that activate once `Cargo.toml` exists:

  ```sh
  cargo fmt --all --check
  cargo clippy --workspace --all-targets --all-features -- -D warnings
  cargo test --workspace --all-features
  ```

- [x] Add Java CI command:

  ```sh
  cd repos/PharmCAT
  ./gradlew test
  ```

- [ ] Remove the temporary "skip Rust if no Cargo.toml" behavior after the Rust workspace skeleton lands.
- [ ] Add a parity CI job after the first Rust CLI/API exists:
  - build the Rust binary or test helper;
  - build or run the Java reference;
  - run both on identical fixtures;
  - diff outputs with documented normalization only.
- [x] Add Java test report upload and keep it green before broad Rust changes.
- [ ] Add coverage later, using Kestrel's pattern as the target:

  ```sh
  cargo install cargo-llvm-cov --locked
  cargo llvm-cov --workspace --all-features --lcov --output-path lcov.info
  ```

## 2. Java Source And Test Inventory

- [ ] Inspect `repos/PharmCAT/build.gradle` and record Java, Gradle, and dependency versions.
- [ ] Run the main Java/JUnit suite locally:

  ```sh
  cd repos/PharmCAT
  ./gradlew test
  ```

- [ ] Inventory Java tests under `repos/PharmCAT/src/test`.
- [ ] Inventory Java fixtures under `repos/PharmCAT/src/test/resources`.
- [ ] Record Gradle report locations:
  - `repos/PharmCAT/build/reports/tests/`
  - `repos/PharmCAT/build/reports/jacoco/`
- [ ] Document tests that require network access, large data, environment variables, or generated resources.
- [ ] Create a test map table: Java test class, Rust test module, fixtures used, parity status.
- [ ] Add a status file later if needed, modeled after `bcftools-rs/docs/test-status.md`.

## 3. Preprocessor Test Inventory

- [ ] Inspect `repos/PharmCAT/preprocessor/tests/README.md`.
- [ ] Install or document Python test dependencies from `repos/PharmCAT/preprocessor/tests/requirements.txt`.
- [ ] Run the preprocessor tests from `repos/PharmCAT/preprocessor/tests`.
- [ ] Inventory all preprocessor fixtures and expected outputs.
- [ ] Decide whether Rust will reimplement preprocessor behavior directly or keep it as a separate compatibility stage.
- [ ] If preprocessor behavior is in scope, add a separate parity matrix for VCF normalization, filters, expected output VCFs, and error cases.

## 4. `pgkb-common` Inventory

- [ ] Inventory every `org.pharmgkb.common` import used by PharmCAT.
- [ ] Classify each imported item as:
  - [ ] standard Rust replacement;
  - [ ] small local adapter;
  - [ ] parity-sensitive behavior;
  - [ ] CLI/test-only behavior;
  - [ ] intentionally out of scope.
- [ ] Port `ChromosomeNameComparator` and `ChromosomePositionComparator` early; sort order can affect byte-for-byte output.
- [ ] Add ordering tests for `chr1`, `chr2`, `chr10`, `chrX`, `chrY`, `chrM`, bare contigs, and nonstandard contigs.
- [ ] Reimplement utility behavior locally as needed:
  - `CliHelper`
  - `AnsiConsole`
  - `ComparisonChain`
  - `ComparatorUtils`
  - `PathUtils`
  - `StreamUtils`
  - `IoUtils`
  - `NoDuplicateMergeFunction`
  - `Throttler`
  - `TimeUtils`

## 5. Rust Workspace Skeleton

- [ ] Initialize a Cargo workspace at the repository root.
- [ ] Start with one crate unless the implementation naturally splits:

  ```text
  pharmcat-rs/
  ├── crates/
  │   └── pharmcat/
  │       ├── src/lib.rs
  │       └── src/bin/pharmcat.rs
  ├── repos/
  │   ├── PharmCAT/
  │   ├── noodles/
  │   └── pgkb-common/
  ├── tests/
  └── fixtures/
  ```

- [ ] Add `rust-toolchain.toml` after choosing the minimum Rust version.
- [ ] Add workspace dependencies conservatively:
  - `thiserror` or `anyhow` for errors;
  - `clap` for CLI compatibility;
  - `serde` / `serde_json` for JSON output;
  - `tracing` for logging;
  - `tempfile` for parity tests;
  - `pretty_assertions` for readable test diffs;
  - `rstest` for Java-style parameterized tests.
- [ ] Add `noodles` dependency through the local submodule:

  ```toml
  noodles = { git = "https://github.com/madhavajay/noodles.git", branch = "madhava/bioscript", features = ["bgzf", "vcf"] }
  ```

- [ ] Keep `.cargo/config.toml` CI patching to `repos/noodles/noodles` so CI validates the checked-out source.

## 6. VCF Adapter

- [ ] Replace Java `org.pharmgkb:vcf-parser` with a small PharmCAT-specific adapter over `noodles-vcf` and `noodles-bgzf`.
- [ ] Match Java `VcfSampleReader` behavior:
  - sample names from header;
  - contig assembly extraction;
  - error on mixed contig assemblies;
  - no need to parse data rows for sample discovery.
- [ ] Match Java `VcfReader` behavior:
  - metadata parsing;
  - `FORMAT/AD` header handling and warnings;
  - selected sample lookup;
  - `GT`, `AD`, and `PS` sample field parsing;
  - REF/ALT/FILTER handling;
  - duplicate position handling;
  - gzip/BGZF input;
  - phased and effectively phased decisions;
  - haploid chromosome handling;
  - warnings and discarded-position behavior.
- [ ] Add fixture-backed Rust tests using Java VCF fixtures first.
- [ ] Escalate to `htslib-rs` only if indexed queries, BCF support, or HTSlib-compatible edge behavior becomes necessary.

## 7. Port Order

- [ ] Port data models with no external behavior first.
- [ ] Port fixture/data loading next.
- [ ] Port VCF parsing and sample handling.
- [ ] Port allele definition data loading.
- [ ] Port named allele matching.
- [ ] Port phasing and diplotype resolution.
- [ ] Port phenotype lookup and activity score handling.
- [ ] Port recommendation/report data loading.
- [ ] Port JSON output generation.
- [ ] Port TSV/report output generation.
- [ ] Port HTML/report output generation if in scope.
- [ ] Port CLI argument compatibility if in scope.
- [ ] Port error messages and exit codes for documented user-facing failures.

## 8. Test Porting Strategy

- [ ] For each Java test class under `repos/PharmCAT/src/test`, create a matching Rust test module or integration test.
- [ ] Preserve test intent and fixture coverage even when the Rust implementation structure differs.
- [ ] Use `rstest` for Java parameterized tests.
- [ ] Prefer golden-file tests for externally visible behavior.
- [ ] Reuse Java fixtures directly where possible.
- [ ] Copy only small immutable fixtures into `fixtures/` when direct reuse would make tests brittle.
- [ ] Add focused Rust unit tests only for internal edge cases discovered while porting.
- [ ] Track known Java quirks or bugs in a "bug decisions" table before changing behavior.

## 9. Java-vs-Rust Parity Harness

- [ ] Define a repeatable command that runs Java PharmCAT for selected fixtures and writes outputs to a temporary baseline directory.
- [ ] Define the matching Rust command/API path and output directory.
- [ ] Add byte-for-byte comparison of Java and Rust outputs.
- [ ] Normalize only unavoidable nondeterminism, and document every rule:
  - timestamps;
  - absolute paths;
  - version strings;
  - temporary directory names;
  - ordering only if Java order is nondeterministic and proven so.
- [ ] Start with one small VCF fixture, then expand by feature area.
- [ ] Keep failing parity cases checked in as ignored or named tests until implemented.
- [ ] Add an env var gate like Kestrel's `KESTREL_RUN_JAVA_PARITY`; for this repo use:

  ```text
  PHARMCAT_RUN_JAVA_PARITY=1
  PHARMCAT_JAVA_DIR=repos/PharmCAT
  ```

- [ ] CI parity job should run once it is stable enough to avoid false greens.

## 10. Upstream-Style Harness Discipline

- [ ] Add scripts for stable local gates:

  ```sh
  ./test-java.sh
  ./test-rust.sh
  ./test-parity.sh
  ```

- [ ] Make scripts fail loudly if required tools are missing. The `samtools-rs` TODO noted a false-green caused by missing `bgzip`; avoid silent skips.
- [ ] Print exact pass/fail/ignored counts for parity groups.
- [ ] Keep a promoted parity subset separate from exploratory tests until the full suite is green.
- [ ] Add artifacts for failed output diffs so CI failures are easy to inspect.

## 11. Byte-For-Byte Parity Gates

- [ ] Full Java test suite passes from a clean checkout.
- [ ] Full Rust test suite passes from a clean checkout.
- [ ] Java-vs-Rust parity tests pass across selected fixtures.
- [ ] Generated JSON/TSV/report outputs are byte-for-byte identical.
- [ ] Sort order, whitespace, number formatting, null handling, and missing-data behavior match Java.
- [ ] Warnings and recoverable errors match where user-visible compatibility matters.
- [ ] Exit codes match where CLI compatibility is in scope.
- [ ] Any non-identical behavior is documented with rationale and tests.

## 12. Dependency Blocker Log

- [ ] Keep dependency-blocked work in this TODO instead of patching dependencies opportunistically.
- [ ] If `noodles` lacks needed VCF behavior, record:
  - fixture;
  - expected Java behavior;
  - current Rust/library behavior;
  - proposed upstream fix;
  - temporary local adapter, if any.
- [ ] If `pgkb-common` behavior is unclear, add a focused Java probe or fixture before porting.
- [ ] Do not bypass `noodles` with hand parsing unless the adapter is explicitly documented as PharmCAT-specific behavior.

## 13. Documentation And Release Readiness

- [ ] Document supported Rust API surface.
- [ ] Document CLI compatibility status.
- [ ] Document data update workflow.
- [ ] Add benchmark fixtures comparing Java and Rust runtime.
- [ ] Add regression fixtures for every parity bug found during porting.
- [ ] Decide when the Java submodules become test-only reference material rather than active implementation references.
