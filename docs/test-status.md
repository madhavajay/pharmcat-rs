# PharmCAT Java Test Status

Initial inventory for porting Java PharmCAT tests into the Rust implementation.

This file is intentionally conservative. `Partial Rust coverage` means there is a Rust owner module with tests for some of the same behavior; it does not mean every Java test method in that class has been ported or that byte-for-byte Java parity is proven. Full parity still requires a runnable Java reference gate plus promoted Java-vs-Rust fixture comparisons.

## Baseline

- Java reference source: `repos/PharmCAT` at `55e3cb30a078537b4bec63b8d2b5035a20bc2fc0`.
- Inventory command:

```sh
find repos/PharmCAT/src/test/java -name '*Test.java' | sort
rg -c '^\s*@(Test|ParameterizedTest)' repos/PharmCAT/src/test/java -g '*Test.java'
```

- Current Java inventory: 63 `*Test.java` files.
- Current JUnit inventory: 62 files with JUnit test annotations, 426 `@Test` / `@ParameterizedTest` annotations.
- `SyntheticBatchTest.java` is a generator utility with `main`, not a JUnit test class.
- Java gate status: blocked locally until Gradle can resolve the wrapper distribution/cache.
- Current Rust broad gate recorded in `TODO.md`: `timeout 750s cargo test -p pharmcat --all-features` passed with 211 tests.

## Status Legend

- `Direct Rust coverage`: a Rust test is explicitly tied to the same Java method intent or fixture.
- `Partial Rust coverage`: Rust tests cover some Java behavior in the listed owner module.
- `Not yet mapped`: Java class is inventoried, but method-level Rust ownership has not been audited.
- `Reference utility`: Java helper/generator, not a direct JUnit test target.
- `Blocked`: needs Java reference execution or promoted Java-vs-Rust fixtures before completion can be claimed.

## Top-Level And Pipeline Tests

| Java test class | Java tests | Rust owner | Status | Notes / next action |
| --- | ---: | --- | --- | --- |
| `BaseConfigTest` | 1 | `crates/pharmcat/src/cli.rs`, `pipeline.rs` | Partial Rust coverage | CLI/base config parsing has Rust tests; map the single Java assertion explicitly. |
| `BatchPharmCATTest` | 13 | `pipeline.rs` | Partial Rust coverage | Rust has single-sample and multi-gene pipeline slices; batch/multisample parity still needs method-level inventory. |
| `Cacna1sTest` | 5 | `report.rs`, `pipeline.rs` | Not yet mapped | Gene-specific end-to-end reporter behavior needs fixture ownership. |
| `CftrTest` | 4 | `matcher.rs`, `report.rs` | Not yet mapped | CFTR matcher fixtures are not yet method-mapped. |
| `Cyp2d6Test` | 9 | `report.rs`, `phenotype.rs` | Partial Rust coverage | CYP2D6 outside-call/suballele and static-message behavior has some Rust coverage; copy-number caller remains separate. |
| `DpydTest` | 34 | `matcher.rs`, `phenotype.rs`, `report.rs`, `pipeline.rs` | Partial Rust coverage | DPYD HapB3, lowest-function, phenotype, reporter, and pipeline slices have Rust tests; full Java method map still needed. |
| `PharmCATTest` | 15 | `cli.rs`, `pipeline.rs`, `report.rs` | Partial Rust coverage | Full CLI/output parity is only partially represented by current CYP3A5 path. |
| `PipelineTest` | 81 | `pipeline.rs`, `cli.rs`, `vcf.rs`, `report.rs` | Partial Rust coverage | Largest Java surface; current Rust tests cover output planning, VCF-to-reporter slices, matcher JSON/HTML intermediates, and warnings, but not all Java cases. |
| `Ryr1Test` | 13 | `matcher.rs`, `phenotype.rs`, `report.rs`, `pipeline.rs` | Partial Rust coverage | RYR1 lowest-function path has real-definition Rust coverage; full fixture sweep remains. |
| `SyntheticBatchTest` | n/a | n/a | Reference utility | Generator utility. Track separately if a Rust synthetic fixture generator is added. |
| `ToxicGenesTest` | 17 | `report.rs`, `phenotype.rs`, `pipeline.rs` | Not yet mapped | Needs explicit toxic-gene fixture inventory and Rust ownership. |

## PipelineTest Method Map

| Java method | Line | Fixture / focus | Closest Rust coverage | Status | Next action |
| --- | ---: | --- | --- | --- | --- |
| `testCounts` | 419 | Full report gene/drug inventory from reference CYP2C9 VCF | broad report/guidance loading tests | Partial Rust coverage | Add a fixture-backed report context count test once full resource loading parity is stable. |
| `testAll` | 481 | All matcher genes plus CYP2D6/HLA outside calls | `run_reporter_from_vcf_handles_standard_and_lowest_function_genes_in_one_run`; outside-call unit coverage | Partial Rust coverage | Add multi-gene pipeline fixture with outside calls after outside-call pipeline wiring is complete. |
| `testNoData` | 539 | Pipeline report with no input data | `pipeline::tests::run_reporter_from_empty_vcf_writes_no_data_report_like_java_pipeline_test` | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testUndocumentedVariation` | 546 | CYP2C19 undocumented variation no-call and compact HTML | `pipeline::tests::run_reporter_from_vcf_marks_cyp2c19_undocumented_variation_uncallable_like_java_pipeline_test`; definition-aware VCF warnings and report HTML markers covered end-to-end | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testUndocumentedVariationExtendedReport` | 564 | CYP2C19 undocumented variation in extended report | `pipeline::tests::run_reporter_from_vcf_marks_cyp2c19_undocumented_variation_uncallable_in_extended_report_like_java_pipeline_test`; compact-vs-extended Section II behavior covered end-to-end | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testUndocumentedVariationsWithTreatAsReference` | 583 | TPMT/RYR1 undocumented-as-reference plus CYP2C19 no-call | `pipeline::tests::run_reporter_from_vcf_treats_toxic_gene_undocumented_variations_as_reference_like_java_pipeline_test`; Java undocumented-as-reference warnings and report markers covered end-to-end | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testUndocumentedVariationsWithTreatAsReferenceAndCombo` | 620 | Combo mode with undocumented-as-reference TPMT and RYR1 | `pipeline::tests::run_reporter_from_vcf_combo_keeps_standard_toxic_gene_undocumented_variation_like_java_pipeline_test`; combo-mode TPMT custom SNP and lowest-function RYR1 treat-as-reference covered end-to-end | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testUncallable` | 657 | CYP2C19/TPMT uncallable genes and HTML no-call state | `pipeline::tests::run_reporter_from_vcf_marks_cyp2c19_and_tpmt_uncallable_like_java_pipeline_test`; exact CYP2C19/TPMT no-call and report output covered end-to-end | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testCyp2c19` | 696 | CYP2C19 `*1/*1` plus CYP2D6 outside call recommendation matching | `pipeline::tests::run_reporter_from_vcf_cyp2c19_with_cyp2d6_and_g6pd_outside_calls_like_java_pipeline_test`; full VCF pipeline now parses outside-call TSVs and asserts CYP2C19, CYP2D6/G6PD outside calls, and amitriptyline/citalopram matched annotations | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testCyp2c19_s1s2rs58973490het` | 717 | CYP2C19 ambiguity message applies with het variant | `pipeline::tests::run_reporter_from_vcf_cyp2c19_s1s2_het_ambiguity_message_like_java_pipeline_test`; full VCF pipeline now loads `messages.json`, applies gene and drug message matching, and asserts CYP2C19 `*1/*2`, heterozygous `rs58973490`, outside CYP2D6, drug annotation counts, and CPIC amitriptyline message count | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testCyp2c19_s1s2` | 766 | CYP2C19 ambiguity message suppressed with hom variant | `pipeline::tests::run_reporter_from_vcf_cyp2c19_s1s2_hom_suppresses_ambiguity_message_like_java_pipeline_test`; paired real-message-catalog fixture asserts CYP2C19 `*1/*2`, homozygous/reference `rs58973490`, unchanged drug annotation matches, zero CYP2C19 ambiguity gene messages, and zero CPIC amitriptyline ambiguity drug messages | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testClomipramineCall` | 813 | CYP2C19 `*2/*2` broad drug recommendation counts | `pipeline::tests::run_reporter_from_vcf_cyp2c19_s2s2_clomipramine_counts_like_java_pipeline_test`; real CYP2C19 VCF + outside-call fixture asserts `*2/*2` and Java matched-annotation counts/source presence for amitriptyline, clomipramine, desipramine, doxepin, imipramine, nortriptyline, trimipramine, clopidogrel, lansoprazole, and voriconazole | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testCyp2c19noCall` | 866 | CYP2C19 no-call suppresses recommendation matches | `pipeline::tests::run_reporter_from_vcf_cyp2c19_no_call_suppresses_recommendation_matches_like_java_pipeline_test`; `matcher::tests::call_standard_gene_returns_no_call_for_cyp2c19_empty_match_set_like_java`; full CYP2C19 VCF + outside-call fixture asserts `NoCall`, `Unknown/Unknown`, outside CYP2D6/G6PD calls, and no CPIC/DPWG matches for citalopram or ivacaftor | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testCyp2c19s4bs17rs28399504missing` | 882 | CYP2C19 missing position yields `*4/*4`, `*4/*17`, `*17/*17` | `pipeline::tests::run_reporter_from_vcf_cyp2c19_s4b_s17_rs28399504_missing_like_java_pipeline_test`; `report::tests::recommendation_genotypes_expand_each_recommendation_diplotype_like_java`; full CYP2C19 VCF fixture omits `rs28399504`, asserts matcher/report calls, missing-position tracking, and citalopram CPIC/DPWG/FDA counts | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testCyp2c19s1s4het` | 902 | CYP2C19 `*4/*17` plus CYP2D6 outside call | `pipeline::tests::run_reporter_from_vcf_cyp2c19_s4_s17_with_cyp2d6_outside_call_like_java_pipeline_test`; exact CYP2C19 VCF + Java outside-call TSV fixture asserts CYP2C19 `*4/*17`, reportable CYP2D6 `*1/*4`, and outside-call marker | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testCyp2c19s1s4missingS1` | 922 | CYP2C19 partial missing `*1/*4`, `*4/*38` ambiguity | `pipeline::tests::run_reporter_from_vcf_cyp2c19_s1_s4_missing_s1_with_cyp2d6_outside_call_like_java_pipeline_test`; exact CYP2C19 VCF omits `rs3758581` like Java `TestVcfBuilder.missing`, uses the Java CYP2D6 outside-call TSV, and asserts matcher/report calls, missing/het variant reports, CYP2C19 ambiguity messages, and CPIC amitriptyline drug messages | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testCyp2c19SingleGeneMatch` | 960 | CYP2C19 single-gene `*1/*38` with outside CYP2D6 | `pipeline::tests::run_reporter_from_vcf_cyp2c19_s1_s38_with_cyp2d6_outside_call_like_java_pipeline_test`; exact CYP2C19 reference VCF fixture sets `rs3758581` `A/G`, omits `rs56337013` like Java `TestVcfBuilder.missing`, and asserts CYP2C19 `*1/*38`, CYP2D6 outside `*1/*4`, and missing-position report surface | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `multipleCalls` | 981 | CYP2C19 partial `*1/*4`, `*1/*17` multiple calls | `pipeline::tests::run_reporter_from_vcf_cyp2c19_multiple_calls_like_java_pipeline_test`; exact CYP2C19 VCF fixture sets `rs12248560` `C/T`, `rs3758581` `G/G`, omits `rs28399504` like Java `TestVcfBuilder.missing`, and asserts matcher/report calls plus missing-position report surface | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testRosuvastatin` | 1001 | ABCG2 + SLCO1B1 rosuvastatin and DPYD no-data HTML | `pipeline::tests::run_reporter_from_vcf_rosuvastatin_with_dpyd_no_data_like_java_pipeline_test`; exact ABCG2/SLCO1B1 VCF fixture sets `rs2231142` `G/T` and `rs56101265` `T/C`, loads DPYD with no rows, asserts ABCG2/SLCO1B1 matcher calls, SLCO1B1 `*1/*2`, rosuvastatin matched annotation count, no capecitabine section, and DPYD no-data HTML marker | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testSlco1b1HomWild` | 1024 | SLCO1B1 `*1/*1`, simvastatin recommendations | `pipeline::tests::run_reporter_from_vcf_slco1b1_hom_wild_simvastatin_like_java_pipeline_test`; exact SLCO1B1 reference VCF fixture asserts matcher/source/recommended `*1/*1`, simvastatin matched annotation count, DPWG `No recommendation` classification, and reporter HTML simvastatin/SLCO1B1 markers | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testSlco1b1HomVar` | 1049 | SLCO1B1 `*5/*15` | `pipeline::tests::run_reporter_from_vcf_slco1b1_hom_var_simvastatin_like_java_pipeline_test`; exact SLCO1B1 VCF fixture sets `rs2306283` `A/G` and `rs4149056` `C/C`, asserts matcher/source/recommended `*5/*15`, CPIC and DPWG simvastatin matched annotation counts, and reporter HTML simvastatin/SLCO1B1 markers | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testSlco1b1Test2` | 1071 | SLCO1B1 `*1/*44` | `pipeline::tests::run_reporter_from_vcf_slco1b1_s1_s44_simvastatin_like_java_pipeline_test`; exact SLCO1B1 VCF fixture sets `rs2306283`, `rs11045852`, and `rs74064213` to `A/G`, asserts matcher/source/recommended `*1/*44`, simvastatin matched annotation counts, and reporter HTML simvastatin/SLCO1B1 markers | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testSlco1b1Test3` | 1093 | SLCO1B1 `*1/*15` | `pipeline::tests::run_reporter_from_vcf_slco1b1_s1_s15_simvastatin_like_java_pipeline_test`; exact SLCO1B1 VCF fixture sets `rs2306283` to `A/G` and `rs4149056` to `T/C`, asserts matcher/source/recommended `*1/*15`, CPIC and DPWG simvastatin matched annotation counts, and reporter HTML simvastatin/SLCO1B1 markers | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testSlco1b1Test4` | 1115 | SLCO1B1 `*5/*45` plus intermediate files | `pipeline::tests::run_reporter_from_vcf_slco1b1_s5_s45_simvastatin_intermediates_like_java_pipeline_test`; exact SLCO1B1 VCF fixture sets `rs4149056` to `T/C` and `rs71581941` to `C/T`, asserts matcher/source/recommended `*5/*45`, simvastatin CPIC/DPWG/FDA source behavior, and matcher/phenotyper/reporter intermediate outputs | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testSlco1b1Test5` | 1141 | SLCO1B1 `*1/*45` warning | `pipeline::tests::run_reporter_from_vcf_slco1b1_s1_s45_warning_intermediates_like_java_pipeline_test`; exact SLCO1B1 VCF fixture sets `rs2306283` to `A/G`, `rs4149056` to `T/C`, and `rs71581941` to `C/T`, asserts matcher/source/recommended `*1/*45`, the `SLCO1B1 *1/*45 warning` gene message, simvastatin CPIC/DPWG/FDA source behavior, and matcher/phenotyper/reporter intermediate outputs | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testSlco1b1UncalledOverride` | 1176 | SLCO1B1 uncalled rs4149056 override | `pipeline::tests::run_reporter_from_vcf_slco1b1_uncalled_override_like_java_pipeline_test`; exact SLCO1B1 VCF fixture sets `rs2306283`, `rs11045853`, and `rs72559748` to homozygous Java variation calls plus heterozygous `rs4149056`, asserts matcher no-call, source `Unknown/Unknown`, inferred recommendation `*1/*5`, print-call source `rs4149056 C/rs4149056 T`, and simvastatin CPIC/DPWG/FDA-association matches | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testUgt1a1Phased` | 1203 | UGT1A1 phased `*1/*80` | `pipeline::tests::run_reporter_from_vcf_ugt1a1_s1_s80_phased_like_java_pipeline_test`; exact UGT1A1 VCF fixture phases all definition rows, sets `rs887829` to `C|T`, asserts matcher/report phased state, matcher/source/recommended `*1/*80`, recommendation alleles `*1` and `*80`, `Indeterminate` lookup, and reporter JSON/HTML diplotype markers | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testUgt1a1Unphased` | 1216 | UGT1A1 unphased `*1/*80` | `pipeline::tests::run_reporter_from_vcf_ugt1a1_s1_s80_unphased_like_java_pipeline_test`; exact UGT1A1 VCF fixture follows Java's true unphased matcher/phenotyper fixture by setting `rs887829` to `C/T` with unphased reference rows, asserts matcher/report unphased state, matcher/source/recommended `*1/*80`, recommendation alleles `*1` and `*80`, `Indeterminate` lookup, and reporter JSON/HTML diplotype markers | Direct Rust coverage | Java `PipelineTest.testUgt1a1Unphased` currently calls `.phased()` despite the test name; promote to Java-vs-Rust fixture later once the reference gate is available. |
| `testUgt1a1s1s1` | 1229 | UGT1A1 `*1/*1` | `pipeline::tests::run_reporter_from_vcf_ugt1a1_s1_s1_reference_like_java_pipeline_test`; exact UGT1A1 reference VCF fixture sets all definition rows to `0/0`, asserts matcher/source/recommended `*1/*1`, recommendation alleles `*1` and `*1`, `Normal Metabolizer` lookup, and reporter JSON/HTML diplotype markers | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testUgt1a1S1S80S28` | 1241 | UGT1A1 `*1/*80+*28` | combination name/report helpers | Partial Rust coverage | Add UGT1A1 combo diplotype fixture. |
| `testUgt1a1S28S37` | 1254 | UGT1A1 `*28/*37` | none specific | Not yet mapped | Add UGT1A1 repeat variant fixture. |
| `testUgt1a1s28s80phased` | 1266 | UGT1A1 phased combo, atazanavir messages | message/report tests | Partial Rust coverage | Add phased atazanavir message fixture. |
| `testUgt1a1s28s80s6s60phased` | 1289 | UGT1A1 phased `*6/*80+*28` | none specific | Not yet mapped | Add UGT1A1 phased multi-variant fixture. |
| `testUgt1a1s28s80s6s60unphased` | 1304 | UGT1A1 unphased `*6/*80+*28` | none specific | Not yet mapped | Add UGT1A1 unphased multi-variant fixture. |
| `testUgt1a1s6s6` | 1318 | UGT1A1 `*6/*6` | none specific | Not yet mapped | Add UGT1A1 homozygous fixture. |
| `testUgt1a1s6s60s80s28MissingPhased` | 1330 | UGT1A1 phased missing repeat position with 3 calls | missing-position matcher/report tests | Partial Rust coverage | Add UGT1A1 phased missing-repeat fixture. |
| `testUgt1a1s6s60s80s28MissingUnphased` | 1347 | UGT1A1 unphased missing repeat position with 3 calls | missing-position matcher/report tests | Partial Rust coverage | Add UGT1A1 unphased missing-repeat fixture. |
| `testUgt1a1s80s28missing` | 1364 | UGT1A1 missing `rs3064744` and atazanavir HTML grouping | report HTML recommendation tests | Partial Rust coverage | Add UGT1A1 atazanavir HTML fixture. |
| `testUgt1a1na12717` | 1391 | UGT1A1 `*80/*80+*28` | none specific | Not yet mapped | Add NA12717-inspired UGT1A1 fixture. |
| `testUgt1a1s28homMissing` | 1404 | UGT1A1 hom repeat with missing `rs887829` | missing-position tests | Partial Rust coverage | Add UGT1A1 hom-missing fixture. |
| `testUgt1a1s28s60Hom` | 1419 | UGT1A1 `*1/*28` | none specific | Not yet mapped | Add UGT1A1 `*28` fixture. |
| `testUgt1a1s27s28unphaseds80s60missing` | 1431 | UGT1A1 `*27/*28`, `*27/*80+*28` | none specific | Not yet mapped | Add UGT1A1 `*27` fixture. |
| `testUgt1a1HG00436` | 1446 | UGT1A1 phased no-call | matcher no-call generic tests | Partial Rust coverage | Add HG00436-inspired no-call fixture. |
| `testUgt1a1s1s80s27s60s28MissingPhased` | 1459 | UGT1A1 phased no-call preserves phased flag | matcher/report phasing helpers | Partial Rust coverage | Add exact phased no-call report fixture. |
| `testUgt1a1s1s60s80s6phased` | 1475 | UGT1A1 phased no-call | matcher/report phasing helpers | Partial Rust coverage | Add exact phased no-call fixture. |
| `testUgt1a1s1s60s80s28s6phased` | 1490 | UGT1A1 phased no-call with repeat | matcher/report phasing helpers | Partial Rust coverage | Add exact phased repeat no-call fixture. |
| `testUgt1a1s1s37s80s60phased` | 1506 | UGT1A1 phased `*1/*80+*37` | none specific | Not yet mapped | Add UGT1A1 `*37` fixture. |
| `testCyp3a5Missing3Message` | 1524 | CYP3A5 missing `*3` position and tacrolimus message | CYP3A5 pipeline fixture; matcher/report missing-position helpers | Partial Rust coverage | Add exact missing-`rs776746` message fixture. |
| `testCyp3a5v1` | 1550 | CYP3A5 `*1/*3` | CYP3A5 CLI/pipeline fixture covers `*1/*2`, not this call | Partial Rust coverage | Add CYP3A5 `*1/*3` fixture. |
| `testCyp3a5v2` | 1563 | CYP3A5 `*3/*9` | CYP3A5 matcher generic coverage | Partial Rust coverage | Add CYP3A5 `*3/*9` fixture. |
| `testCyp3a5v3` | 1577 | CYP3A5 `*3/*3` | CYP3A5 matcher generic coverage | Partial Rust coverage | Add CYP3A5 hom `*3` fixture. |
| `testCyp3a5v4` | 1590 | CYP3A5 `*1/*3` duplicate orientation case | CYP3A5 matcher generic coverage | Partial Rust coverage | Add orientation-specific fixture. |
| `testCyp3a5v5` | 1603 | CYP3A5 `*3/*9` duplicate orientation case | CYP3A5 matcher generic coverage | Partial Rust coverage | Add orientation-specific fixture. |
| `testHlab` | 1617 | HLA-B outside diplotype plus CYP2C9 recommendations | HLA-B phenotype/report tests | Partial Rust coverage | Add HLA-B outside-call pipeline fixture. |
| `testSingleHlabAllele` | 1653 | Single HLA-B allele outside call | HLA-B phenotype tests | Partial Rust coverage | Add single-allele HLA-B pipeline fixture. |
| `testHlabPhenotype` | 1685 | HLA-B phenotype outside call | HLA-B phenotype tests | Partial Rust coverage | Add phenotype-only HLA-B pipeline fixture. |
| `testRecommendationExamples` | 1732 | Multi-drug recommendation example report | broad report HTML/JSON tests | Partial Rust coverage | Promote as a high-value end-to-end parity fixture. |
| `testCyp2c9star61` | 1785 | CYP2C9 `*1/*61` | none specific | Not yet mapped | Add CYP2C9 `*61` fixture. |
| `testCyp2c9star1Hom` | 1798 | CYP2C9 `*1/*1` and NSAID recommendations | report recommendation tests | Partial Rust coverage | Add CYP2C9 reference fixture. |
| `testCyp2b6star1star34` | 1821 | CYP2B6 top-candidate `*1/*34` | CYP2B6 matcher generated fixtures, not this one | Partial Rust coverage | Add CYP2B6 `*1/*34` pipeline fixture. |
| `testCyp2b6star1star34AllMatch` | 1845 | CYP2B6 all-matches `*1/*34`, `*33/*36` | matcher all-candidate helpers | Partial Rust coverage | Add all-matches pipeline fixture. |
| `testMtrnr1` | 1867 | MT-RNR1 outside call with aminoglycoside recommendations | outside-call parser coverage | Partial Rust coverage | Add MT-RNR1 outside-call pipeline fixture. |
| `testIfnl3` | 1889 | IFNL3 reference single-position gene | none specific | Not yet mapped | Add IFNL3 pipeline fixture. |
| `testCyp3a4` | 1905 | CYP3A4 `*8/*8` quetiapine | none specific | Not yet mapped | Add CYP3A4 fixture. |
| `testPartialCall` | 1921 | CYP2C19 partial call with combo/partial mode | matcher partial/combination tests | Partial Rust coverage | Add partial-call pipeline fixture when combo CLI path is ported. |
| `testPartialCallInTwoGene` | 1938 | CYP2C19 partial call with CYP2D6 outside call | matcher partial + outside-call coverage | Partial Rust coverage | Add two-gene partial pipeline fixture. |
| `testOutsideCallCollision` | 1962 | Outside call overrides same-gene VCF call and warning | outside-call parser/report tests | Partial Rust coverage | Add same-gene collision pipeline fixture. |
| `outsideCallCollision2Files` | 1985 | Two outside-call files merge/normalize CYP4F2 and warn once | outside-call parser tests | Partial Rust coverage | Add multi-file outside-call fixture. |
| `testOutsideCallDiplotypeNormalization` | 2019 | Outside-call diplotype order normalization | outside-call/model tests | Partial Rust coverage | Add pipeline-level normalization fixture. |
| `testOutsideCallPhenotypeOverridesDiplotype` | 2037 | Outside phenotype overrides known diplotype phenotype | outside-call phenotype tests | Partial Rust coverage | Add override recommendation fixture. |
| `testOutsideCallActivityScore` | 2072 | Outside activity score without diplotype | outside-call parser tests | Partial Rust coverage | Add activity-score-only pipeline fixture. |
| `testOutsideCallActivityScoreAndPhenotype` | 2101 | Outside diplotype plus phenotype plus activity score override | outside-call parser/report tests | Partial Rust coverage | Add combined outside-call fixture. |
| `testWarfarinDpwg` | 2132 | CYP2C9/VKORC1 warfarin DPWG matching | recommendation/report tests | Partial Rust coverage | Add warfarin DPWG pipeline fixture. |
| `testWarfarinMissingRs12777823` | 2150 | Warfarin missing extra position highlighted in HTML | report extra-position helpers | Partial Rust coverage | Add missing-extra-position HTML fixture. |
| `testOutsideSinglePositionCalls` | 2193 | IFNL3 and CYP4F2 outside calls | outside-call parser tests | Partial Rust coverage | Add single-position outside-call pipeline fixture. |
| `testDiplotypeOverrideRecommendation` | 2224 | CYP2C9 diplotype-specific recommendation overrides phenotype | recommendation matching tests | Partial Rust coverage | Add phenytoin override pipeline fixture. |
| `testDuplicateEntrySecondIsBad` | 2245 | Duplicate VCF row warning when second row invalid | VCF duplicate handling tests | Partial Rust coverage | Add NUDT15 duplicate-row pipeline fixture. |
| `testDuplicateEntryFirstIsBad` | 2266 | REF mismatch plus duplicate warning order | VCF REF mismatch and duplicate tests | Partial Rust coverage | Add exact warning-order pipeline fixture. |
| `phaseSetDpyd` | 2296 | DPYD unphased/phased/phase-set behavior | DPYD matcher phase-set tests; multi-gene pipeline slice | Partial Rust coverage | Add DPYD phase-set end-to-end fixture. |
| `testNoCallOutsideCall` | 2350 | No-call outside calls for HLA-B and NAT2 | outside-call parser/model tests | Partial Rust coverage | Add no-call outside-call pipeline fixture. |
| `phaseSetCyp2C9` | 2368 | CYP2C9 unphased/phased/combo/phase-set behavior and HTML calls | matcher phase/combination tests | Partial Rust coverage | Add CYP2C9 phase-set pipeline fixture. |

## Definition Tests

| Java test class | Java tests | Rust owner | Status | Notes / next action |
| --- | ---: | --- | --- | --- |
| `definition/AssemblyMapTest` | 1 | `definition.rs` | Not yet mapped | Assembly-map behavior needs explicit Rust fixture/test owner. |
| `definition/PhenotypeMapTest` | 5 | `phenotype.rs` | Partial Rust coverage | Rust loads Java phenotype directory and exercises lookup behavior. |
| `definition/model/NamedAlleleTest` | 3 | `definition.rs`, `matcher.rs`, `pipeline.rs` | Partial Rust coverage | Named allele JSON/model and matcher pattern behavior are partially covered; map Java methods. |

## Haplotype And Matcher Tests

| Java test class | Java tests | Rust owner | Status | Notes / next action |
| --- | ---: | --- | --- | --- |
| `haplotype/CombinationUtilTest` | 5 | `matcher.rs` | Partial Rust coverage | Combination viability/name behavior has Rust tests; map all Java utility cases. |
| `haplotype/DefinitionReaderTest` | 3 | `definition.rs` | Partial Rust coverage | Definition loading/indexing has Rust tests. |
| `haplotype/DiplotypeMatcherTest` | 6 | `matcher.rs` | Partial Rust coverage | Rust has diplotype pair tests based on Java intent. |
| `haplotype/IupacTest` | 1 | `pipeline.rs`, `matcher.rs` | Partial Rust coverage | IUPAC display expansion is covered for matcher HTML; core IUPAC API parity still needs a direct owner if required. |
| `haplotype/NamedAlleleMatcherCftrTest` | 2 | `matcher.rs` | Not yet mapped | CFTR gene fixture coverage needs method-level mapping. |
| `haplotype/NamedAlleleMatcherCyp2c19Test` | 18 | `matcher.rs`, `report.rs` | Partial Rust coverage | CYP2C19 appears in matcher/report fixtures, but full Java fixture sweep is not mapped. |
| `haplotype/NamedAlleleMatcherCyp2c9Test` | 8 | `matcher.rs` | Not yet mapped | Needs fixture-by-fixture Rust ownership. |
| `haplotype/NamedAlleleMatcherCyp3a5Test` | 5 | `matcher.rs`, `pipeline.rs` | Partial Rust coverage | Java haplotyper CYP3A5 fixture is used by Rust matcher/pipeline tests. |
| `haplotype/NamedAlleleMatcherCyp4f2Test` | 2 | `matcher.rs` | Not yet mapped | Needs fixture-by-fixture Rust ownership. |
| `haplotype/NamedAlleleMatcherIfnl3Test` | 3 | `matcher.rs` | Not yet mapped | Needs fixture-by-fixture Rust ownership. |
| `haplotype/NamedAlleleMatcherNat2Test` | 2 | `matcher.rs` | Not yet mapped | Needs fixture-by-fixture Rust ownership. |
| `haplotype/NamedAlleleMatcherRyr1Test` | 3 | `matcher.rs`, `report.rs` | Partial Rust coverage | RYR1 lowest-function coverage exists; specific matcher fixture mapping still needed. |
| `haplotype/NamedAlleleMatcherSlco1b1Test` | 6 | `matcher.rs`, `report.rs` | Partial Rust coverage | SLCO1B1 reporter fallback has Rust coverage; matcher fixture sweep not mapped. |
| `haplotype/NamedAlleleMatcherTest` | 42 | `matcher.rs`, `pipeline.rs` | Partial Rust coverage | Core matcher logic has many Rust tests, but this Java class needs the next detailed method-level map. |
| `haplotype/NamedAlleleMatcherTpmtTest` | 6 | `matcher.rs` | Not yet mapped | Needs fixture-by-fixture Rust ownership. |
| `haplotype/NamedAlleleMatcherUgt1a1Test` | 14 | `matcher.rs`, `phenotype.rs` | Not yet mapped | Needs fixture-by-fixture Rust ownership. |
| `haplotype/NamedAlleleMatcherVkorc1Test` | 3 | `matcher.rs` | Not yet mapped | Needs fixture-by-fixture Rust ownership. |
| `haplotype/SampleAlleleTest` | 1 | `matcher.rs`, `vcf.rs` | Partial Rust coverage | `SampleAllele` VCF-call/allele display behavior is partially covered. |
| `haplotype/VcfReaderTest` | 13 | `vcf.rs` | Partial Rust coverage | Rust covers VCF extension checks, filters, phasing, phase sets, allele order, AD warnings, BGZF, and error paths. |
| `haplotype/VcfSampleReaderTest` | 1 | `vcf.rs`, `pipeline.rs` | Partial Rust coverage | Rust covers multisample sample selection and sample metadata ingestion. |
| `haplotype/model/CombinationMatchTest` | 1 | `matcher.rs` | Partial Rust coverage | Combination match behavior is represented in matcher tests. |
| `haplotype/model/DiplotypeMatchTest` | 2 | `matcher.rs`, `report.rs` | Partial Rust coverage | Diplotype ordering/reporting has Rust tests; map Java methods. |
| `haplotype/model/HaplotypeMatchTest` | 1 | `matcher.rs` | Partial Rust coverage | Haplotype match result shape is represented indirectly. |

### `haplotype/NamedAlleleMatcherTest` Method Map

This table maps the 42 Java methods in `NamedAlleleMatcherTest.java`. It is a planning artifact, not a parity claim: `Partial Rust coverage` means the Rust tests cover nearby behavior, but the exact Java fixture or synthetic `TestVcfBuilder` scenario still needs a promoted Rust test or Java-vs-Rust parity case.

| Java method | Line | Java focus / fixture | Rust owner or closest test | Status | Next action |
| --- | ---: | --- | --- | --- | --- |
| `testCall` | 183 | CYP3A5 `haplotyper.vcf` / `haplotyper.json`, save matcher JSON/HTML | `matcher::tests::call_standard_gene_returns_exact_diplotype_for_java_haplotyper_fixture`; `pipeline::tests` matcher output tests | Partial Rust coverage | Add direct save-results parity or promote CLI intermediate comparison. |
| `testCallDiplotypePath` | 208 | CYP3A5 match-data, permutations, diplotype pair | `matcher::tests::builds_match_data_for_java_haplotyper_fixture` | Direct Rust coverage | Promote to parity after Java gate is runnable. |
| `sortException` | 256 | CYP2D6 sort-contract regression and 34 VCF warnings from `NamedAlleleMatcher-sortError.vcf` | `matcher::tests::cyp2d6_sort_exception_fixture_produces_many_diplotypes_like_java_named_allele_matcher` | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testMismatchedRefAlleleWarnings` | 280 | REF mismatch warnings from `NamedAlleleMatcher-mismatchedRefAllele.*` | `matcher::tests::definition_aware_vcf_ref_mismatch_warnings_match_java_named_allele_matcher`; `pipeline::run_reporter_from_vcf` uses the same definition-aware VCF adapter | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testCyp2d6` | 310 | CYP2D6 include/exclude toggle | none | Not yet mapped | Decide Rust CLI/matcher option owner for CYP2D6 inclusion. |
| `testWobbleScoring` | 331 | CYP2D6 wobble score keeps `*1/*2` and `*1/*12` | `matcher::tests::wobble_alleles_do_not_score_when_all_supporting_sequences_are_reference` | Partial Rust coverage | Add Java CYP2D6 synthetic fixture equivalent. |
| `testCombinationBaseline` | 359 | `NamedAlleleMatcher-combinationBaseline.vcf` -> `*1/*1` | `matcher::tests::computes_combination_baseline_from_java_fixture` | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testCombinationHomozygous` | 375 | CYP2B6 synthetic homozygous combination `[*2 + *5]/[*2 + *5]` | `matcher::tests::computes_combination_matches_and_partials_like_java_combination_matcher` | Partial Rust coverage | Add CYP2B6 synthetic fixture or generated VCF helper. |
| `testCombinationPhased` | 397 | `NamedAlleleMatcher-combinationPhased.vcf` | `matcher::tests::computes_combination_phased_from_java_fixture` | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testCombinationUnphased` | 418 | `NamedAlleleMatcher-combinationUnphased.vcf` eight diplotypes | `matcher::tests::computes_combination_unphased_from_java_fixture` | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testPartialWithCombination` | 450 | `NamedAlleleMatcher-partialWithCombination.vcf` | `matcher::tests::computes_partial_with_combination_from_java_fixture`; `call_standard_gene_uses_combination_fallback_when_enabled_like_java` | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testPartialWithCombinationUnphased` | 481 | NAT2 synthetic unphased partial combination | none | Not yet mapped | Add NAT2 generated VCF helper or fixture. |
| `testPartial` | 510 | `NamedAlleleMatcher-partial.vcf` | `matcher::tests::computes_partial_from_java_fixture` | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testPartialUnphasedWithSingleHapMatch` | 538 | CYP2B6 synthetic partial when combination matching yields <= 1 haplotype | none | Not yet mapped | Add generated VCF helper and expected names. |
| `testPartial2Phased` | 572 | CYP2C19 synthetic phased partial and top-candidate reduction | none | Not yet mapped | Add CYP2C19 generated fixture. |
| `testCombo3` | 604 | CYP2B6 phased combo, matcher HTML/report HTML side effects | none | Not yet mapped | Add matcher result fixture first; report HTML parity can follow. |
| `testCombo3WithPhaseSet` | 647 | CYP2B6 phase-set combination ordering | none | Not yet mapped | Add phase-set generated fixture. |
| `testCombinationWithShellPhased` | 689 | CYP2B6 shell allele with combination | none | Not yet mapped | Add shell/combination fixture. |
| `testShellWithPartial` | 727 | CYP2B6 shell allele eliminates child haplotypes with partial | none | Not yet mapped | Add shell/partial fixture. |
| `testCombinationWithShellMissing` | 753 | CYP2B6 shell allele with missing tag position | none | Not yet mapped | Add shell/missing fixture. |
| `testPartialReferenceUnphased` | 784 | `NamedAlleleMatcher-partialReferenceUnphased.vcf` | `matcher::tests::partial_reference_unphased_fixture_matches_java_named_allele_matcher` | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testPartialReferencePhased` | 815 | `NamedAlleleMatcher-partialReferencePhased.vcf` | `matcher::tests::partial_reference_phased_fixture_matches_java_named_allele_matcher` | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testPartialReferenceDouble` | 841 | `NamedAlleleMatcher-partialReferenceDouble.vcf` | `matcher::tests::partial_reference_double_fixture_matches_java_named_allele_matcher` | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testDpydEffectivelyPhased` | 867 | DPYD synthetic effectively phased `Reference/c.62G>A` | `matcher::tests::dpyd_effectively_phased_fixture_matches_java_named_allele_matcher` | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testDpydEffectivelyPhased2` | 898 | DPYD synthetic multi-component effectively phased combination | `matcher::tests::dpyd_effectively_phased_combination_fixture_matches_java_named_allele_matcher` | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testDpydPhased` | 940 | DPYD phased c.62/c.3067 HapB3 merge | `matcher::tests::dpyd_phased_fixture_matches_java_named_allele_matcher` | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testDpydUnPhased` | 973 | DPYD unphased haplotype fallback for c.62/c.3067 | `matcher::tests::dpyd_unphased_fixture_matches_java_named_allele_matcher` | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testDpydPhasedDouble` | 1003 | DPYD phased double combination including HapB3 | `matcher::tests::dpyd_phased_double_fixture_matches_java_named_allele_matcher` | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testDpydHomozygous` | 1042 | DPYD phased homozygous combination | `matcher::tests::dpyd_homozygous_fixture_matches_java_named_allele_matcher` | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testDpydEffectivelyPhasedCombination` | 1075 | DPYD NA18973-inspired effectively phased combination | `matcher::tests::dpyd_effectively_phased_combination_missing_fixture_matches_java_named_allele_matcher` | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testDpydUnphasedHomozygousNoFunction` | 1113 | DPYD unphased no-function haplotype fallback | `matcher::tests::dpyd_unphased_homozygous_no_function_fixture_matches_java_named_allele_matcher` | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testDiplotypeMatcher` | 1146 | DPYD `NamedAlleleMatcher-diplotypeMatcher.vcf` sort stability loop | `matcher::tests::dpyd_diplotype_matcher_fixture_repeats_like_java_named_allele_matcher` | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testUnknownAltMultisample` | 1161 | RYR1 unknown ALT multisample selected sample has no warnings | `matcher::tests::unknown_alt_multisample_selected_sample_has_no_warnings_like_java_named_allele_matcher`; `vcf::tests::accepts_ft_format_number_unknown_like_java_vcf_reader` | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testPermutationGeneration` | 1179 | CYP2C19 phased/unphased/mixed permutations all call `*2/*17` | `matcher::tests::cyp2c19_permutation_generation_fixture_matches_java_named_allele_matcher` | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testPartialMissingAllele` | 1213 | CYP2B6 missing allele without combinations yields no matches | `matcher::tests::call_standard_gene_returns_no_call_for_partial_missing_without_combinations_like_java` | Direct Rust coverage | Add exact VCF warning assertion if needed. |
| `testPartialMissingAllele_combination1` | 1233 | CYP2B6 missing allele with combinations yields `*2/g.40991369?` | `matcher::tests::cyp2b6_partial_missing_allele_combination_fixture_matches_java_named_allele_matcher` | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `testPartialMissingAllele_combination_phased` | 1256 | Phased missing allele with combinations yields `*2/g.40991369?` | `matcher::tests::cyp2b6_phased_partial_missing_allele_combination_fixture_matches_java_named_allele_matcher` | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `unphasedPrioritySameScore` | 1279 | NAT2 same-score unphased priority with exemption warning | `matcher::tests::nat2_unphased_priority_same_score_fixture_matches_java_named_allele_matcher` | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |
| `unphasedPriorityDifferentScore` | 1321 | NAT2 different-score unphased priority and serialized matcher outputs | `matcher::tests::nat2_unphased_priority_different_score_fixture_matches_java_named_allele_matcher`; generic matcher JSON/HTML warning rendering is covered in `pipeline::tests` | Direct Rust matcher coverage | Promote save-results JSON/HTML to Java-vs-Rust fixture later. |
| `requiredPosition` | 1369 | NAT2 required-position exemption gates reference call | `matcher::tests::nat2_required_position_fixture_matches_java_named_allele_matcher`; pipeline/report warning tests cover generic serialization | Direct Rust matcher coverage | Promote save-results JSON/HTML to Java-vs-Rust fixture later. |
| `testNat2Combination` | 1417 | NAT2 phased combination and report HTML side effects | `matcher::tests::nat2_combination_fixture_matches_java_named_allele_matcher`; report HTML side effects remain tracked for Java-vs-Rust promotion | Direct Rust matcher coverage | Promote report HTML output to Java-vs-Rust fixture later. |
| `testWobbleScoringWithMultipleSequenceMatches` | 1470 | CYP2B6 wobble with multiple sequences preserves `*6/*18`, `*9/*18` | `matcher::tests::cyp2b6_wobble_scoring_multiple_sequence_fixture_matches_java_named_allele_matcher` | Direct Rust coverage | Promote to Java-vs-Rust fixture later. |

## Phenotype Tests

| Java test class | Java tests | Rust owner | Status | Notes / next action |
| --- | ---: | --- | --- | --- |
| `phenotype/OutsideCallParserTest` | 7 | `phenotype.rs` | Partial Rust coverage | Outside call parsing, comments, errors, and set behavior have Rust tests. |
| `phenotype/PhenotypeUtilsTest` | 1 | `phenotype.rs` | Partial Rust coverage | Phenotype normalization has Rust tests. |
| `phenotype/PhenotyperTest` | 7 | `phenotype.rs`, `report.rs` | Partial Rust coverage | Diplotype/activity lookup and report feeding have Rust tests; map fixture JSON methods. |
| `phenotype/model/OutsideCallTest` | 5 | `phenotype.rs`, `report.rs` | Partial Rust coverage | Outside call annotation and mismatch behavior have Rust tests. |

## Reporter Tests

| Java test class | Java tests | Rust owner | Status | Notes / next action |
| --- | ---: | --- | --- | --- |
| `reporter/PgkbGuidelineCollectionTest` | 1 | `report.rs` | Partial Rust coverage | Guideline loading/indexing has Rust tests. |
| `reporter/RecommendationUtilsTest` | 1 | `report.rs` | Partial Rust coverage | Recommendation lookup/matching has Rust tests. |
| `reporter/ReporterTest` | 3 | `report.rs`, `pipeline.rs` | Partial Rust coverage | Report writer and end-to-end reporter output have Rust coverage; exact Java fixtures not fully mapped. |
| `reporter/VariantReportFactoryTest` | 1 | `report.rs` | Partial Rust coverage | Variant report construction/sorting/call helpers have Rust tests. |
| `reporter/caller/Cyp2d6CopyNumberCallerTest` | 4 | `report.rs`, `phenotype.rs` | Not yet mapped | Copy-number caller behavior needs a direct Rust owner or documented out-of-scope decision. |
| `reporter/caller/LowestFunctionGeneCallerTest` | 6 | `matcher.rs`, `phenotype.rs`, `report.rs` | Partial Rust coverage | DPYD/RYR1 lowest-function logic has Rust tests. |
| `reporter/format/CallsOnlyFormatTest` | 5 | `report.rs` | Partial Rust coverage | Calls-only TSV header, debug columns, metadata, writing, and validation have Rust tests. |
| `reporter/handlebars/ReportHelpersTest` | 1 | `report.rs` | Partial Rust coverage | Many HTML helper behaviors are covered in Rust report HTML tests; map the Java helper test explicitly. |
| `reporter/model/VariantReportTest` | 1 | `report.rs` | Partial Rust coverage | Variant report ordering/call helpers have Rust tests. |
| `reporter/model/pgkb/RecommendationAnnotationTest` | 1 | `report.rs` | Partial Rust coverage | Recommendation genotype matching has Rust tests. |
| `reporter/model/result/DiplotypeTest` | 5 | `report.rs` | Partial Rust coverage | Report diplotype sorting/serialization/combination/outside-call behavior has Rust tests. |
| `reporter/model/result/DrugLinkTest` | 1 | `report.rs` | Partial Rust coverage | Drug-link/backlink behavior has Rust tests. |
| `reporter/model/result/GeneReportTest` | 4 | `report.rs` | Partial Rust coverage | Gene-report JSON/source/merge behavior has Rust tests. |
| `reporter/model/result/GenotypeTest` | 4 | `report.rs` | Partial Rust coverage | Recommendation genotype ordering and report-gene serialization have Rust tests. |
| `reporter/model/result/HaplotypeTest` | 5 | `report.rs` | Partial Rust coverage | Haplotype serialization/comparison behavior is only indirectly covered; map Java methods. |

## Utility And Stats Tests

| Java test class | Java tests | Rust owner | Status | Notes / next action |
| --- | ---: | --- | --- | --- |
| `stats/CalcAlleleFrequenciesTest` | 2 | n/a | Not yet mapped | Stats command has not been assigned a Rust owner. |
| `util/ActivityUtilsTest` | 1 | `phenotype.rs` | Partial Rust coverage | Activity-score normalization has Rust tests. |
| `util/ChrNameComparatorTest` | 1 | `common/chromosome.rs`, `vcf.rs`, `definition.rs` | Not yet mapped | Chromosome ordering needs a direct Rust test if exposed behavior matters. |
| `util/CliUtilsTest` | 1 | `cli.rs` | Partial Rust coverage | Base filename/path helper behavior has Rust tests. |
| `util/DataSerializerTest` | 2 | `definition.rs`, `report.rs` | Partial Rust coverage | Definition exemption JSON/TSV loading has Rust tests; serializer behavior needs method map. |
| `util/HaplotypeNameComparatorTest` | 5 | `matcher.rs`, `report.rs` | Partial Rust coverage | Haplotype/diplotype sorting has Rust tests; direct comparator parity still needs method map. |
| `util/VariantUtilsTest` | 2 | `report.rs`, `definition.rs` | Partial Rust coverage | Variant call helpers/HGVS-like names have Rust tests. |

## Next Inventory Work

1. Expand `PipelineTest` method-by-method after the Java reference gate can run or after a local Java artifact is supplied.
2. Promote the covered `NamedAlleleMatcherTest` matcher/report output side effects to Java-vs-Rust parity fixtures once the Java reference gate is available.
3. For each `Not yet mapped` gene-specific matcher class, list the VCF/JSON fixtures and decide whether to port fixture assertions into Rust unit tests or promote them into Java-vs-Rust CLI parity.
4. Once Java Gradle is unblocked, record exact Java pass/fail/skip counts here and update `TODO.md`.
