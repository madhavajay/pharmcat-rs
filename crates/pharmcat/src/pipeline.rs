//! Pipeline orchestration helpers.

use std::{
    collections::{BTreeMap, BTreeSet},
    fs, io,
    path::{Path, PathBuf},
    time::{SystemTime, UNIX_EPOCH},
};

use crate::{
    cli::{CliError, PharmcatCliConfig, PipelineOutputPlan, pipeline_output_plan},
    definition::{
        DefinitionFile, DefinitionLoadError, DefinitionReader, NamedAllele, VariantLocus,
    },
    matcher::{
        DiplotypeMatch, DpydHapB3Warning, GeneCallKind, GeneCallResult, GeneCallWarning,
        HaplotypeMatch, MatchData, MatchError, call_dpyd_lowest_function_gene,
        call_ryr1_lowest_function_gene, call_standard_gene_with_exemption,
    },
    phenotype::{
        OutsideCallError, OutsideCallValidation, PhenotypeLoadError, PhenotypeMap,
        parse_outside_calls_file,
    },
    report::{
        CallsOnlyTsvOptions, GuidanceLoadError, HtmlReportOptions, MessageCatalog,
        PgkbGuidelineCollection, ReportContext, ReportContextFromMatcherError, ReportGene,
        ReportGeneFromOutsideCallError, ReportGeneFromStandardCallError, write_calls_only_tsv,
        write_report_html, write_report_json,
    },
    vcf::{
        ReadVcfError, SampleAlleleSummary, VcfWarnings, allele_calls_for_locations_with_genes,
        read_record_summaries,
    },
};

/// Java `Pipeline.Mode`.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum PipelineMode {
    /// Default CLI mode.
    Cli,
    /// Batch mode.
    Batch,
    /// Test mode omits volatile version/timestamp output where possible.
    Test,
}

/// Java `PipelineResult.Status`.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum PipelineStatus {
    /// Nothing was run.
    Noop,
    /// Pipeline completed.
    Success,
    /// Pipeline failed and was converted to a per-sample/batch failure.
    Failure,
}

/// Java `PipelineResult`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PipelineResult {
    /// Result status.
    pub status: PipelineStatus,
    /// Java pipeline basename.
    pub basename: String,
    /// Optional sample ID.
    pub sample_id: Option<String>,
}

impl PipelineResult {
    /// Builds a Java-style pipeline result.
    pub fn new(
        status: PipelineStatus,
        basename: impl Into<String>,
        sample_id: Option<impl Into<String>>,
    ) -> Self {
        Self {
            status,
            basename: basename.into(),
            sample_id: sample_id.map(Into::into),
        }
    }
}

/// Execution metadata derived before running pipeline stages.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PipelineRunPlan {
    /// Parsed output paths and names.
    pub outputs: PipelineOutputPlan,
    /// Java `Pipeline.getInputDescription`.
    pub input_description: String,
    /// Whether Java would use batch display formatting.
    pub batch_display_mode: bool,
    /// Java-style starting line, when batch display mode is active.
    pub starting_line: Option<String>,
    /// Files Java deletes when `-del` is active.
    pub intermediate_files_to_delete: Vec<PathBuf>,
    /// Java batch/per-sample error file path.
    pub error_file: PathBuf,
    /// User-visible save messages Java would print for completed stage writes.
    pub save_messages: Vec<String>,
    /// Java batch finished block without memory-usage text.
    pub finished_block: Option<String>,
}

/// Options for the first Rust VCF-to-reporter pipeline helper.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ReporterPipelineOptions {
    /// Whether standard genes should use combination matching.
    pub include_combinations: bool,
    /// Java phenotyper outside-call TSV inputs.
    pub outside_call_files: Vec<PathBuf>,
    /// Java reporter message catalog.
    pub message_catalog: Option<MessageCatalog>,
    /// HTML rendering options.
    pub html: HtmlReportOptions,
    /// Calls-only TSV rendering options.
    pub calls_only_tsv: CallsOnlyTsvOptions,
}

impl Default for ReporterPipelineOptions {
    fn default() -> Self {
        Self {
            include_combinations: true,
            outside_call_files: Vec::new(),
            message_catalog: None,
            html: HtmlReportOptions::default(),
            calls_only_tsv: CallsOnlyTsvOptions::default(),
        }
    }
}

/// Result of the first Rust VCF-to-reporter pipeline helper.
#[derive(Clone, Debug, PartialEq)]
pub struct ReporterPipelineRun {
    /// Reporter context generated from matcher calls.
    pub context: ReportContext,
    /// Gene call results generated from the VCF and definitions.
    pub gene_call_results: Vec<GeneCallResult>,
    /// VCF warnings collected while reading the input.
    pub vcf_warnings: VcfWarnings,
    /// Reporter output files written by the helper.
    pub written_outputs: Vec<PathBuf>,
}

/// Resource paths required by the Rust CLI pipeline slice.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PharmcatResourcePaths {
    /// Allele definition directory, usually `definition/alleles`.
    pub definitions_dir: PathBuf,
    /// Phenotype JSON directory.
    pub phenotype_dir: PathBuf,
    /// Prescribing guidance JSON file.
    pub prescribing_guidance: PathBuf,
    /// Reporter messages JSON file.
    pub reporter_messages: PathBuf,
}

impl PharmcatResourcePaths {
    /// Builds paths from a Java PharmCAT resource root.
    pub fn from_resource_root(root: impl AsRef<Path>) -> Self {
        let root = root.as_ref();
        Self {
            definitions_dir: root.join("definition").join("alleles"),
            phenotype_dir: root.join("phenotype"),
            prescribing_guidance: root
                .join("reporter")
                .join(PgkbGuidelineCollection::PRESCRIBING_GUIDANCE_FILE_NAME),
            reporter_messages: root
                .join("reporter")
                .join(MessageCatalog::MESSAGES_JSON_FILE_NAME),
        }
    }
}

/// Options for running the currently ported CLI pipeline surface.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CliPipelineOptions {
    /// Resource files used to load definitions, phenotypes, and reporter guidance.
    pub resources: PharmcatResourcePaths,
    /// Pipeline mode used for Java-style plan metadata.
    pub mode: PipelineMode,
}

/// Error from the first Rust VCF-to-reporter pipeline helper.
#[derive(Debug)]
pub enum ReporterPipelineError {
    /// Failed to read VCF records.
    Vcf(ReadVcfError),
    /// Gene matching failed.
    Match {
        /// Gene being matched.
        gene: String,
        /// Matching error.
        source: MatchError,
    },
    /// Failed to build report context.
    ReportContext(ReportContextFromMatcherError),
    /// Failed to parse Java phenotyper outside-call TSV input.
    OutsideCall(OutsideCallError),
    /// Failed to convert an outside call into a reporter gene.
    OutsideReport(ReportGeneFromOutsideCallError),
    /// Failed to write report output.
    ReportOutput(GuidanceLoadError),
    /// Failed to create output directories.
    Io(io::Error),
}

/// Error from the currently ported CLI pipeline surface.
#[derive(Debug)]
pub enum CliPipelineError {
    /// Parsed CLI configuration is not covered by the current Rust pipeline slice.
    Unsupported(String),
    /// Failed to derive Java-style output paths.
    Cli(CliError),
    /// Failed to load allele definitions.
    Definitions(DefinitionLoadError),
    /// Failed to load phenotype maps.
    Phenotypes(PhenotypeLoadError),
    /// Failed to load prescribing guidance.
    Guidance(GuidanceLoadError),
    /// Failed to serialize JSON output.
    Json(serde_json::Error),
    /// Failed while running the reporter pipeline helper.
    Reporter(ReporterPipelineError),
}

impl std::fmt::Display for ReporterPipelineError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Vcf(error) => write!(f, "{error}"),
            Self::Match { gene, source } => write!(f, "Failed to match {gene}: {source}"),
            Self::ReportContext(error) => write!(f, "{error}"),
            Self::OutsideCall(error) => write!(f, "{error}"),
            Self::OutsideReport(error) => write!(f, "{error}"),
            Self::ReportOutput(error) => write!(f, "{error}"),
            Self::Io(error) => write!(f, "{error}"),
        }
    }
}

impl std::error::Error for ReporterPipelineError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Vcf(error) => Some(error),
            Self::Match { source, .. } => Some(source),
            Self::ReportContext(error) => Some(error),
            Self::OutsideCall(error) => Some(error),
            Self::OutsideReport(error) => Some(error),
            Self::ReportOutput(error) => Some(error),
            Self::Io(error) => Some(error),
        }
    }
}

impl std::fmt::Display for CliPipelineError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Unsupported(message) => f.write_str(message),
            Self::Cli(error) => write!(f, "{error}"),
            Self::Definitions(error) => write!(f, "{error}"),
            Self::Phenotypes(error) => write!(f, "{error}"),
            Self::Guidance(error) => write!(f, "{error}"),
            Self::Json(error) => write!(f, "{error}"),
            Self::Reporter(error) => write!(f, "{error}"),
        }
    }
}

impl std::error::Error for CliPipelineError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Unsupported(_) => None,
            Self::Cli(error) => Some(error),
            Self::Definitions(error) => Some(error),
            Self::Phenotypes(error) => Some(error),
            Self::Guidance(error) => Some(error),
            Self::Json(error) => Some(error),
            Self::Reporter(error) => Some(error),
        }
    }
}

impl From<ReadVcfError> for ReporterPipelineError {
    fn from(error: ReadVcfError) -> Self {
        Self::Vcf(error)
    }
}

impl From<ReportContextFromMatcherError> for ReporterPipelineError {
    fn from(error: ReportContextFromMatcherError) -> Self {
        Self::ReportContext(error)
    }
}

impl From<OutsideCallError> for ReporterPipelineError {
    fn from(error: OutsideCallError) -> Self {
        Self::OutsideCall(error)
    }
}

impl From<ReportGeneFromOutsideCallError> for ReporterPipelineError {
    fn from(error: ReportGeneFromOutsideCallError) -> Self {
        Self::OutsideReport(error)
    }
}

impl From<ReportGeneFromStandardCallError> for ReporterPipelineError {
    fn from(error: ReportGeneFromStandardCallError) -> Self {
        Self::ReportContext(ReportContextFromMatcherError::from(error))
    }
}

impl From<GuidanceLoadError> for ReporterPipelineError {
    fn from(error: GuidanceLoadError) -> Self {
        Self::ReportOutput(error)
    }
}

impl From<io::Error> for ReporterPipelineError {
    fn from(error: io::Error) -> Self {
        Self::Io(error)
    }
}

impl From<CliError> for CliPipelineError {
    fn from(error: CliError) -> Self {
        Self::Cli(error)
    }
}

impl From<DefinitionLoadError> for CliPipelineError {
    fn from(error: DefinitionLoadError) -> Self {
        Self::Definitions(error)
    }
}

impl From<PhenotypeLoadError> for CliPipelineError {
    fn from(error: PhenotypeLoadError) -> Self {
        Self::Phenotypes(error)
    }
}

impl From<GuidanceLoadError> for CliPipelineError {
    fn from(error: GuidanceLoadError) -> Self {
        Self::Guidance(error)
    }
}

impl From<serde_json::Error> for CliPipelineError {
    fn from(error: serde_json::Error) -> Self {
        Self::Json(error)
    }
}

impl From<ReporterPipelineError> for CliPipelineError {
    fn from(error: ReporterPipelineError) -> Self {
        Self::Reporter(error)
    }
}

impl PipelineRunPlan {
    /// Builds the Java pipeline pre-run metadata from a parsed CLI config.
    pub fn from_cli(
        config: &PharmcatCliConfig,
        sample_id: Option<&str>,
        single_sample: bool,
        mode: PipelineMode,
        display_count: Option<&str>,
    ) -> Result<Self, CliError> {
        let outputs = pipeline_output_plan(config, sample_id, single_sample)?;
        let input_description = input_description(config);
        let batch_display_mode = !single_sample || mode == PipelineMode::Batch;
        let starting_line = batch_display_mode.then(|| {
            let mut line = String::from("+ ");
            if let Some(display_count) = display_count {
                line.push_str(display_count);
                line.push(' ');
            }
            line.push_str("Starting ");
            line.push_str(&outputs.display_name);
            if config.verbose {
                line.push_str(" (inputs: ");
                line.push_str(&input_description);
                line.push(')');
            }
            line
        });

        let mut intermediate_files_to_delete = Vec::new();
        if config.delete_intermediate_files {
            if let Some(path) = &outputs.matcher_json {
                intermediate_files_to_delete.push(path.clone());
            }
            if let Some(path) = &outputs.phenotyper_json {
                intermediate_files_to_delete.push(path.clone());
            }
        }

        let error_file = outputs
            .base_dir
            .join(format!("{}.ERROR.txt", outputs.basename));
        let save_messages = save_messages(config, &outputs, batch_display_mode);
        let finished_block = batch_display_mode.then(|| {
            let mut block = String::from("- ");
            if let Some(display_count) = display_count {
                block.push_str(display_count);
                block.push(' ');
            }
            block.push_str("Finished processing ");
            block.push_str(&outputs.display_name);
            block.push('\n');
            for message in &save_messages {
                block.push_str("  * ");
                block.push_str(message);
                block.push('\n');
            }
            block
        });

        Ok(Self {
            outputs,
            input_description,
            batch_display_mode,
            starting_line,
            intermediate_files_to_delete,
            error_file,
            save_messages,
            finished_block,
        })
    }
}

/// Runs the first Rust VCF-to-reporter pipeline slice.
pub fn run_reporter_from_vcf(
    vcf_path: &Path,
    sample_name: Option<&str>,
    definitions: &DefinitionReader,
    phenotypes: &PhenotypeMap,
    guidance: &PgkbGuidelineCollection,
    outputs: Option<&PipelineOutputPlan>,
    options: &ReporterPipelineOptions,
) -> Result<ReporterPipelineRun, ReporterPipelineError> {
    let records = read_record_summaries(vcf_path, sample_name)?;
    let sample_name = records.sample_name.clone();
    let mut vcf_warnings = records.warnings.clone();
    let allele_map = allele_map_from_vcf_records(allele_calls_for_locations_with_genes(
        &records,
        definitions.locations_of_interest(),
        definitions.locations_by_gene(),
        options.include_combinations,
        &mut vcf_warnings,
    ));

    let mut gene_call_results = Vec::new();
    let mut report_genes = Vec::new();
    for gene in definitions.genes() {
        let definition = definitions
            .definition_file(gene)
            .expect("definition gene yielded by DefinitionReader::genes");
        let result = match gene {
            "DPYD" => call_dpyd_lowest_function_gene(&sample_name, definition, &allele_map),
            "RYR1" => call_ryr1_lowest_function_gene(&sample_name, definition, &allele_map),
            _ => call_standard_gene_with_exemption(
                &sample_name,
                definition,
                definitions.exemption(gene),
                &allele_map,
                true,
                options.include_combinations,
            ),
        }
        .map_err(|source| ReporterPipelineError::Match {
            gene: gene.to_owned(),
            source,
        })?;
        if let Some(mut report_gene) =
            ReportGene::from_gene_call_result_with_definition_and_messages(
                &result,
                phenotypes.phenotype(gene),
                definition,
                options.message_catalog.as_ref(),
            )?
        {
            report_gene.add_variant_warning_messages(&vcf_warnings);
            report_genes.push(report_gene);
        }
        gene_call_results.push(result);
    }

    let outside_validation = outside_call_validation(definitions, phenotypes);
    for path in &options.outside_call_files {
        for call in parse_outside_calls_file(&outside_validation, path)? {
            let report_gene =
                ReportGene::from_outside_call(&call, phenotypes.phenotype(&call.gene))?;
            report_genes.push(report_gene);
        }
    }

    let mut context = ReportContext::from_gene_reports(
        guidance,
        report_genes,
        outputs.and_then(|outputs| outputs.reporter_title.clone()),
    );
    if let Some(catalog) = options.message_catalog.as_ref() {
        context.apply_report_as_genotype_messages(catalog);
    }
    let written_outputs = if let Some(outputs) = outputs {
        write_reporter_outputs(&context, definitions, outputs, options)?
    } else {
        Vec::new()
    };

    Ok(ReporterPipelineRun {
        context,
        gene_call_results,
        vcf_warnings,
        written_outputs,
    })
}

/// Runs the currently ported CLI path: VCF input through matcher/phenotyper/reporter outputs.
pub fn run_cli_config(
    config: &PharmcatCliConfig,
    options: &CliPipelineOptions,
) -> Result<Vec<ReporterPipelineRun>, CliPipelineError> {
    if !(config.run_matcher && config.run_phenotyper && config.run_reporter) {
        return Err(CliPipelineError::Unsupported(
            "Only the full -vcf matcher + phenotyper + reporter path is currently ported"
                .to_owned(),
        ));
    }
    if config.phenotyper_input.is_some() || config.reporter_input.is_some() {
        return Err(CliPipelineError::Unsupported(
            "Independent phenotyper/reporter input files are not wired to the Rust CLI pipeline yet"
                .to_owned(),
        ));
    }

    let vcf_path = config
        .matcher_vcf
        .as_deref()
        .ok_or_else(|| CliPipelineError::Unsupported("No matcher VCF input file".to_owned()))?;
    let definitions_dir = config
        .definition_dir
        .as_deref()
        .unwrap_or(options.resources.definitions_dir.as_path());
    let definitions = load_cli_definitions(definitions_dir, &config.genes)?;
    let phenotypes = PhenotypeMap::from_dir(&options.resources.phenotype_dir)?;
    let guidance = PgkbGuidelineCollection::from_path(&options.resources.prescribing_guidance)?;
    let message_catalog = MessageCatalog::from_path(&options.resources.reporter_messages)?;
    let reporter_options = ReporterPipelineOptions {
        include_combinations: config.find_combinations,
        outside_call_files: config.phenotyper_outside_call_files.clone(),
        message_catalog: Some(message_catalog),
        html: HtmlReportOptions {
            compact: config.reporter_compact,
            ..HtmlReportOptions::default()
        },
        calls_only_tsv: CallsOnlyTsvOptions {
            show_sample_id: !config.samples.is_empty(),
            ..CallsOnlyTsvOptions::default()
        },
    };

    let sample_names: Vec<Option<&str>> = if config.samples.is_empty() {
        vec![None]
    } else {
        config
            .samples
            .iter()
            .map(|sample| Some(sample.as_str()))
            .collect()
    };
    let single_sample = sample_names.len() <= 1;
    let mut runs = Vec::new();
    for sample_name in sample_names {
        let outputs = pipeline_output_plan(config, sample_name, single_sample)?;
        let run = run_reporter_from_vcf(
            vcf_path,
            sample_name,
            &definitions,
            &phenotypes,
            &guidance,
            Some(&outputs),
            &reporter_options,
        )?;
        write_cli_intermediate_outputs(&run, &outputs, &definitions, vcf_path, config)?;
        runs.push(run);
    }

    Ok(runs)
}

fn write_cli_intermediate_outputs(
    run: &ReporterPipelineRun,
    outputs: &PipelineOutputPlan,
    definitions: &DefinitionReader,
    vcf_path: &Path,
    config: &PharmcatCliConfig,
) -> Result<(), CliPipelineError> {
    if let Some(path) = &outputs.matcher_json {
        let json = matcher_results_json(run, outputs, definitions, vcf_path, config)?;
        fs::write(path, serde_json::to_string_pretty(&json)?).map_err(ReporterPipelineError::Io)?;
    }
    if let Some(path) = &outputs.matcher_warnings {
        fs::write(path, matcher_warnings_text(&run.vcf_warnings))
            .map_err(ReporterPipelineError::Io)?;
    }
    if let Some(path) = &outputs.matcher_html {
        fs::write(path, matcher_html_string(run, definitions, vcf_path))
            .map_err(ReporterPipelineError::Io)?;
    }
    if let Some(path) = &outputs.phenotyper_json {
        write_report_json(&run.context, path)?;
    }
    Ok(())
}

fn matcher_results_json(
    run: &ReporterPipelineRun,
    outputs: &PipelineOutputPlan,
    definitions: &DefinitionReader,
    vcf_path: &Path,
    config: &PharmcatCliConfig,
) -> Result<serde_json::Value, CliPipelineError> {
    let sample_id = matcher_metadata_sample_id(run, outputs);
    let sample_props =
        matcher_sample_props_json(config.sample_metadata_file.as_deref(), sample_id)?;
    Ok(serde_json::json!({
        "metadata": {
            "namedAlleleMatcherVersion": env!("CARGO_PKG_VERSION"),
            "genomeBuild": definitions.genome_build().ok(),
            "inputFilename": vcf_path.file_name().and_then(|name| name.to_str()).unwrap_or_default(),
            "timestamp": current_timestamp_iso8601(),
            "topCandidatesOnly": config.top_candidate_only,
            "findCombinations": config.find_combinations,
            "callCyp2d": config.call_cyp2d6,
            "sampleId": sample_id,
            "sampleProps": sample_props,
        },
        "results": run
            .gene_call_results
            .iter()
            .map(|result| matcher_gene_result_json(result, definitions.definition_file(&result.gene)))
            .collect::<Vec<_>>(),
        "vcfWarnings": run.vcf_warnings,
    }))
}

fn matcher_metadata_sample_id<'a>(
    run: &'a ReporterPipelineRun,
    outputs: &'a PipelineOutputPlan,
) -> &'a str {
    run.gene_call_results
        .first()
        .map(|result| result.match_data.sample_id.as_str())
        .or(run.context.title.as_deref())
        .unwrap_or(&outputs.basename)
}

fn matcher_sample_props_json(
    sample_metadata_file: Option<&Path>,
    sample_id: &str,
) -> Result<serde_json::Value, CliPipelineError> {
    let Some(sample_metadata_file) = sample_metadata_file else {
        return Ok(serde_json::Value::Null);
    };

    let text = fs::read_to_string(sample_metadata_file).map_err(ReporterPipelineError::Io)?;
    let mut sample_props = BTreeMap::new();
    for line in text.lines() {
        let row = line.split('\t').collect::<Vec<_>>();
        if row.len() >= 3 && row[0] == sample_id {
            sample_props.insert(row[1].to_owned(), row[2].to_owned());
        }
    }

    if sample_props.is_empty() {
        Ok(serde_json::Value::Null)
    } else {
        Ok(serde_json::json!(sample_props))
    }
}

fn current_timestamp_iso8601() -> String {
    system_time_to_iso8601(SystemTime::now())
}

fn system_time_to_iso8601(time: SystemTime) -> String {
    let millis = time
        .duration_since(UNIX_EPOCH)
        .map(|duration| duration.as_millis() as i128)
        .unwrap_or(0);
    unix_millis_to_iso8601(millis)
}

fn unix_millis_to_iso8601(millis: i128) -> String {
    let seconds = millis.div_euclid(1000);
    let millis_remainder = millis.rem_euclid(1000);
    let days = seconds.div_euclid(86_400);
    let seconds_of_day = seconds.rem_euclid(86_400);
    let (year, month, day) = civil_from_days(days as i64);
    let hour = seconds_of_day / 3600;
    let minute = (seconds_of_day % 3600) / 60;
    let second = seconds_of_day % 60;

    if millis_remainder == 0 {
        format!("{year:04}-{month:02}-{day:02}T{hour:02}:{minute:02}:{second:02}Z")
    } else {
        format!(
            "{year:04}-{month:02}-{day:02}T{hour:02}:{minute:02}:{second:02}.{millis_remainder:03}Z"
        )
    }
}

fn civil_from_days(days: i64) -> (i64, u32, u32) {
    let z = days + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = mp + if mp < 10 { 3 } else { -9 };
    let year = y + if month <= 2 { 1 } else { 0 };
    (year, month as u32, day as u32)
}

fn matcher_gene_result_json(
    result: &GeneCallResult,
    definition: Option<&DefinitionFile>,
) -> serde_json::Value {
    let (diplotypes, haplotypes, haplotype_matches) = match &result.kind {
        GeneCallKind::NoCall => (Vec::new(), Vec::new(), Vec::new()),
        GeneCallKind::Diplotypes(matches) => (
            matches.iter().map(diplotype_match_json).collect(),
            called_haplotypes_json(matches),
            Vec::new(),
        ),
        GeneCallKind::Haplotypes(matches) => (
            Vec::new(),
            Vec::new(),
            matches.iter().map(haplotype_match_json).collect(),
        ),
    };
    serde_json::json!({
        "source": definition.and_then(|definition| definition.source.as_ref()),
        "version": definition.and_then(|definition| definition.version.clone()),
        "chromosome": definition.map(|definition| definition.chromosome.clone()),
        "gene": result.gene,
        "diplotypes": diplotypes,
        "haplotypes": haplotypes,
        "haplotypeMatches": haplotype_matches,
        "phased": result.match_data.phased,
        "variants": variants_json(&result.match_data),
        "variantsOfInterest": result
            .match_data
            .positions_with_undocumented_variations
            .iter()
            .map(extra_variant_json)
            .collect::<Vec<_>>(),
        "matchData": match_data_json(&result.match_data),
        "uncallableHaplotypes": uncallable_haplotypes_json(definition, &result.match_data),
        "warnings": matcher_warning_messages_json(result),
    })
}

fn uncallable_haplotypes_json(
    definition: Option<&DefinitionFile>,
    match_data: &MatchData,
) -> BTreeSet<String> {
    let matchable_haplotypes = match_data
        .haplotypes()
        .iter()
        .map(|haplotype| haplotype.name.as_str())
        .collect::<BTreeSet<_>>();

    definition
        .into_iter()
        .flat_map(|definition| &definition.named_alleles)
        .map(|haplotype| haplotype.name.clone())
        .filter(|name| !matchable_haplotypes.contains(name.as_str()))
        .collect()
}

fn matcher_warning_messages_json(result: &GeneCallResult) -> Vec<serde_json::Value> {
    result
        .warnings
        .iter()
        .map(|warning| gene_call_warning_json(&result.gene, warning))
        .chain(
            result
                .dpyd_hap_b3_warnings
                .iter()
                .map(dpyd_hap_b3_warning_json),
        )
        .collect()
}

fn gene_call_warning_json(gene: &str, warning: &GeneCallWarning) -> serde_json::Value {
    match warning {
        GeneCallWarning::UnphasedPriority => message_annotation_json(
            "unphased-priority",
            None,
            serde_json::Value::Null,
            "note",
            gene_call_warning_message(gene, warning),
        ),
        GeneCallWarning::MissingRequiredPosition(_) => message_annotation_json(
            "missing-required-position",
            None,
            serde_json::Value::Null,
            "note",
            gene_call_warning_message(gene, warning),
        ),
        GeneCallWarning::MissingAmp1Position(_) => message_annotation_json(
            "missing-amp1-position",
            None,
            serde_json::Value::Null,
            "note",
            gene_call_warning_message(gene, warning),
        ),
    }
}

fn gene_call_warning_message(gene: &str, warning: &GeneCallWarning) -> String {
    match warning {
        GeneCallWarning::UnphasedPriority => format!(
            "Unphased {gene} variants resulted in multiple calls.  PharmCAT is picking a single call based on frequency data.  Please consult the documentation for details."
        ),
        GeneCallWarning::MissingRequiredPosition(positions) => {
            let suffix = if positions.len() > 1 { "s" } else { "" };
            format!(
                "Cannot call {gene} - missing required variant{suffix} ({})",
                positions.join(", ")
            )
        }
        GeneCallWarning::MissingAmp1Position(positions) => {
            let suffix = if positions.len() > 1 { "s" } else { "" };
            format!(
                "Missing variant{suffix} required to meet AMP Tier 1 requirements:  {}. See https://www.clinpgx.org/ampAllelesToTest for details.",
                positions.join(", ")
            )
        }
    }
}

fn dpyd_hap_b3_warning_json(warning: &DpydHapB3Warning) -> serde_json::Value {
    match warning {
        DpydHapB3Warning::IntronicMismatchExonic => message_annotation_json(
            "pcat-dpyd-hapB3-intronic-mismatch-exonic",
            Some("2"),
            empty_match_logic_json(),
            "",
            "The HapB3 haplotype is comprised of an exonic SNP at rs56038477 and an intronic SNP (c.1129-5923C>G) at rs75017182.  HapB3’s decreased function is thought to be due to the intronic SNP. This genotype was assigned based on the c.1129-5923C>G (intron) SNP alone because the VCF input for the exon SNP is either missing or conflicts with the intron SNP.",
        ),
        DpydHapB3Warning::ExonicOnly => message_annotation_json(
            "pcat-dpyd-hapB3-exonic-only",
            Some("2"),
            empty_match_logic_json(),
            "",
            "The HapB3 haplotype is comprised of an exonic SNP at rs56038477 and an intronic SNP (c.1129-5923C>G) at rs75017182.  HapB3’s decreased function is thought to be due to the intronic SNP. This genotype was assigned based on the rs56038477 (exon) SNP alone since the rs75017182 (intron) SNP is missing from the VCF input file. The rs75017182 (intron) SNP should be genotyped to confirm the presence or absence of HapB3 and decreased function.",
        ),
        DpydHapB3Warning::AlleleCount { count, rsid } => message_annotation_json(
            "warn.alleleCount",
            None,
            serde_json::Value::Null,
            "note",
            format!(
                "Only found {count} allele for {}",
                rsid.as_deref().unwrap_or_default()
            ),
        ),
    }
}

fn message_annotation_json(
    rule_name: &str,
    version: Option<&str>,
    matches: serde_json::Value,
    exception_type: &str,
    message: impl Into<String>,
) -> serde_json::Value {
    serde_json::json!({
        "rule_name": rule_name,
        "version": version,
        "matches": matches,
        "exception_type": exception_type,
        "message": message.into(),
    })
}

fn empty_match_logic_json() -> serde_json::Value {
    serde_json::json!({
        "gene": serde_json::Value::Null,
        "hapsCalled": Vec::<String>::new(),
        "hapsMissing": Vec::<String>::new(),
        "variantsMissing": Vec::<String>::new(),
        "variant": Vec::<String>::new(),
        "dips": Vec::<String>::new(),
        "drugs": Vec::<String>::new(),
    })
}

fn called_haplotypes_json(matches: &[DiplotypeMatch]) -> Vec<serde_json::Value> {
    let mut haplotypes = Vec::new();
    let mut seen = BTreeSet::new();
    for match_ in matches {
        push_called_haplotype_json(&mut haplotypes, &mut seen, &match_.haplotype1);
        if let Some(haplotype2) = &match_.haplotype2 {
            push_called_haplotype_json(&mut haplotypes, &mut seen, haplotype2);
        }
    }
    haplotypes
}

fn push_called_haplotype_json(
    haplotypes: &mut Vec<serde_json::Value>,
    seen: &mut BTreeSet<(String, Vec<String>)>,
    haplotype: &HaplotypeMatch,
) {
    let sequences = haplotype.sequences.iter().cloned().collect::<Vec<_>>();
    if seen.insert((haplotype.name.clone(), sequences)) {
        haplotypes.push(haplotype_match_json(haplotype));
    }
}

fn match_data_json(match_data: &MatchData) -> serde_json::Value {
    serde_json::json!({
        "missingPositions": match_data
            .missing_positions
            .iter()
            .map(variant_locus_json)
            .collect::<Vec<_>>(),
        "positionsWithUndocumentedVariations": match_data
            .positions_with_undocumented_variations
            .iter()
            .map(variant_locus_json)
            .collect::<Vec<_>>(),
        "treatUndocumentedVariationsAsReference": match_data.treat_undocumented_variations_as_reference,
        "phased": match_data.phased,
        "phaseSets": phase_sets_json(match_data),
        "posToPhaseSet": position_to_phase_set_json(match_data),
        "effectivelyPhased": match_data.effectively_phased,
        "homozygous": match_data.homozygous,
        "missingRequiredPositions": match_data.missing_required_positions,
        "missingAmp1Positions": match_data.missing_amp1_positions,
    })
}

fn phase_sets_json(match_data: &MatchData) -> BTreeMap<i32, Vec<u64>> {
    let mut phase_sets = BTreeMap::new();
    for position in &match_data.positions {
        let Some(sample) = match_data.sample_allele_at_position(position.position) else {
            continue;
        };
        if sample.phased() {
            let phase_set = sample.phase_set().unwrap_or(i32::MIN);
            phase_sets
                .entry(phase_set)
                .or_insert_with(Vec::new)
                .push(position.position);
        }
    }
    phase_sets
}

fn position_to_phase_set_json(match_data: &MatchData) -> BTreeMap<u64, i32> {
    let phase_sets = phase_sets_json(match_data);
    let uses_phase_sets = if phase_sets.contains_key(&i32::MIN) {
        phase_sets.len() > 1
    } else {
        !phase_sets.is_empty()
    };
    if !uses_phase_sets {
        return BTreeMap::new();
    }

    let mut position_to_phase_set = BTreeMap::new();
    for (phase_set, positions) in phase_sets {
        for position in positions {
            position_to_phase_set.insert(position, phase_set);
        }
    }
    position_to_phase_set
}

fn variants_json(match_data: &MatchData) -> Vec<serde_json::Value> {
    match_data
        .positions
        .iter()
        .map(|position| {
            let sample = match_data.sample_allele_at_position(position.position);
            serde_json::json!({
                "position": position.position,
                "rsid": position.rsid,
                "vcfCall": sample.map(|sample| sample.vcf_call()),
                "phased": sample.is_some_and(|sample| sample.phased()),
                "phaseSet": sample.and_then(|sample| sample.phase_set()),
            })
        })
        .collect()
}

fn extra_variant_json(locus: &VariantLocus) -> serde_json::Value {
    serde_json::json!({
        "position": locus.position,
        "rsid": locus.rsid,
        "vcfCall": serde_json::Value::Null,
        "phased": false,
        "phaseSet": serde_json::Value::Null,
    })
}

fn diplotype_match_json(match_: &DiplotypeMatch) -> serde_json::Value {
    serde_json::json!({
        "name": match_.name,
        "score": match_.score,
        "haplotype1": haplotype_match_json(&match_.haplotype1),
        "haplotype2": match_.haplotype2.as_ref().map(haplotype_match_json),
        "sequencePairs": match_.sequence_pairs,
    })
}

fn haplotype_match_json(match_: &HaplotypeMatch) -> serde_json::Value {
    serde_json::json!({
        "name": match_.name,
        "haplotype": named_allele_json(&match_.haplotype),
        "positions": match_
            .positions
            .iter()
            .map(variant_locus_json)
            .collect::<Vec<_>>(),
        "sequences": match_.sequences,
    })
}

fn named_allele_json(allele: &NamedAllele) -> serde_json::Value {
    serde_json::json!({
        "name": allele.name,
        "id": allele.id,
        "alleles": allele.alleles,
        "cpicAlleles": allele.cpic_alleles,
        "reference": allele.reference,
        "structuralVariant": allele.structural_variant,
        "corePositions": allele.core_positions,
        "score": named_allele_json_score(allele),
        "numCombinations": allele.num_combinations,
        "numPartials": allele.num_partials,
    })
}

fn named_allele_json_score(allele: &NamedAllele) -> i32 {
    allele.score_override.unwrap_or_else(|| {
        allele
            .alleles
            .iter()
            .filter(|allele| allele.is_some())
            .count() as i32
            - allele.num_partials
    })
}

fn variant_locus_json(locus: &VariantLocus) -> serde_json::Value {
    serde_json::json!({
        "chromosome": locus.chromosome,
        "position": locus.position,
        "cpicPosition": locus.cpic_position,
        "rsid": locus.rsid,
        "chromosomeHgvsName": locus.chromosome_hgvs_name,
        "cpicAlleles": locus.cpic_alleles,
        "cpicToVcfAlleleMap": locus.cpic_to_vcf_allele_map,
        "ref": locus.reference,
        "alts": locus.alts,
    })
}

fn matcher_warnings_text(warnings: &VcfWarnings) -> String {
    let mut text = String::new();
    for (position, messages) in warnings {
        for message in messages {
            text.push_str(position);
            text.push('\t');
            text.push_str(message);
            text.push('\n');
        }
    }
    text
}

fn matcher_html_string(
    run: &ReporterPipelineRun,
    definitions: &DefinitionReader,
    vcf_path: &Path,
) -> String {
    matcher_html_string_with_options(run, definitions, vcf_path, MatcherHtmlOptions::default())
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct MatcherHtmlOptions {
    always_show_unmatched_haplotypes: bool,
}

fn matcher_html_string_with_options(
    run: &ReporterPipelineRun,
    definitions: &DefinitionReader,
    vcf_path: &Path,
    options: MatcherHtmlOptions,
) -> String {
    let sample_id = run
        .gene_call_results
        .first()
        .map(|result| result.match_data.sample_id.as_str())
        .unwrap_or("unknown");
    let input_filename = vcf_path
        .file_name()
        .and_then(|name| name.to_str())
        .map(str::to_owned)
        .unwrap_or_else(|| vcf_path.to_string_lossy().into_owned());
    let title = format!("PharmCAT Allele Call Report for {sample_id} in {input_filename}");
    let mut content = String::new();
    for result in &run.gene_call_results {
        matcher_html_gene_section(
            &mut content,
            result,
            definitions.definition_file(&result.gene),
            options,
        );
    }

    render_matcher_html_template(&title, &content, &current_matcher_html_date())
}

const MATCHER_HTML_TEMPLATE: &str = r##"<!DOCTYPE html>
<html class="no-js" lang="en">
<head>
  <meta charset="utf-8" />
  <meta http-equiv="x-ua-compatible" content="ie=edge" />
  <title>${title}</title>
  <meta name="viewport" content="width=device-width, initial-scale=1" />

  <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css" rel="stylesheet"
      integrity="sha384-9ndCyUaIbzAi2FUVXJi0CjmCapSmO7SnpJef0486qhLnuZ2cdeRhO02iuK6FUUVM" crossorigin="anonymous" />
  <link href=" https://cdn.jsdelivr.net/npm/bootswatch@5.3.0/dist/pulse/bootstrap.min.css" rel="stylesheet"
      crossorigin="anonymous"/>
  <style>
    /* Move down content because we have a fixed navbar */
    body {
      padding-top: 3rem;
      padding-bottom: 20px;
      --bs-body-font-size: .8rem;
    }
    .browserUpgrade {
      background: #ccc;
      color: #000;
      padding: 1.5em 1em 1em;
    }
    .navbar {
      --bs-navbar-padding-y: .5rem;
    }
    h3 {
      margin-top: 1em;
    }
    table {
      width: auto !important;
    }
    td, th {
      text-align: center;
      padding-left: 1em !important;
      padding-right: 1em !important;
    }
    th.first {
      text-align: initial;
      margin-right: 1em;
      min-width: 12em;
    }
    .footer {
      margin-top: 4em;
    }
  </style>
</head>
<body>
<!--[if lt IE 9]>
<p class="browserUpgrade">You are using an <strong>outdated</strong> browser. Please <a href="https://browsehappy.com/">upgrade your browser</a> to improve your experience.</p>
<![endif]-->
<nav class="navbar fixed-top navbar-dark bg-dark">
  <div class="container-fluid">
    <a class="navbar-brand" href="#">${title}</a>
  </div>
</nav>

<div class="container-fluid">
${content}
</div>
<div class="footer">
  <hr />
  <footer class="container-fluid">
    <p>
      <small>Generated on ${timestamp}.</small>
    </p>
  </footer>
</div>
</body>
</html>
"##;

fn render_matcher_html_template(title: &str, content: &str, timestamp: &str) -> String {
    let html = MATCHER_HTML_TEMPLATE
        .replace("${title}", title)
        .replace("${content}", content)
        .replace("${timestamp}", timestamp);
    format!("{html}\n")
}

fn matcher_html_gene_section(
    html: &mut String,
    result: &GeneCallResult,
    definition: Option<&DefinitionFile>,
    options: MatcherHtmlOptions,
) {
    html.push_str("<h3>");
    html.push_str(&html_escape_text(&result.gene));
    html.push_str("</h3>\n");

    let mut missing_required = false;
    for warning in &result.warnings {
        match warning {
            GeneCallWarning::UnphasedPriority => {
                html.push_str("<p>");
                html.push_str(&html_escape_text(&gene_call_warning_message(
                    &result.gene,
                    warning,
                )));
                html.push_str("</p>\n");
            }
            GeneCallWarning::MissingRequiredPosition(_) => {
                missing_required = true;
                html.push_str("<p>");
                html.push_str(&html_escape_text(&gene_call_warning_message(
                    &result.gene,
                    warning,
                )));
                html.push_str("</p>\n");
            }
            GeneCallWarning::MissingAmp1Position(_) => {}
        }
    }

    if !missing_required {
        matcher_html_diplotype_list(html, result);
        matcher_html_variant_table(html, result, definition, options);
    }

    if !result.match_data.missing_positions.is_empty() {
        html.push_str("<p>There ");
        if result.match_data.missing_positions.len() > 1 {
            html.push_str("were ");
        } else {
            html.push_str("was ");
        }
        html.push_str(&result.match_data.missing_positions.len().to_string());
        html.push_str(" missing positions from the VCF file:</p>\n<ul>");
        for variant in &result.match_data.missing_positions {
            html.push_str("  <li>");
            html.push_str(&variant.position.to_string());
            html.push_str(" (");
            html.push_str(&html_escape_text(&variant.chromosome_hgvs_name));
            html.push_str(")</li>");
        }
        html.push_str("</ul>\n");

        if !missing_required && let Some(definition) = definition {
            let uncallable = uncallable_haplotypes_json(Some(definition), &result.match_data);
            if !uncallable.is_empty() {
                html.push_str(
                    "<p>The following haplotype(s) were eliminated from consideration:</p><ul>",
                );
                for name in uncallable {
                    html.push_str("  <li>");
                    html.push_str(&html_escape_text(&name));
                    html.push_str("</li>");
                }
                html.push_str("</ul>\n");
            }
        }

        matcher_html_missing_tag_position_notes(html, result);
    }
    html.push('\n');
}

fn matcher_html_diplotype_list(html: &mut String, result: &GeneCallResult) {
    html.push_str("<ul>");
    if let GeneCallKind::Diplotypes(matches) = &result.kind {
        for diplotype in matches {
            html.push_str("  <li>");
            html.push_str(&html_escape_text(&diplotype.name));
            html.push_str(" (");
            html.push_str(&diplotype.score.to_string());
            html.push_str(")</li>");
        }
    }
    html.push_str("</ul>\n");
}

fn matcher_html_variant_table(
    html: &mut String,
    result: &GeneCallResult,
    definition: Option<&DefinitionFile>,
    options: MatcherHtmlOptions,
) {
    let positions = &result.match_data.positions;
    html.push_str("<table class=\"table table-striped table-hover table-sm\">\n");
    matcher_html_header_row(html, "Definition Position", positions, |variant| {
        variant.position.to_string()
    });
    matcher_html_header_row(html, "", positions, |variant| {
        variant.rsid.clone().unwrap_or_default()
    });
    matcher_html_header_row(html, "VCF Position", positions, |variant| {
        variant.position.to_string()
    });
    if result.match_data.is_using_phase_sets() {
        html.push_str("  <tr><th class=\"first\">VCF Phase Set</th>");
        for variant in positions {
            html.push_str("<th>");
            if let Some(phase_set) = result.match_data.phase_set(variant.position) {
                html.push_str(&phase_set.to_string());
            }
            html.push_str("</th>");
        }
        html.push_str("</tr>\n");
    }
    html.push_str("  <tr><th class=\"first\">VCF REF,ALTs</th>");
    for variant in positions {
        html.push_str("<th>");
        if let Some(sample) = result
            .match_data
            .sample_allele_at_position(variant.position)
        {
            html.push_str(&html_escape_text(&sample.vcf_alleles().join(",")));
        }
        html.push_str("</th>");
    }
    html.push_str("</tr>\n");
    html.push_str("  <tr class=\"table-success\"><th class=\"first\">VCF Call</th>");
    let highlight_positions = matcher_html_highlight_positions(result);
    for variant in positions {
        let Some(sample) = result
            .match_data
            .sample_allele_at_position(variant.position)
        else {
            html.push_str("<th />");
            continue;
        };
        if highlight_positions.contains(&variant.position) {
            html.push_str("<th class=\"table-danger\">");
        } else {
            html.push_str("<th>");
        }
        html.push_str(&html_escape_text(sample.vcf_call()));
        html.push_str("</th>");
    }
    html.push_str("</tr>\n");

    match &result.kind {
        GeneCallKind::Diplotypes(matches) => {
            let mut matched_names = BTreeSet::new();
            for diplotype in matches {
                matcher_html_haplotype_row(
                    html,
                    Some(&diplotype.haplotype1.name),
                    &named_allele_position_map(&diplotype.haplotype1.haplotype, positions),
                    Some("table-info"),
                    &highlight_positions,
                    positions,
                );
                matched_names.insert(diplotype.haplotype1.name.clone());
                if let Some(haplotype2) = &diplotype.haplotype2 {
                    matcher_html_haplotype_row(
                        html,
                        Some(&haplotype2.name),
                        &named_allele_position_map(&haplotype2.haplotype, positions),
                        Some("table-info"),
                        &highlight_positions,
                        positions,
                    );
                    matched_names.insert(haplotype2.name.clone());
                }
                for sequence_pair in &diplotype.sequence_pairs {
                    for sequence in sequence_pair {
                        matcher_html_haplotype_row(
                            html,
                            None,
                            &sequence_position_map(sequence),
                            None,
                            &highlight_positions,
                            positions,
                        );
                    }
                }
            }
            matcher_html_unmatched_haplotype_rows(
                html,
                result,
                definition,
                &matched_names,
                &highlight_positions,
                options,
            );
        }
        GeneCallKind::Haplotypes(matches) => {
            let mut matched_names = BTreeSet::new();
            for haplotype in matches {
                matcher_html_haplotype_row(
                    html,
                    Some(&haplotype.name),
                    &named_allele_position_map(&haplotype.haplotype, positions),
                    Some("table-info"),
                    &highlight_positions,
                    positions,
                );
                matched_names.insert(haplotype.name.clone());
                for sequence in &haplotype.sequences {
                    matcher_html_haplotype_row(
                        html,
                        None,
                        &sequence_position_map(sequence),
                        None,
                        &highlight_positions,
                        positions,
                    );
                }
            }
            matcher_html_unmatched_haplotype_rows(
                html,
                result,
                definition,
                &matched_names,
                &highlight_positions,
                options,
            );
        }
        GeneCallKind::NoCall => matcher_html_unmatched_haplotype_rows(
            html,
            result,
            definition,
            &BTreeSet::new(),
            &highlight_positions,
            options,
        ),
    }
    html.push_str("</table>\n");
}

fn matcher_html_unmatched_haplotype_rows(
    html: &mut String,
    result: &GeneCallResult,
    definition: Option<&DefinitionFile>,
    matched_names: &BTreeSet<String>,
    highlight_positions: &BTreeSet<u64>,
    options: MatcherHtmlOptions,
) {
    if !options.always_show_unmatched_haplotypes && !matched_names.is_empty() {
        return;
    }
    let Some(definition) = definition else {
        return;
    };
    for haplotype in result.match_data.haplotypes() {
        if matched_names.contains(&haplotype.name) {
            continue;
        }
        matcher_html_haplotype_row(
            html,
            Some(&haplotype.name),
            &named_allele_definition_position_map(
                haplotype,
                definition,
                &result.match_data.positions,
            ),
            Some("table-danger"),
            highlight_positions,
            &result.match_data.positions,
        );
    }
}

fn matcher_html_header_row<F>(html: &mut String, label: &str, positions: &[VariantLocus], value: F)
where
    F: Fn(&VariantLocus) -> String,
{
    html.push_str("  <tr><th");
    if !label.is_empty() {
        html.push_str(" class=\"first\"");
    }
    html.push('>');
    html.push_str(&html_escape_text(label));
    html.push_str("</th>");
    for variant in positions {
        html.push_str("<th>");
        html.push_str(&html_escape_text(&value(variant)));
        html.push_str("</th>");
    }
    html.push_str("</tr>");
}

fn matcher_html_haplotype_row(
    html: &mut String,
    name: Option<&str>,
    alleles_by_position: &BTreeMap<u64, String>,
    row_class: Option<&str>,
    highlight_positions: &BTreeSet<u64>,
    _positions: &[VariantLocus],
) {
    html.push_str("  <tr");
    if let Some(row_class) = row_class {
        html.push_str(" class=\"");
        html.push_str(row_class);
        html.push('"');
    }
    html.push_str("><th class=\"first\">");
    if let Some(name) = name {
        html.push_str(name);
    }
    html.push_str("</th>");
    for (position, call) in alleles_by_position {
        let call = call.replace('\\', "");
        let call = call.as_str();
        let is_any = call == ".*?";
        html.push_str("<td");
        if highlight_positions.contains(position) && !is_any {
            html.push_str(" class=\"table-danger\"");
        }
        html.push('>');
        if name.is_some() && !is_any {
            html.push_str("<b>");
            html.push_str(call);
            html.push_str("</b>");
        } else {
            html.push_str(call);
        }
        html.push_str("</td>");
    }
    html.push_str("</tr>\n");
}

fn matcher_html_highlight_positions(result: &GeneCallResult) -> BTreeSet<u64> {
    result
        .match_data
        .positions
        .iter()
        .filter_map(|variant| {
            let sample = result
                .match_data
                .sample_allele_at_position(variant.position)?;
            let alleles = sample
                .vcf_call()
                .split(if sample.phased() { '|' } else { '/' })
                .collect::<BTreeSet<_>>();
            let first_allele = alleles.iter().next().copied();
            if alleles.len() > 1
                || first_allele.is_some_and(|allele| allele != variant.reference.as_str())
            {
                Some(variant.position)
            } else {
                None
            }
        })
        .collect()
}

fn named_allele_position_map(
    allele: &NamedAllele,
    positions: &[VariantLocus],
) -> BTreeMap<u64, String> {
    positions
        .iter()
        .enumerate()
        .map(|(index, variant)| {
            (
                variant.position,
                named_allele_pattern_cell(allele.alleles.get(index).and_then(Option::as_deref)),
            )
        })
        .collect()
}

fn named_allele_definition_position_map(
    allele: &NamedAllele,
    definition: &DefinitionFile,
    positions: &[VariantLocus],
) -> BTreeMap<u64, String> {
    positions
        .iter()
        .filter_map(|variant| {
            definition
                .index_for_position(variant.position)
                .map(|index| {
                    (
                        variant.position,
                        named_allele_pattern_cell(
                            allele.alleles.get(index).and_then(Option::as_deref),
                        ),
                    )
                })
        })
        .collect()
}

fn named_allele_pattern_cell(allele: Option<&str>) -> String {
    let Some(allele) = allele else {
        return ".*?".to_owned();
    };
    if allele.chars().count() == 1 {
        return matcher_html_iupac_regex(allele).to_owned();
    }
    if allele.contains(" or ") {
        return format!("({})", allele.split(" or ").collect::<Vec<_>>().join("|"));
    }
    allele.to_owned()
}

fn matcher_html_iupac_regex(allele: &str) -> &str {
    match allele.to_ascii_uppercase().as_str() {
        "A" => "A",
        "C" => "C",
        "G" => "G",
        "T" => "T",
        "R" => "[AG]",
        "Y" => "[CT]",
        "S" => "[GC]",
        "W" => "[AT]",
        "K" => "[GT]",
        "M" => "[AC]",
        "B" => "[CGT]",
        "D" => "[AGT]",
        "H" => "[ACT]",
        "V" => "[ACG]",
        "N" => "[ACGT]",
        "-" => "del",
        _ => allele,
    }
}

fn matcher_html_missing_tag_position_notes(html: &mut String, result: &GeneCallResult) {
    let called_haplotypes = match &result.kind {
        GeneCallKind::Diplotypes(matches) => matches
            .iter()
            .flat_map(|diplotype| {
                std::iter::once(&diplotype.haplotype1).chain(diplotype.haplotype2.iter())
            })
            .collect::<Vec<_>>(),
        GeneCallKind::Haplotypes(matches) => matches.iter().collect::<Vec<_>>(),
        GeneCallKind::NoCall => Vec::new(),
    };
    if called_haplotypes.is_empty()
        || !called_haplotypes
            .iter()
            .any(|match_| !match_.haplotype.missing_positions.is_empty())
    {
        return;
    }

    html.push_str(
        "<p>The following haplotypes were called even though tag positions were missing:</p>\n<ul>",
    );
    for match_ in called_haplotypes {
        if match_.haplotype.missing_positions.is_empty() {
            continue;
        }
        let positions = match_
            .haplotype
            .missing_positions
            .iter()
            .map(|variant| variant.chromosome_hgvs_name.as_str())
            .collect::<Vec<_>>()
            .join(", ");
        html.push_str("  <li>Called ");
        html.push_str(&html_escape_text(&match_.name));
        html.push_str(" without ");
        html.push_str(&html_escape_text(&positions));
        html.push_str("</li>");
    }
    html.push_str("</ul>\n");
}

fn sequence_position_map(sequence: &str) -> BTreeMap<u64, String> {
    sequence
        .split(';')
        .filter_map(|pos_allele| {
            let (position, allele) = pos_allele.split_once(':')?;
            let position = position.parse::<u64>().ok()?;
            let allele = match allele {
                ".?" => "",
                other => other,
            };
            Some((position, allele.replace('\\', "")))
        })
        .collect()
}

fn current_matcher_html_date() -> String {
    let millis = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|duration| duration.as_millis() as i128)
        .unwrap_or(0);
    matcher_html_date_from_millis(millis)
}

fn matcher_html_date_from_millis(millis: i128) -> String {
    let days = millis.div_euclid(1000).div_euclid(86_400);
    let (year, month, day) = civil_from_days(days as i64);
    format!("{month:02}/{day:02}/{:02}", year.rem_euclid(100))
}

fn html_escape_text(value: &str) -> String {
    value
        .replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

fn load_cli_definitions(
    definitions_dir: &Path,
    genes: &std::collections::BTreeSet<String>,
) -> Result<DefinitionReader, DefinitionLoadError> {
    let mut paths = fs::read_dir(definitions_dir)?
        .map(|entry| entry.map(|entry| entry.path()))
        .collect::<Result<Vec<_>, _>>()?;
    paths.retain(|path| {
        let is_json = path
            .extension()
            .is_some_and(|extension| extension == "json");
        let Some(file_name) = path.file_name().and_then(|name| name.to_str()) else {
            return false;
        };
        is_json
            && file_name.ends_with("_translation.json")
            && (genes.is_empty()
                || genes.iter().any(|gene| {
                    file_name
                        .strip_suffix("_translation.json")
                        .is_some_and(|file_gene| file_gene.eq_ignore_ascii_case(gene))
                }))
    });
    paths.sort();

    let exemptions_path = definitions_dir.join("exemptions.json");
    if exemptions_path.is_file() {
        DefinitionReader::from_paths_with_exemptions(paths, exemptions_path)
    } else {
        DefinitionReader::from_paths(paths)
    }
}

fn allele_map_from_vcf_records(
    records: impl IntoIterator<Item = SampleAlleleSummary>,
) -> BTreeMap<String, SampleAlleleSummary> {
    records
        .into_iter()
        .map(|record| (format!("{}:{}", record.chromosome, record.position), record))
        .collect()
}

fn outside_call_validation(
    definitions: &DefinitionReader,
    phenotypes: &PhenotypeMap,
) -> OutsideCallValidation {
    let mut validation =
        OutsideCallValidation::for_supported_genes(definitions.genes().map(str::to_owned));
    validation
        .supported_genes
        .extend(phenotypes.genes().map(str::to_owned));

    for gene in phenotypes.genes() {
        if phenotypes
            .phenotype(gene)
            .is_some_and(|phenotype| phenotype.is_activity_gene())
        {
            validation.activity_score_genes.insert(gene.to_owned());
        }
    }

    for gene in definitions.genes() {
        let definition = definitions
            .definition_file(gene)
            .expect("definition gene yielded by DefinitionReader::genes");
        validation.valid_named_alleles.insert(
            gene.to_owned(),
            definition
                .named_alleles
                .iter()
                .map(|allele| allele.name.clone())
                .collect(),
        );
    }

    validation
}

fn write_reporter_outputs(
    context: &ReportContext,
    definitions: &DefinitionReader,
    outputs: &PipelineOutputPlan,
    options: &ReporterPipelineOptions,
) -> Result<Vec<PathBuf>, ReporterPipelineError> {
    fs::create_dir_all(&outputs.base_dir)?;
    let mut written = Vec::new();
    if let Some(path) = &outputs.reporter_json {
        write_report_json(context, path)?;
        written.push(path.clone());
    }
    if let Some(path) = &outputs.reporter_html {
        let mut html_options = options.html.clone();
        if html_options.definition_genes.is_empty() {
            html_options.definition_genes = definitions.genes().map(str::to_owned).collect();
        }
        write_report_html(context, path, &html_options)?;
        written.push(path.clone());
    }
    if let Some(path) = &outputs.reporter_calls_only_tsv {
        let mut tsv_options = options.calls_only_tsv.clone();
        if tsv_options.sample_id.is_none() {
            tsv_options.sample_id = context.title.clone();
        }
        write_calls_only_tsv(context, path, &tsv_options)?;
        written.push(path.clone());
    }
    Ok(written)
}

/// Java `Pipeline.getInputDescription`.
pub fn input_description(config: &PharmcatCliConfig) -> String {
    let mut inputs = Vec::new();
    if let Some(path) = &config.matcher_vcf {
        inputs.push(file_name(path));
    }
    if let Some(path) = &config.phenotyper_input {
        inputs.push(file_name(path));
    }
    for path in &config.phenotyper_outside_call_files {
        inputs.push(file_name(path));
    }
    if let Some(path) = &config.reporter_input {
        inputs.push(file_name(path));
    }
    inputs.join(", ")
}

/// Java `Pipeline.call` save-output messages for normal stage outputs.
pub fn save_messages(
    config: &PharmcatCliConfig,
    outputs: &PipelineOutputPlan,
    batch_display_mode: bool,
) -> Vec<String> {
    if batch_display_mode {
        return Vec::new();
    }

    let mut messages = Vec::new();
    if config.run_matcher && (!config.delete_intermediate_files || !config.run_phenotyper) {
        if let Some(path) = &outputs.matcher_json {
            messages.push(format!(
                "Saving named allele matcher JSON results to {}",
                path.display()
            ));
        }
        if let Some(path) = &outputs.matcher_html {
            messages.push(format!(
                "Saving named allele matcher HTML results to {}",
                path.display()
            ));
        }
    }

    if config.run_phenotyper
        && (!config.delete_intermediate_files || !config.run_reporter)
        && let Some(path) = &outputs.phenotyper_json
    {
        messages.push(format!(
            "Saving phenotyper JSON results to {}",
            path.display()
        ));
    }

    if config.run_reporter {
        if let Some(path) = &outputs.reporter_html {
            messages.push(format!(
                "Saving reporter HTML results to {}",
                path.display()
            ));
        }
        if let Some(path) = &outputs.reporter_json {
            messages.push(format!(
                "Saving reporter JSON results to {}",
                path.display()
            ));
        }
        if let Some(path) = &outputs.reporter_calls_only_tsv {
            messages.push(format!(
                "Saving calls-only TSV results to {}",
                path.display()
            ));
        }
    }

    messages
}

/// Deletes Java pipeline intermediate files with `Files.deleteIfExists` semantics.
pub fn delete_intermediate_files(plan: &PipelineRunPlan) -> io::Result<Vec<PathBuf>> {
    let mut deleted = Vec::new();
    for path in &plan.intermediate_files_to_delete {
        match fs::remove_file(path) {
            Ok(()) => deleted.push(path.clone()),
            Err(error) if error.kind() == io::ErrorKind::NotFound => {}
            Err(error) => return Err(error),
        }
    }
    Ok(deleted)
}

/// Writes a Java-style pipeline error file.
pub fn write_error_file(path: &Path, error: &(dyn std::error::Error + 'static)) -> io::Result<()> {
    fs::write(path, format!("{error:?}\n"))
}

fn file_name(path: &Path) -> String {
    path.file_name()
        .and_then(|name| name.to_str())
        .unwrap_or_default()
        .to_string()
}

#[cfg(test)]
mod tests {
    use std::{
        collections::{BTreeMap, BTreeSet},
        fs,
        path::PathBuf,
        time::{SystemTime, UNIX_EPOCH},
    };

    use crate::cli::{CliAction, parse_pharmcat_args};
    use crate::definition::{DefinitionFile, DefinitionReader, VariantLocus, read_definition_file};
    use crate::phenotype::PhenotypeMap;
    use crate::report::{MessageAnnotation, PgkbGuidelineCollection, PrescribingGuidanceSource};
    use crate::vcf::ReadVcfError;

    use super::*;

    const GUIDANCE_PATH: &str = "../../repos/PharmCAT/src/main/resources/org/pharmgkb/pharmcat/reporter/prescribing_guidance.json";
    const MESSAGE_PATH: &str =
        "../../repos/PharmCAT/src/main/resources/org/pharmgkb/pharmcat/reporter/messages.json";
    const PHENOTYPE_PATH: &str =
        "../../repos/PharmCAT/src/main/resources/org/pharmgkb/pharmcat/phenotype";
    const CYP3A5_DEFINITION_PATH: &str =
        "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/haplotyper.json";
    const CYP3A5_VCF_PATH: &str =
        "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/haplotyper.vcf";
    const ABCG2_DEFINITION_PATH: &str = "../../repos/PharmCAT/src/main/resources/org/pharmgkb/pharmcat/definition/alleles/ABCG2_translation.json";
    const CYP2C19_DEFINITION_PATH: &str = "../../repos/PharmCAT/src/main/resources/org/pharmgkb/pharmcat/definition/alleles/CYP2C19_translation.json";
    const DPYD_DEFINITION_PATH: &str = "../../repos/PharmCAT/src/main/resources/org/pharmgkb/pharmcat/definition/alleles/DPYD_translation.json";
    const RYR1_DEFINITION_PATH: &str = "../../repos/PharmCAT/src/main/resources/org/pharmgkb/pharmcat/definition/alleles/RYR1_translation.json";
    const SLCO1B1_DEFINITION_PATH: &str = "../../repos/PharmCAT/src/main/resources/org/pharmgkb/pharmcat/definition/alleles/SLCO1B1_translation.json";
    const TPMT_DEFINITION_PATH: &str = "../../repos/PharmCAT/src/main/resources/org/pharmgkb/pharmcat/definition/alleles/TPMT_translation.json";
    const UGT1A1_DEFINITION_PATH: &str = "../../repos/PharmCAT/src/main/resources/org/pharmgkb/pharmcat/definition/alleles/UGT1A1_translation.json";

    #[test]
    fn unix_millis_to_iso8601_matches_java_instant_shape() {
        assert_eq!(unix_millis_to_iso8601(0), "1970-01-01T00:00:00Z");
        assert_eq!(
            unix_millis_to_iso8601(1_700_000_000_123),
            "2023-11-14T22:13:20.123Z"
        );
    }

    #[test]
    fn matcher_html_date_from_millis_matches_java_short_date_shape() {
        assert_eq!(matcher_html_date_from_millis(0), "01/01/70");
        assert_eq!(matcher_html_date_from_millis(1_700_000_000_123), "11/14/23");
    }

    #[test]
    fn matcher_html_template_matches_java_template_shell_and_println_newline() {
        let html = render_matcher_html_template("Title <raw>", "<h3>Gene</h3>\n", "06/03/26");

        assert!(html.starts_with("<!DOCTYPE html>\n<html class=\"no-js\" lang=\"en\">\n<head>\n"));
        assert!(html.contains("    /* Move down content because we have a fixed navbar */\n"));
        assert!(html.contains("<title>Title <raw></title>"));
        assert!(html.contains("<a class=\"navbar-brand\" href=\"#\">Title <raw></a>"));
        assert!(html.contains("<div class=\"container-fluid\">\n<h3>Gene</h3>\n\n</div>"));
        assert!(html.contains("<small>Generated on 06/03/26.</small>"));
        assert!(html.ends_with("</html>\n\n"));
    }

    #[test]
    fn matcher_sample_props_json_reads_java_triplet_metadata_tsv() {
        let metadata = write_temp_named_file(
            "sample-metadata.tsv",
            "Other\tGroup\tignored\nNA12878\tGroup\tcase\nNA12878\tBatch\tB1\nNA12879\tGroup\tcontrol\n",
        );

        let props = matcher_sample_props_json(Some(&metadata), "NA12878").expect("sample props");

        assert_eq!(props["Group"].as_str(), Some("case"));
        assert_eq!(props["Batch"].as_str(), Some("B1"));
        assert!(
            matcher_sample_props_json(Some(&metadata), "missing")
                .expect("missing sample")
                .is_null()
        );
        assert!(
            matcher_sample_props_json(None, "NA12878")
                .expect("no metadata")
                .is_null()
        );
    }

    #[test]
    fn matcher_warning_json_serializes_java_message_annotation_objects() {
        let missing_required = gene_call_warning_json(
            "CYP2B6",
            &GeneCallWarning::MissingRequiredPosition(vec!["chr19:41006936".to_owned()]),
        );
        assert_eq!(
            missing_required,
            serde_json::json!({
                "rule_name": "missing-required-position",
                "version": null,
                "matches": null,
                "exception_type": "note",
                "message": "Cannot call CYP2B6 - missing required variant (chr19:41006936)",
            })
        );

        let dpyd_exonic_only = dpyd_hap_b3_warning_json(&DpydHapB3Warning::ExonicOnly);
        assert_eq!(
            dpyd_exonic_only["rule_name"].as_str(),
            Some("pcat-dpyd-hapB3-exonic-only")
        );
        assert_eq!(dpyd_exonic_only["version"].as_str(), Some("2"));
        assert_eq!(dpyd_exonic_only["exception_type"].as_str(), Some(""));
        assert_eq!(
            dpyd_exonic_only["matches"]["hapsCalled"]
                .as_array()
                .map(Vec::len),
            Some(0)
        );
        assert!(
            dpyd_exonic_only["message"]
                .as_str()
                .is_some_and(|message| message.contains("rs56038477"))
        );

        let allele_count = dpyd_hap_b3_warning_json(&DpydHapB3Warning::AlleleCount {
            count: 1,
            rsid: Some("rs75017182".to_owned()),
        });
        assert_eq!(allele_count["rule_name"].as_str(), Some("warn.alleleCount"));
        assert_eq!(allele_count["exception_type"].as_str(), Some("note"));
        assert_eq!(
            allele_count["message"].as_str(),
            Some("Only found 1 allele for rs75017182")
        );
    }

    #[test]
    fn match_data_json_uses_java_exposed_fields_and_phase_set_maps() {
        let mut definition =
            read_definition_file(Path::new(CYP3A5_DEFINITION_PATH)).expect("CYP3A5 definition");
        definition.variants.truncate(2);
        let first = &definition.variants[0];
        let second = &definition.variants[1];
        let allele_map = allele_map_from_vcf_records([
            sample_summary(first, "C|T", true, Some(10)),
            sample_summary(second, "A|G", true, None),
        ]);
        let match_data = MatchData::new("NA12878", "CYP3A5", &definition, &allele_map);

        let json = match_data_json(&match_data);

        assert_eq!(json["phased"].as_bool(), Some(true));
        assert_eq!(json["effectivelyPhased"].as_bool(), Some(false));
        assert_eq!(json["homozygous"].as_bool(), Some(false));
        assert_eq!(
            json["treatUndocumentedVariationsAsReference"].as_bool(),
            Some(false)
        );
        assert!(json.get("haploid").is_none());
        assert_eq!(json["phaseSets"]["10"][0].as_u64(), Some(first.position));
        assert_eq!(
            json["phaseSets"][i32::MIN.to_string()][0].as_u64(),
            Some(second.position)
        );
        assert_eq!(
            json["posToPhaseSet"][first.position.to_string()].as_i64(),
            Some(10)
        );
        assert_eq!(
            json["posToPhaseSet"][second.position.to_string()].as_i64(),
            Some(i32::MIN as i64)
        );
    }

    #[test]
    fn named_allele_json_uses_java_exposed_fields() {
        let definition =
            read_definition_file(Path::new(CYP3A5_DEFINITION_PATH)).expect("CYP3A5 definition");
        let allele = definition
            .named_alleles
            .iter()
            .find(|allele| allele.name == "*2")
            .expect("CYP3A5 *2");

        let json = named_allele_json(allele);

        assert_eq!(json["name"].as_str(), Some("*2"));
        assert_eq!(json["id"].as_str(), Some(allele.id.as_str()));
        assert!(json["alleles"].is_array());
        assert!(json["cpicAlleles"].is_array());
        assert_eq!(json["reference"].as_bool(), Some(false));
        assert_eq!(json["structuralVariant"].as_bool(), Some(false));
        assert!(json["corePositions"].is_array());
        assert!(json["score"].as_i64().is_some_and(|score| score > 0));
        assert!(json["numCombinations"].is_number());
        assert!(json["numPartials"].is_number());
        assert!(json.get("populationFrequency").is_none());
        assert!(json.get("missingPositions").is_none());
        assert!(json.get("scoreOverride").is_none());
        assert!(json.get("isCombinationOrPartial").is_none());
    }

    #[test]
    fn uncallable_haplotypes_json_matches_java_definition_minus_matchable_haps() {
        let mut definition =
            read_definition_file(Path::new(CYP3A5_DEFINITION_PATH)).expect("CYP3A5 definition");
        definition.variants.truncate(2);
        for haplotype in &mut definition.named_alleles {
            haplotype.alleles.truncate(2);
            haplotype.cpic_alleles.truncate(2);
        }
        let forced_uncallable = definition
            .named_alleles
            .iter_mut()
            .find(|haplotype| haplotype.name == "*2")
            .expect("CYP3A5 *2");
        forced_uncallable.alleles[0] = None;
        forced_uncallable.cpic_alleles[0] = None;
        let locus = &definition.variants[0];
        let allele_map = allele_map_from_vcf_records([sample_summary(locus, "C|T", true, None)]);
        let mut match_data = MatchData::new("NA12878", "CYP3A5", &definition, &allele_map);
        match_data.marshall_haplotypes(&definition);

        let uncallable = uncallable_haplotypes_json(Some(&definition), &match_data);
        let matchable = match_data
            .haplotypes()
            .iter()
            .map(|haplotype| haplotype.name.clone())
            .collect::<BTreeSet<_>>();
        let all = definition
            .named_alleles
            .iter()
            .map(|haplotype| haplotype.name.clone())
            .collect::<BTreeSet<_>>();

        assert_eq!(
            uncallable,
            all.difference(&matchable).cloned().collect::<BTreeSet<_>>()
        );
        assert!(!uncallable.is_empty());
        assert!(uncallable.contains("*2"));
        assert!(uncallable.is_disjoint(&matchable));
    }

    #[test]
    fn matcher_html_renders_unmatched_haplotypes_when_no_haplotypes_matched_like_java() {
        let mut definition =
            read_definition_file(Path::new(CYP3A5_DEFINITION_PATH)).expect("CYP3A5 definition");
        definition.variants.truncate(1);
        definition.named_alleles.truncate(2);
        for haplotype in &mut definition.named_alleles {
            haplotype.alleles.truncate(1);
            haplotype.cpic_alleles.truncate(1);
        }
        let locus = &definition.variants[0];
        let allele_map = allele_map_from_vcf_records([sample_summary(locus, "C|T", true, None)]);
        let mut match_data = MatchData::new("NA12878", "CYP3A5", &definition, &allele_map);
        match_data.marshall_haplotypes(&definition);
        let result = GeneCallResult {
            gene: "CYP3A5".to_owned(),
            match_data,
            kind: GeneCallKind::NoCall,
            dpyd_hap_b3_warnings: BTreeSet::new(),
            warnings: BTreeSet::new(),
        };

        let mut html = String::new();
        matcher_html_gene_section(
            &mut html,
            &result,
            Some(&definition),
            MatcherHtmlOptions::default(),
        );

        assert!(html.contains("<tr class=\"table-danger\"><th class=\"first\">*1</th>"));
        assert!(html.contains("<tr class=\"table-danger\"><th class=\"first\">*2</th>"));
    }

    #[test]
    fn matcher_html_always_show_unmatched_option_matches_java_serializer_mode() {
        let mut definition =
            read_definition_file(Path::new(CYP3A5_DEFINITION_PATH)).expect("CYP3A5 definition");
        definition.variants.truncate(1);
        definition.named_alleles.truncate(2);
        for haplotype in &mut definition.named_alleles {
            haplotype.alleles.truncate(1);
            haplotype.cpic_alleles.truncate(1);
        }
        let locus = &definition.variants[0];
        let allele_map = allele_map_from_vcf_records([sample_summary(locus, "C|C", true, None)]);
        let mut match_data = MatchData::new("NA12878", "CYP3A5", &definition, &allele_map);
        match_data.marshall_haplotypes(&definition);
        let matched = match_data
            .haplotypes()
            .iter()
            .find(|haplotype| haplotype.name == "*1")
            .expect("CYP3A5 *1")
            .clone();
        let result = GeneCallResult {
            gene: "CYP3A5".to_owned(),
            match_data,
            kind: GeneCallKind::Haplotypes(vec![HaplotypeMatch {
                name: matched.name.clone(),
                haplotype: matched,
                positions: definition.variants[..1].to_vec(),
                sequences: BTreeSet::new(),
            }]),
            dpyd_hap_b3_warnings: BTreeSet::new(),
            warnings: BTreeSet::new(),
        };

        let mut default_html = String::new();
        matcher_html_gene_section(
            &mut default_html,
            &result,
            Some(&definition),
            MatcherHtmlOptions::default(),
        );
        assert!(!default_html.contains("<tr class=\"table-danger\"><th class=\"first\">*2</th>"));

        let mut always_show_html = String::new();
        matcher_html_gene_section(
            &mut always_show_html,
            &result,
            Some(&definition),
            MatcherHtmlOptions {
                always_show_unmatched_haplotypes: true,
            },
        );
        assert!(
            always_show_html.contains("<tr class=\"table-danger\"><th class=\"first\">*2</th>")
        );
    }

    #[test]
    fn matcher_html_haplotype_rows_emit_only_sequence_cells_like_java_print_allele() {
        let definition =
            read_definition_file(Path::new(CYP3A5_DEFINITION_PATH)).expect("CYP3A5 definition");
        let positions = definition.variants[..2].to_vec();
        let mut alleles = BTreeMap::new();
        alleles.insert(positions[1].position, "T".to_owned());
        let highlighted = [positions[1].position].into_iter().collect();

        let mut html = String::new();
        matcher_html_haplotype_row(
            &mut html,
            Some("*sparse"),
            &alleles,
            Some("table-info"),
            &highlighted,
            &positions,
        );

        assert_eq!(
            html,
            "  <tr class=\"table-info\"><th class=\"first\">*sparse</th><td class=\"table-danger\"><b>T</b></td></tr>\n"
        );
    }

    #[test]
    fn matcher_html_haplotype_rows_render_name_and_call_raw_like_java_print_allele() {
        let definition =
            read_definition_file(Path::new(CYP3A5_DEFINITION_PATH)).expect("CYP3A5 definition");
        let positions = definition.variants[..1].to_vec();
        let mut alleles = BTreeMap::new();
        alleles.insert(positions[0].position, "A<T".to_owned());

        let mut html = String::new();
        matcher_html_haplotype_row(
            &mut html,
            Some("*raw<name>"),
            &alleles,
            None,
            &BTreeSet::new(),
            &positions,
        );

        assert_eq!(
            html,
            "  <tr><th class=\"first\">*raw<name></th><td><b>A<T</b></td></tr>\n"
        );
    }

    #[test]
    fn matcher_html_named_allele_patterns_match_java_permutation_cells() {
        let definition =
            read_definition_file(Path::new(CYP3A5_DEFINITION_PATH)).expect("CYP3A5 definition");
        let positions = matcher_html_pattern_positions();
        let mut allele = definition.named_alleles[0].clone();
        allele.alleles = vec![
            Some("R".to_owned()),
            None,
            Some("-".to_owned()),
            Some("AT or GC".to_owned()),
            Some("A\\T".to_owned()),
        ];

        let pattern_cells = named_allele_position_map(&allele, &positions);

        assert_eq!(
            pattern_cells
                .get(&positions[0].position)
                .map(String::as_str),
            Some("[AG]")
        );
        assert_eq!(
            pattern_cells
                .get(&positions[1].position)
                .map(String::as_str),
            Some(".*?")
        );
        assert_eq!(
            pattern_cells
                .get(&positions[2].position)
                .map(String::as_str),
            Some("del")
        );
        assert_eq!(
            pattern_cells
                .get(&positions[3].position)
                .map(String::as_str),
            Some("(AT|GC)")
        );
        assert_eq!(
            pattern_cells
                .get(&positions[4].position)
                .map(String::as_str),
            Some("A\\T")
        );
    }

    #[test]
    fn matcher_html_haplotype_rows_render_java_pattern_cells() {
        let definition =
            read_definition_file(Path::new(CYP3A5_DEFINITION_PATH)).expect("CYP3A5 definition");
        let positions = matcher_html_pattern_positions();
        let mut allele = definition.named_alleles[0].clone();
        allele.alleles = vec![
            Some("R".to_owned()),
            None,
            Some("-".to_owned()),
            Some("AT or GC".to_owned()),
            Some("A\\T".to_owned()),
        ];
        let highlighted = positions
            .iter()
            .map(|position| position.position)
            .collect::<BTreeSet<_>>();

        let mut html = String::new();
        matcher_html_haplotype_row(
            &mut html,
            Some("*pattern"),
            &named_allele_position_map(&allele, &positions),
            Some("table-info"),
            &highlighted,
            &positions,
        );

        assert_eq!(
            html,
            concat!(
                "  <tr class=\"table-info\"><th class=\"first\">*pattern</th>",
                "<td class=\"table-danger\"><b>[AG]</b></td>",
                "<td>.*?</td>",
                "<td class=\"table-danger\"><b>del</b></td>",
                "<td class=\"table-danger\"><b>(AT|GC)</b></td>",
                "<td class=\"table-danger\"><b>AT</b></td>",
                "</tr>\n"
            )
        );
    }

    #[test]
    fn matcher_html_warning_paragraphs_use_java_message_annotation_text() {
        let definition =
            read_definition_file(Path::new(CYP3A5_DEFINITION_PATH)).expect("CYP3A5 definition");
        let result = GeneCallResult {
            gene: "CYP3A5".to_owned(),
            match_data: MatchData::new("NA12878", "CYP3A5", &definition, &BTreeMap::new()),
            kind: GeneCallKind::NoCall,
            dpyd_hap_b3_warnings: BTreeSet::new(),
            warnings: [GeneCallWarning::MissingRequiredPosition(vec![
                "chr7:99270539".to_owned(),
            ])]
            .into_iter()
            .collect(),
        };

        let mut html = String::new();
        matcher_html_gene_section(
            &mut html,
            &result,
            Some(&definition),
            MatcherHtmlOptions::default(),
        );

        assert!(
            html.contains("<p>Cannot call CYP3A5 - missing required variant (chr7:99270539)</p>")
        );
        assert!(!html.contains("Missing required position(s):"));
        assert!(!html.contains("<table class=\"table table-striped table-hover table-sm\">"));
    }

    #[test]
    fn matcher_html_renders_missing_tag_position_notes_for_called_haplotypes_like_java() {
        let mut definition =
            read_definition_file(Path::new(CYP3A5_DEFINITION_PATH)).expect("CYP3A5 definition");
        definition.variants.truncate(2);
        for haplotype in &mut definition.named_alleles {
            haplotype.alleles.truncate(2);
            haplotype.cpic_alleles.truncate(2);
        }
        let locus = &definition.variants[0];
        let allele_map = allele_map_from_vcf_records([sample_summary(locus, "C|T", true, None)]);
        let mut match_data = MatchData::new("NA12878", "CYP3A5", &definition, &allele_map);
        match_data.marshall_haplotypes(&definition);
        let haplotype = match_data
            .haplotypes()
            .iter()
            .find(|haplotype| !haplotype.missing_positions.is_empty())
            .expect("haplotype with missing tag")
            .clone();
        let missing_hgvs = haplotype
            .missing_positions
            .iter()
            .map(|variant| variant.chromosome_hgvs_name.as_str())
            .collect::<Vec<_>>()
            .join(", ");
        let result = GeneCallResult {
            gene: "CYP3A5".to_owned(),
            match_data,
            kind: GeneCallKind::Haplotypes(vec![HaplotypeMatch {
                name: haplotype.name.clone(),
                haplotype,
                positions: definition.variants[..1].to_vec(),
                sequences: BTreeSet::new(),
            }]),
            dpyd_hap_b3_warnings: BTreeSet::new(),
            warnings: BTreeSet::new(),
        };

        let mut html = String::new();
        matcher_html_gene_section(
            &mut html,
            &result,
            Some(&definition),
            MatcherHtmlOptions::default(),
        );

        assert!(html.contains(
            "The following haplotypes were called even though tag positions were missing:"
        ));
        assert!(html.contains(&format!(" without {missing_hgvs}</li>")));
    }

    fn sample_summary(
        locus: &VariantLocus,
        vcf_call: &str,
        phased: bool,
        phase_set: Option<i32>,
    ) -> SampleAlleleSummary {
        let alleles = vcf_call
            .split(['|', '/'])
            .map(str::to_owned)
            .collect::<Vec<_>>();
        SampleAlleleSummary {
            chromosome: locus.chromosome.clone(),
            position: locus.position as usize,
            allele1: alleles.first().cloned(),
            allele2: alleles.get(1).cloned(),
            vcf_alleles: alleles,
            genotype: if phased {
                "0|1".to_owned()
            } else {
                "0/1".to_owned()
            },
            vcf_call: vcf_call.to_owned(),
            phased,
            effectively_phased: phased,
            phase_set,
            undocumented_variations: BTreeSet::new(),
            treat_undocumented_variations_as_reference: false,
        }
    }

    #[test]
    fn pipeline_result_models_java_status_basename_and_sample() {
        assert_eq!(
            PipelineResult::new(PipelineStatus::Success, "sample", Some("NA12878")),
            PipelineResult {
                status: PipelineStatus::Success,
                basename: "sample".to_string(),
                sample_id: Some("NA12878".to_string()),
            }
        );
        assert_eq!(
            PipelineResult::new(PipelineStatus::Noop, "sample", None::<String>).sample_id,
            None
        );
    }

    #[test]
    fn run_plan_matches_java_input_description_and_batch_start_line() {
        let vcf_file = write_temp_named_file("sample.vcf", "##fileformat=VCFv4.3\n");
        let outside_file = write_temp_named_file("sample.outside.tsv", "Gene\tDiplotype\n");
        let action = parse_pharmcat_args([
            "-vcf",
            vcf_file.to_str().unwrap(),
            "-po",
            outside_file.to_str().unwrap(),
            "-del",
            "-v",
        ])
        .expect("parse");
        let CliAction::Run(config) = action else {
            panic!("expected run");
        };

        let plan = PipelineRunPlan::from_cli(
            &config,
            Some("NA12878"),
            false,
            PipelineMode::Cli,
            Some("2 / 3 -"),
        )
        .expect("plan");
        assert_eq!(plan.input_description, "sample.vcf, sample.outside.tsv");
        assert!(plan.batch_display_mode);
        assert_eq!(
            plan.starting_line,
            Some(
                "+ 2 / 3 - Starting sample NA12878 in sample.vcf (inputs: sample.vcf, sample.outside.tsv)"
                    .to_string()
            )
        );
        assert_eq!(
            plan.intermediate_files_to_delete,
            vec![
                plan.outputs.base_dir.join("sample.NA12878.match.json"),
                plan.outputs.base_dir.join("sample.NA12878.phenotype.json"),
            ]
        );
        assert_eq!(
            plan.error_file,
            plan.outputs.base_dir.join("sample.NA12878.ERROR.txt")
        );
        assert!(plan.save_messages.is_empty());
        assert_eq!(
            plan.finished_block,
            Some("- 2 / 3 - Finished processing sample NA12878 in sample.vcf\n".to_string())
        );
    }

    #[test]
    fn run_plan_omits_batch_start_line_for_single_sample_cli_like_java() {
        let reporter_input = write_temp_named_file("sample.phenotype.json", "{}\n");
        let action = parse_pharmcat_args([
            "-reporter",
            "-ri",
            reporter_input.to_str().unwrap(),
            "-reporterJson",
        ])
        .expect("parse");
        let CliAction::Run(config) = action else {
            panic!("expected run");
        };

        let plan =
            PipelineRunPlan::from_cli(&config, None, true, PipelineMode::Cli, None).expect("plan");
        assert_eq!(plan.input_description, "sample.phenotype.json");
        assert!(!plan.batch_display_mode);
        assert_eq!(plan.starting_line, None);
        assert!(plan.intermediate_files_to_delete.is_empty());
        assert_eq!(
            plan.error_file,
            plan.outputs.base_dir.join("sample.ERROR.txt")
        );
        assert_eq!(
            plan.save_messages,
            vec![
                "Saving reporter JSON results to ".to_string()
                    + &plan
                        .outputs
                        .base_dir
                        .join("sample.report.json")
                        .display()
                        .to_string()
            ]
        );
        assert_eq!(plan.finished_block, None);
    }

    #[test]
    fn save_messages_match_java_single_sample_cli_output_and_delete_rules() {
        let vcf_file = write_temp_named_file("sample.vcf", "##fileformat=VCFv4.3\n");
        let action = parse_pharmcat_args([
            "-vcf",
            vcf_file.to_str().unwrap(),
            "-matcherHtml",
            "-reporterJson",
            "-reporterHtml",
            "-reporterCallsOnlyTsv",
        ])
        .expect("parse");
        let CliAction::Run(config) = action else {
            panic!("expected run");
        };
        let plan =
            PipelineRunPlan::from_cli(&config, None, true, PipelineMode::Cli, None).expect("plan");

        assert_eq!(
            plan.save_messages,
            vec![
                format!(
                    "Saving named allele matcher JSON results to {}",
                    plan.outputs.base_dir.join("sample.match.json").display()
                ),
                format!(
                    "Saving named allele matcher HTML results to {}",
                    plan.outputs.base_dir.join("sample.match.html").display()
                ),
                format!(
                    "Saving phenotyper JSON results to {}",
                    plan.outputs
                        .base_dir
                        .join("sample.phenotype.json")
                        .display()
                ),
                format!(
                    "Saving reporter HTML results to {}",
                    plan.outputs.base_dir.join("sample.report.html").display()
                ),
                format!(
                    "Saving reporter JSON results to {}",
                    plan.outputs.base_dir.join("sample.report.json").display()
                ),
                format!(
                    "Saving calls-only TSV results to {}",
                    plan.outputs.base_dir.join("sample.report.tsv").display()
                ),
            ]
        );

        let action = parse_pharmcat_args([
            "-vcf",
            vcf_file.to_str().unwrap(),
            "-matcherHtml",
            "-reporterJson",
            "-reporterHtml",
            "-reporterCallsOnlyTsv",
            "-del",
        ])
        .expect("parse");
        let CliAction::Run(config) = action else {
            panic!("expected run");
        };
        let plan =
            PipelineRunPlan::from_cli(&config, None, true, PipelineMode::Cli, None).expect("plan");
        assert_eq!(
            plan.save_messages,
            vec![
                format!(
                    "Saving reporter HTML results to {}",
                    plan.outputs.base_dir.join("sample.report.html").display()
                ),
                format!(
                    "Saving reporter JSON results to {}",
                    plan.outputs.base_dir.join("sample.report.json").display()
                ),
                format!(
                    "Saving calls-only TSV results to {}",
                    plan.outputs.base_dir.join("sample.report.tsv").display()
                ),
            ]
        );
    }

    #[test]
    fn delete_intermediate_files_matches_java_delete_if_exists_semantics() {
        let vcf_file = write_temp_named_file("sample.vcf", "##fileformat=VCFv4.3\n");
        let output_dir = unique_temp_path("pharmcat-cleanup-out");
        let action = parse_pharmcat_args([
            "-vcf",
            vcf_file.to_str().unwrap(),
            "-o",
            output_dir.to_str().unwrap(),
            "-del",
        ])
        .expect("parse");
        let CliAction::Run(config) = action else {
            panic!("expected run");
        };
        let plan =
            PipelineRunPlan::from_cli(&config, None, true, PipelineMode::Cli, None).expect("plan");

        fs::write(plan.outputs.base_dir.join("sample.match.json"), "{}\n").expect("write matcher");
        let deleted = delete_intermediate_files(&plan).expect("delete intermediates");
        assert_eq!(
            deleted,
            vec![plan.outputs.base_dir.join("sample.match.json")]
        );
        assert!(!plan.outputs.base_dir.join("sample.match.json").exists());
        assert!(!plan.outputs.base_dir.join("sample.phenotype.json").exists());

        let deleted = delete_intermediate_files(&plan).expect("second delete is ignored");
        assert!(deleted.is_empty());
    }

    #[test]
    fn write_error_file_uses_java_error_path_from_run_plan() {
        let reporter_input = write_temp_named_file("sample.phenotype.json", "{}\n");
        let action = parse_pharmcat_args([
            "-reporter",
            "-ri",
            reporter_input.to_str().unwrap(),
            "-reporterJson",
        ])
        .expect("parse");
        let CliAction::Run(config) = action else {
            panic!("expected run");
        };
        let plan =
            PipelineRunPlan::from_cli(&config, None, true, PipelineMode::Cli, None).expect("plan");
        let error = io::Error::other("boom");

        write_error_file(&plan.error_file, &error).expect("write error file");
        let contents = fs::read_to_string(&plan.error_file).expect("read error file");
        assert_eq!(
            plan.error_file,
            plan.outputs.base_dir.join("sample.ERROR.txt")
        );
        assert!(contents.contains("boom"));
    }

    #[test]
    fn run_reporter_from_vcf_builds_context_and_writes_reporter_outputs_like_java_pipeline_slice() {
        let definition =
            read_definition_file(Path::new(CYP3A5_DEFINITION_PATH)).expect("definition");
        let definitions = DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition)]
                .into_iter()
                .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");
        let output_dir = unique_temp_path("pharmcat-run-reporter-from-vcf");
        let outputs = PipelineOutputPlan {
            base_dir: output_dir.clone(),
            basename: "haplotyper".to_owned(),
            display_name: "haplotyper.vcf".to_owned(),
            reporter_title: Some("NA12878".to_owned()),
            matcher_json: None,
            matcher_html: None,
            matcher_warnings: None,
            phenotyper_json: None,
            reporter_html: Some(output_dir.join("haplotyper.report.html")),
            reporter_json: Some(output_dir.join("haplotyper.report.json")),
            reporter_calls_only_tsv: Some(output_dir.join("haplotyper.report.tsv")),
        };
        let options = ReporterPipelineOptions {
            include_combinations: false,
            ..ReporterPipelineOptions::default()
        };

        let run = run_reporter_from_vcf(
            Path::new(CYP3A5_VCF_PATH),
            Some("NA12878"),
            &definitions,
            &phenotypes,
            &guidance,
            Some(&outputs),
            &options,
        )
        .expect("pipeline run");

        assert_eq!(run.gene_call_results.len(), 1);
        assert_eq!(run.written_outputs.len(), 3);
        assert!(
            run.written_outputs
                .contains(outputs.reporter_json.as_ref().unwrap())
        );
        assert!(
            run.written_outputs
                .contains(outputs.reporter_html.as_ref().unwrap())
        );
        assert!(
            run.written_outputs
                .contains(outputs.reporter_calls_only_tsv.as_ref().unwrap())
        );

        let cyp3a5 = run.context.gene_report("CYP3A5").expect("CYP3A5 report");
        assert_eq!(cyp3a5.source_diplotype.as_deref(), Some("*1/*2"));
        assert_eq!(cyp3a5.chromosome.as_deref(), Some("chr1"));
        assert!(!cyp3a5.variant_reports.is_empty());

        let json = fs::read_to_string(outputs.reporter_json.as_ref().unwrap()).expect("json");
        assert!(json.contains("\"geneSymbol\": \"CYP3A5\""));
        let html = fs::read_to_string(outputs.reporter_html.as_ref().unwrap()).expect("html");
        assert!(html.contains("PharmCAT Report [NA12878]"));
        let tsv =
            fs::read_to_string(outputs.reporter_calls_only_tsv.as_ref().unwrap()).expect("tsv");
        assert!(tsv.contains("CYP3A5"));
    }

    #[test]
    fn run_reporter_from_empty_vcf_writes_no_data_report_like_java_pipeline_test() {
        let definition =
            read_definition_file(Path::new(CYP3A5_DEFINITION_PATH)).expect("definition");
        let definitions = DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition)]
                .into_iter()
                .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");
        let vcf_file = write_temp_named_file(
            "pipeline-no-data.vcf",
            "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n",
        );
        let output_dir = unique_temp_path("pharmcat-run-no-data");
        let outputs = PipelineOutputPlan {
            base_dir: output_dir.clone(),
            basename: "pipeline-no-data".to_owned(),
            display_name: "pipeline-no-data.vcf".to_owned(),
            reporter_title: Some("PharmCAT".to_owned()),
            matcher_json: None,
            matcher_html: None,
            matcher_warnings: None,
            phenotyper_json: None,
            reporter_html: Some(output_dir.join("pipeline-no-data.report.html")),
            reporter_json: Some(output_dir.join("pipeline-no-data.report.json")),
            reporter_calls_only_tsv: Some(output_dir.join("pipeline-no-data.report.tsv")),
        };

        let run = run_reporter_from_vcf(
            &vcf_file,
            Some("PharmCAT"),
            &definitions,
            &phenotypes,
            &guidance,
            Some(&outputs),
            &ReporterPipelineOptions {
                include_combinations: false,
                ..ReporterPipelineOptions::default()
            },
        )
        .expect("no-data pipeline run");

        assert_eq!(run.gene_call_results.len(), 1);
        assert!(run.vcf_warnings.is_empty());
        assert_eq!(run.written_outputs.len(), 3);

        let gene = run.context.gene_report("CYP3A5").expect("CYP3A5 report");
        assert_eq!(gene.source_diplotype.as_deref(), Some("Unknown/Unknown"));
        assert!(!gene.variant_reports.is_empty());
        assert!(
            gene.variant_reports
                .iter()
                .all(|variant| variant.is_missing())
        );

        let json = fs::read_to_string(outputs.reporter_json.as_ref().unwrap()).expect("json");
        let report_json: serde_json::Value = serde_json::from_str(&json).expect("report json");
        let cyp3a5_json = &report_json["genes"]["CYP3A5"];
        assert_eq!(cyp3a5_json["geneSymbol"].as_str(), Some("CYP3A5"));
        assert_eq!(
            cyp3a5_json["sourceDiplotypes"][0]["label"].as_str(),
            Some("Unknown/Unknown")
        );
        let html = fs::read_to_string(outputs.reporter_html.as_ref().unwrap()).expect("html");
        assert!(html.contains("PharmCAT Report [PharmCAT]"));
        assert!(html.contains("<p>No data provided.</p>"));
        assert!(!html.contains("<table class=\"genotypeSummary\">"));
        let tsv =
            fs::read_to_string(outputs.reporter_calls_only_tsv.as_ref().unwrap()).expect("tsv");
        assert!(tsv.contains("CYP3A5"));
        assert!(tsv.contains("Unknown/Unknown"));
    }

    #[test]
    fn run_reporter_from_vcf_marks_cyp2c19_undocumented_variation_uncallable_like_java_pipeline_test()
     {
        let (definition, outputs, run) =
            run_cyp2c19_undocumented_variation_pipeline(true, "pipeline-cyp2c19-undocumented");

        assert_eq!(run.gene_call_results.len(), 1);
        let warning = run
            .vcf_warnings
            .get(&variant_by_rsid(&definition, "rs3758581").vcf_chr_position())
            .expect("undocumented variation warning")
            .iter()
            .next()
            .expect("warning text");
        assert!(warning.contains("expected G"));
        assert!(warning.contains("found T in VCF"));

        let gene = run.context.gene_report("CYP2C19").expect("CYP2C19 report");
        assert_eq!(gene.source_diplotype.as_deref(), Some("Unknown/Unknown"));
        assert!(gene.has_undocumented_variations);
        assert!(!gene.treat_undocumented_variations_as_reference);
        let undocumented = gene
            .variant_reports
            .iter()
            .find(|variant| variant.db_snp_id.as_deref() == Some("rs3758581"))
            .expect("rs3758581 report");
        assert!(undocumented.has_undocumented_variations);
        assert_eq!(undocumented.call.as_deref(), Some("A/T"));

        let html = fs::read_to_string(outputs.reporter_html.as_ref().unwrap()).expect("html");
        assert!(html.contains("id=\"gs-uncallable-CYP2C19\""));
        assert!(html.contains("<td class=\"top-aligned genotype-result\">Not called</td>"));
        assert!(html.contains("<div class=\"callMessage\">Undocumented variation</div>"));
        assert!(html.contains("CYP2C19 allele match data"));
        assert!(html.contains("<p class=\"rx-no-recs\">No recommendations.</p>"));
        assert!(!html.contains("<span class=\"rx-no-call\">No call data for CYP2C19</span>"));

        let json = fs::read_to_string(outputs.reporter_json.as_ref().unwrap()).expect("json");
        let report_json: serde_json::Value = serde_json::from_str(&json).expect("report json");
        let cyp2c19_json = &report_json["genes"]["CYP2C19"];
        assert_eq!(
            cyp2c19_json["sourceDiplotypes"][0]["label"].as_str(),
            Some("Unknown/Unknown")
        );
        assert_eq!(
            cyp2c19_json["hasUndocumentedVariations"].as_bool(),
            Some(true)
        );
    }

    #[test]
    fn run_reporter_from_vcf_marks_cyp2c19_undocumented_variation_uncallable_in_extended_report_like_java_pipeline_test()
     {
        let (_definition, outputs, run) = run_cyp2c19_undocumented_variation_pipeline(
            false,
            "pipeline-cyp2c19-undocumented-extended",
        );

        assert_eq!(run.gene_call_results.len(), 1);
        let gene = run.context.gene_report("CYP2C19").expect("CYP2C19 report");
        assert_eq!(gene.source_diplotype.as_deref(), Some("Unknown/Unknown"));
        assert!(gene.has_undocumented_variations);

        let html = fs::read_to_string(outputs.reporter_html.as_ref().unwrap()).expect("html");
        assert!(html.contains("id=\"gs-uncallable-CYP2C19\""));
        assert!(html.contains("<td class=\"top-aligned genotype-result\">Not called</td>"));
        assert!(html.contains("<div class=\"callMessage\">Undocumented variation</div>"));
        assert!(html.contains("CYP2C19 allele match data"));
        assert!(!html.contains("<p class=\"rx-no-recs\">No recommendations.</p>"));
        assert!(html.contains("<section class=\"guideline drugReport amitriptyline\">"));
        assert!(html.contains("<span class=\"rx-no-call\">No call data for CYP2C19</span>"));

        let json = fs::read_to_string(outputs.reporter_json.as_ref().unwrap()).expect("json");
        let report_json: serde_json::Value = serde_json::from_str(&json).expect("report json");
        let cyp2c19_json = &report_json["genes"]["CYP2C19"];
        assert_eq!(
            cyp2c19_json["sourceDiplotypes"][0]["label"].as_str(),
            Some("Unknown/Unknown")
        );
        assert_eq!(
            cyp2c19_json["hasUndocumentedVariations"].as_bool(),
            Some(true)
        );
    }

    #[test]
    fn run_reporter_from_vcf_treats_toxic_gene_undocumented_variations_as_reference_like_java_pipeline_test()
     {
        let cyp2c19_definition =
            read_definition_file(Path::new(CYP2C19_DEFINITION_PATH)).expect("CYP2C19 definition");
        let ryr1_definition =
            read_definition_file(Path::new(RYR1_DEFINITION_PATH)).expect("RYR1 definition");
        let tpmt_definition =
            read_definition_file(Path::new(TPMT_DEFINITION_PATH)).expect("TPMT definition");
        let definitions = DefinitionReader::from_definitions(
            [
                (
                    cyp2c19_definition.gene_symbol.clone(),
                    cyp2c19_definition.clone(),
                ),
                (ryr1_definition.gene_symbol.clone(), ryr1_definition.clone()),
                (tpmt_definition.gene_symbol.clone(), tpmt_definition.clone()),
            ]
            .into_iter()
            .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");

        let mut vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n".to_owned();
        append_definition_vcf_rows(
            &mut vcf,
            &cyp2c19_definition,
            &[("rs3758581", "G,T", "0/2")],
            &[],
        );
        append_definition_vcf_rows(
            &mut vcf,
            &ryr1_definition,
            &[("rs193922753", "A,T,C", "0/3")],
            &[],
        );
        append_definition_vcf_rows(
            &mut vcf,
            &tpmt_definition,
            &[("rs1800462", "G,T", "0/2")],
            &[],
        );
        let vcf_file = write_temp_named_file("pipeline-undocumented-treat-as-reference.vcf", &vcf);
        let output_dir = unique_temp_path("pharmcat-run-undocumented-treat-as-reference");
        let outputs = PipelineOutputPlan {
            base_dir: output_dir.clone(),
            basename: "pipeline-undocumented-treat-as-reference".to_owned(),
            display_name: "pipeline-undocumented-treat-as-reference.vcf".to_owned(),
            reporter_title: Some("PharmCAT".to_owned()),
            matcher_json: None,
            matcher_html: None,
            matcher_warnings: None,
            phenotyper_json: None,
            reporter_html: Some(
                output_dir.join("pipeline-undocumented-treat-as-reference.report.html"),
            ),
            reporter_json: Some(
                output_dir.join("pipeline-undocumented-treat-as-reference.report.json"),
            ),
            reporter_calls_only_tsv: Some(
                output_dir.join("pipeline-undocumented-treat-as-reference.report.tsv"),
            ),
        };

        let run = run_reporter_from_vcf(
            &vcf_file,
            Some("PharmCAT"),
            &definitions,
            &phenotypes,
            &guidance,
            Some(&outputs),
            &ReporterPipelineOptions {
                include_combinations: false,
                html: HtmlReportOptions {
                    compact: true,
                    ..HtmlReportOptions::default()
                },
                ..ReporterPipelineOptions::default()
            },
        )
        .expect("undocumented treat-as-reference pipeline run");

        assert_eq!(run.gene_call_results.len(), 3);
        assert_undocumented_warning(
            &run.vcf_warnings,
            &variant_by_rsid(&cyp2c19_definition, "rs3758581").vcf_chr_position(),
            "expected G",
            "found T in VCF",
            false,
        );
        assert_undocumented_warning(
            &run.vcf_warnings,
            &variant_by_rsid(&ryr1_definition, "rs193922753").vcf_chr_position(),
            "expected A/T",
            "found C in VCF",
            true,
        );
        assert_undocumented_warning(
            &run.vcf_warnings,
            &variant_by_rsid(&tpmt_definition, "rs1800462").vcf_chr_position(),
            "expected G",
            "found T in VCF",
            true,
        );

        let cyp2c19 = run.context.gene_report("CYP2C19").expect("CYP2C19 report");
        assert_eq!(cyp2c19.source_diplotype.as_deref(), Some("Unknown/Unknown"));
        assert!(cyp2c19.has_undocumented_variations);
        assert!(!cyp2c19.treat_undocumented_variations_as_reference);

        let ryr1 = run.context.gene_report("RYR1").expect("RYR1 report");
        assert_eq!(
            ryr1.source_diplotype.as_deref(),
            Some("Reference/Reference")
        );
        assert!(ryr1.has_undocumented_variations);
        assert!(ryr1.treat_undocumented_variations_as_reference);
        assert_eq!(
            ryr1.recommendation_diplotypes[0].label,
            "Reference/Reference"
        );

        let tpmt = run.context.gene_report("TPMT").expect("TPMT report");
        assert_eq!(tpmt.source_diplotype.as_deref(), Some("*1/*1"));
        assert!(tpmt.has_undocumented_variations);
        assert!(tpmt.treat_undocumented_variations_as_reference);
        assert_eq!(tpmt.recommendation_diplotypes[0].label, "*1/*1");

        let html = fs::read_to_string(outputs.reporter_html.as_ref().unwrap()).expect("html");
        assert!(html.contains("id=\"gs-uncallable-CYP2C19\""));
        assert!(html.contains("id=\"gs-undocVarAsRef-TPMT\""));
        assert!(html.contains("id=\"gs-undocVarAsRef-RYR1\""));
        assert!(html.contains("<tr class=\"top-aligned gs-TPMT\""));
        assert!(html.contains("<tr class=\"top-aligned gs-RYR1\""));
        assert!(html.contains("<td class=\"top-aligned genotype-result\">*1/*1"));
        assert!(html.contains("<td class=\"top-aligned genotype-result\">Reference/Reference"));
        assert!(html.contains("<td class=\"top-aligned genotype-result\">Not called</td>"));

        let json = fs::read_to_string(outputs.reporter_json.as_ref().unwrap()).expect("json");
        let report_json: serde_json::Value = serde_json::from_str(&json).expect("report json");
        assert_eq!(
            report_json["genes"]["CYP2C19"]["sourceDiplotypes"][0]["label"].as_str(),
            Some("Unknown/Unknown")
        );
        assert_eq!(
            report_json["genes"]["TPMT"]["sourceDiplotypes"][0]["label"].as_str(),
            Some("*1/*1")
        );
        assert_eq!(
            report_json["genes"]["RYR1"]["sourceDiplotypes"][0]["label"].as_str(),
            Some("Reference/Reference")
        );
        assert_eq!(
            report_json["genes"]["TPMT"]["treatUndocumentedVariationsAsReference"].as_bool(),
            Some(true)
        );
        assert_eq!(
            report_json["genes"]["RYR1"]["treatUndocumentedVariationsAsReference"].as_bool(),
            Some(true)
        );
        assert_eq!(
            report_json["genes"]["CYP2C19"]["treatUndocumentedVariationsAsReference"].as_bool(),
            Some(false)
        );
    }

    #[test]
    fn run_reporter_from_vcf_combo_keeps_standard_toxic_gene_undocumented_variation_like_java_pipeline_test()
     {
        let ryr1_definition =
            read_definition_file(Path::new(RYR1_DEFINITION_PATH)).expect("RYR1 definition");
        let tpmt_definition =
            read_definition_file(Path::new(TPMT_DEFINITION_PATH)).expect("TPMT definition");
        let definitions = DefinitionReader::from_definitions(
            [
                (ryr1_definition.gene_symbol.clone(), ryr1_definition.clone()),
                (tpmt_definition.gene_symbol.clone(), tpmt_definition.clone()),
            ]
            .into_iter()
            .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");

        let mut vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n".to_owned();
        append_definition_vcf_rows(
            &mut vcf,
            &ryr1_definition,
            &[("rs193922753", "A,T,C", "0/3")],
            &[],
        );
        append_definition_vcf_rows(
            &mut vcf,
            &tpmt_definition,
            &[("rs1800462", "G,T", "0/2")],
            &[],
        );
        let vcf_file =
            write_temp_named_file("pipeline-undocumented-treat-as-reference-combo.vcf", &vcf);
        let output_dir = unique_temp_path("pharmcat-run-undocumented-treat-as-reference-combo");
        let outputs = PipelineOutputPlan {
            base_dir: output_dir.clone(),
            basename: "pipeline-undocumented-treat-as-reference-combo".to_owned(),
            display_name: "pipeline-undocumented-treat-as-reference-combo.vcf".to_owned(),
            reporter_title: Some("PharmCAT".to_owned()),
            matcher_json: None,
            matcher_html: None,
            matcher_warnings: None,
            phenotyper_json: None,
            reporter_html: Some(
                output_dir.join("pipeline-undocumented-treat-as-reference-combo.report.html"),
            ),
            reporter_json: Some(
                output_dir.join("pipeline-undocumented-treat-as-reference-combo.report.json"),
            ),
            reporter_calls_only_tsv: Some(
                output_dir.join("pipeline-undocumented-treat-as-reference-combo.report.tsv"),
            ),
        };

        let run = run_reporter_from_vcf(
            &vcf_file,
            Some("PharmCAT"),
            &definitions,
            &phenotypes,
            &guidance,
            Some(&outputs),
            &ReporterPipelineOptions {
                include_combinations: true,
                html: HtmlReportOptions {
                    compact: true,
                    ..HtmlReportOptions::default()
                },
                ..ReporterPipelineOptions::default()
            },
        )
        .expect("undocumented treat-as-reference combo pipeline run");

        assert_eq!(run.gene_call_results.len(), 2);
        assert_undocumented_warning(
            &run.vcf_warnings,
            &variant_by_rsid(&ryr1_definition, "rs193922753").vcf_chr_position(),
            "expected A/T",
            "found C in VCF",
            true,
        );
        assert_undocumented_warning(
            &run.vcf_warnings,
            &variant_by_rsid(&tpmt_definition, "rs1800462").vcf_chr_position(),
            "expected G",
            "found T in VCF",
            false,
        );

        let ryr1 = run.context.gene_report("RYR1").expect("RYR1 report");
        assert_eq!(
            ryr1.source_diplotype.as_deref(),
            Some("Reference/Reference")
        );
        assert!(ryr1.has_undocumented_variations);
        assert!(ryr1.treat_undocumented_variations_as_reference);
        assert_eq!(
            ryr1.recommendation_diplotypes[0].label,
            "Reference/Reference"
        );

        let tpmt = run.context.gene_report("TPMT").expect("TPMT report");
        assert_eq!(tpmt.source_diplotype.as_deref(), Some("*1/g.18143724C>T"));
        assert!(tpmt.has_undocumented_variations);
        assert!(!tpmt.treat_undocumented_variations_as_reference);
        assert_eq!(tpmt.recommendation_diplotypes[0].label, "*1/g.18143724C>T");
        assert_eq!(tpmt.phenotypes, ["n/a"]);

        let html = fs::read_to_string(outputs.reporter_html.as_ref().unwrap()).expect("html");
        assert!(!html.contains("id=\"gs-undocVarAsRef-TPMT\""));
        assert!(html.contains("id=\"gs-undocVarAsRef-RYR1\""));
        assert!(html.contains("<tr class=\"top-aligned gs-TPMT\""));
        assert!(html.contains("<tr class=\"top-aligned gs-RYR1\""));
        assert!(html.contains("*1/g.18143724C&gt;T"));
        assert!(html.contains("<td class=\"top-aligned genotype-result\">Reference/Reference"));

        let json = fs::read_to_string(outputs.reporter_json.as_ref().unwrap()).expect("json");
        let report_json: serde_json::Value = serde_json::from_str(&json).expect("report json");
        assert_eq!(
            report_json["genes"]["TPMT"]["sourceDiplotypes"][0]["label"].as_str(),
            Some("*1/g.18143724C>T")
        );
        assert_eq!(
            report_json["genes"]["RYR1"]["sourceDiplotypes"][0]["label"].as_str(),
            Some("Reference/Reference")
        );
        assert_eq!(
            report_json["genes"]["TPMT"]["treatUndocumentedVariationsAsReference"].as_bool(),
            Some(false)
        );
        assert_eq!(
            report_json["genes"]["RYR1"]["treatUndocumentedVariationsAsReference"].as_bool(),
            Some(true)
        );
    }

    #[test]
    fn run_reporter_from_vcf_marks_cyp2c19_and_tpmt_uncallable_like_java_pipeline_test() {
        let abcg2_definition =
            read_definition_file(Path::new(ABCG2_DEFINITION_PATH)).expect("ABCG2 definition");
        let cyp2c19_definition =
            read_definition_file(Path::new(CYP2C19_DEFINITION_PATH)).expect("CYP2C19 definition");
        let tpmt_definition =
            read_definition_file(Path::new(TPMT_DEFINITION_PATH)).expect("TPMT definition");
        let definitions = DefinitionReader::from_definitions(
            [
                (
                    abcg2_definition.gene_symbol.clone(),
                    abcg2_definition.clone(),
                ),
                (
                    cyp2c19_definition.gene_symbol.clone(),
                    cyp2c19_definition.clone(),
                ),
                (tpmt_definition.gene_symbol.clone(), tpmt_definition.clone()),
            ]
            .into_iter()
            .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");

        let mut vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n".to_owned();
        append_definition_vcf_rows(&mut vcf, &abcg2_definition, &[], &[]);
        append_definition_vcf_rows(
            &mut vcf,
            &cyp2c19_definition,
            &[("rs3758581", "G,T", "0/2")],
            &[],
        );
        append_definition_vcf_rows(
            &mut vcf,
            &tpmt_definition,
            &[("rs1256618794", "A", "1/1"), ("rs753545734", "T", "0/0")],
            &[],
        );
        let vcf_file = write_temp_named_file("pipeline-uncallable.vcf", &vcf);
        let output_dir = unique_temp_path("pharmcat-run-uncallable");
        let outputs = PipelineOutputPlan {
            base_dir: output_dir.clone(),
            basename: "pipeline-uncallable".to_owned(),
            display_name: "pipeline-uncallable.vcf".to_owned(),
            reporter_title: Some("PharmCAT".to_owned()),
            matcher_json: None,
            matcher_html: None,
            matcher_warnings: None,
            phenotyper_json: None,
            reporter_html: Some(output_dir.join("pipeline-uncallable.report.html")),
            reporter_json: Some(output_dir.join("pipeline-uncallable.report.json")),
            reporter_calls_only_tsv: Some(output_dir.join("pipeline-uncallable.report.tsv")),
        };

        let run = run_reporter_from_vcf(
            &vcf_file,
            Some("PharmCAT"),
            &definitions,
            &phenotypes,
            &guidance,
            Some(&outputs),
            &ReporterPipelineOptions {
                include_combinations: false,
                html: HtmlReportOptions {
                    compact: true,
                    ..HtmlReportOptions::default()
                },
                ..ReporterPipelineOptions::default()
            },
        )
        .expect("uncallable pipeline run");

        assert_eq!(run.gene_call_results.len(), 3);
        assert_undocumented_warning(
            &run.vcf_warnings,
            &variant_by_rsid(&cyp2c19_definition, "rs3758581").vcf_chr_position(),
            "expected G",
            "found T in VCF",
            false,
        );

        let cyp2c19 = run.context.gene_report("CYP2C19").expect("CYP2C19 report");
        assert_eq!(cyp2c19.source_diplotype.as_deref(), Some("Unknown/Unknown"));
        assert!(cyp2c19.has_undocumented_variations);
        assert!(!cyp2c19.treat_undocumented_variations_as_reference);

        let tpmt = run.context.gene_report("TPMT").expect("TPMT report");
        assert_eq!(tpmt.source_diplotype.as_deref(), Some("Unknown/Unknown"));
        assert!(!tpmt.has_undocumented_variations);
        assert!(!tpmt.treat_undocumented_variations_as_reference);

        let abcg2 = run.context.gene_report("ABCG2").expect("ABCG2 report");
        assert_ne!(
            abcg2.source_diplotype.as_deref(),
            Some("Unknown/Unknown"),
            "ABCG2 reference rows should remain callable"
        );

        let html = fs::read_to_string(outputs.reporter_html.as_ref().unwrap()).expect("html");
        assert!(html.contains("id=\"gs-uncallable-CYP2C19\""));
        assert!(html.contains("id=\"gs-uncallable-TPMT\""));
        assert!(html.contains("<section class=\"gene CYP2C19\">"));
        assert!(html.contains("<section class=\"gene TPMT\">"));
        assert!(html.contains("<td class=\"top-aligned genotype-result\">Not called</td>"));

        let tsv =
            fs::read_to_string(outputs.reporter_calls_only_tsv.as_ref().unwrap()).expect("tsv");
        assert!(tsv.contains("CYP2C19"));
        assert!(tsv.contains("TPMT"));
        assert!(tsv.contains("Unknown/Unknown"));

        let json = fs::read_to_string(outputs.reporter_json.as_ref().unwrap()).expect("json");
        let report_json: serde_json::Value = serde_json::from_str(&json).expect("report json");
        assert_eq!(
            report_json["genes"]["CYP2C19"]["sourceDiplotypes"][0]["label"].as_str(),
            Some("Unknown/Unknown")
        );
        assert_eq!(
            report_json["genes"]["TPMT"]["sourceDiplotypes"][0]["label"].as_str(),
            Some("Unknown/Unknown")
        );
    }

    #[test]
    fn run_reporter_from_vcf_cyp2c19_with_cyp2d6_and_g6pd_outside_calls_like_java_pipeline_test() {
        let definition =
            read_definition_file(Path::new(CYP2C19_DEFINITION_PATH)).expect("CYP2C19 definition");
        let definitions = DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition.clone())]
                .into_iter()
                .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");

        let mut vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n".to_owned();
        append_definition_vcf_rows(&mut vcf, &definition, &[("rs3758581", "G", "1/1")], &[]);
        let vcf_file = write_temp_named_file("pipeline-cyp2c19.vcf", &vcf);
        let outside_file = write_temp_named_file(
            "pipeline-cyp2c19.outside.tsv",
            "CYP2D6\t*3/*4\nG6PD\tB (wildtype)/B (wildtype)\n",
        );
        let output_dir = unique_temp_path("pharmcat-run-cyp2c19");
        let outputs = PipelineOutputPlan {
            base_dir: output_dir.clone(),
            basename: "pipeline-cyp2c19".to_owned(),
            display_name: "pipeline-cyp2c19.vcf".to_owned(),
            reporter_title: Some("PharmCAT".to_owned()),
            matcher_json: None,
            matcher_html: None,
            matcher_warnings: None,
            phenotyper_json: None,
            reporter_html: Some(output_dir.join("pipeline-cyp2c19.report.html")),
            reporter_json: Some(output_dir.join("pipeline-cyp2c19.report.json")),
            reporter_calls_only_tsv: Some(output_dir.join("pipeline-cyp2c19.report.tsv")),
        };

        let run = run_reporter_from_vcf(
            &vcf_file,
            Some("PharmCAT"),
            &definitions,
            &phenotypes,
            &guidance,
            Some(&outputs),
            &ReporterPipelineOptions {
                include_combinations: false,
                outside_call_files: vec![outside_file],
                html: HtmlReportOptions {
                    compact: true,
                    ..HtmlReportOptions::default()
                },
                ..ReporterPipelineOptions::default()
            },
        )
        .expect("CYP2C19 pipeline run with outside calls");

        assert_eq!(run.gene_call_results.len(), 1);
        let cyp2c19 = run.context.gene_report("CYP2C19").expect("CYP2C19 report");
        assert_eq!(cyp2c19.source_diplotype.as_deref(), Some("*1/*1"));
        assert!(
            cyp2c19
                .recommendation_diplotypes
                .iter()
                .any(|diplotype| diplotype.label == "*1/*1")
        );

        let cyp2d6 = run.context.gene_report("CYP2D6").expect("CYP2D6 report");
        assert!(cyp2d6.outside_call);
        assert_eq!(cyp2d6.source_diplotype.as_deref(), Some("*3/*4"));
        let g6pd = run.context.gene_report("G6PD").expect("G6PD report");
        assert!(g6pd.outside_call);
        assert_eq!(
            g6pd.source_diplotype.as_deref(),
            Some("B (wildtype)/B (wildtype)")
        );

        for (drug, source) in [
            ("amitriptyline", PrescribingGuidanceSource::CpicGuideline),
            ("amitriptyline", PrescribingGuidanceSource::DpwgGuideline),
            ("citalopram", PrescribingGuidanceSource::CpicGuideline),
            ("citalopram", PrescribingGuidanceSource::DpwgGuideline),
        ] {
            let report = run
                .context
                .drug_report(source, drug)
                .unwrap_or_else(|| panic!("{drug} {source:?} report"));
            assert_eq!(report.matched_annotation_count(), 1, "{drug} {source:?}");
        }
        for source in [
            PrescribingGuidanceSource::CpicGuideline,
            PrescribingGuidanceSource::DpwgGuideline,
            PrescribingGuidanceSource::FdaLabel,
            PrescribingGuidanceSource::FdaAssoc,
        ] {
            assert!(
                run.context
                    .drug_report(source, "ivacaftor")
                    .is_none_or(|report| report.matched_annotation_count() == 0),
                "ivacaftor {source:?}"
            );
        }

        let tsv =
            fs::read_to_string(outputs.reporter_calls_only_tsv.as_ref().unwrap()).expect("tsv");
        assert!(tsv.contains("CYP2C19"));
        assert!(tsv.contains("*1/*1"));
        assert!(tsv.contains("CYP2D6"));
        assert!(tsv.contains("*3/*4"));
    }

    #[test]
    fn run_reporter_from_vcf_cyp2c19_s1s2_het_ambiguity_message_like_java_pipeline_test() {
        let definition =
            read_definition_file(Path::new(CYP2C19_DEFINITION_PATH)).expect("CYP2C19 definition");
        let definitions = DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition.clone())]
                .into_iter()
                .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");
        let messages = MessageCatalog::from_path(Path::new(MESSAGE_PATH)).expect("messages");

        let mut vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n".to_owned();
        append_definition_vcf_rows(
            &mut vcf,
            &definition,
            &[
                ("rs12769205", "G", "0/1"),
                ("rs58973490", "A", "0/1"),
                ("rs4244285", "A", "0/1"),
                ("rs3758581", "G", "1/1"),
            ],
            &[],
        );
        let vcf_file = write_temp_named_file("pipeline-cyp2c19-s1s2-rs58973490-het.vcf", &vcf);
        let outside_file = write_temp_named_file(
            "pipeline-cyp2c19-s1s2-rs58973490-het.outside.tsv",
            "CYP2D6\t*3/*4\nG6PD\tB (wildtype)/B (wildtype)\n",
        );

        let run = run_reporter_from_vcf(
            &vcf_file,
            Some("PharmCAT"),
            &definitions,
            &phenotypes,
            &guidance,
            None,
            &ReporterPipelineOptions {
                include_combinations: false,
                outside_call_files: vec![outside_file],
                message_catalog: Some(messages),
                html: HtmlReportOptions {
                    compact: true,
                    ..HtmlReportOptions::default()
                },
                ..ReporterPipelineOptions::default()
            },
        )
        .expect("CYP2C19 *1/*2 pipeline run with messages");

        let cyp2c19 = run.context.gene_report("CYP2C19").expect("CYP2C19 report");
        assert_eq!(cyp2c19.source_diplotype.as_deref(), Some("*1/*2"));
        assert!(
            cyp2c19
                .recommendation_diplotypes
                .iter()
                .any(|diplotype| diplotype.label == "*1/*2")
        );
        let rs58973490 = cyp2c19
            .variant_reports
            .iter()
            .find(|variant| variant.db_snp_id.as_deref() == Some("rs58973490"))
            .expect("rs58973490 variant report");
        assert!(rs58973490.is_het_call());
        let ambiguity_messages = cyp2c19
            .messages
            .iter()
            .filter(|message| {
                message.exception_type == "ambiguity"
                    && message
                        .matches
                        .variants
                        .iter()
                        .any(|rsid| rsid == "rs58973490")
            })
            .collect::<Vec<_>>();
        assert_eq!(ambiguity_messages.len(), 1);
        assert_eq!(ambiguity_messages[0].name, "CYP2C19 *1/*2 warning");

        let cyp2d6 = run.context.gene_report("CYP2D6").expect("CYP2D6 report");
        assert!(cyp2d6.outside_call);
        assert_eq!(cyp2d6.source_diplotype.as_deref(), Some("*3/*4"));

        assert_eq!(matched_annotation_count(&run.context, "amitriptyline"), 3);
        assert!(
            run.context
                .drug_report(PrescribingGuidanceSource::CpicGuideline, "amitriptyline")
                .is_some_and(|report| report.matched_annotation_count() > 0)
        );
        assert!(
            run.context
                .drug_report(PrescribingGuidanceSource::DpwgGuideline, "amitriptyline")
                .is_some_and(|report| report.matched_annotation_count() > 0)
        );
        assert!(
            run.context
                .drug_report(PrescribingGuidanceSource::FdaAssoc, "amitriptyline")
                .is_some_and(|report| report.matched_annotation_count() > 0)
        );
        assert_eq!(matched_annotation_count(&run.context, "citalopram"), 2);
        assert!(
            run.context
                .drug_report(PrescribingGuidanceSource::CpicGuideline, "citalopram")
                .is_some_and(|report| report.matched_annotation_count() > 0)
        );
        assert!(
            run.context
                .drug_report(PrescribingGuidanceSource::DpwgGuideline, "citalopram")
                .is_some_and(|report| report.matched_annotation_count() > 0)
        );
        assert_eq!(matched_annotation_count(&run.context, "clomipramine"), 4);
        assert_eq!(
            run.context
                .drug_report(PrescribingGuidanceSource::CpicGuideline, "clomipramine")
                .map(|report| report.matched_annotation_count()),
            Some(1)
        );
        assert_eq!(
            run.context
                .drug_report(PrescribingGuidanceSource::DpwgGuideline, "clomipramine")
                .map(|report| report.matched_annotation_count()),
            Some(2)
        );
        assert_eq!(
            run.context
                .drug_report(PrescribingGuidanceSource::FdaAssoc, "clomipramine")
                .map(|report| report.matched_annotation_count()),
            Some(1)
        );
        assert_eq!(matched_annotation_count(&run.context, "ivacaftor"), 0);

        assert_eq!(
            run.context
                .drug_report(PrescribingGuidanceSource::CpicGuideline, "amitriptyline")
                .expect("amitriptyline CPIC report")
                .messages
                .iter()
                .filter(|message| message.name == "CYP2C19 *1/*2 warning")
                .count(),
            1
        );
    }

    #[test]
    fn run_reporter_from_vcf_cyp2c19_s1s2_hom_suppresses_ambiguity_message_like_java_pipeline_test()
    {
        let definition =
            read_definition_file(Path::new(CYP2C19_DEFINITION_PATH)).expect("CYP2C19 definition");
        let definitions = DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition.clone())]
                .into_iter()
                .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");
        let messages = MessageCatalog::from_path(Path::new(MESSAGE_PATH)).expect("messages");

        let mut vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n".to_owned();
        append_definition_vcf_rows(
            &mut vcf,
            &definition,
            &[
                ("rs12769205", "G", "0/1"),
                ("rs4244285", "A", "0/1"),
                ("rs3758581", "G", "1/1"),
            ],
            &[],
        );
        let vcf_file = write_temp_named_file("pipeline-cyp2c19-s1s2-rs58973490-hom.vcf", &vcf);
        let outside_file = write_temp_named_file(
            "pipeline-cyp2c19-s1s2-rs58973490-hom.outside.tsv",
            "CYP2D6\t*3/*4\nG6PD\tB (wildtype)/B (wildtype)\n",
        );

        let run = run_reporter_from_vcf(
            &vcf_file,
            Some("PharmCAT"),
            &definitions,
            &phenotypes,
            &guidance,
            None,
            &ReporterPipelineOptions {
                include_combinations: false,
                outside_call_files: vec![outside_file],
                message_catalog: Some(messages),
                html: HtmlReportOptions {
                    compact: true,
                    ..HtmlReportOptions::default()
                },
                ..ReporterPipelineOptions::default()
            },
        )
        .expect("CYP2C19 *1/*2 pipeline run with homozygous rs58973490");

        let cyp2c19 = run.context.gene_report("CYP2C19").expect("CYP2C19 report");
        assert_eq!(cyp2c19.source_diplotype.as_deref(), Some("*1/*2"));
        assert!(
            cyp2c19
                .recommendation_diplotypes
                .iter()
                .any(|diplotype| diplotype.label == "*1/*2")
        );
        let rs58973490 = cyp2c19
            .variant_reports
            .iter()
            .find(|variant| variant.db_snp_id.as_deref() == Some("rs58973490"))
            .expect("rs58973490 variant report");
        assert!(!rs58973490.is_het_call());
        assert_eq!(
            cyp2c19
                .messages
                .iter()
                .filter(|message| {
                    message.exception_type == "ambiguity"
                        && message
                            .matches
                            .variants
                            .iter()
                            .any(|rsid| rsid == "rs58973490")
                })
                .count(),
            0
        );
        assert!(cyp2c19.messages.is_empty());

        assert_eq!(matched_annotation_count(&run.context, "amitriptyline"), 3);
        assert!(
            run.context
                .drug_report(PrescribingGuidanceSource::CpicGuideline, "amitriptyline")
                .is_some_and(|report| report.matched_annotation_count() > 0)
        );
        assert!(
            run.context
                .drug_report(PrescribingGuidanceSource::DpwgGuideline, "amitriptyline")
                .is_some_and(|report| report.matched_annotation_count() > 0)
        );
        assert!(
            run.context
                .drug_report(PrescribingGuidanceSource::FdaAssoc, "amitriptyline")
                .is_some_and(|report| report.matched_annotation_count() > 0)
        );
        assert_eq!(matched_annotation_count(&run.context, "citalopram"), 2);
        assert!(
            run.context
                .drug_report(PrescribingGuidanceSource::CpicGuideline, "citalopram")
                .is_some_and(|report| report.matched_annotation_count() > 0)
        );
        assert!(
            run.context
                .drug_report(PrescribingGuidanceSource::DpwgGuideline, "citalopram")
                .is_some_and(|report| report.matched_annotation_count() > 0)
        );
        assert_eq!(matched_annotation_count(&run.context, "clomipramine"), 4);
        assert!(
            run.context
                .drug_report(PrescribingGuidanceSource::CpicGuideline, "clomipramine")
                .is_some_and(|report| report.matched_annotation_count() > 0)
        );
        assert!(
            run.context
                .drug_report(PrescribingGuidanceSource::DpwgGuideline, "clomipramine")
                .is_some_and(|report| report.matched_annotation_count() > 0)
        );
        assert!(
            run.context
                .drug_report(PrescribingGuidanceSource::FdaAssoc, "clomipramine")
                .is_some_and(|report| report.matched_annotation_count() > 0)
        );
        assert_eq!(matched_annotation_count(&run.context, "ivacaftor"), 0);

        assert_eq!(
            run.context
                .drug_report(PrescribingGuidanceSource::CpicGuideline, "amitriptyline")
                .expect("amitriptyline CPIC report")
                .messages
                .iter()
                .filter(|message| message.name == "CYP2C19 *1/*2 warning")
                .count(),
            0
        );
    }

    #[test]
    fn run_reporter_from_vcf_cyp2c19_s2s2_clomipramine_counts_like_java_pipeline_test() {
        let definition =
            read_definition_file(Path::new(CYP2C19_DEFINITION_PATH)).expect("CYP2C19 definition");
        let definitions = DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition.clone())]
                .into_iter()
                .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");
        let messages = MessageCatalog::from_path(Path::new(MESSAGE_PATH)).expect("messages");

        let mut vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n".to_owned();
        append_definition_vcf_rows(
            &mut vcf,
            &definition,
            &[
                ("rs12769205", "G", "1/1"),
                ("rs4244285", "A", "1/1"),
                ("rs3758581", "G", "1/1"),
            ],
            &[],
        );
        let vcf_file = write_temp_named_file("pipeline-cyp2c19-s2s2-clomipramine.vcf", &vcf);
        let outside_file = write_temp_named_file(
            "pipeline-cyp2c19-s2s2-clomipramine.outside.tsv",
            "CYP2D6\t*3/*4\nG6PD\tB (wildtype)/B (wildtype)\n",
        );

        let run = run_reporter_from_vcf(
            &vcf_file,
            Some("PharmCAT"),
            &definitions,
            &phenotypes,
            &guidance,
            None,
            &ReporterPipelineOptions {
                include_combinations: false,
                outside_call_files: vec![outside_file],
                message_catalog: Some(messages),
                html: HtmlReportOptions {
                    compact: true,
                    ..HtmlReportOptions::default()
                },
                ..ReporterPipelineOptions::default()
            },
        )
        .expect("CYP2C19 *2/*2 clomipramine pipeline run");

        let cyp2c19 = run.context.gene_report("CYP2C19").expect("CYP2C19 report");
        assert_eq!(cyp2c19.source_diplotype.as_deref(), Some("*2/*2"));
        assert!(
            cyp2c19
                .recommendation_diplotypes
                .iter()
                .any(|diplotype| diplotype.label == "*2/*2")
        );

        assert_eq!(matched_annotation_count(&run.context, "amitriptyline"), 3);
        assert_drug_has_match_from_source(
            &run.context,
            "amitriptyline",
            PrescribingGuidanceSource::CpicGuideline,
        );
        assert_drug_has_match_from_source(
            &run.context,
            "amitriptyline",
            PrescribingGuidanceSource::DpwgGuideline,
        );
        assert_drug_has_match_from_source(
            &run.context,
            "amitriptyline",
            PrescribingGuidanceSource::FdaAssoc,
        );

        assert_eq!(matched_annotation_count(&run.context, "clomipramine"), 4);
        assert_drug_has_match_from_source(
            &run.context,
            "clomipramine",
            PrescribingGuidanceSource::CpicGuideline,
        );
        assert_drug_has_match_from_source(
            &run.context,
            "clomipramine",
            PrescribingGuidanceSource::DpwgGuideline,
        );

        assert_eq!(matched_annotation_count(&run.context, "desipramine"), 2);
        assert_drug_has_match_from_source(
            &run.context,
            "desipramine",
            PrescribingGuidanceSource::CpicGuideline,
        );
        assert_drug_has_match_from_source(
            &run.context,
            "desipramine",
            PrescribingGuidanceSource::FdaAssoc,
        );

        assert_eq!(matched_annotation_count(&run.context, "doxepin"), 4);
        assert_drug_has_match_from_source(
            &run.context,
            "doxepin",
            PrescribingGuidanceSource::CpicGuideline,
        );
        assert_drug_has_match_from_source(
            &run.context,
            "doxepin",
            PrescribingGuidanceSource::DpwgGuideline,
        );
        assert_drug_has_match_from_source(
            &run.context,
            "doxepin",
            PrescribingGuidanceSource::FdaAssoc,
        );

        assert_eq!(matched_annotation_count(&run.context, "imipramine"), 4);
        assert_drug_has_match_from_source(
            &run.context,
            "imipramine",
            PrescribingGuidanceSource::CpicGuideline,
        );
        assert_drug_has_match_from_source(
            &run.context,
            "imipramine",
            PrescribingGuidanceSource::DpwgGuideline,
        );
        assert_drug_has_match_from_source(
            &run.context,
            "imipramine",
            PrescribingGuidanceSource::FdaAssoc,
        );

        assert_eq!(matched_annotation_count(&run.context, "nortriptyline"), 3);
        assert_drug_has_match_from_source(
            &run.context,
            "nortriptyline",
            PrescribingGuidanceSource::CpicGuideline,
        );
        assert_drug_has_match_from_source(
            &run.context,
            "nortriptyline",
            PrescribingGuidanceSource::DpwgGuideline,
        );
        assert_drug_has_match_from_source(
            &run.context,
            "nortriptyline",
            PrescribingGuidanceSource::FdaAssoc,
        );

        assert_eq!(matched_annotation_count(&run.context, "trimipramine"), 2);

        assert_eq!(matched_annotation_count(&run.context, "clopidogrel"), 6);
        assert_drug_has_match_from_source(
            &run.context,
            "clopidogrel",
            PrescribingGuidanceSource::CpicGuideline,
        );
        assert_drug_has_match_from_source(
            &run.context,
            "clopidogrel",
            PrescribingGuidanceSource::DpwgGuideline,
        );
        assert_drug_has_match_from_source(
            &run.context,
            "clopidogrel",
            PrescribingGuidanceSource::FdaLabel,
        );
        assert_drug_has_match_from_source(
            &run.context,
            "clopidogrel",
            PrescribingGuidanceSource::FdaAssoc,
        );

        assert_eq!(matched_annotation_count(&run.context, "lansoprazole"), 3);
        assert_drug_has_match_from_source(
            &run.context,
            "lansoprazole",
            PrescribingGuidanceSource::CpicGuideline,
        );
        assert_drug_has_match_from_source(
            &run.context,
            "lansoprazole",
            PrescribingGuidanceSource::DpwgGuideline,
        );
        assert_drug_has_match_from_source(
            &run.context,
            "lansoprazole",
            PrescribingGuidanceSource::FdaAssoc,
        );

        assert_eq!(
            run.context
                .drug_report(PrescribingGuidanceSource::CpicGuideline, "voriconazole")
                .map(|report| report.matched_annotation_count()),
            Some(2)
        );
        assert_eq!(
            run.context
                .drug_report(PrescribingGuidanceSource::DpwgGuideline, "voriconazole")
                .map(|report| report.matched_annotation_count()),
            Some(1)
        );
    }

    #[test]
    fn run_reporter_from_vcf_cyp2c19_no_call_suppresses_recommendation_matches_like_java_pipeline_test()
     {
        let definition =
            read_definition_file(Path::new(CYP2C19_DEFINITION_PATH)).expect("CYP2C19 definition");
        let definitions = DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition.clone())]
                .into_iter()
                .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");
        let messages = MessageCatalog::from_path(Path::new(MESSAGE_PATH)).expect("messages");

        let mut vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n".to_owned();
        append_definition_vcf_rows(
            &mut vcf,
            &definition,
            &[("rs12769205", "G", "0/1"), ("rs4244285", "A", "1/1")],
            &[],
        );
        let vcf_file = write_temp_named_file("pipeline-cyp2c19-no-call.vcf", &vcf);
        let outside_file = write_temp_named_file(
            "pipeline-cyp2c19-no-call.outside.tsv",
            "CYP2D6\t*3/*4\nG6PD\tB (wildtype)/B (wildtype)\n",
        );

        let run = run_reporter_from_vcf(
            &vcf_file,
            Some("PharmCAT"),
            &definitions,
            &phenotypes,
            &guidance,
            None,
            &ReporterPipelineOptions {
                include_combinations: false,
                outside_call_files: vec![outside_file],
                message_catalog: Some(messages),
                html: HtmlReportOptions {
                    compact: true,
                    ..HtmlReportOptions::default()
                },
                ..ReporterPipelineOptions::default()
            },
        )
        .expect("CYP2C19 no-call pipeline run");

        let result = run
            .gene_call_results
            .iter()
            .find(|result| result.gene == "CYP2C19")
            .expect("CYP2C19 matcher result");
        assert!(
            matches!(result.kind, GeneCallKind::NoCall),
            "expected CYP2C19 no-call, got {:?}",
            result.kind
        );
        let cyp2c19 = run.context.gene_report("CYP2C19").expect("CYP2C19 report");
        assert_eq!(cyp2c19.source_diplotype.as_deref(), Some("Unknown/Unknown"));
        assert!(
            cyp2c19
                .recommendation_diplotypes
                .iter()
                .all(|diplotype| diplotype.label == "Unknown/Unknown")
        );

        let cyp2d6 = run.context.gene_report("CYP2D6").expect("CYP2D6 report");
        assert!(cyp2d6.outside_call);
        assert_eq!(cyp2d6.source_diplotype.as_deref(), Some("*3/*4"));

        assert_drug_has_no_match_from_source(
            &run.context,
            "citalopram",
            PrescribingGuidanceSource::CpicGuideline,
        );
        assert_drug_has_no_match_from_source(
            &run.context,
            "citalopram",
            PrescribingGuidanceSource::DpwgGuideline,
        );
        assert_drug_has_no_match_from_source(
            &run.context,
            "ivacaftor",
            PrescribingGuidanceSource::CpicGuideline,
        );
        assert_drug_has_no_match_from_source(
            &run.context,
            "ivacaftor",
            PrescribingGuidanceSource::DpwgGuideline,
        );
    }

    #[test]
    fn run_reporter_from_vcf_cyp2c19_s4b_s17_rs28399504_missing_like_java_pipeline_test() {
        let definition =
            read_definition_file(Path::new(CYP2C19_DEFINITION_PATH)).expect("CYP2C19 definition");
        let definitions = DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition.clone())]
                .into_iter()
                .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");

        let mut vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n".to_owned();
        append_definition_vcf_rows_omitting(
            &mut vcf,
            &definition,
            &[("rs12248560", "T", "1/1"), ("rs3758581", "G", "1/1")],
            &["rs28399504"],
        );
        let vcf_file =
            write_temp_named_file("pipeline-cyp2c19-s4b-s17-rs28399504-missing.vcf", &vcf);

        let run = run_reporter_from_vcf(
            &vcf_file,
            Some("PharmCAT"),
            &definitions,
            &phenotypes,
            &guidance,
            None,
            &ReporterPipelineOptions {
                include_combinations: false,
                html: HtmlReportOptions {
                    compact: true,
                    ..HtmlReportOptions::default()
                },
                ..ReporterPipelineOptions::default()
            },
        )
        .expect("CYP2C19 *4B/*17 missing rs28399504 pipeline run");

        let result = run
            .gene_call_results
            .iter()
            .find(|result| result.gene == "CYP2C19")
            .expect("CYP2C19 matcher result");
        let GeneCallKind::Diplotypes(diplotypes) = &result.kind else {
            panic!("expected CYP2C19 diplotype call, got {:?}", result.kind);
        };
        assert_eq!(
            diplotypes
                .iter()
                .map(|diplotype| diplotype.name.as_str())
                .collect::<Vec<_>>(),
            ["*4/*4", "*4/*17", "*17/*17"]
        );
        assert_eq!(
            result
                .match_data
                .missing_positions
                .iter()
                .map(|variant| variant.rsid.as_deref())
                .collect::<Vec<_>>(),
            [Some("rs28399504")]
        );

        let cyp2c19 = run.context.gene_report("CYP2C19").expect("CYP2C19 report");
        assert_eq!(cyp2c19.source_diplotype.as_deref(), Some("*4/*4"));
        assert_eq!(
            cyp2c19
                .recommendation_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["*4/*4", "*4/*17", "*17/*17"]
        );

        assert_eq!(matched_annotation_count(&run.context, "citalopram"), 8);
        assert_eq!(
            run.context
                .drug_report(PrescribingGuidanceSource::CpicGuideline, "citalopram")
                .map(|report| report.matched_annotation_count()),
            Some(3)
        );
        assert_eq!(
            run.context
                .drug_report(PrescribingGuidanceSource::DpwgGuideline, "citalopram")
                .map(|report| report.matched_annotation_count()),
            Some(3)
        );
        assert_eq!(
            run.context
                .drug_report(PrescribingGuidanceSource::FdaLabel, "citalopram")
                .map(|report| report.matched_annotation_count()),
            Some(1)
        );
        assert_eq!(
            run.context
                .drug_report(PrescribingGuidanceSource::FdaAssoc, "citalopram")
                .map(|report| report.matched_annotation_count()),
            Some(1)
        );
    }

    #[test]
    fn run_reporter_from_vcf_cyp2c19_s4_s17_with_cyp2d6_outside_call_like_java_pipeline_test() {
        let definition =
            read_definition_file(Path::new(CYP2C19_DEFINITION_PATH)).expect("CYP2C19 definition");
        let definitions = DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition.clone())]
                .into_iter()
                .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");

        let mut vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n".to_owned();
        append_definition_vcf_rows(
            &mut vcf,
            &definition,
            &[
                ("rs12248560", "T", "1/1"),
                ("rs28399504", "G", "0/1"),
                ("rs3758581", "G", "1/1"),
            ],
            &[],
        );
        let vcf_file = write_temp_named_file("pipeline-cyp2c19-s4-s17.vcf", &vcf);
        let outside_file = write_temp_named_file(
            "pipeline-cyp2c19-s4-s17.outside.tsv",
            "##Test Outside Call Data\nCYP2D6\t*1/*4\t\t\t0.6\t0.75\tp: 0.0\t\t\tv1.9-2017_02_09\n",
        );

        let run = run_reporter_from_vcf(
            &vcf_file,
            Some("PharmCAT"),
            &definitions,
            &phenotypes,
            &guidance,
            None,
            &ReporterPipelineOptions {
                include_combinations: false,
                outside_call_files: vec![outside_file],
                html: HtmlReportOptions {
                    compact: true,
                    ..HtmlReportOptions::default()
                },
                ..ReporterPipelineOptions::default()
            },
        )
        .expect("CYP2C19 *4/*17 with CYP2D6 outside-call pipeline run");

        let result = run
            .gene_call_results
            .iter()
            .find(|result| result.gene == "CYP2C19")
            .expect("CYP2C19 matcher result");
        let GeneCallKind::Diplotypes(diplotypes) = &result.kind else {
            panic!("expected CYP2C19 diplotype call, got {:?}", result.kind);
        };
        assert_eq!(
            diplotypes
                .iter()
                .map(|diplotype| diplotype.name.as_str())
                .collect::<Vec<_>>(),
            ["*4/*17"]
        );

        let cyp2c19 = run.context.gene_report("CYP2C19").expect("CYP2C19 report");
        assert_eq!(cyp2c19.source_diplotype.as_deref(), Some("*4/*17"));
        assert_eq!(
            cyp2c19
                .recommendation_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["*4/*17"]
        );

        let cyp2d6 = run.context.gene_report("CYP2D6").expect("CYP2D6 report");
        assert!(cyp2d6.outside_call);
        assert_eq!(cyp2d6.source_diplotype.as_deref(), Some("*1/*4"));
        assert_eq!(
            cyp2d6
                .recommendation_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["*1/*4"]
        );
    }

    #[test]
    fn run_reporter_from_vcf_cyp2c19_s1_s4_missing_s1_with_cyp2d6_outside_call_like_java_pipeline_test()
     {
        let definition =
            read_definition_file(Path::new(CYP2C19_DEFINITION_PATH)).expect("CYP2C19 definition");
        let definitions = DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition.clone())]
                .into_iter()
                .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");
        let messages = MessageCatalog::from_path(Path::new(MESSAGE_PATH)).expect("messages");

        let mut vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n".to_owned();
        append_definition_vcf_rows_omitting(
            &mut vcf,
            &definition,
            &[("rs12248560", "T", "0/1"), ("rs28399504", "G", "0/1")],
            &["rs3758581"],
        );
        let vcf_file = write_temp_named_file("pipeline-cyp2c19-s1-s4-missing-s1.vcf", &vcf);
        let outside_file = write_temp_named_file(
            "pipeline-cyp2c19-s1-s4-missing-s1.outside.tsv",
            "##Test Outside Call Data\nCYP2D6\t*1/*4\t\t\t0.6\t0.75\tp: 0.0\t\t\tv1.9-2017_02_09\n",
        );

        let run = run_reporter_from_vcf(
            &vcf_file,
            Some("PharmCAT"),
            &definitions,
            &phenotypes,
            &guidance,
            None,
            &ReporterPipelineOptions {
                include_combinations: false,
                outside_call_files: vec![outside_file],
                message_catalog: Some(messages),
                html: HtmlReportOptions {
                    compact: true,
                    ..HtmlReportOptions::default()
                },
                ..ReporterPipelineOptions::default()
            },
        )
        .expect("CYP2C19 *1/*4 partial missing with CYP2D6 outside-call pipeline run");

        let result = run
            .gene_call_results
            .iter()
            .find(|result| result.gene == "CYP2C19")
            .expect("CYP2C19 matcher result");
        let GeneCallKind::Diplotypes(diplotypes) = &result.kind else {
            panic!("expected CYP2C19 diplotype call, got {:?}", result.kind);
        };
        assert_eq!(
            diplotypes
                .iter()
                .map(|diplotype| diplotype.name.as_str())
                .collect::<Vec<_>>(),
            ["*1/*4", "*4/*38"]
        );

        let cyp2d6 = run.context.gene_report("CYP2D6").expect("CYP2D6 report");
        assert!(cyp2d6.outside_call);
        assert_eq!(cyp2d6.source_diplotype.as_deref(), Some("*1/*4"));
        assert_eq!(
            cyp2d6
                .recommendation_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["*1/*4"]
        );

        let cyp2c19 = run.context.gene_report("CYP2C19").expect("CYP2C19 report");
        assert_eq!(cyp2c19.source_diplotype.as_deref(), Some("*1/*4"));
        assert_eq!(
            cyp2c19
                .recommendation_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["*1/*4", "*4/*38"]
        );
        assert!(!cyp2c19.phased);
        assert!(
            cyp2c19
                .variant_reports
                .iter()
                .any(|variant| variant.is_missing())
        );

        let rs12248560 = cyp2c19
            .variant_reports
            .iter()
            .find(|variant| variant.db_snp_id.as_deref() == Some("rs12248560"))
            .expect("rs12248560 variant report");
        assert!(rs12248560.is_het_call());

        let rs3758581 = cyp2c19
            .variant_reports
            .iter()
            .find(|variant| variant.db_snp_id.as_deref() == Some("rs3758581"))
            .expect("rs3758581 variant report");
        assert!(rs3758581.is_missing());
        assert!(
            cyp2c19
                .recommendation_diplotypes
                .iter()
                .any(|diplotype| diplotype.label == "*4/*38")
        );
        assert_eq!(
            cyp2c19
                .messages
                .iter()
                .filter(|message| message.exception_type == MessageAnnotation::TYPE_AMBIGUITY)
                .count(),
            2
        );

        assert_drug_has_match_from_source(
            &run.context,
            "amitriptyline",
            PrescribingGuidanceSource::CpicGuideline,
        );
        let amitriptyline = run
            .context
            .drug_report(PrescribingGuidanceSource::CpicGuideline, "amitriptyline")
            .expect("CPIC amitriptyline report");
        assert_eq!(amitriptyline.messages.len(), 2);
    }

    #[test]
    fn run_reporter_from_vcf_cyp2c19_s1_s38_with_cyp2d6_outside_call_like_java_pipeline_test() {
        let definition =
            read_definition_file(Path::new(CYP2C19_DEFINITION_PATH)).expect("CYP2C19 definition");
        let definitions = DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition.clone())]
                .into_iter()
                .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");

        let mut vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n".to_owned();
        append_definition_vcf_rows_omitting(
            &mut vcf,
            &definition,
            &[("rs3758581", "G", "0/1")],
            &["rs56337013"],
        );
        let vcf_file = write_temp_named_file("pipeline-cyp2c19-s1-s38.vcf", &vcf);
        let outside_file = write_temp_named_file(
            "pipeline-cyp2c19-s1-s38.outside.tsv",
            "##Test Outside Call Data\nCYP2D6\t*1/*4\t\t\t0.6\t0.75\tp: 0.0\t\t\tv1.9-2017_02_09\n",
        );

        let run = run_reporter_from_vcf(
            &vcf_file,
            Some("PharmCAT"),
            &definitions,
            &phenotypes,
            &guidance,
            None,
            &ReporterPipelineOptions {
                include_combinations: false,
                outside_call_files: vec![outside_file],
                html: HtmlReportOptions {
                    compact: true,
                    ..HtmlReportOptions::default()
                },
                ..ReporterPipelineOptions::default()
            },
        )
        .expect("CYP2C19 *1/*38 with CYP2D6 outside-call pipeline run");

        let result = run
            .gene_call_results
            .iter()
            .find(|result| result.gene == "CYP2C19")
            .expect("CYP2C19 matcher result");
        let GeneCallKind::Diplotypes(diplotypes) = &result.kind else {
            panic!("expected CYP2C19 diplotype call, got {:?}", result.kind);
        };
        assert_eq!(
            diplotypes
                .iter()
                .map(|diplotype| diplotype.name.as_str())
                .collect::<Vec<_>>(),
            ["*1/*38"]
        );
        assert_eq!(
            result
                .match_data
                .missing_positions
                .iter()
                .map(|variant| variant.rsid.as_deref())
                .collect::<Vec<_>>(),
            [Some("rs56337013")]
        );

        let cyp2d6 = run.context.gene_report("CYP2D6").expect("CYP2D6 report");
        assert!(cyp2d6.outside_call);
        assert_eq!(cyp2d6.source_diplotype.as_deref(), Some("*1/*4"));
        assert_eq!(
            cyp2d6
                .recommendation_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["*1/*4"]
        );

        let cyp2c19 = run.context.gene_report("CYP2C19").expect("CYP2C19 report");
        assert_eq!(cyp2c19.source_diplotype.as_deref(), Some("*1/*38"));
        assert_eq!(
            cyp2c19
                .recommendation_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["*1/*38"]
        );
        assert!(
            cyp2c19
                .variant_reports
                .iter()
                .find(|variant| variant.db_snp_id.as_deref() == Some("rs56337013"))
                .is_some_and(|variant| variant.is_missing())
        );
    }

    #[test]
    fn run_reporter_from_vcf_cyp2c19_multiple_calls_like_java_pipeline_test() {
        let definition =
            read_definition_file(Path::new(CYP2C19_DEFINITION_PATH)).expect("CYP2C19 definition");
        let definitions = DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition.clone())]
                .into_iter()
                .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");

        let mut vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n".to_owned();
        append_definition_vcf_rows_omitting(
            &mut vcf,
            &definition,
            &[("rs12248560", "T", "0/1"), ("rs3758581", "G", "1/1")],
            &["rs28399504"],
        );
        let vcf_file = write_temp_named_file("pipeline-cyp2c19-multiple-calls.vcf", &vcf);

        let run = run_reporter_from_vcf(
            &vcf_file,
            Some("PharmCAT"),
            &definitions,
            &phenotypes,
            &guidance,
            None,
            &ReporterPipelineOptions {
                include_combinations: false,
                html: HtmlReportOptions {
                    compact: true,
                    ..HtmlReportOptions::default()
                },
                ..ReporterPipelineOptions::default()
            },
        )
        .expect("CYP2C19 multiple-calls pipeline run");

        let result = run
            .gene_call_results
            .iter()
            .find(|result| result.gene == "CYP2C19")
            .expect("CYP2C19 matcher result");
        let GeneCallKind::Diplotypes(diplotypes) = &result.kind else {
            panic!("expected CYP2C19 diplotype call, got {:?}", result.kind);
        };
        assert_eq!(
            diplotypes
                .iter()
                .map(|diplotype| diplotype.name.as_str())
                .collect::<Vec<_>>(),
            ["*1/*4", "*1/*17"]
        );
        assert_eq!(
            result
                .match_data
                .missing_positions
                .iter()
                .map(|variant| variant.rsid.as_deref())
                .collect::<Vec<_>>(),
            [Some("rs28399504")]
        );

        let cyp2c19 = run.context.gene_report("CYP2C19").expect("CYP2C19 report");
        assert_eq!(cyp2c19.source_diplotype.as_deref(), Some("*1/*4"));
        assert_eq!(
            cyp2c19
                .recommendation_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["*1/*4", "*1/*17"]
        );
        assert!(
            cyp2c19
                .variant_reports
                .iter()
                .find(|variant| variant.db_snp_id.as_deref() == Some("rs28399504"))
                .is_some_and(|variant| variant.is_missing())
        );
    }

    #[test]
    fn run_reporter_from_vcf_rosuvastatin_with_dpyd_no_data_like_java_pipeline_test() {
        let abcg2_definition =
            read_definition_file(Path::new(ABCG2_DEFINITION_PATH)).expect("ABCG2 definition");
        let slco1b1_definition =
            read_definition_file(Path::new(SLCO1B1_DEFINITION_PATH)).expect("SLCO1B1 definition");
        let dpyd_definition =
            read_definition_file(Path::new(DPYD_DEFINITION_PATH)).expect("DPYD definition");
        let definitions = DefinitionReader::from_definitions(
            [
                (
                    abcg2_definition.gene_symbol.clone(),
                    abcg2_definition.clone(),
                ),
                (
                    slco1b1_definition.gene_symbol.clone(),
                    slco1b1_definition.clone(),
                ),
                (dpyd_definition.gene_symbol.clone(), dpyd_definition),
            ]
            .into_iter()
            .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");

        let mut vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n".to_owned();
        append_definition_vcf_rows(
            &mut vcf,
            &abcg2_definition,
            &[("rs2231142", "T", "0/1")],
            &[],
        );
        append_definition_vcf_rows(
            &mut vcf,
            &slco1b1_definition,
            &[("rs56101265", "C", "0/1")],
            &[],
        );
        let vcf_file = write_temp_named_file("pipeline-rosuvastatin.vcf", &vcf);
        let output_dir = unique_temp_path("pipeline-rosuvastatin-output");
        let outputs = PipelineOutputPlan {
            base_dir: output_dir.clone(),
            basename: "pipeline-rosuvastatin".to_owned(),
            display_name: "pipeline-rosuvastatin.vcf".to_owned(),
            reporter_title: Some("PharmCAT".to_owned()),
            matcher_json: None,
            matcher_html: None,
            matcher_warnings: None,
            phenotyper_json: None,
            reporter_html: Some(output_dir.join("pipeline-rosuvastatin.report.html")),
            reporter_json: None,
            reporter_calls_only_tsv: None,
        };

        let run = run_reporter_from_vcf(
            &vcf_file,
            Some("PharmCAT"),
            &definitions,
            &phenotypes,
            &guidance,
            Some(&outputs),
            &ReporterPipelineOptions {
                include_combinations: false,
                html: HtmlReportOptions {
                    compact: true,
                    ..HtmlReportOptions::default()
                },
                ..ReporterPipelineOptions::default()
            },
        )
        .expect("rosuvastatin pipeline run");

        for gene in ["ABCG2", "SLCO1B1"] {
            let result = run
                .gene_call_results
                .iter()
                .find(|result| result.gene == gene)
                .unwrap_or_else(|| panic!("{gene} matcher result"));
            assert!(
                matches!(result.kind, GeneCallKind::Diplotypes(_)),
                "{gene} should be called by matcher: {:?}",
                result.kind
            );
        }

        let slco1b1 = run.context.gene_report("SLCO1B1").expect("SLCO1B1 report");
        assert_eq!(slco1b1.source_diplotype.as_deref(), Some("*1/*2"));
        assert_eq!(
            slco1b1
                .recommendation_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["*1/*2"]
        );
        assert_eq!(matched_annotation_count(&run.context, "rosuvastatin"), 2);

        let dpyd = run.context.gene_report("DPYD").expect("DPYD report");
        assert!(!dpyd.variant_reports.is_empty());
        assert!(
            dpyd.variant_reports
                .iter()
                .all(|variant| variant.is_missing())
        );

        let html = fs::read_to_string(outputs.reporter_html.as_ref().unwrap()).expect("html");
        assert!(!html.contains("capecitabine"));
        assert!(
            html.contains("<span class=\"gene dpyd\"><span class=\"no-data\">DPYD</span></span>")
        );
    }

    #[test]
    fn run_reporter_from_vcf_slco1b1_hom_wild_simvastatin_like_java_pipeline_test() {
        let definition =
            read_definition_file(Path::new(SLCO1B1_DEFINITION_PATH)).expect("SLCO1B1 definition");
        let definitions = DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition.clone())]
                .into_iter()
                .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");

        let mut vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n".to_owned();
        append_definition_vcf_rows(&mut vcf, &definition, &[], &[]);
        let vcf_file = write_temp_named_file("pipeline-slco1b1-hom-wild.vcf", &vcf);
        let output_dir = unique_temp_path("pipeline-slco1b1-hom-wild-output");
        let outputs = PipelineOutputPlan {
            base_dir: output_dir.clone(),
            basename: "pipeline-slco1b1-hom-wild".to_owned(),
            display_name: "pipeline-slco1b1-hom-wild.vcf".to_owned(),
            reporter_title: Some("PharmCAT".to_owned()),
            matcher_json: None,
            matcher_html: None,
            matcher_warnings: None,
            phenotyper_json: None,
            reporter_html: Some(output_dir.join("pipeline-slco1b1-hom-wild.report.html")),
            reporter_json: None,
            reporter_calls_only_tsv: None,
        };

        let run = run_reporter_from_vcf(
            &vcf_file,
            Some("PharmCAT"),
            &definitions,
            &phenotypes,
            &guidance,
            Some(&outputs),
            &ReporterPipelineOptions {
                include_combinations: false,
                html: HtmlReportOptions {
                    compact: true,
                    ..HtmlReportOptions::default()
                },
                ..ReporterPipelineOptions::default()
            },
        )
        .expect("SLCO1B1 *1/*1 simvastatin pipeline run");

        let result = run
            .gene_call_results
            .iter()
            .find(|result| result.gene == "SLCO1B1")
            .expect("SLCO1B1 matcher result");
        let GeneCallKind::Diplotypes(diplotypes) = &result.kind else {
            panic!("expected SLCO1B1 diplotype call, got {:?}", result.kind);
        };
        assert_eq!(
            diplotypes
                .iter()
                .map(|diplotype| diplotype.name.as_str())
                .collect::<Vec<_>>(),
            ["*1/*1"]
        );

        let slco1b1 = run.context.gene_report("SLCO1B1").expect("SLCO1B1 report");
        assert_eq!(slco1b1.source_diplotype.as_deref(), Some("*1/*1"));
        assert_eq!(
            slco1b1
                .source_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["*1/*1"]
        );
        assert_eq!(
            slco1b1
                .recommendation_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["*1/*1"]
        );
        assert_eq!(matched_annotation_count(&run.context, "simvastatin"), 2);

        let dpwg_simvastatin = run
            .context
            .drug_report(PrescribingGuidanceSource::DpwgGuideline, "simvastatin")
            .expect("DPWG simvastatin report");
        let first_classification = dpwg_simvastatin
            .guidelines
            .iter()
            .flat_map(|guideline| &guideline.annotations)
            .next()
            .map(|annotation| annotation.classification.as_str());
        assert_eq!(first_classification, Some("No recommendation"));

        let html = fs::read_to_string(outputs.reporter_html.as_ref().unwrap()).expect("html");
        assert!(html.contains("SLCO1B1"));
        assert!(html.contains("*1/*1"));
        assert!(html.contains("simvastatin"));
    }

    #[test]
    fn run_reporter_from_vcf_slco1b1_hom_var_simvastatin_like_java_pipeline_test() {
        let definition =
            read_definition_file(Path::new(SLCO1B1_DEFINITION_PATH)).expect("SLCO1B1 definition");
        let definitions = DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition.clone())]
                .into_iter()
                .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");

        let mut vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n".to_owned();
        append_definition_vcf_rows(
            &mut vcf,
            &definition,
            &[("rs2306283", "G", "0/1"), ("rs4149056", "C", "1/1")],
            &[],
        );
        let vcf_file = write_temp_named_file("pipeline-slco1b1-hom-var.vcf", &vcf);
        let output_dir = unique_temp_path("pipeline-slco1b1-hom-var-output");
        let outputs = PipelineOutputPlan {
            base_dir: output_dir.clone(),
            basename: "pipeline-slco1b1-hom-var".to_owned(),
            display_name: "pipeline-slco1b1-hom-var.vcf".to_owned(),
            reporter_title: Some("PharmCAT".to_owned()),
            matcher_json: None,
            matcher_html: None,
            matcher_warnings: None,
            phenotyper_json: None,
            reporter_html: Some(output_dir.join("pipeline-slco1b1-hom-var.report.html")),
            reporter_json: None,
            reporter_calls_only_tsv: None,
        };

        let run = run_reporter_from_vcf(
            &vcf_file,
            Some("PharmCAT"),
            &definitions,
            &phenotypes,
            &guidance,
            Some(&outputs),
            &ReporterPipelineOptions {
                include_combinations: false,
                html: HtmlReportOptions {
                    compact: true,
                    ..HtmlReportOptions::default()
                },
                ..ReporterPipelineOptions::default()
            },
        )
        .expect("SLCO1B1 *5/*15 simvastatin pipeline run");

        let result = run
            .gene_call_results
            .iter()
            .find(|result| result.gene == "SLCO1B1")
            .expect("SLCO1B1 matcher result");
        let GeneCallKind::Diplotypes(diplotypes) = &result.kind else {
            panic!("expected SLCO1B1 diplotype call, got {:?}", result.kind);
        };
        assert_eq!(
            diplotypes
                .iter()
                .map(|diplotype| diplotype.name.as_str())
                .collect::<Vec<_>>(),
            ["*5/*15"]
        );

        let slco1b1 = run.context.gene_report("SLCO1B1").expect("SLCO1B1 report");
        assert_eq!(slco1b1.source_diplotype.as_deref(), Some("*5/*15"));
        assert_eq!(
            slco1b1
                .source_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["*5/*15"]
        );
        assert_eq!(
            slco1b1
                .recommendation_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["*5/*15"]
        );
        assert_eq!(
            run.context
                .drug_report(PrescribingGuidanceSource::CpicGuideline, "simvastatin")
                .map(|report| report.matched_annotation_count()),
            Some(1)
        );
        assert_eq!(
            run.context
                .drug_report(PrescribingGuidanceSource::DpwgGuideline, "simvastatin")
                .map(|report| report.matched_annotation_count()),
            Some(1)
        );

        let html = fs::read_to_string(outputs.reporter_html.as_ref().unwrap()).expect("html");
        assert!(html.contains("SLCO1B1"));
        assert!(html.contains("*5/*15"));
        assert!(html.contains("simvastatin"));
    }

    #[test]
    fn run_reporter_from_vcf_slco1b1_s1_s44_simvastatin_like_java_pipeline_test() {
        let definition =
            read_definition_file(Path::new(SLCO1B1_DEFINITION_PATH)).expect("SLCO1B1 definition");
        let definitions = DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition.clone())]
                .into_iter()
                .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");

        let mut vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n".to_owned();
        append_definition_vcf_rows(
            &mut vcf,
            &definition,
            &[
                ("rs2306283", "G", "0/1"),
                ("rs11045852", "G", "0/1"),
                ("rs74064213", "G", "0/1"),
            ],
            &[],
        );
        let vcf_file = write_temp_named_file("pipeline-slco1b1-s1-s44.vcf", &vcf);
        let output_dir = unique_temp_path("pipeline-slco1b1-s1-s44-output");
        let outputs = PipelineOutputPlan {
            base_dir: output_dir.clone(),
            basename: "pipeline-slco1b1-s1-s44".to_owned(),
            display_name: "pipeline-slco1b1-s1-s44.vcf".to_owned(),
            reporter_title: Some("PharmCAT".to_owned()),
            matcher_json: None,
            matcher_html: None,
            matcher_warnings: None,
            phenotyper_json: None,
            reporter_html: Some(output_dir.join("pipeline-slco1b1-s1-s44.report.html")),
            reporter_json: None,
            reporter_calls_only_tsv: None,
        };

        let run = run_reporter_from_vcf(
            &vcf_file,
            Some("PharmCAT"),
            &definitions,
            &phenotypes,
            &guidance,
            Some(&outputs),
            &ReporterPipelineOptions {
                include_combinations: false,
                html: HtmlReportOptions {
                    compact: true,
                    ..HtmlReportOptions::default()
                },
                ..ReporterPipelineOptions::default()
            },
        )
        .expect("SLCO1B1 *1/*44 simvastatin pipeline run");

        let result = run
            .gene_call_results
            .iter()
            .find(|result| result.gene == "SLCO1B1")
            .expect("SLCO1B1 matcher result");
        let GeneCallKind::Diplotypes(diplotypes) = &result.kind else {
            panic!("expected SLCO1B1 diplotype call, got {:?}", result.kind);
        };
        assert_eq!(
            diplotypes
                .iter()
                .map(|diplotype| diplotype.name.as_str())
                .collect::<Vec<_>>(),
            ["*1/*44"]
        );

        let slco1b1 = run.context.gene_report("SLCO1B1").expect("SLCO1B1 report");
        assert_eq!(slco1b1.source_diplotype.as_deref(), Some("*1/*44"));
        assert_eq!(
            slco1b1
                .source_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["*1/*44"]
        );
        assert_eq!(
            slco1b1
                .recommendation_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["*1/*44"]
        );
        assert_eq!(matched_annotation_count(&run.context, "simvastatin"), 2);

        assert_eq!(
            run.context
                .drug_report(PrescribingGuidanceSource::CpicGuideline, "simvastatin")
                .map(|report| report.matched_annotation_count()),
            Some(1)
        );
        assert_eq!(
            run.context
                .drug_report(PrescribingGuidanceSource::DpwgGuideline, "simvastatin")
                .map(|report| report.matched_annotation_count()),
            Some(1)
        );

        let html = fs::read_to_string(outputs.reporter_html.as_ref().unwrap()).expect("html");
        assert!(html.contains("SLCO1B1"));
        assert!(html.contains("*1/*44"));
        assert!(html.contains("simvastatin"));
    }

    #[test]
    fn run_reporter_from_vcf_slco1b1_s1_s15_simvastatin_like_java_pipeline_test() {
        let definition =
            read_definition_file(Path::new(SLCO1B1_DEFINITION_PATH)).expect("SLCO1B1 definition");
        let definitions = DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition.clone())]
                .into_iter()
                .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");

        let mut vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n".to_owned();
        append_definition_vcf_rows(
            &mut vcf,
            &definition,
            &[("rs2306283", "G", "0/1"), ("rs4149056", "C", "0/1")],
            &[],
        );
        let vcf_file = write_temp_named_file("pipeline-slco1b1-s1-s15.vcf", &vcf);
        let output_dir = unique_temp_path("pipeline-slco1b1-s1-s15-output");
        let outputs = PipelineOutputPlan {
            base_dir: output_dir.clone(),
            basename: "pipeline-slco1b1-s1-s15".to_owned(),
            display_name: "pipeline-slco1b1-s1-s15.vcf".to_owned(),
            reporter_title: Some("PharmCAT".to_owned()),
            matcher_json: None,
            matcher_html: None,
            matcher_warnings: None,
            phenotyper_json: None,
            reporter_html: Some(output_dir.join("pipeline-slco1b1-s1-s15.report.html")),
            reporter_json: None,
            reporter_calls_only_tsv: None,
        };

        let run = run_reporter_from_vcf(
            &vcf_file,
            Some("PharmCAT"),
            &definitions,
            &phenotypes,
            &guidance,
            Some(&outputs),
            &ReporterPipelineOptions {
                include_combinations: false,
                html: HtmlReportOptions {
                    compact: true,
                    ..HtmlReportOptions::default()
                },
                ..ReporterPipelineOptions::default()
            },
        )
        .expect("SLCO1B1 *1/*15 simvastatin pipeline run");

        let result = run
            .gene_call_results
            .iter()
            .find(|result| result.gene == "SLCO1B1")
            .expect("SLCO1B1 matcher result");
        let GeneCallKind::Diplotypes(diplotypes) = &result.kind else {
            panic!("expected SLCO1B1 diplotype call, got {:?}", result.kind);
        };
        assert_eq!(
            diplotypes
                .iter()
                .map(|diplotype| diplotype.name.as_str())
                .collect::<Vec<_>>(),
            ["*1/*15"]
        );

        let slco1b1 = run.context.gene_report("SLCO1B1").expect("SLCO1B1 report");
        assert_eq!(slco1b1.source_diplotype.as_deref(), Some("*1/*15"));
        assert_eq!(
            slco1b1
                .source_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["*1/*15"]
        );
        assert_eq!(
            slco1b1
                .recommendation_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["*1/*15"]
        );
        assert_eq!(
            run.context
                .drug_report(PrescribingGuidanceSource::CpicGuideline, "simvastatin")
                .map(|report| report.matched_annotation_count()),
            Some(1)
        );
        assert_eq!(
            run.context
                .drug_report(PrescribingGuidanceSource::DpwgGuideline, "simvastatin")
                .map(|report| report.matched_annotation_count()),
            Some(1)
        );

        let html = fs::read_to_string(outputs.reporter_html.as_ref().unwrap()).expect("html");
        assert!(html.contains("SLCO1B1"));
        assert!(html.contains("*1/*15"));
        assert!(html.contains("simvastatin"));
    }

    #[test]
    fn run_reporter_from_vcf_slco1b1_s5_s45_simvastatin_intermediates_like_java_pipeline_test() {
        let definition =
            read_definition_file(Path::new(SLCO1B1_DEFINITION_PATH)).expect("SLCO1B1 definition");
        let definitions = DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition.clone())]
                .into_iter()
                .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");

        let mut vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n".to_owned();
        append_definition_vcf_rows(
            &mut vcf,
            &definition,
            &[("rs4149056", "C", "0/1"), ("rs71581941", "T", "0/1")],
            &[],
        );
        let vcf_file = write_temp_named_file("pipeline-slco1b1-s5-s45.vcf", &vcf);
        let output_dir = unique_temp_path("pipeline-slco1b1-s5-s45-output");
        let outputs = PipelineOutputPlan {
            base_dir: output_dir.clone(),
            basename: "pipeline-slco1b1-s5-s45".to_owned(),
            display_name: "pipeline-slco1b1-s5-s45.vcf".to_owned(),
            reporter_title: Some("PharmCAT".to_owned()),
            matcher_json: Some(output_dir.join("pipeline-slco1b1-s5-s45.match.json")),
            matcher_html: None,
            matcher_warnings: Some(output_dir.join("pipeline-slco1b1-s5-s45.match_warnings.txt")),
            phenotyper_json: Some(output_dir.join("pipeline-slco1b1-s5-s45.phenotype.json")),
            reporter_html: Some(output_dir.join("pipeline-slco1b1-s5-s45.report.html")),
            reporter_json: Some(output_dir.join("pipeline-slco1b1-s5-s45.report.json")),
            reporter_calls_only_tsv: None,
        };

        let run = run_reporter_from_vcf(
            &vcf_file,
            Some("PharmCAT"),
            &definitions,
            &phenotypes,
            &guidance,
            Some(&outputs),
            &ReporterPipelineOptions {
                include_combinations: false,
                html: HtmlReportOptions {
                    compact: true,
                    ..HtmlReportOptions::default()
                },
                ..ReporterPipelineOptions::default()
            },
        )
        .expect("SLCO1B1 *5/*45 simvastatin pipeline run");
        write_cli_intermediate_outputs(
            &run,
            &outputs,
            &definitions,
            &vcf_file,
            &PharmcatCliConfig::default(),
        )
        .expect("write intermediates");

        let result = run
            .gene_call_results
            .iter()
            .find(|result| result.gene == "SLCO1B1")
            .expect("SLCO1B1 matcher result");
        let GeneCallKind::Diplotypes(diplotypes) = &result.kind else {
            panic!("expected SLCO1B1 diplotype call, got {:?}", result.kind);
        };
        assert_eq!(
            diplotypes
                .iter()
                .map(|diplotype| diplotype.name.as_str())
                .collect::<Vec<_>>(),
            ["*5/*45"]
        );

        let slco1b1 = run.context.gene_report("SLCO1B1").expect("SLCO1B1 report");
        assert_eq!(slco1b1.source_diplotype.as_deref(), Some("*5/*45"));
        assert_eq!(
            slco1b1
                .source_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["*5/*45"]
        );
        assert_eq!(
            slco1b1
                .recommendation_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["*5/*45"]
        );
        assert_eq!(matched_annotation_count(&run.context, "simvastatin"), 3);
        assert_drug_has_match_from_source(
            &run.context,
            "simvastatin",
            PrescribingGuidanceSource::CpicGuideline,
        );
        assert_drug_has_match_from_source(
            &run.context,
            "simvastatin",
            PrescribingGuidanceSource::DpwgGuideline,
        );
        assert_drug_has_no_match_from_source(
            &run.context,
            "simvastatin",
            PrescribingGuidanceSource::FdaLabel,
        );
        assert_drug_has_match_from_source(
            &run.context,
            "simvastatin",
            PrescribingGuidanceSource::FdaAssoc,
        );

        for path in [
            outputs.matcher_json.as_ref().expect("matcher JSON"),
            outputs.matcher_warnings.as_ref().expect("matcher warnings"),
            outputs.phenotyper_json.as_ref().expect("phenotyper JSON"),
            outputs.reporter_html.as_ref().expect("reporter HTML"),
            outputs.reporter_json.as_ref().expect("reporter JSON"),
        ] {
            assert!(path.is_file(), "missing output {}", path.display());
        }

        let matcher_json =
            fs::read_to_string(outputs.matcher_json.as_ref().unwrap()).expect("matcher JSON");
        let matcher_json: serde_json::Value =
            serde_json::from_str(&matcher_json).expect("parse matcher JSON");
        assert_eq!(matcher_json["results"][0]["gene"].as_str(), Some("SLCO1B1"));
        assert_eq!(
            matcher_json["results"][0]["diplotypes"][0]["name"].as_str(),
            Some("*5/*45")
        );

        let phenotyper_json =
            fs::read_to_string(outputs.phenotyper_json.as_ref().unwrap()).expect("phenotyper JSON");
        assert!(phenotyper_json.contains("*5/*45"));
        let reporter_json =
            fs::read_to_string(outputs.reporter_json.as_ref().unwrap()).expect("reporter JSON");
        assert!(reporter_json.contains("*5/*45"));
        let html = fs::read_to_string(outputs.reporter_html.as_ref().unwrap()).expect("html");
        assert!(html.contains("SLCO1B1"));
        assert!(html.contains("*5/*45"));
        assert!(html.contains("simvastatin"));
    }

    #[test]
    fn run_reporter_from_vcf_slco1b1_s1_s45_warning_intermediates_like_java_pipeline_test() {
        let definition =
            read_definition_file(Path::new(SLCO1B1_DEFINITION_PATH)).expect("SLCO1B1 definition");
        let definitions = DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition.clone())]
                .into_iter()
                .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");
        let messages = MessageCatalog::from_path(Path::new(MESSAGE_PATH)).expect("messages");

        let mut vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n".to_owned();
        append_definition_vcf_rows(
            &mut vcf,
            &definition,
            &[
                ("rs2306283", "G", "0/1"),
                ("rs4149056", "C", "0/1"),
                ("rs71581941", "T", "0/1"),
            ],
            &[],
        );
        let vcf_file = write_temp_named_file("pipeline-slco1b1-s1-s45.vcf", &vcf);
        let output_dir = unique_temp_path("pipeline-slco1b1-s1-s45-output");
        let outputs = PipelineOutputPlan {
            base_dir: output_dir.clone(),
            basename: "pipeline-slco1b1-s1-s45".to_owned(),
            display_name: "pipeline-slco1b1-s1-s45.vcf".to_owned(),
            reporter_title: Some("PharmCAT".to_owned()),
            matcher_json: Some(output_dir.join("pipeline-slco1b1-s1-s45.match.json")),
            matcher_html: None,
            matcher_warnings: Some(output_dir.join("pipeline-slco1b1-s1-s45.match_warnings.txt")),
            phenotyper_json: Some(output_dir.join("pipeline-slco1b1-s1-s45.phenotype.json")),
            reporter_html: Some(output_dir.join("pipeline-slco1b1-s1-s45.report.html")),
            reporter_json: Some(output_dir.join("pipeline-slco1b1-s1-s45.report.json")),
            reporter_calls_only_tsv: None,
        };

        let run = run_reporter_from_vcf(
            &vcf_file,
            Some("PharmCAT"),
            &definitions,
            &phenotypes,
            &guidance,
            Some(&outputs),
            &ReporterPipelineOptions {
                include_combinations: false,
                message_catalog: Some(messages),
                html: HtmlReportOptions {
                    compact: true,
                    ..HtmlReportOptions::default()
                },
                ..ReporterPipelineOptions::default()
            },
        )
        .expect("SLCO1B1 *1/*45 simvastatin pipeline run with messages");
        write_cli_intermediate_outputs(
            &run,
            &outputs,
            &definitions,
            &vcf_file,
            &PharmcatCliConfig::default(),
        )
        .expect("write intermediates");

        let result = run
            .gene_call_results
            .iter()
            .find(|result| result.gene == "SLCO1B1")
            .expect("SLCO1B1 matcher result");
        let GeneCallKind::Diplotypes(diplotypes) = &result.kind else {
            panic!("expected SLCO1B1 diplotype call, got {:?}", result.kind);
        };
        assert_eq!(
            diplotypes
                .iter()
                .map(|diplotype| diplotype.name.as_str())
                .collect::<Vec<_>>(),
            ["*1/*45"]
        );

        let slco1b1 = run.context.gene_report("SLCO1B1").expect("SLCO1B1 report");
        assert_eq!(slco1b1.source_diplotype.as_deref(), Some("*1/*45"));
        assert_eq!(
            slco1b1
                .source_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["*1/*45"]
        );
        assert_eq!(
            slco1b1
                .recommendation_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["*1/*45"]
        );
        assert!(
            slco1b1
                .messages
                .iter()
                .any(|message| message.name == "SLCO1B1 *1/*45 warning"),
            "missing SLCO1B1 *1/*45 warning message"
        );
        assert_eq!(matched_annotation_count(&run.context, "simvastatin"), 3);
        assert_drug_has_match_from_source(
            &run.context,
            "simvastatin",
            PrescribingGuidanceSource::CpicGuideline,
        );
        assert_drug_has_match_from_source(
            &run.context,
            "simvastatin",
            PrescribingGuidanceSource::DpwgGuideline,
        );
        assert_drug_has_no_match_from_source(
            &run.context,
            "simvastatin",
            PrescribingGuidanceSource::FdaLabel,
        );
        assert_drug_has_match_from_source(
            &run.context,
            "simvastatin",
            PrescribingGuidanceSource::FdaAssoc,
        );

        for path in [
            outputs.matcher_json.as_ref().expect("matcher JSON"),
            outputs.matcher_warnings.as_ref().expect("matcher warnings"),
            outputs.phenotyper_json.as_ref().expect("phenotyper JSON"),
            outputs.reporter_html.as_ref().expect("reporter HTML"),
            outputs.reporter_json.as_ref().expect("reporter JSON"),
        ] {
            assert!(path.is_file(), "missing output {}", path.display());
        }

        let matcher_json =
            fs::read_to_string(outputs.matcher_json.as_ref().unwrap()).expect("matcher JSON");
        let matcher_json: serde_json::Value =
            serde_json::from_str(&matcher_json).expect("parse matcher JSON");
        assert_eq!(matcher_json["results"][0]["gene"].as_str(), Some("SLCO1B1"));
        assert_eq!(
            matcher_json["results"][0]["diplotypes"][0]["name"].as_str(),
            Some("*1/*45")
        );

        let phenotyper_json =
            fs::read_to_string(outputs.phenotyper_json.as_ref().unwrap()).expect("phenotyper JSON");
        assert!(phenotyper_json.contains("*1/*45"));
        let reporter_json =
            fs::read_to_string(outputs.reporter_json.as_ref().unwrap()).expect("reporter JSON");
        assert!(reporter_json.contains("*1/*45"));
        assert!(reporter_json.contains("SLCO1B1 *1/*45 warning"));
        let html = fs::read_to_string(outputs.reporter_html.as_ref().unwrap()).expect("html");
        assert!(html.contains("SLCO1B1"));
        assert!(html.contains("*1/*45"));
        assert!(html.contains("simvastatin"));
    }

    #[test]
    fn run_reporter_from_vcf_slco1b1_uncalled_override_like_java_pipeline_test() {
        let definition =
            read_definition_file(Path::new(SLCO1B1_DEFINITION_PATH)).expect("SLCO1B1 definition");
        let definitions = DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition.clone())]
                .into_iter()
                .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");

        let mut vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n".to_owned();
        append_definition_vcf_rows(
            &mut vcf,
            &definition,
            &[
                ("rs2306283", "G", "1/1"),
                ("rs4149056", "C", "0/1"),
                ("rs11045853", "A", "1/1"),
                ("rs72559748", "G", "1/1"),
            ],
            &[],
        );
        let vcf_file = write_temp_named_file("pipeline-slco1b1-uncalled-override.vcf", &vcf);
        let output_dir = unique_temp_path("pipeline-slco1b1-uncalled-override-output");
        let outputs = PipelineOutputPlan {
            base_dir: output_dir.clone(),
            basename: "pipeline-slco1b1-uncalled-override".to_owned(),
            display_name: "pipeline-slco1b1-uncalled-override.vcf".to_owned(),
            reporter_title: Some("PharmCAT".to_owned()),
            matcher_json: None,
            matcher_html: None,
            matcher_warnings: None,
            phenotyper_json: None,
            reporter_html: Some(output_dir.join("pipeline-slco1b1-uncalled-override.report.html")),
            reporter_json: Some(output_dir.join("pipeline-slco1b1-uncalled-override.report.json")),
            reporter_calls_only_tsv: None,
        };

        let run = run_reporter_from_vcf(
            &vcf_file,
            Some("PharmCAT"),
            &definitions,
            &phenotypes,
            &guidance,
            Some(&outputs),
            &ReporterPipelineOptions {
                include_combinations: false,
                html: HtmlReportOptions {
                    compact: true,
                    ..HtmlReportOptions::default()
                },
                ..ReporterPipelineOptions::default()
            },
        )
        .expect("SLCO1B1 uncalled override simvastatin pipeline run");

        let result = run
            .gene_call_results
            .iter()
            .find(|result| result.gene == "SLCO1B1")
            .expect("SLCO1B1 matcher result");
        assert!(
            matches!(result.kind, GeneCallKind::NoCall),
            "expected SLCO1B1 matcher no-call, got {:?}",
            result.kind
        );

        let slco1b1 = run.context.gene_report("SLCO1B1").expect("SLCO1B1 report");
        assert_eq!(slco1b1.source_diplotype.as_deref(), Some("Unknown/Unknown"));
        assert_eq!(
            slco1b1
                .source_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["Unknown/Unknown"]
        );
        assert_eq!(
            slco1b1
                .recommendation_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["*1/*5"]
        );
        let recommendation = &slco1b1.recommendation_diplotypes[0];
        assert!(recommendation.inferred);
        assert_eq!(
            recommendation
                .inferred_source_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["rs4149056 C/rs4149056 T"]
        );
        assert_eq!(slco1b1.lookup_keys, ["Decreased Function"]);
        assert_eq!(matched_annotation_count(&run.context, "simvastatin"), 3);
        assert_drug_has_match_from_source(
            &run.context,
            "simvastatin",
            PrescribingGuidanceSource::CpicGuideline,
        );
        assert_drug_has_match_from_source(
            &run.context,
            "simvastatin",
            PrescribingGuidanceSource::DpwgGuideline,
        );
        assert_drug_has_match_from_source(
            &run.context,
            "simvastatin",
            PrescribingGuidanceSource::FdaAssoc,
        );

        let reporter_json =
            fs::read_to_string(outputs.reporter_json.as_ref().unwrap()).expect("reporter JSON");
        assert!(reporter_json.contains("Unknown/Unknown"));
        assert!(reporter_json.contains("*1/*5"));
        assert!(reporter_json.contains("rs4149056 C/rs4149056 T"));
        let html = fs::read_to_string(outputs.reporter_html.as_ref().unwrap()).expect("html");
        assert!(html.contains("SLCO1B1"));
        assert!(html.contains("rs4149056 C/rs4149056 T"));
        assert!(html.contains("simvastatin"));
    }

    #[test]
    fn run_reporter_from_vcf_ugt1a1_s1_s80_phased_like_java_pipeline_test() {
        let definition =
            read_definition_file(Path::new(UGT1A1_DEFINITION_PATH)).expect("UGT1A1 definition");
        let definitions = DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition.clone())]
                .into_iter()
                .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");

        let mut vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n".to_owned();
        append_definition_vcf_rows_with_default_genotype(
            &mut vcf,
            &definition,
            &[("rs887829", "T", "0|1")],
            &[],
            "0|0",
        );
        let vcf_file = write_temp_named_file("pipeline-ugt1a1-s1-s80-phased.vcf", &vcf);
        let output_dir = unique_temp_path("pipeline-ugt1a1-s1-s80-phased-output");
        let outputs = PipelineOutputPlan {
            base_dir: output_dir.clone(),
            basename: "pipeline-ugt1a1-s1-s80-phased".to_owned(),
            display_name: "pipeline-ugt1a1-s1-s80-phased.vcf".to_owned(),
            reporter_title: Some("PharmCAT".to_owned()),
            matcher_json: None,
            matcher_html: None,
            matcher_warnings: None,
            phenotyper_json: None,
            reporter_html: Some(output_dir.join("pipeline-ugt1a1-s1-s80-phased.report.html")),
            reporter_json: Some(output_dir.join("pipeline-ugt1a1-s1-s80-phased.report.json")),
            reporter_calls_only_tsv: None,
        };

        let run = run_reporter_from_vcf(
            &vcf_file,
            Some("PharmCAT"),
            &definitions,
            &phenotypes,
            &guidance,
            Some(&outputs),
            &ReporterPipelineOptions {
                include_combinations: false,
                html: HtmlReportOptions {
                    compact: true,
                    ..HtmlReportOptions::default()
                },
                ..ReporterPipelineOptions::default()
            },
        )
        .expect("UGT1A1 *1/*80 phased pipeline run");

        let result = run
            .gene_call_results
            .iter()
            .find(|result| result.gene == "UGT1A1")
            .expect("UGT1A1 matcher result");
        assert!(result.match_data.phased);
        let GeneCallKind::Diplotypes(diplotypes) = &result.kind else {
            panic!("expected UGT1A1 diplotype call, got {:?}", result.kind);
        };
        assert_eq!(
            diplotypes
                .iter()
                .map(|diplotype| diplotype.name.as_str())
                .collect::<Vec<_>>(),
            ["*1/*80"]
        );

        let ugt1a1 = run.context.gene_report("UGT1A1").expect("UGT1A1 report");
        assert!(ugt1a1.phased);
        assert_eq!(ugt1a1.source_diplotype.as_deref(), Some("*1/*80"));
        assert_eq!(
            ugt1a1
                .source_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["*1/*80"]
        );
        assert_eq!(
            ugt1a1
                .recommendation_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["*1/*80"]
        );
        let recommendation = &ugt1a1.recommendation_diplotypes[0];
        assert_eq!(
            recommendation
                .allele1
                .as_ref()
                .map(|haplotype| haplotype.name.as_str()),
            Some("*1")
        );
        assert_eq!(
            recommendation
                .allele2
                .as_ref()
                .map(|haplotype| haplotype.name.as_str()),
            Some("*80")
        );
        assert_eq!(ugt1a1.lookup_keys, ["Indeterminate"]);

        let reporter_json =
            fs::read_to_string(outputs.reporter_json.as_ref().unwrap()).expect("reporter JSON");
        assert!(reporter_json.contains("*1/*80"));
        let html = fs::read_to_string(outputs.reporter_html.as_ref().unwrap()).expect("html");
        assert!(html.contains("UGT1A1"));
        assert!(html.contains("*1/*80"));
    }

    #[test]
    fn run_reporter_from_vcf_ugt1a1_s1_s80_unphased_like_java_pipeline_test() {
        let definition =
            read_definition_file(Path::new(UGT1A1_DEFINITION_PATH)).expect("UGT1A1 definition");
        let definitions = DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition.clone())]
                .into_iter()
                .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");

        let mut vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n".to_owned();
        append_definition_vcf_rows(&mut vcf, &definition, &[("rs887829", "T", "0/1")], &[]);
        let vcf_file = write_temp_named_file("pipeline-ugt1a1-s1-s80-unphased.vcf", &vcf);
        let output_dir = unique_temp_path("pipeline-ugt1a1-s1-s80-unphased-output");
        let outputs = PipelineOutputPlan {
            base_dir: output_dir.clone(),
            basename: "pipeline-ugt1a1-s1-s80-unphased".to_owned(),
            display_name: "pipeline-ugt1a1-s1-s80-unphased.vcf".to_owned(),
            reporter_title: Some("PharmCAT".to_owned()),
            matcher_json: None,
            matcher_html: None,
            matcher_warnings: None,
            phenotyper_json: None,
            reporter_html: Some(output_dir.join("pipeline-ugt1a1-s1-s80-unphased.report.html")),
            reporter_json: Some(output_dir.join("pipeline-ugt1a1-s1-s80-unphased.report.json")),
            reporter_calls_only_tsv: None,
        };

        let run = run_reporter_from_vcf(
            &vcf_file,
            Some("PharmCAT"),
            &definitions,
            &phenotypes,
            &guidance,
            Some(&outputs),
            &ReporterPipelineOptions {
                include_combinations: false,
                html: HtmlReportOptions {
                    compact: true,
                    ..HtmlReportOptions::default()
                },
                ..ReporterPipelineOptions::default()
            },
        )
        .expect("UGT1A1 *1/*80 unphased pipeline run");

        let result = run
            .gene_call_results
            .iter()
            .find(|result| result.gene == "UGT1A1")
            .expect("UGT1A1 matcher result");
        assert!(!result.match_data.phased);
        let GeneCallKind::Diplotypes(diplotypes) = &result.kind else {
            panic!("expected UGT1A1 diplotype call, got {:?}", result.kind);
        };
        assert_eq!(
            diplotypes
                .iter()
                .map(|diplotype| diplotype.name.as_str())
                .collect::<Vec<_>>(),
            ["*1/*80"]
        );

        let ugt1a1 = run.context.gene_report("UGT1A1").expect("UGT1A1 report");
        assert!(!ugt1a1.phased);
        assert_eq!(ugt1a1.source_diplotype.as_deref(), Some("*1/*80"));
        assert_eq!(
            ugt1a1
                .source_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["*1/*80"]
        );
        assert_eq!(
            ugt1a1
                .recommendation_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["*1/*80"]
        );
        let recommendation = &ugt1a1.recommendation_diplotypes[0];
        assert_eq!(
            recommendation
                .allele1
                .as_ref()
                .map(|haplotype| haplotype.name.as_str()),
            Some("*1")
        );
        assert_eq!(
            recommendation
                .allele2
                .as_ref()
                .map(|haplotype| haplotype.name.as_str()),
            Some("*80")
        );
        assert_eq!(ugt1a1.lookup_keys, ["Indeterminate"]);

        let reporter_json =
            fs::read_to_string(outputs.reporter_json.as_ref().unwrap()).expect("reporter JSON");
        assert!(reporter_json.contains("*1/*80"));
        let html = fs::read_to_string(outputs.reporter_html.as_ref().unwrap()).expect("html");
        assert!(html.contains("UGT1A1"));
        assert!(html.contains("*1/*80"));
    }

    #[test]
    fn run_reporter_from_vcf_ugt1a1_s1_s1_reference_like_java_pipeline_test() {
        let definition =
            read_definition_file(Path::new(UGT1A1_DEFINITION_PATH)).expect("UGT1A1 definition");
        let definitions = DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition.clone())]
                .into_iter()
                .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");

        let mut vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n".to_owned();
        append_definition_vcf_rows(&mut vcf, &definition, &[], &[]);
        let vcf_file = write_temp_named_file("pipeline-ugt1a1-s1-s1-reference.vcf", &vcf);
        let output_dir = unique_temp_path("pipeline-ugt1a1-s1-s1-reference-output");
        let outputs = PipelineOutputPlan {
            base_dir: output_dir.clone(),
            basename: "pipeline-ugt1a1-s1-s1-reference".to_owned(),
            display_name: "pipeline-ugt1a1-s1-s1-reference.vcf".to_owned(),
            reporter_title: Some("PharmCAT".to_owned()),
            matcher_json: None,
            matcher_html: None,
            matcher_warnings: None,
            phenotyper_json: None,
            reporter_html: Some(output_dir.join("pipeline-ugt1a1-s1-s1-reference.report.html")),
            reporter_json: Some(output_dir.join("pipeline-ugt1a1-s1-s1-reference.report.json")),
            reporter_calls_only_tsv: None,
        };

        let run = run_reporter_from_vcf(
            &vcf_file,
            Some("PharmCAT"),
            &definitions,
            &phenotypes,
            &guidance,
            Some(&outputs),
            &ReporterPipelineOptions {
                include_combinations: false,
                html: HtmlReportOptions {
                    compact: true,
                    ..HtmlReportOptions::default()
                },
                ..ReporterPipelineOptions::default()
            },
        )
        .expect("UGT1A1 *1/*1 reference pipeline run");

        let result = run
            .gene_call_results
            .iter()
            .find(|result| result.gene == "UGT1A1")
            .expect("UGT1A1 matcher result");
        let GeneCallKind::Diplotypes(diplotypes) = &result.kind else {
            panic!("expected UGT1A1 diplotype call, got {:?}", result.kind);
        };
        assert_eq!(
            diplotypes
                .iter()
                .map(|diplotype| diplotype.name.as_str())
                .collect::<Vec<_>>(),
            ["*1/*1"]
        );

        let ugt1a1 = run.context.gene_report("UGT1A1").expect("UGT1A1 report");
        assert_eq!(ugt1a1.source_diplotype.as_deref(), Some("*1/*1"));
        assert_eq!(
            ugt1a1
                .source_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["*1/*1"]
        );
        assert_eq!(
            ugt1a1
                .recommendation_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["*1/*1"]
        );
        let recommendation = &ugt1a1.recommendation_diplotypes[0];
        assert_eq!(
            recommendation
                .allele1
                .as_ref()
                .map(|haplotype| haplotype.name.as_str()),
            Some("*1")
        );
        assert_eq!(
            recommendation
                .allele2
                .as_ref()
                .map(|haplotype| haplotype.name.as_str()),
            Some("*1")
        );
        assert_eq!(ugt1a1.lookup_keys, ["Normal Metabolizer"]);

        let reporter_json =
            fs::read_to_string(outputs.reporter_json.as_ref().unwrap()).expect("reporter JSON");
        assert!(reporter_json.contains("*1/*1"));
        let html = fs::read_to_string(outputs.reporter_html.as_ref().unwrap()).expect("html");
        assert!(html.contains("UGT1A1"));
        assert!(html.contains("*1/*1"));
    }

    #[test]
    fn run_reporter_from_vcf_ugt1a1_s1_s80_s28_like_java_pipeline_test() {
        let definition =
            read_definition_file(Path::new(UGT1A1_DEFINITION_PATH)).expect("UGT1A1 definition");
        let definitions = DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition.clone())]
                .into_iter()
                .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");

        let mut vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n".to_owned();
        append_definition_vcf_rows(
            &mut vcf,
            &definition,
            &[("rs887829", "T", "0/1"), ("rs3064744", "CATAT", "0/1")],
            &[],
        );
        let vcf_file = write_temp_named_file("pipeline-ugt1a1-s1-s80-s28.vcf", &vcf);
        let output_dir = unique_temp_path("pipeline-ugt1a1-s1-s80-s28-output");
        let outputs = PipelineOutputPlan {
            base_dir: output_dir.clone(),
            basename: "pipeline-ugt1a1-s1-s80-s28".to_owned(),
            display_name: "pipeline-ugt1a1-s1-s80-s28.vcf".to_owned(),
            reporter_title: Some("PharmCAT".to_owned()),
            matcher_json: None,
            matcher_html: None,
            matcher_warnings: None,
            phenotyper_json: None,
            reporter_html: Some(output_dir.join("pipeline-ugt1a1-s1-s80-s28.report.html")),
            reporter_json: Some(output_dir.join("pipeline-ugt1a1-s1-s80-s28.report.json")),
            reporter_calls_only_tsv: None,
        };

        let run = run_reporter_from_vcf(
            &vcf_file,
            Some("PharmCAT"),
            &definitions,
            &phenotypes,
            &guidance,
            Some(&outputs),
            &ReporterPipelineOptions {
                include_combinations: false,
                html: HtmlReportOptions {
                    compact: true,
                    ..HtmlReportOptions::default()
                },
                ..ReporterPipelineOptions::default()
            },
        )
        .expect("UGT1A1 *1/*80+*28 pipeline run");

        let result = run
            .gene_call_results
            .iter()
            .find(|result| result.gene == "UGT1A1")
            .expect("UGT1A1 matcher result");
        let GeneCallKind::Diplotypes(diplotypes) = &result.kind else {
            panic!("expected UGT1A1 diplotype call, got {:?}", result.kind);
        };
        assert_eq!(
            diplotypes
                .iter()
                .map(|diplotype| diplotype.name.as_str())
                .collect::<Vec<_>>(),
            ["*1/*80+*28"]
        );

        let ugt1a1 = run.context.gene_report("UGT1A1").expect("UGT1A1 report");
        assert_eq!(ugt1a1.source_diplotype.as_deref(), Some("*1/*80+*28"));
        assert_eq!(
            ugt1a1
                .source_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["*1/*80+*28"]
        );
        assert_eq!(
            ugt1a1
                .recommendation_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["*1/*80+*28"]
        );
        let recommendation = &ugt1a1.recommendation_diplotypes[0];
        assert_eq!(
            recommendation
                .allele1
                .as_ref()
                .map(|haplotype| haplotype.name.as_str()),
            Some("*1")
        );
        assert_eq!(
            recommendation
                .allele2
                .as_ref()
                .map(|haplotype| haplotype.name.as_str()),
            Some("*80+*28")
        );

        let reporter_json =
            fs::read_to_string(outputs.reporter_json.as_ref().unwrap()).expect("reporter JSON");
        assert!(reporter_json.contains("*1/*80+*28"));
        let html = fs::read_to_string(outputs.reporter_html.as_ref().unwrap()).expect("html");
        assert!(html.contains("UGT1A1"));
        assert!(html.contains("*1/*80+*28"));
    }

    #[test]
    fn run_reporter_from_vcf_ugt1a1_s28_s37_like_java_pipeline_test() {
        let definition =
            read_definition_file(Path::new(UGT1A1_DEFINITION_PATH)).expect("UGT1A1 definition");
        let definitions = DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition.clone())]
                .into_iter()
                .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");

        let mut vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n".to_owned();
        append_definition_vcf_rows(
            &mut vcf,
            &definition,
            &[("rs3064744", "CATAT,CATATAT", "1/2")],
            &[],
        );
        let vcf_file = write_temp_named_file("pipeline-ugt1a1-s28-s37.vcf", &vcf);
        let output_dir = unique_temp_path("pipeline-ugt1a1-s28-s37-output");
        let outputs = PipelineOutputPlan {
            base_dir: output_dir.clone(),
            basename: "pipeline-ugt1a1-s28-s37".to_owned(),
            display_name: "pipeline-ugt1a1-s28-s37.vcf".to_owned(),
            reporter_title: Some("PharmCAT".to_owned()),
            matcher_json: None,
            matcher_html: None,
            matcher_warnings: None,
            phenotyper_json: None,
            reporter_html: Some(output_dir.join("pipeline-ugt1a1-s28-s37.report.html")),
            reporter_json: Some(output_dir.join("pipeline-ugt1a1-s28-s37.report.json")),
            reporter_calls_only_tsv: None,
        };

        let run = run_reporter_from_vcf(
            &vcf_file,
            Some("PharmCAT"),
            &definitions,
            &phenotypes,
            &guidance,
            Some(&outputs),
            &ReporterPipelineOptions {
                include_combinations: false,
                html: HtmlReportOptions {
                    compact: true,
                    ..HtmlReportOptions::default()
                },
                ..ReporterPipelineOptions::default()
            },
        )
        .expect("UGT1A1 *28/*37 pipeline run");

        let result = run
            .gene_call_results
            .iter()
            .find(|result| result.gene == "UGT1A1")
            .expect("UGT1A1 matcher result");
        let GeneCallKind::Diplotypes(diplotypes) = &result.kind else {
            panic!("expected UGT1A1 diplotype call, got {:?}", result.kind);
        };
        assert_eq!(
            diplotypes
                .iter()
                .map(|diplotype| diplotype.name.as_str())
                .collect::<Vec<_>>(),
            ["*28/*37"]
        );

        let ugt1a1 = run.context.gene_report("UGT1A1").expect("UGT1A1 report");
        assert_eq!(ugt1a1.source_diplotype.as_deref(), Some("*28/*37"));
        assert_eq!(
            ugt1a1
                .source_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["*28/*37"]
        );
        // Java testRecommendedDiplotypes compares an unordered count map, so *37/*28 is a set.
        let recommendation = &ugt1a1.recommendation_diplotypes[0];
        let mut alleles = [
            recommendation.allele1.as_ref().map(|h| h.name.as_str()),
            recommendation.allele2.as_ref().map(|h| h.name.as_str()),
        ];
        alleles.sort();
        assert_eq!(alleles, [Some("*28"), Some("*37")]);

        let reporter_json =
            fs::read_to_string(outputs.reporter_json.as_ref().unwrap()).expect("reporter JSON");
        assert!(reporter_json.contains("*28/*37"));
        let html = fs::read_to_string(outputs.reporter_html.as_ref().unwrap()).expect("html");
        assert!(html.contains("UGT1A1"));
        assert!(html.contains("*28/*37"));
    }

    #[test]
    fn run_reporter_from_vcf_ugt1a1_s28_s80_phased_like_java_pipeline_test() {
        let definition =
            read_definition_file(Path::new(UGT1A1_DEFINITION_PATH)).expect("UGT1A1 definition");
        let definitions = DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition.clone())]
                .into_iter()
                .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");
        let messages = MessageCatalog::from_path(Path::new(MESSAGE_PATH)).expect("messages");

        let mut vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n".to_owned();
        append_definition_vcf_rows_with_default_genotype(
            &mut vcf,
            &definition,
            &[("rs887829", "T", "0|1"), ("rs3064744", "CATAT", "0|1")],
            &[],
            "0|0",
        );
        let vcf_file = write_temp_named_file("pipeline-ugt1a1-s28-s80-phased.vcf", &vcf);

        let run = run_reporter_from_vcf(
            &vcf_file,
            Some("PharmCAT"),
            &definitions,
            &phenotypes,
            &guidance,
            None,
            &ReporterPipelineOptions {
                include_combinations: false,
                message_catalog: Some(messages),
                html: HtmlReportOptions {
                    compact: true,
                    ..HtmlReportOptions::default()
                },
                ..ReporterPipelineOptions::default()
            },
        )
        .expect("UGT1A1 *1/*80+*28 phased pipeline run");

        let result = run
            .gene_call_results
            .iter()
            .find(|result| result.gene == "UGT1A1")
            .expect("UGT1A1 matcher result");
        assert!(result.match_data.phased);
        let GeneCallKind::Diplotypes(diplotypes) = &result.kind else {
            panic!("expected UGT1A1 diplotype call, got {:?}", result.kind);
        };
        assert_eq!(
            diplotypes
                .iter()
                .map(|diplotype| diplotype.name.as_str())
                .collect::<Vec<_>>(),
            ["*1/*80+*28"]
        );

        let ugt1a1 = run.context.gene_report("UGT1A1").expect("UGT1A1 report");
        assert_eq!(ugt1a1.source_diplotype.as_deref(), Some("*1/*80+*28"));
        let recommendation = &ugt1a1.recommendation_diplotypes[0];
        let mut alleles = [
            recommendation.allele1.as_ref().map(|h| h.name.as_str()),
            recommendation.allele2.as_ref().map(|h| h.name.as_str()),
        ];
        alleles.sort();
        assert_eq!(alleles, [Some("*1"), Some("*80+*28")]);

        // Phased data suppresses the *1/*80+*28 ambiguity message on atazanavir.
        assert_eq!(matched_annotation_count(&run.context, "atazanavir"), 2);
        assert_eq!(
            run.context
                .drug_report(PrescribingGuidanceSource::CpicGuideline, "atazanavir")
                .map(|report| report.matched_annotation_count()),
            Some(1)
        );
        assert_eq!(
            run.context
                .drug_report(PrescribingGuidanceSource::DpwgGuideline, "atazanavir")
                .map(|report| report.matched_annotation_count()),
            Some(1)
        );
        assert_eq!(
            run.context
                .drug_report(PrescribingGuidanceSource::CpicGuideline, "atazanavir")
                .expect("atazanavir CPIC report")
                .messages
                .len(),
            0
        );

        assert_eq!(ugt1a1.messages.len(), 2);
        assert!(
            ugt1a1
                .messages
                .iter()
                .any(|message| message.name == "reference-allele")
        );
    }

    #[test]
    fn run_reporter_from_vcf_ugt1a1_s6_s80_s28_phased_like_java_pipeline_test() {
        let definition =
            read_definition_file(Path::new(UGT1A1_DEFINITION_PATH)).expect("UGT1A1 definition");
        let definitions = DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition.clone())]
                .into_iter()
                .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");

        let mut vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n".to_owned();
        // Java phases allele1/allele2 left/right: hap1 = *80+*28 (T, CATAT), hap2 = *6 (rs4148323 A).
        append_definition_vcf_rows_with_default_genotype(
            &mut vcf,
            &definition,
            &[
                ("rs887829", "T", "1|0"),
                ("rs3064744", "CATAT", "1|0"),
                ("rs4148323", "A", "0|1"),
            ],
            &[],
            "0|0",
        );
        let vcf_file = write_temp_named_file("pipeline-ugt1a1-s6-s80-s28-phased.vcf", &vcf);

        let run = run_reporter_from_vcf(
            &vcf_file,
            Some("PharmCAT"),
            &definitions,
            &phenotypes,
            &guidance,
            None,
            &ReporterPipelineOptions {
                include_combinations: false,
                html: HtmlReportOptions {
                    compact: true,
                    ..HtmlReportOptions::default()
                },
                ..ReporterPipelineOptions::default()
            },
        )
        .expect("UGT1A1 *6/*80+*28 phased pipeline run");

        let result = run
            .gene_call_results
            .iter()
            .find(|result| result.gene == "UGT1A1")
            .expect("UGT1A1 matcher result");
        assert!(result.match_data.phased);
        let GeneCallKind::Diplotypes(diplotypes) = &result.kind else {
            panic!("expected UGT1A1 diplotype call, got {:?}", result.kind);
        };
        assert_eq!(
            diplotypes
                .iter()
                .map(|diplotype| diplotype.name.as_str())
                .collect::<Vec<_>>(),
            ["*6/*80+*28"]
        );

        let ugt1a1 = run.context.gene_report("UGT1A1").expect("UGT1A1 report");
        assert_eq!(ugt1a1.source_diplotype.as_deref(), Some("*6/*80+*28"));
        let recommendation = &ugt1a1.recommendation_diplotypes[0];
        let mut alleles = [
            recommendation.allele1.as_ref().map(|h| h.name.as_str()),
            recommendation.allele2.as_ref().map(|h| h.name.as_str()),
        ];
        alleles.sort();
        assert_eq!(alleles, [Some("*6"), Some("*80+*28")]);
    }

    #[test]
    fn run_reporter_from_vcf_ugt1a1_s6_s6_like_java_pipeline_test() {
        let definition =
            read_definition_file(Path::new(UGT1A1_DEFINITION_PATH)).expect("UGT1A1 definition");
        let definitions = DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition.clone())]
                .into_iter()
                .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");

        let mut vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n".to_owned();
        append_definition_vcf_rows(&mut vcf, &definition, &[("rs4148323", "A", "1/1")], &[]);
        let vcf_file = write_temp_named_file("pipeline-ugt1a1-s6-s6.vcf", &vcf);

        let run = run_reporter_from_vcf(
            &vcf_file,
            Some("PharmCAT"),
            &definitions,
            &phenotypes,
            &guidance,
            None,
            &ReporterPipelineOptions {
                include_combinations: false,
                html: HtmlReportOptions {
                    compact: true,
                    ..HtmlReportOptions::default()
                },
                ..ReporterPipelineOptions::default()
            },
        )
        .expect("UGT1A1 *6/*6 pipeline run");

        let result = run
            .gene_call_results
            .iter()
            .find(|result| result.gene == "UGT1A1")
            .expect("UGT1A1 matcher result");
        let GeneCallKind::Diplotypes(diplotypes) = &result.kind else {
            panic!("expected UGT1A1 diplotype call, got {:?}", result.kind);
        };
        assert_eq!(
            diplotypes
                .iter()
                .map(|diplotype| diplotype.name.as_str())
                .collect::<Vec<_>>(),
            ["*6/*6"]
        );

        let ugt1a1 = run.context.gene_report("UGT1A1").expect("UGT1A1 report");
        assert_eq!(ugt1a1.source_diplotype.as_deref(), Some("*6/*6"));
        let recommendation = &ugt1a1.recommendation_diplotypes[0];
        assert_eq!(
            recommendation
                .allele1
                .as_ref()
                .map(|haplotype| haplotype.name.as_str()),
            Some("*6")
        );
        assert_eq!(
            recommendation
                .allele2
                .as_ref()
                .map(|haplotype| haplotype.name.as_str()),
            Some("*6")
        );
    }

    #[test]
    fn run_reporter_from_vcf_ugt1a1_s6_s80_s28_missing_phased_like_java_pipeline_test() {
        let definition =
            read_definition_file(Path::new(UGT1A1_DEFINITION_PATH)).expect("UGT1A1 definition");
        let definitions = DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition.clone())]
                .into_iter()
                .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");

        let mut vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n".to_owned();
        // Phased with rs3064744 missing: hap1 = *6 (rs4148323 A), hap2 = *80 (rs887829 T) but the
        // missing repeat leaves *80, *80+*28, and *80+*37 indistinguishable.
        append_definition_vcf_rows_with_default_genotype(
            &mut vcf,
            &definition,
            &[("rs887829", "T", "0|1"), ("rs4148323", "A", "1|0")],
            &["rs3064744"],
            "0|0",
        );
        let vcf_file = write_temp_named_file("pipeline-ugt1a1-s6-s80-s28-missing-phased.vcf", &vcf);

        let run = run_reporter_from_vcf(
            &vcf_file,
            Some("PharmCAT"),
            &definitions,
            &phenotypes,
            &guidance,
            None,
            &ReporterPipelineOptions {
                include_combinations: false,
                html: HtmlReportOptions {
                    compact: true,
                    ..HtmlReportOptions::default()
                },
                ..ReporterPipelineOptions::default()
            },
        )
        .expect("UGT1A1 *6/*80 missing phased pipeline run");

        let result = run
            .gene_call_results
            .iter()
            .find(|result| result.gene == "UGT1A1")
            .expect("UGT1A1 matcher result");
        assert!(result.match_data.phased);
        let GeneCallKind::Diplotypes(diplotypes) = &result.kind else {
            panic!("expected UGT1A1 diplotype call, got {:?}", result.kind);
        };
        assert_eq!(
            diplotypes
                .iter()
                .map(|diplotype| diplotype.name.as_str())
                .collect::<Vec<_>>(),
            ["*6/*80", "*6/*80+*28", "*6/*80+*37"]
        );

        let ugt1a1 = run.context.gene_report("UGT1A1").expect("UGT1A1 report");
        assert_eq!(
            ugt1a1
                .source_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["*6/*80", "*6/*80+*28", "*6/*80+*37"]
        );
        let recommendation_labels = ugt1a1
            .recommendation_diplotypes
            .iter()
            .map(|diplotype| diplotype.label.as_str())
            .collect::<Vec<_>>();
        for expected in ["*6/*80", "*6/*80+*28", "*6/*80+*37"] {
            assert!(
                recommendation_labels.contains(&expected),
                "missing recommendation {expected} in {recommendation_labels:?}"
            );
        }
    }

    #[test]
    fn run_reporter_from_vcf_ugt1a1_s80_s28_missing_grouped_recommendations_like_java_pipeline_test()
     {
        let definition =
            read_definition_file(Path::new(UGT1A1_DEFINITION_PATH)).expect("UGT1A1 definition");
        let definitions = DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition.clone())]
                .into_iter()
                .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");

        let mut vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n".to_owned();
        append_definition_vcf_rows(
            &mut vcf,
            &definition,
            &[("rs887829", "T", "0/1")],
            &["rs3064744"],
        );
        let vcf_file = write_temp_named_file("pipeline-ugt1a1-s80-s28-missing.vcf", &vcf);
        let output_dir = unique_temp_path("pipeline-ugt1a1-s80-s28-missing-output");
        let outputs = PipelineOutputPlan {
            base_dir: output_dir.clone(),
            basename: "pipeline-ugt1a1-s80-s28-missing".to_owned(),
            display_name: "pipeline-ugt1a1-s80-s28-missing.vcf".to_owned(),
            reporter_title: Some("PharmCAT".to_owned()),
            matcher_json: None,
            matcher_html: None,
            matcher_warnings: None,
            phenotyper_json: None,
            reporter_html: Some(output_dir.join("pipeline-ugt1a1-s80-s28-missing.report.html")),
            reporter_json: Some(output_dir.join("pipeline-ugt1a1-s80-s28-missing.report.json")),
            reporter_calls_only_tsv: None,
        };

        let run = run_reporter_from_vcf(
            &vcf_file,
            Some("PharmCAT"),
            &definitions,
            &phenotypes,
            &guidance,
            Some(&outputs),
            &ReporterPipelineOptions {
                include_combinations: false,
                html: HtmlReportOptions {
                    compact: true,
                    ..HtmlReportOptions::default()
                },
                ..ReporterPipelineOptions::default()
            },
        )
        .expect("UGT1A1 *1/*80 missing pipeline run");

        let result = run
            .gene_call_results
            .iter()
            .find(|result| result.gene == "UGT1A1")
            .expect("UGT1A1 matcher result");
        let GeneCallKind::Diplotypes(diplotypes) = &result.kind else {
            panic!("expected UGT1A1 diplotype call, got {:?}", result.kind);
        };
        assert_eq!(
            diplotypes
                .iter()
                .map(|diplotype| diplotype.name.as_str())
                .collect::<Vec<_>>(),
            ["*1/*80", "*1/*80+*28", "*1/*80+*37"]
        );

        let ugt1a1 = run.context.gene_report("UGT1A1").expect("UGT1A1 report");
        let recommendation_labels = ugt1a1
            .recommendation_diplotypes
            .iter()
            .map(|diplotype| diplotype.label.as_str())
            .collect::<Vec<_>>();
        for expected in ["*1/*80", "*1/*80+*28", "*1/*80+*37"] {
            assert!(
                recommendation_labels.contains(&expected),
                "missing recommendation {expected} in {recommendation_labels:?}"
            );
        }

        // Java selects two .cpic-guideline-atazanavir rows: the first lists *1/*80 alone, the
        // second groups *1/*80+*28 and *1/*80+*37 under one recommendation.
        let html = fs::read_to_string(outputs.reporter_html.as_ref().unwrap()).expect("html");
        let rows = atazanavir_ugt1a1_rx_dip_rows(&html);
        assert_eq!(
            rows,
            vec![
                vec!["*1/*80".to_owned()],
                vec!["*1/*80+*28".to_owned(), "*1/*80+*37".to_owned()],
            ]
        );
    }

    /// Extract the ordered UGT1A1 `rx-dip` diplotype labels per `cpic-guideline-atazanavir` row,
    /// mirroring Java's `document.select(".cpic-guideline-atazanavir") .. .select(".rx-dip")`.
    fn atazanavir_ugt1a1_rx_dip_rows(html: &str) -> Vec<Vec<String>> {
        const DIP_PAT: &str = "rx-dip\"><a href=\"#UGT1A1\">UGT1A1</a>:";
        html.split("cpic-guideline-atazanavir")
            .skip(1)
            .map(|section| {
                let row = section.split("</tr>").next().unwrap_or(section);
                row.match_indices(DIP_PAT)
                    .map(|(index, pattern)| {
                        row[index + pattern.len()..]
                            .split('<')
                            .next()
                            .unwrap_or("")
                            .to_owned()
                    })
                    .collect::<Vec<_>>()
            })
            .collect()
    }

    #[test]
    fn run_reporter_from_vcf_ugt1a1_na12717_like_java_pipeline_test() {
        let definition =
            read_definition_file(Path::new(UGT1A1_DEFINITION_PATH)).expect("UGT1A1 definition");
        let definitions = DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition.clone())]
                .into_iter()
                .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");

        let mut vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n".to_owned();
        // Homozygous *80 (rs887829 T/T) with a heterozygous repeat (CAT/CATAT) splits the haplotypes
        // into *80 and *80+*28.
        append_definition_vcf_rows(
            &mut vcf,
            &definition,
            &[("rs887829", "T", "1/1"), ("rs3064744", "CATAT", "0/1")],
            &[],
        );
        let vcf_file = write_temp_named_file("pipeline-ugt1a1-na12717.vcf", &vcf);

        let run = run_reporter_from_vcf(
            &vcf_file,
            Some("PharmCAT"),
            &definitions,
            &phenotypes,
            &guidance,
            None,
            &ReporterPipelineOptions {
                include_combinations: false,
                html: HtmlReportOptions {
                    compact: true,
                    ..HtmlReportOptions::default()
                },
                ..ReporterPipelineOptions::default()
            },
        )
        .expect("UGT1A1 *80/*80+*28 pipeline run");

        let result = run
            .gene_call_results
            .iter()
            .find(|result| result.gene == "UGT1A1")
            .expect("UGT1A1 matcher result");
        let GeneCallKind::Diplotypes(diplotypes) = &result.kind else {
            panic!("expected UGT1A1 diplotype call, got {:?}", result.kind);
        };
        assert_eq!(
            diplotypes
                .iter()
                .map(|diplotype| diplotype.name.as_str())
                .collect::<Vec<_>>(),
            ["*80/*80+*28"]
        );

        let ugt1a1 = run.context.gene_report("UGT1A1").expect("UGT1A1 report");
        assert_eq!(ugt1a1.source_diplotype.as_deref(), Some("*80/*80+*28"));
        let recommendation = &ugt1a1.recommendation_diplotypes[0];
        let mut alleles = [
            recommendation.allele1.as_ref().map(|h| h.name.as_str()),
            recommendation.allele2.as_ref().map(|h| h.name.as_str()),
        ];
        alleles.sort();
        assert_eq!(alleles, [Some("*80"), Some("*80+*28")]);
    }

    #[test]
    fn run_reporter_from_vcf_ugt1a1_s28_hom_missing_like_java_pipeline_test() {
        let definition =
            read_definition_file(Path::new(UGT1A1_DEFINITION_PATH)).expect("UGT1A1 definition");
        let definitions = DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition.clone())]
                .into_iter()
                .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");

        let mut vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n".to_owned();
        // Homozygous *28 (rs3064744 CATAT/CATAT) with rs887829 missing leaves *80 status unknown, so
        // each *28 could also be *80+*28.
        append_definition_vcf_rows(
            &mut vcf,
            &definition,
            &[("rs3064744", "CATAT", "1/1")],
            &["rs887829"],
        );
        let vcf_file = write_temp_named_file("pipeline-ugt1a1-s28-hom-missing.vcf", &vcf);

        let run = run_reporter_from_vcf(
            &vcf_file,
            Some("PharmCAT"),
            &definitions,
            &phenotypes,
            &guidance,
            None,
            &ReporterPipelineOptions {
                include_combinations: false,
                html: HtmlReportOptions {
                    compact: true,
                    ..HtmlReportOptions::default()
                },
                ..ReporterPipelineOptions::default()
            },
        )
        .expect("UGT1A1 *28 hom missing pipeline run");

        let result = run
            .gene_call_results
            .iter()
            .find(|result| result.gene == "UGT1A1")
            .expect("UGT1A1 matcher result");
        let GeneCallKind::Diplotypes(diplotypes) = &result.kind else {
            panic!("expected UGT1A1 diplotype call, got {:?}", result.kind);
        };
        assert_eq!(
            diplotypes
                .iter()
                .map(|diplotype| diplotype.name.as_str())
                .collect::<Vec<_>>(),
            ["*28/*28", "*28/*80+*28", "*80+*28/*80+*28"]
        );

        let ugt1a1 = run.context.gene_report("UGT1A1").expect("UGT1A1 report");
        assert_eq!(
            ugt1a1
                .source_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["*28/*28", "*28/*80+*28", "*80+*28/*80+*28"]
        );
        let recommendation_labels = ugt1a1
            .recommendation_diplotypes
            .iter()
            .map(|diplotype| diplotype.label.as_str())
            .collect::<Vec<_>>();
        for expected in ["*28/*28", "*28/*80+*28", "*80+*28/*80+*28"] {
            assert!(
                recommendation_labels.contains(&expected),
                "missing recommendation {expected} in {recommendation_labels:?}"
            );
        }
    }

    #[test]
    fn run_reporter_from_vcf_ugt1a1_s1_s28_like_java_pipeline_test() {
        let definition =
            read_definition_file(Path::new(UGT1A1_DEFINITION_PATH)).expect("UGT1A1 definition");
        let definitions = DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition.clone())]
                .into_iter()
                .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");

        let mut vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n".to_owned();
        append_definition_vcf_rows(&mut vcf, &definition, &[("rs3064744", "CATAT", "0/1")], &[]);
        let vcf_file = write_temp_named_file("pipeline-ugt1a1-s1-s28.vcf", &vcf);

        let run = run_reporter_from_vcf(
            &vcf_file,
            Some("PharmCAT"),
            &definitions,
            &phenotypes,
            &guidance,
            None,
            &ReporterPipelineOptions {
                include_combinations: false,
                html: HtmlReportOptions {
                    compact: true,
                    ..HtmlReportOptions::default()
                },
                ..ReporterPipelineOptions::default()
            },
        )
        .expect("UGT1A1 *1/*28 pipeline run");

        let result = run
            .gene_call_results
            .iter()
            .find(|result| result.gene == "UGT1A1")
            .expect("UGT1A1 matcher result");
        let GeneCallKind::Diplotypes(diplotypes) = &result.kind else {
            panic!("expected UGT1A1 diplotype call, got {:?}", result.kind);
        };
        assert_eq!(
            diplotypes
                .iter()
                .map(|diplotype| diplotype.name.as_str())
                .collect::<Vec<_>>(),
            ["*1/*28"]
        );

        let ugt1a1 = run.context.gene_report("UGT1A1").expect("UGT1A1 report");
        assert_eq!(ugt1a1.source_diplotype.as_deref(), Some("*1/*28"));
        let recommendation = &ugt1a1.recommendation_diplotypes[0];
        let mut alleles = [
            recommendation.allele1.as_ref().map(|h| h.name.as_str()),
            recommendation.allele2.as_ref().map(|h| h.name.as_str()),
        ];
        alleles.sort();
        assert_eq!(alleles, [Some("*1"), Some("*28")]);
    }

    #[test]
    fn run_reporter_from_vcf_ugt1a1_s27_s28_unphased_s80_missing_like_java_pipeline_test() {
        let definition =
            read_definition_file(Path::new(UGT1A1_DEFINITION_PATH)).expect("UGT1A1 definition");
        let definitions = DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition.clone())]
                .into_iter()
                .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");

        let mut vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n".to_owned();
        // *27 (rs35350960 A) over a het *28 repeat, with rs887829 missing so *80 status is unknown.
        append_definition_vcf_rows(
            &mut vcf,
            &definition,
            &[("rs3064744", "CATAT", "0/1"), ("rs35350960", "A", "0/1")],
            &["rs887829"],
        );
        let vcf_file = write_temp_named_file("pipeline-ugt1a1-s27-s28-s80-missing.vcf", &vcf);

        let run = run_reporter_from_vcf(
            &vcf_file,
            Some("PharmCAT"),
            &definitions,
            &phenotypes,
            &guidance,
            None,
            &ReporterPipelineOptions {
                include_combinations: false,
                html: HtmlReportOptions {
                    compact: true,
                    ..HtmlReportOptions::default()
                },
                ..ReporterPipelineOptions::default()
            },
        )
        .expect("UGT1A1 *27/*28 missing pipeline run");

        let result = run
            .gene_call_results
            .iter()
            .find(|result| result.gene == "UGT1A1")
            .expect("UGT1A1 matcher result");
        let GeneCallKind::Diplotypes(diplotypes) = &result.kind else {
            panic!("expected UGT1A1 diplotype call, got {:?}", result.kind);
        };
        assert_eq!(
            diplotypes
                .iter()
                .map(|diplotype| diplotype.name.as_str())
                .collect::<Vec<_>>(),
            ["*27/*28", "*27/*80+*28"]
        );

        let ugt1a1 = run.context.gene_report("UGT1A1").expect("UGT1A1 report");
        assert_eq!(
            ugt1a1
                .source_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["*27/*28", "*27/*80+*28"]
        );
        let recommendation_labels = ugt1a1
            .recommendation_diplotypes
            .iter()
            .map(|diplotype| diplotype.label.as_str())
            .collect::<Vec<_>>();
        for expected in ["*27/*28", "*27/*80+*28"] {
            assert!(
                recommendation_labels.contains(&expected),
                "missing recommendation {expected} in {recommendation_labels:?}"
            );
        }
    }

    #[test]
    fn run_reporter_from_vcf_ugt1a1_s6_s80_s28_missing_unphased_like_java_pipeline_test() {
        let definition =
            read_definition_file(Path::new(UGT1A1_DEFINITION_PATH)).expect("UGT1A1 definition");
        let definitions = DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition.clone())]
                .into_iter()
                .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");

        let mut vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n".to_owned();
        // Unphased twin of the missing-phased case: rs3064744 missing leaves *80/*80+*28/*80+*37.
        append_definition_vcf_rows(
            &mut vcf,
            &definition,
            &[("rs887829", "T", "0/1"), ("rs4148323", "A", "0/1")],
            &["rs3064744"],
        );
        let vcf_file =
            write_temp_named_file("pipeline-ugt1a1-s6-s80-s28-missing-unphased.vcf", &vcf);

        let run = run_reporter_from_vcf(
            &vcf_file,
            Some("PharmCAT"),
            &definitions,
            &phenotypes,
            &guidance,
            None,
            &ReporterPipelineOptions {
                include_combinations: false,
                html: HtmlReportOptions {
                    compact: true,
                    ..HtmlReportOptions::default()
                },
                ..ReporterPipelineOptions::default()
            },
        )
        .expect("UGT1A1 *6/*80 missing unphased pipeline run");

        let result = run
            .gene_call_results
            .iter()
            .find(|result| result.gene == "UGT1A1")
            .expect("UGT1A1 matcher result");
        assert!(!result.match_data.phased);
        let GeneCallKind::Diplotypes(diplotypes) = &result.kind else {
            panic!("expected UGT1A1 diplotype call, got {:?}", result.kind);
        };
        assert_eq!(
            diplotypes
                .iter()
                .map(|diplotype| diplotype.name.as_str())
                .collect::<Vec<_>>(),
            ["*6/*80", "*6/*80+*28", "*6/*80+*37"]
        );

        let ugt1a1 = run.context.gene_report("UGT1A1").expect("UGT1A1 report");
        assert_eq!(
            ugt1a1
                .source_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["*6/*80", "*6/*80+*28", "*6/*80+*37"]
        );
        let recommendation_labels = ugt1a1
            .recommendation_diplotypes
            .iter()
            .map(|diplotype| diplotype.label.as_str())
            .collect::<Vec<_>>();
        for expected in ["*6/*80", "*6/*80+*28", "*6/*80+*37"] {
            assert!(
                recommendation_labels.contains(&expected),
                "missing recommendation {expected} in {recommendation_labels:?}"
            );
        }
    }

    #[test]
    fn run_reporter_from_vcf_ugt1a1_s6_s80_s28_unphased_like_java_pipeline_test() {
        let definition =
            read_definition_file(Path::new(UGT1A1_DEFINITION_PATH)).expect("UGT1A1 definition");
        let definitions = DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition.clone())]
                .into_iter()
                .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");

        let mut vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n".to_owned();
        append_definition_vcf_rows(
            &mut vcf,
            &definition,
            &[
                ("rs887829", "T", "0/1"),
                ("rs3064744", "CATAT", "0/1"),
                ("rs4148323", "A", "0/1"),
            ],
            &[],
        );
        let vcf_file = write_temp_named_file("pipeline-ugt1a1-s6-s80-s28-unphased.vcf", &vcf);

        let run = run_reporter_from_vcf(
            &vcf_file,
            Some("PharmCAT"),
            &definitions,
            &phenotypes,
            &guidance,
            None,
            &ReporterPipelineOptions {
                include_combinations: false,
                html: HtmlReportOptions {
                    compact: true,
                    ..HtmlReportOptions::default()
                },
                ..ReporterPipelineOptions::default()
            },
        )
        .expect("UGT1A1 *6/*80+*28 unphased pipeline run");

        let result = run
            .gene_call_results
            .iter()
            .find(|result| result.gene == "UGT1A1")
            .expect("UGT1A1 matcher result");
        assert!(!result.match_data.phased);
        let GeneCallKind::Diplotypes(diplotypes) = &result.kind else {
            panic!("expected UGT1A1 diplotype call, got {:?}", result.kind);
        };
        assert_eq!(
            diplotypes
                .iter()
                .map(|diplotype| diplotype.name.as_str())
                .collect::<Vec<_>>(),
            ["*6/*80+*28"]
        );

        let ugt1a1 = run.context.gene_report("UGT1A1").expect("UGT1A1 report");
        assert_eq!(ugt1a1.source_diplotype.as_deref(), Some("*6/*80+*28"));
        let recommendation = &ugt1a1.recommendation_diplotypes[0];
        let mut alleles = [
            recommendation.allele1.as_ref().map(|h| h.name.as_str()),
            recommendation.allele2.as_ref().map(|h| h.name.as_str()),
        ];
        alleles.sort();
        assert_eq!(alleles, [Some("*6"), Some("*80+*28")]);
    }

    #[test]
    fn run_cli_config_wires_parsed_vcf_args_to_reporter_outputs() {
        let definition_dir = unique_temp_path("pharmcat-cli-definitions");
        fs::create_dir_all(&definition_dir).expect("definition dir");
        fs::copy(
            Path::new(CYP3A5_DEFINITION_PATH),
            definition_dir.join("CYP3A5_translation.json"),
        )
        .expect("copy definition");
        let output_dir = unique_temp_path("pharmcat-cli-output");
        let sample_metadata = write_temp_named_file(
            "pharmcat-cli-sample-metadata.tsv",
            "NA12878\tStudy\tPGx\nNA12878\tBatch\tB1\nOther\tStudy\tignored\n",
        );
        let outside_file = write_temp_named_file("pharmcat-cli.outside.tsv", "CYP2D6\t*3/*4\n");
        let action = parse_pharmcat_args([
            "-vcf",
            CYP3A5_VCF_PATH,
            "-def",
            definition_dir.to_str().unwrap(),
            "-sm",
            sample_metadata.to_str().unwrap(),
            "-po",
            outside_file.to_str().unwrap(),
            "-o",
            output_dir.to_str().unwrap(),
            "-bf",
            "cli-cyp3a5",
            "-matcherHtml",
            "-reporterHtml",
            "-reporterJson",
            "-reporterCallsOnlyTsv",
        ])
        .expect("parse");
        let CliAction::Run(config) = action else {
            panic!("expected run");
        };
        let options = CliPipelineOptions {
            resources: PharmcatResourcePaths {
                definitions_dir: definition_dir,
                phenotype_dir: PathBuf::from(PHENOTYPE_PATH),
                prescribing_guidance: PathBuf::from(GUIDANCE_PATH),
                reporter_messages: PathBuf::from(MESSAGE_PATH),
            },
            mode: PipelineMode::Cli,
        };

        let runs = run_cli_config(&config, &options).expect("CLI pipeline");

        assert_eq!(runs.len(), 1);
        let run = &runs[0];
        assert_eq!(run.gene_call_results.len(), 1);
        assert_eq!(run.context.title.as_deref(), Some("cli-cyp3a5"));
        assert_eq!(
            run.context
                .gene_report("CYP3A5")
                .expect("CYP3A5 report")
                .source_diplotype
                .as_deref(),
            Some("*1/*2")
        );
        let cyp2d6 = run.context.gene_report("CYP2D6").expect("CYP2D6 report");
        assert!(cyp2d6.outside_call);
        assert_eq!(cyp2d6.source_diplotype.as_deref(), Some("*3/*4"));
        assert!(output_dir.join("cli-cyp3a5.match.json").is_file());
        assert!(output_dir.join("cli-cyp3a5.match.html").is_file());
        assert!(output_dir.join("cli-cyp3a5.match_warnings.txt").is_file());
        assert!(output_dir.join("cli-cyp3a5.phenotype.json").is_file());
        assert!(output_dir.join("cli-cyp3a5.report.html").is_file());
        assert!(output_dir.join("cli-cyp3a5.report.json").is_file());
        assert!(output_dir.join("cli-cyp3a5.report.tsv").is_file());

        let matcher_json =
            fs::read_to_string(output_dir.join("cli-cyp3a5.match.json")).expect("matcher JSON");
        let matcher_json: serde_json::Value =
            serde_json::from_str(&matcher_json).expect("parse matcher JSON");
        assert_eq!(
            matcher_json["metadata"]["inputFilename"].as_str(),
            Some("haplotyper.vcf")
        );
        assert_eq!(
            matcher_json["metadata"]["genomeBuild"].as_str(),
            Some("GRCh38.p13")
        );
        assert_eq!(
            matcher_json["metadata"]["sampleId"].as_str(),
            Some("NA12878")
        );
        assert!(
            matcher_json["metadata"]["timestamp"]
                .as_str()
                .is_some_and(|timestamp| timestamp.contains('T') && timestamp.ends_with('Z'))
        );
        assert_eq!(
            matcher_json["metadata"]["sampleProps"]["Study"].as_str(),
            Some("PGx")
        );
        assert_eq!(
            matcher_json["metadata"]["sampleProps"]["Batch"].as_str(),
            Some("B1")
        );
        assert!(matcher_json["metadata"].get("sample").is_none());
        assert!(matcher_json["vcfWarnings"].is_object());
        let gene_call = &matcher_json["results"][0];
        assert_eq!(gene_call["gene"].as_str(), Some("CYP3A5"));
        assert!(gene_call.get("dpydHapB3Warnings").is_none());
        assert_eq!(gene_call["diplotypes"][0]["name"].as_str(), Some("*1/*2"));
        assert_eq!(gene_call["haplotypes"][0]["name"].as_str(), Some("*1"));
        let haplotype = &gene_call["haplotypes"][0]["haplotype"];
        assert!(haplotype["corePositions"].is_array());
        assert!(haplotype["score"].is_number());
        assert_eq!(haplotype["structuralVariant"].as_bool(), Some(false));
        assert!(haplotype.get("populationFrequency").is_none());
        assert!(
            gene_call["haplotypeMatches"]
                .as_array()
                .is_some_and(Vec::is_empty)
        );
        assert_eq!(gene_call["variants"][0]["vcfCall"].as_str(), Some("C|T"));
        assert!(gene_call["variantsOfInterest"].as_array().is_some());
        assert!(gene_call["uncallableHaplotypes"].as_array().is_some());
        assert!(gene_call["warnings"].as_array().is_some_and(Vec::is_empty));
        let match_data = &gene_call["matchData"];
        assert_eq!(match_data["phased"].as_bool(), Some(true));
        assert_eq!(match_data["effectivelyPhased"].as_bool(), Some(false));
        assert_eq!(
            match_data["treatUndocumentedVariationsAsReference"].as_bool(),
            Some(false)
        );
        assert!(match_data["phaseSets"].is_object());
        assert!(match_data["posToPhaseSet"].is_object());
        assert!(match_data.get("haploid").is_none());
        let matcher_html =
            fs::read_to_string(output_dir.join("cli-cyp3a5.match.html")).expect("matcher HTML");
        assert!(matcher_html.contains("<html class=\"no-js\" lang=\"en\">"));
        assert!(matcher_html.contains("PharmCAT Allele Call Report for NA12878 in haplotyper.vcf"));
        assert!(matcher_html.contains("<h3>CYP3A5</h3>"));
        assert!(
            matcher_html.contains("<table class=\"table table-striped table-hover table-sm\">")
        );
        assert!(matcher_html.contains("<th class=\"first\">Definition Position</th>"));
        assert!(matcher_html.contains("<th class=\"first\">VCF REF,ALTs</th>"));
        assert!(matcher_html.contains("<th>C,T</th>"));
        assert!(matcher_html.contains("<th class=\"first\">VCF Call</th>"));
        assert!(matcher_html.contains("<li>*1/*2 ("));
        assert!(matcher_html.contains("<th class=\"table-danger\">C|T</th>"));
        assert!(matcher_html.contains("<th class=\"first\">*1</th>"));
        let phenotype_json = fs::read_to_string(output_dir.join("cli-cyp3a5.phenotype.json"))
            .expect("phenotype JSON");
        assert!(phenotype_json.contains("\"geneSymbol\": \"CYP3A5\""));
    }

    #[test]
    fn run_cli_config_rejects_unported_stage_combinations_explicitly() {
        let config = PharmcatCliConfig {
            run_matcher: false,
            run_phenotyper: false,
            run_reporter: true,
            reporter_input: Some(PathBuf::from("sample.phenotype.json")),
            ..PharmcatCliConfig::default()
        };
        let options = CliPipelineOptions {
            resources: PharmcatResourcePaths::from_resource_root("missing-resources"),
            mode: PipelineMode::Cli,
        };

        let error = run_cli_config(&config, &options).expect_err("unsupported path");

        assert!(
            matches!(error, CliPipelineError::Unsupported(message) if message.contains("full -vcf"))
        );
    }

    #[test]
    fn run_reporter_from_vcf_attaches_vcf_warnings_to_report_variant_reports() {
        let definition =
            read_definition_file(Path::new(CYP3A5_DEFINITION_PATH)).expect("definition");
        let warning_variant = definition
            .variants
            .first()
            .expect("definition variant")
            .clone();
        let alt = if warning_variant.reference == "A" {
            "C"
        } else {
            "A"
        };
        let vcf = format!(
            "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878\n{}\t{}\t{}\t{}\t{}\t.\tPASS\t.\tGT\t./.\n",
            warning_variant.chromosome,
            warning_variant.position,
            warning_variant.rsid.as_deref().unwrap_or("."),
            warning_variant.reference,
            alt
        );
        let vcf_path = write_temp_named_file("warning.vcf", &vcf);
        let definitions = DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition)]
                .into_iter()
                .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");

        let run = run_reporter_from_vcf(
            &vcf_path,
            Some("NA12878"),
            &definitions,
            &phenotypes,
            &guidance,
            None,
            &ReporterPipelineOptions::default(),
        )
        .expect("pipeline run");

        let chr_position = warning_variant.vcf_chr_position();
        assert_eq!(
            run.vcf_warnings
                .get(&chr_position)
                .expect("pipeline warning")
                .iter()
                .next()
                .map(String::as_str),
            Some("Ignoring: no call (./.)")
        );
        let report_gene = run.context.gene_report("CYP3A5").expect("CYP3A5 report");
        let report_variant = report_gene
            .variant_reports
            .iter()
            .find(|variant| variant.position == Some(warning_variant.position as i64))
            .expect("warning variant report");
        assert_eq!(report_variant.warnings, ["Ignoring: no call (./.)"]);
    }

    #[test]
    fn run_reporter_from_vcf_handles_standard_and_lowest_function_genes_in_one_run() {
        let cyp3a5_definition =
            read_definition_file(Path::new(CYP3A5_DEFINITION_PATH)).expect("CYP3A5 definition");
        let dpyd_definition =
            read_definition_file(Path::new(DPYD_DEFINITION_PATH)).expect("DPYD definition");
        let ryr1_definition =
            read_definition_file(Path::new(RYR1_DEFINITION_PATH)).expect("RYR1 definition");
        let mut vcf = fs::read_to_string(Path::new(CYP3A5_VCF_PATH)).expect("CYP3A5 VCF");
        append_definition_vcf_rows(
            &mut vcf,
            &dpyd_definition,
            &[("rs67376798", "A", "0|1")],
            &[],
        );
        append_definition_vcf_rows(
            &mut vcf,
            &ryr1_definition,
            &[("rs193922746", "G", "0/1")],
            &["rs193922753"],
        );
        let vcf_path = write_temp_named_file("multi-gene.vcf", &vcf);
        let definitions = DefinitionReader::from_definitions(
            [
                (cyp3a5_definition.gene_symbol.clone(), cyp3a5_definition),
                (dpyd_definition.gene_symbol.clone(), dpyd_definition),
                (ryr1_definition.gene_symbol.clone(), ryr1_definition),
            ]
            .into_iter()
            .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");
        let options = ReporterPipelineOptions {
            include_combinations: false,
            ..ReporterPipelineOptions::default()
        };

        let run = run_reporter_from_vcf(
            &vcf_path,
            Some("NA12878"),
            &definitions,
            &phenotypes,
            &guidance,
            None,
            &options,
        )
        .expect("pipeline run");

        assert_eq!(run.gene_call_results.len(), 3);

        let cyp3a5 = run.context.gene_report("CYP3A5").expect("CYP3A5 report");
        assert_eq!(cyp3a5.source_diplotype.as_deref(), Some("*1/*2"));

        let dpyd = run.context.gene_report("DPYD").expect("DPYD report");
        assert_eq!(
            dpyd.source_diplotype.as_deref(),
            Some("c.2846A>T (heterozygous)")
        );
        assert_eq!(dpyd.lookup_keys, ["1.5"]);
        assert!(
            dpyd.variant_reports
                .iter()
                .any(|variant| variant.db_snp_id.as_deref() == Some("rs67376798"))
        );

        let ryr1 = run.context.gene_report("RYR1").expect("RYR1 report");
        assert_eq!(
            ryr1.source_diplotype.as_deref(),
            Some("c.97A>G (heterozygous)")
        );
        assert_eq!(ryr1.lookup_keys, ["Malignant Hyperthermia Susceptibility"]);
        assert!(
            ryr1.variant_reports
                .iter()
                .any(|variant| variant.db_snp_id.as_deref() == Some("rs193922746"))
        );

        let warning_variant = variant_by_rsid(
            definitions
                .definition_file("RYR1")
                .expect("RYR1 definition"),
            "rs193922753",
        );
        assert_eq!(
            run.vcf_warnings
                .get(&warning_variant.vcf_chr_position())
                .expect("pipeline warning")
                .iter()
                .next()
                .map(String::as_str),
            Some("Ignoring: no call (./.)")
        );
        let report_variant = ryr1
            .variant_reports
            .iter()
            .find(|variant| variant.db_snp_id.as_deref() == Some("rs193922753"))
            .expect("warning variant report");
        assert_eq!(report_variant.warnings, ["Ignoring: no call (./.)"]);
    }

    #[test]
    fn run_reporter_from_vcf_surfaces_vcf_sample_errors() {
        let definitions = DefinitionReader::from_definitions(BTreeMap::new());
        let phenotypes = PhenotypeMap::default();
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");

        let error = run_reporter_from_vcf(
            Path::new(CYP3A5_VCF_PATH),
            Some("missing-sample"),
            &definitions,
            &phenotypes,
            &guidance,
            None,
            &ReporterPipelineOptions::default(),
        )
        .expect_err("missing sample");

        assert!(matches!(
            error,
            ReporterPipelineError::Vcf(ReadVcfError::SampleNotFound(sample))
                if sample == "missing-sample"
        ));
    }

    fn run_cyp2c19_undocumented_variation_pipeline(
        compact: bool,
        basename: &str,
    ) -> (DefinitionFile, PipelineOutputPlan, ReporterPipelineRun) {
        let definition =
            read_definition_file(Path::new(CYP2C19_DEFINITION_PATH)).expect("CYP2C19 definition");
        let definitions = DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition.clone())]
                .into_iter()
                .collect(),
        );
        let phenotypes = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let guidance =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");

        let mut vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n".to_owned();
        append_definition_vcf_rows(&mut vcf, &definition, &[("rs3758581", "G,T", "0/2")], &[]);
        let vcf_file = write_temp_named_file(&format!("{basename}.vcf"), &vcf);
        let output_dir = unique_temp_path(&format!("pharmcat-run-{basename}"));
        let outputs = PipelineOutputPlan {
            base_dir: output_dir.clone(),
            basename: basename.to_owned(),
            display_name: format!("{basename}.vcf"),
            reporter_title: Some("PharmCAT".to_owned()),
            matcher_json: None,
            matcher_html: None,
            matcher_warnings: None,
            phenotyper_json: None,
            reporter_html: Some(output_dir.join(format!("{basename}.report.html"))),
            reporter_json: Some(output_dir.join(format!("{basename}.report.json"))),
            reporter_calls_only_tsv: Some(output_dir.join(format!("{basename}.report.tsv"))),
        };

        let run = run_reporter_from_vcf(
            &vcf_file,
            Some("PharmCAT"),
            &definitions,
            &phenotypes,
            &guidance,
            Some(&outputs),
            &ReporterPipelineOptions {
                include_combinations: false,
                html: HtmlReportOptions {
                    compact,
                    ..HtmlReportOptions::default()
                },
                ..ReporterPipelineOptions::default()
            },
        )
        .expect("undocumented CYP2C19 pipeline run");

        (definition, outputs, run)
    }

    fn assert_undocumented_warning(
        warnings: &VcfWarnings,
        chr_position: &str,
        expected_text: &str,
        found_text: &str,
        treat_as_reference: bool,
    ) {
        let warning = warnings
            .get(chr_position)
            .unwrap_or_else(|| panic!("undocumented warning at {chr_position}"))
            .iter()
            .next()
            .unwrap_or_else(|| panic!("warning text at {chr_position}"));
        assert!(warning.contains(expected_text), "{warning}");
        assert!(warning.contains(found_text), "{warning}");
        assert_eq!(
            warning.contains("Undocumented variations will be replaced with reference."),
            treat_as_reference,
            "{warning}"
        );
    }

    fn write_temp_named_file(filename: &str, contents: &str) -> PathBuf {
        let dir = unique_temp_path("pharmcat-pipeline-test-dir");
        fs::create_dir_all(&dir).expect("create temp dir");
        let path = dir.join(filename);
        fs::write(&path, contents).expect("write temp file");
        path
    }

    fn append_definition_vcf_rows(
        vcf: &mut String,
        definition: &DefinitionFile,
        overrides: &[(&str, &str, &str)],
        missing_rsids: &[&str],
    ) {
        append_definition_vcf_rows_with_default_genotype(
            vcf,
            definition,
            overrides,
            missing_rsids,
            "0/0",
        );
    }

    fn append_definition_vcf_rows_with_default_genotype(
        vcf: &mut String,
        definition: &DefinitionFile,
        overrides: &[(&str, &str, &str)],
        missing_rsids: &[&str],
        default_genotype: &str,
    ) {
        if !vcf.ends_with('\n') {
            vcf.push('\n');
        }
        for variant in &definition.variants {
            let rsid = variant.rsid.as_deref().unwrap_or(".");
            if missing_rsids.contains(&rsid) {
                append_vcf_row(vcf, variant, &alternate_for(variant), "./.");
            } else if let Some((_, alt, genotype)) = overrides
                .iter()
                .find(|(override_rsid, _, _)| *override_rsid == rsid)
            {
                append_vcf_row(vcf, variant, alt, genotype);
            } else {
                append_vcf_row(vcf, variant, &alternate_for(variant), default_genotype);
            }
        }
    }

    fn append_definition_vcf_rows_omitting(
        vcf: &mut String,
        definition: &DefinitionFile,
        overrides: &[(&str, &str, &str)],
        omitted_rsids: &[&str],
    ) {
        if !vcf.ends_with('\n') {
            vcf.push('\n');
        }
        for variant in &definition.variants {
            let rsid = variant.rsid.as_deref().unwrap_or(".");
            if omitted_rsids.contains(&rsid) {
                continue;
            }
            if let Some((_, alt, genotype)) = overrides
                .iter()
                .find(|(override_rsid, _, _)| *override_rsid == rsid)
            {
                append_vcf_row(vcf, variant, alt, genotype);
            } else {
                append_vcf_row(vcf, variant, &alternate_for(variant), "0/0");
            }
        }
    }

    fn append_vcf_row(vcf: &mut String, variant: &VariantLocus, alt: &str, genotype: &str) {
        vcf.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t.\tPASS\t.\tGT\t{}\n",
            variant.chromosome,
            variant.position,
            variant.rsid.as_deref().unwrap_or("."),
            variant.reference,
            alt,
            genotype
        ));
    }

    fn alternate_for(variant: &VariantLocus) -> String {
        variant
            .alts
            .first()
            .cloned()
            .unwrap_or_else(|| match variant.reference.as_str() {
                "A" => "C".to_owned(),
                _ => "A".to_owned(),
            })
    }

    fn variant_by_rsid<'a>(definition: &'a DefinitionFile, rsid: &str) -> &'a VariantLocus {
        definition
            .variants
            .iter()
            .find(|variant| variant.rsid.as_deref() == Some(rsid))
            .unwrap_or_else(|| panic!("missing variant {rsid}"))
    }

    fn matcher_html_pattern_positions() -> Vec<VariantLocus> {
        (0..5)
            .map(|index| VariantLocus {
                chromosome: "chr1".to_owned(),
                position: 100 + index,
                cpic_position: 100 + index,
                rsid: Some(format!("rs-pattern-{index}")),
                chromosome_hgvs_name: format!("NC_000001.11:g.{}A>T", 100 + index),
                cpic_alleles: Default::default(),
                cpic_to_vcf_allele_map: Default::default(),
                reference: "A".to_owned(),
                alts: vec!["T".to_owned()],
            })
            .collect()
    }

    fn matched_annotation_count(context: &ReportContext, drug: &str) -> usize {
        PrescribingGuidanceSource::list_values()
            .iter()
            .filter_map(|source| context.drug_report(*source, drug))
            .map(|report| report.matched_annotation_count())
            .sum()
    }

    fn assert_drug_has_match_from_source(
        context: &ReportContext,
        drug: &str,
        source: PrescribingGuidanceSource,
    ) {
        assert!(
            context
                .drug_report(source, drug)
                .is_some_and(|report| report.matched_annotation_count() > 0),
            "{drug} {source:?}"
        );
    }

    fn assert_drug_has_no_match_from_source(
        context: &ReportContext,
        drug: &str,
        source: PrescribingGuidanceSource,
    ) {
        assert!(
            context
                .drug_report(source, drug)
                .is_none_or(|report| report.matched_annotation_count() == 0),
            "{drug} {source:?}"
        );
    }

    fn unique_temp_path(prefix: &str) -> PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("time")
            .as_nanos();
        std::env::temp_dir().join(format!("{prefix}-{nanos}"))
    }
}
