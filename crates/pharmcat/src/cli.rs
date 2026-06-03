//! Command-line parsing helpers.

use std::{
    collections::{BTreeMap, BTreeSet},
    fmt, fs, io,
    path::{Path, PathBuf},
};

use crate::report::PrescribingGuidanceSource;

/// PharmCAT CLI version string prefix used by Java.
pub const PHARMCAT_VERSION_PREFIX: &str = "PharmCAT";

const VCF_PREPROCESSED_SUFFIX: &str = ".preprocessed";
const MATCHER_SUFFIX: &str = ".match";
const PHENOTYPER_SUFFIX: &str = ".phenotype";
const REPORTER_SUFFIX: &str = ".report";

/// Parsed action for the top-level PharmCAT CLI.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum CliAction {
    /// Run PharmCAT with parsed configuration.
    Run(Box<PharmcatCliConfig>),
    /// Print help and exit successfully.
    Help,
    /// Print version and exit successfully.
    Version,
}

/// First Rust port of Java `BaseConfig`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PharmcatCliConfig {
    /// Run the named allele matcher stage.
    pub run_matcher: bool,
    /// Run the phenotyper stage.
    pub run_phenotyper: bool,
    /// Run the reporter stage.
    pub run_reporter: bool,
    /// Optional definitions directory.
    pub definition_dir: Option<PathBuf>,
    /// Return only top candidate matcher results.
    pub top_candidate_only: bool,
    /// Enable research-mode combination calling.
    pub find_combinations: bool,
    /// Enable research-mode CYP2D6 calling.
    pub call_cyp2d6: bool,
    /// Save matcher HTML.
    pub matcher_html: bool,
    /// Reporter title.
    pub reporter_title: Option<String>,
    /// Use compact reporter output.
    pub reporter_compact: bool,
    /// Reporter source filters.
    pub reporter_sources: Option<Vec<PrescribingGuidanceSource>>,
    /// Save reporter JSON.
    pub reporter_json: bool,
    /// Save reporter HTML.
    pub reporter_html: bool,
    /// Save reporter calls-only TSV.
    pub reporter_calls_only_tsv: bool,
    /// Output directory.
    pub output_dir: Option<PathBuf>,
    /// Base output filename.
    pub base_filename: Option<String>,
    /// Delete intermediate files.
    pub delete_intermediate_files: bool,
    /// Verbose output.
    pub verbose: bool,
    /// Selected samples, sorted and deduplicated like Java `TreeSet`.
    pub samples: BTreeSet<String>,
    /// Sample metadata file.
    pub sample_metadata_file: Option<PathBuf>,
    /// Selected genes, uppercased and sorted like Java `TreeSet`.
    pub genes: BTreeSet<String>,
    /// Matcher VCF input.
    pub matcher_vcf: Option<PathBuf>,
    /// Phenotyper JSON input.
    pub phenotyper_input: Option<PathBuf>,
    /// Phenotyper outside-call TSV inputs.
    pub phenotyper_outside_call_files: Vec<PathBuf>,
    /// Reporter JSON input.
    pub reporter_input: Option<PathBuf>,
}

/// Output paths derived the same way Java `Pipeline` derives them before running.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PipelineOutputPlan {
    /// Directory where generated files will be written.
    pub base_dir: PathBuf,
    /// Java `Pipeline.m_basename`.
    pub basename: String,
    /// Java `Pipeline.m_displayName`.
    pub display_name: String,
    /// Reporter title, defaulting to `basename`.
    pub reporter_title: Option<String>,
    /// Named allele matcher JSON output path.
    pub matcher_json: Option<PathBuf>,
    /// Named allele matcher HTML output path.
    pub matcher_html: Option<PathBuf>,
    /// Named allele matcher VCF warning output path.
    pub matcher_warnings: Option<PathBuf>,
    /// Phenotyper JSON output path.
    pub phenotyper_json: Option<PathBuf>,
    /// Reporter HTML output path.
    pub reporter_html: Option<PathBuf>,
    /// Reporter JSON output path.
    pub reporter_json: Option<PathBuf>,
    /// Reporter calls-only TSV output path.
    pub reporter_calls_only_tsv: Option<PathBuf>,
}

impl Default for PharmcatCliConfig {
    fn default() -> Self {
        Self {
            run_matcher: true,
            run_phenotyper: true,
            run_reporter: true,
            definition_dir: None,
            top_candidate_only: true,
            find_combinations: false,
            call_cyp2d6: false,
            matcher_html: false,
            reporter_title: None,
            reporter_compact: true,
            reporter_sources: None,
            reporter_json: false,
            reporter_html: true,
            reporter_calls_only_tsv: false,
            output_dir: None,
            base_filename: None,
            delete_intermediate_files: false,
            verbose: false,
            samples: BTreeSet::new(),
            sample_metadata_file: None,
            genes: BTreeSet::new(),
            matcher_vcf: None,
            phenotyper_input: None,
            phenotyper_outside_call_files: Vec::new(),
            reporter_input: None,
        }
    }
}

/// CLI parse error.
#[derive(Debug)]
pub enum CliError {
    /// Unknown argument.
    UnknownOption(String),
    /// Missing value for an option.
    MissingValue(&'static str),
    /// Invalid argument combination or value.
    Invalid(String),
    /// File I/O while reading a referenced input.
    Io { path: PathBuf, source: io::Error },
}

impl fmt::Display for CliError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnknownOption(option) => write!(f, "Unknown option: {option}"),
            Self::MissingValue(option) => write!(f, "Missing value for {option}"),
            Self::Invalid(message) => f.write_str(message),
            Self::Io { path, source } => write!(f, "{}: {source}", path.display()),
        }
    }
}

impl std::error::Error for CliError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io { source, .. } => Some(source),
            _ => None,
        }
    }
}

#[derive(Default)]
struct ParsedArgs {
    values: BTreeMap<&'static str, Vec<String>>,
    flags: BTreeSet<&'static str>,
}

impl ParsedArgs {
    fn has(&self, key: &'static str) -> bool {
        self.flags.contains(key) || self.values.contains_key(key)
    }

    fn value(&self, key: &'static str) -> Option<&str> {
        self.values
            .get(key)
            .and_then(|values| values.last())
            .map(String::as_str)
    }

    fn values(&self, key: &'static str) -> impl Iterator<Item = &str> {
        self.values
            .get(key)
            .into_iter()
            .flat_map(|values| values.iter().map(String::as_str))
    }
}

/// Parse Java-compatible PharmCAT CLI args.
pub fn parse_pharmcat_args<I, S>(args: I) -> Result<CliAction, CliError>
where
    I: IntoIterator<Item = S>,
    S: Into<String>,
{
    let parsed = parse_raw_args(args)?;
    if parsed.has("help") {
        return Ok(CliAction::Help);
    }
    if parsed.has("version") {
        return Ok(CliAction::Version);
    }
    config_from_parsed(parsed).map(|config| CliAction::Run(Box::new(config)))
}

/// Minimal help text for the currently ported CLI surface.
pub fn help_text() -> &'static str {
    concat!(
        "Usage: pharmcat -vcf <file> [options]\n\n",
        "Inputs:\n",
        "  -vcf, --matcher-vcf <file>\n",
        "  -s, --samples <samples>\n",
        "  -S, --sample-file <file>\n",
        "  -g, --genes <genes>\n\n",
        "Stages:\n",
        "  -matcher, --matcher\n",
        "  -phenotyper, --phenotyper\n",
        "  -reporter, --reporter\n\n",
        "Reporter:\n",
        "  -reporterHtml, --reporter-save-html\n",
        "  -reporterJson, --reporter-save-json\n",
        "  -reporterCallsOnlyTsv, --reporter-save-calls-only-tsv\n"
    )
}

/// Builds Java `Pipeline`-style output paths for a parsed CLI config.
pub fn pipeline_output_plan(
    config: &PharmcatCliConfig,
    sample_id: Option<&str>,
    single_sample: bool,
) -> Result<PipelineOutputPlan, CliError> {
    let mut builder = PipelineOutputPlanBuilder::new(config.output_dir.clone());

    if config.run_matcher {
        let vcf_file = config
            .matcher_vcf
            .as_deref()
            .ok_or_else(|| CliError::Invalid("No matcher VCF input file".to_string()))?;
        builder.generate_basename(
            config.base_filename.as_deref(),
            vcf_file,
            sample_id,
            single_sample,
        )?;
        let base_dir = builder
            .base_dir
            .as_ref()
            .expect("basename generation sets base dir")
            .clone();
        let basename = builder.basename.as_ref().expect("basename set").clone();
        builder.matcher_json = Some(base_dir.join(format!("{basename}{MATCHER_SUFFIX}.json")));
        if config.matcher_html {
            builder.matcher_html = Some(base_dir.join(format!("{basename}{MATCHER_SUFFIX}.html")));
        }
        builder.matcher_warnings =
            Some(base_dir.join(format!("{basename}{MATCHER_SUFFIX}_warnings.txt")));
    }

    if config.run_phenotyper {
        let input_file = builder
            .matcher_json
            .clone()
            .or_else(|| config.phenotyper_input.clone())
            .or_else(|| config.phenotyper_outside_call_files.first().cloned())
            .ok_or_else(|| CliError::Invalid("No phenotyper input file".to_string()))?;
        builder.generate_basename(
            config.base_filename.as_deref(),
            &input_file,
            sample_id,
            single_sample,
        )?;
        let base_dir = builder
            .base_dir
            .as_ref()
            .expect("basename generation sets base dir")
            .clone();
        let basename = builder.basename.as_ref().expect("basename set").clone();
        builder.phenotyper_json =
            Some(base_dir.join(format!("{basename}{PHENOTYPER_SUFFIX}.json")));
    }

    if config.run_reporter {
        let input_file = builder
            .phenotyper_json
            .clone()
            .or_else(|| config.reporter_input.clone())
            .ok_or_else(|| CliError::Invalid("No reporter input file".to_string()))?;
        builder.generate_basename(
            config.base_filename.as_deref(),
            &input_file,
            sample_id,
            single_sample,
        )?;
        let base_dir = builder
            .base_dir
            .as_ref()
            .expect("basename generation sets base dir")
            .clone();
        let basename = builder.basename.as_ref().expect("basename set").clone();
        if config.reporter_html {
            builder.reporter_html =
                Some(base_dir.join(format!("{basename}{REPORTER_SUFFIX}.html")));
        }
        if config.reporter_json {
            builder.reporter_json =
                Some(base_dir.join(format!("{basename}{REPORTER_SUFFIX}.json")));
        }
        if config.reporter_calls_only_tsv {
            builder.reporter_calls_only_tsv =
                Some(base_dir.join(format!("{basename}{REPORTER_SUFFIX}.tsv")));
        }
        builder.reporter_title = Some(
            config
                .reporter_title
                .clone()
                .unwrap_or_else(|| basename.to_string()),
        );
    }

    builder.finish()
}

#[derive(Debug)]
struct PipelineOutputPlanBuilder {
    base_dir: Option<PathBuf>,
    basename: Option<String>,
    display_name: Option<String>,
    reporter_title: Option<String>,
    matcher_json: Option<PathBuf>,
    matcher_html: Option<PathBuf>,
    matcher_warnings: Option<PathBuf>,
    phenotyper_json: Option<PathBuf>,
    reporter_html: Option<PathBuf>,
    reporter_json: Option<PathBuf>,
    reporter_calls_only_tsv: Option<PathBuf>,
}

impl PipelineOutputPlanBuilder {
    fn new(base_dir: Option<PathBuf>) -> Self {
        Self {
            base_dir,
            basename: None,
            display_name: None,
            reporter_title: None,
            matcher_json: None,
            matcher_html: None,
            matcher_warnings: None,
            phenotyper_json: None,
            reporter_html: None,
            reporter_json: None,
            reporter_calls_only_tsv: None,
        }
    }

    fn generate_basename(
        &mut self,
        base_filename: Option<&str>,
        input_file: &Path,
        sample_id: Option<&str>,
        single_sample: bool,
    ) -> Result<(), CliError> {
        if self.base_dir.is_none() {
            self.base_dir = Some(base_dir(input_file)?);
        }
        if self.basename.is_some() {
            return Ok(());
        }

        let mut basename = base_filename
            .map(str::to_string)
            .unwrap_or_else(|| base_filename_from_path(input_file));
        let mut display_name = basename.clone();
        if let Some(sample_id) = sample_id {
            if !single_sample
                && basename != sample_id
                && !basename.starts_with(&format!("{sample_id}."))
                && !basename.contains(&format!(".{sample_id}."))
            {
                basename.push('.');
                basename.push_str(sample_id);
            }
            let input_name = input_file
                .file_name()
                .and_then(|name| name.to_str())
                .unwrap_or_default();
            display_name = format!("sample {sample_id} in {input_name}");
        }
        self.basename = Some(basename);
        self.display_name = Some(display_name);
        Ok(())
    }

    fn finish(self) -> Result<PipelineOutputPlan, CliError> {
        Ok(PipelineOutputPlan {
            base_dir: self.base_dir.ok_or_else(|| {
                CliError::Invalid(
                    "Cannot determine directory to save results to.  Please specify output directory."
                        .to_string(),
                )
            })?,
            basename: self.basename.ok_or_else(|| {
                CliError::Invalid("Cannot determine base filename for output files".to_string())
            })?,
            display_name: self.display_name.unwrap_or_default(),
            reporter_title: self.reporter_title,
            matcher_json: self.matcher_json,
            matcher_html: self.matcher_html,
            matcher_warnings: self.matcher_warnings,
            phenotyper_json: self.phenotyper_json,
            reporter_html: self.reporter_html,
            reporter_json: self.reporter_json,
            reporter_calls_only_tsv: self.reporter_calls_only_tsv,
        })
    }
}

fn parse_raw_args<I, S>(args: I) -> Result<ParsedArgs, CliError>
where
    I: IntoIterator<Item = S>,
    S: Into<String>,
{
    let mut parsed = ParsedArgs::default();
    let mut iter = args.into_iter().map(Into::into).peekable();
    while let Some(arg) = iter.next() {
        let Some(option) = normalize_option(&arg) else {
            return Err(CliError::UnknownOption(arg));
        };
        if option_takes_value(option) {
            let Some(value) = iter.next() else {
                return Err(CliError::MissingValue(display_option(option)));
            };
            if normalize_option(&value).is_some() {
                return Err(CliError::MissingValue(display_option(option)));
            }
            parsed.values.entry(option).or_default().push(value);
        } else {
            parsed.flags.insert(option);
        }
    }
    Ok(parsed)
}

fn normalize_option(arg: &str) -> Option<&'static str> {
    match arg {
        "-h" | "--help" | "-help" => Some("help"),
        "-version" | "--version" | "-V" => Some("version"),
        "-v" | "--verbose" | "-verbose" => Some("verbose"),
        "-s" | "--samples" => Some("s"),
        "-S" | "--sample-file" => Some("S"),
        "-sm" | "--sample-metadata" => Some("sm"),
        "-g" | "--genes" => Some("g"),
        "-matcher" | "--matcher" => Some("matcher"),
        "-vcf" | "--matcher-vcf" => Some("vcf"),
        "-ma" | "--matcher-all-results" => Some("ma"),
        "-matcherHtml" | "--matcher-save-html" => Some("matcherHtml"),
        "-phenotyper" | "--phenotyper" => Some("phenotyper"),
        "-pi" | "--phenotyper-input" => Some("pi"),
        "-po" | "--phenotyper-outside-call-file" => Some("po"),
        "-reporter" | "--reporter" => Some("reporter"),
        "-ri" | "--reporter-input" => Some("ri"),
        "-rt" | "--reporter-title" => Some("rt"),
        "-rs" | "--reporter-sources" => Some("rs"),
        "-re" | "--reporter-extended" => Some("re"),
        "-reporterHtml" | "--reporter-save-html" => Some("reporterHtml"),
        "-reporterJson" | "--reporter-save-json" => Some("reporterJson"),
        "-reporterCallsOnlyTsv" | "--reporter-save-calls-only-tsv" => Some("reporterCallsOnlyTsv"),
        "-o" | "--output-dir" => Some("o"),
        "-bf" | "--base-filename" => Some("bf"),
        "-del" | "--delete-intermediate-files" => Some("del"),
        "-def" | "--definitions-dir" => Some("def"),
        "-research" | "--research-mode" => Some("research"),
        _ => None,
    }
}

fn option_takes_value(option: &str) -> bool {
    matches!(
        option,
        "s" | "S"
            | "sm"
            | "g"
            | "vcf"
            | "pi"
            | "po"
            | "ri"
            | "rt"
            | "rs"
            | "o"
            | "bf"
            | "def"
            | "research"
    )
}

fn display_option(option: &'static str) -> &'static str {
    match option {
        "s" => "-s",
        "S" => "-S",
        "sm" => "-sm",
        "g" => "-g",
        "vcf" => "-vcf",
        "pi" => "-pi",
        "po" => "-po",
        "ri" => "-ri",
        "rt" => "-rt",
        "rs" => "-rs",
        "o" => "-o",
        "bf" => "-bf",
        "def" => "-def",
        "research" => "-research",
        _ => "-<option>",
    }
}

fn config_from_parsed(parsed: ParsedArgs) -> Result<PharmcatCliConfig, CliError> {
    let mut config = PharmcatCliConfig::default();

    if parsed.has("matcher") || parsed.has("phenotyper") || parsed.has("reporter") {
        config.run_matcher = parsed.has("matcher");
        config.run_phenotyper = parsed.has("phenotyper");
        config.run_reporter = parsed.has("reporter");
    }
    if config.run_matcher && !config.run_phenotyper && config.run_reporter {
        return Err(CliError::Invalid(
            "Cannot run matcher and reporter without also running phenotyper.".to_string(),
        ));
    }

    config.definition_dir = parsed
        .value("def")
        .map(|path| valid_directory(path, false))
        .transpose()?;

    if parsed.has("s") && parsed.has("S") {
        return Err(CliError::Invalid(
            "Cannot specify both -s and -S".to_string(),
        ));
    }
    for value in parsed.values("s") {
        config
            .samples
            .extend(split_comma(value).map(str::to_string));
    }
    if let Some(sample_file) = parsed.value("S").map(PathBuf::from) {
        read_sample_file(&sample_file, &mut config.samples)?;
    }
    config.sample_metadata_file = parsed
        .value("sm")
        .map(|path| valid_file(path, true))
        .transpose()?;
    for value in parsed.values("g") {
        config
            .genes
            .extend(split_comma(value).map(|gene| gene.to_uppercase()));
    }

    let mut research_mode = false;
    if config.run_matcher {
        config.top_candidate_only = !parsed.has("ma");
        if let Some(research) = parsed.value("research") {
            let mut unknown = Vec::new();
            for option in split_comma(research).map(str::to_lowercase) {
                match option.as_str() {
                    "cyp2d6" => {
                        config.call_cyp2d6 = true;
                        research_mode = true;
                    }
                    "combinations" | "combination" => {
                        config.find_combinations = true;
                        research_mode = true;
                    }
                    _ => unknown.push(option),
                }
            }
            if !unknown.is_empty() {
                return Err(CliError::Invalid(format!(
                    "Unrecognized research option: {}",
                    unknown.join(",")
                )));
            }
        }
        config.matcher_html = parsed.has("matcherHtml");
    }

    if config.run_reporter {
        config.reporter_title = parsed.value("rt").map(str::to_string);
        config.reporter_compact = !parsed.has("re");
        config.reporter_json = parsed.has("reporterJson");
        config.reporter_calls_only_tsv = parsed.has("reporterCallsOnlyTsv");
        if config.reporter_json || config.reporter_calls_only_tsv {
            config.reporter_html = parsed.has("reporterHtml");
        }
        if research_mode {
            if !config.reporter_calls_only_tsv {
                config.reporter_calls_only_tsv = true;
            }
            config.reporter_html = false;
            config.reporter_json = false;
        }
        if let Some(sources) = parsed.value("rs") {
            config.reporter_sources = Some(parse_reporter_sources(sources)?);
        }
    }

    config.output_dir = parsed
        .value("o")
        .map(|path| valid_directory(path, true))
        .transpose()?;
    config.base_filename = parsed.value("bf").map(str::to_string);
    config.delete_intermediate_files = parsed.has("del");
    config.verbose = parsed.has("verbose");
    config.matcher_vcf = parsed
        .value("vcf")
        .map(|path| valid_file(path, true))
        .transpose()?;
    config.phenotyper_input = parsed
        .value("pi")
        .map(|path| valid_file(path, true))
        .transpose()?;
    config.phenotyper_outside_call_files = parsed
        .values("po")
        .map(|path| {
            let path_buf = PathBuf::from(path);
            if path_buf.is_file() {
                Ok(path_buf)
            } else {
                Err(CliError::Invalid(format!("Not a valid file: '{path}")))
            }
        })
        .collect::<Result<_, _>>()?;
    config.reporter_input = parsed
        .value("ri")
        .map(|path| valid_file(path, true))
        .transpose()?;

    if config.run_matcher {
        if config.matcher_vcf.is_none() {
            return Err(CliError::Invalid(
                "No input for Named Allele Matcher!\n\nPlease specify a VCF file (-vcf)"
                    .to_string(),
            ));
        }
        if config.phenotyper_input.is_some() {
            return Err(CliError::Invalid(
                "Cannot specify phenotyper-input (-pi) if running named allele matcher".to_string(),
            ));
        }
    } else if !config.samples.is_empty() {
        return Err(CliError::Invalid(
            "Cannot specify samples unless running matcher.".to_string(),
        ));
    }

    if config.run_phenotyper
        && config.matcher_vcf.is_none()
        && config.phenotyper_input.is_none()
        && config.phenotyper_outside_call_files.is_empty()
    {
        return Err(CliError::Invalid(
            "No input for Phenotyper!\n\nEither:\n  1. Run named allele matcher with VCF input, or\n  2. Specify phenotyper-input (-pi) and/or phenotyper-outside-call-file (-po)"
                .to_string(),
        ));
    }

    if config.run_reporter
        && config.matcher_vcf.is_none()
        && config.phenotyper_input.is_none()
        && config.phenotyper_outside_call_files.is_empty()
        && config.reporter_input.is_none()
    {
        return Err(CliError::Invalid(
            "No input for Reporter!\n\nEither:\n  1. Run phenotyper, or\n  2. Specify reporter-input (-ri)"
                .to_string(),
        ));
    }

    Ok(config)
}

fn read_sample_file(path: &PathBuf, samples: &mut BTreeSet<String>) -> Result<(), CliError> {
    let path = valid_file(path, true)?;
    let contents = fs::read_to_string(&path).map_err(|source| CliError::Io { path, source })?;
    for line in contents.lines() {
        let sample = line.trim();
        if sample.is_empty() {
            continue;
        }
        if sample.contains(',') {
            return Err(CliError::Invalid(
                "Error: Please remove comma ',' from sample names".to_string(),
            ));
        }
        samples.insert(sample.to_string());
    }
    Ok(())
}

fn valid_file(path: impl AsRef<Path>, must_exist: bool) -> Result<PathBuf, CliError> {
    let path = path.as_ref();
    if !path.exists() {
        if must_exist {
            return Err(CliError::Invalid(format!(
                "File '{}' does not exist",
                path.display()
            )));
        }
    } else if !path.is_file() {
        return Err(CliError::Invalid(format!(
            "Not a file: '{}",
            path.display()
        )));
    }
    Ok(absolute_if_no_parent(path))
}

fn valid_directory(path: impl AsRef<Path>, create_if_not_exist: bool) -> Result<PathBuf, CliError> {
    let path = path.as_ref();
    if path.exists() {
        if path.is_dir() {
            return Ok(path.to_path_buf());
        }
        return Err(CliError::Invalid(format!(
            "Not a valid directory: {}",
            path.display()
        )));
    }
    if create_if_not_exist {
        fs::create_dir_all(path).map_err(|source| CliError::Io {
            path: path.to_path_buf(),
            source,
        })?;
        return Ok(path.to_path_buf());
    }
    Err(CliError::Invalid(format!(
        "No such directory: {}",
        path.display()
    )))
}

fn absolute_if_no_parent(path: &Path) -> PathBuf {
    if path.parent().is_none() {
        path.to_path_buf()
            .canonicalize()
            .unwrap_or_else(|_| path.to_path_buf())
    } else {
        path.to_path_buf()
    }
}

fn split_comma(value: &str) -> impl Iterator<Item = &str> {
    value
        .split(',')
        .map(str::trim)
        .filter(|item| !item.is_empty())
}

fn parse_reporter_sources(value: &str) -> Result<Vec<PrescribingGuidanceSource>, CliError> {
    let mut sources = Vec::new();
    for source in split_comma(value) {
        match source {
            "CPIC" => sources.push(PrescribingGuidanceSource::CpicGuideline),
            "DPWG" => sources.push(PrescribingGuidanceSource::DpwgGuideline),
            "FDA" => {
                sources.push(PrescribingGuidanceSource::FdaLabel);
                sources.push(PrescribingGuidanceSource::FdaAssoc);
            }
            "ClinPGx" | "Unknown" => {
                return Err(CliError::Invalid(format!("Unsupported source: {source}")));
            }
            _ => return Err(CliError::Invalid(format!("Unknown source: {source}"))),
        }
    }
    Ok(sources)
}

/// Java `BaseConfig.getBaseFilename`.
pub fn base_filename(input_file: impl Into<PathBuf>) -> String {
    let input_file = input_file.into();
    base_filename_from_path(&input_file)
}

fn base_filename_from_path(input_file: &Path) -> String {
    let file_name = input_file
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or_default();
    let mut filename = strip_extension(file_name);
    if filename.ends_with(".vcf") {
        filename = strip_extension(&filename);
    }
    for suffix in [
        VCF_PREPROCESSED_SUFFIX,
        MATCHER_SUFFIX,
        PHENOTYPER_SUFFIX,
        REPORTER_SUFFIX,
    ] {
        if filename.ends_with(suffix) {
            filename.truncate(filename.len() - suffix.len());
        }
    }
    if let Some(outside_base) = strip_outside_suffix(&filename) {
        filename = outside_base.to_string();
    }
    filename
}

fn base_dir(input_file: &Path) -> Result<PathBuf, CliError> {
    let abs_path = if input_file.is_absolute() {
        input_file.to_path_buf()
    } else {
        std::env::current_dir()
            .map_err(|source| CliError::Io {
                path: input_file.to_path_buf(),
                source,
            })?
            .join(input_file)
    };
    abs_path.parent().map(Path::to_path_buf).ok_or_else(|| {
        CliError::Invalid(
            "Cannot determine directory to save results to.  Please specify output directory."
                .to_string(),
        )
    })
}

fn strip_extension(value: &str) -> String {
    value
        .rsplit_once('.')
        .map_or_else(|| value.to_string(), |(base, _)| base.to_string())
}

fn strip_outside_suffix(value: &str) -> Option<&str> {
    let (base, suffix) = value.rsplit_once(".outside")?;
    if suffix.chars().all(|ch| ch.is_ascii_digit()) {
        Some(base)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use std::{
        fs,
        path::PathBuf,
        time::{SystemTime, UNIX_EPOCH},
    };

    use super::*;

    #[test]
    fn parses_default_full_run_with_vcf_like_java_pharmcat() {
        let vcf_file = write_temp_file("pharmcat-sample", "##fileformat=VCFv4.3\n");
        let output_dir = unique_temp_path("pharmcat-out");
        let action = parse_pharmcat_args([
            "-vcf",
            vcf_file.to_str().unwrap(),
            "-o",
            output_dir.to_str().unwrap(),
            "-g",
            "vkorc1, cyp2c9",
            "-s",
            "S2,S1,S2",
        ])
        .expect("parse");
        let CliAction::Run(config) = action else {
            panic!("expected run");
        };
        assert!(config.run_matcher);
        assert!(config.run_phenotyper);
        assert!(config.run_reporter);
        assert_eq!(config.matcher_vcf, Some(vcf_file));
        assert_eq!(config.output_dir, Some(output_dir.clone()));
        assert!(output_dir.is_dir());
        assert_eq!(
            config.genes.into_iter().collect::<Vec<_>>(),
            ["CYP2C9", "VKORC1"]
        );
        assert_eq!(config.samples.into_iter().collect::<Vec<_>>(), ["S1", "S2"]);
        assert!(config.reporter_html);
        assert!(!config.reporter_json);
    }

    #[test]
    fn mirrors_java_sample_file_and_duplicate_sample_rules() {
        let vcf_file = write_temp_file("pharmcat-sample", "##fileformat=VCFv4.3\n");
        let sample_file = write_temp_file("pharmcat-samples", "S3\n\n S1 \nS2\nS1\n");
        let action = parse_pharmcat_args([
            "-vcf",
            vcf_file.to_str().unwrap(),
            "-S",
            sample_file.to_str().unwrap(),
        ])
        .expect("parse");
        let CliAction::Run(config) = action else {
            panic!("expected run");
        };
        assert_eq!(
            config.samples.into_iter().collect::<Vec<_>>(),
            ["S1", "S2", "S3"]
        );

        let err = parse_pharmcat_args([
            "-vcf",
            "sample.vcf",
            "-s",
            "S1",
            "-S",
            sample_file.to_str().unwrap(),
        ])
        .expect_err("sample args should conflict");
        assert_eq!(err.to_string(), "Cannot specify both -s and -S");

        let bad_sample_file = write_temp_file("pharmcat-samples-bad", "S1,S2\n");
        let err = parse_pharmcat_args([
            "-vcf",
            "sample.vcf",
            "-S",
            bad_sample_file.to_str().unwrap(),
        ])
        .expect_err("comma in sample file should fail");
        assert!(err.to_string().contains("Please remove comma"));
    }

    #[test]
    fn mirrors_java_stage_input_validation() {
        let err = parse_pharmcat_args(Vec::<String>::new()).expect_err("empty CLI needs VCF");
        assert!(
            err.to_string()
                .contains("No input for Named Allele Matcher")
        );
        assert!(err.to_string().contains("-vcf"));

        let err = parse_pharmcat_args(["-matcher", "-reporter", "-vcf", "sample.vcf"])
            .expect_err("matcher reporter without phenotyper");
        assert_eq!(
            err.to_string(),
            "Cannot run matcher and reporter without also running phenotyper."
        );

        let err = parse_pharmcat_args(["-phenotyper"]).expect_err("phenotyper needs input");
        assert!(err.to_string().contains("No input for Phenotyper"));

        let err = parse_pharmcat_args(["-reporter"]).expect_err("reporter needs input");
        assert!(err.to_string().contains("No input for Reporter"));
    }

    #[test]
    fn parses_reporter_sources_and_output_format_defaults_like_java() {
        let reporter_input = write_temp_file("pharmcat-phenotype", "{}\n");
        let action = parse_pharmcat_args([
            "-reporter",
            "-ri",
            reporter_input.to_str().unwrap(),
            "-rs",
            "CPIC,FDA",
            "-reporterJson",
        ])
        .expect("parse");
        let CliAction::Run(config) = action else {
            panic!("expected run");
        };
        assert_eq!(
            config.reporter_sources,
            Some(vec![
                PrescribingGuidanceSource::CpicGuideline,
                PrescribingGuidanceSource::FdaLabel,
                PrescribingGuidanceSource::FdaAssoc,
            ])
        );
        assert!(config.reporter_json);
        assert!(!config.reporter_html);

        let err = parse_pharmcat_args(["-reporter", "-ri", "x.json", "-rs", "ClinPGx"])
            .expect_err("unsupported source");
        assert_eq!(err.to_string(), "Unsupported source: ClinPGx");
    }

    #[test]
    fn research_mode_matches_java_reporter_output_restriction() {
        let vcf_file = write_temp_file("pharmcat-sample", "##fileformat=VCFv4.3\n");
        let action = parse_pharmcat_args([
            "-vcf",
            vcf_file.to_str().unwrap(),
            "-research",
            "cyp2d6,combination",
            "-reporterJson",
            "-reporterHtml",
        ])
        .expect("parse");
        let CliAction::Run(config) = action else {
            panic!("expected run");
        };
        assert!(config.call_cyp2d6);
        assert!(config.find_combinations);
        assert!(config.reporter_calls_only_tsv);
        assert!(!config.reporter_json);
        assert!(!config.reporter_html);

        let err = parse_pharmcat_args(["-vcf", "sample.vcf", "-research", "bogus"])
            .expect_err("unknown research");
        assert_eq!(err.to_string(), "Unrecognized research option: bogus");
    }

    #[test]
    fn validates_cli_paths_like_pgkb_common_cli_helper() {
        let missing = unique_temp_path("pharmcat-missing-vcf");
        let err = parse_pharmcat_args(["-vcf", missing.to_str().unwrap()])
            .expect_err("missing matcher VCF should fail");
        assert_eq!(
            err.to_string(),
            format!("File '{}' does not exist", missing.display())
        );

        let dir = unique_temp_path("pharmcat-not-file");
        fs::create_dir_all(&dir).expect("create temp dir");
        let err = parse_pharmcat_args(["-vcf", dir.to_str().unwrap()])
            .expect_err("directory as VCF should fail");
        assert_eq!(err.to_string(), format!("Not a file: '{}", dir.display()));

        let file = write_temp_file("pharmcat-not-dir", "");
        let err = parse_pharmcat_args([
            "-vcf",
            file.to_str().unwrap(),
            "-def",
            file.to_str().unwrap(),
        ])
        .expect_err("file as definitions dir should fail");
        assert_eq!(
            err.to_string(),
            format!("Not a valid directory: {}", file.display())
        );

        let missing_def = unique_temp_path("pharmcat-missing-def");
        let err = parse_pharmcat_args([
            "-vcf",
            file.to_str().unwrap(),
            "-def",
            missing_def.to_str().unwrap(),
        ])
        .expect_err("missing definitions dir should fail");
        assert_eq!(
            err.to_string(),
            format!("No such directory: {}", missing_def.display())
        );

        let missing_outside = unique_temp_path("pharmcat-missing-outside");
        let err = parse_pharmcat_args(["-phenotyper", "-po", missing_outside.to_str().unwrap()])
            .expect_err("missing outside call file should fail");
        assert_eq!(
            err.to_string(),
            format!("Not a valid file: '{}", missing_outside.display())
        );
    }

    #[test]
    fn plans_full_pipeline_output_paths_like_java_pipeline() {
        let vcf_file = write_temp_named_file("pharmcat-reference.vcf", "##fileformat=VCFv4.3\n");
        let output_dir = unique_temp_path("pharmcat-plan-out");
        let action = parse_pharmcat_args([
            "-vcf",
            vcf_file.to_str().unwrap(),
            "-o",
            output_dir.to_str().unwrap(),
            "-matcherHtml",
            "-reporterJson",
            "-reporterHtml",
            "-reporterCallsOnlyTsv",
        ])
        .expect("parse");
        let CliAction::Run(config) = action else {
            panic!("expected run");
        };

        let plan = pipeline_output_plan(&config, Some("NA12878"), false).expect("plan");
        assert_eq!(plan.base_dir, output_dir);
        assert_eq!(plan.basename, "pharmcat-reference.NA12878");
        assert_eq!(
            plan.display_name,
            format!(
                "sample NA12878 in {}",
                vcf_file.file_name().unwrap().to_string_lossy()
            )
        );
        assert_eq!(
            plan.matcher_json,
            Some(plan.base_dir.join("pharmcat-reference.NA12878.match.json"))
        );
        assert_eq!(
            plan.matcher_html,
            Some(plan.base_dir.join("pharmcat-reference.NA12878.match.html"))
        );
        assert_eq!(
            plan.matcher_warnings,
            Some(
                plan.base_dir
                    .join("pharmcat-reference.NA12878.match_warnings.txt")
            )
        );
        assert_eq!(
            plan.phenotyper_json,
            Some(
                plan.base_dir
                    .join("pharmcat-reference.NA12878.phenotype.json")
            )
        );
        assert_eq!(
            plan.reporter_html,
            Some(plan.base_dir.join("pharmcat-reference.NA12878.report.html"))
        );
        assert_eq!(
            plan.reporter_json,
            Some(plan.base_dir.join("pharmcat-reference.NA12878.report.json"))
        );
        assert_eq!(
            plan.reporter_calls_only_tsv,
            Some(plan.base_dir.join("pharmcat-reference.NA12878.report.tsv"))
        );
        assert_eq!(
            plan.reporter_title,
            Some("pharmcat-reference.NA12878".to_string())
        );
    }

    #[test]
    fn plans_independent_reporter_and_base_filename_override_like_java_pipeline() {
        let reporter_input = write_temp_named_file("sample.phenotype.json", "{}\n");
        let output_dir = unique_temp_path("pharmcat-reporter-out");
        let action = parse_pharmcat_args([
            "-reporter",
            "-ri",
            reporter_input.to_str().unwrap(),
            "-o",
            output_dir.to_str().unwrap(),
            "-bf",
            "custom",
            "-rt",
            "Custom report",
            "-reporterJson",
        ])
        .expect("parse");
        let CliAction::Run(config) = action else {
            panic!("expected run");
        };

        let plan = pipeline_output_plan(&config, None, true).expect("plan");
        assert_eq!(plan.basename, "custom");
        assert_eq!(plan.display_name, "custom");
        assert_eq!(plan.reporter_title, Some("Custom report".to_string()));
        assert_eq!(plan.matcher_json, None);
        assert_eq!(plan.phenotyper_json, None);
        assert_eq!(plan.reporter_html, None);
        assert_eq!(
            plan.reporter_json,
            Some(output_dir.join("custom.report.json"))
        );
    }

    #[test]
    fn avoids_duplicate_sample_suffix_like_java_pipeline() {
        let vcf_file = write_temp_named_file("NA12878.reference.vcf", "##fileformat=VCFv4.3\n");
        let action = parse_pharmcat_args(["-vcf", vcf_file.to_str().unwrap()]).expect("parse");
        let CliAction::Run(config) = action else {
            panic!("expected run");
        };

        let plan = pipeline_output_plan(&config, Some("NA12878"), false).expect("plan");
        assert_eq!(plan.basename, "NA12878.reference");
        assert_eq!(plan.reporter_title, Some("NA12878.reference".to_string()));
    }

    #[test]
    fn base_filename_matches_java_suffix_cleanup() {
        assert_eq!(base_filename("sample.vcf"), "sample");
        assert_eq!(base_filename("sample.vcf.gz"), "sample");
        assert_eq!(base_filename("sample.preprocessed.vcf"), "sample");
        assert_eq!(base_filename("sample.match.json"), "sample");
        assert_eq!(base_filename("sample.outside.tsv"), "sample");
        assert_eq!(base_filename("sample.outside2.tsv"), "sample");
        assert_eq!(base_filename("sample.phenotype.json"), "sample");
        assert_eq!(base_filename("sample.report.html"), "sample");
    }

    fn write_temp_file(prefix: &str, contents: &str) -> PathBuf {
        let path = unique_temp_path(prefix).with_extension("txt");
        fs::write(&path, contents).expect("write temp file");
        path
    }

    fn write_temp_named_file(filename: &str, contents: &str) -> PathBuf {
        let dir = unique_temp_path("pharmcat-cli-test-dir");
        fs::create_dir_all(&dir).expect("create temp dir");
        let path = dir.join(filename);
        fs::write(&path, contents).expect("write temp file");
        path
    }

    fn unique_temp_path(prefix: &str) -> PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("time")
            .as_nanos();
        std::env::temp_dir().join(format!("{prefix}-{nanos}"))
    }
}
