//! Report and prescribing-guidance data helpers.

use std::{
    cmp::Ordering,
    collections::{BTreeMap, BTreeSet},
    fs, io,
    path::{Path, PathBuf},
};

use serde::{Deserialize, Serialize, ser::SerializeStruct};
use serde_json::Value;

use crate::{
    definition::{DefinitionFile, DefinitionReader, VariantLocus},
    matcher::{GeneCallKind, GeneCallResult, GeneCallWarning, MatchData, compare_haplotype_names},
    phenotype::{
        AnnotatedDiplotype, DiplotypeAnnotationInput, GenePhenotype, OutsideCall, OutsideCallError,
        PhenotypeLookupError, PhenotypeMap,
    },
};

/// Java `DataSource`.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd, Serialize)]
pub enum DataSource {
    /// CPIC source.
    #[serde(rename = "CPIC")]
    Cpic,
    /// DPWG source.
    #[serde(rename = "DPWG")]
    Dpwg,
    /// ClinPGx source.
    #[serde(rename = "CLINPGX")]
    ClinPgx,
    /// FDA source.
    #[serde(rename = "FDA")]
    Fda,
    /// Unknown source.
    #[serde(rename = "UNKNOWN")]
    Unknown,
}

impl Default for DataSource {
    fn default() -> Self {
        Self::Unknown
    }
}

impl DataSource {
    /// Java `DataSource.getPharmgkbName`.
    pub fn pharmgkb_name(self) -> &'static str {
        match self {
            Self::Cpic => "CPIC",
            Self::Dpwg => "DPWG",
            Self::ClinPgx => "ClinPGx",
            Self::Fda => "FDA",
            Self::Unknown => "Unknown",
        }
    }
}

fn data_source_from_definition(definition: &DefinitionFile) -> DataSource {
    definition
        .source
        .as_ref()
        .and_then(Value::as_str)
        .map(|source| match source {
            "CPIC" => DataSource::Cpic,
            "DPWG" => DataSource::Dpwg,
            "CLINPGX" | "ClinPGx" => DataSource::ClinPgx,
            "FDA" => DataSource::Fda,
            _ => DataSource::Unknown,
        })
        .unwrap_or(DataSource::Unknown)
}

/// Java `DrugLink` serialized in `GeneReport.relatedDrugs`.
#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd, Serialize)]
pub struct DrugLink {
    /// Drug display name.
    #[serde(rename = "name")]
    pub name: String,
    /// Drug PharmGKB id.
    #[serde(rename = "id")]
    pub id: String,
}

impl DrugLink {
    fn new(name: impl Into<String>, id: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            id: id.into(),
        }
    }
}

/// Java `PrescribingGuidanceSource`.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd, Serialize)]
pub enum PrescribingGuidanceSource {
    /// CPIC guideline annotation.
    #[serde(rename = "CPIC_GUIDELINE")]
    CpicGuideline,
    /// DPWG guideline annotation.
    #[serde(rename = "DPWG_GUIDELINE")]
    DpwgGuideline,
    /// FDA label annotation.
    #[serde(rename = "FDA_LABEL")]
    FdaLabel,
    /// FDA PGx association.
    #[serde(rename = "FDA_ASSOC")]
    FdaAssoc,
}

impl PrescribingGuidanceSource {
    /// Lists sources in Java enum order.
    pub fn list_values() -> [Self; 4] {
        [
            Self::CpicGuideline,
            Self::DpwgGuideline,
            Self::FdaLabel,
            Self::FdaAssoc,
        ]
    }

    /// Java `getPgkbSource`.
    pub fn pgkb_source(self) -> DataSource {
        match self {
            Self::CpicGuideline => DataSource::Cpic,
            Self::DpwgGuideline => DataSource::Dpwg,
            Self::FdaLabel | Self::FdaAssoc => DataSource::Fda,
        }
    }

    /// Java `getPgkbObjectType`.
    pub fn pgkb_object_type(self) -> &'static str {
        match self {
            Self::CpicGuideline | Self::DpwgGuideline => "Guideline Annotation",
            Self::FdaLabel => "Label Annotation",
            Self::FdaAssoc => "PGx Association",
        }
    }

    /// Java `getDisplayName`.
    pub fn display_name(self) -> &'static str {
        match self {
            Self::CpicGuideline => "CPIC Guideline Annotation",
            Self::DpwgGuideline => "DPWG Guideline Annotation",
            Self::FdaLabel => "FDA Label Annotation",
            Self::FdaAssoc => "FDA PGx Association",
        }
    }

    /// Java `getCodeName`.
    pub fn code_name(self) -> &'static str {
        match self {
            Self::CpicGuideline => "cpic-guideline",
            Self::DpwgGuideline => "dpwg-guideline",
            Self::FdaLabel => "fda-label",
            Self::FdaAssoc => "fda-assoc",
        }
    }

    /// Java `matches`.
    pub fn matches(self, guideline: &DosingGuideline) -> bool {
        guideline.source == self.pgkb_source().pharmgkb_name()
            && guideline.obj_cls == self.pgkb_object_type()
    }
}

/// Java `CallSource` ordering for preserved per-source report genes.
#[derive(Clone, Copy, Debug, Default, Eq, Ord, PartialEq, PartialOrd, Serialize)]
pub enum ReportCallSource {
    /// Java `CallSource.OUTSIDE`.
    #[serde(rename = "OUTSIDE")]
    Outside,
    /// Java `CallSource.MATCHER`.
    #[serde(rename = "MATCHER")]
    Matcher,
    /// Java `CallSource.NONE`.
    #[serde(rename = "NONE")]
    #[default]
    None,
}

/// Java `PgkbGuidelineCollection`.
#[derive(Clone, Debug, PartialEq)]
pub struct PgkbGuidelineCollection {
    version: Option<String>,
    guideline_packages: Vec<GuidelinePackage>,
    guideline_map: BTreeMap<String, BTreeSet<usize>>,
}

impl PgkbGuidelineCollection {
    /// Java `PRESCRIBING_GUIDANCE_FILE_NAME`.
    pub const PRESCRIBING_GUIDANCE_FILE_NAME: &'static str = "prescribing_guidance.json";

    /// Loads prescribing guidance JSON from `path`.
    pub fn from_path(path: &Path) -> Result<Self, GuidanceLoadError> {
        let data = fs::read(path)?;
        let dataset: PrescribingGuidanceDataset = serde_json::from_slice(&data)?;

        let mut guideline_map = BTreeMap::<String, BTreeSet<usize>>::new();
        for (index, package) in dataset.guideline_packages.iter().enumerate() {
            for chemical in &package.guideline.related_chemicals {
                guideline_map
                    .entry(chemical.name.clone())
                    .or_default()
                    .insert(index);
            }
        }

        Ok(Self {
            version: dataset.version,
            guideline_packages: dataset.guideline_packages,
            guideline_map,
        })
    }

    /// Java `getVersion`.
    pub fn version(&self) -> Option<&str> {
        self.version.as_deref()
    }

    /// Java `getGuidelinePackages`.
    pub fn guideline_packages(&self) -> &[GuidelinePackage] {
        &self.guideline_packages
    }

    /// Java `findGuidelinePackages`.
    pub fn find_guideline_packages(
        &self,
        chemical_name: &str,
        source: PrescribingGuidanceSource,
    ) -> Vec<&GuidelinePackage> {
        self.guideline_map
            .get(chemical_name)
            .into_iter()
            .flat_map(|indices| indices.iter())
            .filter_map(|index| self.guideline_packages.get(*index))
            .filter(|package| package.is_data_source_type(source))
            .collect()
    }

    /// Chemical names indexed in Java `getGuidelineMap().keys()` order.
    pub fn chemical_names(&self) -> BTreeSet<String> {
        self.guideline_map.keys().cloned().collect()
    }

    /// Java `getGuidelinesFromSource(PrescribingGuidanceSource)`.
    pub fn guidelines_from_source(
        &self,
        source: PrescribingGuidanceSource,
    ) -> Vec<&GuidelinePackage> {
        self.guideline_packages
            .iter()
            .filter(|package| package.is_data_source_type(source))
            .collect()
    }

    /// Java `getGenesWithRecommendations`.
    pub fn genes_with_recommendations(&self) -> BTreeSet<String> {
        self.guideline_packages
            .iter()
            .flat_map(|package| package.recommendations.iter())
            .flat_map(|recommendation| recommendation.lookup_genes())
            .collect()
    }

    /// Java `getGenesUsedInSource`.
    pub fn genes_used_in_source(&self, source: DataSource) -> BTreeSet<String> {
        self.guideline_packages
            .iter()
            .filter(|package| package.guideline.source.eq(source.pharmgkb_name()))
            .flat_map(GuidelinePackage::genes)
            .collect()
    }
}

/// Java `PrescribingGuidanceDataset`.
#[derive(Clone, Debug, Deserialize, PartialEq)]
pub struct PrescribingGuidanceDataset {
    /// Data version.
    #[serde(default)]
    pub version: Option<String>,
    /// Guideline packages.
    #[serde(rename = "guidelines", default)]
    pub guideline_packages: Vec<GuidelinePackage>,
}

/// Java `GuidelinePackage`.
#[derive(Clone, Debug, Deserialize, PartialEq)]
pub struct GuidelinePackage {
    /// Guideline annotation.
    pub guideline: DosingGuideline,
    /// Recommendation annotations.
    #[serde(default)]
    pub recommendations: Vec<RecommendationAnnotation>,
    /// Citation publications.
    #[serde(default)]
    pub citations: Vec<Publication>,
    /// PharmGKB URL.
    #[serde(default)]
    pub url: Option<String>,
}

impl GuidelinePackage {
    /// Java `getGenes`.
    pub fn genes(&self) -> BTreeSet<String> {
        self.guideline
            .related_genes
            .iter()
            .filter_map(|gene| gene.symbol.clone())
            .collect()
    }

    /// Java `getDrugs`.
    pub fn drugs(&self) -> BTreeSet<String> {
        self.guideline
            .related_chemicals
            .iter()
            .map(|chemical| chemical.name.clone())
            .collect()
    }

    /// Java `isDataSourceType`.
    pub fn is_data_source_type(&self, source: PrescribingGuidanceSource) -> bool {
        source.matches(&self.guideline)
    }
}

/// Java `DosingGuideline`.
#[derive(Clone, Debug, Deserialize, PartialEq)]
#[serde(rename_all = "camelCase")]
pub struct DosingGuideline {
    /// PharmGKB id.
    pub id: String,
    /// Guideline name.
    pub name: String,
    /// PharmGKB object class.
    #[serde(rename = "objCls")]
    pub obj_cls: String,
    /// PharmGKB source.
    pub source: String,
    /// Guideline version.
    #[serde(default)]
    pub version: Option<i64>,
    /// URL.
    #[serde(default)]
    pub url: Option<String>,
    /// Related chemicals.
    #[serde(default)]
    pub related_chemicals: Vec<AccessionObject>,
    /// Related genes.
    #[serde(default)]
    pub related_genes: Vec<AccessionObject>,
    /// Related alleles.
    #[serde(default)]
    pub related_alleles: Vec<AccessionObject>,
    /// Whether this is a recommendation guideline.
    #[serde(default)]
    pub recommendation: bool,
}

/// Java `RecommendationAnnotation`.
#[derive(Clone, Debug, Deserialize, PartialEq)]
#[serde(rename_all = "camelCase")]
pub struct RecommendationAnnotation {
    /// PharmGKB id.
    pub id: String,
    /// Recommendation name.
    pub name: String,
    /// Population.
    #[serde(default)]
    pub population: Option<String>,
    /// Classification.
    #[serde(default)]
    pub classification: Option<OntologyTerm>,
    /// Related chemicals.
    #[serde(default)]
    pub related_chemicals: Vec<AccessionObject>,
    /// Recommendation text.
    #[serde(default)]
    pub text: Option<Markdown>,
    /// Implication strings.
    #[serde(default)]
    pub implications: Vec<String>,
    /// Recommendation lookup keys.
    #[serde(rename = "lookupKey", default)]
    pub lookup_key: Vec<BTreeMap<String, Value>>,
    /// Whether this is dosing information.
    #[serde(default)]
    pub dosing_information: bool,
    /// Whether an alternate drug is available.
    #[serde(default)]
    pub alternate_drug_available: bool,
    /// Whether other prescribing guidance is available.
    #[serde(default)]
    pub other_prescribing_guidance: bool,
}

impl RecommendationAnnotation {
    /// Java `getLookupGenes`.
    pub fn lookup_genes(&self) -> BTreeSet<String> {
        self.lookup_key
            .iter()
            .flat_map(|lookup| lookup.keys().cloned())
            .collect()
    }

    /// Java `matchesGenotype`.
    pub fn matches_genotype(&self, genotype: &RecommendationGenotype) -> bool {
        !self.lookup_key.is_empty()
            && genotype.lookup_keys.iter().any(|key| {
                self.lookup_key
                    .iter()
                    .any(|lookup| map_contains(key, lookup))
            })
    }

    /// Java `matchesDiplotype`.
    pub fn matches_diplotype(&self, genotype: &RecommendationGenotype) -> bool {
        !self.lookup_key.is_empty()
            && self
                .lookup_key
                .iter()
                .any(|lookup| map_contains(&genotype.diplotype_key, lookup))
    }

    /// Java `appliesToDrug`.
    pub fn applies_to_drug(&self, drug_name: &str) -> bool {
        self.related_chemicals
            .iter()
            .any(|chemical| chemical.name.eq_ignore_ascii_case(drug_name))
    }
}

/// Minimal Java `Genotype` lookup-key model for recommendation matching.
#[derive(Clone, Debug, Default, Eq, PartialEq, Serialize)]
pub struct RecommendationGenotype {
    /// Java-style genotype lookup-key combinations.
    #[serde(skip)]
    lookup_keys: Vec<BTreeMap<String, Value>>,
    /// Java-style genotype diplotype key.
    #[serde(skip)]
    diplotype_key: BTreeMap<String, Value>,
    /// Report genes represented by this genotype.
    #[serde(rename = "diplotypes")]
    report_genes: Vec<ReportGene>,
}

impl RecommendationGenotype {
    /// Builds genotype lookup-key combinations like Java `Genotype.addDiplotype`.
    pub fn from_gene_lookup_keys<I, K>(entries: I) -> Self
    where
        I: IntoIterator<Item = (String, K)>,
        K: IntoIterator<Item = String>,
    {
        Self::from_report_genes(entries.into_iter().map(|(gene, lookup_keys)| ReportGene {
            allele_definition_version: None,
            allele_definition_source: DataSource::Unknown,
            phenotype_version: None,
            gene,
            chromosome: None,
            phased: false,
            effectively_phased: false,
            lookup_keys: lookup_keys.into_iter().collect(),
            diplotype_key: serde_json::Map::new(),
            phenotypes: Vec::new(),
            activity_score: None,
            is_activity_score_type: false,
            is_allele_presence_type: false,
            source_diplotype: None,
            source_diplotypes: Vec::new(),
            matcher_component_haplotypes: BTreeSet::new(),
            matcher_component_diplotypes: Vec::new(),
            matcher_homozygous_component_haplotypes: BTreeSet::new(),
            recommendation_diplotypes: Vec::new(),
            allele_function_map: BTreeMap::new(),
            related_drugs: BTreeSet::new(),
            match_score: None,
            outside_call: false,
            call_source: ReportCallSource::None,
            guidance_source: None,
            variant_reports: Vec::new(),
            variant_of_interest_reports: Vec::new(),
            uncalled_haplotypes: BTreeSet::new(),
            has_undocumented_variations: false,
            treat_undocumented_variations_as_reference: false,
            messages: BTreeSet::new(),
            outside_phenotype_mismatch: None,
            outside_activity_score_mismatch: None,
        }))
    }

    /// Builds genotype lookup-key and diplotype-key combinations from report genes.
    pub fn from_report_genes<I>(entries: I) -> Self
    where
        I: IntoIterator<Item = ReportGene>,
    {
        let mut lookup_keys: Option<Vec<BTreeMap<String, Value>>> = None;
        let mut diplotype_key = BTreeMap::new();
        let mut entries = entries.into_iter().collect::<Vec<_>>();
        entries.sort_by(report_gene_source_cmp_like_java);
        let mut report_genes = Vec::new();

        for report_gene in entries {
            let gene = report_gene.gene.clone();
            diplotype_key.insert(
                gene.clone(),
                Value::Object(report_gene.diplotype_key.clone()),
            );
            let gene_lookup_keys = report_gene.lookup_keys.clone();
            lookup_keys = Some(match lookup_keys {
                None => gene_lookup_keys
                    .into_iter()
                    .map(|lookup_key| {
                        [(gene.clone(), Value::String(lookup_key))]
                            .into_iter()
                            .collect()
                    })
                    .collect(),
                Some(existing) => gene_lookup_keys
                    .into_iter()
                    .flat_map(|lookup_key| {
                        existing.iter().map({
                            let gene = gene.clone();
                            move |original| {
                                let mut map = original.clone();
                                map.insert(gene.clone(), Value::String(lookup_key.clone()));
                                map
                            }
                        })
                    })
                    .collect(),
            });
            report_genes.push(report_gene);
        }

        Self {
            lookup_keys: lookup_keys.unwrap_or_default(),
            diplotype_key,
            report_genes,
        }
    }

    /// Returns Java-style genotype lookup-key combinations.
    pub fn lookup_keys(&self) -> &[BTreeMap<String, Value>] {
        &self.lookup_keys
    }

    /// Returns Java-style genotype diplotype key.
    pub fn diplotype_key(&self) -> &BTreeMap<String, Value> {
        &self.diplotype_key
    }

    /// Returns the report genes represented by this genotype.
    pub fn report_genes(&self) -> &[ReportGene] {
        &self.report_genes
    }

    fn uses_activity_score(&self) -> bool {
        self.report_genes
            .iter()
            .any(|report_gene| report_gene.is_activity_score_type)
    }
}

/// Minimal gene-report input used to build recommendation genotypes.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct ReportGene {
    /// Allele definition version used for matcher-backed calls.
    pub allele_definition_version: Option<String>,
    /// Allele definition source used for matcher-backed calls.
    pub allele_definition_source: DataSource,
    /// Phenotype data version used for this gene.
    pub phenotype_version: Option<String>,
    /// Gene symbol.
    pub gene: String,
    /// Chromosome this gene appears on.
    pub chromosome: Option<String>,
    /// Whether the matcher call was entirely phased.
    pub phased: bool,
    /// Whether Java treats the matcher call as effectively phased.
    pub effectively_phased: bool,
    /// Recommendation lookup keys for this gene.
    pub lookup_keys: Vec<String>,
    /// Diplotype key for exact diplotype matching.
    pub diplotype_key: serde_json::Map<String, Value>,
    /// Phenotypes assigned to this gene report.
    pub phenotypes: Vec<String>,
    /// Activity score assigned to this gene report.
    pub activity_score: Option<String>,
    /// Whether this gene uses activity-score lookup.
    pub is_activity_score_type: bool,
    /// Whether this gene uses allele-presence lookup.
    pub is_allele_presence_type: bool,
    /// Source diplotype label for calls-only TSV output.
    pub source_diplotype: Option<String>,
    /// Java-style source diplotypes based on matcher or outside-call input.
    pub source_diplotypes: Vec<ReportDiplotype>,
    /// Lowest-function matcher component haplotypes, matching Java `matcherComponentHaplotypes` names for helper parity.
    pub matcher_component_haplotypes: BTreeSet<String>,
    /// Java `matcherComponentHaplotypes` diplotypes for report JSON.
    pub matcher_component_diplotypes: Vec<ReportDiplotype>,
    /// Lowest-function matcher component haplotypes that Java marks homozygous.
    pub matcher_homozygous_component_haplotypes: BTreeSet<String>,
    /// Java-style recommendation diplotypes used to match prescribing guidance.
    pub recommendation_diplotypes: Vec<ReportDiplotype>,
    /// Java `functionMap` entries used by HTML allele-function helpers.
    pub allele_function_map: BTreeMap<String, String>,
    /// Drugs linked back from guideline reports, matching Java `GeneReport.relatedDrugs`.
    pub related_drugs: BTreeSet<DrugLink>,
    /// Match score for calls-only TSV output.
    pub match_score: Option<String>,
    /// Whether this report came from an outside call.
    pub outside_call: bool,
    /// Java `CallSource` for source ordering; skipped because Java JSON exposes `outsideCall`.
    pub call_source: ReportCallSource,
    /// Prescribing-guidance owner for source-level summary ordering, when known.
    pub guidance_source: Option<PrescribingGuidanceSource>,
    /// Variant reports used by message rules and output.
    pub variant_reports: Vec<VariantReport>,
    /// Variant-of-interest reports used by Java report-as-genotype message rules.
    pub variant_of_interest_reports: Vec<VariantReport>,
    /// Haplotypes the matcher could not evaluate.
    pub uncalled_haplotypes: BTreeSet<String>,
    /// Whether the matcher found definition positions with undocumented variations.
    pub has_undocumented_variations: bool,
    /// Whether undocumented variation calls were replaced with reference during matching.
    pub treat_undocumented_variations_as_reference: bool,
    /// Java gene-report messages.
    pub messages: BTreeSet<MessageAnnotation>,
    /// Expected phenotype when an outside phenotype conflicts with annotation data.
    pub outside_phenotype_mismatch: Option<String>,
    /// Expected activity score when an outside score conflicts with annotation data.
    pub outside_activity_score_mismatch: Option<String>,
}

impl Serialize for ReportGene {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        let mut state = serializer.serialize_struct("ReportGene", 19)?;
        state.serialize_field("alleleDefinitionVersion", &self.allele_definition_version)?;
        state.serialize_field("alleleDefinitionSource", &self.allele_definition_source)?;
        state.serialize_field("phenotypeVersion", &self.phenotype_version)?;
        state.serialize_field("geneSymbol", &self.gene)?;
        state.serialize_field("chr", &self.chromosome)?;
        state.serialize_field("phased", &self.phased)?;
        state.serialize_field("effectivelyPhased", &self.effectively_phased)?;
        state.serialize_field("callSource", &self.call_source)?;
        state.serialize_field("uncalledHaplotypes", &self.uncalled_haplotypes)?;
        state.serialize_field("messages", &self.messages)?;
        state.serialize_field("relatedDrugs", &self.related_drugs)?;
        state.serialize_field("sourceDiplotypes", &self.source_diplotypes)?;
        state.serialize_field(
            "matcherComponentHaplotypes",
            &self.matcher_component_diplotypes,
        )?;
        state.serialize_field(
            "matcherHomozygousComponentHaplotypes",
            &self.matcher_homozygous_component_haplotypes,
        )?;
        state.serialize_field("recommendationDiplotypes", &self.recommendation_diplotypes)?;
        state.serialize_field("variants", &self.variant_reports)?;
        state.serialize_field("variantsOfInterest", &self.variant_of_interest_reports)?;
        state.serialize_field(
            "hasUndocumentedVariations",
            &self.has_undocumented_variations,
        )?;
        state.serialize_field(
            "treatUndocumentedVariationsAsReference",
            &self.treat_undocumented_variations_as_reference,
        )?;
        state.end()
    }
}

/// Error converting a Java-style outside call to a report gene.
#[derive(Debug)]
pub enum ReportGeneFromOutsideCallError {
    /// Outside-call diplotype parsing or validation failed.
    OutsideCall(OutsideCallError),
    /// Phenotype annotation failed.
    Phenotype(PhenotypeLookupError),
}

impl std::fmt::Display for ReportGeneFromOutsideCallError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::OutsideCall(error) => write!(f, "{error}"),
            Self::Phenotype(error) => write!(f, "{error}"),
        }
    }
}

impl std::error::Error for ReportGeneFromOutsideCallError {}

impl From<OutsideCallError> for ReportGeneFromOutsideCallError {
    fn from(error: OutsideCallError) -> Self {
        Self::OutsideCall(error)
    }
}

impl From<PhenotypeLookupError> for ReportGeneFromOutsideCallError {
    fn from(error: PhenotypeLookupError) -> Self {
        Self::Phenotype(error)
    }
}

/// Error converting a Java matcher result to a report gene.
#[derive(Debug)]
pub enum ReportGeneFromStandardCallError {
    /// Phenotype annotation failed.
    Phenotype(PhenotypeLookupError),
    /// SLCO1B1 custom fallback failed.
    Slco1b1(Slco1b1CustomCallError),
}

impl std::fmt::Display for ReportGeneFromStandardCallError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Phenotype(error) => write!(f, "{error}"),
            Self::Slco1b1(error) => write!(f, "{error}"),
        }
    }
}

impl std::error::Error for ReportGeneFromStandardCallError {}

impl From<PhenotypeLookupError> for ReportGeneFromStandardCallError {
    fn from(error: PhenotypeLookupError) -> Self {
        Self::Phenotype(error)
    }
}

impl From<Slco1b1CustomCallError> for ReportGeneFromStandardCallError {
    fn from(error: Slco1b1CustomCallError) -> Self {
        Self::Slco1b1(error)
    }
}

/// Error building a report context from Java matcher results.
#[derive(Debug)]
pub enum ReportContextFromMatcherError {
    /// A matcher result referenced a gene without a loaded definition.
    MissingDefinition(String),
    /// A matcher result could not be converted to a report gene.
    ReportGene(ReportGeneFromStandardCallError),
}

impl std::fmt::Display for ReportContextFromMatcherError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::MissingDefinition(gene) => {
                write!(f, "No allele definition loaded for gene {gene}")
            }
            Self::ReportGene(error) => write!(f, "{error}"),
        }
    }
}

impl std::error::Error for ReportContextFromMatcherError {}

impl From<ReportGeneFromStandardCallError> for ReportContextFromMatcherError {
    fn from(error: ReportGeneFromStandardCallError) -> Self {
        Self::ReportGene(error)
    }
}

/// Error applying Java's SLCO1B1 rs4149056 fallback caller.
#[derive(Debug, Eq, PartialEq)]
pub enum Slco1b1CustomCallError {
    /// More than one rs4149056 variant report was attached to the gene.
    MultipleRs4149056Reports,
    /// Phenotype annotation failed for the inferred diplotype.
    Phenotype(PhenotypeLookupError),
}

impl std::fmt::Display for Slco1b1CustomCallError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::MultipleRs4149056Reports => {
                write!(f, "More than one report found for rs4149056")
            }
            Self::Phenotype(error) => write!(f, "{error}"),
        }
    }
}

impl std::error::Error for Slco1b1CustomCallError {}

impl From<PhenotypeLookupError> for Slco1b1CustomCallError {
    fn from(error: PhenotypeLookupError) -> Self {
        Self::Phenotype(error)
    }
}

impl ReportGene {
    /// Creates a report gene with lookup keys.
    pub fn new(gene: impl Into<String>, lookup_keys: impl IntoIterator<Item = String>) -> Self {
        let gene = gene.into();
        let lookup_keys = lookup_keys.into_iter().collect::<Vec<_>>();
        Self {
            allele_definition_version: None,
            allele_definition_source: DataSource::Unknown,
            phenotype_version: None,
            is_allele_presence_type: is_allele_presence_gene(&gene),
            gene,
            chromosome: None,
            phased: false,
            effectively_phased: false,
            phenotypes: lookup_keys.clone(),
            lookup_keys,
            diplotype_key: serde_json::Map::new(),
            activity_score: None,
            is_activity_score_type: false,
            source_diplotype: None,
            source_diplotypes: Vec::new(),
            matcher_component_haplotypes: BTreeSet::new(),
            matcher_component_diplotypes: Vec::new(),
            matcher_homozygous_component_haplotypes: BTreeSet::new(),
            recommendation_diplotypes: Vec::new(),
            allele_function_map: BTreeMap::new(),
            related_drugs: BTreeSet::new(),
            match_score: None,
            outside_call: false,
            call_source: ReportCallSource::None,
            guidance_source: None,
            variant_reports: Vec::new(),
            variant_of_interest_reports: Vec::new(),
            uncalled_haplotypes: BTreeSet::new(),
            has_undocumented_variations: false,
            treat_undocumented_variations_as_reference: false,
            messages: BTreeSet::new(),
            outside_phenotype_mismatch: None,
            outside_activity_score_mismatch: None,
        }
    }

    /// Creates a report gene with Java `Diplotype.getDiplotypeKey`-style allele counts.
    pub fn with_diplotype_counts(
        gene: impl Into<String>,
        lookup_keys: impl IntoIterator<Item = String>,
        diplotype_key: impl IntoIterator<Item = (String, i32)>,
    ) -> Self {
        let gene = gene.into();
        let lookup_keys = lookup_keys.into_iter().collect::<Vec<_>>();
        Self {
            allele_definition_version: None,
            allele_definition_source: DataSource::Unknown,
            phenotype_version: None,
            is_allele_presence_type: is_allele_presence_gene(&gene),
            gene,
            chromosome: None,
            phased: false,
            effectively_phased: false,
            phenotypes: lookup_keys.clone(),
            lookup_keys,
            diplotype_key: diplotype_key
                .into_iter()
                .map(|(allele, count)| {
                    (
                        allele,
                        Value::Number(
                            serde_json::Number::from_f64(f64::from(count))
                                .expect("integer diplotype count converts to finite JSON number"),
                        ),
                    )
                })
                .collect(),
            activity_score: None,
            is_activity_score_type: false,
            source_diplotype: None,
            source_diplotypes: Vec::new(),
            matcher_component_haplotypes: BTreeSet::new(),
            matcher_component_diplotypes: Vec::new(),
            matcher_homozygous_component_haplotypes: BTreeSet::new(),
            recommendation_diplotypes: Vec::new(),
            allele_function_map: BTreeMap::new(),
            related_drugs: BTreeSet::new(),
            match_score: None,
            outside_call: false,
            call_source: ReportCallSource::None,
            guidance_source: None,
            variant_reports: Vec::new(),
            variant_of_interest_reports: Vec::new(),
            uncalled_haplotypes: BTreeSet::new(),
            has_undocumented_variations: false,
            treat_undocumented_variations_as_reference: false,
            messages: BTreeSet::new(),
            outside_phenotype_mismatch: None,
            outside_activity_score_mismatch: None,
        }
    }

    /// Creates a Java-style unknown/no-call report gene.
    pub fn unknown(
        gene: impl Into<String>,
        phenotype: Option<&GenePhenotype>,
    ) -> Result<Self, PhenotypeLookupError> {
        let gene = gene.into();
        let allele2 = (!is_single_ploidy_report_gene(&gene)).then_some("Unknown");
        let annotated =
            DiplotypeAnnotationInput::from_alleles(gene, "Unknown", allele2).annotate(phenotype)?;
        Ok(Self::from_annotated_diplotype(annotated))
    }

    /// Creates a report gene from a Java `OutsideCall`.
    pub fn from_outside_call(
        call: &OutsideCall,
        phenotype: Option<&GenePhenotype>,
    ) -> Result<Self, ReportGeneFromOutsideCallError> {
        let phenotype = phenotype.filter(|phenotype| phenotype.gene == call.gene);
        let mut report_gene = if call.is_no_call() {
            Self::unknown(&call.gene, phenotype)?
        } else {
            let annotated = call.to_annotation_input()?.annotate(phenotype)?;
            Self::from_annotated_diplotype_with_phenotype(annotated, phenotype)
        };
        report_gene.outside_call = true;
        report_gene.call_source = ReportCallSource::Outside;
        Ok(report_gene)
    }

    /// Creates a report gene from a standard, non-lowest-function matcher result.
    pub fn from_standard_gene_call_result(
        result: &GeneCallResult,
        phenotype: Option<&GenePhenotype>,
    ) -> Result<Option<Self>, PhenotypeLookupError> {
        if let Some(phenotype) = phenotype
            && result.gene != phenotype.gene
        {
            return Ok(None);
        }
        if matches!(result.kind, GeneCallKind::NoCall) {
            let mut report_gene = Self::unknown(&result.gene, phenotype)?;
            report_gene.call_source = ReportCallSource::Matcher;
            report_gene.messages = gene_call_warning_messages(result);
            return Ok(Some(report_gene));
        }
        if let GeneCallKind::Diplotypes(diplotypes) = &result.kind {
            let mut report_diplotypes = Vec::new();
            for diplotype in diplotypes {
                let (allele1, allele2) = split_report_diplotype_label(&diplotype.name);
                let Some(allele1) = allele1 else {
                    continue;
                };
                let annotated =
                    DiplotypeAnnotationInput::from_alleles(&result.gene, allele1, allele2)
                        .annotate(phenotype)?;
                let mut report_diplotype = ReportDiplotype::from_annotated(&annotated, phenotype);
                report_diplotype.match_score = Some(diplotype.score.to_string());
                report_diplotype.combination = diplotype.name.split('/').any(is_combination_label);
                report_diplotypes.push(report_diplotype);
            }

            if report_diplotypes.is_empty() {
                return Self::unknown(&result.gene, phenotype).map(Some);
            }

            let first = report_diplotypes
                .first()
                .expect("checked non-empty report diplotypes")
                .clone();
            let mut lookup_keys = Vec::new();
            let mut phenotypes = Vec::new();
            for diplotype in &report_diplotypes {
                append_unique(&mut lookup_keys, diplotype.lookup_keys.clone());
                append_unique(&mut phenotypes, diplotype.phenotypes.clone());
            }
            let mut report_gene = Self {
                allele_definition_version: None,
                allele_definition_source: DataSource::Unknown,
                phenotype_version: phenotype.and_then(|phenotype| phenotype.version.clone()),
                is_allele_presence_type: is_allele_presence_gene(&result.gene),
                gene: result.gene.clone(),
                chromosome: None,
                phased: false,
                effectively_phased: false,
                lookup_keys,
                diplotype_key: first.diplotype_key.clone(),
                phenotypes,
                activity_score: first.activity_score.clone(),
                is_activity_score_type: phenotype.is_some_and(GenePhenotype::is_activity_gene),
                source_diplotype: Some(first.label.clone()),
                source_diplotypes: report_diplotypes.clone(),
                matcher_component_haplotypes: BTreeSet::new(),
                matcher_component_diplotypes: Vec::new(),
                matcher_homozygous_component_haplotypes: BTreeSet::new(),
                recommendation_diplotypes: report_diplotypes,
                allele_function_map: allele_function_map_from_phenotype(phenotype),
                related_drugs: BTreeSet::new(),
                match_score: first.match_score.clone(),
                outside_call: false,
                call_source: ReportCallSource::Matcher,
                guidance_source: None,
                variant_reports: Vec::new(),
                variant_of_interest_reports: Vec::new(),
                uncalled_haplotypes: BTreeSet::new(),
                has_undocumented_variations: false,
                treat_undocumented_variations_as_reference: false,
                messages: gene_call_warning_messages(result),
                outside_phenotype_mismatch: None,
                outside_activity_score_mismatch: None,
            };
            sort_report_diplotypes(&mut report_gene.source_diplotypes);
            sort_report_diplotypes(&mut report_gene.recommendation_diplotypes);
            return Ok(Some(report_gene));
        }
        Ok(None)
    }

    /// Creates a report gene from a standard matcher result and attaches Java-style variant reports.
    pub fn from_standard_gene_call_result_with_definition(
        result: &GeneCallResult,
        phenotype: Option<&GenePhenotype>,
        definition: &DefinitionFile,
    ) -> Result<Option<Self>, ReportGeneFromStandardCallError> {
        Self::from_standard_gene_call_result_with_definition_and_messages(
            result, phenotype, definition, None,
        )
    }

    /// Creates a report gene from a Java matcher result and definition, selecting special
    /// lowest-function handling for genes that Java does not route through the standard caller.
    pub fn from_gene_call_result_with_definition(
        result: &GeneCallResult,
        phenotype: Option<&GenePhenotype>,
        definition: &DefinitionFile,
    ) -> Result<Option<Self>, ReportGeneFromStandardCallError> {
        Self::from_gene_call_result_with_definition_and_messages(
            result, phenotype, definition, None,
        )
    }

    /// Creates a report gene from a Java matcher result and definition, applying Java reporter
    /// messages when a catalog is available for the standard matcher path.
    pub fn from_gene_call_result_with_definition_and_messages(
        result: &GeneCallResult,
        phenotype: Option<&GenePhenotype>,
        definition: &DefinitionFile,
        catalog: Option<&MessageCatalog>,
    ) -> Result<Option<Self>, ReportGeneFromStandardCallError> {
        match result.gene.as_str() {
            "DPYD" => phenotype
                .map(|phenotype| {
                    Self::from_dpyd_gene_call_result_with_definition(result, phenotype, definition)
                })
                .transpose()
                .map(Option::flatten)
                .map_err(ReportGeneFromStandardCallError::from),
            "RYR1" => phenotype
                .map(|phenotype| {
                    Self::from_ryr1_gene_call_result_with_definition(result, phenotype, definition)
                })
                .transpose()
                .map(Option::flatten)
                .map_err(ReportGeneFromStandardCallError::from),
            _ => Self::from_standard_gene_call_result_with_definition_and_messages(
                result, phenotype, definition, catalog,
            ),
        }
    }

    /// Creates a report gene from a standard matcher result, attaches Java-style variant reports,
    /// and applies static Java matcher messages when a message catalog is available.
    pub fn from_standard_gene_call_result_with_definition_and_messages(
        result: &GeneCallResult,
        phenotype: Option<&GenePhenotype>,
        definition: &DefinitionFile,
        catalog: Option<&MessageCatalog>,
    ) -> Result<Option<Self>, ReportGeneFromStandardCallError> {
        let Some(mut report_gene) = Self::from_standard_gene_call_result(result, phenotype)? else {
            return Ok(None);
        };
        report_gene
            .messages
            .extend(gene_call_warning_messages(result));
        report_gene.attach_matcher_variant_reports(result, definition);
        report_gene.apply_definition_haplotype_metadata(definition);
        report_gene.apply_definition_report_metadata(definition, phenotype);
        report_gene.add_reference_allele_message(result, definition);
        if let Some(catalog) = catalog {
            report_gene.apply_matcher_static_messages(result, catalog);
            report_gene.apply_matching_gene_messages(catalog);
        }
        report_gene.apply_slco1b1_custom_recommendation(phenotype)?;
        Ok(Some(report_gene))
    }

    /// Creates a report gene from a phenotype-side annotated diplotype.
    pub fn from_annotated_diplotype(diplotype: AnnotatedDiplotype) -> Self {
        Self::from_annotated_diplotype_with_phenotype(diplotype, None)
    }

    /// Creates a report gene from a phenotype-side annotated diplotype with Java haplotype annotation context.
    pub fn from_annotated_diplotype_with_phenotype(
        diplotype: AnnotatedDiplotype,
        phenotype: Option<&GenePhenotype>,
    ) -> Self {
        let diplotype_key =
            diplotype_key_from_alleles(diplotype.allele1.as_deref(), diplotype.allele2.as_deref());
        let source_diplotype = source_diplotype_label(&diplotype);
        let report_diplotype = ReportDiplotype::from_annotated(&diplotype, phenotype);
        let lookup_keys = if diplotype.lookup_keys.is_empty() {
            diplotype.phenotypes.clone()
        } else {
            diplotype.lookup_keys.clone()
        };
        let mut report_gene = Self {
            allele_definition_version: None,
            allele_definition_source: DataSource::Unknown,
            phenotype_version: phenotype.and_then(|phenotype| phenotype.version.clone()),
            is_allele_presence_type: is_allele_presence_gene(&diplotype.gene),
            gene: diplotype.gene,
            chromosome: None,
            phased: false,
            effectively_phased: false,
            lookup_keys,
            diplotype_key,
            phenotypes: diplotype.phenotypes,
            activity_score: diplotype.activity_score,
            is_activity_score_type: diplotype.is_activity_score_type,
            source_diplotype,
            source_diplotypes: vec![report_diplotype.clone()],
            matcher_component_haplotypes: BTreeSet::new(),
            matcher_component_diplotypes: Vec::new(),
            matcher_homozygous_component_haplotypes: BTreeSet::new(),
            recommendation_diplotypes: vec![report_diplotype],
            allele_function_map: allele_function_map_from_phenotype(phenotype),
            related_drugs: BTreeSet::new(),
            match_score: None,
            outside_call: diplotype.outside_phenotype || diplotype.outside_activity_score,
            call_source: if diplotype.outside_phenotype || diplotype.outside_activity_score {
                ReportCallSource::Outside
            } else {
                ReportCallSource::None
            },
            guidance_source: None,
            variant_reports: Vec::new(),
            variant_of_interest_reports: Vec::new(),
            uncalled_haplotypes: BTreeSet::new(),
            has_undocumented_variations: false,
            treat_undocumented_variations_as_reference: false,
            messages: BTreeSet::new(),
            outside_phenotype_mismatch: diplotype.outside_phenotype_mismatch,
            outside_activity_score_mismatch: diplotype.outside_activity_score_mismatch,
        };
        sort_report_diplotypes(&mut report_gene.source_diplotypes);
        sort_report_diplotypes(&mut report_gene.recommendation_diplotypes);
        report_gene
    }

    /// Creates a report gene from an RYR1 lowest-function matcher result.
    pub fn from_ryr1_gene_call_result(
        result: &GeneCallResult,
        phenotype: &GenePhenotype,
    ) -> Result<Option<Self>, PhenotypeLookupError> {
        Self::from_lowest_function_gene_call_result(result, phenotype)
    }

    /// Creates a report gene from an RYR1 lowest-function matcher result with definition metadata.
    pub fn from_ryr1_gene_call_result_with_definition(
        result: &GeneCallResult,
        phenotype: &GenePhenotype,
        definition: &DefinitionFile,
    ) -> Result<Option<Self>, PhenotypeLookupError> {
        Self::from_lowest_function_gene_call_result_with_definition(result, phenotype, definition)
    }

    /// Creates a report gene from a DPYD lowest-function matcher result.
    pub fn from_dpyd_gene_call_result(
        result: &GeneCallResult,
        phenotype: &GenePhenotype,
    ) -> Result<Option<Self>, PhenotypeLookupError> {
        Self::from_lowest_function_gene_call_result(result, phenotype)
    }

    /// Creates a report gene from a DPYD lowest-function matcher result with definition metadata.
    pub fn from_dpyd_gene_call_result_with_definition(
        result: &GeneCallResult,
        phenotype: &GenePhenotype,
        definition: &DefinitionFile,
    ) -> Result<Option<Self>, PhenotypeLookupError> {
        Self::from_lowest_function_gene_call_result_with_definition(result, phenotype, definition)
    }

    /// Creates a report gene from a lowest-function matcher result.
    pub fn from_lowest_function_gene_call_result(
        result: &GeneCallResult,
        phenotype: &GenePhenotype,
    ) -> Result<Option<Self>, PhenotypeLookupError> {
        if result.gene != phenotype.gene {
            return Ok(None);
        }
        if matches!(result.kind, GeneCallKind::NoCall) {
            let mut report_gene = Self::unknown(&result.gene, Some(phenotype))?;
            report_gene.call_source = ReportCallSource::Matcher;
            report_gene.messages = gene_call_warning_messages(result);
            return Ok(Some(report_gene));
        }

        let input = lowest_function_annotation_input(result, phenotype);
        input
            .map(|input| {
                let annotated = input.annotate(Some(phenotype))?;
                let mut recommendation_diplotype =
                    ReportDiplotype::from_annotated(&annotated, Some(phenotype));
                let mut report_gene =
                    Self::from_annotated_diplotype_with_phenotype(annotated, Some(phenotype));
                report_gene.call_source = ReportCallSource::Matcher;
                report_gene.messages = gene_call_warning_messages(result);
                report_gene.source_diplotypes = source_report_diplotypes(result, phenotype);
                report_gene.matcher_component_haplotypes = matcher_component_haplotypes(result);
                report_gene.matcher_component_diplotypes =
                    matcher_component_report_diplotypes(result, phenotype);
                report_gene.matcher_homozygous_component_haplotypes =
                    matcher_homozygous_component_haplotypes(result);
                recommendation_diplotype.inferred =
                    inferred_lowest_function_diplotype(result, &recommendation_diplotype);
                recommendation_diplotype.inferred_source_diplotypes =
                    report_gene.source_diplotypes.clone();
                report_gene.recommendation_diplotypes = vec![recommendation_diplotype];
                sort_report_diplotypes(&mut report_gene.source_diplotypes);
                sort_report_diplotypes(&mut report_gene.matcher_component_diplotypes);
                sort_report_diplotypes(&mut report_gene.recommendation_diplotypes);
                Ok(report_gene)
            })
            .transpose()
    }

    /// Creates a report gene from a lowest-function matcher result with definition metadata.
    pub fn from_lowest_function_gene_call_result_with_definition(
        result: &GeneCallResult,
        phenotype: &GenePhenotype,
        definition: &DefinitionFile,
    ) -> Result<Option<Self>, PhenotypeLookupError> {
        let Some(mut report_gene) = Self::from_lowest_function_gene_call_result(result, phenotype)?
        else {
            return Ok(None);
        };
        report_gene.attach_matcher_variant_reports(result, definition);
        report_gene.apply_definition_haplotype_metadata(definition);
        report_gene.apply_definition_report_metadata(definition, Some(phenotype));
        report_gene.add_reference_allele_message(result, definition);
        Ok(Some(report_gene))
    }

    /// Sets phenotypes assigned to this gene.
    pub fn with_phenotypes(mut self, phenotypes: impl IntoIterator<Item = String>) -> Self {
        self.phenotypes = phenotypes.into_iter().collect();
        self
    }

    /// Sets activity score assigned to this gene and marks it as activity-score based.
    pub fn with_activity_score(mut self, activity_score: impl Into<String>) -> Self {
        self.activity_score = Some(activity_score.into());
        self.is_activity_score_type = true;
        self
    }

    /// Overrides activity-score lookup type.
    pub fn with_activity_score_type(mut self, is_activity_score_type: bool) -> Self {
        self.is_activity_score_type = is_activity_score_type;
        self
    }

    /// Overrides allele-presence lookup type.
    pub fn with_allele_presence_type(mut self, is_allele_presence_type: bool) -> Self {
        self.is_allele_presence_type = is_allele_presence_type;
        self
    }

    /// Sets source diplotype label for calls-only TSV output.
    pub fn with_source_diplotype(mut self, source_diplotype: impl Into<String>) -> Self {
        self.source_diplotype = Some(source_diplotype.into());
        self
    }

    /// Sets match score for calls-only TSV output.
    pub fn with_match_score(mut self, match_score: impl Into<String>) -> Self {
        self.match_score = Some(match_score.into());
        self
    }

    /// Sets outside-call flag for calls-only TSV output.
    pub fn with_outside_call(mut self, outside_call: bool) -> Self {
        self.outside_call = outside_call;
        self.call_source = if outside_call {
            ReportCallSource::Outside
        } else {
            ReportCallSource::None
        };
        self
    }

    /// Sets Java `CallSource` metadata used by source-level report ordering.
    pub fn with_call_source(mut self, call_source: ReportCallSource) -> Self {
        self.call_source = call_source;
        self.outside_call = call_source == ReportCallSource::Outside;
        self
    }

    /// Sets prescribing-guidance ownership metadata used by source-level report ordering.
    pub fn with_guidance_source(mut self, guidance_source: PrescribingGuidanceSource) -> Self {
        self.guidance_source = Some(guidance_source);
        self
    }

    /// Adds variant reports used by Java message/highlight rules.
    pub fn with_variant_reports(
        mut self,
        variants: impl IntoIterator<Item = VariantReport>,
    ) -> Self {
        self.variant_reports = variants.into_iter().collect();
        self
    }

    /// Adds Java variants-of-interest used by report-as-genotype message rules.
    pub fn with_variant_of_interest_reports(
        mut self,
        variants: impl IntoIterator<Item = VariantReport>,
    ) -> Self {
        self.variant_of_interest_reports = variants.into_iter().collect();
        self
    }

    /// Adds Java gene-report messages.
    pub fn with_messages(mut self, messages: impl IntoIterator<Item = MessageAnnotation>) -> Self {
        self.messages = messages.into_iter().collect();
        self
    }

    fn add_related_drug(&mut self, drug: DrugLink) {
        self.related_drugs.insert(drug);
    }

    /// Attaches Java `VariantReportFactory`-style variant reports from matcher data.
    pub fn attach_matcher_variant_reports(
        &mut self,
        result: &GeneCallResult,
        definition: &DefinitionFile,
    ) {
        self.chromosome = Some(definition.chromosome.clone());
        self.phased = result.match_data.phased;
        self.effectively_phased = result.match_data.effectively_phased;
        self.uncalled_haplotypes =
            uncalled_haplotypes_for_match_data(&result.match_data, definition);

        let mut variants = Vec::new();
        for locus in &result.match_data.positions {
            let mut variant = VariantReport::from_matcher_locus(
                &result.gene,
                locus,
                result.match_data.sample_allele_at_position(locus.position),
                definition,
            );
            if result
                .match_data
                .positions_with_undocumented_variations
                .contains(locus)
            {
                variant.has_undocumented_variations = true;
            }
            variants.push(variant);
        }
        for locus in &result.match_data.missing_positions {
            variants.push(VariantReport::from_matcher_locus(
                &result.gene,
                locus,
                None,
                definition,
            ));
        }
        variants.sort();
        self.variant_reports = variants;
        self.has_undocumented_variations = !result
            .match_data
            .positions_with_undocumented_variations
            .is_empty();
        self.treat_undocumented_variations_as_reference =
            result.match_data.treat_undocumented_variations_as_reference;
    }

    /// Adds VCF warning messages to matching variant reports, like Java `GeneReport.addVariantWarningMessages`.
    pub fn add_variant_warning_messages(
        &mut self,
        variant_warnings: &BTreeMap<String, BTreeSet<String>>,
    ) {
        if variant_warnings.is_empty() {
            return;
        }
        for variant in &mut self.variant_reports {
            let Some(chr_position) = variant.chr_position() else {
                continue;
            };
            if let Some(warnings) = variant_warnings.get(&chr_position) {
                variant.warnings = warnings.iter().cloned().collect();
            }
        }
    }

    /// Applies Java `GeneReport.applyMatcherMessages` static message rules ported so far.
    pub fn apply_matcher_static_messages(
        &mut self,
        result: &GeneCallResult,
        catalog: &MessageCatalog,
    ) {
        if result_has_combination_or_partial_call(result) && !is_lowest_function_gene(&result.gene)
        {
            self.add_static_message(catalog, "pcat-combo-naming");
            if !result.match_data.phased {
                self.add_static_message(catalog, "pcat-combo-unphased");
            }
        }

        if result.gene == "CYP2D6" {
            self.add_static_message(catalog, "pcat-cyp2d6-research-mode");
            self.add_static_message(catalog, "pcat-cyp2d6-gene-note");
        }
    }

    fn add_static_message(&mut self, catalog: &MessageCatalog, key: &str) {
        if let Some(message) = catalog.message(key) {
            self.messages.insert(message.clone());
        }
    }

    /// Applies Java `MessageHelper.addMatchingMessagesTo` gene-level matcher messages.
    pub fn apply_matching_gene_messages(&mut self, catalog: &MessageCatalog) {
        let candidate_messages = catalog.messages_for_gene(&self.gene).to_vec();
        if !self.is_reportable_like_java() {
            if self.is_no_data_like_java() {
                return;
            }
            for message in candidate_messages {
                if message.exception_type == MessageAnnotation::TYPE_NONMATCH
                    && self.matches_gene_message(&message)
                {
                    self.messages.insert(message);
                }
            }
            return;
        }

        if self.outside_call {
            return;
        }
        for message in candidate_messages {
            if message.exception_type != MessageAnnotation::TYPE_NONMATCH
                && self.matches_gene_message(&message)
            {
                self.messages.insert(message);
            }
        }
    }

    fn matches_gene_message(&self, message: &MessageAnnotation) -> bool {
        let matches = &message.matches;
        if matches.gene.as_deref() != Some(self.gene.as_str()) {
            return false;
        }
        if !matches
            .haps_called
            .iter()
            .all(|haplotype| self.has_haplotype_like_java(haplotype))
        {
            return false;
        }
        if !matches
            .haps_missing
            .iter()
            .all(|haplotype| self.uncalled_haplotypes.contains(haplotype))
        {
            return false;
        }
        if !matches.variants.iter().all(|rsid| {
            self.find_variant_report(rsid)
                .is_some_and(|variant| !variant.is_missing())
        }) {
            return false;
        }
        if !matches.variants_missing.iter().all(|rsid| {
            self.find_variant_report(rsid)
                .is_some_and(VariantReport::is_missing)
        }) {
            return false;
        }
        if !matches
            .dips
            .iter()
            .all(|diplotype| self.has_source_diplotype_label(diplotype))
        {
            return false;
        }
        if message.exception_type == MessageAnnotation::TYPE_AMBIGUITY {
            if !matches.dips.is_empty() && self.phased {
                return false;
            }
            if !matches.variants.is_empty()
                && !matches.variants.iter().all(|rsid| {
                    self.find_variant_report(rsid)
                        .is_some_and(VariantReport::is_het_call)
                })
            {
                return false;
            }
        }
        true
    }

    fn is_reportable_like_java(&self) -> bool {
        !self.recommendation_diplotypes.is_empty()
            && self
                .recommendation_diplotypes
                .iter()
                .all(|diplotype| !diplotype.is_unknown())
    }

    fn is_no_data_like_java(&self) -> bool {
        !self.outside_call
            && (self.variant_reports.is_empty()
                || self.variant_reports.iter().all(VariantReport::is_missing))
    }

    fn has_haplotype_like_java(&self, haplotype: &str) -> bool {
        self.recommendation_diplotypes
            .iter()
            .any(|diplotype| diplotype.has_allele(haplotype))
    }

    fn find_variant_report(&self, rsid: &str) -> Option<&VariantReport> {
        self.variant_reports.iter().find(|variant| {
            variant
                .db_snp_id
                .as_deref()
                .is_some_and(|db_snp_id| db_snp_id.contains(rsid))
        })
    }

    fn has_source_diplotype_label(&self, label: &str) -> bool {
        self.source_diplotypes
            .iter()
            .any(|diplotype| diplotype.label == label)
    }

    fn apply_definition_haplotype_metadata(&mut self, definition: &DefinitionFile) {
        for diplotype in &mut self.source_diplotypes {
            diplotype.apply_definition_haplotype_metadata(definition);
        }
        for diplotype in &mut self.matcher_component_diplotypes {
            diplotype.apply_definition_haplotype_metadata(definition);
        }
        for diplotype in &mut self.recommendation_diplotypes {
            diplotype.apply_definition_haplotype_metadata(definition);
        }
        sort_report_diplotypes(&mut self.source_diplotypes);
        sort_report_diplotypes(&mut self.matcher_component_diplotypes);
        sort_report_diplotypes(&mut self.recommendation_diplotypes);
        self.source_diplotype = self
            .source_diplotypes
            .first()
            .map(|diplotype| diplotype.label.clone());
    }

    fn apply_definition_report_metadata(
        &mut self,
        definition: &DefinitionFile,
        phenotype: Option<&GenePhenotype>,
    ) {
        if definition.gene_symbol != self.gene {
            return;
        }
        self.allele_definition_version = definition
            .version
            .clone()
            .or_else(|| definition.data_version.clone());
        self.allele_definition_source = data_source_from_definition(definition);
        if self.phenotype_version.is_none() {
            self.phenotype_version = phenotype.and_then(|phenotype| phenotype.version.clone());
        }
    }

    fn add_reference_allele_message(
        &mut self,
        result: &GeneCallResult,
        definition: &DefinitionFile,
    ) {
        if result.gene == "CFTR" || definition.variants.len() <= 1 {
            return;
        }
        let Some(reference_haplotype) = first_reference_haplotype_name(result) else {
            return;
        };

        let mut message = format!(
            "The {} {} allele assignment is characterized by the absence of variants at the positions that are included in the underlying allele definitions",
            result.gene, reference_haplotype
        );
        if self.is_missing_variants_like_java() {
            message.push_str(" either because the position is reference or missing");
        }
        message.push('.');
        self.messages
            .insert(MessageAnnotation::new_note("reference-allele", message));
    }

    fn is_missing_variants_like_java(&self) -> bool {
        !self.outside_call
            && (self.variant_reports.is_empty()
                || self.variant_reports.iter().any(VariantReport::is_missing))
    }

    /// Applies Java's SLCO1B1 rs4149056 custom recommendation fallback.
    pub fn apply_slco1b1_custom_recommendation(
        &mut self,
        phenotype: Option<&GenePhenotype>,
    ) -> Result<(), Slco1b1CustomCallError> {
        if self.gene != "SLCO1B1"
            || self.outside_call
            || !self
                .source_diplotypes
                .iter()
                .all(report_diplotype_is_unknown)
        {
            return Ok(());
        }

        let variants = self
            .variant_reports
            .iter()
            .filter(|variant| variant.db_snp_id.as_deref() == Some("rs4149056"))
            .collect::<Vec<_>>();
        if variants.is_empty() {
            return Ok(());
        }
        if variants.len() > 1 {
            return Err(Slco1b1CustomCallError::MultipleRs4149056Reports);
        }

        let variant = variants[0].clone();
        let Some(alleles) = split_slco1b1_rs4149056_variant_call(&variant) else {
            return Ok(());
        };
        let Some(haplotype1) = slco1b1_rs4149056_allele_to_haplotype(&alleles[0]) else {
            return Ok(());
        };
        let Some(haplotype2) = slco1b1_rs4149056_allele_to_haplotype(&alleles[1]) else {
            return Ok(());
        };

        let annotated =
            DiplotypeAnnotationInput::from_alleles("SLCO1B1", haplotype1, Some(haplotype2))
                .annotate(phenotype)?;
        let mut recommendation = ReportDiplotype::from_annotated(&annotated, phenotype);
        recommendation.match_score = Some("0".to_owned());
        recommendation.variant = Some(variant.clone());
        recommendation.inferred = true;
        recommendation.inferred_source_diplotypes = vec![ReportDiplotype::from_match_label(
            "SLCO1B1",
            &format!("rs4149056 {}/rs4149056 {}", alleles[0], alleles[1]),
            Some(0),
            phenotype,
        )];

        if let Some(source) = self.source_diplotypes.first_mut() {
            source.variant = Some(variant);
        }
        self.recommendation_diplotypes = vec![recommendation.clone()];
        self.lookup_keys = recommendation.lookup_keys.clone();
        self.phenotypes = recommendation.phenotypes.clone();
        self.activity_score = recommendation.activity_score.clone();
        self.diplotype_key = recommendation.diplotype_key.clone();
        self.match_score = recommendation.match_score.clone();
        sort_report_diplotypes(&mut self.recommendation_diplotypes);
        Ok(())
    }
}

/// Minimal Java `Diplotype` report surface.
#[derive(Clone, Debug, Default, Eq, PartialEq, Serialize)]
pub struct ReportDiplotype {
    /// Gene symbol.
    pub gene: String,
    /// Display label.
    pub label: String,
    /// First allele, if known.
    #[serde(rename = "allele1", skip_serializing_if = "Option::is_none")]
    pub allele1: Option<ReportHaplotype>,
    /// Second allele, if known.
    #[serde(rename = "allele2", skip_serializing_if = "Option::is_none")]
    pub allele2: Option<ReportHaplotype>,
    /// Recommendation lookup keys.
    #[serde(rename = "lookupKey", skip_serializing_if = "Vec::is_empty")]
    pub lookup_keys: Vec<String>,
    /// Phenotypes assigned to this diplotype.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub phenotypes: Vec<String>,
    /// Activity score assigned to this diplotype.
    #[serde(rename = "activityScore", skip_serializing_if = "Option::is_none")]
    pub activity_score: Option<String>,
    /// Java `Diplotype.getDiplotypeKey`-style allele counts.
    #[serde(
        rename = "diplotypeKey",
        skip_serializing_if = "serde_json::Map::is_empty"
    )]
    pub diplotype_key: serde_json::Map<String, Value>,
    /// Matcher score, when this diplotype came directly from matcher output.
    #[serde(rename = "matchScore", skip_serializing_if = "Option::is_none")]
    pub match_score: Option<String>,
    /// Whether the phenotype came from an outside call.
    #[serde(rename = "outsidePhenotype")]
    pub outside_phenotype: bool,
    /// Expected phenotype when an outside phenotype mismatches.
    #[serde(
        rename = "outsidePhenotypeMismatch",
        skip_serializing_if = "Option::is_none"
    )]
    pub outside_phenotype_mismatch: Option<String>,
    /// Whether the activity score came from an outside call.
    #[serde(rename = "outsideActivityScore")]
    pub outside_activity_score: bool,
    /// Expected activity score when an outside activity score mismatches.
    #[serde(
        rename = "outsideActivityScoreMismatch",
        skip_serializing_if = "Option::is_none"
    )]
    pub outside_activity_score_mismatch: Option<String>,
    /// Variant used to make this diplotype call.
    #[serde(rename = "variant", skip_serializing_if = "Option::is_none")]
    pub variant: Option<VariantReport>,
    /// Whether this diplotype was inferred from a lower-level source call.
    pub inferred: bool,
    /// Source diplotypes used to infer this diplotype.
    #[serde(
        rename = "inferredSourceDiplotypes",
        skip_serializing_if = "Vec::is_empty"
    )]
    pub inferred_source_diplotypes: Vec<ReportDiplotype>,
    /// Whether this diplotype contains a Java combination haplotype name.
    pub combination: bool,
}

impl ReportDiplotype {
    fn from_annotated(diplotype: &AnnotatedDiplotype, phenotype: Option<&GenePhenotype>) -> Self {
        let allele1 = diplotype
            .allele1
            .as_deref()
            .map(|allele| ReportHaplotype::new(&diplotype.gene, allele, phenotype));
        let allele2 = diplotype
            .allele2
            .as_deref()
            .map(|allele| ReportHaplotype::new(&diplotype.gene, allele, phenotype));
        let label = report_diplotype_label(
            &diplotype.gene,
            allele1.as_ref(),
            allele2.as_ref(),
            &diplotype.phenotypes,
            diplotype.activity_score.as_deref(),
        );
        let lookup_keys = if diplotype.lookup_keys.is_empty() {
            diplotype.phenotypes.clone()
        } else {
            diplotype.lookup_keys.clone()
        };
        let combination = allele1
            .as_ref()
            .is_some_and(ReportHaplotype::is_combination)
            || allele2
                .as_ref()
                .is_some_and(ReportHaplotype::is_combination);
        Self {
            gene: diplotype.gene.clone(),
            label,
            allele1,
            allele2,
            lookup_keys,
            phenotypes: diplotype.phenotypes.clone(),
            activity_score: diplotype.activity_score.clone(),
            diplotype_key: diplotype_key_from_alleles(
                diplotype.allele1.as_deref(),
                diplotype.allele2.as_deref(),
            ),
            match_score: None,
            outside_phenotype: diplotype.outside_phenotype,
            outside_phenotype_mismatch: diplotype.outside_phenotype_mismatch.clone(),
            outside_activity_score: diplotype.outside_activity_score,
            outside_activity_score_mismatch: diplotype.outside_activity_score_mismatch.clone(),
            combination,
            ..Self::default()
        }
    }

    fn from_match_label(
        gene: &str,
        label: &str,
        match_score: Option<i32>,
        phenotype: Option<&GenePhenotype>,
    ) -> Self {
        let (allele1, allele2) = split_report_diplotype_label(label);
        let haplotype1 = allele1
            .as_deref()
            .map(|allele| ReportHaplotype::new(gene, allele, phenotype));
        let haplotype2 = allele2
            .as_deref()
            .map(|allele| ReportHaplotype::new(gene, allele, phenotype));
        let label = if haplotype1.is_some() || haplotype2.is_some() {
            report_diplotype_label(gene, haplotype1.as_ref(), haplotype2.as_ref(), &[], None)
        } else {
            label.to_owned()
        };
        let combination = haplotype1
            .as_ref()
            .is_some_and(ReportHaplotype::is_combination)
            || haplotype2
                .as_ref()
                .is_some_and(ReportHaplotype::is_combination)
            || label.contains(" + ");
        Self {
            gene: gene.to_owned(),
            label,
            diplotype_key: diplotype_key_from_alleles(allele1.as_deref(), allele2.as_deref()),
            allele1: haplotype1,
            allele2: haplotype2,
            match_score: match_score.map(|score| score.to_string()),
            combination,
            ..Self::default()
        }
    }

    fn apply_definition_haplotype_metadata(&mut self, definition: &DefinitionFile) {
        if let Some(allele1) = &mut self.allele1 {
            allele1.apply_definition_metadata(definition);
        }
        if let Some(allele2) = &mut self.allele2 {
            allele2.apply_definition_metadata(definition);
        }
        for diplotype in &mut self.inferred_source_diplotypes {
            diplotype.apply_definition_haplotype_metadata(definition);
        }
        self.label = report_diplotype_label(
            &self.gene,
            self.allele1.as_ref(),
            self.allele2.as_ref(),
            &self.phenotypes,
            self.activity_score.as_deref(),
        );
    }

    fn has_allele(&self, haplotype: &str) -> bool {
        self.allele1
            .as_ref()
            .is_some_and(|allele| allele.name == haplotype)
            || self
                .allele2
                .as_ref()
                .is_some_and(|allele| allele.name == haplotype)
    }

    fn is_unknown(&self) -> bool {
        self.has_allele("Unknown") || self.label.contains("Unknown")
    }
}

/// Minimal Java `Haplotype` report surface nested under report diplotypes.
#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct ReportHaplotype {
    /// Gene symbol.
    pub gene: String,
    /// Haplotype name.
    pub name: String,
    /// Haplotype function.
    pub function: String,
    /// Whether this is the reference haplotype.
    pub reference: bool,
    /// Haplotype activity value, when this gene uses activity scores.
    #[serde(rename = "activityValue", skip_serializing_if = "Option::is_none")]
    pub activity_value: Option<String>,
}

impl ReportHaplotype {
    fn new(gene: &str, name: &str, phenotype: Option<&GenePhenotype>) -> Self {
        Self {
            gene: gene.to_owned(),
            name: name.to_owned(),
            function: phenotype
                .map(|phenotype| phenotype.haplotype_function(name).to_owned())
                .unwrap_or_else(|| GenePhenotype::UNASSIGNED_FUNCTION.to_owned()),
            reference: name == "Reference",
            activity_value: phenotype
                .and_then(|phenotype| phenotype.haplotype_activity(name))
                .map(str::to_owned),
        }
    }

    fn apply_definition_metadata(&mut self, definition: &DefinitionFile) {
        if definition.gene_symbol != self.gene {
            return;
        }
        if let Some(named_allele) = definition.named_allele(&self.name) {
            self.reference = named_allele.reference;
        }
    }

    /// Java `Haplotype.toString`.
    pub fn display_name(&self) -> String {
        if self.name.starts_with('*') {
            format!("{}{}", self.gene, self.name)
        } else {
            format!("{} {}", self.gene, self.name)
        }
    }

    fn is_combination(&self) -> bool {
        self.name.starts_with('[') && self.name.contains(" + ") && self.name.ends_with(']')
    }
}

/// Minimal Java `VariantReport` surface used by report message rules.
#[derive(Clone, Debug, Default, Eq, PartialEq, Serialize)]
pub struct VariantReport {
    /// Gene symbol.
    #[serde(rename = "gene", skip_serializing_if = "Option::is_none")]
    pub gene: Option<String>,
    /// Chromosome or contig name.
    #[serde(rename = "chromosome", skip_serializing_if = "Option::is_none")]
    pub chromosome: Option<String>,
    /// dbSNP identifier.
    #[serde(rename = "dbSnpId")]
    pub db_snp_id: Option<String>,
    /// CPIC/VCF position.
    #[serde(rename = "position", skip_serializing_if = "Option::is_none")]
    pub position: Option<i64>,
    /// VCF call string.
    #[serde(rename = "call", skip_serializing_if = "Option::is_none")]
    pub call: Option<String>,
    /// Reference allele.
    #[serde(rename = "referenceAllele", skip_serializing_if = "Option::is_none")]
    pub reference_allele: Option<String>,
    /// Haplotype alleles associated with this variant.
    #[serde(rename = "alleles", skip_serializing_if = "Vec::is_empty")]
    pub alleles: Vec<String>,
    /// Whether the variant call was phased.
    #[serde(rename = "phased")]
    pub phased: bool,
    /// Java-retained phase set.
    #[serde(rename = "phaseSet", skip_serializing_if = "Option::is_none")]
    pub phase_set: Option<i32>,
    /// Whether the call has undocumented variation.
    #[serde(rename = "hasUndocumentedVariations")]
    pub has_undocumented_variations: bool,
    /// VCF warnings attached to this variant report.
    #[serde(rename = "warnings", skip_serializing_if = "Vec::is_empty")]
    pub warnings: Vec<String>,
}

impl VariantReport {
    /// Creates a variant report with a dbSNP id and optional call.
    pub fn new(db_snp_id: impl Into<String>, call: Option<impl Into<String>>) -> Self {
        Self {
            gene: None,
            chromosome: None,
            db_snp_id: Some(db_snp_id.into()),
            position: None,
            call: call.map(Into::into),
            reference_allele: None,
            alleles: Vec::new(),
            phased: false,
            phase_set: None,
            has_undocumented_variations: false,
            warnings: Vec::new(),
        }
    }

    /// Creates a Java `VariantReportFactory`-style report from matcher data.
    pub fn from_matcher_locus(
        gene: &str,
        locus: &VariantLocus,
        sample_allele: Option<&crate::matcher::SampleAllele>,
        definition: &DefinitionFile,
    ) -> Self {
        Self {
            gene: Some(gene.to_owned()),
            chromosome: Some(locus.chromosome.clone()),
            db_snp_id: locus.rsid.clone(),
            position: Some(locus.position as i64),
            call: sample_allele
                .map(crate::matcher::SampleAllele::vcf_call)
                .filter(|call| is_valid_variant_call(call))
                .map(str::to_owned),
            reference_allele: Some(locus.reference.clone()),
            alleles: variant_report_alleles(definition, locus),
            phased: sample_allele.is_some_and(crate::matcher::SampleAllele::phased),
            phase_set: sample_allele.and_then(crate::matcher::SampleAllele::phase_set),
            has_undocumented_variations: false,
            warnings: Vec::new(),
        }
    }

    /// Sets CPIC/VCF position.
    pub fn with_position(mut self, position: i64) -> Self {
        self.position = Some(position);
        self
    }

    /// Sets reference allele.
    pub fn with_reference_allele(mut self, reference_allele: impl Into<String>) -> Self {
        self.reference_allele = Some(reference_allele.into());
        self
    }

    /// Sets Java `VariantReport.getAlleles`.
    pub fn with_alleles(mut self, alleles: impl IntoIterator<Item = String>) -> Self {
        self.alleles = alleles.into_iter().collect();
        self.alleles
            .sort_by(|left, right| compare_haplotype_names(left, right));
        self.alleles.dedup();
        self
    }

    /// Sets phasing metadata.
    pub fn with_phasing(mut self, phased: bool, phase_set: Option<i32>) -> Self {
        self.phased = phased;
        self.phase_set = phase_set;
        self
    }

    /// Java `VariantReport.toChrPosition`.
    pub fn chr_position(&self) -> Option<String> {
        Some(format!(
            "{}:{}",
            self.chromosome.as_deref()?,
            self.position?
        ))
    }

    /// Sets undocumented-variation flag.
    pub fn with_undocumented_variations(mut self, has_undocumented_variations: bool) -> Self {
        self.has_undocumented_variations = has_undocumented_variations;
        self
    }

    /// Java `VariantReport.isMissing`.
    pub fn is_missing(&self) -> bool {
        self.call
            .as_deref()
            .is_none_or(|call| call.trim().is_empty())
    }

    /// Java `VariantReport.isHetCall`.
    pub fn is_het_call(&self) -> bool {
        let Some(call) = self.call.as_deref() else {
            return false;
        };
        let mut alleles = call.split(['|', '/']);
        let Some(allele1) = alleles.next() else {
            return false;
        };
        let Some(allele2) = alleles.next() else {
            return false;
        };
        alleles.next().is_none()
            && !allele1.trim().is_empty()
            && !allele2.trim().is_empty()
            && allele1
                .chars()
                .all(|character| character.is_ascii_alphabetic())
            && allele2
                .chars()
                .all(|character| character.is_ascii_alphabetic())
            && allele1 != allele2
    }

    /// Java `VariantReport.isNonReference`.
    pub fn is_non_reference(&self) -> bool {
        let Some(call) = self.call.as_deref() else {
            return false;
        };
        let Some(reference_allele) = self.reference_allele.as_deref() else {
            return false;
        };
        !call.trim().is_empty()
            && call
                .split(['|', '/'])
                .any(|allele| allele != reference_allele)
    }

    fn highlighted_call(&self) -> Option<String> {
        (!self.is_missing()).then(|| self.call.as_deref().unwrap_or("").replace('|', "/"))
    }
}

impl Ord for VariantReport {
    fn cmp(&self, other: &Self) -> Ordering {
        self.is_missing()
            .cmp(&other.is_missing())
            .then_with(|| self.chromosome.cmp(&other.chromosome))
            .then_with(|| self.position.cmp(&other.position))
    }
}

impl PartialOrd for VariantReport {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

fn diplotype_key_from_alleles(
    allele1: Option<&str>,
    allele2: Option<&str>,
) -> serde_json::Map<String, Value> {
    let mut counts = BTreeMap::<String, i32>::new();
    if let Some(allele1) = allele1 {
        *counts.entry(allele1.to_owned()).or_default() += 1;
    }
    if let Some(allele2) = allele2 {
        *counts.entry(allele2.to_owned()).or_default() += 1;
    }

    counts
        .into_iter()
        .map(|(allele, count)| {
            (
                allele,
                Value::Number(
                    serde_json::Number::from_f64(f64::from(count))
                        .expect("integer diplotype count converts to finite JSON number"),
                ),
            )
        })
        .collect()
}

fn sort_report_diplotypes(diplotypes: &mut [ReportDiplotype]) {
    diplotypes.sort_by(compare_report_diplotypes_like_java);
    for diplotype in diplotypes {
        sort_report_diplotypes(&mut diplotype.inferred_source_diplotypes);
    }
}

fn compare_report_diplotypes_like_java(
    left: &ReportDiplotype,
    right: &ReportDiplotype,
) -> Ordering {
    left.gene
        .cmp(&right.gene)
        .then_with(|| {
            if left.label != right.label {
                compare_report_haplotypes_like_java(left.allele1.as_ref(), right.allele1.as_ref())
                    .then_with(|| {
                        compare_report_haplotypes_like_java(
                            left.allele2.as_ref(),
                            right.allele2.as_ref(),
                        )
                    })
            } else {
                Ordering::Equal
            }
        })
        .then_with(|| left.inferred.cmp(&right.inferred))
}

fn compare_report_haplotypes_like_java(
    left: Option<&ReportHaplotype>,
    right: Option<&ReportHaplotype>,
) -> Ordering {
    match (left, right) {
        (None, None) => Ordering::Equal,
        (None, Some(_)) => Ordering::Less,
        (Some(_), None) => Ordering::Greater,
        (Some(left), Some(right)) => compare_haplotype_names(&left.name, &right.name),
    }
}

fn report_diplotype_label(
    gene: &str,
    allele1: Option<&ReportHaplotype>,
    allele2: Option<&ReportHaplotype>,
    phenotypes: &[String],
    activity_score: Option<&str>,
) -> String {
    if cpic_style_diplotype_gene(gene)
        && let Some(allele1) = allele1
    {
        let allele1_reference = allele1.reference;
        let allele2_reference = allele2.is_some_and(|allele| allele.reference);
        if allele1_reference && allele2_reference {
            return "Reference/Reference".to_owned();
        }
        if allele1_reference || allele2_reference {
            if let Some(allele2) = allele2 {
                let allele = if allele1_reference {
                    &allele2.name
                } else {
                    &allele1.name
                };
                return format!("{allele} (heterozygous)");
            }
            return allele1.name.clone();
        }
    }

    match (allele1, allele2) {
        (Some(allele1), Some(allele2)) => {
            let mut alleles = [allele1.name.as_str(), allele2.name.as_str()];
            alleles.sort_by(|left, right| compare_haplotype_names(left, right));
            alleles.join("/")
        }
        (Some(allele1), None) => allele1.name.clone(),
        (None, Some(_)) => crate::phenotype::NA.to_owned(),
        (None, None) => {
            if let Some(activity_score) = activity_score.filter(|score| !is_unspecified(score)) {
                if phenotypes.is_empty() || phenotypes.iter().any(|value| value == "No Result") {
                    activity_score.to_owned()
                } else {
                    format!(
                        "{phenotype} ({activity_score})",
                        phenotype = phenotypes.join("/")
                    )
                }
            } else if !phenotypes.is_empty() && phenotypes.iter().all(|value| value != "No Result")
            {
                phenotypes.join("/")
            } else {
                crate::phenotype::NA.to_owned()
            }
        }
    }
}

fn cpic_style_diplotype_gene(gene: &str) -> bool {
    matches!(gene, "CACNA1S" | "CFTR" | "DPYD" | "RYR1")
}

fn is_unspecified(value: &str) -> bool {
    value.trim().is_empty() || value.trim().eq_ignore_ascii_case(crate::phenotype::NA)
}

fn allele_function_map_from_phenotype(
    phenotype: Option<&GenePhenotype>,
) -> BTreeMap<String, String> {
    phenotype
        .map(GenePhenotype::formatted_function_score_map)
        .unwrap_or_default()
}

fn is_single_ploidy_report_gene(gene: &str) -> bool {
    matches!(gene, "G6PD" | "MT-RNR1")
}

fn source_diplotype_label(diplotype: &AnnotatedDiplotype) -> Option<String> {
    match (&diplotype.allele1, &diplotype.allele2) {
        (Some(allele1), Some(allele2)) => Some(format!("{allele1}/{allele2}")),
        (Some(allele1), None) => Some(allele1.clone()),
        (None, Some(allele2)) => Some(allele2.clone()),
        (None, None) => {
            // Match report_diplotype_label: an activity-score-only call prints "phenotype (score)".
            if let Some(score) = diplotype
                .activity_score
                .as_deref()
                .filter(|score| !is_unspecified(score))
            {
                if diplotype.phenotypes.is_empty()
                    || diplotype
                        .phenotypes
                        .iter()
                        .any(|value| value == "No Result")
                {
                    Some(score.to_owned())
                } else {
                    Some(format!("{} ({score})", diplotype.phenotypes.join("/")))
                }
            } else {
                diplotype
                    .phenotypes
                    .first()
                    .cloned()
                    .or_else(|| diplotype.activity_score.clone())
            }
        }
    }
}

fn lowest_function_annotation_input(
    result: &GeneCallResult,
    phenotype: &GenePhenotype,
) -> Option<crate::phenotype::DiplotypeAnnotationInput> {
    if result.gene != phenotype.gene {
        return None;
    }

    match (result.gene.as_str(), &result.kind) {
        (_, GeneCallKind::NoCall) => None,
        ("DPYD", GeneCallKind::Diplotypes(diplotypes)) => phenotype
            .infer_dpyd_lowest_function_from_diplotypes(
                diplotypes.iter().map(|diplotype| diplotype.name.as_str()),
            ),
        ("DPYD", GeneCallKind::Haplotypes(haplotypes)) => phenotype
            .infer_dpyd_lowest_function_from_haplotypes(
                haplotypes.iter().map(|haplotype| haplotype.name.as_str()),
            ),
        ("RYR1", GeneCallKind::Diplotypes(diplotypes)) => phenotype
            .infer_ryr1_lowest_function_from_diplotypes(
                diplotypes.iter().map(|diplotype| diplotype.name.as_str()),
            ),
        ("RYR1", GeneCallKind::Haplotypes(haplotypes)) => phenotype
            .infer_ryr1_lowest_function_from_haplotypes(
                haplotypes.iter().map(|haplotype| haplotype.name.as_str()),
            ),
        _ => None,
    }
}

fn source_report_diplotypes(
    result: &GeneCallResult,
    phenotype: &GenePhenotype,
) -> Vec<ReportDiplotype> {
    match &result.kind {
        GeneCallKind::NoCall => Vec::new(),
        GeneCallKind::Diplotypes(diplotypes) => diplotypes
            .iter()
            .map(|diplotype| {
                ReportDiplotype::from_match_label(
                    &result.gene,
                    &diplotype.name,
                    Some(diplotype.score),
                    Some(phenotype),
                )
            })
            .collect(),
        GeneCallKind::Haplotypes(haplotypes) => haplotypes
            .iter()
            .map(|haplotype| {
                ReportDiplotype::from_match_label(
                    &result.gene,
                    &haplotype.name,
                    None,
                    Some(phenotype),
                )
            })
            .collect(),
    }
}

fn matcher_homozygous_component_haplotypes(result: &GeneCallResult) -> BTreeSet<String> {
    let GeneCallKind::Haplotypes(haplotypes) = &result.kind else {
        return BTreeSet::new();
    };
    let mut counts = BTreeMap::<String, usize>::new();
    for haplotype in haplotypes {
        for name in report_component_names_for_match(haplotype) {
            *counts.entry(name).or_default() += 1;
        }
    }
    counts
        .into_iter()
        .filter_map(|(name, count)| (count > 1).then_some(name))
        .collect()
}

fn matcher_component_haplotypes(result: &GeneCallResult) -> BTreeSet<String> {
    match &result.kind {
        GeneCallKind::NoCall => BTreeSet::new(),
        GeneCallKind::Diplotypes(diplotypes) => diplotypes
            .iter()
            .flat_map(|diplotype| {
                report_component_names_for_match(&diplotype.haplotype1)
                    .into_iter()
                    .chain(
                        diplotype
                            .haplotype2
                            .as_ref()
                            .into_iter()
                            .flat_map(report_component_names_for_match),
                    )
            })
            .collect(),
        GeneCallKind::Haplotypes(haplotypes) => haplotypes
            .iter()
            .flat_map(report_component_names_for_match)
            .collect(),
    }
}

fn matcher_component_report_diplotypes(
    result: &GeneCallResult,
    phenotype: &GenePhenotype,
) -> Vec<ReportDiplotype> {
    matcher_component_haplotypes(result)
        .into_iter()
        .map(|haplotype| {
            ReportDiplotype::from_match_label(&result.gene, &haplotype, None, Some(phenotype))
        })
        .collect()
}

fn report_component_names_for_match(haplotype: &crate::matcher::HaplotypeMatch) -> Vec<String> {
    if haplotype.haplotype.is_combination_or_partial {
        let components = haplotype
            .haplotype
            .id
            .split(" + ")
            .map(str::to_owned)
            .filter(|name| !report_is_partial_name(name))
            .collect::<Vec<_>>();
        if !components.is_empty() {
            return components;
        }
    }

    if is_combination_label(&haplotype.name) {
        report_split_combination_name(&haplotype.name)
            .into_iter()
            .filter(|name| !report_is_partial_name(name))
            .collect()
    } else {
        vec![haplotype.name.clone()]
    }
}

fn report_is_partial_name(name: &str) -> bool {
    name.starts_with('(') && name.ends_with(')')
}

fn report_split_combination_name(name: &str) -> Vec<String> {
    name.trim_start_matches('[')
        .trim_end_matches(']')
        .split(" + ")
        .map(str::to_owned)
        .collect()
}

fn inferred_lowest_function_diplotype(
    result: &GeneCallResult,
    recommendation: &ReportDiplotype,
) -> bool {
    if recommendation.combination {
        return true;
    }

    match &result.kind {
        GeneCallKind::NoCall => false,
        GeneCallKind::Diplotypes(diplotypes) => diplotypes
            .iter()
            .any(|diplotype| diplotype.name.split('/').any(is_combination_label)),
        GeneCallKind::Haplotypes(haplotypes) => {
            haplotypes.len() > 2
                || haplotypes
                    .iter()
                    .any(|haplotype| is_combination_label(&haplotype.name))
        }
    }
}

fn result_has_combination_or_partial_call(result: &GeneCallResult) -> bool {
    match &result.kind {
        GeneCallKind::NoCall => false,
        GeneCallKind::Diplotypes(diplotypes) => diplotypes.iter().any(|diplotype| {
            haplotype_match_is_combination_or_partial(&diplotype.haplotype1)
                || diplotype
                    .haplotype2
                    .as_ref()
                    .is_some_and(haplotype_match_is_combination_or_partial)
                || diplotype.name.split('/').any(is_combination_label)
        }),
        GeneCallKind::Haplotypes(haplotypes) => haplotypes
            .iter()
            .any(haplotype_match_is_combination_or_partial),
    }
}

fn uncalled_haplotypes_for_match_data(
    match_data: &MatchData,
    definition: &DefinitionFile,
) -> BTreeSet<String> {
    let matchable_haplotypes = match_data
        .haplotypes()
        .iter()
        .map(|haplotype| haplotype.name.as_str())
        .collect::<BTreeSet<_>>();
    definition
        .named_alleles
        .iter()
        .map(|haplotype| haplotype.name.clone())
        .filter(|name| !matchable_haplotypes.contains(name.as_str()))
        .collect()
}

fn first_reference_haplotype_name(result: &GeneCallResult) -> Option<&str> {
    match &result.kind {
        GeneCallKind::NoCall => None,
        GeneCallKind::Diplotypes(diplotypes) => diplotypes.iter().find_map(|diplotype| {
            reference_haplotype_name(&diplotype.haplotype1).or_else(|| {
                diplotype
                    .haplotype2
                    .as_ref()
                    .and_then(reference_haplotype_name)
            })
        }),
        GeneCallKind::Haplotypes(haplotypes) => {
            haplotypes.iter().find_map(reference_haplotype_name)
        }
    }
}

fn reference_haplotype_name(haplotype: &crate::matcher::HaplotypeMatch) -> Option<&str> {
    haplotype
        .haplotype
        .reference
        .then_some(haplotype.name.as_str())
}

fn haplotype_match_is_combination_or_partial(haplotype: &crate::matcher::HaplotypeMatch) -> bool {
    haplotype.haplotype.is_combination_or_partial || is_combination_label(&haplotype.name)
}

fn is_lowest_function_gene(gene: &str) -> bool {
    matches!(gene, "DPYD" | "RYR1")
}

fn split_report_diplotype_label(label: &str) -> (Option<String>, Option<String>) {
    let Some((allele1, allele2)) = label.split_once('/') else {
        return (clean_report_allele_label(label), None);
    };
    (
        clean_report_allele_label(allele1),
        clean_report_allele_label(allele2),
    )
}

fn clean_report_allele_label(label: &str) -> Option<String> {
    let trimmed = label.trim();
    if trimmed.is_empty() {
        None
    } else {
        Some(trimmed.to_owned())
    }
}

fn is_combination_label(label: &str) -> bool {
    let trimmed = label.trim();
    trimmed.starts_with('[') && trimmed.contains(" + ") && trimmed.ends_with(']')
}

fn merge_report_gene(existing: &mut ReportGene, incoming: ReportGene) {
    if report_gene_summary_precedes_like_java(&incoming, existing) {
        merge_report_gene_prefer_incoming_summary(existing, incoming);
        return;
    }

    let incoming_call_source = incoming.call_source;
    let incoming_guidance_source = incoming.guidance_source;
    append_unique(&mut existing.lookup_keys, incoming.lookup_keys);
    append_unique(&mut existing.phenotypes, incoming.phenotypes);
    for (allele, count) in incoming.diplotype_key {
        existing.diplotype_key.entry(allele).or_insert(count);
    }

    if existing.activity_score.is_none() {
        existing.activity_score = incoming.activity_score;
    }
    existing.is_activity_score_type |= incoming.is_activity_score_type;
    existing.is_allele_presence_type |= incoming.is_allele_presence_type;
    if existing.chromosome.is_none() {
        existing.chromosome = incoming.chromosome;
    }
    existing.phased |= incoming.phased;
    existing.effectively_phased |= incoming.effectively_phased;
    existing.outside_call |= incoming.outside_call;
    merge_report_gene_source_metadata(existing, incoming_call_source, incoming_guidance_source);
    append_unique_items(&mut existing.source_diplotypes, incoming.source_diplotypes);
    existing
        .matcher_component_haplotypes
        .extend(incoming.matcher_component_haplotypes);
    append_unique_items(
        &mut existing.matcher_component_diplotypes,
        incoming.matcher_component_diplotypes,
    );
    existing
        .matcher_homozygous_component_haplotypes
        .extend(incoming.matcher_homozygous_component_haplotypes);
    append_unique_items(
        &mut existing.recommendation_diplotypes,
        incoming.recommendation_diplotypes,
    );
    existing
        .allele_function_map
        .extend(incoming.allele_function_map);
    existing.related_drugs.extend(incoming.related_drugs);
    sort_report_diplotypes(&mut existing.source_diplotypes);
    sort_report_diplotypes(&mut existing.matcher_component_diplotypes);
    sort_report_diplotypes(&mut existing.recommendation_diplotypes);
    existing.variant_reports.extend(incoming.variant_reports);
    existing
        .variant_of_interest_reports
        .extend(incoming.variant_of_interest_reports);
    existing.has_undocumented_variations |= incoming.has_undocumented_variations;
    existing.treat_undocumented_variations_as_reference |=
        incoming.treat_undocumented_variations_as_reference;
    existing.messages.extend(incoming.messages);
    if existing.outside_phenotype_mismatch.is_none() {
        existing.outside_phenotype_mismatch = incoming.outside_phenotype_mismatch;
    }
    if existing.outside_activity_score_mismatch.is_none() {
        existing.outside_activity_score_mismatch = incoming.outside_activity_score_mismatch;
    }

    match (&mut existing.source_diplotype, incoming.source_diplotype) {
        (Some(existing_label), Some(incoming_label)) if !existing_label.is_empty() => {
            if !existing_label
                .split("; ")
                .any(|label| label == incoming_label)
            {
                existing_label.push_str("; ");
                existing_label.push_str(&incoming_label);
            }
        }
        (None, Some(incoming_label)) => existing.source_diplotype = Some(incoming_label),
        _ => {}
    }
    if existing.match_score.is_none() {
        existing.match_score = incoming.match_score;
    }
}

fn report_gene_summary_precedes_like_java(incoming: &ReportGene, existing: &ReportGene) -> bool {
    incoming.gene == existing.gene && report_gene_source_cmp_like_java(incoming, existing).is_lt()
}

fn merge_report_gene_prefer_incoming_summary(existing: &mut ReportGene, mut incoming: ReportGene) {
    incoming
        .related_drugs
        .extend(std::mem::take(&mut existing.related_drugs));
    incoming
        .messages
        .extend(std::mem::take(&mut existing.messages));
    incoming
        .variant_reports
        .extend(std::mem::take(&mut existing.variant_reports));
    incoming
        .variant_of_interest_reports
        .extend(std::mem::take(&mut existing.variant_of_interest_reports));
    incoming
        .allele_function_map
        .extend(std::mem::take(&mut existing.allele_function_map));
    append_unique_items(
        &mut incoming.matcher_component_diplotypes,
        std::mem::take(&mut existing.matcher_component_diplotypes),
    );
    incoming
        .uncalled_haplotypes
        .extend(std::mem::take(&mut existing.uncalled_haplotypes));
    incoming.has_undocumented_variations |= existing.has_undocumented_variations;
    incoming.treat_undocumented_variations_as_reference |=
        existing.treat_undocumented_variations_as_reference;
    sort_report_diplotypes(&mut incoming.source_diplotypes);
    sort_report_diplotypes(&mut incoming.matcher_component_diplotypes);
    sort_report_diplotypes(&mut incoming.recommendation_diplotypes);
    *existing = incoming;
}

fn insert_report_gene_source(
    report_gene_sources: &mut BTreeMap<String, Vec<ReportGene>>,
    report_gene: ReportGene,
) {
    let entries = report_gene_sources
        .entry(report_gene.gene.clone())
        .or_default();
    if let Some(existing) = entries.iter_mut().find(|existing| {
        existing.call_source == report_gene.call_source
            && existing.guidance_source == report_gene.guidance_source
    }) {
        merge_report_gene(existing, report_gene);
    } else {
        entries.push(report_gene);
    }
}

fn sort_report_gene_sources_like_java(report_gene_sources: &mut BTreeMap<String, Vec<ReportGene>>) {
    for report_genes in report_gene_sources.values_mut() {
        report_genes.sort_by(report_gene_source_cmp_like_java);
    }
}

fn report_gene_source_cmp_like_java(a: &ReportGene, b: &ReportGene) -> std::cmp::Ordering {
    a.gene
        .to_lowercase()
        .cmp(&b.gene.to_lowercase())
        .then_with(|| {
            report_gene_call_source_rank_like_java(a)
                .cmp(&report_gene_call_source_rank_like_java(b))
        })
        .then_with(|| {
            report_gene_guidance_source_rank_like_java(a)
                .cmp(&report_gene_guidance_source_rank_like_java(b))
        })
}

fn report_gene_call_source_rank_like_java(report_gene: &ReportGene) -> u8 {
    match report_gene.call_source {
        ReportCallSource::Outside => 0,
        ReportCallSource::Matcher => 1,
        ReportCallSource::None => 2,
    }
}

fn report_gene_guidance_source_rank_like_java(report_gene: &ReportGene) -> u8 {
    match report_gene.guidance_source {
        Some(PrescribingGuidanceSource::CpicGuideline) => 0,
        Some(PrescribingGuidanceSource::DpwgGuideline) => 1,
        Some(PrescribingGuidanceSource::FdaLabel) => 2,
        Some(PrescribingGuidanceSource::FdaAssoc) => 3,
        None => 4,
    }
}

fn merge_report_gene_source_metadata(
    existing: &mut ReportGene,
    incoming_call_source: ReportCallSource,
    incoming_guidance_source: Option<PrescribingGuidanceSource>,
) {
    let incoming = ReportGene {
        gene: existing.gene.clone(),
        call_source: incoming_call_source,
        guidance_source: incoming_guidance_source,
        ..ReportGene::default()
    };
    if report_gene_source_cmp_like_java(&incoming, existing).is_lt() {
        existing.call_source = incoming_call_source;
        existing.guidance_source = incoming_guidance_source;
    }
}

fn append_unique(values: &mut Vec<String>, incoming: Vec<String>) {
    for value in incoming {
        if !values.contains(&value) {
            values.push(value);
        }
    }
}

fn append_unique_items<T: PartialEq>(values: &mut Vec<T>, incoming: Vec<T>) {
    for value in incoming {
        if !values.contains(&value) {
            values.push(value);
        }
    }
}

/// First Rust `ReportContext` slice.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct ReportContext {
    /// Optional title.
    #[serde(rename = "title")]
    pub title: Option<String>,
    /// Data version.
    #[serde(rename = "dataVersion")]
    pub data_version: String,
    /// Gene reports keyed by gene.
    #[serde(rename = "genes")]
    gene_reports: BTreeMap<String, ReportGene>,
    /// Java-style per-source gene reports keyed by gene, preserving same-gene reports before merge.
    #[serde(skip)]
    report_gene_sources: BTreeMap<String, Vec<ReportGene>>,
    /// Drug reports keyed by prescribing guidance source and drug name.
    #[serde(rename = "drugs")]
    drug_reports: BTreeMap<PrescribingGuidanceSource, BTreeMap<String, DrugReport>>,
    /// Global report messages.
    #[serde(rename = "messages")]
    messages: Vec<MessageAnnotation>,
    /// Unannotated gene calls.
    #[serde(rename = "unannotatedGeneCalls")]
    unannotated_gene_calls: Vec<ReportGene>,
}

impl ReportContext {
    /// Builds a report context directly from Java matcher results, definitions, and phenotypes.
    pub fn from_gene_call_results_with_definitions<'a>(
        guidance: &PgkbGuidelineCollection,
        gene_call_results: impl IntoIterator<Item = &'a GeneCallResult>,
        definitions: &DefinitionReader,
        phenotypes: &PhenotypeMap,
        title: Option<String>,
    ) -> Result<Self, ReportContextFromMatcherError> {
        let mut report_genes = Vec::new();
        for result in gene_call_results {
            let definition = definitions.definition_file(&result.gene).ok_or_else(|| {
                ReportContextFromMatcherError::MissingDefinition(result.gene.clone())
            })?;
            if let Some(report_gene) = ReportGene::from_gene_call_result_with_definition(
                result,
                phenotypes.phenotype(&result.gene),
                definition,
            )? {
                report_genes.push(report_gene);
            }
        }
        Ok(Self::from_gene_reports(guidance, report_genes, title))
    }

    /// Builds drug reports from prescribing guidance and report genes.
    pub fn from_gene_reports(
        guidance: &PgkbGuidelineCollection,
        gene_reports: impl IntoIterator<Item = ReportGene>,
        title: Option<String>,
    ) -> Self {
        let mut gene_report_map = BTreeMap::new();
        let mut report_gene_sources: BTreeMap<String, Vec<ReportGene>> = BTreeMap::new();
        for gene_report in gene_reports {
            report_gene_sources
                .entry(gene_report.gene.clone())
                .or_default()
                .push(gene_report.clone());
            match gene_report_map.get_mut(&gene_report.gene) {
                Some(existing) => merge_report_gene(existing, gene_report),
                None => {
                    gene_report_map.insert(gene_report.gene.clone(), gene_report);
                }
            }
        }
        sort_report_gene_sources_like_java(&mut report_gene_sources);
        let data_version = guidance
            .version()
            .unwrap_or(crate::phenotype::NA)
            .to_owned();
        let mut context = Self {
            title,
            data_version,
            gene_reports: gene_report_map,
            report_gene_sources,
            drug_reports: BTreeMap::new(),
            messages: Vec::new(),
            unannotated_gene_calls: Vec::new(),
        };

        for source in PrescribingGuidanceSource::list_values() {
            let mut source_reports = BTreeMap::new();
            for drug_name in guidance.chemical_names() {
                let guideline_packages = guidance.find_guideline_packages(&drug_name, source);
                if !guideline_packages.is_empty() {
                    let drug_report =
                        DrugReport::new(drug_name.clone(), guideline_packages, &context);
                    source_reports.insert(drug_name.to_lowercase(), drug_report);
                }
            }
            context.drug_reports.insert(source, source_reports);
        }
        context.rebuild_report_gene_sources_from_guidelines();
        context.apply_related_drug_backlinks();

        context
    }

    /// Adds Java report-context summary messages.
    pub fn with_messages(mut self, messages: impl IntoIterator<Item = MessageAnnotation>) -> Self {
        self.messages = messages.into_iter().collect();
        self
    }

    /// Java `getMessages`.
    pub fn messages(&self) -> &[MessageAnnotation] {
        &self.messages
    }

    /// Returns all drug reports by prescribing-guidance source.
    pub fn drug_reports(
        &self,
    ) -> &BTreeMap<PrescribingGuidanceSource, BTreeMap<String, DrugReport>> {
        &self.drug_reports
    }

    /// Java `getDrugReport`.
    pub fn drug_report(
        &self,
        source: PrescribingGuidanceSource,
        drug: &str,
    ) -> Option<&DrugReport> {
        self.drug_reports.get(&source)?.get(drug)
    }

    /// Java `getGeneReport`.
    pub fn gene_report(&self, gene: &str) -> Option<&ReportGene> {
        self.gene_reports.get(gene)
    }

    fn report_genes_for_guideline_source(
        &self,
        gene: &str,
        source: PrescribingGuidanceSource,
    ) -> Vec<ReportGene> {
        if let Some(report_genes) = self.report_gene_sources.get(gene) {
            return report_genes
                .iter()
                .cloned()
                .map(|mut report_gene| {
                    report_gene.guidance_source = Some(source);
                    report_gene
                })
                .collect();
        }
        self.gene_report(gene)
            .cloned()
            .map(|mut report_gene| {
                report_gene.guidance_source = Some(source);
                vec![report_gene]
            })
            .unwrap_or_default()
    }

    fn rebuild_report_gene_sources_from_guidelines(&mut self) {
        let mut source_reports = BTreeMap::<String, Vec<ReportGene>>::new();
        for report_gene in self
            .drug_reports
            .values()
            .flat_map(|source_reports| source_reports.values())
            .flat_map(|drug_report| drug_report.guidelines.iter())
            .flat_map(|guideline| guideline.report_genes.iter().cloned())
        {
            insert_report_gene_source(&mut source_reports, report_gene);
        }
        for (gene, report_genes) in source_reports {
            self.report_gene_sources.insert(gene, report_genes);
        }
        sort_report_gene_sources_like_java(&mut self.report_gene_sources);
    }

    fn apply_related_drug_backlinks(&mut self) {
        let mut links = Vec::new();
        for source_reports in self.drug_reports.values() {
            for drug_report in source_reports.values() {
                let drug = DrugLink::new(drug_report.name.clone(), drug_report.id.clone());
                for gene in drug_report.related_gene_symbols() {
                    links.push((gene, drug.clone()));
                }
            }
        }
        for (gene, drug) in links {
            if let Some(report_gene) = self.gene_reports.get_mut(&gene) {
                report_gene.add_related_drug(drug.clone());
            }
        }
    }

    /// Applies Java drug-level `MessageHelper.addMatchingMessagesTo` behavior.
    pub fn apply_matching_drug_messages(&mut self, catalog: &MessageCatalog) {
        let gene_reports = &self.gene_reports;
        for (source, source_reports) in &mut self.drug_reports {
            for drug_report in source_reports.values_mut() {
                drug_report.apply_matching_messages(catalog, gene_reports, *source);
                drug_report.add_missing_variant_messages(gene_reports);
            }
        }
    }

    /// Applies Java `report-as-genotype` drug messages to matching annotation reports.
    pub fn apply_report_as_genotype_messages(&mut self, catalog: &MessageCatalog) {
        self.apply_matching_drug_messages(catalog);
    }

    /// Serializes this report context as Java-style pretty JSON.
    pub fn to_json_string(&self) -> Result<String, serde_json::Error> {
        serde_json::to_string_pretty(self)
    }

    /// Writes this report context to a `.json` file.
    pub fn write_json(&self, path: &Path) -> Result<(), GuidanceLoadError> {
        write_report_json(self, path)
    }

    /// Writes the first Rust calls-only TSV surface to a `.tsv` file.
    pub fn write_calls_only_tsv(&self, path: &Path) -> Result<(), GuidanceLoadError> {
        write_calls_only_tsv(self, path, &CallsOnlyTsvOptions::default())
    }

    /// Writes the first Rust HTML report surface to a `.html` file.
    pub fn write_html(
        &self,
        path: &Path,
        options: &HtmlReportOptions,
    ) -> Result<(), GuidanceLoadError> {
        write_report_html(self, path, options)
    }
}

/// First Rust `DrugReport` slice.
#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct DrugReport {
    /// Drug name.
    #[serde(rename = "name")]
    pub name: String,
    /// Drug PharmGKB id.
    #[serde(rename = "id")]
    pub id: String,
    /// Prescribing-guidance source.
    #[serde(rename = "source")]
    pub source: PrescribingGuidanceSource,
    /// Data version.
    #[serde(rename = "version")]
    pub version: String,
    /// Drug-level messages.
    #[serde(rename = "messages")]
    pub messages: BTreeSet<MessageAnnotation>,
    /// Variants to display as part of genotype recommendations.
    #[serde(rename = "variants")]
    pub report_variants: BTreeSet<String>,
    /// URLs for guideline packages.
    #[serde(rename = "urls")]
    pub urls: Vec<String>,
    /// Literature citations for this guideline/drug report.
    #[serde(rename = "citations")]
    pub citations: BTreeSet<Publication>,
    /// Guideline reports.
    #[serde(rename = "guidelines")]
    pub guidelines: BTreeSet<GuidelineReport>,
}

impl DrugReport {
    fn new(
        name: String,
        guideline_packages: Vec<&GuidelinePackage>,
        report_context: &ReportContext,
    ) -> Self {
        assert!(!guideline_packages.is_empty());
        let first_guideline = &guideline_packages[0].guideline;
        let id = first_guideline
            .related_chemicals
            .iter()
            .find(|chemical| chemical.name == name)
            .unwrap_or_else(|| {
                panic!(
                    "DPWG guideline {} is supposed to be related to {name} but is not",
                    first_guideline.id
                )
            })
            .id
            .clone();
        let source = PrescribingGuidanceSource::type_for(first_guideline)
            .expect("prescribing guidance source");
        let mut urls = Vec::new();
        let mut citations = BTreeSet::new();
        let mut guidelines = BTreeSet::new();
        for guideline_package in guideline_packages {
            if let Some(url) = &guideline_package.url {
                urls.push(url.clone());
            }
            citations.extend(guideline_package.citations.iter().cloned());
            guidelines.insert(GuidelineReport::new(
                guideline_package,
                report_context,
                &name,
            ));
        }

        let mut drug_report = Self {
            name,
            id,
            source,
            version: report_context.data_version.clone(),
            messages: BTreeSet::new(),
            report_variants: BTreeSet::new(),
            urls,
            citations,
            guidelines,
        };
        drug_report.apply_related_drug_backlinks_to_guideline_genes();
        drug_report
    }

    fn apply_related_drug_backlinks_to_guideline_genes(&mut self) {
        let drug = DrugLink::new(self.name.clone(), self.id.clone());
        let guidelines = std::mem::take(&mut self.guidelines);
        self.guidelines = guidelines
            .into_iter()
            .map(|mut guideline| {
                for report_gene in &mut guideline.report_genes {
                    report_gene.add_related_drug(drug.clone());
                }
                guideline
            })
            .collect();
    }

    /// Java `isMatched`.
    pub fn is_matched(&self) -> bool {
        (self.source == PrescribingGuidanceSource::CpicGuideline && self.id == "RxNorm:11289")
            || self.guidelines.iter().any(GuidelineReport::is_matched)
    }

    /// Java `getMatchedAnnotationCount`.
    pub fn matched_annotation_count(&self) -> usize {
        self.guidelines
            .iter()
            .map(|guideline| guideline.annotations.len())
            .sum()
    }

    fn apply_matching_messages(
        &mut self,
        catalog: &MessageCatalog,
        gene_reports: &BTreeMap<String, ReportGene>,
        source: PrescribingGuidanceSource,
    ) {
        let messages = catalog.messages_for_drug(&self.name, source);
        if messages.is_empty() {
            return;
        }

        let mut report_as_genotype = Vec::new();
        for message in messages {
            if message.exception_type == MessageAnnotation::TYPE_REPORT_AS_GENOTYPE {
                self.report_variants
                    .extend(message.matches.variants.iter().cloned());
                report_as_genotype.push(message);
            } else if match_drug_report_message(message, gene_reports) {
                self.messages.insert(message.clone());
            }
        }

        for message in report_as_genotype {
            self.apply_report_as_genotype_message(message, gene_reports);
        }
    }

    fn apply_report_as_genotype_message(
        &mut self,
        message: &MessageAnnotation,
        gene_reports: &BTreeMap<String, ReportGene>,
    ) {
        let gene = message.matches.gene.as_deref();
        let genotype = compute_report_as_genotype(message, gene_reports);
        let guidelines = std::mem::take(&mut self.guidelines);
        self.guidelines = guidelines
            .into_iter()
            .map(|mut guideline| {
                if gene.is_none_or(|gene| guideline.genes.contains(gene)) {
                    let annotations = std::mem::take(&mut guideline.annotations);
                    guideline.annotations = annotations
                        .into_iter()
                        .map(|mut annotation| {
                            annotation.add_highlighted_variant(genotype.clone());
                            annotation
                        })
                        .collect();
                }
                guideline
            })
            .collect();
    }

    fn add_missing_variant_messages(&mut self, gene_reports: &BTreeMap<String, ReportGene>) {
        for gene in self.related_gene_symbols() {
            let Some(gene_report) = gene_reports.get(&gene) else {
                continue;
            };
            if gene_report.outside_call
                || !gene_report.is_missing_variants_like_java()
                || gene_report.is_no_data_like_java()
            {
                continue;
            }
            let gene_display = gene_display_name(&gene_report.gene);
            self.messages.insert(MessageAnnotation::new_note(
                "missing-variants",
                format!(
                    "Some position data used to define {gene_display} alleles is missing which may change the matched genotype. See <a href=\"#{gene_display}\">{gene_display}</a> in Section III for for more information."
                ),
            ));
        }
    }

    fn related_gene_symbols(&self) -> BTreeSet<String> {
        self.guidelines
            .iter()
            .flat_map(|guideline| guideline.genes.iter().cloned())
            .collect()
    }
}

impl Ord for DrugReport {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.name
            .cmp(&other.name)
            .then_with(|| self.source.cmp(&other.source))
            .then_with(|| self.version.cmp(&other.version))
            .then_with(|| self.id.cmp(&other.id))
            .then_with(|| self.guidelines.cmp(&other.guidelines))
            .then_with(|| self.messages.cmp(&other.messages))
    }
}

impl PartialOrd for DrugReport {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

/// First Rust `GuidelineReport` slice.
#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct GuidelineReport {
    /// Guideline id.
    #[serde(rename = "id")]
    pub id: String,
    /// Guideline name.
    #[serde(rename = "name")]
    pub name: String,
    /// Prescribing-guidance source.
    #[serde(rename = "source")]
    pub source: PrescribingGuidanceSource,
    /// Data version.
    #[serde(rename = "version")]
    pub version: String,
    /// URL.
    #[serde(rename = "url")]
    pub url: Option<String>,
    /// Genes from this guideline that are present in the report context.
    #[serde(skip)]
    pub genes: BTreeSet<String>,
    /// Java transient `GeneReport` links used by HTML helpers.
    #[serde(skip)]
    pub report_genes: Vec<ReportGene>,
    /// Matched annotation reports.
    #[serde(rename = "annotations")]
    pub annotations: BTreeSet<AnnotationReport>,
}

impl GuidelineReport {
    fn new(
        guideline_package: &GuidelinePackage,
        report_context: &ReportContext,
        drug_name: &str,
    ) -> Self {
        let source = PrescribingGuidanceSource::type_for(&guideline_package.guideline)
            .expect("prescribing guidance source");
        let report_genes = guideline_package
            .genes()
            .into_iter()
            .flat_map(|gene| report_context.report_genes_for_guideline_source(&gene, source))
            .collect::<Vec<_>>();
        let genes = report_genes
            .iter()
            .map(|report_gene| report_gene.gene.clone())
            .collect::<BTreeSet<_>>();
        let recommendation_genotypes = make_recommendation_genotypes(&genes, report_context);
        let mut annotations = BTreeSet::new();

        for recommendation in
            matching_recommendations(guideline_package, drug_name, &recommendation_genotypes)
        {
            let local_id = format!(
                "{}-{}",
                guideline_package.guideline.source, recommendation.id
            );
            annotations.insert(AnnotationReport::new(
                recommendation,
                local_id,
                &recommendation_genotypes,
                &guideline_package.guideline.related_alleles,
            ));
        }
        if drug_name == "warfarin" && source == PrescribingGuidanceSource::CpicGuideline {
            annotations.insert(AnnotationReport::for_cpic_warfarin(
                &recommendation_genotypes,
            ));
        }

        Self {
            id: guideline_package.guideline.id.clone(),
            name: guideline_package.guideline.name.clone(),
            source,
            version: report_context.data_version.clone(),
            url: guideline_package.url.clone(),
            genes,
            report_genes,
            annotations,
        }
    }

    /// Java `isMatched`.
    pub fn is_matched(&self) -> bool {
        !self.annotations.is_empty()
    }

    /// Java `isReportable`.
    pub fn is_reportable(&self) -> bool {
        self.report_genes
            .iter()
            .any(ReportGene::is_reportable_like_java)
    }

    fn homozygous_component_haplotypes(&self) -> BTreeSet<&str> {
        self.report_genes
            .iter()
            .flat_map(|report_gene| {
                report_gene
                    .matcher_homozygous_component_haplotypes
                    .iter()
                    .map(String::as_str)
            })
            .collect()
    }
}

impl Ord for GuidelineReport {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.name
            .cmp(&other.name)
            .then_with(|| self.source.cmp(&other.source))
            .then_with(|| self.version.cmp(&other.version))
            .then_with(|| self.id.cmp(&other.id))
            .then_with(|| self.annotations.cmp(&other.annotations))
    }
}

impl PartialOrd for GuidelineReport {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

/// First Rust `AnnotationReport` slice.
#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct AnnotationReport {
    /// Local report id.
    #[serde(skip)]
    pub local_id: String,
    /// Recommendation text HTML without wrapping paragraph tags.
    #[serde(rename = "drugRecommendation")]
    pub drug_recommendation: Option<String>,
    /// Recommendation classification.
    #[serde(rename = "classification")]
    pub classification: String,
    /// Population.
    #[serde(rename = "population")]
    pub population: String,
    /// Matched recommendation genotypes.
    #[serde(rename = "genotypes")]
    pub genotypes: Vec<RecommendationGenotype>,
    /// Implications.
    #[serde(rename = "implications")]
    pub implications: Vec<String>,
    /// Matched phenotype values by gene.
    #[serde(rename = "phenotypes")]
    pub phenotypes: BTreeMap<String, String>,
    /// Activity-score values by gene.
    #[serde(rename = "activityScore")]
    pub activity_scores: BTreeMap<String, String>,
    /// Highlighted variants that should be displayed as genotype text.
    #[serde(rename = "highlightedVariants")]
    pub highlighted_variants: BTreeSet<String>,
    /// Java `dosingInformation` flag used by genotype-summary drug tags.
    #[serde(rename = "dosingInformation")]
    pub dosing_information: bool,
    /// Java `alternateDrugAvailable` flag used by genotype-summary drug tags.
    #[serde(rename = "alternateDrugAvailable")]
    pub alternate_drug_available: bool,
    /// Java `otherPrescribingGuidance` flag used by genotype-summary drug tags.
    #[serde(rename = "otherPrescribingGuidance")]
    pub other_prescribing_guidance: bool,
    /// Annotation-level messages.
    #[serde(rename = "messages")]
    pub messages: BTreeSet<MessageAnnotation>,
    /// Lookup key maps from the recommendation.
    #[serde(rename = "lookupKey")]
    pub lookup_key: Vec<BTreeMap<String, Value>>,
}

impl AnnotationReport {
    fn new(
        recommendation: &RecommendationAnnotation,
        local_id: String,
        genotypes: &[RecommendationGenotype],
        alleles: &[AccessionObject],
    ) -> Self {
        let mut phenotypes = BTreeMap::new();
        let mut activity_scores = BTreeMap::new();
        let mut messages = BTreeSet::new();
        let mut matched_genotypes = Vec::new();
        for genotype in genotypes {
            if recommendation.matches_genotype(genotype)
                || recommendation.matches_diplotype(genotype)
            {
                matched_genotypes.push(genotype.clone());
                add_genotype_annotation(
                    genotype,
                    &recommendation.lookup_key,
                    alleles,
                    &mut phenotypes,
                    &mut activity_scores,
                    &mut messages,
                );
            }
        }

        Self {
            local_id,
            drug_recommendation: recommendation
                .text
                .as_ref()
                .map(|text| strip_paragraph_tags(&text.html)),
            classification: recommendation
                .classification
                .as_ref()
                .map(|classification| classification.term.clone())
                .filter(|term| {
                    !term.trim().is_empty()
                        && !term.trim().eq_ignore_ascii_case(crate::phenotype::NA)
                })
                .unwrap_or_else(|| "Unspecified".to_owned()),
            population: recommendation
                .population
                .clone()
                .unwrap_or_else(|| crate::phenotype::NA.to_owned()),
            genotypes: matched_genotypes,
            implications: recommendation.implications.clone(),
            phenotypes,
            activity_scores,
            highlighted_variants: BTreeSet::new(),
            dosing_information: recommendation.dosing_information,
            alternate_drug_available: recommendation.alternate_drug_available,
            other_prescribing_guidance: recommendation.other_prescribing_guidance,
            messages,
            lookup_key: recommendation.lookup_key.clone(),
        }
    }

    fn for_cpic_warfarin(genotypes: &[RecommendationGenotype]) -> Self {
        Self {
            local_id: "warfarin-cpic-1-1".to_owned(),
            drug_recommendation: None,
            classification: "Unspecified".to_owned(),
            population: crate::phenotype::NA.to_owned(),
            genotypes: genotypes.to_vec(),
            implications: Vec::new(),
            phenotypes: BTreeMap::new(),
            activity_scores: BTreeMap::new(),
            highlighted_variants: BTreeSet::new(),
            dosing_information: false,
            alternate_drug_available: false,
            other_prescribing_guidance: false,
            messages: BTreeSet::new(),
            lookup_key: Vec::new(),
        }
    }

    /// Java `addHighlightedVariant`.
    pub fn add_highlighted_variant(&mut self, variant: impl Into<String>) {
        self.highlighted_variants.insert(variant.into());
    }
}

impl Ord for AnnotationReport {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        compare_recommendation_genotypes_like_java(&self.genotypes, &other.genotypes)
            .then_with(|| self.population.cmp(&other.population))
            .then_with(|| self.highlighted_variants.cmp(&other.highlighted_variants))
            .then_with(|| self.activity_scores.cmp(&other.activity_scores))
            .then_with(|| self.classification.cmp(&other.classification))
            .then_with(|| self.drug_recommendation.cmp(&other.drug_recommendation))
            .then_with(|| self.implications.cmp(&other.implications))
            .then_with(|| self.messages.cmp(&other.messages))
            .then_with(|| self.local_id.cmp(&other.local_id))
    }
}

fn compare_recommendation_genotypes_like_java(
    left: &[RecommendationGenotype],
    right: &[RecommendationGenotype],
) -> std::cmp::Ordering {
    let mut left = left.iter().collect::<Vec<_>>();
    let mut right = right.iter().collect::<Vec<_>>();
    left.sort_by(|a, b| compare_recommendation_genotype_like_java(a, b));
    right.sort_by(|a, b| compare_recommendation_genotype_like_java(a, b));

    for (left, right) in left.iter().zip(&right) {
        let ordering = compare_recommendation_genotype_like_java(left, right);
        if ordering != std::cmp::Ordering::Equal {
            return ordering;
        }
    }
    left.len().cmp(&right.len())
}

fn compare_recommendation_genotype_like_java(
    left: &RecommendationGenotype,
    right: &RecommendationGenotype,
) -> std::cmp::Ordering {
    let mut left = left.report_genes.iter().collect::<Vec<_>>();
    let mut right = right.report_genes.iter().collect::<Vec<_>>();
    left.sort_by(|a, b| report_gene_source_cmp_like_java(a, b));
    right.sort_by(|a, b| report_gene_source_cmp_like_java(a, b));

    for (left, right) in left.iter().zip(&right) {
        let ordering = report_gene_source_cmp_like_java(left, right);
        if ordering != std::cmp::Ordering::Equal {
            return ordering;
        }
    }
    left.len().cmp(&right.len())
}

impl PartialOrd for AnnotationReport {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl PrescribingGuidanceSource {
    /// Java `typeFor`.
    pub fn type_for(guideline: &DosingGuideline) -> Option<Self> {
        Self::list_values()
            .into_iter()
            .find(|source| source.matches(guideline))
    }
}

fn make_recommendation_genotypes(
    genes: &BTreeSet<String>,
    report_context: &ReportContext,
) -> Vec<RecommendationGenotype> {
    let report_genes = genes
        .iter()
        .filter_map(|gene| report_context.gene_report(gene).cloned())
        .collect::<Vec<_>>();
    let mut possible_genes = Vec::<Vec<ReportGene>>::new();
    for report_gene in report_genes {
        let candidate_genes = if report_gene.recommendation_diplotypes.is_empty() {
            if report_gene.lookup_keys.is_empty() && !report_gene.diplotype_key.is_empty() {
                vec![report_gene.clone()]
            } else {
                report_gene
                    .lookup_keys
                    .iter()
                    .map(|lookup_key| {
                        report_gene_for_recommendation_lookup_key(&report_gene, lookup_key)
                    })
                    .collect::<Vec<_>>()
            }
        } else {
            report_gene
                .recommendation_diplotypes
                .iter()
                .map(|diplotype| report_gene_for_recommendation_genotype(&report_gene, diplotype))
                .collect::<Vec<_>>()
        };
        if possible_genes.is_empty() {
            possible_genes = candidate_genes
                .into_iter()
                .map(|report_gene| vec![report_gene])
                .collect();
        } else {
            let old_genes = std::mem::take(&mut possible_genes);
            possible_genes = candidate_genes
                .into_iter()
                .flat_map(|candidate| {
                    old_genes.iter().map(move |existing| {
                        let mut genes = existing.clone();
                        genes.push(candidate.clone());
                        genes
                    })
                })
                .collect();
        }
    }
    possible_genes
        .into_iter()
        .map(RecommendationGenotype::from_report_genes)
        .collect()
}

fn report_gene_for_recommendation_genotype(
    report_gene: &ReportGene,
    diplotype: &ReportDiplotype,
) -> ReportGene {
    let mut report_gene = report_gene.clone();
    report_gene.lookup_keys = if diplotype.lookup_keys.is_empty() {
        diplotype.phenotypes.clone()
    } else {
        diplotype.lookup_keys.clone()
    };
    report_gene.diplotype_key = diplotype.diplotype_key.clone();
    report_gene.phenotypes = diplotype.phenotypes.clone();
    report_gene.activity_score = diplotype.activity_score.clone();
    report_gene.source_diplotype = Some(diplotype.label.clone());
    report_gene.source_diplotypes = vec![diplotype.clone()];
    report_gene.recommendation_diplotypes = vec![diplotype.clone()];
    report_gene
}

fn report_gene_for_recommendation_lookup_key(
    report_gene: &ReportGene,
    lookup_key: &str,
) -> ReportGene {
    let mut report_gene = report_gene.clone();
    report_gene.lookup_keys = vec![lookup_key.to_owned()];
    if report_gene.is_activity_score_type {
        if report_gene.activity_score.is_none() {
            report_gene.activity_score = Some(lookup_key.to_owned());
        }
    } else {
        report_gene.phenotypes = vec![lookup_key.to_owned()];
    }
    report_gene
}

fn matching_recommendations<'a>(
    guideline_package: &'a GuidelinePackage,
    drug_name: &str,
    genotypes: &[RecommendationGenotype],
) -> Vec<&'a RecommendationAnnotation> {
    let mut matched = Vec::new();
    for genotype in genotypes {
        let mut matched_diplotype = false;
        for recommendation in &guideline_package.recommendations {
            if recommendation.applies_to_drug(drug_name)
                && recommendation.matches_diplotype(genotype)
            {
                if !matched
                    .iter()
                    .any(|existing: &&RecommendationAnnotation| existing.id == recommendation.id)
                {
                    matched.push(recommendation);
                }
                matched_diplotype = true;
            }
        }
        if !matched_diplotype {
            for recommendation in &guideline_package.recommendations {
                if recommendation.applies_to_drug(drug_name)
                    && recommendation.matches_genotype(genotype)
                    && !matched.iter().any(|existing: &&RecommendationAnnotation| {
                        existing.id == recommendation.id
                    })
                {
                    matched.push(recommendation);
                }
            }
        }
    }
    matched
}

fn strip_paragraph_tags(html: &str) -> String {
    let trimmed = html.trim();
    let stripped = trimmed.strip_prefix("<p>").unwrap_or(trimmed);
    stripped.strip_suffix("</p>").unwrap_or(stripped).to_owned()
}

fn add_genotype_annotation(
    genotype: &RecommendationGenotype,
    lookup_key: &[BTreeMap<String, Value>],
    alleles: &[AccessionObject],
    phenotypes: &mut BTreeMap<String, String>,
    activity_scores: &mut BTreeMap<String, String>,
    messages: &mut BTreeSet<MessageAnnotation>,
) {
    for report_gene in genotype.report_genes() {
        let gene = &report_gene.gene;
        if lookup_key.iter().all(|lookup| !lookup.contains_key(gene)) {
            continue;
        }

        add_outside_call_mismatch_message(report_gene, messages);

        if report_gene.is_allele_presence_type {
            let relevant_alleles = alleles
                .iter()
                .filter(|allele| {
                    allele
                        .symbol
                        .as_deref()
                        .is_some_and(|symbol| symbol.starts_with(gene))
                })
                .map(|allele| allele.name.as_str())
                .collect::<Vec<_>>();
            for phenotype in &report_gene.phenotypes {
                if relevant_alleles
                    .iter()
                    .any(|allele| phenotype.starts_with(*allele))
                {
                    phenotypes.insert(gene.clone(), phenotype.clone());
                }
            }
        } else {
            for phenotype in &report_gene.phenotypes {
                let old_phenotype = phenotypes.insert(gene.clone(), phenotype.clone());
                if old_phenotype
                    .as_deref()
                    .is_some_and(|old_phenotype| old_phenotype != phenotype)
                {
                    panic!("Multiple phenotypes for gene {gene}");
                }
            }
        }

        if genotype.uses_activity_score() {
            if report_gene.is_activity_score_type {
                let activity_score = report_gene
                    .activity_score
                    .clone()
                    .unwrap_or_else(|| crate::phenotype::NA.to_owned());
                let old_activity = activity_scores.insert(gene.clone(), activity_score.clone());
                if old_activity
                    .as_deref()
                    .is_some_and(|old_activity| old_activity != activity_score)
                {
                    panic!("Multiple activity scores for gene {gene}");
                }
            } else {
                activity_scores.insert(gene.clone(), crate::phenotype::NA.to_owned());
            }
        }
    }
}

fn add_outside_call_mismatch_message(
    report_gene: &ReportGene,
    messages: &mut BTreeSet<MessageAnnotation>,
) {
    if report_gene.outside_phenotype_mismatch.is_none()
        && report_gene.outside_activity_score_mismatch.is_none()
    {
        return;
    }

    let data_type = if report_gene.activity_score.is_some()
        || report_gene.outside_activity_score_mismatch.is_some()
    {
        "activity score"
    } else {
        "phenotype"
    };
    messages.insert(MessageAnnotation::new_note(
        "warn.mismatch.outsideCall",
        format!(
            "Conflicting outside call data was provided for {}.  PharmCAT will use provided {data_type} to match recommendations.",
            report_gene.gene
        ),
    ));
}

fn gene_call_warning_messages(result: &GeneCallResult) -> BTreeSet<MessageAnnotation> {
    result
        .warnings
        .iter()
        .map(|warning| gene_call_warning_message(&result.gene, warning))
        .collect()
}

fn gene_call_warning_message(gene: &str, warning: &GeneCallWarning) -> MessageAnnotation {
    match warning {
        GeneCallWarning::UnphasedPriority => MessageAnnotation::new_note(
            "unphased-priority",
            format!(
                "Unphased {gene} variants resulted in multiple calls.  PharmCAT is picking a single call based on frequency data.  Please consult the documentation for details."
            ),
        ),
        GeneCallWarning::MissingRequiredPosition(positions) => {
            let suffix = if positions.len() > 1 { "s" } else { "" };
            MessageAnnotation::new_note(
                "missing-required-position",
                format!(
                    "{gene} - missing required variant{suffix} ({})",
                    positions.join(", ")
                ),
            )
        }
        GeneCallWarning::MissingAmp1Position(positions) => {
            let suffix = if positions.len() > 1 { "s" } else { "" };
            MessageAnnotation::new_note(
                "missing-amp1-position",
                format!(
                    "Missing variant{suffix} required to meet AMP Tier 1 requirements:  {}. See https://www.clinpgx.org/ampAllelesToTest for details.",
                    positions.join(", ")
                ),
            )
        }
    }
}

fn compute_report_as_genotype(
    message: &MessageAnnotation,
    gene_reports: &BTreeMap<String, ReportGene>,
) -> String {
    let Some(gene) = message.matches.gene.as_deref() else {
        return String::new();
    };
    let gene_report = gene_reports.get(gene);
    message
        .matches
        .variants
        .iter()
        .map(|rsid| compute_single_report_as_genotype(gene_report, rsid))
        .collect::<Vec<_>>()
        .join(", ")
}

fn compute_single_report_as_genotype(gene_report: Option<&ReportGene>, rsid: &str) -> String {
    let call = gene_report
        .and_then(|gene_report| {
            gene_report
                .variant_reports
                .iter()
                .chain(gene_report.variant_of_interest_reports.iter())
                .find(|variant| variant.db_snp_id.as_deref() == Some(rsid) && !variant.is_missing())
        })
        .and_then(VariantReport::highlighted_call);

    format!(
        "{rsid}:{}",
        call.filter(|call| !call.trim().is_empty())
            .unwrap_or_else(|| "Unknown".to_owned())
    )
}

fn match_drug_report_message(
    message: &MessageAnnotation,
    gene_reports: &BTreeMap<String, ReportGene>,
) -> bool {
    let Some(gene) = message.matches.gene.as_deref() else {
        return true;
    };
    if gene.trim().is_empty() {
        return true;
    }
    gene_reports.get(gene).is_some_and(|gene_report| {
        !gene_report.is_no_data_like_java()
            && gene_report
                .messages
                .iter()
                .any(|gene_message| gene_message.name == message.name)
    })
}

fn gene_display_name(gene: &str) -> &str {
    match gene {
        "IFNL3" | "IFNL4" => "IFNL3/4",
        _ => gene,
    }
}

fn is_allele_presence_gene(gene: &str) -> bool {
    matches!(gene, "HLA-A" | "HLA-B")
}

fn is_variant_gene(gene: &str) -> bool {
    matches!(
        gene,
        "ABCG2" | "CACNA1S" | "CFTR" | "DPYD" | "G6PD" | "INFL3" | "MT-RNR1" | "RYR1" | "VKORC1"
    )
}

/// Java `RecommendationUtils.mapContains`.
pub fn map_contains(
    super_set: &BTreeMap<String, Value>,
    sub_set: &BTreeMap<String, Value>,
) -> bool {
    !super_set.is_empty()
        && !sub_set.is_empty()
        && super_set.len() >= sub_set.len()
        && sub_set
            .iter()
            .all(|(key, value)| super_set.get(key) == Some(value))
}

/// Java `AccessionObject`.
#[derive(Clone, Debug, Deserialize, Eq, PartialEq)]
pub struct AccessionObject {
    /// PharmGKB id.
    pub id: String,
    /// Name.
    pub name: String,
    /// Symbol.
    #[serde(default)]
    pub symbol: Option<String>,
}

/// Java `Publication`.
#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct Publication {
    /// PubMed id.
    #[serde(rename = "pmid", default)]
    pub pmid: Option<String>,
    /// Title.
    #[serde(default)]
    pub title: Option<String>,
    /// Journal.
    #[serde(default)]
    pub journal: Option<String>,
    /// Publication year.
    #[serde(default)]
    pub year: Option<i64>,
    /// Same-as URL.
    #[serde(rename = "_sameAs", default)]
    pub same_as: Option<String>,
}

impl Ord for Publication {
    fn cmp(&self, other: &Self) -> Ordering {
        self.year
            .cmp(&other.year)
            .then_with(|| self.pmid.cmp(&other.pmid))
            .then_with(|| self.title.cmp(&other.title))
    }
}

impl PartialOrd for Publication {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Java `OntologyTerm`.
#[derive(Clone, Debug, Deserialize, Eq, PartialEq)]
pub struct OntologyTerm {
    /// Term.
    pub term: String,
    /// Term id.
    #[serde(rename = "termId")]
    pub term_id: String,
}

/// Java `Markdown`.
#[derive(Clone, Debug, Deserialize, Eq, PartialEq)]
pub struct Markdown {
    /// Text id.
    pub id: i64,
    /// HTML.
    pub html: String,
}

/// Java `MatchLogic`.
#[derive(Clone, Debug, Default, Deserialize, Eq, Ord, PartialEq, PartialOrd, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct MatchLogic {
    /// Gene symbol.
    #[serde(default)]
    pub gene: Option<String>,
    /// Called haplotypes required by this rule.
    #[serde(default)]
    pub haps_called: Vec<String>,
    /// Missing haplotypes required by this rule.
    #[serde(default)]
    pub haps_missing: Vec<String>,
    /// Missing variants required by this rule.
    #[serde(default)]
    pub variants_missing: Vec<String>,
    /// Variants required by this rule.
    #[serde(rename = "variant", default)]
    pub variants: Vec<String>,
    /// Diplotypes required by this rule.
    #[serde(default)]
    pub dips: Vec<String>,
    /// Drugs this rule applies to.
    #[serde(default)]
    pub drugs: Vec<String>,
}

/// Java `MessageAnnotation`.
#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct MessageAnnotation {
    /// Rule name.
    #[serde(rename = "rule_name")]
    pub name: String,
    /// Rule version.
    #[serde(default)]
    pub version: Option<String>,
    /// Match logic.
    pub matches: MatchLogic,
    /// Exception type.
    #[serde(rename = "exception_type")]
    pub exception_type: String,
    /// Message HTML/text.
    pub message: String,
}

impl MessageAnnotation {
    /// Java `TYPE_AMBIGUITY`.
    pub const TYPE_AMBIGUITY: &'static str = "ambiguity";
    /// Java `TYPE_COMBO`.
    pub const TYPE_COMBO: &'static str = "combo-partial";
    /// Java `TYPE_EXTRA_POSITION`.
    pub const TYPE_EXTRA_POSITION: &'static str = "extra-position-notes";
    /// Java `TYPE_FOOTNOTE`.
    pub const TYPE_FOOTNOTE: &'static str = "footnote";
    /// Java `TYPE_NOTE`.
    pub const TYPE_NOTE: &'static str = "note";
    /// Java `TYPE_REPORT_AS_GENOTYPE`.
    pub const TYPE_REPORT_AS_GENOTYPE: &'static str = "report-as-genotype";
    /// Java `TYPE_NONMATCH`.
    pub const TYPE_NONMATCH: &'static str = "non-match";

    /// Creates a Java-style note message.
    pub fn new_note(name: impl Into<String>, message: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            version: None,
            matches: MatchLogic::default(),
            exception_type: Self::TYPE_NOTE.to_owned(),
            message: message.into(),
        }
    }

    /// Java `isFootnote`.
    pub fn is_footnote(&self) -> bool {
        self.exception_type == Self::TYPE_FOOTNOTE
    }

    /// Java `isExtraPositionNote`.
    pub fn is_extra_position_note(&self) -> bool {
        self.exception_type == Self::TYPE_EXTRA_POSITION
    }

    /// Java `isMessage`.
    pub fn is_message(&self) -> bool {
        self.exception_type != Self::TYPE_FOOTNOTE
            && self.exception_type != Self::TYPE_EXTRA_POSITION
            && self.exception_type != Self::TYPE_REPORT_AS_GENOTYPE
    }
}

impl Ord for MessageAnnotation {
    fn cmp(&self, other: &Self) -> Ordering {
        self.name
            .cmp(&other.name)
            .then_with(|| self.version.cmp(&other.version))
            .then_with(|| self.exception_type.cmp(&other.exception_type))
            .then_with(|| self.message.cmp(&other.message))
            .then_with(|| self.matches.cmp(&other.matches))
    }
}

impl PartialOrd for MessageAnnotation {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Java `MessageHelper` resource indexes.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct MessageCatalog {
    messages: Vec<MessageAnnotation>,
    gene_map: BTreeMap<String, Vec<MessageAnnotation>>,
    drug_map: BTreeMap<String, Vec<MessageAnnotation>>,
    static_map: BTreeMap<String, MessageAnnotation>,
}

impl MessageCatalog {
    /// Java `MESSAGES_JSON_FILE_NAME`.
    pub const MESSAGES_JSON_FILE_NAME: &'static str = "messages.json";

    /// Loads Java reporter `messages.json`.
    pub fn from_path(path: &Path) -> Result<Self, GuidanceLoadError> {
        let data = fs::read(path)?;
        let messages = serde_json::from_slice::<Vec<MessageAnnotation>>(&data)?;
        Ok(Self::from_messages(messages))
    }

    fn from_messages(messages: Vec<MessageAnnotation>) -> Self {
        let mut gene_map = BTreeMap::<String, Vec<MessageAnnotation>>::new();
        let mut drug_map = BTreeMap::<String, Vec<MessageAnnotation>>::new();
        let mut static_map = BTreeMap::<String, MessageAnnotation>::new();

        for message in &messages {
            if let Some(gene) = &message.matches.gene {
                gene_map
                    .entry(gene.clone())
                    .or_default()
                    .push(message.clone());
            }
            for drug in &message.matches.drugs {
                drug_map
                    .entry(drug.clone())
                    .or_default()
                    .push(message.clone());
            }
            if message.name.starts_with("pcat-") {
                static_map.insert(message.name.clone(), message.clone());
            }
        }

        Self {
            messages,
            gene_map,
            drug_map,
            static_map,
        }
    }

    /// All message annotations in resource order.
    pub fn messages(&self) -> &[MessageAnnotation] {
        &self.messages
    }

    /// Java static message lookup by key.
    pub fn message(&self, key: &str) -> Option<&MessageAnnotation> {
        self.static_map.get(key)
    }

    /// Messages indexed by gene.
    pub fn messages_for_gene(&self, gene: &str) -> &[MessageAnnotation] {
        self.gene_map.get(gene).map(Vec::as_slice).unwrap_or(&[])
    }

    /// Messages indexed by drug and allowed for a prescribing guidance source.
    pub fn messages_for_drug(
        &self,
        drug: &str,
        source: PrescribingGuidanceSource,
    ) -> Vec<&MessageAnnotation> {
        self.drug_map
            .get(drug)
            .into_iter()
            .flat_map(|messages| messages.iter())
            .filter(|message| message_allowed_for_source(message, source))
            .collect()
    }

    /// Java report-as-genotype messages indexed by drug and source.
    pub fn report_as_genotype_messages_for_drug(
        &self,
        drug: &str,
        source: PrescribingGuidanceSource,
    ) -> Vec<&MessageAnnotation> {
        self.messages_for_drug(drug, source)
            .into_iter()
            .filter(|message| message.exception_type == MessageAnnotation::TYPE_REPORT_AS_GENOTYPE)
            .collect()
    }

    /// Static `pcat-` message keys.
    pub fn static_message_keys(&self) -> BTreeSet<String> {
        self.static_map.keys().cloned().collect()
    }
}

fn message_allowed_for_source(
    message: &MessageAnnotation,
    source: PrescribingGuidanceSource,
) -> bool {
    let key = &message.name;
    if key.contains("cpic-") && source != PrescribingGuidanceSource::CpicGuideline {
        return false;
    }
    if key.contains("dpwg-") && source != PrescribingGuidanceSource::DpwgGuideline {
        return false;
    }
    !key.contains("fda-") || source == PrescribingGuidanceSource::FdaLabel
}

/// Loads Java reporter `disclaimers.hbs` as a template string.
pub fn load_disclaimers_template(path: &Path) -> Result<String, GuidanceLoadError> {
    Ok(fs::read_to_string(path)?)
}

/// Java reporter Handlebars template resources.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct HtmlTemplateSet {
    /// Java `report.hbs`.
    pub report: String,
    /// Java `header.hbs`.
    pub header: String,
    /// Java `uncallableGenesNote.hbs`.
    pub uncallable_genes_note: String,
    /// Java `disclaimers.hbs`.
    pub disclaimers: String,
}

impl HtmlTemplateSet {
    /// Loads Java reporter template resources from a reporter resource directory.
    pub fn from_reporter_dir(path: &Path) -> Result<Self, GuidanceLoadError> {
        Ok(Self {
            report: fs::read_to_string(path.join("report.hbs"))?,
            header: fs::read_to_string(path.join("header.hbs"))?,
            uncallable_genes_note: fs::read_to_string(path.join("uncallableGenesNote.hbs"))?,
            disclaimers: load_disclaimers_template(&path.join("disclaimers.hbs"))?,
        })
    }
}

/// Options for the first Rust HTML report writer.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct HtmlReportOptions {
    /// Optional PharmCAT version text.
    pub pharmcat_version: Option<String>,
    /// Optional data version text.
    pub data_version: Option<String>,
    /// Optional date-created text.
    pub timestamp: Option<String>,
    /// Whether to render Java debug-only report helpers.
    pub debug: bool,
    /// Java `HtmlFormat.compact` flag controlling compact report filtering.
    pub compact: bool,
    /// Java `DefinitionReader.getGenes()` set for deriving `HtmlFormat` noDataGenes.
    pub definition_genes: BTreeSet<String>,
    /// Additional Java `HtmlFormat` noDataGenes entries supplied to recommendation helpers.
    pub no_data_genes: BTreeSet<String>,
}

/// Generates a minimal HTML report for the currently ported report model.
pub fn report_html_string(context: &ReportContext, options: &HtmlReportOptions) -> String {
    let title_suffix = context
        .title
        .as_deref()
        .map(|title| format!(" [{}]", html_escape(title)))
        .unwrap_or_default();
    let data_version = options
        .data_version
        .as_deref()
        .unwrap_or(&context.data_version);
    let no_data_genes = html_no_data_genes(context, options);

    let mut html = String::new();
    html.push_str("<!DOCTYPE html>\n<html class=\"no-js\" lang=\"en\">\n<head>\n");
    html.push_str("  <meta charset=\"utf-8\" />\n");
    html.push_str("  <meta http-equiv=\"x-ua-compatible\" content=\"ie=edge\" />\n");
    html.push_str("  <title>PharmCAT Report");
    html.push_str(&title_suffix);
    html.push_str("</title>\n");
    html.push_str("</head>\n<body>\n<header class=\"pageHeader\">\n");
    html.push_str("  <div class=\"title\">PharmCAT Report");
    if let Some(title) = &context.title {
        html.push_str("<div class=\"subtitle\">");
        html.push_str(&html_escape(title));
        html.push_str("</div>");
    }
    html.push_str("</div>\n  <div class=\"metadata\"><table>\n");
    if let Some(timestamp) = &options.timestamp {
        html.push_str("    <tr><th>Date created</th><td>");
        html.push_str(&html_escape(timestamp));
        html.push_str("</td></tr>\n");
    }
    if let Some(version) = &options.pharmcat_version {
        html.push_str("    <tr><th>PharmCAT Version</th><td>");
        html.push_str(&html_escape(version));
        html.push_str("</td></tr>\n");
    }
    html.push_str("    <tr><th>Data Version</th><td>");
    html.push_str(&html_escape(data_version));
    html.push_str("</td></tr>\n  </table></div>\n</header>\n<main>\n");
    html.push_str(&html_genotype_summary_section(context, options));
    html.push_str(
        "  <section id=\"section-ii\">\n    <h2>Section II: Prescribing Recommendations</h2>\n",
    );
    let recommendation_drugs = html_recommendation_drugs(context, options);
    let drugs_without_recommendations = html_drugs_without_recommendations(context, options);
    if recommendation_drugs.is_empty() {
        html.push_str("    <p class=\"rx-no-recs\">No recommendations.</p>\n");
    }
    for drug in recommendation_drugs {
        html.push_str("    <section class=\"guideline drugReport ");
        html.push_str(&html_css_selector(&drug));
        html.push_str("\">\n      <h3 id=\"");
        html.push_str(&html_css_selector(&drug));
        html.push_str("\">");
        html.push_str(&html_escape(&drug));
        html.push_str("</h3>\n");

        let reports = html_reports_for_drug(context, &drug, options);
        let has_fda_recommendation = html_reports_have_matched_fda(&reports);
        for drug_report in &reports {
            for message in drug_report
                .messages
                .iter()
                .filter(|message| message.is_message())
            {
                html.push_str("      <div class=\"alert alert-info ");
                html.push_str(&html_escape(&html_message_class(message)));
                html.push_str("\">");
                html.push_str(&message.message);
                html.push_str("</div>\n");
            }
        }

        html.push_str("      <table><thead><tr><th>Source</th><th>Genes</th><th>Implications</th><th>Recommendation</th><th>Classification</th></tr></thead><tbody>\n");
        for drug_report in &reports {
            if drug_report.is_matched() {
                html.push_str(&html_matched_drug_report_rows(
                    drug_report,
                    options.debug,
                    &no_data_genes,
                    has_fda_recommendation,
                ));
            } else {
                html.push_str(&html_unmatched_drug_report_row(drug_report, &no_data_genes));
            }
        }
        html.push_str("      </tbody></table>\n");

        if html_reports_have_non_dpyd_inferred(&reports) {
            html.push_str("      <div class=\"footnote\" id=\"rx-dagger-");
            html.push_str(&html_css_selector(&drug));
            html.push_str("\"><sup>&dagger;</sup> Inferred genotype used to look up phenotype. Consult the <a href=\"https://pharmcat.org/methods/Gene-Definition-Exceptions/\" target=\"_blank\" rel=\"noopener noreferrer\">PharmCAT documentation</a> for details.</div>\n");
        }
        if html_reports_have_dpyd_inferred(&reports) {
            html.push_str("      <div class=\"footnote\" id=\"rx-ddagger-");
            html.push_str(&html_css_selector(&drug));
            html.push_str("\"><sup>&ddagger;</sup> The DPYD genotype used to look up phenotype is inferred from the two lowest function haplotypes. Consult the <a href=\"https://pharmcat.org/methods/Gene-Definition-Exceptions/#dpyd\" target=\"_blank\" rel=\"noopener noreferrer\">PharmCAT documentation</a> for details.</div>\n");
        }
        if has_fda_recommendation {
            html.push_str("      <div class=\"footnote\" id=\"rx-ast-");
            html.push_str(&html_css_selector(&drug));
            html.push_str("\"><sup>&ast;</sup> Text in quotation is taken directly from the FDA Label or FDA PGx Association table. For a label PDF with highlighted PGx content use the link in the Source column or go to <a href=\"https://www.clinpgx.org/fda\">ClinPGx</a>.</div>\n");
        }
        for footnote in html_report_footnotes(&reports) {
            html.push_str("      <div class=\"footnote ");
            html.push_str(&html_escape(&html_message_class(footnote)));
            html.push_str("\">");
            html.push_str(&footnote.message);
            html.push_str("</div>\n");
        }
        for drug_report in &reports {
            if !drug_report.citations.is_empty() {
                html.push_str("      <div class=\"citations\"><p>Citations:</p><ul>");
                for citation in &drug_report.citations {
                    html.push_str("<li>");
                    html.push_str(&html_publication_citation(citation));
                    html.push_str("</li>");
                }
                html.push_str("</ul></div>\n");
            }
        }
        html.push_str("    </section>\n");
    }
    if !drugs_without_recommendations.is_empty() {
        html.push_str("    <div>\n      <h3>Drugs With No Guidance</h3>\n");
        html.push_str("      <p>The following drugs are known to be associated with genes in this report but have no guidance for the specific genotypes in this report. For more information, see the <a href=\"https://www.clinpgx.org/prescribingInfo\">&quot;Prescribing Info&quot; page on ClinPGx</a>.</p>\n");
        html.push_str("      <ul>\n");
        for drug in drugs_without_recommendations {
            html.push_str("        <li>");
            html.push_str(&html_escape(&drug));
            html.push_str("</li>\n");
        }
        html.push_str("      </ul>\n    </div>\n");
    }
    html.push_str("  </section>\n");
    let amd_section = html_amd_section(context, options, &no_data_genes);
    if !amd_section.is_empty() {
        html.push_str(&amd_section);
    }
    html.push_str("  <section id=\"disclaimer\"><h2>Section IV: Disclaimers and Other Information</h2></section>\n");
    html.push_str("</main>\n</body>\n</html>\n");
    html
}

fn html_no_data_genes(context: &ReportContext, options: &HtmlReportOptions) -> BTreeSet<String> {
    let mut genes = options.no_data_genes.clone();
    for report_gene in context.gene_reports.values() {
        if report_gene.is_no_data_like_java()
            && options.definition_genes.contains(&report_gene.gene)
        {
            genes.insert(gene_display_name(&report_gene.gene).to_owned());
        }
    }
    genes
}

fn html_genotype_summary_section(context: &ReportContext, options: &HtmlReportOptions) -> String {
    let drugs_with_recommendations = html_drugs_with_recommendations(context, options);
    let drug_tags = html_drug_tags(context, options);
    let summary_genes = html_genotype_summary_report_genes(context, options);
    let uncallable_genes = html_uncallable_genes(context);
    let mut html = String::new();
    html.push_str("  <section id=\"section-i\">\n    <h2>Section I: Genotype Summary</h2>\n");
    if summary_genes.is_empty() {
        if uncallable_genes.is_empty() {
            html.push_str("    <p>No data provided.</p>\n");
        } else {
            html.push_str("    <p>No genotypes called.</p>\n");
            html.push_str(&html_uncallable_genes_note(&uncallable_genes));
        }
        html.push_str("    <p>For a full list of disclaimers and limitations, see <a href=\"#disclaimer\">Section IV</a>.</p>\n");
        html.push_str("  </section>\n");
        return html;
    }

    html.push_str("    <p>Genotypes called: ");
    html.push_str(&html_genotype_summary_called_genes(context).to_string());
    html.push_str(" / ");
    html.push_str(&html_genotype_summary_total_genes(context, options).to_string());
    html.push_str(" </p>\n");
    html.push_str("    <table class=\"genotypeSummary\"><thead><tr><th>Drugs</th><th>Gene</th><th>Genotypes<table class=\"diplotype table-small\"><tbody><tr><td>Genotype</td><td>Allele Functionality</td><td>Phenotype</td></tr></tbody></table></th></tr></thead><tbody>\n");
    for report_gene in &summary_genes {
        let gene_display = gene_display_name(&report_gene.gene);
        html.push_str("      <tr class=\"top-aligned gs-");
        html.push_str(&html_java_css_selector(gene_display));
        html.push_str("\"><td>");
        html.push_str(&html_genotype_summary_drugs(
            report_gene,
            &drugs_with_recommendations,
            &drug_tags,
        ));
        html.push_str("</td><td><span class=\"noWrap\"><a href=\"#");
        html.push_str(&html_java_css_selector(gene_display));
        html.push_str("\" class=\"normalWrap\">");
        html.push_str(&html_escape(gene_display));
        html.push_str("</a>");
        if html_gene_summary_has_messages(report_gene) {
            html.push_str("<sup><a href=\"#genotypes-dagger\">&dagger;</a></sup>");
        }
        if html_amd_show_unphased_note(report_gene) {
            html.push_str("<sup><a href=\"#genotypes-ddagger\">&ddagger;</a></sup>");
        }
        html.push_str("</span></td><td>");
        html.push_str(&html_genotype_summary_diplotypes(report_gene));
        html.push_str("</td></tr>\n");
    }
    html.push_str("    </tbody></table>\n");
    html.push_str(&html_uncallable_genes_warning(&uncallable_genes));
    html.push_str(&html_genotype_summary_footnotes(context));
    html.push_str("    <div class=\"footnote\">CPIC terms for allele function and phenotype are used for all CPIC genes. For non-CPIC genes, DPWG terms are used.</div>\n");
    html.push_str("    <div class=\"footnote\">For a full list of disclaimers and limitations see <a href=\"#disclaimer\">Section IV</a>.</div>\n");
    html.push_str(&html_genotype_summary_combo_alert(context));
    html.push_str(&html_genotype_summary_messages(context));
    html.push_str("  </section>\n");
    html
}

fn html_genotype_summary_report_genes(
    context: &ReportContext,
    options: &HtmlReportOptions,
) -> Vec<ReportGene> {
    context
        .report_gene_sources
        .values()
        .filter_map(|report_genes| html_genotype_summary_report_gene(report_genes))
        .filter(|report_gene| !options.compact || !report_gene.related_drugs.is_empty())
        .collect()
}

fn html_genotype_summary_report_gene(report_genes: &[ReportGene]) -> Option<ReportGene> {
    let mut summary = None;
    let mut related_drugs = BTreeSet::new();
    let mut messages = BTreeSet::new();
    for report_gene in report_genes
        .iter()
        .filter(|report_gene| report_gene.is_reportable_like_java())
    {
        related_drugs.extend(report_gene.related_drugs.iter().cloned());
        messages.extend(report_gene.messages.iter().cloned());
        summary.get_or_insert_with(|| report_gene.clone());
    }

    summary.map(|mut report_gene: ReportGene| {
        report_gene.related_drugs = related_drugs;
        report_gene.messages = messages;
        report_gene
    })
}

fn html_genotype_summary_total_genes(
    context: &ReportContext,
    options: &HtmlReportOptions,
) -> usize {
    let mut genes = context
        .gene_reports
        .values()
        .filter(|report_gene| !report_gene.is_no_data_like_java())
        .map(|report_gene| gene_display_name(&report_gene.gene).to_owned())
        .collect::<BTreeSet<_>>();
    genes.extend(html_no_data_genes(context, options));
    genes.len()
}

fn html_genotype_summary_called_genes(context: &ReportContext) -> usize {
    let related_genes = html_related_gene_symbols(context);
    context
        .report_gene_sources
        .values()
        .filter_map(|report_genes| html_genotype_summary_report_gene(report_genes))
        .filter(|report_gene| {
            !report_gene.is_no_data_like_java()
                && (!report_gene.related_drugs.is_empty()
                    || related_genes.contains(&report_gene.gene))
        })
        .map(|report_gene| gene_display_name(&report_gene.gene).to_owned())
        .collect::<BTreeSet<_>>()
        .len()
}

fn html_uncallable_genes(context: &ReportContext) -> BTreeSet<String> {
    context
        .gene_reports
        .values()
        .filter(|report_gene| !report_gene.is_no_data_like_java())
        .filter(|report_gene| !report_gene.related_drugs.is_empty())
        .filter(|report_gene| !report_gene.is_reportable_like_java())
        .map(|report_gene| gene_display_name(&report_gene.gene).to_owned())
        .collect()
}

fn html_uncallable_genes_warning(genes: &BTreeSet<String>) -> String {
    if genes.is_empty() {
        return String::new();
    }

    let mut html = String::new();
    html.push_str("    <div class=\"alert alert-warning\">");
    html.push_str(&html_uncallable_genes_note(genes));
    html.push_str("</div>\n");
    html
}

fn html_uncallable_genes_note(genes: &BTreeSet<String>) -> String {
    let mut html = String::new();
    html.push_str("<p>The following ");
    html.push_str(if genes.len() == 1 { "gene" } else { "genes" });
    html.push_str(" could not be called because there were genetic variations that do not match the allele definition. There could still be actionable variants in these genes. See <a href=\"#section-iii\">Section III</a> for details.</p><ul class=\"mb-2\">");
    for gene in genes {
        html.push_str("<li><a id=\"gs-uncallable-");
        html.push_str(&html_java_css_selector(gene));
        html.push_str("\" href=\"#");
        html.push_str(&html_java_css_selector(gene));
        html.push_str("\">");
        html.push_str(&html_escape(gene));
        html.push_str("</a></li>");
    }
    html.push_str("</ul>");
    html
}

fn html_genotype_summary_combo_alert(context: &ReportContext) -> String {
    let has_combo = html_genotype_summary_report_genes(context, &HtmlReportOptions::default())
        .into_iter()
        .any(|report_gene| {
            report_gene
                .messages
                .iter()
                .any(|message| message.exception_type == MessageAnnotation::TYPE_COMBO)
        });
    if !has_combo {
        return String::new();
    }

    "    <div class=\"alert alert-info\">Partial and combination allele calls are based on the variants identified in the VCF file. Matches to different star allele definitions or star allele with additional position combinations are connected by a '+', which represents combinations in a single gene copy per allele and does NOT indicate gene duplications.</div>\n".to_owned()
}

fn html_genotype_summary_messages(context: &ReportContext) -> String {
    let mut html = String::new();
    for message in context
        .messages
        .iter()
        .filter(|message| message.is_message())
    {
        html.push_str("    <div class=\"alert alert-info ");
        html.push_str(&html_message_class(message));
        html.push_str("\">");
        html.push_str(&html_escape(&message.message));
        html.push_str("</div>\n");
    }
    html
}

fn html_genotype_summary_diplotypes(report_gene: &ReportGene) -> String {
    let diplotypes = html_genotype_summary_report_diplotypes(report_gene);
    let show_components = is_lowest_function_gene(&report_gene.gene)
        && !report_gene.matcher_component_haplotypes.is_empty();
    let mut html = String::new();
    html.push_str("<table class=\"diplotype\"><tbody>");
    if diplotypes.is_empty() {
        if show_components && !report_gene.source_diplotypes.is_empty() {
            let homozygous_component_haplotypes = report_gene
                .matcher_homozygous_component_haplotypes
                .iter()
                .map(String::as_str)
                .collect::<BTreeSet<_>>();
            for component in &report_gene.source_diplotypes {
                html.push_str(&html_genotype_summary_component_row(
                    component,
                    &homozygous_component_haplotypes,
                    None,
                    false,
                ));
            }
        } else {
            html.push_str("<tr class=\"top-aligned gs-dip\"><td>");
            html.push_str(&html_escape(
                report_gene
                    .source_diplotype
                    .as_deref()
                    .or_else(|| report_gene.lookup_keys.first().map(String::as_str))
                    .unwrap_or(""),
            ));
            html.push_str("</td><td>");
            html.push_str(&crate::phenotype::NA.to_uppercase());
            html.push_str("</td><td>");
            html.push_str(&html_genotype_summary_fallback_phenotype(report_gene));
            html.push_str("</td></tr>");
        }
    } else {
        let homozygous_component_haplotypes = report_gene
            .matcher_homozygous_component_haplotypes
            .iter()
            .map(String::as_str)
            .collect::<BTreeSet<_>>();
        if show_components {
            for diplotype in diplotypes {
                html.push_str("<tr class=\"lowestFunctionDiplotype\"><td colspan=\"3\"><b class=\"gs-dip_lowestFunction\">");
                html.push_str(&html_gs_call(diplotype, &homozygous_component_haplotypes));
                html.push_str("</b><br /></td></tr>");
                html.push_str(&html_gs_lowest_function_components(
                    diplotype,
                    &report_gene.source_diplotypes,
                    &homozygous_component_haplotypes,
                ));
            }
        } else {
            for diplotype in diplotypes {
                html.push_str("<tr class=\"top-aligned gs-dip\"><td>");
                html.push_str(&html_gs_call(diplotype, &homozygous_component_haplotypes));
                html.push_str("</td><td>");
                html.push_str(&html_escape(&html_gs_function(diplotype)));
                html.push_str("</td><td>");
                html.push_str(&html_escape(&html_gs_phenotype(diplotype)));
                html.push_str("</td></tr>");
            }
        }
    }
    html.push_str("</tbody></table>");
    if report_gene.is_missing_variants_like_java() {
        html.push_str("<p class=\"tdNote\">Genotype");
        if html_genotype_summary_report_diplotypes(report_gene).len() > 1 {
            html.push('s');
        }
        html.push_str(
            " based on missing variant input <sup><a href=\"#genotypes-star\">*</a></sup>.</p>",
        );
    }
    if report_gene.treat_undocumented_variations_as_reference {
        let gene_display = gene_display_name(&report_gene.gene);
        html.push_str("<p class=\"tdNote\" id=\"gs-undocVarAsRef-");
        html.push_str(&html_java_css_selector(gene_display));
        html.push_str("\">There are genetic variations in this gene that do not match what is in the allele definition. <b>These undocumented variations were replaced with reference.</b> See <a href=\"#");
        html.push_str(&html_java_css_selector(gene_display));
        html.push_str("\">Section III</a> for details.</p>");
    }
    html
}

fn html_genotype_summary_footnotes(context: &ReportContext) -> String {
    let mut html = String::new();
    let summary_genes = html_genotype_summary_report_genes(context, &HtmlReportOptions::default());
    if summary_genes
        .iter()
        .any(ReportGene::is_missing_variants_like_java)
    {
        html.push_str("    <div class=\"footnote\" id=\"genotypes-star\"><sup>*</sup> Some alleles were not considered for the genotype calls due to missing variant information. Please see <a href=\"#section-iii\">Section III</a> for details. Alleles that could not be considered due to missing input might change the metabolizer phenotype and possible recommendation.</div>\n");
    }
    if summary_genes.iter().any(html_gene_summary_has_messages) {
        html.push_str("    <div class=\"footnote\" id=\"genotypes-dagger\"><sup>&dagger;</sup> Check <a href=\"#section-iii\">Section III</a> for more details about this call.</div>\n");
    }
    if summary_genes.iter().any(html_amd_show_unphased_note) {
        html.push_str("    <div class=\"footnote\" id=\"genotypes-ddagger\"><sup>&ddagger;</sup> PharmCAT reports the genotype(s) that receive the highest score during the matcher process. In case of unphased data, additional genotypes might be possible and cannot be ruled out.</div>\n");
    }
    html
}

fn html_genotype_summary_report_diplotypes(report_gene: &ReportGene) -> Vec<&ReportDiplotype> {
    if !report_gene.source_diplotypes.is_empty()
        && !is_lowest_function_gene(&report_gene.gene)
        && report_gene.gene != "SLCO1B1"
    {
        report_gene.source_diplotypes.iter().collect()
    } else {
        report_gene.recommendation_diplotypes.iter().collect()
    }
}

fn html_genotype_summary_fallback_phenotype(report_gene: &ReportGene) -> String {
    if report_gene.phenotypes.is_empty() {
        crate::phenotype::NA.to_uppercase()
    } else {
        html_escape(&report_gene.phenotypes.join("; "))
    }
}

fn html_gs_call(
    diplotype: &ReportDiplotype,
    homozygous_component_haplotypes: &BTreeSet<&str>,
) -> String {
    if diplotype.is_unknown() {
        let call = if diplotype.outside_phenotype || diplotype.outside_activity_score {
            "Not provided"
        } else {
            "Not called"
        };
        return format!(
            "<span class=\"gs-uncalled-{}\">{}</span>",
            html_escape(&diplotype.gene),
            call
        );
    }
    let display_diplotype = match diplotype.inferred_source_diplotypes.as_slice() {
        [source] => source,
        _ => diplotype,
    };
    let mut label = display_diplotype.label.clone();
    if homozygous_component_haplotypes.contains(label.as_str()) {
        label.push_str(" (homozygous)");
    }
    html_escape(&label)
}

fn html_gs_lowest_function_components(
    diplotype: &ReportDiplotype,
    components: &[ReportDiplotype],
    homozygous_component_haplotypes: &BTreeSet<&str>,
) -> String {
    let haps = html_gs_diplotype_haplotype_names(diplotype);
    let matching_components = components
        .iter()
        .filter(|component| {
            component
                .allele1
                .as_ref()
                .is_some_and(|allele| haps.contains(&allele.name))
        })
        .collect::<Vec<_>>();
    let mut html = String::new();
    for (index, component) in matching_components.iter().enumerate() {
        html.push_str(&html_genotype_summary_component_row(
            component,
            homozygous_component_haplotypes,
            (index == 0).then_some(components.len()),
            index != 0,
        ));
    }
    html
}

fn html_genotype_summary_component_row(
    component: &ReportDiplotype,
    homozygous_component_haplotypes: &BTreeSet<&str>,
    phenotype_rowspan: Option<usize>,
    omit_phenotype: bool,
) -> String {
    let mut html = String::new();
    html.push_str("<tr class=\"top-aligned gs-dip_component\"><td>");
    html.push_str(&html_gs_call(component, homozygous_component_haplotypes));
    html.push_str("</td><td>");
    html.push_str(&html_escape(&html_gs_function(component)));
    html.push_str("</td>");
    if omit_phenotype {
        // Java's first component row spans the phenotype cell; later component rows omit it.
    } else if let Some(rowspan) = phenotype_rowspan {
        html.push_str("<td rowspan=\"");
        html.push_str(&rowspan.to_string());
        html.push_str("\" class=\"center\">See Drug Recommendation</td>");
    } else {
        html.push_str("<td>");
        html.push_str(&html_escape(&html_gs_phenotype(component)));
        html.push_str("</td>");
    }
    html.push_str("</tr>");
    html
}

fn html_gs_diplotype_haplotype_names(diplotype: &ReportDiplotype) -> BTreeSet<String> {
    let mut names = BTreeSet::new();
    for allele in [&diplotype.allele1, &diplotype.allele2]
        .into_iter()
        .flatten()
    {
        names.extend(html_gs_parse_haplotype_name(&allele.name));
    }
    names
}

fn html_gs_parse_haplotype_name(name: &str) -> Vec<String> {
    if is_combination_label(name) {
        name.trim_start_matches('[')
            .trim_end_matches(']')
            .split(" + ")
            .map(str::to_owned)
            .collect()
    } else {
        vec![name.to_owned()]
    }
}

fn html_gs_function(diplotype: &ReportDiplotype) -> String {
    if diplotype.combination && is_lowest_function_gene(&diplotype.gene) {
        return "See Drug Recommendation".to_owned();
    }
    let function1 = diplotype
        .allele1
        .as_ref()
        .map(|allele| allele.function.as_str())
        .filter(|function| !function.trim().is_empty());
    let function2 = diplotype
        .allele2
        .as_ref()
        .map(|allele| allele.function.as_str())
        .filter(|function| !function.trim().is_empty());

    match (function1, function2) {
        (None, None) => crate::phenotype::NA.to_uppercase(),
        (Some(function), None) | (None, Some(function)) => function.to_owned(),
        (Some(function1), Some(function2)) if function1 == function2 => {
            format!("Two {function1} alleles")
        }
        (Some(function1), Some(function2)) => {
            let mut functions = [function1, function2];
            functions.sort();
            format!(
                "One {} allele and one {} allele",
                functions[0], functions[1]
            )
        }
    }
}

fn html_gs_phenotype(diplotype: &ReportDiplotype) -> String {
    if is_lowest_function_gene(&diplotype.gene) {
        return "See Drug Recommendation".to_owned();
    }
    if diplotype.phenotypes.is_empty() {
        crate::phenotype::NA.to_uppercase()
    } else {
        diplotype.phenotypes.join("; ")
    }
}

fn html_genotype_summary_drugs(
    report_gene: &ReportGene,
    drugs_with_recommendations: &BTreeSet<String>,
    drug_tags: &BTreeMap<String, BTreeSet<&'static str>>,
) -> String {
    let mut html = String::new();
    for drug in report_gene
        .related_drugs
        .iter()
        .map(|drug| drug.name.as_str())
        .filter(|drug| drugs_with_recommendations.contains(*drug))
    {
        html.push_str("<div class=\"gsDrugs\"><span class=\"drugName\"><a href=\"#");
        html.push_str(&html_java_css_selector(drug));
        html.push_str("\">");
        html.push_str(&html_escape(drug));
        html.push_str("</a></span><span class=\"drugTags\">");
        html.push_str(&html_drug_tag_markup(drug, drug_tags));
        html.push_str("</span></div>");
    }
    html
}

fn html_drugs_with_recommendations(
    context: &ReportContext,
    options: &HtmlReportOptions,
) -> BTreeSet<String> {
    context
        .drug_reports
        .values()
        .flat_map(|source_reports| source_reports.values())
        .filter(|drug_report| html_drug_report_has_recommendations(drug_report, options))
        .map(|drug_report| drug_report.name.clone())
        .collect()
}

fn html_drug_report_visible_for_recommendations(
    drug_report: &DrugReport,
    options: &HtmlReportOptions,
) -> bool {
    !options.compact
        || drug_report
            .guidelines
            .iter()
            .any(GuidelineReport::is_reportable)
}

fn html_drug_tags(
    context: &ReportContext,
    options: &HtmlReportOptions,
) -> BTreeMap<String, BTreeSet<&'static str>> {
    if !options.compact {
        return BTreeMap::new();
    }

    let mut drug_tags = BTreeMap::new();
    for source_reports in context.drug_reports.values() {
        for drug_report in source_reports.values() {
            if !drug_report.is_matched()
                || !html_drug_report_visible_for_recommendations(drug_report, options)
            {
                continue;
            }
            if drug_report.name == "warfarin"
                && drug_report.source == PrescribingGuidanceSource::CpicGuideline
            {
                continue;
            }
            let tags = drug_tags
                .entry(drug_report.name.clone())
                .or_insert_with(BTreeSet::new);
            for annotation in drug_report
                .guidelines
                .iter()
                .flat_map(|guideline| guideline.annotations.iter())
            {
                if annotation.alternate_drug_available {
                    tags.insert("Alternate Drug");
                }
                if annotation.dosing_information {
                    tags.insert("Dosing Info");
                }
                if annotation.other_prescribing_guidance {
                    tags.insert("Other Guidance");
                }
            }
        }
    }
    drug_tags
}

fn html_drug_tag_markup(
    drug: &str,
    drug_tags: &BTreeMap<String, BTreeSet<&'static str>>,
) -> String {
    let Some(tags) = drug_tags.get(drug) else {
        return String::new();
    };
    tags.iter()
        .map(|tag| format!("<div class=\"tag\">{}</div>", html_escape(tag)))
        .collect::<Vec<_>>()
        .join("\n")
}

fn html_gene_summary_has_messages(report_gene: &ReportGene) -> bool {
    report_gene
        .messages
        .iter()
        .any(|message| message.is_message())
}

fn html_amd_no_data_message(no_data_genes: &BTreeSet<String>) -> String {
    let genes = no_data_genes
        .iter()
        .map(|gene| {
            format!(
                "<span class=\"gene {}\"><span class=\"no-data\">{}</span></span>",
                html_css_selector(&gene.to_lowercase()),
                html_escape(gene)
            )
        })
        .collect::<Vec<_>>()
        .join(", ");
    format!("<p class=\"noGeneData\">No data provided for {genes}.</p>")
}

fn html_amd_section(
    context: &ReportContext,
    options: &HtmlReportOptions,
    no_data_genes: &BTreeSet<String>,
) -> String {
    let gene_reports = html_amd_gene_reports(context, options, no_data_genes);
    if gene_reports.is_empty() && no_data_genes.is_empty() {
        return String::new();
    }

    let mut html = String::new();
    html.push_str(
        "  <section id=\"section-iii\">\n    <h2>Section III: Allele Matching Details</h2>\n",
    );
    if !gene_reports.is_empty() {
        html.push_str("\n    <ol>\n");
        for report_gene in &gene_reports {
            let gene_display = gene_display_name(&report_gene.gene);
            html.push_str("      <li><a href=\"#");
            html.push_str(&html_java_css_selector(gene_display));
            html.push_str("\">");
            html.push_str(&html_escape(gene_display));
            html.push_str(" allele match data</a></li>\n");
        }
        html.push_str("    </ol>\n");
    }
    if !no_data_genes.is_empty() {
        html.push_str("    ");
        html.push_str(&html_amd_no_data_message(no_data_genes));
        html.push('\n');
    }
    for report_gene in gene_reports {
        html.push_str(&html_amd_gene_report(report_gene));
    }
    html.push_str("  </section>\n");
    html
}

fn html_amd_gene_reports<'a>(
    context: &'a ReportContext,
    options: &HtmlReportOptions,
    no_data_genes: &BTreeSet<String>,
) -> Vec<&'a ReportGene> {
    let related_genes = html_related_gene_symbols(context);
    context
        .gene_reports
        .values()
        .filter(|report_gene| !no_data_genes.contains(gene_display_name(&report_gene.gene)))
        .filter(|report_gene| !report_gene.is_no_data_like_java())
        .filter(|report_gene| {
            if !options.compact {
                return true;
            }
            related_genes.contains(&report_gene.gene)
        })
        .collect()
}

fn html_related_gene_symbols(context: &ReportContext) -> BTreeSet<String> {
    context
        .drug_reports
        .values()
        .flat_map(|source_reports| source_reports.values())
        .flat_map(|drug_report| drug_report.guidelines.iter())
        .flat_map(|guideline| guideline.report_genes.iter())
        .map(|report_gene| report_gene.gene.clone())
        .collect()
}

fn html_amd_gene_report(report_gene: &ReportGene) -> String {
    let gene_display = gene_display_name(&report_gene.gene);
    let gene_selector = html_java_css_selector(gene_display);
    let mut html = String::new();
    html.push_str("    <section class=\"gene ");
    html.push_str(&gene_selector);
    html.push_str("\">\n      <h3 id=\"");
    html.push_str(&gene_selector);
    html.push_str("\">");
    html.push_str(&html_escape(gene_display));
    html.push_str(" allele match data</h3>\n\n");

    let no_call = html_amd_no_call(report_gene);
    if no_call {
        html.push_str("      <div class=\"alert alert-warning no-data\">\n        ");
        html.push_str(&html_escape(&html_amd_gene_call(report_gene)));
        html.push_str(".\n      </div>\n");
    } else {
        html.push_str(
            "      <table>\n        <tbody>\n        <tr>\n          <th style=\"width: 12em;\">",
        );
        html.push_str(&html_escape(&html_amd_subtitle(report_gene)));
        html.push_str(":</th>\n          <td class=\"top-aligned genotype-result\">");
        let calls = html_amd_gene_calls(report_gene);
        if calls.len() == 1 {
            html.push_str(&html_escape(&calls[0]));
        } else {
            html.push_str("\n            <ul class=\"noPadding mt-0\">\n");
            for call in calls {
                html.push_str("              <li>");
                html.push_str(&html_escape(&call));
                html.push_str("</li>\n");
            }
            html.push_str("            </ul>\n          ");
        }
        if report_gene.treat_undocumented_variations_as_reference {
            html.push_str("\n            <p class=\"tdNote\">\n              There are genetic variations in this gene that do not match what is in the allele definition.\n              <b>These undocumented variations were replaced with reference.</b>  See below for details.\n            </p>");
        }
        html.push_str("</td>\n        </tr>\n");

        if !report_gene.variant_reports.is_empty() {
            html.push_str("        <tr>\n          <th>Phasing Status:</th>\n          <td class=\"top-aligned\">\n            <p>");
            html.push_str(&html_escape(&html_amd_phase_status(report_gene)));
            html.push_str("</p>\n");
            if html_amd_show_unphased_note(report_gene) {
                html.push_str("            <p>PharmCAT reports the genotype(s) that receive the highest score during the matcher process. In case of unphased data, additional genotypes might be possible and cannot be ruled out.</p>\n");
            }
            html.push_str("          </td>\n        </tr>\n");
        }

        if !report_gene.uncalled_haplotypes.is_empty() {
            html.push_str("        <tr>\n          <th>Alleles Not Considered:</th>\n          <td class=\"top-aligned\">\n            <p>The following alleles are not considered due to ");
            html.push_str(&html_amd_total_missing_variants(report_gene).to_string());
            html.push_str(" missing positions of the total ");
            html.push_str(&report_gene.variant_reports.len().to_string());
            html.push_str(" positions: ");
            html.push_str(&html_escape(&html_amd_uncalled_haps(report_gene)));
            html.push_str("</p>\n            <p>Carriage of these alleles might result in a different phenotype and different guideline recommendations.</p>\n          </td>\n        </tr>\n");
        }
        html.push_str("        </tbody>\n      </table>\n");
    }

    for message in report_gene
        .messages
        .iter()
        .filter(|message| message.is_message())
    {
        html.push_str("\n      <div class=\"alert alert-warning ");
        html.push_str(&html_message_class(message));
        html.push_str("\">");
        html.push_str(&html_message_message(message));
        html.push_str("</div>\n");
    }

    if !no_call {
        if !report_gene.variant_reports.is_empty() {
            html.push_str(&html_amd_calls_at_positions(report_gene));
        }
        if !report_gene.variant_of_interest_reports.is_empty() {
            html.push_str(&html_amd_other_positions_of_interest(report_gene));
        }
    }

    html.push_str("    </section>\n");
    html
}

fn html_amd_subtitle(report_gene: &ReportGene) -> String {
    let mut title = if is_variant_gene(&report_gene.gene)
        && report_gene.matcher_component_haplotypes.is_empty()
    {
        String::from("Variant")
    } else {
        String::from("Allele")
    };
    if report_gene.source_diplotypes.len() > 1 {
        title.push('s');
    }
    if report_gene.outside_call {
        title.push_str(" Reported");
    } else {
        title.push_str(" Matched");
    }
    title
}

fn html_amd_no_call(report_gene: &ReportGene) -> bool {
    match report_gene.call_source {
        ReportCallSource::None => true,
        ReportCallSource::Matcher => {
            report_gene.variant_reports.is_empty()
                || report_gene
                    .variant_reports
                    .iter()
                    .all(VariantReport::is_missing)
        }
        ReportCallSource::Outside => false,
    }
}

fn html_amd_gene_calls(report_gene: &ReportGene) -> Vec<String> {
    if !report_gene.outside_call {
        if report_gene.is_no_data_like_java() {
            return vec!["Not called - no variant data provided".to_owned()];
        }
        if !report_gene.is_reportable_like_java() {
            return vec!["Not called".to_owned()];
        }
    }

    let mut calls = BTreeSet::new();
    for diplotype in &report_gene.recommendation_diplotypes {
        let display_diplotypes = if diplotype.inferred_source_diplotypes.is_empty() {
            std::slice::from_ref(diplotype)
        } else {
            diplotype.inferred_source_diplotypes.as_slice()
        };
        for display_diplotype in display_diplotypes {
            calls.insert(display_diplotype.label.clone());
        }
    }
    calls.into_iter().collect()
}

fn html_amd_gene_call(report_gene: &ReportGene) -> String {
    html_amd_gene_calls(report_gene)
        .into_iter()
        .next()
        .unwrap_or_default()
}

fn html_amd_phase_status(report_gene: &ReportGene) -> String {
    if report_gene.outside_call {
        return "Unavailable for calls made outside PharmCAT".to_owned();
    }
    if report_gene.phased {
        if report_gene
            .variant_reports
            .iter()
            .any(|variant| variant.phase_set.is_some())
        {
            return "Phased, with phase sets (PS)".to_owned();
        }
        return "Phased".to_owned();
    }
    "Unphased".to_owned()
}

fn html_amd_show_unphased_note(report_gene: &ReportGene) -> bool {
    !report_gene.phased && !is_lowest_function_gene(gene_display_name(&report_gene.gene))
}

fn html_amd_total_missing_variants(report_gene: &ReportGene) -> usize {
    report_gene
        .variant_reports
        .iter()
        .filter(|variant| variant.is_missing())
        .count()
}

fn html_amd_uncalled_haps(report_gene: &ReportGene) -> String {
    let mut haplotypes = report_gene
        .uncalled_haplotypes
        .iter()
        .cloned()
        .collect::<Vec<_>>();
    haplotypes.sort_by(|left, right| compare_haplotype_names(left, right));
    haplotypes.join(", ")
}

fn html_amd_calls_at_positions(report_gene: &ReportGene) -> String {
    let mut html = String::new();
    let allele_functions = html_amd_allele_functions(report_gene);
    html.push_str("\n      <h4>Calls at Positions</h4>\n      <table>\n        <thead>\n        <tr>\n          <th>Position in VCF</th>\n          <th>RSID</th>\n          <th>Call in VCF</th>\n          <th>Reference</th>\n          <th>Related Alleles and Function</th>\n          <th>Warnings</th>\n        </tr>\n        </thead>\n        <tbody>\n");
    for variant in &report_gene.variant_reports {
        html.push_str("        <tr id=\"");
        html.push_str(&html_variant_position_id(variant));
        html.push_str("\" style=\"vertical-align: initial;\">\n          <td>");
        html.push_str(&html_variant_position(variant));
        html.push_str("</td>\n          <td");
        if let Some(db_snp_id) = &variant.db_snp_id {
            html.push_str(" id=\"");
            html.push_str(&html_escape(db_snp_id));
            html.push('"');
        }
        html.push('>');
        html.push_str(&html_escape(variant.db_snp_id.as_deref().unwrap_or("")));
        html.push_str("</td>\n");
        if variant.is_missing() {
            html.push_str(
                "          <td class=\"missingVariant\"><div class=\"callMessage\">Missing</div></td>\n",
            );
        } else {
            html.push_str("          ");
            html.push_str(&html_variant_alleles(variant));
            html.push('\n');
        }
        html.push_str("          <td>\n            ");
        html.push_str(&html_format_variant_call(
            variant.reference_allele.as_deref().unwrap_or(""),
        ));
        html.push_str("\n          </td>\n          <td>\n");
        if !variant.alleles.is_empty() {
            html.push_str("            <ul class=\"noBullet mt-0 mb-0\">\n");
            for allele in &variant.alleles {
                if let Some(allele_function) =
                    html_amd_allele_function(report_gene, variant, allele, &allele_functions)
                {
                    html.push_str("              ");
                    html.push_str(&allele_function);
                    html.push('\n');
                }
            }
            html.push_str("            </ul>\n");
        }
        html.push_str("          </td>\n          <td>\n");
        if !variant.warnings.is_empty() {
            html.push_str("            <ul class=\"warningList\">\n");
            for warning in &variant.warnings {
                html.push_str("              <li>");
                html.push_str(&html_escape(warning));
                html.push_str("</li>\n");
            }
            html.push_str("            </ul>\n");
        }
        html.push_str("          </td>\n        </tr>\n");
    }
    html.push_str("        </tbody>\n      </table>\n");
    html
}

fn html_amd_other_positions_of_interest(report_gene: &ReportGene) -> String {
    let mut html = String::new();
    html.push_str("\n      <h4>Other Positions of Interest</h4>\n");
    for message in report_gene
        .messages
        .iter()
        .filter(|message| message.is_extra_position_note())
    {
        html.push_str("      <div class=\"alert alert-warning\">");
        html.push_str(&html_message_message(message));
        html.push_str("</div>\n");
    }
    html.push_str("      <table>\n        <thead>\n        <tr>\n          <th>Position in VCF</th>\n          <th>RSID</th>\n          <th>Call in VCF</th>\n        </tr>\n        </thead>\n        <tbody>\n");
    for variant in &report_gene.variant_of_interest_reports {
        html.push_str("        <tr>\n          <td id=\"");
        html.push_str(&html_variant_position_id(variant));
        html.push_str("\">");
        html.push_str(&html_variant_position(variant));
        html.push_str("</td>\n          <td id=\"");
        html.push_str(&html_escape(variant.db_snp_id.as_deref().unwrap_or("")));
        html.push_str("\">");
        html.push_str(&html_escape(variant.db_snp_id.as_deref().unwrap_or("")));
        html.push_str("</td>\n");
        if variant.is_missing() {
            html.push_str("          <td class=\"missingVariant\"><em>missing</em></td>\n");
        } else {
            html.push_str("          ");
            html.push_str(&html_variant_alleles(variant));
            html.push('\n');
        }
        html.push_str("        </tr>\n");
    }
    html.push_str("        </tbody>\n      </table>\n");
    html
}

fn html_variant_position_id(variant: &VariantReport) -> String {
    format!(
        "{}_{}",
        html_escape(variant.chromosome.as_deref().unwrap_or("")),
        variant
            .position
            .map(|position| position.to_string())
            .unwrap_or_default()
    )
}

fn html_variant_position(variant: &VariantReport) -> String {
    format!(
        "{}:{}",
        html_escape(variant.chromosome.as_deref().unwrap_or("")),
        variant
            .position
            .map(|position| position.to_string())
            .unwrap_or_default()
    )
}

fn html_variant_alleles(variant: &VariantReport) -> String {
    let mut cell_class = if variant.is_non_reference() {
        "nonwild".to_owned()
    } else {
        String::new()
    };
    let mismatch = if variant.has_undocumented_variations {
        if cell_class.is_empty() {
            cell_class.push_str("mismatch");
        } else {
            cell_class.push_str(" mismatch");
        }
        "<div class=\"callMessage\">Undocumented variation</div>"
    } else {
        ""
    };
    let mut call = html_format_variant_call(variant.call.as_deref().unwrap_or(""));
    if let Some(phase_set) = variant.phase_set {
        call.push_str(" (PS:");
        call.push_str(&phase_set.to_string());
        call.push(')');
    }
    format!(
        "<td class=\"{}\">{}{}</td>",
        html_escape(&cell_class),
        call,
        mismatch
    )
}

fn html_format_variant_call(call: &str) -> String {
    if call.len() <= 9 {
        return html_escape(call);
    }
    call.split('/')
        .map(|allele| {
            if allele.len() <= 8 {
                return html_escape(allele);
            }
            allele
                .as_bytes()
                .chunks(9)
                .map(|chunk| html_escape(std::str::from_utf8(chunk).unwrap_or("")))
                .collect::<Vec<_>>()
                .join("<br />")
        })
        .collect::<Vec<_>>()
        .join("/<br />")
}

fn html_amd_allele_functions(report_gene: &ReportGene) -> BTreeMap<String, String> {
    if !report_gene.allele_function_map.is_empty() {
        return report_gene.allele_function_map.clone();
    }

    let mut functions = BTreeMap::new();
    for diplotype in report_gene
        .source_diplotypes
        .iter()
        .chain(report_gene.recommendation_diplotypes.iter())
        .flat_map(|diplotype| {
            std::iter::once(diplotype).chain(diplotype.inferred_source_diplotypes.iter())
        })
    {
        for haplotype in diplotype.allele1.iter().chain(diplotype.allele2.iter()) {
            functions
                .entry(haplotype.name.clone())
                .or_insert_with(|| haplotype.function.clone());
        }
    }
    functions
}

fn html_amd_allele_function(
    report_gene: &ReportGene,
    variant: &VariantReport,
    allele: &str,
    allele_functions: &BTreeMap<String, String>,
) -> Option<String> {
    if report_gene.gene == "CYP2C19"
        && allele.eq_ignore_ascii_case("*1")
        && !variant
            .db_snp_id
            .as_deref()
            .is_some_and(|db_snp_id| db_snp_id.eq_ignore_ascii_case("rs3758581"))
    {
        return None;
    }

    let function = allele_functions
        .get(allele)
        .filter(|function| !is_unspecified(function))
        .map(String::as_str)
        .unwrap_or("Unassigned");
    Some(format!(
        "<li>{} - {}</li>",
        html_escape(allele),
        html_escape(function)
    ))
}

fn html_message_class(message: &MessageAnnotation) -> String {
    let mut class = String::new();
    let mut last_dash = false;
    for character in message.name.chars() {
        if character.is_ascii_punctuation() || character.is_whitespace() {
            if !last_dash {
                class.push('-');
                last_dash = true;
            }
        } else {
            class.push(character);
            last_dash = false;
        }
    }
    class
}

fn html_message_message(message: &MessageAnnotation) -> String {
    let mut html = String::new();
    let mut remainder = message.message.as_str();
    while let Some(index) = remainder.find("PMID:") {
        let before = &remainder[..index];
        html.push_str(before);
        let after_prefix = &remainder[index + "PMID:".len()..];
        let digits_len = after_prefix
            .chars()
            .take_while(|character| character.is_ascii_digit())
            .map(char::len_utf8)
            .sum::<usize>();
        if digits_len == 0 {
            html.push_str("PMID:");
            remainder = after_prefix;
            continue;
        }
        let pmid = &after_prefix[..digits_len];
        html.push_str("<a href=\"https://pubmed.ncbi.nlm.nih.gov/");
        html.push_str(pmid);
        html.push_str("\" target=\"_blank\" rel=\"noopener noreferrer\">PMID:");
        html.push_str(pmid);
        html.push_str("</a>");
        remainder = &after_prefix[digits_len..];
    }
    html.push_str(remainder);
    html
}

fn html_recommendation_drugs(
    context: &ReportContext,
    options: &HtmlReportOptions,
) -> BTreeSet<String> {
    let mut drugs = BTreeSet::new();
    for source_reports in context.drug_reports.values() {
        for drug_report in source_reports.values() {
            if html_drug_report_has_recommendations(drug_report, options) {
                drugs.insert(drug_report.name.clone());
            }
        }
    }
    drugs
}

fn html_drugs_without_recommendations(
    context: &ReportContext,
    options: &HtmlReportOptions,
) -> BTreeSet<String> {
    let mut drugs = BTreeSet::new();
    if !options.compact {
        return drugs;
    }
    for source_reports in context.drug_reports.values() {
        for drug_report in source_reports.values() {
            if html_drug_report_visible_for_recommendations(drug_report, options)
                && !drug_report.is_matched()
            {
                drugs.insert(drug_report.name.clone());
            }
        }
    }
    drugs
}

fn html_drug_report_has_recommendations(
    drug_report: &DrugReport,
    options: &HtmlReportOptions,
) -> bool {
    html_drug_report_visible_for_recommendations(drug_report, options)
        && (!options.compact || drug_report.is_matched())
}

fn html_reports_for_drug<'a>(
    context: &'a ReportContext,
    drug: &str,
    options: &HtmlReportOptions,
) -> Vec<&'a DrugReport> {
    let mut reports = Vec::new();
    for source in PrescribingGuidanceSource::list_values() {
        if let Some(drug_report) = context
            .drug_reports
            .get(&source)
            .and_then(|source_reports| source_reports.get(drug))
            && html_drug_report_visible_for_recommendations(drug_report, options)
        {
            reports.push(drug_report);
        }
    }
    reports
}

fn html_matched_drug_report_rows(
    drug_report: &DrugReport,
    debug: bool,
    no_data_genes: &BTreeSet<String>,
    has_fda_recommendation: bool,
) -> String {
    let mut html = String::new();
    for guideline in &drug_report.guidelines {
        for annotation in &guideline.annotations {
            html.push_str("        <tr class=\"top-aligned ");
            html.push_str(&html_rx_annotation_class(
                drug_report.source,
                &drug_report.name,
            ));
            html.push_str("\"><td><p><b>");
            if let Some(url) = &guideline.url {
                html.push_str("<a href=\"");
                html.push_str(&html_escape(url));
                html.push_str("\" target=\"_blank\">");
                html.push_str(&html_source_value(drug_report.source, &guideline.name));
                html.push_str("</a>");
            } else {
                html.push_str(&html_source_value(drug_report.source, &guideline.name));
            }
            html.push_str("</b></p><p>");
            html.push_str(&html_population_value(
                drug_report.source,
                &annotation.population,
            ));
            html.push_str("</p>");
            if !html_rx_is_cpic_warfarin(&drug_report.name, drug_report.source) {
                html.push_str(&html_annotation_tags(annotation));
            }
            html.push_str("</td><td>");
            html.push_str("<div class=\"hint\">");
            html.push_str(&html_pluralize("Genotype", annotation.genotypes.len()));
            html.push_str("</div>");
            html.push_str(&html_annotation_genotypes(
                annotation,
                &drug_report.name,
                guideline,
                no_data_genes,
                debug,
            ));
            if !annotation.phenotypes.is_empty() {
                html.push_str("<div class=\"hint\">");
                html.push_str(&html_pluralize("Phenotype", annotation.phenotypes.len()));
                html.push_str("</div>");
                html.push_str(&html_print_rec_map(
                    &annotation.phenotypes,
                    Some("rx-phenotype"),
                ));
            }
            if !annotation.activity_scores.is_empty() {
                html.push_str("<div class=\"hint\">");
                html.push_str(&html_pluralize(
                    "Activity Score",
                    annotation.activity_scores.len(),
                ));
                html.push_str("</div>");
                html.push_str(&html_print_rec_map(
                    &annotation.activity_scores,
                    Some("rx-activity"),
                ));
            }
            html.push_str("</td>");
            if html_rx_is_cpic_warfarin(&drug_report.name, drug_report.source) {
                html.push_str(&html_cpic_warfarin_recommendation_cell(annotation));
            } else {
                html.push_str("<td>");
                html.push_str(&html_annotation_implications(annotation));
                html.push_str("</td><td class=\"drugRecommendation\">");
                if let Some(text) = &annotation.drug_recommendation {
                    html.push_str(text);
                }
                if has_fda_recommendation {
                    html.push_str("<a href=\"#rx-ast-");
                    html.push_str(&html_css_selector(&drug_report.name));
                    html.push_str("\" style=\"text-decoration: none\">&ast;</a>");
                }
                html.push_str("</td><td class=\"drugRecClass\">");
                html.push_str(&html_capitalize_na(&annotation.classification));
                html.push_str("</td>");
            }
            html.push_str("</tr>\n");
        }
    }
    html
}

fn html_rx_annotation_class(source: PrescribingGuidanceSource, drug: &str) -> String {
    format!("{}-{}", source.code_name(), html_java_css_selector(drug))
}

fn html_source_value(source: PrescribingGuidanceSource, guideline_name: &str) -> String {
    if source == PrescribingGuidanceSource::FdaAssoc {
        html_escape(
            guideline_name
                .split_once(':')
                .map_or(guideline_name, |(name, _)| name),
        )
    } else {
        html_escape(source.display_name())
    }
}

fn html_rx_is_cpic_warfarin(drug: &str, source: PrescribingGuidanceSource) -> bool {
    drug == "warfarin" && source == PrescribingGuidanceSource::CpicGuideline
}

fn html_annotation_tags(annotation: &AnnotationReport) -> String {
    let mut tags = Vec::new();
    if annotation.alternate_drug_available {
        tags.push("Alternate Drug");
    }
    if annotation.dosing_information {
        tags.push("Dosing Info");
    }
    if annotation.other_prescribing_guidance {
        tags.push("Other Guidance");
    }
    if tags.is_empty() {
        return "<div class=\"tag noAction\">No Action</div>".to_owned();
    }
    tags.into_iter()
        .map(|tag| format!("<div class=\"tag\">{tag}</div>"))
        .collect::<Vec<_>>()
        .join("\n")
}

fn html_cpic_warfarin_recommendation_cell(annotation: &AnnotationReport) -> String {
    let mut html = String::from("<td colspan=\"4\">");
    for message in annotation
        .messages
        .iter()
        .filter(|message| message.is_message())
    {
        html.push_str("<div class=\"alert alert-info ");
        html.push_str(&html_escape(&html_message_class(message)));
        html.push_str("\">");
        html.push_str(&message.message);
        html.push_str("</div>");
    }
    html.push_str("<div class=\"warfarinFlowchart\"><img src=\"https://files.cpicpgx.org/images/warfarin/warfarin_recommendation_diagram.png\" alt=\"Figure 2 from the CPIC guideline for warfarin\"/></div>");
    html.push_str("</td>");
    html
}

fn html_unmatched_drug_report_row(
    drug_report: &DrugReport,
    no_data_genes: &BTreeSet<String>,
) -> String {
    let mut html = String::new();
    let not_called = !drug_report
        .guidelines
        .iter()
        .any(GuidelineReport::is_reportable);
    html.push_str("        <tr class=\"top-aligned ");
    html.push_str(&html_rx_annotation_class(
        drug_report.source,
        &drug_report.name,
    ));
    html.push_str("\"><td class=\"top-aligned\"><p><b>");
    html.push_str(&html_escape(drug_report.source.display_name()));
    if !not_called && !drug_report.urls.is_empty() {
        html.push_str("</b><sup class=\"sources\">");
        html.push_str(&html_list_sources(&drug_report.urls));
        html.push_str("</sup>");
    } else {
        html.push_str("</b>");
    }
    html.push_str("</p></td>");
    if not_called {
        html.push_str("<td colspan=\"5\" class=\"top-aligned\"><span class=\"rx-no-call\">");
        let uncalled_genes = html_uncalled_genes(drug_report);
        if uncalled_genes.is_empty() {
            html.push_str("No call data provided");
        } else {
            html.push_str("No call data for ");
            html.push_str(&html_escape(&uncalled_genes));
        }
        html.push_str("</span>.</td>");
    } else {
        html.push_str("<td class=\"top-aligned\">");
        html.push_str("<div class=\"hint\">");
        html.push_str(&html_pluralize(
            "Genotype",
            html_unmatched_diplotypes(drug_report).len(),
        ));
        html.push_str("</div>");
        html.push_str(&html_rx_unmatched_diplotypes(drug_report, no_data_genes));
        if html_unmatched_has_non_dpyd_inferred(drug_report) {
            html.push_str("<sup><a href=\"#rx-dagger-");
            html.push_str(&html_css_selector(&drug_report.name));
            html.push_str("\" title=\"Inferred\">&dagger;</a></sup>");
        }
        if html_unmatched_has_dpyd_inferred(drug_report) {
            html.push_str("<sup><a href=\"#rx-ddagger-");
            html.push_str(&html_css_selector(&drug_report.name));
            html.push_str("\" title=\"Inferred\">&ddagger;</a></sup>");
        }
        html.push_str(
            "</td><td colspan=\"4\" class=\"top-aligned\"><div class=\"hint\">&nbsp;</div>",
        );
        html.push_str(&html_escape(drug_report.source.display_name()));
        html.push_str(
            " provides no genotype-based recommendations for this genotype, after evaluating the evidence.",
        );
        html.push_str("</td>");
    }
    html.push_str("</tr>\n");
    html
}

fn html_annotation_genotypes(
    annotation: &AnnotationReport,
    drug: &str,
    guideline: &GuidelineReport,
    no_data_genes: &BTreeSet<String>,
    debug: bool,
) -> String {
    if annotation.genotypes.is_empty() {
        return String::new();
    }
    let list_class = if annotation.genotypes.len() > 1 {
        "noPadding"
    } else {
        "noBullet"
    };
    let mut html = format!("<ul class=\"{list_class} mt-0\">");
    for genotype in &annotation.genotypes {
        html.push_str("<li>");
        html.push_str("<span class=\"noWrap\">");
        html.push_str(&html_rx_genotype(
            genotype,
            annotation,
            guideline,
            no_data_genes,
        ));
        if html_genotype_has_non_dpyd_inferred(genotype) {
            html.push_str("<sup><a href=\"#rx-dagger-");
            html.push_str(&html_css_selector(drug));
            html.push_str("\" title=\"Inferred\">&dagger;</a></sup>");
        }
        if html_genotype_has_dpyd_inferred(genotype) {
            html.push_str("<sup><a href=\"#rx-ddagger-");
            html.push_str(&html_css_selector(drug));
            html.push_str("\" title=\"Inferred\">&ddagger;</a></sup>");
        }
        html.push_str("</span>");
        if debug {
            html.push_str(&html_rx_genotype_debug(genotype));
        }
        html.push_str("</li>");
    }
    html.push_str("</ul>");
    html
}

fn html_rx_genotype(
    genotype: &RecommendationGenotype,
    annotation: &AnnotationReport,
    guideline: &GuidelineReport,
    no_data_genes: &BTreeSet<String>,
) -> String {
    let diplotypes = genotype
        .report_genes
        .iter()
        .flat_map(|report_gene| report_gene.recommendation_diplotypes.iter())
        .collect::<Vec<_>>();
    let mut html = if diplotypes.is_empty() {
        "Unknown genotype".to_owned()
    } else {
        html_render_rx_diplotypes(
            &diplotypes,
            "rx-dip",
            no_data_genes,
            &guideline.homozygous_component_haplotypes(),
        )
    };
    for variant in &annotation.highlighted_variants {
        if !html.is_empty() {
            html.push_str(";<br />");
        }
        html.push_str("<span class=\"rx-hl-var\">");
        html.push_str(&html_escape(variant));
        html.push_str("</span>");
    }
    html
}

fn html_rx_genotype_debug(genotype: &RecommendationGenotype) -> String {
    let diplotypes = genotype
        .report_genes
        .iter()
        .flat_map(|report_gene| report_gene.recommendation_diplotypes.iter())
        .filter(|diplotype| !diplotype.inferred_source_diplotypes.is_empty())
        .collect::<Vec<_>>();
    if diplotypes.is_empty() {
        return String::new();
    }
    let debug = html_render_rx_diplotypes_debug(&diplotypes);
    format!(
        "<div class=\"alert alert-debug\"><div class=\"hint\">Inferred:</div><span class=\"nowrap\">{debug}</span></div>"
    )
}

fn html_render_rx_diplotypes(
    diplotypes: &[&ReportDiplotype],
    dip_class: &str,
    no_data_genes: &BTreeSet<String>,
    homozygous_component_haplotypes: &BTreeSet<&str>,
) -> String {
    html_render_rx_diplotypes_with_options(
        diplotypes,
        30,
        false,
        dip_class,
        no_data_genes,
        homozygous_component_haplotypes,
    )
}

#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct HtmlRxDiplotypeDisplay {
    gene: String,
    label: String,
    allele1_missing: bool,
    outside_phenotype: bool,
    outside_activity_score: bool,
}

fn html_render_rx_diplotypes_with_options(
    diplotypes: &[&ReportDiplotype],
    length_limit: usize,
    for_debug: bool,
    dip_class: &str,
    no_data_genes: &BTreeSet<String>,
    homozygous_component_haplotypes: &BTreeSet<&str>,
) -> String {
    let mut display = BTreeSet::<HtmlRxDiplotypeDisplay>::new();
    for diplotype in diplotypes {
        let display_diplotypes = if for_debug || diplotype.inferred_source_diplotypes.is_empty() {
            vec![*diplotype]
        } else {
            diplotype
                .inferred_source_diplotypes
                .iter()
                .collect::<Vec<_>>()
        };
        for display_diplotype in display_diplotypes {
            display.insert(HtmlRxDiplotypeDisplay {
                gene: display_diplotype.gene.clone(),
                label: display_diplotype.label.clone(),
                allele1_missing: display_diplotype.allele1.is_none(),
                outside_phenotype: display_diplotype.outside_phenotype,
                outside_activity_score: display_diplotype.outside_activity_score,
            });
        }
    }

    let mut html = String::new();
    for diplotype in display {
        if !html.is_empty() {
            html.push_str(";<br />");
        }
        html.push_str("<span");
        if !for_debug {
            html.push_str(" class=\"");
            html.push_str(&html_escape(dip_class));
            html.push('"');
        }
        html.push('>');
        if no_data_genes.contains(&diplotype.gene) {
            html.push_str(&html_escape(&diplotype.gene));
            html.push(':');
            html.push_str("Uncalled - no variant data provided");
        } else {
            html.push_str("<a href=\"#");
            html.push_str(&html_escape(&diplotype.gene));
            html.push_str("\">");
            html.push_str(&html_escape(&diplotype.gene));
            html.push_str("</a>:");
            let call = html_rx_diplotype_call(&diplotype);
            html.push_str(&html_break_rx_call(&call, length_limit));
            if homozygous_component_haplotypes.contains(call.as_str()) {
                html.push_str(" (homozygous)");
            }
        }
        html.push_str("</span>");
    }
    html
}

fn html_render_rx_diplotypes_debug(diplotypes: &[&ReportDiplotype]) -> String {
    html_render_rx_diplotypes_with_options(
        diplotypes,
        30,
        true,
        "",
        &BTreeSet::new(),
        &BTreeSet::new(),
    )
}

fn html_rx_unmatched_diplotypes(
    drug_report: &DrugReport,
    no_data_genes: &BTreeSet<String>,
) -> String {
    let diplotypes = html_unmatched_diplotypes(drug_report);
    let homozygous_component_haplotypes = drug_report
        .guidelines
        .iter()
        .next()
        .map(GuidelineReport::homozygous_component_haplotypes)
        .unwrap_or_default();
    html_render_rx_diplotypes_with_options(
        &diplotypes,
        30,
        false,
        "rx-unmatched-dip",
        no_data_genes,
        &homozygous_component_haplotypes,
    )
}

fn html_rx_diplotype_call(diplotype: &HtmlRxDiplotypeDisplay) -> String {
    if diplotype.allele1_missing
        && (diplotype.outside_phenotype || diplotype.outside_activity_score)
    {
        "Not provided".to_owned()
    } else {
        diplotype.label.clone()
    }
}

fn html_list_sources(urls: &[String]) -> String {
    let mut html = String::new();
    for (index, url) in urls.iter().enumerate() {
        if index > 0 {
            html.push_str(", ");
        }
        html.push_str("<a href=\"");
        html.push_str(&html_escape(url));
        html.push_str("\" target=\"_blank\" rel=\"noopener noreferrer\">");
        html.push_str(&(index + 1).to_string());
        html.push_str("</a>");
    }
    html
}

fn html_unmatched_diplotypes(drug_report: &DrugReport) -> Vec<&ReportDiplotype> {
    drug_report
        .guidelines
        .iter()
        .flat_map(|guideline| guideline.report_genes.iter())
        .flat_map(|report_gene| report_gene.recommendation_diplotypes.iter())
        .collect()
}

fn html_unmatched_has_non_dpyd_inferred(drug_report: &DrugReport) -> bool {
    html_unmatched_inferred_genes(drug_report)
        .into_iter()
        .any(|gene| gene != "DPYD")
}

fn html_unmatched_has_dpyd_inferred(drug_report: &DrugReport) -> bool {
    html_unmatched_inferred_genes(drug_report)
        .into_iter()
        .any(|gene| gene == "DPYD")
}

fn html_unmatched_inferred_genes(drug_report: &DrugReport) -> BTreeSet<&str> {
    drug_report
        .guidelines
        .iter()
        .flat_map(|guideline| guideline.report_genes.iter())
        .flat_map(|report_gene| report_gene.recommendation_diplotypes.iter())
        .filter(|diplotype| diplotype.inferred)
        .map(|diplotype| diplotype.gene.as_str())
        .collect()
}

fn html_pluralize(word: &str, count: usize) -> String {
    if count == 1 {
        word.to_owned()
    } else {
        format!("{word}s")
    }
}

fn html_break_rx_call(call: &str, length_limit: usize) -> String {
    let escaped = html_escape(call);
    if length_limit > 0 && call.len() <= length_limit {
        return escaped;
    }
    if let Some(idx) = call.find('/') {
        let (first, second_with_delimiter) = call.split_at(idx + 1);
        let mut html = String::new();
        if first.len() > length_limit {
            html.push_str("<br />&nbsp;");
        }
        html.push_str(&html_escape(first));
        html.push_str("<br />&nbsp;");
        html.push_str(&html_escape(second_with_delimiter));
        return html;
    }
    let mut idx = call.rfind(" + ");
    while let Some(current) = idx {
        if current <= length_limit * 2 {
            break;
        }
        let next = call[..current].rfind(" + ");
        if next.is_none() {
            break;
        }
        idx = next;
    }
    if idx.is_none() {
        idx = call.rfind(" (");
    }
    let Some(idx) = idx else {
        return format!("<br />&nbsp;{escaped}");
    };
    format!(
        "<br />&nbsp;{}<br />&nbsp;{}",
        html_escape(&call[..idx]),
        html_escape(&call[idx..])
    )
}

fn html_annotation_implications(annotation: &AnnotationReport) -> String {
    match annotation.implications.as_slice() {
        [] => String::new(),
        [implication] => html_escape(implication),
        implications => {
            let mut html = String::from("<ul class=\"noPadding mt-0\">");
            for implication in implications {
                html.push_str("<li>");
                html.push_str(&html_capitalize_na(implication));
                html.push_str("</li>");
            }
            html.push_str("</ul>");
            html
        }
    }
}

fn html_reports_have_non_dpyd_inferred(reports: &[&DrugReport]) -> bool {
    reports.iter().any(|report| {
        report
            .guidelines
            .iter()
            .flat_map(|guideline| guideline.annotations.iter())
            .flat_map(|annotation| annotation.genotypes.iter())
            .any(html_genotype_has_non_dpyd_inferred)
    })
}

fn html_reports_have_dpyd_inferred(reports: &[&DrugReport]) -> bool {
    reports.iter().any(|report| {
        report
            .guidelines
            .iter()
            .flat_map(|guideline| guideline.annotations.iter())
            .flat_map(|annotation| annotation.genotypes.iter())
            .any(html_genotype_has_dpyd_inferred)
    })
}

fn html_reports_have_matched_fda(reports: &[&DrugReport]) -> bool {
    reports.iter().any(|report| {
        report.is_matched()
            && matches!(
                report.source,
                PrescribingGuidanceSource::FdaLabel | PrescribingGuidanceSource::FdaAssoc
            )
    })
}

fn html_report_footnotes<'a>(reports: &[&'a DrugReport]) -> BTreeSet<&'a MessageAnnotation> {
    reports
        .iter()
        .flat_map(|report| report.messages.iter())
        .filter(|message| message.is_footnote())
        .collect()
}

fn html_genotype_has_non_dpyd_inferred(genotype: &RecommendationGenotype) -> bool {
    html_genotype_inferred_genes(genotype)
        .into_iter()
        .any(|gene| gene != "DPYD")
}

fn html_genotype_has_dpyd_inferred(genotype: &RecommendationGenotype) -> bool {
    html_genotype_inferred_genes(genotype)
        .into_iter()
        .any(|gene| gene == "DPYD")
}

fn html_genotype_inferred_genes(genotype: &RecommendationGenotype) -> BTreeSet<&str> {
    genotype
        .report_genes
        .iter()
        .flat_map(|report_gene| report_gene.recommendation_diplotypes.iter())
        .filter(|diplotype| diplotype.inferred)
        .map(|diplotype| diplotype.gene.as_str())
        .collect()
}

fn html_uncalled_genes(drug_report: &DrugReport) -> String {
    drug_report
        .guidelines
        .iter()
        .flat_map(|guideline| {
            guideline
                .genes
                .iter()
                .map(|gene| gene_display_name(gene).to_owned())
        })
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect::<Vec<_>>()
        .join(", ")
}

fn html_print_rec_map(data: &BTreeMap<String, String>, css_class: Option<&str>) -> String {
    if data.len() == 1 {
        let class_attr = css_class
            .map(|class| format!(" class=\"{}\"", html_escape(class)))
            .unwrap_or_default();
        let value = data.values().next().expect("single map value");
        return format!("<p{}>{}</p>", class_attr, html_capitalize_na(value));
    }

    let mut html = String::from("<dl class=\"compact mt-0\">");
    for (key, value) in data {
        if let Some(class) = css_class {
            html.push_str("<div class=\"");
            html.push_str(&html_escape(class));
            html.push(' ');
            html.push_str(&html_escape(class));
            html.push_str("--");
            html.push_str(&html_escape(key));
            html.push_str("\">");
        }
        html.push_str("<dt>");
        html.push_str(&html_escape(key));
        html.push_str(":</dt><dd>");
        html.push_str(&html_capitalize_na(value));
        html.push_str("</dd>");
        if css_class.is_some() {
            html.push_str("</div>");
        }
    }
    html.push_str("</dl>");
    html
}

fn html_capitalize_na(text: &str) -> String {
    if text.trim().is_empty() {
        "Unspecified".to_owned()
    } else if text.eq_ignore_ascii_case(crate::phenotype::NA) {
        crate::phenotype::NA.to_uppercase()
    } else {
        html_escape(text)
    }
}

fn html_population_value(source: PrescribingGuidanceSource, population: &str) -> String {
    if source == PrescribingGuidanceSource::FdaAssoc {
        html_escape(population)
    } else {
        format!("Population:<br/>{}", html_capitalize_na(population))
    }
}

fn html_css_selector(text: &str) -> String {
    let mut selector = String::new();
    let mut last_dash = false;
    for character in text.chars().flat_map(char::to_lowercase) {
        if character.is_ascii_alphanumeric() || character == '_' || character == '-' {
            selector.push(character);
            last_dash = false;
        } else if !last_dash {
            selector.push('-');
            last_dash = true;
        }
    }
    selector.trim_matches('-').to_owned()
}

fn html_java_css_selector(text: &str) -> String {
    let mut selector = String::new();
    let mut last_underscore = false;
    for character in text.chars() {
        if character.is_ascii_alphanumeric() || character == '_' || character == '-' {
            selector.push(character);
            last_underscore = false;
        } else if !last_underscore {
            selector.push('_');
            last_underscore = true;
        }
    }
    selector.replace("_-_", "-").trim_matches('_').to_owned()
}

/// Writes a minimal HTML report to `path`.
pub fn write_report_html(
    context: &ReportContext,
    path: &Path,
    options: &HtmlReportOptions,
) -> Result<(), GuidanceLoadError> {
    if path.extension().and_then(|extension| extension.to_str()) != Some("html") {
        return Err(GuidanceLoadError::InvalidHtmlPath(path.to_path_buf()));
    }
    fs::write(path, report_html_string(context, options))?;
    Ok(())
}

fn html_escape(text: &str) -> String {
    text.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&#39;")
}

fn html_publication_citation(publication: &Publication) -> String {
    let title = publication.title.as_deref().unwrap_or("null");
    let url = publication.same_as.as_deref().unwrap_or("null");
    let mut citation = format!(
        "<a href=\"{}\" target=\"_blank\" rel=\"noopener noreferrer\">{}</a>",
        html_escape(url),
        html_escape(title)
    );
    if !title
        .chars()
        .last()
        .is_some_and(|character| character.is_ascii_punctuation())
    {
        citation.push('.');
    }
    if let Some(pmid) = publication.pmid.as_deref() {
        let journal = publication.journal.as_deref().unwrap_or("null");
        let year = publication
            .year
            .map(|year| year.to_string())
            .unwrap_or_else(|| "null".to_owned());
        citation.push(' ');
        citation.push_str("<i>");
        citation.push_str(&html_escape(journal));
        citation.push_str("</i>. ");
        citation.push_str(&year);
        citation.push_str(". ");
        citation.push_str("PMID:");
        citation.push_str(&html_escape(pmid));
    }
    citation
}

/// Writes Java-style pretty reporter JSON to `path`.
pub fn write_report_json(context: &ReportContext, path: &Path) -> Result<(), GuidanceLoadError> {
    if path.extension().and_then(|extension| extension.to_str()) != Some("json") {
        return Err(GuidanceLoadError::InvalidJsonPath(path.to_path_buf()));
    }
    let json = context.to_json_string()?;
    fs::write(path, json)?;
    Ok(())
}

/// Java `CallsOnlyFormat.NO_CALL_TAG`.
pub const CALLS_ONLY_NO_CALL_TAG: &str = "no call";
/// Java `CallsOnlyFormat.HEADER_SAMPLE_ID`.
pub const CALLS_ONLY_HEADER_SAMPLE_ID: &str = "Sample ID";
/// Java `CallsOnlyFormat.HEADER_VARIANTS`.
pub const CALLS_ONLY_HEADER_VARIANTS: &str = "Variants";
/// Java `CallsOnlyFormat.HEADER_UNDOCUMENTED_VARIANTS`.
pub const CALLS_ONLY_HEADER_UNDOCUMENTED_VARIANTS: &str = "Undocumented variants";

/// Options for the initial calls-only TSV writer.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CallsOnlyTsvOptions {
    /// PharmCAT version label used in the first line.
    pub pharmcat_version: String,
    /// Optional sample id to include as the first data column.
    pub sample_id: Option<String>,
    /// Whether to include the sample id column.
    pub show_sample_id: bool,
    /// Whether to include variants column.
    pub show_variants: bool,
    /// Whether to list missing variants instead of yes/no.
    pub show_missing_variants: bool,
    /// Whether to include undocumented variants column.
    pub show_undocumented_variants: bool,
    /// Whether undocumented variants are treated as reference.
    pub treat_undocumented_variations_as_reference: bool,
    /// Whether writes append to an existing single calls-only file.
    pub single_file_mode: bool,
    /// Sample metadata columns appended after recommendation lookup fields.
    pub sample_properties: BTreeMap<String, String>,
}

impl Default for CallsOnlyTsvOptions {
    fn default() -> Self {
        Self {
            pharmcat_version: "unknown".to_owned(),
            sample_id: None,
            show_sample_id: false,
            show_variants: false,
            show_missing_variants: false,
            show_undocumented_variants: false,
            treat_undocumented_variations_as_reference: false,
            single_file_mode: false,
            sample_properties: BTreeMap::new(),
        }
    }
}

/// Generates the first Rust calls-only TSV surface for the currently ported report model.
pub fn calls_only_tsv_string(context: &ReportContext, options: &CallsOnlyTsvOptions) -> String {
    let mut tsv = String::new();
    tsv.push_str("PharmCAT ");
    tsv.push_str(&options.pharmcat_version);
    tsv.push('\n');
    if options.show_sample_id {
        tsv.push_str(CALLS_ONLY_HEADER_SAMPLE_ID);
        tsv.push('\t');
    }
    let mut header = CALLS_ONLY_BASE_HEADER.to_owned();
    if options.show_variants {
        header = header.replacen(
            "\tMissing positions",
            &format!("\t{CALLS_ONLY_HEADER_VARIANTS}\tMissing positions"),
            1,
        );
    }
    if options.show_undocumented_variants {
        header = header.replacen(
            "\tRecommendation Lookup Diplotype",
            &format!(
                "\t{CALLS_ONLY_HEADER_UNDOCUMENTED_VARIANTS}\tRecommendation Lookup Diplotype"
            ),
            1,
        );
    }
    for key in options.sample_properties.keys() {
        header.push('\t');
        header.push_str(key);
    }
    tsv.push_str(&header);
    tsv.push('\n');

    tsv.push_str(&calls_only_tsv_rows(context, options));
    tsv
}

fn calls_only_tsv_rows(context: &ReportContext, options: &CallsOnlyTsvOptions) -> String {
    let mut tsv = String::new();
    for report_gene in context.gene_reports.values() {
        append_calls_only_row(&mut tsv, report_gene, options);
    }
    tsv
}

/// Writes calls-only TSV to `path`.
pub fn write_calls_only_tsv(
    context: &ReportContext,
    path: &Path,
    options: &CallsOnlyTsvOptions,
) -> Result<(), GuidanceLoadError> {
    if path.extension().and_then(|extension| extension.to_str()) != Some("tsv") {
        return Err(GuidanceLoadError::InvalidTsvPath(path.to_path_buf()));
    }
    if options.single_file_mode && path.exists() {
        use std::io::Write;
        let mut file = fs::OpenOptions::new()
            .append(true)
            .create(true)
            .open(path)?;
        file.write_all(calls_only_tsv_rows(context, options).as_bytes())?;
    } else {
        fs::write(path, calls_only_tsv_string(context, options))?;
    }
    Ok(())
}

const CALLS_ONLY_BASE_HEADER: &str = "Gene\tSource Diplotype\tPhenotype\tActivity Score\tHaplotype 1\tHaplotype 1 Function\tHaplotype 1 Activity Value\tHaplotype 2\tHaplotype 2 Function\tHaplotype 2 Activity Value\tOutside Call\tMatch Score\tMissing positions\tRecommendation Lookup Diplotype\tRecommendation Lookup Phenotype\tRecommendation Lookup Activity Score";

fn append_calls_only_row(
    tsv: &mut String,
    report_gene: &ReportGene,
    options: &CallsOnlyTsvOptions,
) {
    if options.show_sample_id {
        if let Some(sample_id) = &options.sample_id {
            tsv.push_str(sample_id);
        }
        tsv.push('\t');
    }

    let source_diplotype = report_gene
        .source_diplotype
        .as_deref()
        .map(str::to_owned)
        .unwrap_or_else(|| report_gene.lookup_keys.join(" OR "));
    let phenotype = calls_only_value_list(&report_gene.phenotypes);
    let activity_score = if report_gene.is_activity_score_type {
        report_gene
            .activity_score
            .as_deref()
            .map(calls_only_value)
            .unwrap_or_default()
    } else {
        String::new()
    };
    let recommendation_phenotype = calls_only_value_list(&report_gene.lookup_keys);
    let recommendation_activity_score = if report_gene.is_activity_score_type {
        activity_score.clone()
    } else {
        String::new()
    };

    let mut fields = vec![
        report_gene.gene.clone(),
        source_diplotype,
        phenotype,
        activity_score,
        String::new(),
        String::new(),
        String::new(),
        String::new(),
        String::new(),
        String::new(),
        if report_gene.outside_call {
            "yes"
        } else {
            "no"
        }
        .to_owned(),
        report_gene.match_score.clone().unwrap_or_default(),
    ];
    if options.show_variants {
        fields.push(calls_only_variants(report_gene));
    }
    fields.push(calls_only_missing_variants(report_gene, options));
    if options.show_undocumented_variants {
        fields.push(calls_only_undocumented_variants(report_gene, options));
    }
    fields.extend([
        String::new(),
        recommendation_phenotype,
        recommendation_activity_score,
    ]);
    fields.extend(options.sample_properties.values().cloned());
    tsv.push_str(&fields.join("\t"));
    tsv.push('\n');
}

fn calls_only_variants(report_gene: &ReportGene) -> String {
    report_gene
        .variant_reports
        .iter()
        .filter(|variant| !variant.is_missing() && variant.is_non_reference())
        .filter_map(|variant| {
            Some(format!(
                "{}:{}",
                variant.position?,
                variant.call.as_deref().unwrap_or("")
            ))
        })
        .collect::<Vec<_>>()
        .join(", ")
}

fn calls_only_missing_variants(report_gene: &ReportGene, options: &CallsOnlyTsvOptions) -> String {
    if options.show_missing_variants {
        return report_gene
            .variant_reports
            .iter()
            .filter(|variant| variant.is_missing())
            .filter_map(|variant| variant.position.map(|position| position.to_string()))
            .collect::<Vec<_>>()
            .join(", ");
    }

    if !report_gene.variant_reports.is_empty()
        && report_gene
            .variant_reports
            .iter()
            .any(VariantReport::is_missing)
    {
        "yes".to_owned()
    } else {
        "no".to_owned()
    }
}

fn calls_only_undocumented_variants(
    report_gene: &ReportGene,
    options: &CallsOnlyTsvOptions,
) -> String {
    let mut value = report_gene
        .variant_reports
        .iter()
        .filter(|variant| variant.has_undocumented_variations)
        .filter_map(|variant| variant.position.map(|position| position.to_string()))
        .collect::<Vec<_>>()
        .join(", ");
    if !value.is_empty() && options.treat_undocumented_variations_as_reference {
        value.push_str(" treat as reference");
    }
    value
}

fn calls_only_value_list(values: &[String]) -> String {
    values
        .iter()
        .map(|value| calls_only_value(value))
        .collect::<Vec<_>>()
        .join(", ")
}

fn calls_only_value(value: &str) -> String {
    if value.trim().is_empty()
        || value == crate::phenotype::NA
        || value == crate::phenotype::NO_RESULT
    {
        " ".to_owned()
    } else {
        value.to_owned()
    }
}

fn report_diplotype_is_unknown(diplotype: &ReportDiplotype) -> bool {
    fn unknown_allele(haplotype: Option<&ReportHaplotype>) -> bool {
        haplotype.is_none_or(|haplotype| haplotype.name == "Unknown")
    }

    unknown_allele(diplotype.allele1.as_ref())
        && unknown_allele(diplotype.allele2.as_ref())
        && (diplotype.phenotypes.is_empty()
            || diplotype
                .phenotypes
                .iter()
                .any(|phenotype| phenotype == crate::phenotype::NO_RESULT))
}

fn split_slco1b1_rs4149056_variant_call(variant: &VariantReport) -> Option<[String; 2]> {
    let call = variant.call.as_deref()?.trim();
    if call.is_empty() {
        return None;
    }
    let mut alleles = call
        .split(['|', '/'])
        .map(str::to_owned)
        .collect::<Vec<_>>();
    if alleles.len() != 2 {
        return None;
    }
    alleles.sort_by(|left, right| right.cmp(left));
    Some([
        alleles
            .first()
            .expect("checked two rs4149056 alleles")
            .clone(),
        alleles
            .get(1)
            .expect("checked two rs4149056 alleles")
            .clone(),
    ])
}

fn slco1b1_rs4149056_allele_to_haplotype(allele: &str) -> Option<&'static str> {
    match allele {
        "T" => Some("*1"),
        "C" => Some("*5"),
        _ => None,
    }
}

fn variant_report_alleles(definition: &DefinitionFile, locus: &VariantLocus) -> Vec<String> {
    let Some(index) = definition.index_for_position(locus.position) else {
        return Vec::new();
    };

    let mut alleles = definition
        .named_alleles
        .iter()
        .filter(|named_allele| !named_allele.reference)
        .filter(|named_allele| named_allele.alleles.get(index).is_some_and(Option::is_some))
        .map(|named_allele| {
            definition
                .suballeles_map
                .get(&named_allele.name)
                .cloned()
                .unwrap_or_else(|| named_allele.name.clone())
        })
        .collect::<Vec<_>>();
    alleles.sort_by(|left, right| compare_haplotype_names(left, right));
    alleles.dedup();
    alleles
}

fn is_valid_variant_call(call: &str) -> bool {
    let mut alleles = call.split(['|', '/']);
    let Some(allele1) = alleles.next() else {
        return false;
    };
    if allele1.trim().is_empty()
        || !allele1
            .chars()
            .all(|character| character.is_ascii_alphabetic())
    {
        return false;
    }
    let Some(allele2) = alleles.next() else {
        return true;
    };
    alleles.next().is_none()
        && !allele2.trim().is_empty()
        && allele2
            .chars()
            .all(|character| character.is_ascii_alphabetic())
}

/// Prescribing-guidance loading error.
#[derive(Debug)]
pub enum GuidanceLoadError {
    /// I/O error.
    Io(io::Error),
    /// JSON parse error.
    Json(serde_json::Error),
    /// Output path does not end in `.json`.
    InvalidJsonPath(PathBuf),
    /// Output path does not end in `.tsv`.
    InvalidTsvPath(PathBuf),
    /// Output path does not end in `.html`.
    InvalidHtmlPath(PathBuf),
}

impl std::fmt::Display for GuidanceLoadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io(error) => error.fmt(f),
            Self::Json(error) => error.fmt(f),
            Self::InvalidJsonPath(path) => {
                write!(
                    f,
                    "Invalid format: {} does not end with .json",
                    path.display()
                )
            }
            Self::InvalidTsvPath(path) => {
                write!(
                    f,
                    "Invalid format: {} does not end with .tsv",
                    path.display()
                )
            }
            Self::InvalidHtmlPath(path) => {
                write!(
                    f,
                    "Invalid format: {} does not end with .html",
                    path.display()
                )
            }
        }
    }
}

impl std::error::Error for GuidanceLoadError {}

impl From<io::Error> for GuidanceLoadError {
    fn from(error: io::Error) -> Self {
        Self::Io(error)
    }
}

impl From<serde_json::Error> for GuidanceLoadError {
    fn from(error: serde_json::Error) -> Self {
        Self::Json(error)
    }
}

#[cfg(test)]
mod tests {
    use std::{
        collections::{BTreeMap, BTreeSet},
        fs,
        path::Path,
    };

    use serde_json::Value;

    use crate::{
        definition::{DefinitionFile, DefinitionReader, VariantLocus, read_definition_file},
        matcher::{
            GeneCallKind, GeneCallResult, GeneCallWarning, MatchData,
            call_dpyd_lowest_function_gene, call_ryr1_lowest_function_gene, call_standard_gene,
        },
        phenotype::{
            DiplotypeAnnotationInput, OutsideCallValidation, PhenotypeMap, parse_outside_call_line,
            parse_outside_calls_str, read_gene_phenotype_file,
        },
        vcf::{SampleAlleleSummary, read_record_summaries},
    };

    use super::{
        AnnotationReport, CALLS_ONLY_HEADER_UNDOCUMENTED_VARIANTS, CALLS_ONLY_HEADER_VARIANTS,
        CallsOnlyTsvOptions, DataSource, DrugReport, GuidelineReport, HtmlReportOptions,
        HtmlTemplateSet, MatchLogic, MessageAnnotation, MessageCatalog, PgkbGuidelineCollection,
        PrescribingGuidanceSource, Publication, RecommendationAnnotation, RecommendationGenotype,
        ReportCallSource, ReportContext, ReportContextFromMatcherError, ReportDiplotype,
        ReportGene, Slco1b1CustomCallError, VariantReport, calls_only_tsv_string,
        data_source_from_definition, load_disclaimers_template, make_recommendation_genotypes,
        map_contains, report_html_string, sort_report_diplotypes, write_calls_only_tsv,
        write_report_html, write_report_json,
    };

    const GUIDANCE_PATH: &str = "../../repos/PharmCAT/src/main/resources/org/pharmgkb/pharmcat/reporter/prescribing_guidance.json";
    const PHENOTYPE_PATH: &str =
        "../../repos/PharmCAT/src/main/resources/org/pharmgkb/pharmcat/phenotype";
    const DPYD_DEFINITION_PATH: &str = "../../repos/PharmCAT/src/main/resources/org/pharmgkb/pharmcat/definition/alleles/DPYD_translation.json";
    const DPYD_PHENOTYPE_PATH: &str =
        "../../repos/PharmCAT/src/main/resources/org/pharmgkb/pharmcat/phenotype/DPYD.json";
    const CYP2D6_PHENOTYPE_PATH: &str =
        "../../repos/PharmCAT/src/main/resources/org/pharmgkb/pharmcat/phenotype/CYP2D6.json";
    const CYP3A5_DEFINITION_PATH: &str =
        "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/haplotyper.json";
    const CYP3A5_VCF_PATH: &str =
        "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/haplotyper.vcf";
    const CYP3A5_PHENOTYPE_PATH: &str =
        "../../repos/PharmCAT/src/main/resources/org/pharmgkb/pharmcat/phenotype/CYP3A5.json";
    const UGT1A1_COMBINATION_DEFINITION_PATH: &str = "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-combination.json";
    const UGT1A1_PARTIAL_WITH_COMBINATION_VCF_PATH: &str = "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-partialWithCombination.vcf";
    const UGT1A1_PHENOTYPE_PATH: &str =
        "../../repos/PharmCAT/src/main/resources/org/pharmgkb/pharmcat/phenotype/UGT1A1.json";
    const SLCO1B1_DEFINITION_PATH: &str = "../../repos/PharmCAT/src/main/resources/org/pharmgkb/pharmcat/definition/alleles/SLCO1B1_translation.json";
    const SLCO1B1_PHENOTYPE_PATH: &str =
        "../../repos/PharmCAT/src/main/resources/org/pharmgkb/pharmcat/phenotype/SLCO1B1.json";
    const RYR1_DEFINITION_PATH: &str = "../../repos/PharmCAT/src/main/resources/org/pharmgkb/pharmcat/definition/alleles/RYR1_translation.json";
    const RYR1_PHENOTYPE_PATH: &str =
        "../../repos/PharmCAT/src/main/resources/org/pharmgkb/pharmcat/phenotype/RYR1.json";
    const MESSAGES_PATH: &str =
        "../../repos/PharmCAT/src/main/resources/org/pharmgkb/pharmcat/reporter/messages.json";
    const DISCLAIMERS_PATH: &str =
        "../../repos/PharmCAT/src/main/resources/org/pharmgkb/pharmcat/reporter/disclaimers.hbs";
    const REPORTER_RESOURCE_DIR: &str =
        "../../repos/PharmCAT/src/main/resources/org/pharmgkb/pharmcat/reporter";

    #[test]
    fn loads_guideline_collection_like_java_pgkb_guideline_collection() {
        let collection = PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH))
            .expect("prescribing guidance");
        let guideline_count = collection.guideline_packages().len();

        assert!(guideline_count > 50);

        let cpic = collection
            .guidelines_from_source(PrescribingGuidanceSource::CpicGuideline)
            .len();
        let dpwg = collection
            .guidelines_from_source(PrescribingGuidanceSource::DpwgGuideline)
            .len();
        let fda_labels = collection
            .guidelines_from_source(PrescribingGuidanceSource::FdaLabel)
            .len();
        let fda_assocs = collection
            .guidelines_from_source(PrescribingGuidanceSource::FdaAssoc)
            .len();

        assert_eq!(cpic + dpwg + fda_labels + fda_assocs, guideline_count);
        assert!(
            !collection
                .genes_used_in_source(DataSource::Dpwg)
                .contains("CACNA1S")
        );
    }

    #[test]
    fn indexes_guidelines_by_related_chemical_and_source_like_java() {
        let collection = PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH))
            .expect("prescribing guidance");

        let abacavir_cpic = collection
            .find_guideline_packages("abacavir", PrescribingGuidanceSource::CpicGuideline);
        assert_eq!(abacavir_cpic.len(), 1);
        assert_eq!(
            abacavir_cpic[0].guideline.name,
            "Annotation of CPIC Guideline for abacavir and HLA-B"
        );
        assert_eq!(
            abacavir_cpic[0].genes(),
            ["HLA-B".to_owned()].into_iter().collect()
        );
        assert_eq!(
            abacavir_cpic[0].drugs(),
            ["abacavir".to_owned()].into_iter().collect()
        );
        assert!(collection.genes_with_recommendations().contains("HLA-B"));
    }

    #[test]
    fn loads_message_catalog_like_java_message_helper() {
        let catalog = MessageCatalog::from_path(Path::new(MESSAGES_PATH)).expect("messages");

        assert_eq!(catalog.messages().len(), 64);
        assert_eq!(catalog.messages_for_gene("CYP2B6").len(), 6);

        let combo_naming = catalog
            .message("pcat-combo-naming")
            .expect("pcat combo naming static message");
        assert_eq!(combo_naming.exception_type, "note");
        assert!(combo_naming.message.contains("combination allele calls"));

        assert_eq!(
            catalog.static_message_keys(),
            [
                "pcat-call-multimatch".to_owned(),
                "pcat-combo-naming".to_owned(),
                "pcat-combo-unphased".to_owned(),
                "pcat-cpic-warfarin-1-flowchart".to_owned(),
                "pcat-cpic-warfarin-2-vkorc1".to_owned(),
                "pcat-cyp2d6-gene-note".to_owned(),
                "pcat-cyp2d6-research-mode".to_owned(),
                "pcat-dpyd-hapB3-exonic-only".to_owned(),
                "pcat-dpyd-hapB3-intronic-mismatch-exonic".to_owned(),
                "pcat-outside-call".to_owned(),
                "pcat-score-multimatch".to_owned(),
            ]
            .into_iter()
            .collect()
        );
    }

    #[test]
    fn message_catalog_filters_drug_messages_by_source_like_java() {
        let catalog = MessageCatalog::from_path(Path::new(MESSAGES_PATH)).expect("messages");

        let cpic_warfarin =
            catalog.messages_for_drug("warfarin", PrescribingGuidanceSource::CpicGuideline);
        assert_eq!(cpic_warfarin.len(), 4);
        assert!(
            cpic_warfarin
                .iter()
                .any(|message| message.name == "pcat-cpic-warfarin-1-flowchart")
        );

        let dpwg_warfarin =
            catalog.messages_for_drug("warfarin", PrescribingGuidanceSource::DpwgGuideline);
        assert_eq!(
            dpwg_warfarin
                .iter()
                .map(|message| message.name.as_str())
                .collect::<Vec<_>>(),
            vec!["CYP4F2 *1/*4 warning"]
        );
    }

    #[test]
    fn drug_messages_apply_blank_gene_and_matching_gene_rules_like_java_message_helper() {
        let collection = PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH))
            .expect("prescribing guidance");
        let catalog = MessageCatalog::from_messages(vec![
            gene_message(
                "blank-gene-drug-note",
                MessageAnnotation::TYPE_NOTE,
                MatchLogic {
                    drugs: vec!["warfarin".to_owned()],
                    ..MatchLogic::default()
                },
            ),
            gene_message(
                "gene-backed-drug-note",
                MessageAnnotation::TYPE_NOTE,
                MatchLogic {
                    gene: Some("CYP2C9".to_owned()),
                    drugs: vec!["warfarin".to_owned()],
                    ..MatchLogic::default()
                },
            ),
            gene_message(
                "gene-message-not-present",
                MessageAnnotation::TYPE_NOTE,
                MatchLogic {
                    gene: Some("CYP2C9".to_owned()),
                    drugs: vec!["warfarin".to_owned()],
                    ..MatchLogic::default()
                },
            ),
        ]);
        let mut context = ReportContext::from_gene_reports(
            &collection,
            [ReportGene::new("CYP2C9", ["*1/*1".to_owned()])
                .with_variant_reports([VariantReport::new("rs1057910", Some("A/G"))])
                .with_messages([gene_message(
                    "gene-backed-drug-note",
                    MessageAnnotation::TYPE_NOTE,
                    MatchLogic {
                        gene: Some("CYP2C9".to_owned()),
                        ..MatchLogic::default()
                    },
                )])],
            None,
        );

        context.apply_matching_drug_messages(&catalog);

        let warfarin = context
            .drug_report(PrescribingGuidanceSource::CpicGuideline, "warfarin")
            .expect("warfarin CPIC report");
        let message_names = warfarin
            .messages
            .iter()
            .map(|message| message.name.as_str())
            .collect::<BTreeSet<_>>();
        assert_eq!(
            message_names,
            BTreeSet::from(["blank-gene-drug-note", "gene-backed-drug-note"])
        );
    }

    #[test]
    fn drug_messages_add_missing_variants_note_like_java_report_context() {
        let collection = PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH))
            .expect("prescribing guidance");
        let catalog = MessageCatalog::from_messages(Vec::new());
        let mut context = ReportContext::from_gene_reports(
            &collection,
            [
                ReportGene::new("CYP2C9", ["*1/*1".to_owned()]).with_variant_reports([
                    VariantReport::new("rsMissing", None::<String>),
                    VariantReport::new("rsPresent", Some("A/G")),
                ]),
            ],
            None,
        );

        context.apply_matching_drug_messages(&catalog);

        let warfarin = context
            .drug_report(PrescribingGuidanceSource::CpicGuideline, "warfarin")
            .expect("warfarin CPIC report");
        let message = warfarin
            .messages
            .iter()
            .find(|message| message.name == "missing-variants")
            .expect("missing variants drug message");
        assert_eq!(message.exception_type, MessageAnnotation::TYPE_NOTE);
        assert!(
            message
                .message
                .contains("Some position data used to define CYP2C9 alleles is missing")
        );
    }

    #[test]
    fn report_as_genotype_messages_add_highlighted_variants_like_java_message_helper() {
        let collection = PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH))
            .expect("prescribing guidance");
        let catalog = MessageCatalog::from_path(Path::new(MESSAGES_PATH)).expect("messages");
        let mut context = ReportContext::from_gene_reports(
            &collection,
            [ReportGene::new("CYP2C9", ["*1/*1".to_owned()])
                .with_variant_reports([VariantReport::new("rs12777823", Some("C|T"))])],
            None,
        );

        context.apply_report_as_genotype_messages(&catalog);

        let warfarin = context
            .drug_report(PrescribingGuidanceSource::CpicGuideline, "warfarin")
            .expect("warfarin CPIC report");
        let highlighted = warfarin
            .guidelines
            .iter()
            .flat_map(|guideline| &guideline.annotations)
            .flat_map(|annotation| &annotation.highlighted_variants)
            .cloned()
            .collect::<BTreeSet<_>>();
        assert_eq!(highlighted, BTreeSet::from(["rs12777823:C/T".to_owned()]));
    }

    #[test]
    fn report_as_genotype_uses_variants_of_interest_like_java_message_helper() {
        let collection = PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH))
            .expect("prescribing guidance");
        let catalog = MessageCatalog::from_path(Path::new(MESSAGES_PATH)).expect("messages");
        let mut context = ReportContext::from_gene_reports(
            &collection,
            [ReportGene::new("CYP2C9", ["*1/*1".to_owned()])
                .with_variant_of_interest_reports([VariantReport::new("rs12777823", Some("C|T"))])],
            None,
        );

        context.apply_report_as_genotype_messages(&catalog);

        let warfarin = context
            .drug_report(PrescribingGuidanceSource::CpicGuideline, "warfarin")
            .expect("warfarin CPIC report");
        let highlighted = warfarin
            .guidelines
            .iter()
            .flat_map(|guideline| &guideline.annotations)
            .flat_map(|annotation| &annotation.highlighted_variants)
            .cloned()
            .collect::<BTreeSet<_>>();
        assert_eq!(highlighted, BTreeSet::from(["rs12777823:C/T".to_owned()]));
    }

    #[test]
    fn loads_disclaimer_template_resource() {
        let disclaimers =
            load_disclaimers_template(Path::new(DISCLAIMERS_PATH)).expect("disclaimers");

        assert!(disclaimers.starts_with("<section id=\"disclaimer\">"));
        assert!(disclaimers.contains("Section IV: Disclaimers and Other Information"));
        assert!(disclaimers.contains("CPIC Guideline Disclaimers and Caveats"));
    }

    #[test]
    fn recommendation_lookup_maps_match_java_recommendation_utils() {
        let super_set = [
            ("GENEX".to_owned(), Value::String("TEST".to_owned())),
            ("GENEY".to_owned(), Value::String("TOAST".to_owned())),
        ]
        .into_iter()
        .collect::<BTreeMap<_, _>>();
        let sub_set = [("GENEX".to_owned(), Value::String("TEST".to_owned()))]
            .into_iter()
            .collect::<BTreeMap<_, _>>();

        assert!(map_contains(&super_set, &sub_set));
        assert!(!map_contains(&sub_set, &super_set));
        assert!(!map_contains(&BTreeMap::new(), &sub_set));
    }

    #[test]
    fn recommendation_matches_genotype_like_java_recommendation_annotation_test() {
        let recommendation = RecommendationAnnotation {
            id: "PA-test".to_owned(),
            name: "Recommendation PA-test".to_owned(),
            population: None,
            classification: None,
            related_chemicals: Vec::new(),
            text: None,
            implications: Vec::new(),
            lookup_key: vec![
                [("GENEX".to_owned(), Value::String("TEST".to_owned()))]
                    .into_iter()
                    .collect(),
            ],
            dosing_information: false,
            alternate_drug_available: false,
            other_prescribing_guidance: false,
        };
        let genotype = RecommendationGenotype::from_gene_lookup_keys([
            ("GENEX".to_owned(), vec!["TEST".to_owned()]),
            ("GENEY".to_owned(), vec!["TOAST".to_owned()]),
        ]);

        assert!(recommendation.matches_genotype(&genotype));
    }

    #[test]
    fn report_context_builds_matched_drug_report_from_real_guidance() {
        let collection = PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH))
            .expect("prescribing guidance");
        let context = ReportContext::from_gene_reports(
            &collection,
            [ReportGene::new("HLA-B", ["*57:01 positive".to_owned()])],
            Some("test report".to_owned()),
        );

        let abacavir = context
            .drug_report(PrescribingGuidanceSource::CpicGuideline, "abacavir")
            .expect("abacavir CPIC report");

        assert_eq!(abacavir.name, "abacavir");
        assert_eq!(abacavir.id, "PA448004");
        assert_eq!(abacavir.source, PrescribingGuidanceSource::CpicGuideline);
        assert!(abacavir.is_matched());
        assert_eq!(abacavir.matched_annotation_count(), 1);

        let guideline = abacavir.guidelines.iter().next().expect("guideline");
        assert_eq!(
            guideline.name,
            "Annotation of CPIC Guideline for abacavir and HLA-B"
        );
        assert_eq!(guideline.genes, ["HLA-B".to_owned()].into_iter().collect());
        assert!(guideline.is_matched());

        let annotation = guideline.annotations.iter().next().expect("annotation");
        assert_eq!(annotation.local_id, "CPIC-PA166296759");
        assert_eq!(annotation.phenotypes.get("HLA-B"), None);
        assert_eq!(
            annotation.drug_recommendation.as_deref(),
            Some("Abacavir is not recommended")
        );
    }

    #[test]
    fn report_context_serializes_java_style_report_json_fields() {
        let collection = PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH))
            .expect("prescribing guidance");
        let context = ReportContext::from_gene_reports(
            &collection,
            [ReportGene::new("HLA-B", ["*58:01 positive".to_owned()])],
            Some("test report".to_owned()),
        );

        let json = context.to_json_string().expect("report JSON");
        assert!(json.starts_with("{\n  \"title\""));

        let value = serde_json::from_str::<Value>(&json).expect("report JSON value");
        assert_eq!(value["title"], "test report");
        assert_eq!(
            value["dataVersion"],
            collection.version().unwrap_or(crate::phenotype::NA)
        );
        assert!(value["messages"].as_array().is_some_and(Vec::is_empty));
        assert!(
            value["unannotatedGeneCalls"]
                .as_array()
                .is_some_and(Vec::is_empty)
        );

        let hla_b = &value["genes"]["HLA-B"];
        assert_eq!(hla_b["geneSymbol"], "HLA-B");
        assert_eq!(hla_b["alleleDefinitionVersion"], Value::Null);
        assert_eq!(hla_b["alleleDefinitionSource"], "UNKNOWN");
        assert_eq!(hla_b["phenotypeVersion"], Value::Null);
        assert_eq!(hla_b["chr"], Value::Null);
        assert_eq!(hla_b["callSource"], "NONE");
        assert_eq!(hla_b["phased"], false);
        assert_eq!(hla_b["effectivelyPhased"], false);
        assert!(
            hla_b["uncalledHaplotypes"]
                .as_array()
                .is_some_and(Vec::is_empty)
        );
        assert!(hla_b["messages"].as_array().is_some_and(Vec::is_empty));
        assert!(
            hla_b["sourceDiplotypes"]
                .as_array()
                .is_some_and(Vec::is_empty)
        );
        assert!(
            hla_b["matcherComponentHaplotypes"]
                .as_array()
                .is_some_and(Vec::is_empty)
        );
        assert!(
            hla_b["matcherHomozygousComponentHaplotypes"]
                .as_array()
                .is_some_and(Vec::is_empty)
        );
        assert!(
            hla_b["recommendationDiplotypes"]
                .as_array()
                .is_some_and(Vec::is_empty)
        );
        assert!(hla_b["variants"].as_array().is_some_and(Vec::is_empty));
        assert!(
            hla_b["variantsOfInterest"]
                .as_array()
                .is_some_and(Vec::is_empty)
        );
        assert_eq!(hla_b["hasUndocumentedVariations"], false);
        assert_eq!(hla_b["treatUndocumentedVariationsAsReference"], false);
        let gene_fields = hla_b.as_object().expect("HLA-B gene object");
        for rust_only_field in [
            "lookupKey",
            "diplotypeKey",
            "phenotypes",
            "activityScore",
            "activityScoreType",
            "allelePresenceType",
            "sourceDiplotype",
            "matchScore",
            "outsideCall",
            "outsidePhenotypeMismatch",
            "outsideActivityScoreMismatch",
        ] {
            assert!(
                !gene_fields.contains_key(rust_only_field),
                "Rust-only field {rust_only_field} leaked into Java GeneReport JSON"
            );
        }
        assert!(
            context
                .gene_report("HLA-B")
                .expect("HLA-B gene report")
                .related_drugs
                .contains(&super::DrugLink::new("allopurinol", "PA448320"))
        );
        let related_drugs = hla_b["relatedDrugs"]
            .as_array()
            .expect("HLA-B relatedDrugs");
        assert!(related_drugs.contains(&serde_json::json!({
            "name": "abacavir",
            "id": "PA448004"
        })));
        assert!(related_drugs.contains(&serde_json::json!({
            "name": "allopurinol",
            "id": "PA448320"
        })));

        let allopurinol = &value["drugs"]["CPIC_GUIDELINE"]["allopurinol"];
        assert_eq!(allopurinol["name"], "allopurinol");
        assert_eq!(allopurinol["source"], "CPIC_GUIDELINE");
        assert!(
            allopurinol["messages"]
                .as_array()
                .is_some_and(Vec::is_empty)
        );
        assert!(
            allopurinol["variants"]
                .as_array()
                .is_some_and(Vec::is_empty)
        );
        let citations = allopurinol["citations"]
            .as_array()
            .expect("allopurinol citations");
        assert_eq!(citations.len(), 2);
        assert_eq!(citations[0]["year"], 2013);
        assert_eq!(citations[0]["pmid"], "23232549");
        assert_eq!(citations[1]["year"], 2016);
        assert_eq!(citations[1]["pmid"], "26094938");
        assert!(allopurinol["guidelines"].is_array());

        let guideline = &allopurinol["guidelines"][0];
        assert_eq!(
            guideline["name"],
            "Annotation of CPIC Guideline for allopurinol and HLA-B"
        );
        let annotation = &guideline["annotations"][0];
        assert_eq!(
            annotation["drugRecommendation"],
            "Allopurinol is contraindicated"
        );
        assert_eq!(annotation["phenotypes"]["HLA-B"], "*58:01 positive");
        assert_eq!(
            annotation["lookupKey"][0]["HLA-B"],
            serde_json::json!("*58:01 positive")
        );
    }

    #[test]
    fn report_context_serializes_message_sets_in_java_compare_to_order() {
        fn message(
            exception_type: &str,
            text: &str,
            gene_match: Option<&str>,
        ) -> MessageAnnotation {
            MessageAnnotation {
                exception_type: exception_type.to_owned(),
                name: "same-json-message".to_owned(),
                version: Some("v1".to_owned()),
                message: text.to_owned(),
                matches: MatchLogic {
                    gene: gene_match.map(str::to_owned),
                    ..MatchLogic::default()
                },
            }
        }

        fn assert_message_order(messages: &Value) {
            let messages = messages.as_array().expect("messages array");
            let order = messages
                .iter()
                .map(|message| {
                    (
                        message["exception_type"].as_str().expect("exception_type"),
                        message["message"].as_str().expect("message"),
                    )
                })
                .collect::<Vec<_>>();
            assert_eq!(
                order,
                [
                    ("ambiguity", "Ambiguity JSON message"),
                    ("footnote", "Footnote JSON message"),
                    ("note", "Note JSON message"),
                ]
            );
        }

        let ordered_messages = [
            message(MessageAnnotation::TYPE_NOTE, "Note JSON message", None),
            message(
                MessageAnnotation::TYPE_AMBIGUITY,
                "Ambiguity JSON message",
                Some("ZZZ"),
            ),
            message(
                MessageAnnotation::TYPE_FOOTNOTE,
                "Footnote JSON message",
                Some("AAA"),
            ),
        ];

        let mut gene = ReportGene::new("GENE1", ["Normal".to_owned()]);
        gene.messages = ordered_messages.clone().into_iter().collect();

        let annotation = AnnotationReport {
            local_id: "annotation-json-order".to_owned(),
            drug_recommendation: Some("Use recommendation".to_owned()),
            classification: "Strong".to_owned(),
            population: crate::phenotype::NA.to_owned(),
            genotypes: vec![RecommendationGenotype::from_report_genes([gene.clone()])],
            implications: Vec::new(),
            phenotypes: BTreeMap::new(),
            activity_scores: BTreeMap::new(),
            highlighted_variants: BTreeSet::new(),
            dosing_information: false,
            alternate_drug_available: false,
            other_prescribing_guidance: false,
            messages: ordered_messages.clone().into_iter().collect(),
            lookup_key: Vec::new(),
        };

        let drug_report = DrugReport {
            name: "json drug".to_owned(),
            id: "RxNorm:json".to_owned(),
            source: PrescribingGuidanceSource::CpicGuideline,
            version: "v-test".to_owned(),
            messages: ordered_messages.into_iter().collect(),
            report_variants: BTreeSet::new(),
            urls: Vec::new(),
            citations: BTreeSet::new(),
            guidelines: [GuidelineReport {
                id: "guideline-json-order".to_owned(),
                name: "Guideline JSON order".to_owned(),
                source: PrescribingGuidanceSource::CpicGuideline,
                version: "v-test".to_owned(),
                url: None,
                genes: ["GENE1".to_owned()].into_iter().collect(),
                report_genes: vec![gene.clone()],
                annotations: [annotation].into_iter().collect(),
            }]
            .into_iter()
            .collect(),
        };

        let context = ReportContext {
            title: None,
            data_version: "v-test".to_owned(),
            gene_reports: [("GENE1".to_owned(), gene)].into_iter().collect(),
            report_gene_sources: BTreeMap::new(),
            drug_reports: [(
                PrescribingGuidanceSource::CpicGuideline,
                [("json drug".to_owned(), drug_report)]
                    .into_iter()
                    .collect(),
            )]
            .into_iter()
            .collect(),
            messages: Vec::new(),
            unannotated_gene_calls: Vec::new(),
        };

        let json = context.to_json_string().expect("report JSON");
        let value = serde_json::from_str::<Value>(&json).expect("report JSON value");

        assert_message_order(&value["genes"]["GENE1"]["messages"]);
        let drug = &value["drugs"]["CPIC_GUIDELINE"]["json drug"];
        assert_message_order(&drug["messages"]);
        assert_message_order(&drug["guidelines"][0]["annotations"][0]["messages"]);
    }

    #[test]
    fn report_context_serializes_annotations_in_java_compare_to_order() {
        fn annotation(
            local_id: &str,
            gene: &ReportGene,
            message_text: &str,
            recommendation: &str,
        ) -> AnnotationReport {
            AnnotationReport {
                local_id: local_id.to_owned(),
                drug_recommendation: Some(recommendation.to_owned()),
                classification: "Strong".to_owned(),
                population: crate::phenotype::NA.to_owned(),
                genotypes: vec![RecommendationGenotype::from_report_genes([gene.clone()])],
                implications: Vec::new(),
                phenotypes: BTreeMap::new(),
                activity_scores: BTreeMap::new(),
                highlighted_variants: BTreeSet::new(),
                dosing_information: false,
                alternate_drug_available: false,
                other_prescribing_guidance: false,
                messages: [MessageAnnotation::new_note(
                    "same-annotation-message",
                    message_text.to_owned(),
                )]
                .into_iter()
                .collect(),
                lookup_key: Vec::new(),
            }
        }

        let gene_a = ReportGene::new("GENEA", ["Normal".to_owned()]);
        let gene_b = ReportGene::new("GENEB", ["Normal".to_owned()]);
        let annotations = [
            annotation("00-b-gene", &gene_b, "Alpha annotation message", "B gene"),
            annotation(
                "zz-alpha-message",
                &gene_a,
                "Alpha annotation message",
                "A alpha",
            ),
            annotation(
                "aa-beta-message",
                &gene_a,
                "Beta annotation message",
                "A beta",
            ),
        ]
        .into_iter()
        .collect();

        let drug_report = DrugReport {
            name: "annotation order drug".to_owned(),
            id: "RxNorm:annotation-order".to_owned(),
            source: PrescribingGuidanceSource::CpicGuideline,
            version: "v-test".to_owned(),
            messages: BTreeSet::new(),
            report_variants: BTreeSet::new(),
            urls: Vec::new(),
            citations: BTreeSet::new(),
            guidelines: [GuidelineReport {
                id: "guideline-annotation-order".to_owned(),
                name: "Guideline annotation order".to_owned(),
                source: PrescribingGuidanceSource::CpicGuideline,
                version: "v-test".to_owned(),
                url: None,
                genes: ["GENEA".to_owned(), "GENEB".to_owned()]
                    .into_iter()
                    .collect(),
                report_genes: vec![gene_a.clone(), gene_b.clone()],
                annotations,
            }]
            .into_iter()
            .collect(),
        };
        let context = ReportContext {
            title: None,
            data_version: "v-test".to_owned(),
            gene_reports: [("GENEA".to_owned(), gene_a), ("GENEB".to_owned(), gene_b)]
                .into_iter()
                .collect(),
            report_gene_sources: BTreeMap::new(),
            drug_reports: [(
                PrescribingGuidanceSource::CpicGuideline,
                [("annotation order drug".to_owned(), drug_report)]
                    .into_iter()
                    .collect(),
            )]
            .into_iter()
            .collect(),
            messages: Vec::new(),
            unannotated_gene_calls: Vec::new(),
        };

        let json = context.to_json_string().expect("report JSON");
        let value = serde_json::from_str::<Value>(&json).expect("report JSON value");
        let annotations = value["drugs"]["CPIC_GUIDELINE"]["annotation order drug"]["guidelines"]
            [0]["annotations"]
            .as_array()
            .expect("annotations");
        let ordered = annotations
            .iter()
            .map(|annotation| {
                (
                    annotation["genotypes"][0]["diplotypes"][0]["geneSymbol"]
                        .as_str()
                        .expect("geneSymbol"),
                    annotation["messages"][0]["message"]
                        .as_str()
                        .expect("message"),
                    annotation["drugRecommendation"]
                        .as_str()
                        .expect("drugRecommendation"),
                )
            })
            .collect::<Vec<_>>();

        assert_eq!(
            ordered,
            [
                ("GENEA", "Alpha annotation message", "A alpha"),
                ("GENEA", "Beta annotation message", "A beta"),
                ("GENEB", "Alpha annotation message", "B gene"),
            ]
        );
    }

    #[test]
    fn recommendation_genotype_serializes_report_genes_in_java_order() {
        let genotype = RecommendationGenotype::from_report_genes([
            ReportGene::new("ZZZ", ["z1".to_owned(), "z2".to_owned()]),
            ReportGene::new("AAA", ["a1".to_owned(), "a2".to_owned()]),
        ]);

        assert_eq!(
            genotype
                .report_genes()
                .iter()
                .map(|gene| gene.gene.as_str())
                .collect::<Vec<_>>(),
            ["AAA", "ZZZ"]
        );
        assert_eq!(
            genotype.lookup_keys(),
            &[
                BTreeMap::from([
                    ("AAA".to_owned(), Value::String("a1".to_owned())),
                    ("ZZZ".to_owned(), Value::String("z1".to_owned())),
                ]),
                BTreeMap::from([
                    ("AAA".to_owned(), Value::String("a2".to_owned())),
                    ("ZZZ".to_owned(), Value::String("z1".to_owned())),
                ]),
                BTreeMap::from([
                    ("AAA".to_owned(), Value::String("a1".to_owned())),
                    ("ZZZ".to_owned(), Value::String("z2".to_owned())),
                ]),
                BTreeMap::from([
                    ("AAA".to_owned(), Value::String("a2".to_owned())),
                    ("ZZZ".to_owned(), Value::String("z2".to_owned())),
                ]),
            ]
        );

        let value = serde_json::to_value(&genotype).expect("genotype JSON");
        let genes = value["diplotypes"]
            .as_array()
            .expect("diplotypes")
            .iter()
            .map(|diplotype| diplotype["geneSymbol"].as_str().expect("geneSymbol"))
            .collect::<Vec<_>>();
        assert_eq!(genes, ["AAA", "ZZZ"]);
    }

    #[test]
    fn publication_sorting_matches_java_publication_compare_to() {
        let publications = BTreeSet::from([
            publication(Some("333"), Some("Later PMID"), Some(2020)),
            publication(Some("111"), Some("Earlier PMID"), Some(2020)),
            publication(None, Some("Missing PMID sorts first"), Some(2020)),
            publication(Some("999"), Some("Older year"), Some(2019)),
        ]);

        let ordered = publications
            .iter()
            .map(|publication| {
                (
                    publication.year,
                    publication.pmid.as_deref(),
                    publication.title.as_deref(),
                )
            })
            .collect::<Vec<_>>();
        assert_eq!(
            ordered,
            vec![
                (Some(2019), Some("999"), Some("Older year")),
                (Some(2020), None, Some("Missing PMID sorts first")),
                (Some(2020), Some("111"), Some("Earlier PMID")),
                (Some(2020), Some("333"), Some("Later PMID")),
            ]
        );
    }

    #[test]
    fn write_report_json_writes_pretty_json_and_rejects_non_json_path_like_java() {
        let collection = PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH))
            .expect("prescribing guidance");
        let context = ReportContext::from_gene_reports(
            &collection,
            [ReportGene::new("HLA-B", ["*58:01 positive".to_owned()])],
            Some("test report".to_owned()),
        );
        let base =
            std::env::temp_dir().join(format!("pharmcat-report-json-{}", std::process::id()));
        fs::create_dir_all(&base).expect("temp dir");
        let json_path = base.join("report.json");
        let txt_path = base.join("report.txt");

        write_report_json(&context, &json_path).expect("write report JSON");
        let json = fs::read_to_string(&json_path).expect("report JSON contents");
        assert!(json.starts_with("{\n  \"title\""));
        assert!(json.contains("\"CPIC_GUIDELINE\""));

        let error = write_report_json(&context, &txt_path).expect_err("invalid extension");
        assert!(error.to_string().contains("Invalid format: "));

        fs::remove_file(json_path).ok();
        fs::remove_dir(base).ok();
    }

    #[test]
    fn calls_only_tsv_uses_java_header_shape_for_current_report_gene_surface() {
        let collection = PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH))
            .expect("prescribing guidance");
        let context = ReportContext::from_gene_reports(
            &collection,
            [ReportGene::new("CYP2D6", ["1.0".to_owned()])
                .with_source_diplotype("*1/*4")
                .with_phenotypes(["Intermediate Metabolizer".to_owned()])
                .with_activity_score("1.0")
                .with_match_score("42")],
            None,
        );
        let options = CallsOnlyTsvOptions {
            pharmcat_version: "test-version".to_owned(),
            sample_id: Some("sample-1".to_owned()),
            show_sample_id: true,
            ..CallsOnlyTsvOptions::default()
        };

        let tsv = calls_only_tsv_string(&context, &options);
        let lines = tsv.lines().collect::<Vec<_>>();
        assert_eq!(lines[0], "PharmCAT test-version");
        assert!(
            lines[1].starts_with("Sample ID\tGene\tSource Diplotype\tPhenotype\tActivity Score")
        );

        let header_fields = lines[1].split('\t').collect::<Vec<_>>();
        assert_eq!(header_fields.len(), 17);
        assert_eq!(header_fields[0], "Sample ID");
        assert_eq!(header_fields[13], "Missing positions");
        assert_eq!(header_fields[16], "Recommendation Lookup Activity Score");

        let row_fields = lines[2].split('\t').collect::<Vec<_>>();
        assert_eq!(row_fields.len(), 17);
        assert_eq!(row_fields[0], "sample-1");
        assert_eq!(row_fields[1], "CYP2D6");
        assert_eq!(row_fields[2], "*1/*4");
        assert_eq!(row_fields[3], "Intermediate Metabolizer");
        assert_eq!(row_fields[4], "1.0");
        assert_eq!(row_fields[11], "no");
        assert_eq!(row_fields[12], "42");
        assert_eq!(row_fields[13], "no");
        assert_eq!(row_fields[15], "1.0");
        assert_eq!(row_fields[16], "1.0");
    }

    #[test]
    fn calls_only_tsv_debug_columns_match_java_variant_flags() {
        let collection = PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH))
            .expect("prescribing guidance");
        let context =
            ReportContext::from_gene_reports(
                &collection,
                [ReportGene::new("NUDT15", ["Normal Metabolizer".to_owned()])
                    .with_variant_reports([
                        VariantReport::new("rs1", Some("A|T"))
                            .with_position(101)
                            .with_reference_allele("A"),
                        VariantReport::new("rs2", None::<String>).with_position(202),
                        VariantReport::new("rs3", Some("G/G"))
                            .with_position(303)
                            .with_reference_allele("G")
                            .with_undocumented_variations(true),
                    ])],
                None,
            );
        let options = CallsOnlyTsvOptions {
            show_variants: true,
            show_missing_variants: true,
            show_undocumented_variants: true,
            treat_undocumented_variations_as_reference: true,
            ..CallsOnlyTsvOptions::default()
        };

        let tsv = calls_only_tsv_string(&context, &options);
        let lines = tsv.lines().collect::<Vec<_>>();
        let header_fields = lines[1].split('\t').collect::<Vec<_>>();
        let row_fields = lines[2].split('\t').collect::<Vec<_>>();

        assert_eq!(header_fields.len(), 18);
        assert_eq!(header_fields[12], CALLS_ONLY_HEADER_VARIANTS);
        assert_eq!(header_fields[13], "Missing positions");
        assert_eq!(header_fields[14], CALLS_ONLY_HEADER_UNDOCUMENTED_VARIANTS);
        assert_eq!(row_fields[12], "101:A|T");
        assert_eq!(row_fields[13], "202");
        assert_eq!(row_fields[14], "303 treat as reference");
    }

    #[test]
    fn calls_only_tsv_appends_sample_metadata_and_single_file_rows_like_java() {
        let collection = PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH))
            .expect("prescribing guidance");
        let context = ReportContext::from_gene_reports(
            &collection,
            [ReportGene::new("CYP2C9", ["Normal Metabolizer".to_owned()])],
            None,
        );
        let sample_properties = BTreeMap::from([
            ("State".to_owned(), "CA".to_owned()),
            ("Town".to_owned(), "Stanford".to_owned()),
        ]);
        let options = CallsOnlyTsvOptions {
            pharmcat_version: "test-version".to_owned(),
            sample_id: Some("sample-1".to_owned()),
            show_sample_id: true,
            single_file_mode: true,
            sample_properties,
            ..CallsOnlyTsvOptions::default()
        };

        let tsv = calls_only_tsv_string(&context, &options);
        let lines = tsv.lines().collect::<Vec<_>>();
        let header_fields = lines[1].split('\t').collect::<Vec<_>>();
        let row_fields = lines[2].split('\t').collect::<Vec<_>>();
        assert_eq!(header_fields[17], "State");
        assert_eq!(header_fields[18], "Town");
        assert_eq!(row_fields[0], "sample-1");
        assert_eq!(row_fields[17], "CA");
        assert_eq!(row_fields[18], "Stanford");

        let base =
            std::env::temp_dir().join(format!("pharmcat-calls-only-single-{}", std::process::id()));
        fs::create_dir_all(&base).expect("temp dir");
        let path = base.join("calls.tsv");
        write_calls_only_tsv(&context, &path, &options).expect("first write");
        write_calls_only_tsv(&context, &path, &options).expect("append write");
        let contents = fs::read_to_string(&path).expect("calls-only contents");
        assert_eq!(contents.matches("PharmCAT test-version").count(), 1);
        assert_eq!(contents.lines().count(), 4);

        fs::remove_file(path).ok();
        fs::remove_dir(base).ok();
    }

    #[test]
    fn write_calls_only_tsv_writes_file_and_rejects_non_tsv_path_like_java() {
        let collection = PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH))
            .expect("prescribing guidance");
        let context = ReportContext::from_gene_reports(
            &collection,
            [ReportGene::new("HLA-B", ["*58:01 positive".to_owned()])
                .with_source_diplotype("*58:01 positive")],
            None,
        );
        let base = std::env::temp_dir().join(format!("pharmcat-report-tsv-{}", std::process::id()));
        fs::create_dir_all(&base).expect("temp dir");
        let tsv_path = base.join("report.tsv");
        let txt_path = base.join("report.txt");

        write_calls_only_tsv(&context, &tsv_path, &CallsOnlyTsvOptions::default())
            .expect("write calls-only TSV");
        let tsv = fs::read_to_string(&tsv_path).expect("calls-only TSV contents");
        assert!(tsv.starts_with("PharmCAT unknown\nGene\tSource Diplotype"));
        assert!(tsv.contains("HLA-B\t*58:01 positive"));

        let error = write_calls_only_tsv(&context, &txt_path, &CallsOnlyTsvOptions::default())
            .expect_err("invalid extension");
        assert!(error.to_string().contains("does not end with .tsv"));

        fs::remove_file(tsv_path).ok();
        fs::remove_dir(base).ok();
    }

    #[test]
    fn loads_html_template_resources_like_java_html_format() {
        let templates = HtmlTemplateSet::from_reporter_dir(Path::new(REPORTER_RESOURCE_DIR))
            .expect("templates");

        assert!(templates.report.starts_with("<!DOCTYPE html>"));
        assert!(templates.report.contains("{{> \"header\"}}"));
        assert!(templates.header.contains("<meta name=\"viewport\""));
        assert!(
            templates
                .uncallable_genes_note
                .contains("could not be called")
        );
        assert!(
            templates
                .disclaimers
                .contains("Disclaimers and Other Information")
        );
    }

    #[test]
    fn html_java_css_selector_matches_report_helpers_sanitize_css_selector() {
        assert_eq!(super::html_java_css_selector("__ABC__"), "ABC");
        assert_eq!(super::html_java_css_selector("(ABC)!"), "ABC");
        assert_eq!(super::html_java_css_selector("A B - C"), "A_B-C");
        assert_eq!(super::html_java_css_selector("MT-RNR1"), "MT-RNR1");
        assert_eq!(super::html_java_css_selector("IFNL3/4"), "IFNL3_4");
    }

    #[test]
    fn html_amd_subtitle_matches_java_variant_gene_rules() {
        let mut abcg2 = ReportGene::new("ABCG2", ["Normal Function".to_owned()]);
        assert_eq!(super::html_amd_subtitle(&abcg2), "Variant Matched");

        abcg2.source_diplotypes = vec![
            ReportDiplotype {
                gene: "ABCG2".to_owned(),
                label: "421C>A".to_owned(),
                ..ReportDiplotype::default()
            },
            ReportDiplotype {
                gene: "ABCG2".to_owned(),
                label: "34G>A".to_owned(),
                ..ReportDiplotype::default()
            },
        ];
        abcg2.outside_call = true;
        assert_eq!(super::html_amd_subtitle(&abcg2), "Variants Reported");

        abcg2.outside_call = false;
        abcg2
            .matcher_component_haplotypes
            .insert("421C>A".to_owned());
        assert_eq!(super::html_amd_subtitle(&abcg2), "Alleles Matched");

        let cyp2d6 = ReportGene::new("CYP2D6", ["Normal Metabolizer".to_owned()]);
        assert_eq!(super::html_amd_subtitle(&cyp2d6), "Allele Matched");
    }

    #[test]
    fn html_amd_allele_function_matches_java_helper_edge_cases() {
        let cyp2c19 = ReportGene::new("CYP2C19", ["*1/*38".to_owned()]);
        let mut functions = BTreeMap::new();
        functions.insert("*1".to_owned(), "Normal function".to_owned());
        functions.insert("*2".to_owned(), crate::phenotype::NA.to_owned());
        functions.insert("*38".to_owned(), "Increased function".to_owned());

        let suppressed_variant = VariantReport::new("rsOther", Some("C/T"));
        assert_eq!(
            super::html_amd_allele_function(&cyp2c19, &suppressed_variant, "*1", &functions),
            None
        );

        let allowed_variant = VariantReport::new("RS3758581", Some("C/T"));
        assert_eq!(
            super::html_amd_allele_function(&cyp2c19, &allowed_variant, "*1", &functions)
                .as_deref(),
            Some("<li>*1 - Normal function</li>")
        );
        assert_eq!(
            super::html_amd_allele_function(&cyp2c19, &allowed_variant, "*2", &functions)
                .as_deref(),
            Some("<li>*2 - Unassigned</li>")
        );
        assert_eq!(
            super::html_amd_allele_function(&cyp2c19, &allowed_variant, "*38", &functions)
                .as_deref(),
            Some("<li>*38 - Increased function</li>")
        );
        assert_eq!(
            super::html_amd_allele_function(&cyp2c19, &allowed_variant, "*99", &functions)
                .as_deref(),
            Some("<li>*99 - Unassigned</li>")
        );
    }

    #[test]
    fn html_amd_no_call_and_phase_status_match_java_helper_edges() {
        let none_with_called_variant = ReportGene::new("CYP2D6", Vec::<String>::new())
            .with_call_source(ReportCallSource::None)
            .with_variant_reports([VariantReport::new("rsCalled", Some("A/G"))]);
        assert!(super::html_amd_no_call(&none_with_called_variant));

        let matcher_empty = ReportGene::new("CYP2D6", ["*1/*1".to_owned()])
            .with_call_source(ReportCallSource::Matcher);
        assert!(super::html_amd_no_call(&matcher_empty));

        let matcher_all_missing = ReportGene::new("CYP2D6", ["*1/*1".to_owned()])
            .with_call_source(ReportCallSource::Matcher)
            .with_variant_reports([VariantReport::new("rsMissing", None::<String>)]);
        assert!(super::html_amd_no_call(&matcher_all_missing));

        let matcher_called = ReportGene::new("CYP2D6", ["*1/*1".to_owned()])
            .with_call_source(ReportCallSource::Matcher)
            .with_variant_reports([VariantReport::new("rsCalled", Some("A/G"))]);
        assert!(!super::html_amd_no_call(&matcher_called));
        assert_eq!(super::html_amd_phase_status(&matcher_called), "Unphased");

        let outside = ReportGene::new("CYP2D6", ["*1/*1".to_owned()])
            .with_call_source(ReportCallSource::Outside);
        assert!(!super::html_amd_no_call(&outside));
        assert_eq!(
            super::html_amd_phase_status(&outside),
            "Unavailable for calls made outside PharmCAT"
        );

        let mut phased = matcher_called.clone();
        phased.phased = true;
        assert_eq!(super::html_amd_phase_status(&phased), "Phased");

        let mut phased_with_phase_set = matcher_called.clone();
        phased_with_phase_set.phased = true;
        phased_with_phase_set.variant_reports =
            vec![VariantReport::new("rsPhaseSet", Some("A|G")).with_phasing(true, Some(11))];
        assert_eq!(
            super::html_amd_phase_status(&phased_with_phase_set),
            "Phased, with phase sets (PS)"
        );

        let mut effectively_phased = matcher_called;
        effectively_phased.effectively_phased = true;
        assert_eq!(
            super::html_amd_phase_status(&effectively_phased),
            "Unphased"
        );
    }

    #[test]
    fn html_amd_uncalled_haps_match_java_helper_edges() {
        let mut no_uncalled = ReportGene::new("CYP2D6", ["*1/*1".to_owned()])
            .with_call_source(ReportCallSource::Matcher)
            .with_variant_reports([VariantReport::new("rsCalled", Some("A/G"))]);
        assert!(no_uncalled.uncalled_haplotypes.is_empty());
        assert_eq!(super::html_amd_total_missing_variants(&no_uncalled), 0);
        assert_eq!(super::html_amd_uncalled_haps(&no_uncalled), "");

        no_uncalled
            .variant_reports
            .push(VariantReport::new("rsMissing", None::<String>));
        no_uncalled.uncalled_haplotypes = ["Unknown", "*10", "*2", "Any"]
            .into_iter()
            .map(str::to_owned)
            .collect();

        assert!(!no_uncalled.uncalled_haplotypes.is_empty());
        assert_eq!(super::html_amd_total_missing_variants(&no_uncalled), 1);
        assert_eq!(
            super::html_amd_uncalled_haps(&no_uncalled),
            "Any, *2, *10, Unknown"
        );
    }

    #[test]
    fn html_amd_messages_match_java_filtering_order_and_message_rendering() {
        let mut ambiguity = MessageAnnotation::new_note("same-rule", "Ambiguity PMID:12345");
        ambiguity.exception_type = MessageAnnotation::TYPE_AMBIGUITY.to_owned();
        ambiguity.version = Some("v1".to_owned());
        ambiguity.matches.gene = Some("Z".to_owned());

        let mut note = MessageAnnotation::new_note("same-rule", "Note PMID:67890");
        note.version = Some("v1".to_owned());

        let mut extra_position =
            MessageAnnotation::new_note("same-rule", "Extra position PMID:11111");
        extra_position.version = Some("v1".to_owned());
        extra_position.exception_type = MessageAnnotation::TYPE_EXTRA_POSITION.to_owned();

        let mut report_as_genotype =
            MessageAnnotation::new_note("same-rule", "Report-as-genotype PMID:22222");
        report_as_genotype.version = Some("v1".to_owned());
        report_as_genotype.exception_type = MessageAnnotation::TYPE_REPORT_AS_GENOTYPE.to_owned();

        let mut footnote = MessageAnnotation::new_note("same-rule", "Footnote PMID:33333");
        footnote.version = Some("v1".to_owned());
        footnote.exception_type = MessageAnnotation::TYPE_FOOTNOTE.to_owned();

        let messages = [
            note,
            extra_position.clone(),
            ambiguity,
            report_as_genotype,
            footnote,
        ]
        .into_iter()
        .collect::<BTreeSet<_>>();

        let amd_messages = messages
            .iter()
            .filter(|message| message.is_message())
            .map(|message| message.message.as_str())
            .collect::<Vec<_>>();
        assert_eq!(amd_messages, ["Ambiguity PMID:12345", "Note PMID:67890"]);

        let extra_position_notes = messages
            .iter()
            .filter(|message| message.is_extra_position_note())
            .map(|message| message.message.as_str())
            .collect::<Vec<_>>();
        assert_eq!(extra_position_notes, ["Extra position PMID:11111"]);

        assert_eq!(super::html_message_class(&extra_position), "same-rule");
        assert_eq!(
            super::html_message_message(
                messages
                    .iter()
                    .find(|message| message.message == "Ambiguity PMID:12345")
                    .expect("ambiguity message")
            ),
            "Ambiguity <a href=\"https://pubmed.ncbi.nlm.nih.gov/12345\" target=\"_blank\" rel=\"noopener noreferrer\">PMID:12345</a>"
        );
    }

    #[test]
    fn html_amd_calls_at_positions_uses_activity_value_function_map() {
        let phenotype =
            read_gene_phenotype_file(Path::new(CYP2D6_PHENOTYPE_PATH)).expect("CYP2D6 phenotype");
        let mut cyp2d6 = ReportGene::new("CYP2D6", ["*1/*10".to_owned()]);
        cyp2d6.allele_function_map = phenotype.formatted_function_score_map();
        cyp2d6.variant_reports =
            vec![VariantReport::new("rs1065852", Some("C/T")).with_alleles(["*10".to_owned()])];

        let html = super::html_amd_calls_at_positions(&cyp2d6);

        assert!(html.contains("<li>*10 - Activity Value 0.25 (Decreased function)</li>"));
        assert!(!html.contains("<li>*10 - Decreased function</li>"));
        assert!(!html.contains("<li>*10 - Unassigned</li>"));
    }

    #[test]
    fn report_html_string_writes_minimal_report_from_current_report_context() {
        let collection = PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH))
            .expect("prescribing guidance");
        let mut hla_b = ReportGene::new("HLA-B", ["*58:01 positive".to_owned()])
            .with_source_diplotype("*58:01 positive");
        hla_b.variant_reports =
            vec![VariantReport::new("rsHlaB", Some("A/T")).with_alleles(["*58:01".to_owned()])];
        hla_b.recommendation_diplotypes = vec![ReportDiplotype {
            gene: "HLA-B".to_owned(),
            label: "*58:01 positive".to_owned(),
            phenotypes: vec!["*58:01 positive".to_owned()],
            ..ReportDiplotype::default()
        }];
        let context =
            ReportContext::from_gene_reports(&collection, [hla_b], Some("sample <one>".to_owned()));
        let html = report_html_string(
            &context,
            &HtmlReportOptions {
                pharmcat_version: Some("v-test".to_owned()),
                data_version: Some("data-test".to_owned()),
                timestamp: Some("June 01, 2026".to_owned()),
                debug: false,
                compact: false,
                definition_genes: BTreeSet::new(),
                no_data_genes: BTreeSet::new(),
            },
        );

        assert!(html.starts_with("<!DOCTYPE html>"));
        assert!(html.contains("<title>PharmCAT Report [sample &lt;one&gt;]</title>"));
        assert!(html.contains("<th>PharmCAT Version</th><td>v-test</td>"));
        assert!(html.contains("<th>Data Version</th><td>data-test</td>"));
        assert!(html.contains("<h2>Section I: Genotype Summary</h2>"));
        assert_eq!(
            html.lines()
                .find(|line| line.contains("Genotypes called:"))
                .map(str::trim),
            Some("<p>Genotypes called: 1 / 1 </p>")
        );
        assert!(html.contains(
            "<th>Drugs</th><th>Gene</th><th>Genotypes<table class=\"diplotype table-small\">"
        ));
        assert!(html.contains("<td>Genotype</td><td>Allele Functionality</td><td>Phenotype</td>"));
        assert!(html.contains("<tr class=\"top-aligned gs-HLA-B\">"));
        assert!(html.contains("<a href=\"#HLA-B\" class=\"normalWrap\">HLA-B</a>"));
        assert!(html.contains("<a href=\"#allopurinol\">allopurinol</a>"));
        assert!(html.contains("<tr class=\"top-aligned gs-dip\"><td>*58:01 positive</td><td>N/A</td><td>*58:01 positive</td></tr>"));
        assert!(html.contains("<div class=\"footnote\">CPIC terms for allele function and phenotype are used for all CPIC genes. For non-CPIC genes, DPWG terms are used.</div>"));
        assert!(html.contains("<div class=\"footnote\">For a full list of disclaimers and limitations see <a href=\"#disclaimer\">Section IV</a>.</div>"));
        assert!(html.contains("Allopurinol is contraindicated"));
        assert!(html.contains("<div class=\"citations\"><p>Citations:</p><ul>"));
        assert!(html.contains("PMID:23232549"));
        assert!(html.contains("PMID:26094938"));
        assert!(html.contains("Section IV: Disclaimers and Other Information"));

        let compact_html = report_html_string(
            &context,
            &HtmlReportOptions {
                compact: true,
                ..HtmlReportOptions::default()
            },
        );
        assert!(compact_html.contains(
            "<a href=\"#allopurinol\">allopurinol</a></span><span class=\"drugTags\"><div class=\"tag\">Alternate Drug</div>"
        ));
    }

    #[test]
    fn report_html_string_renders_genotype_summary_combo_and_context_messages_like_java() {
        let collection = PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH))
            .expect("prescribing guidance");
        let mut combo_message =
            MessageAnnotation::new_note("combo-test", "gene combo summary marker");
        combo_message.exception_type = MessageAnnotation::TYPE_COMBO.to_owned();
        let context_message = MessageAnnotation::new_note(
            "summary test/message",
            "Summary context message <b>escaped</b> PMID:12345",
        );
        let mut cyp2b6 = ReportGene::new("CYP2B6", ["Normal Metabolizer".to_owned()])
            .with_source_diplotype("*1/*1")
            .with_messages([combo_message]);
        cyp2b6.variant_reports =
            vec![VariantReport::new("rsCyp2b6", Some("C/T")).with_alleles(["*1".to_owned()])];
        cyp2b6.recommendation_diplotypes = vec![ReportDiplotype {
            gene: "CYP2B6".to_owned(),
            label: "*1/*1".to_owned(),
            phenotypes: vec!["Normal Metabolizer".to_owned()],
            ..ReportDiplotype::default()
        }];
        let context = ReportContext::from_gene_reports(&collection, [cyp2b6], None)
            .with_messages([context_message]);

        let html = report_html_string(&context, &HtmlReportOptions::default());

        let disclaimer_index = html
            .find("For a full list of disclaimers and limitations see")
            .expect("disclaimer footnote");
        let combo_index = html
            .find("<div class=\"alert alert-info\">Partial and combination allele calls are based on the variants identified in the VCF file.")
            .expect("combo summary alert");
        let summary_index = html
            .find("<div class=\"alert alert-info summary-test-message\">Summary context message &lt;b&gt;escaped&lt;/b&gt; PMID:12345</div>")
            .expect("context summary message");

        assert!(disclaimer_index < combo_index);
        assert!(combo_index < summary_index);
        assert!(!html.contains("https://pubmed.ncbi.nlm.nih.gov/12345"));
    }

    #[test]
    fn report_html_string_renders_no_summary_fallbacks_like_java_report_template() {
        let collection = PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH))
            .expect("prescribing guidance");
        let no_data_context = ReportContext::from_gene_reports(
            &collection,
            [ReportGene::new("CYP2D6", Vec::<String>::new())],
            None,
        );
        let no_data_html = report_html_string(&no_data_context, &HtmlReportOptions::default());

        assert!(no_data_html.contains("<h2>Section I: Genotype Summary</h2>"));
        assert!(no_data_html.contains("<p>No data provided.</p>"));
        assert!(no_data_html.contains("For a full list of disclaimers and limitations, see <a href=\"#disclaimer\">Section IV</a>."));
        assert!(!no_data_html.contains("<table class=\"genotypeSummary\">"));
        assert!(!no_data_html.contains("CPIC terms for allele function and phenotype"));

        let mut uncallable = ReportGene::new("CYP2C19", Vec::<String>::new());
        uncallable.variant_reports =
            vec![VariantReport::new("rsCalled", Some("C/T")).with_alleles(["*1".to_owned()])];
        let uncallable_context = ReportContext::from_gene_reports(&collection, [uncallable], None);
        let uncallable_html =
            report_html_string(&uncallable_context, &HtmlReportOptions::default());

        assert!(uncallable_html.contains("<p>No genotypes called.</p>"));
        assert!(uncallable_html.contains("<p>The following gene could not be called because there were genetic variations that do not match the allele definition. There could still be actionable variants in these genes. See <a href=\"#section-iii\">Section III</a> for details.</p><ul class=\"mb-2\"><li><a id=\"gs-uncallable-CYP2C19\" href=\"#CYP2C19\">CYP2C19</a></li></ul>"));
        assert!(!uncallable_html.contains("<table class=\"genotypeSummary\">"));
        assert!(!uncallable_html.contains("<div class=\"alert alert-warning\"><p>No genotypes"));
    }

    #[test]
    fn report_html_string_renders_genotypes_called_count_like_java_html_format() {
        let collection = PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH))
            .expect("prescribing guidance");
        let mut called = ReportGene::new("HLA-B", ["*58:01 positive".to_owned()])
            .with_source_diplotype("*58:01 positive");
        called.variant_reports =
            vec![VariantReport::new("rsHlaB", Some("A/T")).with_alleles(["*58:01".to_owned()])];
        called.recommendation_diplotypes = vec![ReportDiplotype {
            gene: "HLA-B".to_owned(),
            label: "*58:01 positive".to_owned(),
            phenotypes: vec!["*58:01 positive".to_owned()],
            ..ReportDiplotype::default()
        }];
        let no_data_definition_gene = ReportGene::new("CYP2D6", Vec::<String>::new());
        let subset_removed_gene = ReportGene::new("RYR1", Vec::<String>::new());
        let context = ReportContext::from_gene_reports(
            &collection,
            [called, no_data_definition_gene, subset_removed_gene],
            None,
        );
        let options = HtmlReportOptions {
            definition_genes: ["CYP2D6".to_owned()].into_iter().collect(),
            ..HtmlReportOptions::default()
        };

        assert!(context.gene_report("CYP2D6").is_some());
        assert_eq!(super::html_no_data_genes(&context, &options).len(), 1);
        assert!(super::html_related_gene_symbols(&context).contains("HLA-B"));
        assert!(
            context
                .gene_report("HLA-B")
                .expect("HLA-B")
                .is_reportable_like_java()
        );
        assert_eq!(super::html_genotype_summary_called_genes(&context), 1);
        assert_eq!(
            super::html_genotype_summary_total_genes(&context, &options),
            2
        );

        let html = report_html_string(&context, &options);

        assert_eq!(
            html.lines()
                .find(|line| line.contains("Genotypes called:"))
                .map(str::trim),
            Some("<p>Genotypes called: 1 / 2 </p>")
        );
        assert!(html.contains("<tr class=\"top-aligned gs-HLA-B\">"));
        assert!(!html.contains("<tr class=\"top-aligned gs-CYP2D6\">"));
        assert!(!html.contains("<tr class=\"top-aligned gs-RYR1\">"));
    }

    #[test]
    fn report_html_string_omits_non_reportable_summary_rows_like_java_html_format() {
        let collection = PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH))
            .expect("prescribing guidance");
        let non_reportable = ReportGene::new("CYP2B6", ["Normal Metabolizer".to_owned()])
            .with_source_diplotype("*1/*1");
        let context = ReportContext::from_gene_reports(&collection, [non_reportable], None);

        let html = report_html_string(&context, &HtmlReportOptions::default());

        assert!(html.contains("<p>No data provided.</p>"));
        assert!(!html.contains("<table class=\"genotypeSummary\">"));
        assert!(!html.contains("<tr class=\"top-aligned gs-CYP2B6\">"));
        assert!(!html.contains("<p>Genotypes called:"));
    }

    #[test]
    fn report_html_string_derives_no_data_genes_like_java_html_format() {
        let collection = PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH))
            .expect("prescribing guidance");
        let no_data = ReportGene::new("CYP2D6", Vec::<String>::new());
        let subset_removed = ReportGene::new("RYR1", Vec::<String>::new());
        let context =
            ReportContext::from_gene_reports(&collection, [no_data, subset_removed], None);

        let html = report_html_string(
            &context,
            &HtmlReportOptions {
                definition_genes: ["CYP2D6".to_owned()].into_iter().collect(),
                ..HtmlReportOptions::default()
            },
        );

        assert!(html.contains("<section id=\"section-iii\">"));
        assert!(html.contains("<p class=\"noGeneData\">No data provided for <span class=\"gene cyp2d6\"><span class=\"no-data\">CYP2D6</span></span>.</p>"));
        assert!(!html.contains("gene ryr1"));
    }

    #[test]
    fn report_html_string_renders_initial_section_iii_allele_match_data_like_java() {
        let collection = PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH))
            .expect("prescribing guidance");
        let mut cyp2d6 = ReportGene::new("CYP2D6", ["Normal Metabolizer".to_owned()]);
        cyp2d6.source_diplotypes = vec![ReportDiplotype {
            gene: "CYP2D6".to_owned(),
            label: "*1/*2".to_owned(),
            allele1: Some(super::ReportHaplotype {
                gene: "CYP2D6".to_owned(),
                name: "*1".to_owned(),
                function: "Normal function".to_owned(),
                reference: true,
                activity_value: Some("1.0".to_owned()),
            }),
            allele2: Some(super::ReportHaplotype {
                gene: "CYP2D6".to_owned(),
                name: "*2".to_owned(),
                function: "Decreased function".to_owned(),
                reference: false,
                activity_value: Some("0.5".to_owned()),
            }),
            phenotypes: vec!["Normal Metabolizer".to_owned()],
            ..ReportDiplotype::default()
        }];
        cyp2d6.recommendation_diplotypes = cyp2d6.source_diplotypes.clone();
        cyp2d6.call_source = ReportCallSource::Matcher;
        cyp2d6.phased = true;
        cyp2d6.treat_undocumented_variations_as_reference = true;
        let mut cyp2d6_variant = VariantReport::new("rsTest", Some("A/G"))
            .with_position(123)
            .with_reference_allele("A")
            .with_alleles(["*2".to_owned()])
            .with_phasing(true, Some(7))
            .with_undocumented_variations(true);
        cyp2d6_variant.chromosome = Some("chr22".to_owned());
        cyp2d6_variant.warnings = vec!["Low depth".to_owned()];
        let mut missing_cyp2d6_variant = VariantReport::new("rsMissingCall", None::<String>)
            .with_position(124)
            .with_reference_allele("C");
        missing_cyp2d6_variant.chromosome = Some("chr22".to_owned());
        let mut long_call_variant =
            VariantReport::new("rsLongCall", Some("ACGTACGTAC/TTTTTTTTTTT"))
                .with_position(125)
                .with_reference_allele("GGGGGGGGGG")
                .with_phasing(true, Some(99));
        long_call_variant.chromosome = Some("chr22".to_owned());
        cyp2d6.variant_reports = vec![cyp2d6_variant, missing_cyp2d6_variant, long_call_variant];
        let mut extra_position_note =
            MessageAnnotation::new_note("extra-position-note", "Extra position note PMID:67890");
        extra_position_note.exception_type = MessageAnnotation::TYPE_EXTRA_POSITION.to_owned();
        let mut variant_of_interest = VariantReport::new("rsInterest", Some("T/C"))
            .with_position(456)
            .with_reference_allele("T");
        variant_of_interest.chromosome = Some("chr22".to_owned());
        let mut missing_variant_of_interest =
            VariantReport::new("rsInterestMissing", None::<String>).with_position(457);
        missing_variant_of_interest.chromosome = Some("chr22".to_owned());
        cyp2d6.variant_of_interest_reports = vec![variant_of_interest, missing_variant_of_interest];
        cyp2d6.messages = [
            MessageAnnotation::new_note("pcat-test-note", "Gene-level AMD note PMID:12345"),
            extra_position_note,
        ]
        .into_iter()
        .collect();

        let mut cyp2c19 = ReportGene::new("CYP2C19", Vec::<String>::new());
        cyp2c19.call_source = ReportCallSource::Matcher;
        cyp2c19.variant_reports = vec![
            VariantReport::new("rsCalled", Some("C/T")).with_alleles(["*1".to_owned()]),
            VariantReport::new("rsMissing", None::<String>),
        ];
        cyp2c19.uncalled_haplotypes = ["*2".to_owned()].into_iter().collect();

        let context = ReportContext::from_gene_reports(&collection, [cyp2d6, cyp2c19], None);
        let html = report_html_string(&context, &HtmlReportOptions::default());

        assert!(html.contains("<section id=\"section-iii\">"));
        assert!(html.contains("<h2>Section III: Allele Matching Details</h2>"));
        assert!(html.contains("<li><a href=\"#CYP2C19\">CYP2C19 allele match data</a></li>"));
        assert!(html.contains("<li><a href=\"#CYP2D6\">CYP2D6 allele match data</a></li>"));
        assert!(html.contains("<tr class=\"top-aligned gs-CYP2D6\">"));
        assert!(!html.contains("<tr class=\"top-aligned gs-CYP2C19\">"));
        assert!(html.contains("<tr class=\"top-aligned gs-dip\"><td>*1/*2</td><td>One Decreased function allele and one Normal function allele</td><td>Normal Metabolizer</td></tr>"));
        assert!(html.contains(
            "<p class=\"tdNote\" id=\"gs-undocVarAsRef-CYP2D6\">There are genetic variations in this gene that do not match what is in the allele definition. <b>These undocumented variations were replaced with reference.</b> See <a href=\"#CYP2D6\">Section III</a> for details.</p>"
        ));
        assert!(html.contains("<section class=\"gene CYP2D6\">"));
        assert!(html.contains("<h3 id=\"CYP2D6\">CYP2D6 allele match data</h3>"));
        assert!(html.contains("<th style=\"width: 12em;\">Allele Matched:</th>"));
        assert!(html.contains("<td class=\"top-aligned genotype-result\">*1/*2"));
        assert!(html.contains("<b>These undocumented variations were replaced with reference.</b>  See below for details."));
        assert!(html.contains("<th>Phasing Status:</th>"));
        assert!(html.contains("<p>Phased, with phase sets (PS)</p>"));
        assert!(html.contains("<div class=\"alert alert-warning pcat-test-note\">Gene-level AMD note <a href=\"https://pubmed.ncbi.nlm.nih.gov/12345\" target=\"_blank\" rel=\"noopener noreferrer\">PMID:12345</a></div>"));
        assert!(html.contains("<h4>Calls at Positions</h4>"));
        assert!(html.contains("<tr id=\"chr22_123\" style=\"vertical-align: initial;\">"));
        assert!(html.contains("<td>chr22:123</td>"));
        assert!(html.contains("<td id=\"rsTest\">rsTest</td>"));
        assert!(html.contains(
            "<td class=\"nonwild mismatch\">A/G (PS:7)<div class=\"callMessage\">Undocumented variation</div></td>"
        ));
        assert!(html.contains("<td>\n            A\n          </td>"));
        assert!(html.contains("<li>*2 - Decreased function</li>"));
        assert!(html.contains("<ul class=\"warningList\">"));
        assert!(html.contains("<li>Low depth</li>"));
        assert!(html.contains("<tr id=\"chr22_124\" style=\"vertical-align: initial;\">"));
        assert!(html.contains("<td id=\"rsMissingCall\">rsMissingCall</td>"));
        assert!(html.contains(
            "<td class=\"missingVariant\"><div class=\"callMessage\">Missing</div></td>"
        ));
        assert!(html.contains("<td>\n            C\n          </td>"));
        assert!(html.contains("<tr id=\"chr22_125\" style=\"vertical-align: initial;\">"));
        assert!(html.contains("<td id=\"rsLongCall\">rsLongCall</td>"));
        assert!(html.contains(
            "<td class=\"nonwild\">ACGTACGTA<br />C/<br />TTTTTTTTT<br />TT (PS:99)</td>"
        ));
        assert!(html.contains("<td>\n            GGGGGGGGG<br />G\n          </td>"));
        assert!(html.contains("<h4>Other Positions of Interest</h4>"));
        assert!(html.contains("<div class=\"alert alert-warning\">Extra position note <a href=\"https://pubmed.ncbi.nlm.nih.gov/67890\" target=\"_blank\" rel=\"noopener noreferrer\">PMID:67890</a></div>"));
        assert!(html.contains("<td id=\"chr22_456\">chr22:456</td>"));
        assert!(html.contains("<td id=\"rsInterest\">rsInterest</td>"));
        assert!(html.contains("<td class=\"nonwild\">T/C</td>"));
        assert!(html.contains("<td id=\"chr22_457\">chr22:457</td>"));
        assert!(html.contains("<td id=\"rsInterestMissing\">rsInterestMissing</td>"));
        assert!(html.contains("<td class=\"missingVariant\"><em>missing</em></td>"));
        assert!(html.contains("<section class=\"gene CYP2C19\">"));
        assert!(html.contains("<div class=\"alert alert-warning\"><p>The following gene could not be called because there were genetic variations that do not match the allele definition. There could still be actionable variants in these genes. See <a href=\"#section-iii\">Section III</a> for details.</p><ul class=\"mb-2\"><li><a id=\"gs-uncallable-CYP2C19\" href=\"#CYP2C19\">CYP2C19</a></li></ul></div>"));
        assert!(html.contains("<div class=\"footnote\" id=\"genotypes-dagger\"><sup>&dagger;</sup> Check <a href=\"#section-iii\">Section III</a> for more details about this call.</div>"));
        assert!(html.contains("<td class=\"top-aligned genotype-result\">Not called</td>"));
        assert!(!html.contains("<li>*1 - "));
        assert!(html.contains("<th>Alleles Not Considered:</th>"));
        assert!(html.contains(
            "The following alleles are not considered due to 1 missing positions of the total 2 positions: *2"
        ));
        assert!(html.contains("PharmCAT reports the genotype(s) that receive the highest score"));
    }

    #[test]
    fn report_html_string_filters_section_iii_gene_reports_like_java_compact_mode() {
        let collection = PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH))
            .expect("prescribing guidance");
        let related_uncallable = ReportGene::new("CYP2C19", Vec::<String>::new())
            .with_variant_reports([VariantReport::new("rsRelated", Some("C/T"))]);
        let unrelated = ReportGene::new("GENEX", Vec::<String>::new())
            .with_variant_reports([VariantReport::new("rsUnrelated", Some("A/G"))]);
        let context =
            ReportContext::from_gene_reports(&collection, [related_uncallable, unrelated], None);

        let extended_html = report_html_string(&context, &HtmlReportOptions::default());
        let compact_html = report_html_string(
            &context,
            &HtmlReportOptions {
                compact: true,
                ..HtmlReportOptions::default()
            },
        );

        assert!(extended_html.contains("CYP2C19 allele match data"));
        assert!(extended_html.contains("GENEX allele match data"));
        assert!(compact_html.contains("CYP2C19 allele match data"));
        assert!(!compact_html.contains("GENEX allele match data"));
    }

    #[test]
    fn report_html_string_renders_lowest_function_genotype_summary_components_like_java() {
        let mut dpyd = ReportGene::new("DPYD", ["1.0".to_owned()]);
        dpyd.matcher_component_haplotypes = ["c.1129-5923C>G".to_owned(), "c.1236G>A".to_owned()]
            .into_iter()
            .collect();
        dpyd.matcher_component_diplotypes = vec![
            ReportDiplotype {
                gene: "DPYD".to_owned(),
                label: "c.1236G>A".to_owned(),
                allele1: Some(super::ReportHaplotype {
                    gene: "DPYD".to_owned(),
                    name: "c.1236G>A".to_owned(),
                    function: "Normal function".to_owned(),
                    reference: false,
                    activity_value: Some("1.0".to_owned()),
                }),
                ..ReportDiplotype::default()
            },
            ReportDiplotype {
                gene: "DPYD".to_owned(),
                label: "c.1129-5923C>G".to_owned(),
                allele1: Some(super::ReportHaplotype {
                    gene: "DPYD".to_owned(),
                    name: "c.1129-5923C>G".to_owned(),
                    function: "Decreased function".to_owned(),
                    reference: false,
                    activity_value: Some("0.5".to_owned()),
                }),
                ..ReportDiplotype::default()
            },
        ];
        sort_report_diplotypes(&mut dpyd.matcher_component_diplotypes);
        dpyd.recommendation_diplotypes = vec![ReportDiplotype {
            gene: "DPYD".to_owned(),
            label: "[c.1129-5923C>G + c.1236G>A]/Reference".to_owned(),
            allele1: Some(super::ReportHaplotype {
                gene: "DPYD".to_owned(),
                name: "[c.1129-5923C>G + c.1236G>A]".to_owned(),
                function: "Decreased function".to_owned(),
                reference: false,
                activity_value: Some("0.5".to_owned()),
            }),
            allele2: Some(super::ReportHaplotype {
                gene: "DPYD".to_owned(),
                name: "Reference".to_owned(),
                function: "Normal function".to_owned(),
                reference: true,
                activity_value: Some("1.0".to_owned()),
            }),
            combination: true,
            ..ReportDiplotype::default()
        }];
        dpyd.source_diplotypes = vec![
            ReportDiplotype {
                gene: "DPYD".to_owned(),
                label: "c.1129-5923C>G".to_owned(),
                allele1: Some(super::ReportHaplotype {
                    gene: "DPYD".to_owned(),
                    name: "c.1129-5923C>G".to_owned(),
                    function: "Decreased function".to_owned(),
                    reference: false,
                    activity_value: Some("0.5".to_owned()),
                }),
                ..ReportDiplotype::default()
            },
            ReportDiplotype {
                gene: "DPYD".to_owned(),
                label: "c.1236G>A".to_owned(),
                allele1: Some(super::ReportHaplotype {
                    gene: "DPYD".to_owned(),
                    name: "c.1236G>A".to_owned(),
                    function: "Normal function".to_owned(),
                    reference: false,
                    activity_value: Some("1.0".to_owned()),
                }),
                ..ReportDiplotype::default()
            },
        ];

        let context = ReportContext {
            title: None,
            data_version: "v-test".to_owned(),
            gene_reports: [("DPYD".to_owned(), dpyd.clone())].into_iter().collect(),
            report_gene_sources: [("DPYD".to_owned(), vec![dpyd])].into_iter().collect(),
            drug_reports: BTreeMap::new(),
            messages: Vec::new(),
            unannotated_gene_calls: Vec::new(),
        };
        let html = report_html_string(&context, &HtmlReportOptions::default());

        assert!(html.contains("<tr class=\"lowestFunctionDiplotype\"><td colspan=\"3\"><b class=\"gs-dip_lowestFunction\">[c.1129-5923C&gt;G + c.1236G&gt;A]/Reference</b><br /></td></tr>"));
        assert!(html.contains("<tr class=\"top-aligned gs-dip_component\"><td>c.1129-5923C&gt;G</td><td>Decreased function</td><td rowspan=\"2\" class=\"center\">See Drug Recommendation</td></tr>"));
        assert!(html.contains("<tr class=\"top-aligned gs-dip_component\"><td>c.1236G&gt;A</td><td>Normal function</td></tr>"));

        let json = context.to_json_string().expect("report JSON");
        let value = serde_json::from_str::<Value>(&json).expect("report JSON value");
        let components = value["genes"]["DPYD"]["matcherComponentHaplotypes"]
            .as_array()
            .expect("matcherComponentHaplotypes");
        assert_eq!(components.len(), 2);
        assert_eq!(components[0]["label"], "c.1129-5923C>G");
        assert_eq!(components[0]["allele1"]["name"], "c.1129-5923C>G");
        assert_eq!(components[0]["allele1"]["function"], "Decreased function");
        assert_eq!(components[1]["label"], "c.1236G>A");
        assert_eq!(components[1]["allele1"]["name"], "c.1236G>A");
        assert_eq!(components[1]["allele1"]["function"], "Normal function");
    }

    #[test]
    fn report_html_string_renders_drug_messages_and_report_as_genotype_like_java_report_template() {
        let collection = PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH))
            .expect("prescribing guidance");
        let catalog = MessageCatalog::from_path(Path::new(MESSAGES_PATH)).expect("messages");
        let mut context = ReportContext::from_gene_reports(
            &collection,
            [ReportGene::new("CYP2C9", ["*1/*1".to_owned()])
                .with_variant_reports([VariantReport::new("rs12777823", Some("C|T"))])],
            Some("warfarin sample".to_owned()),
        );
        context.apply_matching_drug_messages(&catalog);

        let html = report_html_string(&context, &HtmlReportOptions::default());

        assert!(html.contains("Please follow the flow chart in figure 2"));
        assert!(html.contains("<td colspan=\"4\"><div class=\"warfarinFlowchart\"><img src=\"https://files.cpicpgx.org/images/warfarin/warfarin_recommendation_diagram.png\" alt=\"Figure 2 from the CPIC guideline for warfarin\"/></div></td>"));
        assert!(!html.contains("<div class=\"reportVariants\">"));
        assert!(!html.contains("<li>rs12777823</li>"));
        assert!(html.contains("<span class=\"rx-hl-var\">rs12777823:C/T</span>"));
        assert!(html.contains("PMID:21900891"));
        assert!(html.contains("PMID:28198005"));

        let mut annotation = AnnotationReport::for_cpic_warfarin(&[]);
        annotation.messages = [
            MessageAnnotation {
                exception_type: MessageAnnotation::TYPE_NOTE.to_owned(),
                name: "same-warfarin-message".to_owned(),
                version: Some("v1".to_owned()),
                message: "Warfarin note message".to_owned(),
                matches: MatchLogic::default(),
            },
            MessageAnnotation {
                exception_type: MessageAnnotation::TYPE_AMBIGUITY.to_owned(),
                name: "same-warfarin-message".to_owned(),
                version: Some("v1".to_owned()),
                message: "Warfarin ambiguity message".to_owned(),
                matches: MatchLogic {
                    gene: Some("ZZZ".to_owned()),
                    ..MatchLogic::default()
                },
            },
            MessageAnnotation {
                exception_type: MessageAnnotation::TYPE_FOOTNOTE.to_owned(),
                name: "footnote.warfarin.test".to_owned(),
                version: None,
                message: "Warfarin footnote".to_owned(),
                matches: MatchLogic::default(),
            },
        ]
        .into_iter()
        .collect();
        let warfarin_cell = super::html_cpic_warfarin_recommendation_cell(&annotation);
        assert!(warfarin_cell.contains(
            "<div class=\"alert alert-info same-warfarin-message\">Warfarin ambiguity message</div>"
        ));
        assert!(warfarin_cell.contains(
            "<div class=\"alert alert-info same-warfarin-message\">Warfarin note message</div>"
        ));
        assert!(
            warfarin_cell
                .find("Warfarin ambiguity message")
                .expect("warfarin ambiguity message")
                < warfarin_cell
                    .find("Warfarin note message")
                    .expect("warfarin note message")
        );
        assert!(!warfarin_cell.contains("Warfarin footnote"));
        assert!(warfarin_cell.contains("warfarin_recommendation_diagram.png"));
    }

    #[test]
    fn report_html_string_groups_source_rows_like_java_recommendation_model() {
        let mut drug_reports =
            BTreeMap::<PrescribingGuidanceSource, BTreeMap<String, DrugReport>>::new();
        let mut cyp2d6 = ReportGene::new("CYP2D6", ["Normal Metabolizer".to_owned()]);
        cyp2d6.recommendation_diplotypes = vec![ReportDiplotype {
            gene: "CYP2D6".to_owned(),
            label: "*2/*1".to_owned(),
            inferred: true,
            inferred_source_diplotypes: vec![ReportDiplotype {
                gene: "CYP2D6".to_owned(),
                label: "*4/*1".to_owned(),
                ..ReportDiplotype::default()
            }],
            ..ReportDiplotype::default()
        }];
        let mut dpyd = ReportGene::new("DPYD", ["1.5".to_owned()]);
        dpyd.recommendation_diplotypes = vec![ReportDiplotype {
            gene: "DPYD".to_owned(),
            label: "*1/*2A".to_owned(),
            inferred: true,
            ..ReportDiplotype::default()
        }];
        let annotation = AnnotationReport {
            local_id: "CPIC-test-1".to_owned(),
            drug_recommendation: Some("Use the CPIC recommendation".to_owned()),
            classification: "Strong".to_owned(),
            population: "general".to_owned(),
            genotypes: vec![
                RecommendationGenotype::from_report_genes([cyp2d6.clone(), dpyd.clone()]),
                RecommendationGenotype::from_report_genes([ReportGene::new(
                    "GENE5",
                    ["Normal Function".to_owned()],
                )]),
            ],
            implications: vec![
                "increased exposure".to_owned(),
                crate::phenotype::NA.to_owned(),
            ],
            phenotypes: [
                ("GENE1".to_owned(), "intermediate metabolizer".to_owned()),
                ("GENE2".to_owned(), crate::phenotype::NA.to_owned()),
            ]
            .into_iter()
            .collect(),
            activity_scores: [("GENE1".to_owned(), String::new())].into_iter().collect(),
            highlighted_variants: ["rsTest:A/T".to_owned()].into_iter().collect(),
            dosing_information: true,
            alternate_drug_available: true,
            other_prescribing_guidance: true,
            messages: BTreeSet::new(),
            lookup_key: Vec::new(),
        };
        let mut no_action_annotation = annotation.clone();
        no_action_annotation.dosing_information = false;
        no_action_annotation.alternate_drug_available = false;
        no_action_annotation.other_prescribing_guidance = false;
        assert_eq!(
            super::html_annotation_tags(&no_action_annotation),
            "<div class=\"tag noAction\">No Action</div>"
        );
        let cpic_report = DrugReport {
            name: "test drug".to_owned(),
            id: "RxNorm:test".to_owned(),
            source: PrescribingGuidanceSource::CpicGuideline,
            version: "v-test".to_owned(),
            messages: [
                MessageAnnotation {
                    exception_type: MessageAnnotation::TYPE_NOTE.to_owned(),
                    name: "same-source-message".to_owned(),
                    version: Some("v1".to_owned()),
                    message: "Source-level note message".to_owned(),
                    matches: MatchLogic::default(),
                },
                MessageAnnotation {
                    exception_type: MessageAnnotation::TYPE_AMBIGUITY.to_owned(),
                    name: "same-source-message".to_owned(),
                    version: Some("v1".to_owned()),
                    message: "Source-level ambiguity message".to_owned(),
                    matches: MatchLogic {
                        gene: Some("ZZZ".to_owned()),
                        ..MatchLogic::default()
                    },
                },
                MessageAnnotation {
                    exception_type: MessageAnnotation::TYPE_FOOTNOTE.to_owned(),
                    name: "same-source-footnote".to_owned(),
                    version: Some("v1".to_owned()),
                    message: "Beta source-level footnote".to_owned(),
                    matches: MatchLogic::default(),
                },
                MessageAnnotation {
                    exception_type: MessageAnnotation::TYPE_FOOTNOTE.to_owned(),
                    name: "same-source-footnote".to_owned(),
                    version: Some("v1".to_owned()),
                    message: "Alpha source-level footnote".to_owned(),
                    matches: MatchLogic {
                        gene: Some("ZZZ".to_owned()),
                        ..MatchLogic::default()
                    },
                },
            ]
            .into_iter()
            .collect(),
            report_variants: BTreeSet::new(),
            urls: vec!["https://example.test/cpic".to_owned()],
            citations: [Publication {
                pmid: Some("12345".to_owned()),
                title: Some("Synthetic citation".to_owned()),
                journal: Some("Journal".to_owned()),
                year: Some(2026),
                same_as: Some("https://example.test/citation".to_owned()),
            }]
            .into_iter()
            .collect(),
            guidelines: [GuidelineReport {
                id: "guideline-cpic".to_owned(),
                name: "CPIC guideline".to_owned(),
                source: PrescribingGuidanceSource::CpicGuideline,
                version: "v-test".to_owned(),
                url: Some("https://example.test/cpic".to_owned()),
                genes: ["GENE1".to_owned()].into_iter().collect(),
                report_genes: vec![cyp2d6, dpyd],
                annotations: [annotation].into_iter().collect(),
            }]
            .into_iter()
            .collect(),
        };
        let mut unmatched_cyp2d6 = ReportGene::new("CYP2D6", ["Normal Metabolizer".to_owned()]);
        unmatched_cyp2d6.recommendation_diplotypes = vec![ReportDiplotype {
            gene: "CYP2D6".to_owned(),
            label: "*36 + *10/*1".to_owned(),
            inferred: true,
            inferred_source_diplotypes: vec![ReportDiplotype {
                gene: "CYP2D6".to_owned(),
                label: "*36 + *10 + *2 + *3 + *4/*1 + *5".to_owned(),
                ..ReportDiplotype::default()
            }],
            ..ReportDiplotype::default()
        }];
        let mut unmatched_dpyd = ReportGene::new("DPYD", ["1.0".to_owned()]);
        unmatched_dpyd.matcher_homozygous_component_haplotypes =
            ["*2A".to_owned()].into_iter().collect();
        unmatched_dpyd.recommendation_diplotypes = vec![ReportDiplotype {
            gene: "DPYD".to_owned(),
            label: "*2A/*1".to_owned(),
            inferred: true,
            inferred_source_diplotypes: vec![ReportDiplotype {
                gene: "DPYD".to_owned(),
                label: "*2A".to_owned(),
                ..ReportDiplotype::default()
            }],
            ..ReportDiplotype::default()
        }];
        let mut outside_no_genotype = ReportGene::new("GENE3", ["Normal Function".to_owned()]);
        outside_no_genotype.recommendation_diplotypes = vec![ReportDiplotype {
            gene: "GENE3".to_owned(),
            label: "Normal Function".to_owned(),
            outside_phenotype: true,
            ..ReportDiplotype::default()
        }];
        let mut no_data_gene = ReportGene::new("NODATA", ["Unknown".to_owned()]);
        no_data_gene.recommendation_diplotypes = vec![ReportDiplotype {
            gene: "NODATA".to_owned(),
            label: "Unknown".to_owned(),
            ..ReportDiplotype::default()
        }];
        let fda_report = DrugReport {
            name: "test drug".to_owned(),
            id: "RxNorm:test".to_owned(),
            source: PrescribingGuidanceSource::FdaLabel,
            version: "v-test".to_owned(),
            messages: BTreeSet::new(),
            report_variants: BTreeSet::new(),
            urls: vec!["https://example.test/fda".to_owned()],
            citations: BTreeSet::new(),
            guidelines: [GuidelineReport {
                id: "guideline-fda".to_owned(),
                name: "FDA label".to_owned(),
                source: PrescribingGuidanceSource::FdaLabel,
                version: "v-test".to_owned(),
                url: Some("https://example.test/fda".to_owned()),
                genes: [
                    "CYP2D6".to_owned(),
                    "DPYD".to_owned(),
                    "GENE3".to_owned(),
                    "NODATA".to_owned(),
                ]
                .into_iter()
                .collect(),
                report_genes: vec![
                    unmatched_cyp2d6,
                    unmatched_dpyd,
                    outside_no_genotype,
                    no_data_gene,
                ],
                annotations: BTreeSet::new(),
            }]
            .into_iter()
            .collect(),
        };
        let mut gene4 = ReportGene::new("GENE4", ["Normal Function".to_owned()]);
        gene4.recommendation_diplotypes = vec![ReportDiplotype {
            gene: "GENE4".to_owned(),
            label: "Normal Function".to_owned(),
            ..ReportDiplotype::default()
        }];
        let fda_assoc_annotation = AnnotationReport {
            local_id: "FDA-test-1".to_owned(),
            drug_recommendation: Some("\"Use the FDA recommendation\"".to_owned()),
            classification: "n/a".to_owned(),
            population: "general".to_owned(),
            genotypes: vec![RecommendationGenotype::from_report_genes([gene4])],
            implications: vec!["fda implication".to_owned()],
            phenotypes: BTreeMap::new(),
            activity_scores: BTreeMap::new(),
            highlighted_variants: BTreeSet::new(),
            dosing_information: false,
            alternate_drug_available: false,
            other_prescribing_guidance: false,
            messages: BTreeSet::new(),
            lookup_key: Vec::new(),
        };
        let fda_assoc_report = DrugReport {
            name: "test drug".to_owned(),
            id: "RxNorm:test".to_owned(),
            source: PrescribingGuidanceSource::FdaAssoc,
            version: "v-test".to_owned(),
            messages: BTreeSet::new(),
            report_variants: BTreeSet::new(),
            urls: vec!["https://example.test/fda-assoc".to_owned()],
            citations: BTreeSet::new(),
            guidelines: [GuidelineReport {
                id: "guideline-fda-assoc".to_owned(),
                name: "FDA association: PGx table entry".to_owned(),
                source: PrescribingGuidanceSource::FdaAssoc,
                version: "v-test".to_owned(),
                url: Some("https://example.test/fda-assoc".to_owned()),
                genes: ["GENE4".to_owned()].into_iter().collect(),
                report_genes: Vec::new(),
                annotations: [fda_assoc_annotation].into_iter().collect(),
            }]
            .into_iter()
            .collect(),
        };
        drug_reports.insert(
            PrescribingGuidanceSource::CpicGuideline,
            [("test drug".to_owned(), cpic_report)]
                .into_iter()
                .collect(),
        );
        drug_reports.insert(
            PrescribingGuidanceSource::FdaLabel,
            [("test drug".to_owned(), fda_report)].into_iter().collect(),
        );
        drug_reports.insert(
            PrescribingGuidanceSource::FdaAssoc,
            [("test drug".to_owned(), fda_assoc_report)]
                .into_iter()
                .collect(),
        );
        let context = ReportContext {
            title: Some("grouping".to_owned()),
            data_version: "v-test".to_owned(),
            gene_reports: BTreeMap::new(),
            report_gene_sources: BTreeMap::new(),
            drug_reports,
            messages: Vec::new(),
            unannotated_gene_calls: Vec::new(),
        };

        let html = report_html_string(
            &context,
            &HtmlReportOptions {
                debug: true,
                no_data_genes: ["NODATA".to_owned()].into_iter().collect(),
                ..HtmlReportOptions::default()
            },
        );

        assert!(html.contains("<section class=\"guideline drugReport test-drug\">"));
        assert!(html.contains("<th>Source</th><th>Genes</th><th>Implications</th><th>Recommendation</th><th>Classification</th>"));
        assert!(html.contains("<tr class=\"top-aligned cpic-guideline-test_drug\">"));
        assert!(html.contains("<tr class=\"top-aligned fda-label-test_drug\">"));
        assert!(html.contains("<tr class=\"top-aligned fda-assoc-test_drug\">"));
        assert_eq!(html.matches("<p>Population:<br/>general</p>").count(), 1);
        assert!(html.contains(
            "<a href=\"https://example.test/fda-assoc\" target=\"_blank\">FDA association</a></b></p><p>general</p>"
        ));
        assert!(!html.contains(
            "<a href=\"https://example.test/fda-assoc\" target=\"_blank\">FDA association</a></b></p><p>Population:<br/>general</p>"
        ));
        assert!(html.contains("Source-level ambiguity message"));
        assert!(html.contains("Source-level note message"));
        assert!(
            html.find("Source-level ambiguity message")
                .expect("ambiguity source message")
                < html
                    .find("Source-level note message")
                    .expect("note source message")
        );
        assert!(html.contains("Use the CPIC recommendation"));
        assert!(html.contains("<div class=\"hint\">Genotypes</div><ul class=\"noPadding mt-0\"><li><span class=\"noWrap\"><span class=\"rx-dip\"><a href=\"#CYP2D6\">CYP2D6</a>:*4/*1</span>"));
        assert!(html.contains("<div class=\"alert alert-debug\"><div class=\"hint\">Inferred:</div><span class=\"nowrap\"><span><a href=\"#CYP2D6\">CYP2D6</a>:*2/*1</span></span></div>"));
        assert!(html.contains("<span class=\"rx-dip\"><a href=\"#DPYD\">DPYD</a>:*1/*2A</span>"));
        assert!(html.contains("<div class=\"hint\">Genotype</div><ul class=\"noBullet mt-0\"><li><span class=\"noWrap\"><span class=\"rx-dip\"><a href=\"#GENE4\">GENE4</a>:Normal Function</span>"));
        assert!(html.contains("<span class=\"rx-hl-var\">rsTest:A/T</span>"));
        assert!(html.contains("<div class=\"tag\">Alternate Drug</div>\n<div class=\"tag\">Dosing Info</div>\n<div class=\"tag\">Other Guidance</div>"));
        assert!(html.contains("<div class=\"tag noAction\">No Action</div>"));
        assert!(html.contains("<span class=\"rx-unmatched-dip\"><a href=\"#CYP2D6\">CYP2D6</a>:*36 + *10 + *2 + *3 + *4/<br />&nbsp;*1 + *5</span>"));
        assert!(html.contains(
            "<span class=\"rx-unmatched-dip\"><a href=\"#DPYD\">DPYD</a>:*2A (homozygous)</span>"
        ));
        assert!(html.contains(
            "<span class=\"rx-unmatched-dip\"><a href=\"#GENE3\">GENE3</a>:Not provided</span>"
        ));
        assert!(html.contains(
            "<span class=\"rx-unmatched-dip\">NODATA:Uncalled - no variant data provided</span>"
        ));
        assert!(html.contains("<tr class=\"top-aligned fda-label-test_drug\"><td class=\"top-aligned\"><p><b>FDA Label Annotation</b><sup class=\"sources\"><a href=\"https://example.test/fda\" target=\"_blank\" rel=\"noopener noreferrer\">1</a></sup></p></td><td class=\"top-aligned\"><div class=\"hint\">Genotypes</div>"));
        assert!(html.contains("</td><td colspan=\"4\" class=\"top-aligned\"><div class=\"hint\">&nbsp;</div>FDA Label Annotation provides no genotype-based recommendations for this genotype, after evaluating the evidence.</td>"));
        assert!(html.contains("<p class=\"noGeneData\">No data provided for <span class=\"gene nodata\"><span class=\"no-data\">NODATA</span></span>.</p>"));
        assert!(html.contains("<b>FDA Label Annotation</b><sup class=\"sources\"><a href=\"https://example.test/fda\" target=\"_blank\" rel=\"noopener noreferrer\">1</a></sup>"));
        assert!(html.contains(
            "<a href=\"https://example.test/fda-assoc\" target=\"_blank\">FDA association</a>"
        ));
        assert!(html.contains("*2/*1"));
        assert!(html.contains("*1/*2A"));
        assert!(html.contains(
            "<sup><a href=\"#rx-dagger-test-drug\" title=\"Inferred\">&dagger;</a></sup>"
        ));
        assert!(html.contains(
            "<sup><a href=\"#rx-ddagger-test-drug\" title=\"Inferred\">&ddagger;</a></sup>"
        ));
        assert!(
            html.contains(
                "<ul class=\"noPadding mt-0\"><li>increased exposure</li><li>N/A</li></ul>"
            )
        );
        assert!(html.contains("<td>fda implication</td><td class=\"drugRecommendation\">\"Use the FDA recommendation\""));
        assert!(html.contains("<div class=\"hint\">Phenotypes</div><dl class=\"compact mt-0\"><div class=\"rx-phenotype rx-phenotype--GENE1\"><dt>GENE1:</dt><dd>intermediate metabolizer</dd></div><div class=\"rx-phenotype rx-phenotype--GENE2\"><dt>GENE2:</dt><dd>N/A</dd></div></dl>"));
        assert!(html.contains(
            "<div class=\"hint\">Activity Score</div><p class=\"rx-activity\">Unspecified</p>"
        ));
        assert!(html.contains("<td class=\"drugRecommendation\">\"Use the FDA recommendation\"<a href=\"#rx-ast-test-drug\" style=\"text-decoration: none\">&ast;</a></td>"));
        assert!(html.contains("<td class=\"drugRecClass\">N/A</td>"));
        assert!(
            !html.contains("<div class=\"recommendation\">\"Use the FDA recommendation\"</div>")
        );
        assert!(html.contains("FDA Label Annotation provides no genotype-based recommendations for this genotype, after evaluating the evidence."));
        assert!(html.contains("<div class=\"footnote\" id=\"rx-dagger-test-drug\"><sup>&dagger;</sup> Inferred genotype used to look up phenotype."));
        assert!(html.contains("<a href=\"https://pharmcat.org/methods/Gene-Definition-Exceptions/\" target=\"_blank\" rel=\"noopener noreferrer\">PharmCAT documentation</a>"));
        assert!(html.contains("<div class=\"footnote\" id=\"rx-ddagger-test-drug\"><sup>&ddagger;</sup> The DPYD genotype used to look up phenotype is inferred from the two lowest function haplotypes."));
        assert!(html.contains("<a href=\"https://pharmcat.org/methods/Gene-Definition-Exceptions/#dpyd\" target=\"_blank\" rel=\"noopener noreferrer\">PharmCAT documentation</a>"));
        assert!(html.contains("<div class=\"footnote\" id=\"rx-ast-test-drug\"><sup>&ast;</sup> Text in quotation is taken directly from the FDA Label or FDA PGx Association table."));
        assert!(html.contains(
            "<div class=\"footnote same-source-footnote\">Alpha source-level footnote</div>"
        ));
        assert!(html.contains(
            "<div class=\"footnote same-source-footnote\">Beta source-level footnote</div>"
        ));
        assert!(!html.contains("alert alert-info same-source-footnote"));
        assert!(
            html.find("Alpha source-level footnote")
                .expect("alpha drug report footnote")
                < html
                    .find("Beta source-level footnote")
                    .expect("beta drug report footnote")
        );
        assert!(html.contains("PMID:12345"));
        assert!(html.contains("<li><a href=\"https://example.test/citation\" target=\"_blank\" rel=\"noopener noreferrer\">Synthetic citation</a>. <i>Journal</i>. 2026. PMID:12345</li>"));
        assert!(
            html.find("<div class=\"footnote\" id=\"rx-dagger-test-drug\"")
                .expect("inferred footnote")
                < html
                    .find("<div class=\"citations\"><p>Citations:</p><ul>")
                    .expect("citations")
        );
        assert!(
            html.find(
                "<div class=\"footnote same-source-footnote\">Alpha source-level footnote</div>"
            )
            .expect("drug report footnote")
                < html
                    .find("<div class=\"citations\"><p>Citations:</p><ul>")
                    .expect("citations")
        );
        assert!(
            html.find("CPIC Guideline Annotation").expect("cpic row")
                < html.find("FDA Label Annotation").expect("fda row")
        );
    }

    #[test]
    fn report_html_string_includes_unmatched_recommendations_only_in_extended_mode_like_java() {
        let mut cyp2c19 = ReportGene::new("CYP2C19", ["Normal Metabolizer".to_owned()]);
        cyp2c19.recommendation_diplotypes = vec![ReportDiplotype {
            gene: "CYP2C19".to_owned(),
            label: "*1/*1".to_owned(),
            ..ReportDiplotype::default()
        }];
        let drug_report = DrugReport {
            name: "unmatched drug".to_owned(),
            id: "RxNorm:unmatched".to_owned(),
            source: PrescribingGuidanceSource::CpicGuideline,
            version: "v-test".to_owned(),
            messages: BTreeSet::new(),
            report_variants: BTreeSet::new(),
            urls: vec!["https://example.test/cpic".to_owned()],
            citations: BTreeSet::new(),
            guidelines: [GuidelineReport {
                id: "guideline-unmatched".to_owned(),
                name: "CPIC unmatched guideline".to_owned(),
                source: PrescribingGuidanceSource::CpicGuideline,
                version: "v-test".to_owned(),
                url: Some("https://example.test/cpic".to_owned()),
                genes: ["CYP2C19".to_owned()].into_iter().collect(),
                report_genes: vec![cyp2c19],
                annotations: BTreeSet::new(),
            }]
            .into_iter()
            .collect(),
        };
        assert!(!drug_report.is_matched());
        assert!(
            drug_report
                .guidelines
                .iter()
                .any(GuidelineReport::is_reportable)
        );
        let not_called_drug_report = DrugReport {
            name: "not called drug".to_owned(),
            id: "RxNorm:not-called".to_owned(),
            source: PrescribingGuidanceSource::CpicGuideline,
            version: "v-test".to_owned(),
            messages: BTreeSet::new(),
            report_variants: BTreeSet::new(),
            urls: vec!["https://example.test/not-called".to_owned()],
            citations: BTreeSet::new(),
            guidelines: [GuidelineReport {
                id: "guideline-not-called".to_owned(),
                name: "CPIC not called guideline".to_owned(),
                source: PrescribingGuidanceSource::CpicGuideline,
                version: "v-test".to_owned(),
                url: Some("https://example.test/not-called".to_owned()),
                genes: ["CYP2D6".to_owned()].into_iter().collect(),
                report_genes: Vec::new(),
                annotations: BTreeSet::new(),
            }]
            .into_iter()
            .collect(),
        };
        assert!(!not_called_drug_report.is_matched());
        assert!(
            !not_called_drug_report
                .guidelines
                .iter()
                .any(GuidelineReport::is_reportable)
        );

        let context = ReportContext {
            title: None,
            data_version: "v-test".to_owned(),
            gene_reports: BTreeMap::new(),
            report_gene_sources: BTreeMap::new(),
            drug_reports: [(
                PrescribingGuidanceSource::CpicGuideline,
                [
                    ("unmatched drug".to_owned(), drug_report),
                    ("not called drug".to_owned(), not_called_drug_report),
                ]
                .into_iter()
                .collect(),
            )]
            .into_iter()
            .collect(),
            messages: Vec::new(),
            unannotated_gene_calls: Vec::new(),
        };
        let compact_options = HtmlReportOptions {
            compact: true,
            ..HtmlReportOptions::default()
        };

        assert!(
            super::html_recommendation_drugs(&context, &HtmlReportOptions::default())
                .contains("unmatched drug")
        );
        assert!(
            super::html_recommendation_drugs(&context, &HtmlReportOptions::default())
                .contains("not called drug")
        );
        assert!(
            super::html_drugs_with_recommendations(&context, &HtmlReportOptions::default())
                .contains("unmatched drug")
        );
        assert!(
            super::html_drugs_with_recommendations(&context, &HtmlReportOptions::default())
                .contains("not called drug")
        );
        assert!(
            !super::html_recommendation_drugs(&context, &compact_options)
                .contains("unmatched drug")
        );
        assert!(
            !super::html_recommendation_drugs(&context, &compact_options)
                .contains("not called drug")
        );
        assert!(
            !super::html_drugs_with_recommendations(&context, &compact_options)
                .contains("unmatched drug")
        );
        assert!(
            !super::html_drugs_with_recommendations(&context, &compact_options)
                .contains("not called drug")
        );
        assert!(
            super::html_drugs_without_recommendations(&context, &compact_options)
                .contains("unmatched drug")
        );
        assert!(
            super::html_drugs_without_recommendations(&context, &HtmlReportOptions::default())
                .is_empty()
        );

        let extended_html = report_html_string(&context, &HtmlReportOptions::default());
        assert!(extended_html.contains("<section class=\"guideline drugReport unmatched-drug\">"));
        assert!(extended_html.contains(
            "CPIC Guideline Annotation provides no genotype-based recommendations for this genotype, after evaluating the evidence."
        ));
        assert!(extended_html.contains("<section class=\"guideline drugReport not-called-drug\">"));
        assert!(extended_html.contains("<tr class=\"top-aligned cpic-guideline-not_called_drug\"><td class=\"top-aligned\"><p><b>CPIC Guideline Annotation</b></p></td><td colspan=\"5\" class=\"top-aligned\"><span class=\"rx-no-call\">No call data for CYP2D6</span>.</td></tr>"));
        assert!(!extended_html.contains(
            "https://example.test/not-called\" target=\"_blank\" rel=\"noopener noreferrer\">1</a>"
        ));
        assert!(!extended_html.contains("Drugs With No Guidance"));

        let compact_html = report_html_string(&context, &compact_options);
        assert!(!compact_html.contains("<section class=\"guideline drugReport unmatched-drug\">"));
        assert!(!compact_html.contains("<section class=\"guideline drugReport not-called-drug\">"));
        assert!(compact_html.contains("<p class=\"rx-no-recs\">No recommendations.</p>"));
        assert!(compact_html.contains("<h3>Drugs With No Guidance</h3>"));
        assert!(compact_html.contains("The following drugs are known to be associated with genes in this report but have no guidance for the specific genotypes in this report."));
        assert!(compact_html.contains(
            "<a href=\"https://www.clinpgx.org/prescribingInfo\">&quot;Prescribing Info&quot; page on ClinPGx</a>"
        ));
        assert!(compact_html.contains("<li>unmatched drug</li>"));
    }

    #[test]
    fn write_report_html_writes_file_and_rejects_non_html_path_like_java() {
        let collection = PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH))
            .expect("prescribing guidance");
        let context = ReportContext::from_gene_reports(
            &collection,
            [ReportGene::new("HLA-B", ["*58:01 positive".to_owned()])
                .with_source_diplotype("*58:01 positive")],
            None,
        );
        let base =
            std::env::temp_dir().join(format!("pharmcat-report-html-{}", std::process::id()));
        fs::create_dir_all(&base).expect("temp dir");
        let html_path = base.join("report.html");
        let txt_path = base.join("report.txt");

        write_report_html(&context, &html_path, &HtmlReportOptions::default())
            .expect("write report HTML");
        let html = fs::read_to_string(&html_path).expect("report HTML contents");
        assert!(html.starts_with("<!DOCTYPE html>"));
        assert!(html.contains("<main>"));

        let error = write_report_html(&context, &txt_path, &HtmlReportOptions::default())
            .expect_err("invalid extension");
        assert!(error.to_string().contains("does not end with .html"));

        fs::remove_file(html_path).ok();
        fs::remove_dir(base).ok();
    }

    #[test]
    fn annotation_report_filters_allele_presence_phenotypes_like_java() {
        let collection = PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH))
            .expect("prescribing guidance");
        let context = ReportContext::from_gene_reports(
            &collection,
            [ReportGene::new("HLA-B", ["*58:01 positive".to_owned()])],
            None,
        );

        let allopurinol = context
            .drug_report(PrescribingGuidanceSource::CpicGuideline, "allopurinol")
            .expect("allopurinol CPIC report");

        assert!(allopurinol.is_matched());

        let guideline = allopurinol.guidelines.iter().next().expect("guideline");
        assert_eq!(
            guideline.name,
            "Annotation of CPIC Guideline for allopurinol and HLA-B"
        );

        let annotation = guideline.annotations.iter().next().expect("annotation");
        assert_eq!(annotation.local_id, "CPIC-PA166296961");
        assert_eq!(
            annotation.phenotypes.get("HLA-B").map(String::as_str),
            Some("*58:01 positive")
        );
    }

    #[test]
    fn report_gene_from_annotated_outside_call_feeds_recommendation_matching_like_java() {
        let validation = OutsideCallValidation::for_supported_genes(["CYP2D6"]);
        let call = parse_outside_call_line(
            &validation,
            "CYP2D6\t*1/*3\tIntermediate Metabolizer\t1.0",
            1,
        )
        .expect("outside call");
        let phenotype_map = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let cyp2d6 = phenotype_map.phenotype("CYP2D6").expect("CYP2D6");

        let report_gene =
            ReportGene::from_outside_call(&call, Some(cyp2d6)).expect("outside report gene");
        let genotype = RecommendationGenotype::from_report_genes([report_gene.clone()]);

        assert_eq!(report_gene.gene, "CYP2D6");
        assert_eq!(report_gene.lookup_keys, ["1.0"]);
        assert_eq!(report_gene.phenotypes, ["Intermediate Metabolizer"]);
        assert_eq!(report_gene.activity_score.as_deref(), Some("1.0"));
        assert_eq!(report_gene.source_diplotype.as_deref(), Some("*1/*3"));
        assert_eq!(report_gene.source_diplotypes.len(), 1);
        assert_eq!(report_gene.source_diplotypes[0].label, "*1/*3");
        assert_eq!(report_gene.recommendation_diplotypes.len(), 1);
        assert_eq!(report_gene.recommendation_diplotypes[0].label, "*1/*3");
        let source_allele1 = report_gene.source_diplotypes[0]
            .allele1
            .as_ref()
            .expect("source allele1");
        let source_allele2 = report_gene.source_diplotypes[0]
            .allele2
            .as_ref()
            .expect("source allele2");
        assert_eq!(source_allele1.name, "*1");
        assert_eq!(source_allele1.function, "Normal function");
        assert_eq!(source_allele1.activity_value.as_deref(), Some("1.0"));
        assert_eq!(source_allele2.name, "*3");
        assert_eq!(source_allele2.function, "No function");
        assert_eq!(source_allele2.activity_value.as_deref(), Some("0.0"));
        assert_eq!(
            report_gene.recommendation_diplotypes[0]
                .allele2
                .as_ref()
                .expect("recommendation allele2")
                .function,
            "No function"
        );
        assert!(report_gene.source_diplotypes[0].outside_phenotype);
        assert!(report_gene.source_diplotypes[0].outside_activity_score);
        assert!(report_gene.recommendation_diplotypes[0].outside_phenotype);
        assert!(report_gene.recommendation_diplotypes[0].outside_activity_score);
        assert_eq!(report_gene.diplotype_key["*1"], Value::from(1.0));
        assert_eq!(report_gene.diplotype_key["*3"], Value::from(1.0));
        assert!(report_gene.is_activity_score_type);
        assert!(report_gene.outside_call);
        assert_eq!(genotype.lookup_keys()[0]["CYP2D6"], Value::from("1.0"));
    }

    #[test]
    fn report_gene_from_dpyd_gene_call_result_feeds_recommendation_matching_like_java() {
        let definition =
            read_definition_file(Path::new(DPYD_DEFINITION_PATH)).expect("DPYD definition");
        let phenotype =
            read_gene_phenotype_file(Path::new(DPYD_PHENOTYPE_PATH)).expect("DPYD phenotype");
        let mut allele_map = reference_allele_map(&definition);
        let variant = variant_by_rsid(&definition, "rs67376798");
        allele_map.insert(
            variant.vcf_chr_position(),
            sample_call(
                &variant.chromosome,
                variant.position as usize,
                Some(&variant.reference),
                Some("A"),
                true,
                true,
            ),
        );

        let result =
            call_dpyd_lowest_function_gene("Sample_1", &definition, &allele_map).expect("DPYD");
        let GeneCallKind::Diplotypes(diplotypes) = &result.kind else {
            panic!("expected DPYD diplotype call");
        };
        assert_eq!(diplotypes[0].name, "Reference/c.2846A>T");

        let report_gene = ReportGene::from_gene_call_result_with_definition(
            &result,
            Some(&phenotype),
            &definition,
        )
        .expect("report gene")
        .expect("DPYD report gene");
        let genotype = RecommendationGenotype::from_report_genes([report_gene.clone()]);
        let report_variant = report_gene
            .variant_reports
            .iter()
            .find(|report| report.db_snp_id.as_deref() == Some("rs67376798"))
            .expect("DPYD variant report");
        let reference_message = report_gene
            .messages
            .iter()
            .find(|message| message.name == "reference-allele")
            .expect("reference allele message");
        let expected_variant_call = format!("{}|A", variant.reference);

        assert_eq!(report_gene.gene, "DPYD");
        assert_eq!(report_gene.chromosome.as_deref(), Some("chr1"));
        assert_eq!(
            report_gene.allele_definition_version,
            definition
                .version
                .clone()
                .or_else(|| definition.data_version.clone())
        );
        assert_eq!(
            report_gene.allele_definition_source,
            data_source_from_definition(&definition)
        );
        assert_eq!(report_gene.phenotype_version, phenotype.version);
        assert_eq!(report_gene.phased, result.match_data.phased);
        assert_eq!(
            report_gene.effectively_phased,
            result.match_data.effectively_phased
        );
        assert_eq!(report_gene.variant_reports.len(), definition.variants.len());
        assert_eq!(report_variant.gene.as_deref(), Some("DPYD"));
        assert_eq!(
            report_variant.chromosome.as_deref(),
            Some(variant.chromosome.as_str())
        );
        assert_eq!(report_variant.position, Some(variant.position as i64));
        assert_eq!(
            report_variant.call.as_deref(),
            Some(expected_variant_call.as_str())
        );
        assert_eq!(
            report_variant.reference_allele.as_deref(),
            Some(variant.reference.as_str())
        );
        assert!(report_variant.phased);
        assert_eq!(
            reference_message.exception_type,
            MessageAnnotation::TYPE_NOTE
        );
        assert_eq!(
            reference_message.message,
            "The DPYD Reference allele assignment is characterized by the absence of variants at the positions that are included in the underlying allele definitions."
        );
        assert_eq!(report_gene.lookup_keys, ["1.5"]);
        assert_eq!(report_gene.phenotypes, ["Intermediate Metabolizer"]);
        assert_eq!(report_gene.activity_score.as_deref(), Some("1.5"));
        assert_eq!(
            report_gene.source_diplotype.as_deref(),
            Some("c.2846A>T (heterozygous)")
        );
        assert_eq!(report_gene.source_diplotypes.len(), 1);
        assert_eq!(
            report_gene.source_diplotypes[0].label,
            "c.2846A>T (heterozygous)"
        );
        let source_allele1 = report_gene.source_diplotypes[0]
            .allele1
            .as_ref()
            .expect("source allele1");
        let source_allele2 = report_gene.source_diplotypes[0]
            .allele2
            .as_ref()
            .expect("source allele2");
        assert_eq!(source_allele1.name, "Reference");
        assert_eq!(source_allele1.function, "Normal function");
        assert!(source_allele1.reference);
        assert_eq!(source_allele1.activity_value.as_deref(), Some("1.0"));
        assert_eq!(source_allele2.name, "c.2846A>T");
        assert_eq!(source_allele2.function, "Decreased function");
        assert!(!source_allele2.reference);
        assert_eq!(source_allele2.activity_value.as_deref(), Some("0.5"));
        // Java DiplotypeMatch score sums the two NamedAllele scores from the definition JSON:
        // Reference (83) + c.2846A>T (1) = 84. (The earlier 166 was a wrong 2xReference guess.)
        assert_eq!(
            report_gene.source_diplotypes[0].match_score.as_deref(),
            Some("84")
        );
        assert!(!report_gene.source_diplotypes[0].inferred);
        assert!(
            report_gene.source_diplotypes[0]
                .inferred_source_diplotypes
                .is_empty()
        );
        assert!(!report_gene.source_diplotypes[0].combination);
        assert_eq!(
            report_gene.matcher_component_haplotypes,
            ["Reference".to_owned(), "c.2846A>T".to_owned()]
                .into_iter()
                .collect()
        );
        assert_eq!(report_gene.matcher_component_diplotypes.len(), 2);
        assert_eq!(
            report_gene.matcher_component_diplotypes[0].label,
            "Reference"
        );
        assert_eq!(
            report_gene.matcher_component_diplotypes[1].label,
            "c.2846A>T"
        );
        assert_eq!(
            report_gene.matcher_component_diplotypes[1]
                .allele1
                .as_ref()
                .expect("component allele1")
                .activity_value
                .as_deref(),
            Some("0.5")
        );
        assert_eq!(report_gene.recommendation_diplotypes.len(), 1);
        assert_eq!(
            report_gene.recommendation_diplotypes[0].label,
            "c.2846A>T (heterozygous)"
        );
        assert!(!report_gene.recommendation_diplotypes[0].inferred);
        assert_eq!(
            report_gene.recommendation_diplotypes[0]
                .inferred_source_diplotypes
                .len(),
            1
        );
        assert_eq!(
            report_gene.recommendation_diplotypes[0].inferred_source_diplotypes[0].label,
            "c.2846A>T (heterozygous)"
        );
        assert!(!report_gene.recommendation_diplotypes[0].combination);
        assert_eq!(
            report_gene.recommendation_diplotypes[0]
                .activity_score
                .as_deref(),
            Some("1.5")
        );
        assert_eq!(
            report_gene.recommendation_diplotypes[0]
                .allele2
                .as_ref()
                .expect("recommendation allele2")
                .display_name(),
            "DPYD c.2846A>T"
        );
        assert_eq!(report_gene.diplotype_key["Reference"], Value::from(1.0));
        assert_eq!(report_gene.diplotype_key["c.2846A>T"], Value::from(1.0));
        assert!(report_gene.is_activity_score_type);
        assert_eq!(genotype.lookup_keys()[0]["DPYD"], Value::from("1.5"));
    }

    #[test]
    fn report_gene_from_dpyd_no_call_result_emits_unknown_diplotypes_like_java() {
        let definition =
            read_definition_file(Path::new(DPYD_DEFINITION_PATH)).expect("DPYD definition");
        let phenotype =
            read_gene_phenotype_file(Path::new(DPYD_PHENOTYPE_PATH)).expect("DPYD phenotype");
        let allele_map = BTreeMap::new();

        let result =
            call_dpyd_lowest_function_gene("Sample_1", &definition, &allele_map).expect("DPYD");
        assert!(matches!(result.kind, GeneCallKind::NoCall));

        let report_gene = ReportGene::from_dpyd_gene_call_result(&result, &phenotype)
            .expect("report gene")
            .expect("DPYD no-call report gene");
        let genotype = RecommendationGenotype::from_report_genes([report_gene.clone()]);

        assert_eq!(report_gene.gene, "DPYD");
        assert_eq!(report_gene.lookup_keys, ["No Result"]);
        assert_eq!(report_gene.phenotypes, ["No Result"]);
        assert_eq!(report_gene.activity_score.as_deref(), Some("No Result"));
        assert_eq!(
            report_gene.source_diplotype.as_deref(),
            Some("Unknown/Unknown")
        );
        assert_eq!(report_gene.source_diplotypes.len(), 1);
        assert_eq!(report_gene.recommendation_diplotypes.len(), 1);
        assert_eq!(report_gene.source_diplotypes[0].label, "Unknown/Unknown");
        assert_eq!(
            report_gene.recommendation_diplotypes[0].label,
            "Unknown/Unknown"
        );
        assert_eq!(
            report_gene.source_diplotypes[0]
                .diplotype_key
                .get("Unknown"),
            Some(&Value::from(2.0))
        );
        assert_eq!(
            report_gene.recommendation_diplotypes[0]
                .diplotype_key
                .get("Unknown"),
            Some(&Value::from(2.0))
        );
        assert_eq!(genotype.lookup_keys()[0]["DPYD"], Value::from("No Result"));
    }

    #[test]
    fn report_gene_from_standard_no_call_result_emits_unknown_diplotypes_like_java() {
        let definition =
            read_definition_file(Path::new(CYP3A5_DEFINITION_PATH)).expect("CYP3A5 definition");
        let phenotype =
            read_gene_phenotype_file(Path::new(CYP3A5_PHENOTYPE_PATH)).expect("CYP3A5 phenotype");
        let allele_map = BTreeMap::new();

        let result = call_standard_gene("Sample_1", &definition, &allele_map, false, false)
            .expect("standard gene call");
        assert_eq!(result.gene, "CYP3A5");
        assert!(matches!(result.kind, GeneCallKind::NoCall));

        let report_gene = ReportGene::from_standard_gene_call_result(&result, Some(&phenotype))
            .expect("report gene")
            .expect("CYP3A5 no-call report gene");
        let genotype = RecommendationGenotype::from_report_genes([report_gene.clone()]);

        assert_eq!(report_gene.lookup_keys, ["No Result"]);
        assert_eq!(report_gene.phenotypes, ["No Result"]);
        assert_eq!(report_gene.activity_score, None);
        assert_eq!(
            report_gene.source_diplotype.as_deref(),
            Some("Unknown/Unknown")
        );
        assert_eq!(report_gene.source_diplotypes.len(), 1);
        assert_eq!(report_gene.recommendation_diplotypes.len(), 1);
        assert_eq!(report_gene.source_diplotypes[0].label, "Unknown/Unknown");
        assert_eq!(
            report_gene.recommendation_diplotypes[0].label,
            "Unknown/Unknown"
        );
        assert_eq!(
            genotype.lookup_keys()[0]["CYP3A5"],
            Value::from("No Result")
        );
    }

    #[test]
    fn report_gene_from_standard_no_call_carries_gene_call_warnings_like_java() {
        let definition =
            read_definition_file(Path::new(CYP3A5_DEFINITION_PATH)).expect("CYP3A5 definition");
        let phenotype =
            read_gene_phenotype_file(Path::new(CYP3A5_PHENOTYPE_PATH)).expect("CYP3A5 phenotype");
        let allele_map = BTreeMap::new();
        let mut result = call_standard_gene("Sample_1", &definition, &allele_map, false, false)
            .expect("standard gene call");
        assert!(matches!(result.kind, GeneCallKind::NoCall));
        result
            .warnings
            .insert(GeneCallWarning::MissingRequiredPosition(vec![
                "chr7:99672916".to_owned(),
            ]));

        let report_gene = ReportGene::from_standard_gene_call_result(&result, Some(&phenotype))
            .expect("report gene")
            .expect("CYP3A5 no-call report gene");
        let message = report_gene
            .messages
            .iter()
            .next()
            .expect("gene warning message");

        assert_eq!(message.name, "missing-required-position");
        assert_eq!(message.exception_type, MessageAnnotation::TYPE_NOTE);
        assert_eq!(
            message.message,
            "CYP3A5 - missing required variant (chr7:99672916)"
        );
    }

    #[test]
    fn report_gene_from_standard_diplotype_carries_gene_call_warnings_like_java() {
        let definition =
            read_definition_file(Path::new(CYP3A5_DEFINITION_PATH)).expect("CYP3A5 definition");
        let phenotype =
            read_gene_phenotype_file(Path::new(CYP3A5_PHENOTYPE_PATH)).expect("CYP3A5 phenotype");
        let records =
            read_record_summaries(Path::new(CYP3A5_VCF_PATH), Some("NA12878")).expect("VCF");
        let allele_map = allele_map_from_record_summaries(
            records
                .records
                .into_iter()
                .filter_map(|record| record.allele_call),
        );

        let mut result = call_standard_gene("NA12878", &definition, &allele_map, false, false)
            .expect("standard gene call");
        result.warnings.insert(GeneCallWarning::UnphasedPriority);
        result
            .warnings
            .insert(GeneCallWarning::MissingAmp1Position(vec![
                "chr7:99672916".to_owned(),
                "chr7:99672917".to_owned(),
            ]));

        let report_gene = ReportGene::from_standard_gene_call_result(&result, Some(&phenotype))
            .expect("report gene")
            .expect("CYP3A5 report gene");
        let messages = report_gene
            .messages
            .iter()
            .map(|message| (message.name.as_str(), message.message.as_str()))
            .collect::<Vec<_>>();

        assert_eq!(report_gene.messages.len(), 2);
        assert!(messages.iter().any(|(name, message)| {
            *name == "unphased-priority"
                && *message
                    == "Unphased CYP3A5 variants resulted in multiple calls.  PharmCAT is picking a single call based on frequency data.  Please consult the documentation for details."
        }));
        assert!(messages.iter().any(|(name, message)| {
            *name == "missing-amp1-position"
                && *message
                    == "Missing variants required to meet AMP Tier 1 requirements:  chr7:99672916, chr7:99672917. See https://www.clinpgx.org/ampAllelesToTest for details."
        }));
    }

    #[test]
    fn slco1b1_custom_fallback_infers_recommendation_from_rs4149056_like_java() {
        let phenotype =
            read_gene_phenotype_file(Path::new(SLCO1B1_PHENOTYPE_PATH)).expect("SLCO1B1 phenotype");

        let mut report_gene = ReportGene::unknown("SLCO1B1", Some(&phenotype))
            .expect("SLCO1B1 unknown report")
            .with_variant_reports([VariantReport::new("rs4149056", Some("T/C"))]);
        report_gene
            .apply_slco1b1_custom_recommendation(Some(&phenotype))
            .expect("SLCO1B1 custom caller");
        let genotype = RecommendationGenotype::from_report_genes([report_gene.clone()]);

        assert_eq!(
            report_gene.source_diplotype.as_deref(),
            Some("Unknown/Unknown")
        );
        assert_eq!(report_gene.source_diplotypes.len(), 1);
        assert_eq!(report_gene.source_diplotypes[0].label, "Unknown/Unknown");
        assert_eq!(
            report_gene.source_diplotypes[0]
                .variant
                .as_ref()
                .and_then(|variant| variant.db_snp_id.as_deref()),
            Some("rs4149056")
        );
        assert_eq!(report_gene.recommendation_diplotypes.len(), 1);
        let recommendation = &report_gene.recommendation_diplotypes[0];
        assert_eq!(recommendation.label, "*1/*5");
        assert!(recommendation.inferred);
        assert_eq!(recommendation.match_score.as_deref(), Some("0"));
        assert_eq!(
            recommendation
                .variant
                .as_ref()
                .and_then(|variant| variant.call.as_deref()),
            Some("T/C")
        );
        assert_eq!(recommendation.inferred_source_diplotypes.len(), 1);
        assert_eq!(
            recommendation.inferred_source_diplotypes[0].label,
            "rs4149056 C/rs4149056 T"
        );
        assert_eq!(report_gene.lookup_keys, recommendation.lookup_keys);
        assert_eq!(report_gene.phenotypes, recommendation.phenotypes);
        assert_eq!(report_gene.diplotype_key["*1"], Value::from(1.0));
        assert_eq!(report_gene.diplotype_key["*5"], Value::from(1.0));
        assert_eq!(
            genotype.lookup_keys()[0]["SLCO1B1"],
            Value::from("Decreased Function")
        );
    }

    #[test]
    fn slco1b1_custom_fallback_rejects_duplicate_rs4149056_reports_like_java() {
        let phenotype =
            read_gene_phenotype_file(Path::new(SLCO1B1_PHENOTYPE_PATH)).expect("SLCO1B1 phenotype");
        let mut report_gene = ReportGene::unknown("SLCO1B1", Some(&phenotype))
            .expect("SLCO1B1 unknown report")
            .with_variant_reports([
                VariantReport::new("rs4149056", Some("T/C")),
                VariantReport::new("rs4149056", Some("T/T")),
            ]);

        let error = report_gene
            .apply_slco1b1_custom_recommendation(Some(&phenotype))
            .expect_err("duplicate rs4149056 reports should fail");

        assert_eq!(error, Slco1b1CustomCallError::MultipleRs4149056Reports);
    }

    #[test]
    fn slco1b1_custom_fallback_ignores_unsupported_rs4149056_call_like_java() {
        let phenotype =
            read_gene_phenotype_file(Path::new(SLCO1B1_PHENOTYPE_PATH)).expect("SLCO1B1 phenotype");
        let mut report_gene = ReportGene::unknown("SLCO1B1", Some(&phenotype))
            .expect("SLCO1B1 unknown report")
            .with_variant_reports([VariantReport::new("rs4149056", Some("A/C"))]);

        report_gene
            .apply_slco1b1_custom_recommendation(Some(&phenotype))
            .expect("unsupported rs4149056 call should be ignored");

        assert_eq!(report_gene.source_diplotypes[0].label, "Unknown/Unknown");
        assert_eq!(
            report_gene.recommendation_diplotypes[0].label,
            "Unknown/Unknown"
        );
        assert_eq!(report_gene.lookup_keys, ["No Result"]);
        assert!(!report_gene.recommendation_diplotypes[0].inferred);
    }

    #[test]
    fn standard_report_with_definition_invokes_slco1b1_fallback_from_matcher_variants() {
        let definition =
            read_definition_file(Path::new(SLCO1B1_DEFINITION_PATH)).expect("SLCO1B1 definition");
        let phenotype =
            read_gene_phenotype_file(Path::new(SLCO1B1_PHENOTYPE_PATH)).expect("SLCO1B1 phenotype");
        let rs2306283 = variant_by_rsid(&definition, "rs2306283");
        let rs4149056 = variant_by_rsid(&definition, "rs4149056");
        let rs11045853 = variant_by_rsid(&definition, "rs11045853");
        let rs72559748 = variant_by_rsid(&definition, "rs72559748");
        let allele_map = allele_map_from_record_summaries([
            sample_call(
                &rs2306283.chromosome,
                rs2306283.position as usize,
                Some("G"),
                Some("G"),
                false,
                true,
            ),
            sample_call(
                &rs4149056.chromosome,
                rs4149056.position as usize,
                Some("T"),
                Some("C"),
                false,
                true,
            ),
            sample_call(
                &rs11045853.chromosome,
                rs11045853.position as usize,
                Some("A"),
                Some("A"),
                false,
                true,
            ),
            sample_call(
                &rs72559748.chromosome,
                rs72559748.position as usize,
                Some("G"),
                Some("G"),
                false,
                true,
            ),
        ]);

        let result = GeneCallResult {
            gene: "SLCO1B1".to_owned(),
            match_data: MatchData::new("Sample_1", "SLCO1B1", &definition, &allele_map),
            kind: GeneCallKind::NoCall,
            dpyd_hap_b3_warnings: BTreeSet::new(),
            warnings: BTreeSet::new(),
        };

        let report_gene = ReportGene::from_standard_gene_call_result_with_definition(
            &result,
            Some(&phenotype),
            &definition,
        )
        .expect("report gene")
        .expect("SLCO1B1 report gene");
        let genotype = RecommendationGenotype::from_report_genes([report_gene.clone()]);
        let rs4149056_report = report_gene
            .variant_reports
            .iter()
            .find(|variant| variant.db_snp_id.as_deref() == Some("rs4149056"))
            .expect("rs4149056 report");

        assert_eq!(report_gene.source_diplotypes[0].label, "Unknown/Unknown");
        assert_eq!(report_gene.recommendation_diplotypes[0].label, "*1/*5");
        assert!(report_gene.recommendation_diplotypes[0].inferred);
        assert_eq!(rs4149056_report.call.as_deref(), Some("T/C"));
        assert_eq!(rs4149056_report.reference_allele.as_deref(), Some("T"));
        assert_eq!(rs4149056_report.chromosome.as_deref(), Some("chr12"));
        assert_eq!(rs4149056_report.position, Some(21178615));
        assert!(!rs4149056_report.phased);
        assert!(rs4149056_report.alleles.contains(&"*5".to_owned()));
        assert!(
            report_gene
                .variant_reports
                .iter()
                .any(VariantReport::is_missing)
        );
        assert_eq!(
            genotype.lookup_keys()[0]["SLCO1B1"],
            Value::from("Decreased Function")
        );
    }

    #[test]
    fn report_gene_adds_vcf_warnings_to_matching_variant_reports_like_java() {
        let definition =
            read_definition_file(Path::new(SLCO1B1_DEFINITION_PATH)).expect("SLCO1B1 definition");
        let rs4149056 = variant_by_rsid(&definition, "rs4149056");
        let allele_map = allele_map_from_record_summaries([sample_call(
            &rs4149056.chromosome,
            rs4149056.position as usize,
            Some("T"),
            Some("C"),
            false,
            true,
        )]);
        let result = GeneCallResult {
            gene: "SLCO1B1".to_owned(),
            match_data: MatchData::new("Sample_1", "SLCO1B1", &definition, &allele_map),
            kind: GeneCallKind::NoCall,
            dpyd_hap_b3_warnings: BTreeSet::new(),
            warnings: BTreeSet::new(),
        };
        let mut report_gene = ReportGene::unknown("SLCO1B1", None).expect("unknown report gene");
        report_gene.attach_matcher_variant_reports(&result, &definition);

        report_gene.add_variant_warning_messages(&BTreeMap::from([(
            "chr12:21178615".to_owned(),
            BTreeSet::from(["Ignoring: GT and AD disagree".to_owned()]),
        )]));

        let variant = report_gene
            .variant_reports
            .iter()
            .find(|variant| variant.db_snp_id.as_deref() == Some("rs4149056"))
            .expect("rs4149056 report");
        assert_eq!(
            variant.warnings,
            ["Ignoring: GT and AD disagree".to_owned()]
        );
        assert_eq!(variant.chr_position().as_deref(), Some("chr12:21178615"));
    }

    #[test]
    fn standard_report_marks_undocumented_variation_flags_like_java() {
        let definition =
            read_definition_file(Path::new(SLCO1B1_DEFINITION_PATH)).expect("SLCO1B1 definition");
        let rs2306283 = variant_by_rsid(&definition, "rs2306283");
        let rs4149056 = variant_by_rsid(&definition, "rs4149056");
        let mut rs2306283_call = sample_call(
            &rs2306283.chromosome,
            rs2306283.position as usize,
            Some(&rs2306283.reference),
            Some("G"),
            false,
            true,
        );
        rs2306283_call.vcf_alleles = vec![rs2306283.reference.clone(), "G".to_owned()];
        rs2306283_call.undocumented_variations = BTreeSet::from(["G".to_owned()]);
        rs2306283_call.treat_undocumented_variations_as_reference = true;
        let allele_map = allele_map_from_record_summaries([
            rs2306283_call,
            sample_call(
                &rs4149056.chromosome,
                rs4149056.position as usize,
                Some(&rs4149056.reference),
                Some(&rs4149056.reference),
                false,
                true,
            ),
        ]);
        let result = GeneCallResult {
            gene: "SLCO1B1".to_owned(),
            match_data: MatchData::new("Sample_1", "SLCO1B1", &definition, &allele_map),
            kind: GeneCallKind::NoCall,
            dpyd_hap_b3_warnings: BTreeSet::new(),
            warnings: BTreeSet::new(),
        };

        let report_gene =
            ReportGene::from_standard_gene_call_result_with_definition(&result, None, &definition)
                .expect("report gene")
                .expect("SLCO1B1 report gene");
        let undocumented_report = report_gene
            .variant_reports
            .iter()
            .find(|variant| variant.db_snp_id.as_deref() == Some("rs2306283"))
            .expect("rs2306283 report");
        let documented_report = report_gene
            .variant_reports
            .iter()
            .find(|variant| variant.db_snp_id.as_deref() == Some("rs4149056"))
            .expect("rs4149056 report");

        assert!(result.match_data.treat_undocumented_variations_as_reference);
        assert_eq!(
            result.match_data.positions_with_undocumented_variations,
            BTreeSet::from([rs2306283.clone()])
        );
        assert!(report_gene.has_undocumented_variations);
        assert!(report_gene.treat_undocumented_variations_as_reference);
        assert!(undocumented_report.has_undocumented_variations);
        assert!(!documented_report.has_undocumented_variations);
        assert_eq!(undocumented_report.call.as_deref(), Some("A/G"));
    }

    #[test]
    fn variant_report_sorting_matches_java_compare_to() {
        let mut variants = vec![
            VariantReport {
                chromosome: Some("chr1".to_owned()),
                position: Some(10),
                call: None,
                ..VariantReport::new("missing", None::<String>)
            },
            VariantReport {
                chromosome: Some("chr2".to_owned()),
                position: Some(1),
                call: Some("A/T".to_owned()),
                ..VariantReport::new("chr2", Some("A/T"))
            },
            VariantReport {
                chromosome: Some("chr1".to_owned()),
                position: Some(20),
                call: Some("A/T".to_owned()),
                ..VariantReport::new("chr1-20", Some("A/T"))
            },
            VariantReport {
                chromosome: Some("chr1".to_owned()),
                position: Some(5),
                call: Some("A/T".to_owned()),
                ..VariantReport::new("chr1-5", Some("A/T"))
            },
        ];

        variants.sort();

        assert_eq!(
            variants
                .iter()
                .map(|variant| variant.db_snp_id.as_deref().expect("rsid"))
                .collect::<Vec<_>>(),
            ["chr1-5", "chr1-20", "chr2", "missing"]
        );
    }

    #[test]
    fn variant_report_call_helpers_match_java_variant_utils() {
        assert!(VariantReport::new("het", Some("A/T")).is_het_call());
        assert!(VariantReport::new("het", Some("A|T")).is_het_call());
        assert!(!VariantReport::new("hom", Some("A/A")).is_het_call());
        assert!(!VariantReport::new("haploid", Some("A")).is_het_call());
        assert!(!VariantReport::new("partial", Some("./T")).is_het_call());
        assert!(VariantReport::new("blank", Some(" ")).is_missing());
    }

    #[test]
    fn report_gene_from_standard_diplotype_result_populates_java_diplotypes() {
        let definition =
            read_definition_file(Path::new(CYP3A5_DEFINITION_PATH)).expect("CYP3A5 definition");
        let phenotype =
            read_gene_phenotype_file(Path::new(CYP3A5_PHENOTYPE_PATH)).expect("CYP3A5 phenotype");
        let records =
            read_record_summaries(Path::new(CYP3A5_VCF_PATH), Some("NA12878")).expect("VCF");
        let allele_map = allele_map_from_record_summaries(
            records
                .records
                .into_iter()
                .filter_map(|record| record.allele_call),
        );

        let result = call_standard_gene("NA12878", &definition, &allele_map, false, false)
            .expect("standard gene call");
        let GeneCallKind::Diplotypes(diplotypes) = &result.kind else {
            panic!("expected CYP3A5 diplotype call");
        };
        assert_eq!(diplotypes[0].name, "*1/*2");

        let report_gene = ReportGene::from_standard_gene_call_result(&result, Some(&phenotype))
            .expect("report gene")
            .expect("CYP3A5 report gene");
        let genotype = RecommendationGenotype::from_report_genes([report_gene.clone()]);

        assert_eq!(report_gene.gene, "CYP3A5");
        assert_eq!(report_gene.source_diplotype.as_deref(), Some("*1/*2"));
        assert_eq!(report_gene.source_diplotypes.len(), 1);
        assert_eq!(report_gene.recommendation_diplotypes.len(), 1);
        assert_eq!(report_gene.source_diplotypes[0].label, "*1/*2");
        assert_eq!(report_gene.recommendation_diplotypes[0].label, "*1/*2");
        assert_eq!(
            report_gene.source_diplotypes[0].match_score.as_deref(),
            Some(diplotypes[0].score.to_string()).as_deref()
        );
        let allele1 = report_gene.source_diplotypes[0]
            .allele1
            .as_ref()
            .expect("allele1");
        let allele2 = report_gene.source_diplotypes[0]
            .allele2
            .as_ref()
            .expect("allele2");
        assert_eq!(allele1.name, "*1");
        assert_eq!(allele1.function, "Normal function");
        assert_eq!(allele2.name, "*2");
        assert_eq!(allele2.function, "Unassigned function");
        assert_eq!(report_gene.diplotype_key["*1"], Value::from(1.0));
        assert_eq!(report_gene.diplotype_key["*2"], Value::from(1.0));
        assert_eq!(genotype.lookup_keys()[0]["CYP3A5"], Value::from("n/a"));
    }

    #[test]
    fn report_gene_from_standard_with_definition_sets_matcher_metadata_and_reference_note_like_java()
     {
        let definition =
            read_definition_file(Path::new(CYP3A5_DEFINITION_PATH)).expect("CYP3A5 definition");
        let phenotype =
            read_gene_phenotype_file(Path::new(CYP3A5_PHENOTYPE_PATH)).expect("CYP3A5 phenotype");
        let records =
            read_record_summaries(Path::new(CYP3A5_VCF_PATH), Some("NA12878")).expect("VCF");
        let allele_map = allele_map_from_record_summaries(
            records
                .records
                .into_iter()
                .filter_map(|record| record.allele_call),
        );
        let result = call_standard_gene("NA12878", &definition, &allele_map, false, false)
            .expect("standard gene call");

        let report_gene = ReportGene::from_standard_gene_call_result_with_definition(
            &result,
            Some(&phenotype),
            &definition,
        )
        .expect("report gene")
        .expect("CYP3A5 report gene");
        let reference_message = report_gene
            .messages
            .iter()
            .find(|message| message.name == "reference-allele")
            .expect("reference allele message");

        assert_eq!(report_gene.chromosome.as_deref(), Some("chr1"));
        assert_eq!(
            report_gene.allele_definition_version,
            definition.version.clone()
        );
        assert_eq!(
            report_gene.allele_definition_source,
            super::data_source_from_definition(&definition)
        );
        assert_eq!(report_gene.phenotype_version, phenotype.version.clone());
        assert_eq!(report_gene.phased, result.match_data.phased);
        assert_eq!(
            report_gene.effectively_phased,
            result.match_data.effectively_phased
        );
        assert_eq!(report_gene.source_diplotypes[0].label, "*1/*2");
        assert!(
            report_gene.source_diplotypes[0]
                .allele1
                .as_ref()
                .is_some_and(|haplotype| haplotype.reference)
        );
        assert!(
            !report_gene.source_diplotypes[0]
                .allele2
                .as_ref()
                .is_some_and(|haplotype| haplotype.reference)
        );
        assert_eq!(
            reference_message.exception_type,
            MessageAnnotation::TYPE_NOTE
        );
        assert_eq!(
            reference_message.message,
            "The CYP3A5 *1 allele assignment is characterized by the absence of variants at the positions that are included in the underlying allele definitions."
        );
    }

    #[test]
    fn report_gene_from_standard_combination_result_preserves_all_diplotypes_like_java() {
        let definition = read_definition_file(Path::new(UGT1A1_COMBINATION_DEFINITION_PATH))
            .expect("UGT1A1 definition");
        let phenotype =
            read_gene_phenotype_file(Path::new(UGT1A1_PHENOTYPE_PATH)).expect("UGT1A1 phenotype");
        let records = read_record_summaries(
            Path::new(UGT1A1_PARTIAL_WITH_COMBINATION_VCF_PATH),
            Some("PharmCAT"),
        )
        .expect("VCF");
        let allele_map = allele_map_from_record_summaries(
            records
                .records
                .into_iter()
                .filter_map(|record| record.allele_call),
        );

        let result = call_standard_gene("PharmCAT", &definition, &allele_map, false, true)
            .expect("standard gene call");
        let GeneCallKind::Diplotypes(diplotypes) = &result.kind else {
            panic!("expected UGT1A1 diplotype call");
        };
        assert_eq!(diplotypes.len(), 4);

        let report_gene = ReportGene::from_standard_gene_call_result(&result, Some(&phenotype))
            .expect("report gene")
            .expect("UGT1A1 report gene");
        let source_labels = report_gene
            .source_diplotypes
            .iter()
            .map(|diplotype| diplotype.label.as_str())
            .collect::<Vec<_>>();
        let recommendation_labels = report_gene
            .recommendation_diplotypes
            .iter()
            .map(|diplotype| diplotype.label.as_str())
            .collect::<Vec<_>>();

        assert_eq!(report_gene.gene, "UGT1A1");
        assert_eq!(report_gene.source_diplotypes.len(), 4);
        assert_eq!(report_gene.recommendation_diplotypes.len(), 4);
        assert_eq!(source_labels, recommendation_labels);
        assert!(
            report_gene
                .source_diplotypes
                .iter()
                .all(|diplotype| diplotype.combination)
        );
        assert!(
            report_gene
                .recommendation_diplotypes
                .iter()
                .all(|diplotype| diplotype.combination)
        );
        assert!(source_labels.contains(&"*1/[*6 + *28 + g.233760973C>T]"));
        assert!(source_labels.contains(&"g.233760973C>T/[*6 + *28]"));
        assert!(source_labels.contains(&"*6/[*28 + g.233760973C>T]"));
        assert!(source_labels.contains(&"*28/[*6 + g.233760973C>T]"));
        assert_eq!(
            report_gene
                .source_diplotypes
                .iter()
                .find(|diplotype| diplotype.label == "*1/[*6 + *28 + g.233760973C>T]")
                .and_then(|diplotype| diplotype.match_score.as_deref()),
            diplotypes
                .iter()
                .find(|diplotype| diplotype.name == "*1/[*6 + *28 + g.233760973C>T]")
                .map(|diplotype| diplotype.score.to_string())
                .as_deref()
        );
    }

    #[test]
    fn html_amd_calls_at_positions_uses_java_function_map_for_combination_alleles() {
        let definition = read_definition_file(Path::new(UGT1A1_COMBINATION_DEFINITION_PATH))
            .expect("UGT1A1 definition");
        let phenotype =
            read_gene_phenotype_file(Path::new(UGT1A1_PHENOTYPE_PATH)).expect("UGT1A1 phenotype");
        let records = read_record_summaries(
            Path::new(UGT1A1_PARTIAL_WITH_COMBINATION_VCF_PATH),
            Some("PharmCAT"),
        )
        .expect("VCF");
        let allele_map = allele_map_from_record_summaries(
            records
                .records
                .into_iter()
                .filter_map(|record| record.allele_call),
        );
        let result = call_standard_gene("PharmCAT", &definition, &allele_map, false, true)
            .expect("standard gene call");
        let report_gene = ReportGene::from_standard_gene_call_result_with_definition(
            &result,
            Some(&phenotype),
            &definition,
        )
        .expect("report gene")
        .expect("UGT1A1 report gene");

        let html = super::html_amd_calls_at_positions(&report_gene);

        assert_eq!(
            report_gene
                .allele_function_map
                .get("*80+*28")
                .map(String::as_str),
            Some("Decreased function")
        );
        assert_eq!(
            report_gene
                .allele_function_map
                .get("*80+*37")
                .map(String::as_str),
            Some("Decreased function")
        );
        assert!(html.contains("<li>*80 - Unknown function</li>"));
        assert!(html.contains("<li>*6 - Decreased function</li>"));
        assert!(html.contains("<li>*27 - Decreased function</li>"));
        assert!(html.contains("<li>*28 - Decreased function</li>"));
        assert!(html.contains("<li>*36 - Increased function</li>"));
        assert!(html.contains("<li>*37 - Decreased function</li>"), "{html}");
        assert!(!html.contains("<li>*80 - Unassigned</li>"));
        assert!(!html.contains("<li>*27 - Unassigned</li>"));
        assert!(!html.contains("<li>*36 - Unassigned</li>"));
        assert!(!html.contains("<li>*37 - Unassigned</li>"));
        assert!(!html.contains("<li>*80+*28 - Unassigned</li>"));
        assert!(!html.contains("<li>*80+*37 - Unassigned</li>"));
    }

    #[test]
    fn report_gene_from_standard_combination_result_adds_static_combo_messages_like_java() {
        let definition = read_definition_file(Path::new(UGT1A1_COMBINATION_DEFINITION_PATH))
            .expect("UGT1A1 definition");
        let phenotype =
            read_gene_phenotype_file(Path::new(UGT1A1_PHENOTYPE_PATH)).expect("UGT1A1 phenotype");
        let catalog = MessageCatalog::from_path(Path::new(MESSAGES_PATH)).expect("messages");
        let records = read_record_summaries(
            Path::new(UGT1A1_PARTIAL_WITH_COMBINATION_VCF_PATH),
            Some("PharmCAT"),
        )
        .expect("VCF");
        let allele_map = allele_map_from_record_summaries(
            records
                .records
                .into_iter()
                .filter_map(|record| record.allele_call),
        );
        let result = call_standard_gene("PharmCAT", &definition, &allele_map, false, true)
            .expect("standard gene call");
        assert!(!result.match_data.phased);

        let report_gene = ReportGene::from_standard_gene_call_result_with_definition_and_messages(
            &result,
            Some(&phenotype),
            &definition,
            Some(&catalog),
        )
        .expect("report gene")
        .expect("UGT1A1 report gene");
        let message_names = report_gene
            .messages
            .iter()
            .map(|message| message.name.as_str())
            .collect::<BTreeSet<_>>();

        assert!(message_names.contains("pcat-combo-naming"));
        assert!(message_names.contains("pcat-combo-unphased"));
        assert_eq!(
            report_gene
                .messages
                .iter()
                .filter(|message| message.name.starts_with("pcat-combo-"))
                .count(),
            2
        );
    }

    #[test]
    fn report_gene_from_cyp2d6_matcher_result_adds_static_cyp2d6_messages_like_java() {
        let definition =
            read_definition_file(Path::new(CYP3A5_DEFINITION_PATH)).expect("fixture definition");
        let catalog = MessageCatalog::from_path(Path::new(MESSAGES_PATH)).expect("messages");
        let allele_map = BTreeMap::new();
        let result = GeneCallResult {
            gene: "CYP2D6".to_owned(),
            match_data: MatchData::new("Sample_1", "CYP2D6", &definition, &allele_map),
            kind: GeneCallKind::NoCall,
            dpyd_hap_b3_warnings: BTreeSet::new(),
            warnings: BTreeSet::new(),
        };

        let report_gene = ReportGene::from_standard_gene_call_result_with_definition_and_messages(
            &result,
            None,
            &definition,
            Some(&catalog),
        )
        .expect("report gene")
        .expect("CYP2D6 report gene");
        let message_names = report_gene
            .messages
            .iter()
            .map(|message| message.name.as_str())
            .collect::<BTreeSet<_>>();

        assert!(message_names.contains("pcat-cyp2d6-research-mode"));
        assert!(message_names.contains("pcat-cyp2d6-gene-note"));
        assert_eq!(report_gene.messages.len(), 2);
    }

    #[test]
    fn gene_message_matching_applies_haplotype_variant_missing_and_diplotype_rules_like_java() {
        let catalog = MessageCatalog::from_messages(vec![
            gene_message(
                "called-haplotype",
                MessageAnnotation::TYPE_NOTE,
                MatchLogic {
                    gene: Some("TEST".to_owned()),
                    haps_called: vec!["*2".to_owned()],
                    ..MatchLogic::default()
                },
            ),
            gene_message(
                "missing-haplotype",
                MessageAnnotation::TYPE_NOTE,
                MatchLogic {
                    gene: Some("TEST".to_owned()),
                    haps_missing: vec!["*4".to_owned()],
                    ..MatchLogic::default()
                },
            ),
            gene_message(
                "variant-present",
                MessageAnnotation::TYPE_NOTE,
                MatchLogic {
                    gene: Some("TEST".to_owned()),
                    variants: vec!["rsPresent".to_owned()],
                    ..MatchLogic::default()
                },
            ),
            gene_message(
                "variant-missing",
                MessageAnnotation::TYPE_NOTE,
                MatchLogic {
                    gene: Some("TEST".to_owned()),
                    variants_missing: vec!["rsMissing".to_owned()],
                    ..MatchLogic::default()
                },
            ),
            gene_message(
                "source-diplotype",
                MessageAnnotation::TYPE_NOTE,
                MatchLogic {
                    gene: Some("TEST".to_owned()),
                    dips: vec!["*1/*2".to_owned()],
                    ..MatchLogic::default()
                },
            ),
        ]);
        let mut report_gene = report_gene_for_gene_messages();
        report_gene.uncalled_haplotypes.insert("*4".to_owned());
        report_gene.variant_reports = vec![
            VariantReport::new("rsPresent", Some("A/G")),
            VariantReport::new("rsMissing", None::<String>),
        ];

        report_gene.apply_matching_gene_messages(&catalog);

        let message_names = report_gene
            .messages
            .iter()
            .map(|message| message.name.as_str())
            .collect::<BTreeSet<_>>();
        assert_eq!(
            message_names,
            BTreeSet::from([
                "called-haplotype",
                "missing-haplotype",
                "source-diplotype",
                "variant-missing",
                "variant-present",
            ])
        );
    }

    #[test]
    fn gene_message_matching_handles_ambiguity_and_non_match_gates_like_java() {
        let catalog = MessageCatalog::from_messages(vec![
            gene_message(
                "dip-ambiguity",
                MessageAnnotation::TYPE_AMBIGUITY,
                MatchLogic {
                    gene: Some("TEST".to_owned()),
                    dips: vec!["*1/*2".to_owned()],
                    ..MatchLogic::default()
                },
            ),
            gene_message(
                "variant-ambiguity",
                MessageAnnotation::TYPE_AMBIGUITY,
                MatchLogic {
                    gene: Some("TEST".to_owned()),
                    variants: vec!["rsAmbiguous".to_owned()],
                    ..MatchLogic::default()
                },
            ),
            gene_message(
                "non-match",
                MessageAnnotation::TYPE_NONMATCH,
                MatchLogic {
                    gene: Some("TEST".to_owned()),
                    variants: vec!["rsAmbiguous".to_owned()],
                    ..MatchLogic::default()
                },
            ),
        ]);

        let mut unphased = report_gene_for_gene_messages();
        unphased.variant_reports = vec![VariantReport::new("rsAmbiguous", Some("A/G"))];
        unphased.apply_matching_gene_messages(&catalog);
        assert!(
            unphased
                .messages
                .iter()
                .any(|message| message.name == "dip-ambiguity")
        );
        assert!(
            unphased
                .messages
                .iter()
                .any(|message| message.name == "variant-ambiguity")
        );
        assert!(
            !unphased
                .messages
                .iter()
                .any(|message| message.name == "non-match")
        );

        let mut phased = report_gene_for_gene_messages();
        phased.phased = true;
        phased.variant_reports = vec![VariantReport::new("rsAmbiguous", Some("A/G"))];
        phased.apply_matching_gene_messages(&catalog);
        assert!(
            !phased
                .messages
                .iter()
                .any(|message| message.name == "dip-ambiguity")
        );
        assert!(
            phased
                .messages
                .iter()
                .any(|message| message.name == "variant-ambiguity")
        );

        let mut homozygous = report_gene_for_gene_messages();
        homozygous.variant_reports = vec![VariantReport::new("rsAmbiguous", Some("A/A"))];
        homozygous.apply_matching_gene_messages(&catalog);
        assert!(
            !homozygous
                .messages
                .iter()
                .any(|message| message.name == "variant-ambiguity")
        );

        let mut non_reportable_with_data = ReportGene::new("TEST", Vec::<String>::new());
        non_reportable_with_data.variant_reports =
            vec![VariantReport::new("rsAmbiguous", Some("A/G"))];
        non_reportable_with_data.apply_matching_gene_messages(&catalog);
        assert!(
            non_reportable_with_data
                .messages
                .iter()
                .any(|message| message.name == "non-match")
        );

        let mut non_reportable_no_data = ReportGene::new("TEST", Vec::<String>::new());
        non_reportable_no_data.variant_reports =
            vec![VariantReport::new("rsAmbiguous", None::<String>)];
        non_reportable_no_data.apply_matching_gene_messages(&catalog);
        assert!(non_reportable_no_data.messages.is_empty());
    }

    #[test]
    fn report_gene_from_no_call_outside_call_emits_unknown_diplotypes_like_java() {
        let validation = OutsideCallValidation::for_supported_genes(["HLA-B"]);
        let outside_call =
            parse_outside_call_line(&validation, "HLA-B\t.", 1).expect("outside no-call");
        assert!(outside_call.is_no_call());

        let report_gene =
            ReportGene::from_outside_call(&outside_call, None).expect("outside report gene");
        let genotype = RecommendationGenotype::from_report_genes([report_gene.clone()]);

        assert_eq!(report_gene.gene, "HLA-B");
        assert_eq!(report_gene.lookup_keys, ["No Result"]);
        assert!(report_gene.phenotypes.is_empty());
        assert!(report_gene.outside_call);
        assert_eq!(
            report_gene.source_diplotype.as_deref(),
            Some("Unknown/Unknown")
        );
        assert_eq!(report_gene.source_diplotypes.len(), 1);
        assert_eq!(report_gene.recommendation_diplotypes.len(), 1);
        assert_eq!(report_gene.source_diplotypes[0].label, "Unknown/Unknown");
        assert_eq!(
            report_gene.recommendation_diplotypes[0].label,
            "Unknown/Unknown"
        );
        assert_eq!(genotype.lookup_keys()[0]["HLA-B"], Value::from("No Result"));
    }

    #[test]
    fn report_diplotype_marks_combination_labels_like_java_diplotype_factory() {
        let phenotype =
            read_gene_phenotype_file(Path::new(DPYD_PHENOTYPE_PATH)).expect("DPYD phenotype");

        let diplotype = ReportDiplotype::from_match_label(
            "DPYD",
            "[c.1129-5923C>G + c.1236G>A]/Reference",
            Some(42),
            Some(&phenotype),
        );

        assert_eq!(
            diplotype.label,
            "[c.1129-5923C>G + c.1236G>A] (heterozygous)"
        );
        assert!(diplotype.combination);
        assert!(!diplotype.inferred);
        assert!(diplotype.inferred_source_diplotypes.is_empty());
        assert_eq!(
            diplotype.allele1.as_ref().expect("allele1").name,
            "[c.1129-5923C>G + c.1236G>A]"
        );
        assert_eq!(
            diplotype.diplotype_key["[c.1129-5923C>G + c.1236G>A]"],
            Value::from(1.0)
        );
        assert_eq!(diplotype.match_score.as_deref(), Some("42"));
    }

    #[test]
    fn report_diplotype_sorting_matches_java_comparable_shape() {
        let mut inferred = ReportDiplotype::from_match_label("CYP2D6", "*2/*1", None, None);
        inferred.inferred = true;
        inferred.inferred_source_diplotypes = vec![
            ReportDiplotype::from_match_label("CYP2D6", "*4/*1", None, None),
            ReportDiplotype::from_match_label("CYP2D6", "*2/*1", None, None),
        ];
        let non_inferred_same_label =
            ReportDiplotype::from_match_label("CYP2D6", "*2/*1", None, None);

        let mut diplotypes = vec![
            ReportDiplotype::from_match_label("CYP2D6", "*4/*1", None, None),
            ReportDiplotype::from_match_label("CYP2D6", "*3/*1", None, None),
            inferred,
            non_inferred_same_label,
        ];

        sort_report_diplotypes(&mut diplotypes);

        assert_eq!(diplotypes[0].label, "*1/*2");
        assert!(!diplotypes[0].inferred);
        assert_eq!(diplotypes[1].label, "*1/*2");
        assert!(diplotypes[1].inferred);
        assert_eq!(diplotypes[2].label, "*1/*3");
        assert_eq!(diplotypes[3].label, "*1/*4");
        assert_eq!(diplotypes[1].inferred_source_diplotypes[0].label, "*1/*2");
        assert_eq!(diplotypes[1].inferred_source_diplotypes[1].label, "*1/*4");

        let inferred_json = serde_json::to_value(&diplotypes[1]).expect("inferred diplotype JSON");
        let inferred_sources = inferred_json["inferredSourceDiplotypes"]
            .as_array()
            .expect("inferredSourceDiplotypes")
            .iter()
            .map(|diplotype| diplotype["label"].as_str().expect("label"))
            .collect::<Vec<_>>();
        assert_eq!(inferred_sources, ["*1/*2", "*1/*4"]);
    }

    #[test]
    fn report_diplotype_serializes_outside_call_and_variant_fields_like_java() {
        let phenotype_map = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let cyp2d6 = phenotype_map.phenotype("CYP2D6").expect("CYP2D6");
        let annotated = DiplotypeAnnotationInput::outside_call(
            "CYP2D6",
            Some("*1"),
            Some("*1"),
            Some("Normal Metabolizer"),
            Some("4.0"),
        )
        .annotate(Some(cyp2d6))
        .expect("outside call annotation");

        let mut diplotype = ReportDiplotype::from_annotated(&annotated, Some(cyp2d6));
        diplotype.variant = Some(VariantReport::new("rs-test", Some("A/T")).with_position(101));

        assert_eq!(diplotype.label, "*1/*1");
        assert!(diplotype.outside_phenotype);
        assert!(diplotype.outside_phenotype_mismatch.is_some());
        assert!(diplotype.outside_activity_score);
        assert!(diplotype.outside_activity_score_mismatch.is_some());
        assert_eq!(
            diplotype
                .variant
                .as_ref()
                .and_then(|variant| variant.position),
            Some(101)
        );

        let value = serde_json::to_value(&diplotype).expect("diplotype JSON");
        assert_eq!(value["outsidePhenotype"], Value::from(true));
        assert!(value["outsidePhenotypeMismatch"].is_string());
        assert_eq!(value["outsideActivityScore"], Value::from(true));
        assert!(value["outsideActivityScoreMismatch"].is_string());
        assert_eq!(value["variant"]["position"], Value::from(101));
    }

    #[test]
    fn report_gene_from_ryr1_gene_call_result_feeds_recommendation_matching_like_java() {
        let definition =
            read_definition_file(Path::new(RYR1_DEFINITION_PATH)).expect("RYR1 definition");
        let phenotype =
            read_gene_phenotype_file(Path::new(RYR1_PHENOTYPE_PATH)).expect("RYR1 phenotype");
        let mut allele_map = reference_allele_map(&definition);
        let variant = variant_by_rsid(&definition, "rs193922746");
        allele_map.insert(
            variant.vcf_chr_position(),
            sample_call(
                &variant.chromosome,
                variant.position as usize,
                Some(&variant.reference),
                Some("G"),
                false,
                true,
            ),
        );

        let result =
            call_ryr1_lowest_function_gene("Sample_1", &definition, &allele_map).expect("RYR1");
        let GeneCallKind::Diplotypes(diplotypes) = &result.kind else {
            panic!("expected RYR1 diplotype call");
        };
        assert_eq!(diplotypes[0].name, "Reference/c.97A>G");

        let report_gene = ReportGene::from_gene_call_result_with_definition(
            &result,
            Some(&phenotype),
            &definition,
        )
        .expect("report gene")
        .expect("RYR1 report gene");
        let genotype = RecommendationGenotype::from_report_genes([report_gene.clone()]);
        let report_variant = report_gene
            .variant_reports
            .iter()
            .find(|report| report.db_snp_id.as_deref() == Some("rs193922746"))
            .expect("RYR1 variant report");
        let reference_message = report_gene
            .messages
            .iter()
            .find(|message| message.name == "reference-allele")
            .expect("reference allele message");
        let expected_variant_call = format!("{}/G", variant.reference);

        assert_eq!(report_gene.gene, "RYR1");
        assert_eq!(report_gene.chromosome.as_deref(), Some("chr19"));
        assert_eq!(
            report_gene.allele_definition_version,
            definition
                .version
                .clone()
                .or_else(|| definition.data_version.clone())
        );
        assert_eq!(
            report_gene.allele_definition_source,
            data_source_from_definition(&definition)
        );
        assert_eq!(report_gene.phenotype_version, phenotype.version);
        assert_eq!(report_gene.phased, result.match_data.phased);
        assert_eq!(
            report_gene.effectively_phased,
            result.match_data.effectively_phased
        );
        assert_eq!(report_gene.variant_reports.len(), definition.variants.len());
        assert_eq!(report_variant.gene.as_deref(), Some("RYR1"));
        assert_eq!(
            report_variant.chromosome.as_deref(),
            Some(variant.chromosome.as_str())
        );
        assert_eq!(report_variant.position, Some(variant.position as i64));
        assert_eq!(
            report_variant.call.as_deref(),
            Some(expected_variant_call.as_str())
        );
        assert_eq!(
            report_variant.reference_allele.as_deref(),
            Some(variant.reference.as_str())
        );
        assert!(!report_variant.phased);
        assert_eq!(
            reference_message.exception_type,
            MessageAnnotation::TYPE_NOTE
        );
        assert_eq!(
            reference_message.message,
            "The RYR1 Reference allele assignment is characterized by the absence of variants at the positions that are included in the underlying allele definitions."
        );
        assert_eq!(
            report_gene.lookup_keys,
            ["Malignant Hyperthermia Susceptibility"]
        );
        assert_eq!(
            report_gene.phenotypes,
            ["Malignant Hyperthermia Susceptibility"]
        );
        assert_eq!(
            report_gene.source_diplotype.as_deref(),
            Some("c.97A>G (heterozygous)")
        );
        assert_eq!(report_gene.source_diplotypes.len(), 1);
        assert_eq!(
            report_gene.source_diplotypes[0].label,
            "c.97A>G (heterozygous)"
        );
        let source_allele1 = report_gene.source_diplotypes[0]
            .allele1
            .as_ref()
            .expect("source allele1");
        let source_allele2 = report_gene.source_diplotypes[0]
            .allele2
            .as_ref()
            .expect("source allele2");
        assert_eq!(source_allele1.name, "Reference");
        assert_eq!(source_allele1.function, "Normal function");
        assert!(source_allele1.reference);
        assert_eq!(source_allele1.activity_value, None);
        assert_eq!(source_allele2.name, "c.97A>G");
        assert_eq!(source_allele2.function, "Malignant Hyperthermia associated");
        assert!(!source_allele2.reference);
        assert_eq!(source_allele2.activity_value, None);
        // Java DiplotypeMatch score sums the two NamedAllele scores from the definition JSON:
        // Reference (313) + c.97A>G (1) = 314. (The earlier 626 was a wrong 2xReference guess.)
        assert_eq!(
            report_gene.source_diplotypes[0].match_score.as_deref(),
            Some("314")
        );
        assert!(!report_gene.source_diplotypes[0].inferred);
        assert!(
            report_gene.source_diplotypes[0]
                .inferred_source_diplotypes
                .is_empty()
        );
        assert!(!report_gene.source_diplotypes[0].combination);
        assert_eq!(
            report_gene.matcher_component_haplotypes,
            ["Reference".to_owned(), "c.97A>G".to_owned()]
                .into_iter()
                .collect()
        );
        assert_eq!(report_gene.matcher_component_diplotypes.len(), 2);
        assert_eq!(
            report_gene.matcher_component_diplotypes[0].label,
            "Reference"
        );
        assert_eq!(report_gene.matcher_component_diplotypes[1].label, "c.97A>G");
        assert_eq!(
            report_gene.matcher_component_diplotypes[1]
                .allele1
                .as_ref()
                .expect("component allele1")
                .function,
            "Malignant Hyperthermia associated"
        );
        assert_eq!(report_gene.recommendation_diplotypes.len(), 1);
        assert_eq!(
            report_gene.recommendation_diplotypes[0].label,
            "c.97A>G (heterozygous)"
        );
        assert!(!report_gene.recommendation_diplotypes[0].inferred);
        assert_eq!(
            report_gene.recommendation_diplotypes[0]
                .inferred_source_diplotypes
                .len(),
            1
        );
        assert_eq!(
            report_gene.recommendation_diplotypes[0].inferred_source_diplotypes[0].label,
            "c.97A>G (heterozygous)"
        );
        assert!(!report_gene.recommendation_diplotypes[0].combination);
        assert_eq!(
            report_gene.recommendation_diplotypes[0].phenotypes,
            ["Malignant Hyperthermia Susceptibility"]
        );
        assert_eq!(
            report_gene.recommendation_diplotypes[0]
                .allele2
                .as_ref()
                .expect("recommendation allele2")
                .display_name(),
            "RYR1 c.97A>G"
        );
        assert_eq!(report_gene.diplotype_key["Reference"], Value::from(1.0));
        assert_eq!(report_gene.diplotype_key["c.97A>G"], Value::from(1.0));
        assert_eq!(
            genotype.lookup_keys()[0]["RYR1"],
            Value::from("Malignant Hyperthermia Susceptibility")
        );
    }

    #[test]
    fn report_context_from_gene_call_results_uses_gene_specific_matcher_report_paths() {
        let collection =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");
        let phenotype_map = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let cyp3a5_definition =
            read_definition_file(Path::new(CYP3A5_DEFINITION_PATH)).expect("CYP3A5 definition");
        let dpyd_definition =
            read_definition_file(Path::new(DPYD_DEFINITION_PATH)).expect("DPYD definition");
        let ryr1_definition =
            read_definition_file(Path::new(RYR1_DEFINITION_PATH)).expect("RYR1 definition");
        let definition_reader = DefinitionReader::from_definitions(
            [
                (
                    cyp3a5_definition.gene_symbol.clone(),
                    cyp3a5_definition.clone(),
                ),
                (dpyd_definition.gene_symbol.clone(), dpyd_definition.clone()),
                (ryr1_definition.gene_symbol.clone(), ryr1_definition.clone()),
            ]
            .into_iter()
            .collect(),
        );

        let cyp3a5_records =
            read_record_summaries(Path::new(CYP3A5_VCF_PATH), Some("NA12878")).expect("VCF");
        let cyp3a5_allele_map = allele_map_from_record_summaries(
            cyp3a5_records
                .records
                .into_iter()
                .filter_map(|record| record.allele_call),
        );
        let cyp3a5_result = call_standard_gene(
            "NA12878",
            &cyp3a5_definition,
            &cyp3a5_allele_map,
            false,
            false,
        )
        .expect("CYP3A5 call");

        let mut dpyd_allele_map = reference_allele_map(&dpyd_definition);
        let dpyd_variant = variant_by_rsid(&dpyd_definition, "rs67376798");
        dpyd_allele_map.insert(
            dpyd_variant.vcf_chr_position(),
            sample_call(
                &dpyd_variant.chromosome,
                dpyd_variant.position as usize,
                Some(&dpyd_variant.reference),
                Some("A"),
                true,
                true,
            ),
        );
        let dpyd_result =
            call_dpyd_lowest_function_gene("Sample_1", &dpyd_definition, &dpyd_allele_map)
                .expect("DPYD call");

        let mut ryr1_allele_map = reference_allele_map(&ryr1_definition);
        let ryr1_variant = variant_by_rsid(&ryr1_definition, "rs193922746");
        ryr1_allele_map.insert(
            ryr1_variant.vcf_chr_position(),
            sample_call(
                &ryr1_variant.chromosome,
                ryr1_variant.position as usize,
                Some(&ryr1_variant.reference),
                Some("G"),
                false,
                true,
            ),
        );
        let ryr1_result =
            call_ryr1_lowest_function_gene("Sample_1", &ryr1_definition, &ryr1_allele_map)
                .expect("RYR1 call");

        let context = ReportContext::from_gene_call_results_with_definitions(
            &collection,
            [&cyp3a5_result, &dpyd_result, &ryr1_result],
            &definition_reader,
            &phenotype_map,
            Some("Sample_1".to_owned()),
        )
        .expect("report context");

        assert_eq!(context.title.as_deref(), Some("Sample_1"));
        let cyp3a5 = context.gene_report("CYP3A5").expect("CYP3A5 report");
        assert_eq!(cyp3a5.source_diplotype.as_deref(), Some("*1/*2"));
        assert_eq!(cyp3a5.chromosome.as_deref(), Some("chr1"));
        assert!(
            cyp3a5
                .messages
                .iter()
                .any(|message| message.name == "reference-allele")
        );

        let dpyd = context.gene_report("DPYD").expect("DPYD report");
        assert_eq!(
            dpyd.source_diplotype.as_deref(),
            Some("c.2846A>T (heterozygous)")
        );
        assert_eq!(dpyd.lookup_keys, ["1.5"]);
        assert_eq!(dpyd.chromosome.as_deref(), Some("chr1"));
        assert_eq!(dpyd.variant_reports.len(), dpyd_definition.variants.len());
        assert!(
            dpyd.variant_reports
                .iter()
                .any(|report| report.db_snp_id.as_deref() == Some("rs67376798"))
        );
        assert!(
            dpyd.messages
                .iter()
                .any(|message| message.name == "reference-allele")
        );

        let ryr1 = context.gene_report("RYR1").expect("RYR1 report");
        assert_eq!(
            ryr1.source_diplotype.as_deref(),
            Some("c.97A>G (heterozygous)")
        );
        assert_eq!(ryr1.lookup_keys, ["Malignant Hyperthermia Susceptibility"]);
        assert_eq!(ryr1.chromosome.as_deref(), Some("chr19"));
        assert_eq!(ryr1.variant_reports.len(), ryr1_definition.variants.len());
        assert!(
            ryr1.variant_reports
                .iter()
                .any(|report| report.db_snp_id.as_deref() == Some("rs193922746"))
        );
        assert!(
            ryr1.messages
                .iter()
                .any(|message| message.name == "reference-allele")
        );
    }

    #[test]
    fn report_context_from_gene_call_results_fails_when_definition_is_missing() {
        let collection =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");
        let definition =
            read_definition_file(Path::new(DPYD_DEFINITION_PATH)).expect("DPYD definition");
        let result = call_dpyd_lowest_function_gene("Sample_1", &definition, &BTreeMap::new())
            .expect("DPYD no-call");
        let empty_definitions = DefinitionReader::from_definitions(BTreeMap::new());
        let empty_phenotypes = PhenotypeMap::default();

        let error = ReportContext::from_gene_call_results_with_definitions(
            &collection,
            [&result],
            &empty_definitions,
            &empty_phenotypes,
            None,
        )
        .expect_err("missing definition");

        match error {
            ReportContextFromMatcherError::MissingDefinition(gene) => {
                assert_eq!(gene, "DPYD");
            }
            ReportContextFromMatcherError::ReportGene(error) => {
                panic!("expected missing definition, got {error}");
            }
        }
    }

    #[test]
    fn annotation_report_adds_outside_call_mismatch_note_like_java_html_recommendation() {
        let validation = OutsideCallValidation::for_supported_genes(["CYP2D6"]);
        let call =
            parse_outside_call_line(&validation, "CYP2D6\t*1/*1\tNormal Metabolizer\t4.0", 1)
                .expect("outside call");
        let phenotype_map = PhenotypeMap::from_dir(Path::new(PHENOTYPE_PATH)).expect("phenotypes");
        let cyp2d6 = phenotype_map.phenotype("CYP2D6").expect("CYP2D6");
        let report_gene =
            ReportGene::from_outside_call(&call, Some(cyp2d6)).expect("outside report gene");
        let recommendation = RecommendationAnnotation {
            id: "PA-test".to_owned(),
            name: "Recommendation PA-test".to_owned(),
            population: None,
            classification: None,
            related_chemicals: Vec::new(),
            text: None,
            implications: Vec::new(),
            lookup_key: vec![
                [("CYP2D6".to_owned(), Value::String("4.0".to_owned()))]
                    .into_iter()
                    .collect(),
            ],
            dosing_information: false,
            alternate_drug_available: false,
            other_prescribing_guidance: false,
        };
        let annotation = AnnotationReport::new(
            &recommendation,
            "test-local-id".to_owned(),
            &[RecommendationGenotype::from_report_genes([report_gene])],
            &[],
        );

        let message = annotation.messages.iter().next().expect("mismatch message");
        assert_eq!(message.name, "warn.mismatch.outsideCall");
        assert_eq!(message.exception_type, MessageAnnotation::TYPE_NOTE);
        assert_eq!(
            message.message,
            "Conflicting outside call data was provided for CYP2D6.  PharmCAT will use provided activity score to match recommendations."
        );
    }

    #[test]
    fn report_context_merges_annotated_outside_calls_for_same_gene_like_java_add_outside_call() {
        let validation = OutsideCallValidation::for_supported_genes(["HLA-B"]);
        let calls = parse_outside_calls_str(
            &validation,
            "HLA-B\t*57:01\t*57:01 positive\nHLA-B\t*58:01\t*58:01 positive\n",
        )
        .expect("outside calls");
        let collection =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");
        let report_genes = calls
            .iter()
            .enumerate()
            .map(|(index, call)| {
                ReportGene::from_annotated_diplotype(
                    call.to_annotation_input()
                        .expect("annotation input")
                        .annotate(None)
                        .expect("annotated HLA outside call"),
                )
                .with_messages([MessageAnnotation::new_note(
                    format!("hla-b-message-{index}"),
                    format!("message {index}"),
                )])
            })
            .collect::<Vec<_>>();

        let context = ReportContext::from_gene_reports(&collection, report_genes, None);

        assert_eq!(
            context
                .gene_report("HLA-B")
                .expect("HLA-B")
                .phenotypes
                .as_slice(),
            ["*57:01 positive", "*58:01 positive"]
        );
        assert_eq!(
            context
                .gene_report("HLA-B")
                .expect("HLA-B")
                .source_diplotype
                .as_deref(),
            Some("*57:01; *58:01")
        );
        assert_eq!(
            context
                .gene_report("HLA-B")
                .expect("HLA-B")
                .source_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["*57:01", "*58:01"]
        );
        assert_eq!(
            context
                .gene_report("HLA-B")
                .expect("HLA-B")
                .recommendation_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["*57:01", "*58:01"]
        );
        assert!(
            context
                .drug_report(PrescribingGuidanceSource::CpicGuideline, "allopurinol")
                .expect("allopurinol")
                .is_matched()
        );
        assert_eq!(
            context
                .gene_report("HLA-B")
                .expect("HLA-B")
                .messages
                .iter()
                .map(|message| message.name.as_str())
                .collect::<Vec<_>>(),
            ["hla-b-message-0", "hla-b-message-1"]
        );
    }

    #[test]
    fn merge_report_gene_prefers_outside_summary_fields_like_java_gene_report_order() {
        let mut matcher = ReportGene::new("CYP2D6", ["1.0".to_owned()])
            .with_phenotypes(["Intermediate Metabolizer".to_owned()])
            .with_activity_score("1.0")
            .with_source_diplotype("*1/*4")
            .with_match_score("42")
            .with_variant_reports([VariantReport::new("rsMatcher", Some("A/T"))])
            .with_messages([MessageAnnotation::new_note(
                "matcher-message",
                "matcher message",
            )]);
        matcher.source_diplotypes = vec![ReportDiplotype {
            gene: "CYP2D6".to_owned(),
            label: "*1/*4".to_owned(),
            phenotypes: vec!["Intermediate Metabolizer".to_owned()],
            activity_score: Some("1.0".to_owned()),
            match_score: Some("42".to_owned()),
            ..ReportDiplotype::default()
        }];
        matcher.recommendation_diplotypes = matcher.source_diplotypes.clone();
        matcher
            .related_drugs
            .insert(super::DrugLink::new("matcher drug", "PA-matcher"));

        let mut outside = ReportGene::new("CYP2D6", ["2.0".to_owned()])
            .with_phenotypes(["Ultrarapid Metabolizer".to_owned()])
            .with_activity_score("2.0")
            .with_source_diplotype("*1x2/*1")
            .with_outside_call(true)
            .with_variant_reports([VariantReport::new("rsOutside", Some("G/G"))])
            .with_messages([MessageAnnotation::new_note(
                "outside-message",
                "outside message",
            )]);
        outside.source_diplotypes = vec![ReportDiplotype {
            gene: "CYP2D6".to_owned(),
            label: "*1x2/*1".to_owned(),
            phenotypes: vec!["Ultrarapid Metabolizer".to_owned()],
            activity_score: Some("2.0".to_owned()),
            outside_activity_score: true,
            ..ReportDiplotype::default()
        }];
        outside.recommendation_diplotypes = outside.source_diplotypes.clone();
        outside
            .related_drugs
            .insert(super::DrugLink::new("outside drug", "PA-outside"));

        super::merge_report_gene(&mut matcher, outside);

        assert!(matcher.outside_call);
        assert_eq!(matcher.lookup_keys, ["2.0"]);
        assert_eq!(matcher.phenotypes, ["Ultrarapid Metabolizer"]);
        assert_eq!(matcher.activity_score.as_deref(), Some("2.0"));
        assert_eq!(matcher.source_diplotype.as_deref(), Some("*1x2/*1"));
        assert_eq!(matcher.match_score, None);
        assert_eq!(
            matcher
                .source_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["*1x2/*1"]
        );
        assert_eq!(
            matcher
                .recommendation_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["*1x2/*1"]
        );
        assert_eq!(
            matcher
                .related_drugs
                .iter()
                .map(|drug| drug.name.as_str())
                .collect::<Vec<_>>(),
            ["matcher drug", "outside drug"]
        );
        assert_eq!(
            matcher
                .messages
                .iter()
                .map(|message| message.name.as_str())
                .collect::<Vec<_>>(),
            ["matcher-message", "outside-message"]
        );
        assert_eq!(
            matcher
                .variant_reports
                .iter()
                .map(|variant| variant.db_snp_id.as_deref())
                .collect::<Vec<_>>(),
            [Some("rsOutside"), Some("rsMatcher")]
        );
    }

    #[test]
    fn genotype_summary_uses_preserved_same_gene_source_reports_like_java_build_genotype_summary() {
        let collection =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");
        let mut matcher = ReportGene::new("GENEX", ["matcher-key".to_owned()])
            .with_source_diplotype("matcher diplotype")
            .with_messages([MessageAnnotation::new_note(
                "matcher-summary-message",
                "matcher summary message",
            )]);
        matcher.source_diplotypes = vec![ReportDiplotype {
            gene: "GENEX".to_owned(),
            label: "matcher diplotype".to_owned(),
            lookup_keys: vec!["matcher-key".to_owned()],
            phenotypes: vec!["Matcher Phenotype".to_owned()],
            ..ReportDiplotype::default()
        }];
        matcher.recommendation_diplotypes = matcher.source_diplotypes.clone();
        matcher
            .related_drugs
            .insert(super::DrugLink::new("matcher drug", "PA-matcher"));

        let mut outside = ReportGene::new("GENEX", ["outside-key".to_owned()])
            .with_source_diplotype("outside diplotype")
            .with_outside_call(true)
            .with_messages([MessageAnnotation::new_note(
                "outside-summary-message",
                "outside summary message",
            )]);
        outside.source_diplotypes = vec![ReportDiplotype {
            gene: "GENEX".to_owned(),
            label: "outside diplotype".to_owned(),
            lookup_keys: vec!["outside-key".to_owned()],
            phenotypes: vec!["Outside Phenotype".to_owned()],
            outside_phenotype: true,
            ..ReportDiplotype::default()
        }];
        outside.recommendation_diplotypes = outside.source_diplotypes.clone();
        outside
            .related_drugs
            .insert(super::DrugLink::new("outside drug", "PA-outside"));

        let context = ReportContext::from_gene_reports(&collection, [matcher, outside], None);
        let source_reports = context
            .report_gene_sources
            .get("GENEX")
            .expect("GENEX source reports");

        assert_eq!(source_reports.len(), 2);
        assert!(source_reports[0].outside_call);
        assert!(!source_reports[1].outside_call);

        let summaries =
            super::html_genotype_summary_report_genes(&context, &HtmlReportOptions::default());
        assert_eq!(summaries.len(), 1);
        let summary = &summaries[0];
        assert_eq!(summary.lookup_keys, ["outside-key"]);
        assert_eq!(
            summary.source_diplotype.as_deref(),
            Some("outside diplotype")
        );
        assert_eq!(
            summary
                .recommendation_diplotypes
                .iter()
                .map(|diplotype| diplotype.label.as_str())
                .collect::<Vec<_>>(),
            ["outside diplotype"]
        );
        assert_eq!(
            summary
                .related_drugs
                .iter()
                .map(|drug| drug.name.as_str())
                .collect::<Vec<_>>(),
            ["matcher drug", "outside drug"]
        );
        assert_eq!(
            summary
                .messages
                .iter()
                .map(|message| message.name.as_str())
                .collect::<Vec<_>>(),
            ["matcher-summary-message", "outside-summary-message"]
        );
    }

    #[test]
    fn genotype_summary_orders_same_gene_matcher_sources_cpic_before_dpwg() {
        let collection =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");
        let mut dpwg = ReportGene::new("GENEY", ["dpwg-key".to_owned()])
            .with_source_diplotype("dpwg diplotype")
            .with_call_source(ReportCallSource::Matcher)
            .with_guidance_source(PrescribingGuidanceSource::DpwgGuideline);
        dpwg.source_diplotypes = vec![ReportDiplotype {
            gene: "GENEY".to_owned(),
            label: "dpwg diplotype".to_owned(),
            lookup_keys: vec!["dpwg-key".to_owned()],
            phenotypes: vec!["DPWG Phenotype".to_owned()],
            ..ReportDiplotype::default()
        }];
        dpwg.recommendation_diplotypes = dpwg.source_diplotypes.clone();
        dpwg.related_drugs
            .insert(super::DrugLink::new("dpwg drug", "PA-dpwg"));

        let mut cpic = ReportGene::new("GENEY", ["cpic-key".to_owned()])
            .with_source_diplotype("cpic diplotype")
            .with_call_source(ReportCallSource::Matcher)
            .with_guidance_source(PrescribingGuidanceSource::CpicGuideline);
        cpic.source_diplotypes = vec![ReportDiplotype {
            gene: "GENEY".to_owned(),
            label: "cpic diplotype".to_owned(),
            lookup_keys: vec!["cpic-key".to_owned()],
            phenotypes: vec!["CPIC Phenotype".to_owned()],
            ..ReportDiplotype::default()
        }];
        cpic.recommendation_diplotypes = cpic.source_diplotypes.clone();
        cpic.related_drugs
            .insert(super::DrugLink::new("cpic drug", "PA-cpic"));

        let context = ReportContext::from_gene_reports(&collection, [dpwg, cpic], None);
        let source_reports = context
            .report_gene_sources
            .get("GENEY")
            .expect("GENEY source reports");

        assert_eq!(
            source_reports
                .iter()
                .map(|report_gene| report_gene.guidance_source)
                .collect::<Vec<_>>(),
            [
                Some(PrescribingGuidanceSource::CpicGuideline),
                Some(PrescribingGuidanceSource::DpwgGuideline)
            ]
        );

        let summaries =
            super::html_genotype_summary_report_genes(&context, &HtmlReportOptions::default());
        assert_eq!(summaries.len(), 1);
        let summary = &summaries[0];
        assert_eq!(summary.lookup_keys, ["cpic-key"]);
        assert_eq!(summary.source_diplotype.as_deref(), Some("cpic diplotype"));
        assert_eq!(
            summary
                .related_drugs
                .iter()
                .map(|drug| drug.name.as_str())
                .collect::<Vec<_>>(),
            ["cpic drug", "dpwg drug"]
        );
    }

    #[test]
    fn report_context_populates_guidance_source_from_real_guideline_packages() {
        let collection =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");
        let mut abcg2 = ReportGene::new("ABCG2", ["Normal Function".to_owned()])
            .with_source_diplotype("Normal Function")
            .with_call_source(ReportCallSource::Matcher);
        abcg2.source_diplotypes = vec![ReportDiplotype {
            gene: "ABCG2".to_owned(),
            label: "Normal Function".to_owned(),
            lookup_keys: vec!["Normal Function".to_owned()],
            phenotypes: vec!["Normal Function".to_owned()],
            ..ReportDiplotype::default()
        }];
        abcg2.recommendation_diplotypes = abcg2.source_diplotypes.clone();

        let context = ReportContext::from_gene_reports(&collection, [abcg2], None);
        let source_reports = context
            .report_gene_sources
            .get("ABCG2")
            .expect("ABCG2 source reports");
        let guidance_sources = source_reports
            .iter()
            .filter_map(|report_gene| report_gene.guidance_source)
            .collect::<Vec<_>>();

        assert_eq!(
            guidance_sources.first(),
            Some(&PrescribingGuidanceSource::CpicGuideline)
        );
        assert!(guidance_sources.contains(&PrescribingGuidanceSource::DpwgGuideline));
        assert_eq!(
            context
                .gene_report("ABCG2")
                .expect("ABCG2 merged report")
                .related_drugs,
            [
                super::DrugLink::new("allopurinol", "PA448320"),
                super::DrugLink::new("rosuvastatin", "PA134308647"),
            ]
            .into_iter()
            .collect()
        );
        assert_eq!(
            source_reports[0].related_drugs,
            [super::DrugLink::new("rosuvastatin", "PA134308647")]
                .into_iter()
                .collect()
        );
        assert_eq!(
            source_reports[1].related_drugs,
            [super::DrugLink::new("allopurinol", "PA448320")]
                .into_iter()
                .collect()
        );

        let summaries =
            super::html_genotype_summary_report_genes(&context, &HtmlReportOptions::default());
        assert_eq!(summaries.len(), 1);
        assert_eq!(
            summaries[0].guidance_source,
            Some(PrescribingGuidanceSource::CpicGuideline)
        );
        let compact_drugs_with_recommendations = super::html_drugs_with_recommendations(
            &context,
            &HtmlReportOptions {
                compact: true,
                ..HtmlReportOptions::default()
            },
        );
        let compact_drugs = super::html_genotype_summary_drugs(
            &summaries[0],
            &compact_drugs_with_recommendations,
            &BTreeMap::new(),
        );
        assert!(!compact_drugs.contains("rosuvastatin"));
        assert!(compact_drugs.contains("allopurinol"));
    }

    #[test]
    fn compact_recommendation_drugs_ignore_non_reportable_warfarin_escape() {
        let collection =
            PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH)).expect("guidance");
        let mut abcg2 = ReportGene::new("ABCG2", ["Normal Function".to_owned()])
            .with_source_diplotype("Normal Function")
            .with_call_source(ReportCallSource::Matcher);
        abcg2.source_diplotypes = vec![ReportDiplotype {
            gene: "ABCG2".to_owned(),
            label: "Normal Function".to_owned(),
            lookup_keys: vec!["Normal Function".to_owned()],
            phenotypes: vec!["Normal Function".to_owned()],
            ..ReportDiplotype::default()
        }];
        abcg2.recommendation_diplotypes = abcg2.source_diplotypes.clone();

        let context = ReportContext::from_gene_reports(&collection, [abcg2], None);
        let compact_options = HtmlReportOptions {
            compact: true,
            ..HtmlReportOptions::default()
        };
        let warfarin = context
            .drug_report(PrescribingGuidanceSource::CpicGuideline, "warfarin")
            .expect("CPIC warfarin report");

        assert!(warfarin.is_matched());
        assert!(
            !warfarin
                .guidelines
                .iter()
                .any(GuidelineReport::is_reportable)
        );

        let compact_drugs = super::html_recommendation_drugs(&context, &compact_options);
        assert!(!compact_drugs.contains("warfarin"));
        assert!(super::html_reports_for_drug(&context, "warfarin", &compact_options).is_empty());
    }

    #[test]
    fn annotation_report_extracts_activity_scores_from_matched_genotype_like_java() {
        let recommendation = RecommendationAnnotation {
            id: "PA-test".to_owned(),
            name: "Recommendation PA-test".to_owned(),
            population: None,
            classification: None,
            related_chemicals: Vec::new(),
            text: None,
            implications: Vec::new(),
            lookup_key: vec![
                [("CYP2D6".to_owned(), Value::String("1.0".to_owned()))]
                    .into_iter()
                    .collect(),
            ],
            dosing_information: false,
            alternate_drug_available: false,
            other_prescribing_guidance: false,
        };
        let genotype = RecommendationGenotype::from_report_genes([ReportGene::new(
            "CYP2D6",
            ["1.0".to_owned()],
        )
        .with_phenotypes(["Intermediate Metabolizer".to_owned()])
        .with_activity_score("1.0")]);

        let annotation = AnnotationReport::new(
            &recommendation,
            "test-local-id".to_owned(),
            &[genotype],
            &[],
        );

        assert_eq!(
            annotation.phenotypes.get("CYP2D6").map(String::as_str),
            Some("Intermediate Metabolizer")
        );
        assert_eq!(
            annotation.activity_scores.get("CYP2D6").map(String::as_str),
            Some("1.0")
        );
    }

    #[test]
    fn recommendation_genotypes_expand_each_recommendation_diplotype_like_java() {
        let mut cyp2c19 = ReportGene::new(
            "CYP2C19",
            [
                "Poor Metabolizer".to_owned(),
                "Intermediate Metabolizer".to_owned(),
                "Ultrarapid Metabolizer".to_owned(),
            ],
        );
        cyp2c19.recommendation_diplotypes = [
            ("*4/*4", "Poor Metabolizer"),
            ("*4/*17", "Intermediate Metabolizer"),
            ("*17/*17", "Ultrarapid Metabolizer"),
        ]
        .into_iter()
        .map(|(label, phenotype)| ReportDiplotype {
            gene: "CYP2C19".to_owned(),
            label: label.to_owned(),
            lookup_keys: vec![phenotype.to_owned()],
            phenotypes: vec![phenotype.to_owned()],
            ..ReportDiplotype::default()
        })
        .collect();
        let context = ReportContext {
            title: None,
            data_version: "v-test".to_owned(),
            gene_reports: [("CYP2C19".to_owned(), cyp2c19)].into_iter().collect(),
            report_gene_sources: BTreeMap::new(),
            drug_reports: BTreeMap::new(),
            messages: Vec::new(),
            unannotated_gene_calls: Vec::new(),
        };

        let genotypes =
            make_recommendation_genotypes(&["CYP2C19".to_owned()].into_iter().collect(), &context);

        assert_eq!(genotypes.len(), 3);
        assert_eq!(
            genotypes
                .iter()
                .map(|genotype| genotype.lookup_keys()[0]["CYP2C19"]
                    .as_str()
                    .expect("lookup key"))
                .collect::<Vec<_>>(),
            [
                "Poor Metabolizer",
                "Intermediate Metabolizer",
                "Ultrarapid Metabolizer"
            ]
        );
        assert_eq!(
            genotypes
                .iter()
                .map(|genotype| genotype.report_genes()[0]
                    .phenotypes
                    .iter()
                    .map(String::as_str)
                    .collect::<Vec<_>>())
                .collect::<Vec<_>>(),
            [
                vec!["Poor Metabolizer"],
                vec!["Intermediate Metabolizer"],
                vec!["Ultrarapid Metabolizer"],
            ]
        );
    }

    #[test]
    fn annotation_report_marks_non_activity_score_gene_na_when_genotype_uses_activity_score() {
        let recommendation = RecommendationAnnotation {
            id: "PA-test".to_owned(),
            name: "Recommendation PA-test".to_owned(),
            population: None,
            classification: None,
            related_chemicals: Vec::new(),
            text: None,
            implications: Vec::new(),
            lookup_key: vec![
                [
                    ("CYP2D6".to_owned(), Value::String("1.0".to_owned())),
                    ("GENEX".to_owned(), Value::String("Normal".to_owned())),
                ]
                .into_iter()
                .collect(),
            ],
            dosing_information: false,
            alternate_drug_available: false,
            other_prescribing_guidance: false,
        };
        let genotype = RecommendationGenotype::from_report_genes([
            ReportGene::new("CYP2D6", ["1.0".to_owned()])
                .with_phenotypes(["Intermediate Metabolizer".to_owned()])
                .with_activity_score("1.0"),
            ReportGene::new("GENEX", ["Normal".to_owned()]),
        ]);

        let annotation = AnnotationReport::new(
            &recommendation,
            "test-local-id".to_owned(),
            &[genotype],
            &[],
        );

        assert_eq!(
            annotation.activity_scores.get("GENEX").map(String::as_str),
            Some(crate::phenotype::NA)
        );
    }

    #[test]
    #[should_panic(expected = "Multiple phenotypes for gene GENEX")]
    fn annotation_report_rejects_conflicting_non_allele_presence_phenotypes_like_java() {
        let recommendation = RecommendationAnnotation {
            id: "PA-test".to_owned(),
            name: "Recommendation PA-test".to_owned(),
            population: None,
            classification: None,
            related_chemicals: Vec::new(),
            text: None,
            implications: Vec::new(),
            lookup_key: vec![
                [("GENEX".to_owned(), Value::String("Normal".to_owned()))]
                    .into_iter()
                    .collect(),
            ],
            dosing_information: false,
            alternate_drug_available: false,
            other_prescribing_guidance: false,
        };
        let genotype = RecommendationGenotype::from_report_genes([ReportGene::new(
            "GENEX",
            ["Normal".to_owned()],
        )
        .with_phenotypes(["Normal".to_owned(), "Abnormal".to_owned()])]);

        let _ = AnnotationReport::new(
            &recommendation,
            "test-local-id".to_owned(),
            &[genotype],
            &[],
        );
    }

    #[test]
    fn report_context_adds_synthetic_cpic_warfarin_annotation_like_java() {
        let collection = PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH))
            .expect("prescribing guidance");
        let context = ReportContext::from_gene_reports(
            &collection,
            [ReportGene::new("CYP2C9", ["1.0".to_owned()])
                .with_phenotypes(["Intermediate Metabolizer".to_owned()])
                .with_activity_score("1.0")],
            None,
        );

        let warfarin = context
            .drug_report(PrescribingGuidanceSource::CpicGuideline, "warfarin")
            .expect("warfarin CPIC report");

        assert_eq!(warfarin.name, "warfarin");
        assert_eq!(warfarin.id, "PA451906");
        assert_eq!(warfarin.source, PrescribingGuidanceSource::CpicGuideline);
        assert!(warfarin.is_matched());
        assert_eq!(warfarin.matched_annotation_count(), 1);

        let guideline = warfarin.guidelines.iter().next().expect("guideline");
        assert_eq!(
            guideline.name,
            "Annotation of CPIC Guideline for warfarin and CYP2C9, CYP4F2, VKORC1"
        );
        assert_eq!(guideline.genes, ["CYP2C9".to_owned()].into_iter().collect());

        let annotation = guideline.annotations.iter().next().expect("annotation");
        assert_eq!(annotation.local_id, "warfarin-cpic-1-1");
        assert_eq!(annotation.genotypes.len(), 1);
        assert_eq!(annotation.drug_recommendation, None);
        assert!(annotation.lookup_key.is_empty());
        assert!(annotation.phenotypes.is_empty());
        assert!(annotation.activity_scores.is_empty());
    }

    #[test]
    fn report_context_matches_exact_diplotype_key_before_phenotype_lookup_like_java() {
        let collection = PgkbGuidelineCollection::from_path(Path::new(GUIDANCE_PATH))
            .expect("prescribing guidance");
        let context = ReportContext::from_gene_reports(
            &collection,
            [ReportGene::with_diplotype_counts(
                "UGT1A1",
                Vec::<String>::new(),
                [("*28".to_owned(), 2)],
            )],
            None,
        );

        let atazanavir = context
            .drug_report(PrescribingGuidanceSource::DpwgGuideline, "atazanavir")
            .expect("atazanavir DPWG report");

        assert_eq!(atazanavir.name, "atazanavir");
        assert_eq!(atazanavir.source, PrescribingGuidanceSource::DpwgGuideline);
        assert!(atazanavir.is_matched());
        assert_eq!(atazanavir.matched_annotation_count(), 1);

        let guideline = atazanavir.guidelines.iter().next().expect("guideline");
        assert_eq!(
            guideline.name,
            "Annotation of DPWG Guideline for atazanavir and UGT1A1"
        );
        assert_eq!(guideline.genes, ["UGT1A1".to_owned()].into_iter().collect());

        let annotation = guideline.annotations.iter().next().expect("annotation");
        assert_eq!(annotation.local_id, "DPWG-PA166411721");
        assert_eq!(annotation.phenotypes.get("UGT1A1"), None);
        assert_eq!(annotation.activity_scores.get("UGT1A1"), None);
        assert!(
            annotation
                .drug_recommendation
                .as_deref()
                .is_some_and(|text| text.contains("avoid atazanavir"))
        );
    }

    fn gene_message(name: &str, exception_type: &str, matches: MatchLogic) -> MessageAnnotation {
        MessageAnnotation {
            name: name.to_owned(),
            version: None,
            matches,
            exception_type: exception_type.to_owned(),
            message: format!("{name} message"),
        }
    }

    fn report_gene_for_gene_messages() -> ReportGene {
        let diplotype = ReportDiplotype::from_match_label("TEST", "*1/*2", Some(0), None);
        let mut report_gene = ReportGene::new("TEST", ["Normal".to_owned()]);
        report_gene.source_diplotype = Some("*1/*2".to_owned());
        report_gene.source_diplotypes = vec![diplotype.clone()];
        report_gene.recommendation_diplotypes = vec![diplotype];
        report_gene
    }

    fn publication(pmid: Option<&str>, title: Option<&str>, year: Option<i64>) -> Publication {
        Publication {
            pmid: pmid.map(str::to_owned),
            title: title.map(str::to_owned),
            journal: None,
            year,
            same_as: None,
        }
    }

    fn reference_allele_map(definition: &DefinitionFile) -> BTreeMap<String, SampleAlleleSummary> {
        definition
            .variants
            .iter()
            .map(|variant| {
                let call = sample_call(
                    &variant.chromosome,
                    variant.position as usize,
                    Some(&variant.reference),
                    Some(&variant.reference),
                    false,
                    true,
                );
                (variant.vcf_chr_position(), call)
            })
            .collect()
    }

    fn allele_map_from_record_summaries(
        records: impl IntoIterator<Item = SampleAlleleSummary>,
    ) -> BTreeMap<String, SampleAlleleSummary> {
        records
            .into_iter()
            .map(|record| (format!("{}:{}", record.chromosome, record.position), record))
            .collect()
    }

    fn variant_by_rsid<'a>(definition: &'a DefinitionFile, rsid: &str) -> &'a VariantLocus {
        definition
            .variants
            .iter()
            .find(|variant| variant.rsid.as_deref() == Some(rsid))
            .unwrap_or_else(|| panic!("variant {rsid}"))
    }

    fn sample_call(
        chromosome: &str,
        position: usize,
        allele1: Option<&str>,
        allele2: Option<&str>,
        phased: bool,
        effectively_phased: bool,
    ) -> SampleAlleleSummary {
        SampleAlleleSummary {
            chromosome: chromosome.to_owned(),
            position,
            allele1: allele1.map(str::to_owned),
            allele2: allele2.map(str::to_owned),
            vcf_alleles: Vec::new(),
            genotype: if phased { "0|1" } else { "0/1" }.to_owned(),
            vcf_call: format!(
                "{}{}{}",
                allele1.unwrap_or("."),
                if phased { "|" } else { "/" },
                allele2.unwrap_or(".")
            ),
            phased,
            effectively_phased,
            phase_set: None,
            undocumented_variations: BTreeSet::new(),
            treat_undocumented_variations_as_reference: false,
        }
    }
}
