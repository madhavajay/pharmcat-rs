//! Allele definition JSON loading.

use std::{
    cmp::Ordering,
    collections::{BTreeMap, BTreeSet},
    fs, io,
    path::Path,
};

use crate::common::chromosome;
use serde::{Deserialize, Deserializer};

/// One PharmCAT allele definition file.
#[derive(Clone, Debug, Deserialize, PartialEq)]
#[serde(rename_all = "camelCase")]
pub struct DefinitionFile {
    /// Definition format version.
    pub format_version: String,
    /// Data version, when present in newer test fixtures.
    #[serde(default)]
    pub data_version: Option<String>,
    /// Source object, preserved for later reporter parity.
    #[serde(default)]
    pub source: Option<serde_json::Value>,
    /// Source version, when present.
    #[serde(default)]
    pub version: Option<String>,
    /// Modification date string from JSON.
    #[serde(default)]
    pub modification_date: Option<String>,
    /// Gene symbol.
    #[serde(rename = "gene")]
    pub gene_symbol: String,
    /// Gene orientation.
    #[serde(default)]
    pub orientation: Option<String>,
    /// Chromosome name.
    pub chromosome: String,
    /// Genome build.
    pub genome_build: String,
    /// RefSeq chromosome identifier.
    pub ref_seq_chromosome_id: String,
    /// Variant loci in definition order.
    pub variants: Vec<VariantLocus>,
    /// Named alleles in JSON order.
    pub named_alleles: Vec<NamedAllele>,
    /// Positions only used by a single allele that only has one core position.
    #[serde(default)]
    pub singular_variants: BTreeSet<u64>,
    /// Maps VCF positions to named allele names for which the position is core.
    #[serde(default)]
    pub position_to_alleles: BTreeMap<u64, Vec<String>>,
    /// Maps VCF positions to indexes in `variants`.
    #[serde(default)]
    pub position_to_locus: BTreeMap<u64, usize>,
    /// Ambiguous core alleles suppressed in favor of suballeles.
    #[serde(default)]
    pub hidden_core_alleles: Vec<NamedAllele>,
    /// Maps suballele names to parent core allele names.
    #[serde(default)]
    pub suballeles_map: BTreeMap<String, String>,
}

impl DefinitionFile {
    /// Initializes derived definition indexes used by matching.
    pub fn initialize_derived_fields(&mut self) {
        self.position_to_alleles.clear();
        self.position_to_locus.clear();
        self.singular_variants.clear();

        for (index, variant) in self.variants.iter().enumerate() {
            self.position_to_locus.insert(variant.position, index);
        }

        let mut alleles_with_one_position = BTreeSet::new();

        for named_allele in &mut self.named_alleles {
            named_allele.core_positions.clear();

            if named_allele.reference {
                continue;
            }

            for (index, variant) in self.variants.iter().enumerate() {
                if named_allele
                    .cpic_alleles
                    .get(index)
                    .is_some_and(Option::is_some)
                {
                    self.position_to_alleles
                        .entry(variant.position)
                        .or_default()
                        .push(named_allele.name.clone());
                    named_allele.core_positions.insert(variant.position);
                }
            }

            if named_allele.core_positions.len() == 1 {
                alleles_with_one_position.insert(named_allele.name.clone());
            }
        }

        for (position, allele_names) in &self.position_to_alleles {
            if allele_names.len() == 1 && alleles_with_one_position.contains(&allele_names[0]) {
                self.singular_variants.insert(*position);
            }
        }
    }

    /// Returns the reference named allele, matching Java's first reference lookup.
    pub fn reference_named_allele(&self) -> Option<&NamedAllele> {
        self.named_alleles
            .iter()
            .find(|named_allele| named_allele.reference)
    }

    /// Returns a named allele by name.
    pub fn named_allele(&self, name: &str) -> Option<&NamedAllele> {
        self.named_alleles
            .iter()
            .find(|named_allele| named_allele.name == name)
    }

    /// Returns the variant at a VCF position.
    pub fn variant_for_position(&self, position: u64) -> Option<&VariantLocus> {
        self.position_to_locus
            .get(&position)
            .and_then(|index| self.variants.get(*index))
    }

    /// Returns the variant index for a VCF position.
    pub fn index_for_position(&self, position: u64) -> Option<usize> {
        self.position_to_locus.get(&position).copied()
    }

    /// Removes structural named alleles, matching PharmCAT data ingestion.
    pub fn remove_structural_variants(&mut self) {
        self.named_alleles
            .retain(|named_allele| !named_allele.structural_variant);
        self.initialize_derived_fields();
    }

    /// Removes named alleles listed in a definition exemption.
    pub fn remove_ignored_named_alleles(&mut self, exemption: &DefinitionExemption) {
        self.named_alleles
            .retain(|named_allele| !exemption.should_ignore_allele(&named_allele.name));
        self.initialize_derived_fields();
    }

    /// Removes positions listed in a definition exemption and positions no longer used by any
    /// named allele after named-allele filtering.
    pub fn remove_ignored_positions(
        &mut self,
        exemption: &DefinitionExemption,
    ) -> Result<DefinitionPositionRemoval, DefinitionTransformError> {
        let ignored_positions = self
            .variants
            .iter()
            .filter(|variant| exemption.should_ignore_position(variant))
            .cloned()
            .collect::<BTreeSet<_>>();

        if ignored_positions.len() != exemption.ignored_positions.len() {
            return Err(DefinitionTransformError::IgnoredPositionMismatch {
                expected: exemption.ignored_positions.len(),
                found: ignored_positions.len(),
            });
        }

        self.remove_positions(&ignored_positions)
    }

    fn remove_positions(
        &mut self,
        ignored_positions: &BTreeSet<VariantLocus>,
    ) -> Result<DefinitionPositionRemoval, DefinitionTransformError> {
        let mut unused_positions = BTreeSet::new();

        for (index, variant) in self.variants.iter().enumerate() {
            let in_use = self.named_alleles.iter().any(|named_allele| {
                named_allele
                    .cpic_alleles
                    .get(index)
                    .is_some_and(Option::is_some)
            });

            if !in_use && !ignored_positions.contains(variant) {
                unused_positions.insert(variant.clone());
            }
        }

        let mut skipped_indexes = BTreeSet::new();
        let mut new_variants = Vec::new();

        for (index, variant) in self.variants.iter().enumerate() {
            if ignored_positions.contains(variant) || unused_positions.contains(variant) {
                skipped_indexes.insert(index);
            } else {
                new_variants.push(variant.clone());
            }
        }

        let mut updated_named_alleles = Vec::new();
        for named_allele in &self.named_alleles {
            let cpic_alleles = retain_unskipped(&named_allele.cpic_alleles, &skipped_indexes);
            if cpic_alleles.iter().all(Option::is_none) {
                continue;
            }

            let alleles = retain_unskipped(&named_allele.alleles, &skipped_indexes);
            let mut updated = named_allele.clone();
            updated.alleles = alleles;
            updated.cpic_alleles = cpic_alleles;
            updated_named_alleles.push(updated);
        }

        self.variants = new_variants;
        self.named_alleles = updated_named_alleles;
        self.initialize_derived_fields();

        Ok(DefinitionPositionRemoval {
            ignored: ignored_positions.len(),
            unused: unused_positions.len(),
        })
    }
}

fn retain_unskipped<T: Clone>(values: &[T], skipped_indexes: &BTreeSet<usize>) -> Vec<T> {
    values
        .iter()
        .enumerate()
        .filter(|(index, _)| !skipped_indexes.contains(index))
        .map(|(_, value)| value.clone())
        .collect()
}

/// One variant locus from a definition file.
#[derive(Clone, Debug, Deserialize, Eq, PartialEq)]
#[serde(rename_all = "camelCase")]
pub struct VariantLocus {
    /// Chromosome name.
    pub chromosome: String,
    /// VCF position.
    pub position: u64,
    /// CPIC/source position.
    pub cpic_position: u64,
    /// dbSNP identifier.
    #[serde(default)]
    pub rsid: Option<String>,
    /// Chromosomal HGVS name.
    pub chromosome_hgvs_name: String,
    /// CPIC alleles.
    pub cpic_alleles: BTreeSet<String>,
    /// CPIC to VCF allele mapping.
    pub cpic_to_vcf_allele_map: BTreeMap<String, String>,
    /// VCF reference allele.
    #[serde(rename = "ref")]
    pub reference: String,
    /// VCF alternate alleles.
    #[serde(default)]
    pub alts: Vec<String>,
}

impl VariantLocus {
    /// Returns `<chromosome>:<position>`.
    pub fn vcf_chr_position(&self) -> String {
        format!("{}:{}", self.chromosome, self.position)
    }

    /// Returns whether `allele` is present in the CPIC-to-VCF allele map.
    pub fn has_vcf_allele(&self, allele: &str) -> bool {
        self.cpic_to_vcf_allele_map
            .values()
            .any(|vcf_allele| vcf_allele == allele)
    }

    /// Returns a simple HGVS-like name for a VCF allele, matching Java `VariantLocus#getHgvsForVcfAllele`.
    pub fn hgvs_for_vcf_allele(&self, vcf_allele: &str) -> String {
        if vcf_allele == "." {
            return format!("g.{}?", self.position);
        }
        if vcf_allele == self.reference {
            return format!("g.{}=", self.position);
        }

        for (index, alt) in self.alts.iter().enumerate() {
            if vcf_allele == alt {
                return self
                    .chromosome_hgvs_name
                    .split(';')
                    .nth(index)
                    .unwrap_or(&self.chromosome_hgvs_name)
                    .to_owned();
            }
        }

        format!("g.{}{}>{}", self.position, self.reference, vcf_allele)
    }
}

impl Ord for VariantLocus {
    fn cmp(&self, other: &Self) -> Ordering {
        chromosome::compare_names(Some(&self.chromosome), Some(&other.chromosome))
            .then_with(|| self.position.cmp(&other.position))
            .then_with(|| self.chromosome_hgvs_name.cmp(&other.chromosome_hgvs_name))
            .then_with(|| self.rsid.cmp(&other.rsid))
    }
}

impl PartialOrd for VariantLocus {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// One named allele from a definition file.
#[derive(Clone, Debug, Deserialize, PartialEq)]
#[serde(rename_all = "camelCase")]
pub struct NamedAllele {
    /// Stable allele ID.
    pub id: String,
    /// Star allele or other allele name.
    pub name: String,
    /// VCF alleles by variant index.
    #[serde(default)]
    pub alleles: Vec<Option<String>>,
    /// CPIC alleles by variant index.
    #[serde(default)]
    pub cpic_alleles: Vec<Option<String>>,
    /// Population frequency from source data, if present.
    #[serde(default)]
    pub population_frequency: Option<serde_json::Value>,
    /// Whether this is the reference named allele.
    #[serde(default)]
    pub reference: bool,
    /// Whether this JSON marks the allele as a combination or partial.
    #[serde(default)]
    pub is_combination_or_partial: bool,
    /// Whether this allele is a structural variant removed during data ingestion.
    #[serde(default)]
    pub structural_variant: bool,
    /// Java matcher score initialized from the definition file.
    #[serde(default)]
    pub score: Option<i32>,
    /// Core VCF positions defining this named allele.
    #[serde(default, deserialize_with = "null_to_default")]
    pub core_positions: BTreeSet<u64>,
    /// Definition positions missing from the sample after match-data marshalling.
    #[serde(default, skip)]
    pub missing_positions: BTreeSet<VariantLocus>,
    /// Matcher score override used after missing-position marshalling.
    #[serde(default, skip)]
    pub score_override: Option<i32>,
    /// Combination count used by matcher ranking.
    #[serde(default)]
    pub num_combinations: i32,
    /// Partial count used by matcher ranking.
    #[serde(default)]
    pub num_partials: i32,
}

impl NamedAllele {
    /// Java `NamedAllele.isCombination`.
    pub fn is_combination(&self) -> bool {
        self.num_combinations > 1
    }

    /// Java `NamedAllele.isPartial`.
    pub fn is_partial(&self) -> bool {
        self.num_partials > 0
    }
}

fn null_to_default<'de, D, T>(deserializer: D) -> Result<T, D::Error>
where
    D: Deserializer<'de>,
    T: Default + Deserialize<'de>,
{
    Ok(Option::<T>::deserialize(deserializer)?.unwrap_or_default())
}

/// Special handling applied to one gene's definition.
#[derive(Clone, Debug, Deserialize, PartialEq)]
#[serde(rename_all = "camelCase")]
pub struct DefinitionExemption {
    /// Gene symbol.
    pub gene: String,
    /// Positions required before making a call.
    #[serde(default)]
    pub required_positions: BTreeSet<u64>,
    /// Definition positions ignored during data preparation.
    #[serde(default)]
    pub ignored_positions: BTreeSet<VariantLocus>,
    /// Additional positions of interest not in the definition variants.
    #[serde(default)]
    pub extra_positions: BTreeSet<VariantLocus>,
    /// Ignored named alleles.
    #[serde(default)]
    pub ignored_alleles: BTreeSet<String>,
    /// Lowercase ignored named alleles as serialized by Java.
    #[serde(default)]
    pub ignored_alleles_lc: BTreeSet<String>,
    /// Unphased diplotype priority overrides.
    #[serde(default)]
    pub unphased_diplotype_priorities: BTreeSet<UnphasedDiplotypePriority>,
    /// AMP1 alleles.
    #[serde(default)]
    pub amp1_alleles: Vec<String>,
    /// AMP1 positions.
    #[serde(default)]
    pub amp1_positions: BTreeSet<u64>,
}

impl DefinitionExemption {
    /// Returns whether `allele` should be ignored, case-insensitively.
    pub fn should_ignore_allele(&self, allele: &str) -> bool {
        if self.ignored_alleles_lc.is_empty() {
            self.ignored_alleles
                .iter()
                .any(|ignored| ignored.eq_ignore_ascii_case(allele))
        } else {
            self.ignored_alleles_lc.contains(&allele.to_lowercase())
        }
    }

    /// Returns whether `position` should be ignored, matching Java's RSID/HGVS logic.
    pub fn should_ignore_position(&self, position: &VariantLocus) -> bool {
        self.ignored_positions.iter().any(|ignored| {
            if let Some(rsid) = &ignored.rsid {
                position.rsid.as_ref() == Some(rsid)
            } else {
                ignored.chromosome_hgvs_name == position.chromosome_hgvs_name
            }
        })
    }
}

/// Unphased diplotype priority override.
#[derive(Clone, Debug, Deserialize, Eq, Ord, PartialEq, PartialOrd)]
#[serde(rename_all = "camelCase")]
pub struct UnphasedDiplotypePriority {
    /// Serialized priority key.
    pub id: String,
    /// Candidate diplotypes.
    pub list: BTreeSet<String>,
    /// Diplotype to pick.
    pub pick: String,
}

/// Counts from definition position removal.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct DefinitionPositionRemoval {
    /// Explicitly ignored positions removed.
    pub ignored: usize,
    /// Unused positions removed after ignored named alleles were removed.
    pub unused: usize,
}

/// Definition ingestion transform error.
#[derive(Debug, Eq, PartialEq)]
pub enum DefinitionTransformError {
    /// An exemption referenced positions that were not found in the definition.
    IgnoredPositionMismatch {
        /// Number of ignored positions in the exemption.
        expected: usize,
        /// Number of ignored positions found in the definition.
        found: usize,
    },
}

impl std::fmt::Display for DefinitionTransformError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::IgnoredPositionMismatch { expected, found } => write!(
                f,
                "Should have {expected} ignored positions, but only found {found}"
            ),
        }
    }
}

impl std::error::Error for DefinitionTransformError {}

/// Multiple definition files indexed like Java `DefinitionReader`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct DefinitionReader {
    definitions: BTreeMap<String, DefinitionFile>,
    exemptions: BTreeMap<String, DefinitionExemption>,
    locations_of_interest: BTreeMap<String, VariantLocus>,
    locations_by_gene: BTreeMap<String, String>,
}

impl DefinitionReader {
    /// Loads definition files.
    pub fn from_paths<I, P>(paths: I) -> Result<Self, DefinitionLoadError>
    where
        I: IntoIterator<Item = P>,
        P: AsRef<Path>,
    {
        let mut definitions = BTreeMap::new();

        for path in paths {
            let definition = read_definition_file(path.as_ref())?;
            definitions.insert(definition.gene_symbol.clone(), definition);
        }

        Ok(Self::from_definitions(definitions))
    }

    /// Loads definition files and exemptions.
    pub fn from_paths_with_exemptions<I, P, E>(
        paths: I,
        exemptions_path: E,
    ) -> Result<Self, DefinitionLoadError>
    where
        I: IntoIterator<Item = P>,
        P: AsRef<Path>,
        E: AsRef<Path>,
    {
        let mut reader = Self::from_paths(paths)?;
        reader.exemptions = read_exemptions_file(exemptions_path.as_ref())?;
        reader.rebuild_location_indexes();
        Ok(reader)
    }

    /// Builds a reader from already-loaded definitions.
    pub fn from_definitions(definitions: BTreeMap<String, DefinitionFile>) -> Self {
        Self::from_definitions_and_exemptions(definitions, BTreeMap::new())
    }

    /// Builds a reader from already-loaded definitions and exemptions.
    pub fn from_definitions_and_exemptions(
        definitions: BTreeMap<String, DefinitionFile>,
        exemptions: BTreeMap<String, DefinitionExemption>,
    ) -> Self {
        let mut reader = Self {
            definitions,
            exemptions,
            locations_of_interest: BTreeMap::new(),
            locations_by_gene: BTreeMap::new(),
        };
        reader.rebuild_location_indexes();
        reader
    }

    fn rebuild_location_indexes(&mut self) {
        let mut locations_of_interest = BTreeMap::new();
        let mut locations_by_gene = BTreeMap::new();

        for (gene, definition) in &self.definitions {
            for variant in &definition.variants {
                let chr_position = variant.vcf_chr_position();
                locations_of_interest.insert(chr_position.clone(), variant.clone());
                locations_by_gene.insert(chr_position, gene.clone());
            }

            if let Some(exemption) = self.exemptions.get(gene) {
                for variant in &exemption.extra_positions {
                    let chr_position = variant.vcf_chr_position();
                    if !locations_of_interest.contains_key(&chr_position) {
                        locations_of_interest.insert(chr_position.clone(), variant.clone());
                        locations_by_gene.insert(chr_position, gene.clone());
                    }
                }
            }
        }

        self.locations_of_interest = locations_of_interest;
        self.locations_by_gene = locations_by_gene;
    }

    /// Returns sorted gene names.
    pub fn genes(&self) -> impl Iterator<Item = &str> {
        self.definitions.keys().map(String::as_str)
    }

    /// Returns a definition file for `gene`.
    pub fn definition_file(&self, gene: &str) -> Option<&DefinitionFile> {
        self.definitions.get(gene)
    }

    /// Returns an exemption for `gene`.
    pub fn exemption(&self, gene: &str) -> Option<&DefinitionExemption> {
        self.exemptions.get(gene)
    }

    /// Returns the common genome build, or an error if files disagree.
    pub fn genome_build(&self) -> Result<&str, DefinitionLoadError> {
        let mut builds = self.definitions.values().map(|definition| {
            (
                definition.gene_symbol.as_str(),
                definition.genome_build.as_str(),
            )
        });

        let Some((_, first_build)) = builds.next() else {
            return Err(DefinitionLoadError::NoDefinitions);
        };

        for (gene, build) in builds {
            if !first_build.eq_ignore_ascii_case(build) {
                return Err(DefinitionLoadError::MixedGenomeBuilds {
                    first: first_build.to_owned(),
                    second: build.to_owned(),
                    gene: gene.to_owned(),
                });
            }
        }

        Ok(first_build)
    }

    /// Returns `<chr:position>` to locus mappings.
    pub fn locations_of_interest(&self) -> &BTreeMap<String, VariantLocus> {
        &self.locations_of_interest
    }

    /// Returns `<chr:position>` to gene mappings.
    pub fn locations_by_gene(&self) -> &BTreeMap<String, String> {
        &self.locations_by_gene
    }
}

/// Reads one definition file.
pub fn read_definition_file(path: &Path) -> Result<DefinitionFile, DefinitionLoadError> {
    // `serde_json::from_reader` iterates `Read::bytes`; using it with a bare `File` makes large
    // PharmCAT resources syscall-bound. Read in bulk and let the slice reader parse in memory.
    let data = fs::read(path)?;
    let mut definition: DefinitionFile =
        serde_json::from_slice(&data).map_err(DefinitionLoadError::Parse)?;
    definition.initialize_derived_fields();
    Ok(definition)
}

/// Reads Java PharmCAT exemptions JSON.
pub fn read_exemptions_file(
    path: &Path,
) -> Result<BTreeMap<String, DefinitionExemption>, DefinitionLoadError> {
    let data = fs::read(path)?;
    serde_json::from_slice(&data).map_err(DefinitionLoadError::Parse)
}

/// Definition loading error.
#[derive(Debug)]
pub enum DefinitionLoadError {
    /// File I/O failed.
    Io(io::Error),
    /// JSON parsing failed.
    Parse(serde_json::Error),
    /// No definition files were loaded.
    NoDefinitions,
    /// Definition files disagree on genome build.
    MixedGenomeBuilds {
        /// First genome build encountered.
        first: String,
        /// Conflicting genome build.
        second: String,
        /// Gene carrying the conflicting build.
        gene: String,
    },
}

impl std::fmt::Display for DefinitionLoadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io(err) => err.fmt(f),
            Self::Parse(err) => write!(f, "failed to parse definition JSON: {err}"),
            Self::NoDefinitions => write!(f, "no definition files loaded"),
            Self::MixedGenomeBuilds {
                first,
                second,
                gene,
            } => write!(
                f,
                "Definition files use different genome builds ({first} vs {second} for {gene})"
            ),
        }
    }
}

impl std::error::Error for DefinitionLoadError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(err) => Some(err),
            Self::Parse(err) => Some(err),
            Self::NoDefinitions | Self::MixedGenomeBuilds { .. } => None,
        }
    }
}

impl From<io::Error> for DefinitionLoadError {
    fn from(err: io::Error) -> Self {
        Self::Io(err)
    }
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use super::{
        DefinitionExemption, DefinitionFile, DefinitionLoadError, DefinitionReader,
        DefinitionTransformError, NamedAllele, VariantLocus, read_definition_file,
        read_exemptions_file,
    };

    #[test]
    fn reads_definition_file_metadata_variants_and_named_alleles() {
        let definition = read_definition_file(Path::new(
            "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-cyp2b6.json",
        ))
        .expect("definition JSON");

        assert_eq!(definition.format_version, "2");
        assert_eq!(definition.gene_symbol, "CYP2B6");
        assert_eq!(definition.chromosome, "chr19");
        assert_eq!(definition.genome_build, "GRCh38.p13");
        assert_eq!(definition.variants.len(), 35);
        assert_eq!(definition.named_alleles.len(), 35);

        let first_variant = &definition.variants[0];
        assert_eq!(first_variant.vcf_chr_position(), "chr19:40991224");
        assert_eq!(first_variant.rsid.as_deref(), Some("rs34223104"));
        assert_eq!(first_variant.reference, "T");
        assert_eq!(first_variant.alts, ["C"]);
        assert!(first_variant.has_vcf_allele("C"));

        let reference = definition
            .reference_named_allele()
            .expect("reference named allele");
        assert_eq!(reference.name, "*1");
        assert_eq!(reference.alleles.len(), definition.variants.len());
        assert!(reference.reference);

        assert_eq!(
            definition.index_for_position(40991224).expect("position"),
            0
        );
        assert_eq!(
            definition
                .variant_for_position(40991224)
                .expect("variant")
                .rsid
                .as_deref(),
            Some("rs34223104")
        );
        assert_eq!(
            definition
                .named_allele("*2")
                .expect("named allele")
                .core_positions
                .iter()
                .copied()
                .collect::<Vec<_>>(),
            [40991369]
        );
        assert_eq!(
            definition
                .named_allele("*6")
                .expect("named allele")
                .core_positions
                .iter()
                .copied()
                .collect::<Vec<_>>(),
            [41006936, 41009358]
        );
        assert_eq!(definition.position_to_alleles.len(), 35);
        assert_eq!(
            definition
                .position_to_alleles
                .get(&40991369)
                .expect("alleles"),
            &vec!["*2".to_owned(), "*10".to_owned()]
        );
        assert_eq!(
            definition
                .singular_variants
                .iter()
                .copied()
                .collect::<Vec<_>>(),
            [
                40991441, 41004125, 41004381, 41009350, 41010108, 41012693, 41012740, 41012803,
                41016726, 41016778, 41016805
            ]
        );
    }

    #[test]
    fn indexes_multiple_definitions_like_java_definition_reader() {
        let reader = DefinitionReader::from_paths([
            "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-cyp2b6.json",
            "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-dpyd.json",
        ])
        .expect("definitions");

        assert_eq!(reader.genes().collect::<Vec<_>>(), ["CYP2B6", "DPYD"]);
        assert_eq!(reader.genome_build().expect("build"), "GRCh38.p13");
        assert_eq!(
            reader
                .definition_file("DPYD")
                .expect("DPYD definition")
                .reference_named_allele()
                .expect("reference")
                .name,
            "Reference"
        );
        assert_eq!(
            reader
                .locations_by_gene()
                .get("chr19:40991224")
                .map(String::as_str),
            Some("CYP2B6")
        );
        assert!(reader.locations_of_interest().contains_key("chr1:97078987"));
    }

    #[test]
    fn rejects_empty_definition_reader_for_genome_build() {
        let reader = DefinitionReader::default();
        assert!(matches!(
            reader.genome_build(),
            Err(DefinitionLoadError::NoDefinitions)
        ));
    }

    #[test]
    fn reads_exemptions_json_like_java_dataserializer() {
        let exemptions = read_exemptions_file(Path::new(
            "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/definition/DataSerializerTest-exemptions.json",
        ))
        .expect("exemptions JSON");

        assert_eq!(exemptions.len(), 4);
        assert!(!exemptions.contains_key("TPMT"));

        let cyp2c9 = exemptions.get("CYP2C9").expect("CYP2C9 exemption");
        assert_eq!(cyp2c9.required_positions.len(), 0);
        assert_eq!(cyp2c9.ignored_positions.len(), 0);
        assert_eq!(cyp2c9.extra_positions.len(), 1);
        assert_eq!(cyp2c9.ignored_alleles.len(), 0);
        let extra = cyp2c9.extra_positions.iter().next().expect("extra locus");
        assert_eq!(extra.vcf_chr_position(), "chr10:94645745");
        assert_eq!(extra.rsid.as_deref(), Some("rs12777823"));

        let g6pd = exemptions.get("G6PD").expect("G6PD exemption");
        assert_eq!(g6pd.required_positions.len(), 0);
        assert_eq!(g6pd.ignored_positions.len(), 2);
        assert_eq!(g6pd.extra_positions.len(), 0);
        assert_eq!(g6pd.ignored_alleles.len(), 1);
        assert!(g6pd.should_ignore_allele("mediterranean haplotype"));
        assert!(g6pd.should_ignore_allele("Mediterranean Haplotype"));
        assert!(g6pd.should_ignore_position(&variant_probe(
            "chrX",
            154532439,
            Some("rs2230037"),
            "different-hgvs",
        )));
        assert!(g6pd.should_ignore_position(&variant_probe(
            "chrX",
            154532990,
            None,
            "g.154532991_154532993delGGT",
        )));
        assert!(!g6pd.should_ignore_position(&variant_probe(
            "chrX",
            154532439,
            None,
            "g.154532439A>G",
        )));

        let nat2 = exemptions.get("NAT2").expect("NAT2 exemption");
        assert_eq!(
            nat2.required_positions.iter().copied().collect::<Vec<_>>(),
            [18400344, 18400593, 18400806, 18400860]
        );
        assert_eq!(nat2.ignored_positions.len(), 0);
        assert_eq!(nat2.extra_positions.len(), 0);
        assert_eq!(nat2.ignored_alleles.len(), 0);

        let ryr1 = exemptions.get("RYR1").expect("RYR1 exemption");
        assert_eq!(ryr1.required_positions.len(), 0);
        assert_eq!(ryr1.ignored_positions.len(), 1);
    }

    #[test]
    fn reads_unphased_diplotype_priorities_from_exemptions() {
        let exemptions = read_exemptions_file(Path::new(
            "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/definition/DataSerializerTest-exemptionsWithUnphasedPriorities.json",
        ))
        .expect("exemptions JSON");

        let nat2 = exemptions.get("NAT2").expect("NAT2 exemption");
        assert_eq!(nat2.unphased_diplotype_priorities.len(), 4);
        let priority = nat2
            .unphased_diplotype_priorities
            .iter()
            .find(|priority| priority.id == "*1/*15|*14/*34|*6/*46")
            .expect("priority");
        assert_eq!(priority.pick, "*1/*15");
        assert_eq!(
            priority.list.iter().map(String::as_str).collect::<Vec<_>>(),
            ["*1/*15", "*14/*34", "*6/*46"]
        );
    }

    #[test]
    fn definition_reader_indexes_extra_positions_from_exemptions() {
        let reader = DefinitionReader::from_paths_with_exemptions(
            ["../../repos/PharmCAT/src/main/resources/org/pharmgkb/pharmcat/definition/alleles/CYP2C9_translation.json"],
            "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/definition/DataSerializerTest-exemptions.json",
        )
        .expect("definition reader");

        assert_eq!(
            reader
                .locations_by_gene()
                .get("chr10:94645745")
                .map(String::as_str),
            Some("CYP2C9")
        );
        assert_eq!(
            reader
                .locations_of_interest()
                .get("chr10:94645745")
                .and_then(|locus| locus.rsid.as_deref()),
            Some("rs12777823")
        );
    }

    #[test]
    fn removes_structural_variants_like_java_ingestion() {
        let mut definition = synthetic_definition();

        definition.remove_structural_variants();

        assert!(definition.named_allele("*structural").is_none());
        assert_eq!(
            definition
                .named_alleles
                .iter()
                .map(|allele| allele.name.as_str())
                .collect::<Vec<_>>(),
            ["*1", "*ignore", "*keep"]
        );
        assert_eq!(definition.position_to_alleles.len(), 2);
    }

    #[test]
    fn removes_ignored_alleles_positions_and_unused_positions_like_java_ingestion() {
        let mut definition = synthetic_definition();
        let exemption = synthetic_exemption();

        definition.remove_ignored_named_alleles(&exemption);
        assert!(definition.named_allele("*ignore").is_none());

        let removal = definition
            .remove_ignored_positions(&exemption)
            .expect("remove positions");

        assert_eq!(removal.ignored, 1);
        assert_eq!(removal.unused, 1);
        assert_eq!(
            definition
                .variants
                .iter()
                .map(|variant| variant.position)
                .collect::<Vec<_>>(),
            [200]
        );
        assert_eq!(
            definition
                .named_allele("*keep")
                .expect("keep allele")
                .cpic_alleles,
            [Some("T".to_owned())]
        );
        assert_eq!(
            definition.position_to_alleles.get(&200).expect("alleles"),
            &vec!["*keep".to_owned()]
        );
    }

    #[test]
    fn errors_when_exemption_ignored_positions_are_not_found_like_java() {
        let mut definition = synthetic_definition();
        let mut exemption = synthetic_exemption();
        exemption.ignored_positions.insert(variant_probe(
            "chr1",
            999,
            Some("rs-missing"),
            "g.999A>G",
        ));

        assert_eq!(
            definition.remove_ignored_positions(&exemption),
            Err(DefinitionTransformError::IgnoredPositionMismatch {
                expected: 2,
                found: 1
            })
        );
    }

    fn synthetic_definition() -> DefinitionFile {
        let variants = vec![
            variant_probe("chr1", 100, Some("rs-ignore"), "g.100A>G"),
            variant_probe("chr1", 200, Some("rs-keep"), "g.200A>T"),
            variant_probe("chr1", 300, None, "g.300A>C"),
        ];

        let mut definition = DefinitionFile {
            format_version: "2".to_owned(),
            data_version: None,
            source: None,
            version: None,
            modification_date: None,
            gene_symbol: "SYN".to_owned(),
            orientation: None,
            chromosome: "chr1".to_owned(),
            genome_build: "GRCh38.p13".to_owned(),
            ref_seq_chromosome_id: "NC_000001.11".to_owned(),
            variants,
            named_alleles: vec![
                synthetic_allele("*1", false, [Some("A"), None, None]),
                synthetic_allele("*ignore", false, [None, None, Some("C")]),
                synthetic_allele("*keep", false, [None, Some("T"), None]),
                synthetic_allele("*structural", true, [Some("G"), None, None]),
            ],
            singular_variants: Default::default(),
            position_to_alleles: Default::default(),
            position_to_locus: Default::default(),
            hidden_core_alleles: Vec::new(),
            suballeles_map: Default::default(),
        };
        definition.initialize_derived_fields();
        definition
    }

    fn synthetic_allele(
        name: &str,
        structural_variant: bool,
        cpic_alleles: [Option<&str>; 3],
    ) -> NamedAllele {
        let cpic_alleles = cpic_alleles
            .into_iter()
            .map(|allele| allele.map(str::to_owned))
            .collect::<Vec<_>>();
        NamedAllele {
            id: format!("SYN{name}"),
            name: name.to_owned(),
            alleles: cpic_alleles.clone(),
            cpic_alleles,
            population_frequency: None,
            reference: name == "*1",
            is_combination_or_partial: false,
            structural_variant,
            score: None,
            core_positions: Default::default(),
            missing_positions: Default::default(),
            score_override: None,
            num_combinations: 0,
            num_partials: 0,
        }
    }

    fn synthetic_exemption() -> DefinitionExemption {
        DefinitionExemption {
            gene: "SYN".to_owned(),
            required_positions: Default::default(),
            ignored_positions: [variant_probe(
                "chr1",
                100,
                Some("rs-ignore"),
                "ignored-hgvs",
            )]
            .into_iter()
            .collect(),
            extra_positions: Default::default(),
            ignored_alleles: ["*ignore".to_owned()].into_iter().collect(),
            ignored_alleles_lc: ["*ignore".to_owned()].into_iter().collect(),
            unphased_diplotype_priorities: Default::default(),
            amp1_alleles: Vec::new(),
            amp1_positions: Default::default(),
        }
    }

    fn variant_probe(
        chromosome: &str,
        position: u64,
        rsid: Option<&str>,
        chromosome_hgvs_name: &str,
    ) -> VariantLocus {
        VariantLocus {
            chromosome: chromosome.to_owned(),
            position,
            cpic_position: position,
            rsid: rsid.map(str::to_owned),
            chromosome_hgvs_name: chromosome_hgvs_name.to_owned(),
            cpic_alleles: Default::default(),
            cpic_to_vcf_allele_map: Default::default(),
            reference: String::new(),
            alts: Vec::new(),
        }
    }
}
