//! Phenotype and activity-score helpers.

use std::{
    collections::{BTreeMap, BTreeSet},
    fs, io,
    path::Path,
};

use serde::Deserialize;

use crate::matcher::compare_haplotype_names;

/// Java `TextConstants.NA`.
pub const NA: &str = "n/a";

/// Java `TextConstants.GTE`.
pub const GTE: &str = "\u{2265}";

/// Java `TextConstants.INDETERMINATE`.
pub const INDETERMINATE: &str = "Indeterminate";

/// Java `TextConstants.NO_RESULT`.
pub const NO_RESULT: &str = "No Result";

/// Loaded phenotype mappings keyed by gene.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PhenotypeMap {
    genes: BTreeMap<String, GenePhenotype>,
}

impl PhenotypeMap {
    /// Loads all Java phenotype JSON files from `dir`.
    pub fn from_dir(dir: &Path) -> Result<Self, PhenotypeLoadError> {
        if !dir.is_dir() {
            return Err(PhenotypeLoadError::NotDirectory(dir.to_path_buf()));
        }

        let mut json_files = fs::read_dir(dir)?
            .map(|entry| entry.map(|entry| entry.path()))
            .collect::<Result<Vec<_>, _>>()?;
        json_files.retain(|path| {
            path.extension()
                .is_some_and(|extension| extension == "json")
        });
        json_files.sort();
        if json_files.is_empty() {
            return Err(PhenotypeLoadError::NoPhenotypeFiles(dir.to_path_buf()));
        }

        let mut genes = BTreeMap::new();
        for path in json_files {
            let phenotype = read_gene_phenotype_file(&path)?;
            let gene = phenotype.gene.clone();
            if genes.insert(gene.clone(), phenotype).is_some() {
                return Err(PhenotypeLoadError::DuplicateGene(gene));
            }
        }

        Ok(Self { genes })
    }

    /// Returns sorted gene names.
    pub fn genes(&self) -> impl Iterator<Item = &str> {
        self.genes.keys().map(String::as_str)
    }

    /// Returns a phenotype map for `gene`.
    pub fn phenotype(&self, gene: &str) -> Option<&GenePhenotype> {
        self.genes.get(gene)
    }

    /// Returns the version for `gene`, if known.
    pub fn version(&self, gene: &str) -> Option<&str> {
        self.phenotype(gene)
            .and_then(|phenotype| phenotype.version.as_deref())
    }
}

/// One Java `GenePhenotype` JSON file.
#[derive(Clone, Debug, Deserialize, PartialEq)]
#[serde(rename_all = "camelCase")]
pub struct GenePhenotype {
    /// HGNC gene symbol.
    pub gene: String,
    /// Haplotype-to-function or haplotype-to-activity lookup map.
    #[serde(default)]
    pub haplotypes: BTreeMap<String, String>,
    /// Haplotype activity values for activity-score genes.
    #[serde(default)]
    pub activity_values: BTreeMap<String, String>,
    /// Diplotype phenotype records.
    #[serde(default)]
    pub diplotypes: Vec<DiplotypeRecord>,
    /// Haplotype records.
    #[serde(default)]
    pub named_alleles: Vec<HaplotypeRecord>,
    /// Phenotype data version.
    #[serde(default)]
    pub version: Option<String>,
}

impl GenePhenotype {
    /// Java `GenePhenotype.UNASSIGNED_FUNCTION`.
    pub const UNASSIGNED_FUNCTION: &'static str = "Unassigned function";

    /// Returns whether this gene uses activity scores.
    pub fn is_activity_gene(&self) -> bool {
        !self.activity_values.is_empty()
    }

    /// Looks up a haplotype function like Java `getHaplotypeFunction`.
    pub fn haplotype_function(&self, haplotype: &str) -> &str {
        if haplotype.trim().is_empty() {
            return Self::UNASSIGNED_FUNCTION;
        }

        self.named_alleles
            .iter()
            .find(|record| record.name == haplotype)
            .and_then(|record| record.function_value.as_deref())
            .or_else(|| self.haplotypes.get(haplotype).map(String::as_str))
            .unwrap_or(Self::UNASSIGNED_FUNCTION)
    }

    /// Looks up a haplotype activity value.
    pub fn haplotype_activity(&self, haplotype: &str) -> Option<&str> {
        self.activity_values.get(haplotype).map(String::as_str)
    }

    /// Java `GenePhenotype.makeFormattedFunctionScoreMap`.
    pub fn formatted_function_score_map(&self) -> BTreeMap<String, String> {
        self.named_alleles
            .iter()
            .map(|record| {
                let function = if is_unspecified(record.activity_value.as_deref()) {
                    record.function_value.clone().unwrap_or_default()
                } else {
                    format!(
                        "Activity Value {} ({})",
                        record.activity_value.as_deref().unwrap_or_default(),
                        record.function_value.as_deref().unwrap_or_default()
                    )
                };
                (record.name.clone(), function)
            })
            .collect()
    }

    /// Finds a diplotype record by diplotype key.
    pub fn find_diplotype(
        &self,
        diplotype_key: &BTreeMap<String, i32>,
    ) -> Option<&DiplotypeRecord> {
        if diplotype_key.is_empty() {
            return None;
        }
        self.diplotypes
            .iter()
            .find(|record| record.matches_key(diplotype_key))
    }

    /// Looks up phenotype by diplotype key like Java `lookupPhenotypesByDiplotype`.
    pub fn lookup_phenotype_by_diplotype_key(
        &self,
        diplotype_key: &BTreeMap<String, i32>,
        has_unknown_alleles: bool,
        is_activity_score_type: bool,
    ) -> Result<String, PhenotypeLookupError> {
        if has_unknown_alleles {
            return Ok(NO_RESULT.to_owned());
        }

        let keys = self
            .diplotypes
            .iter()
            .filter(|record| record.matches_key(diplotype_key))
            .map(|record| record.gene_result.clone())
            .collect::<BTreeSet<_>>();

        unique_lookup(
            "phenotype",
            keys,
            if is_activity_score_type {
                INDETERMINATE
            } else {
                NA
            },
        )
    }

    /// Looks up activity score by diplotype key like Java `lookupActivityByDiplotype`.
    pub fn lookup_activity_by_diplotype_key(
        &self,
        diplotype_key: &BTreeMap<String, i32>,
        has_unknown_alleles: bool,
    ) -> Result<String, PhenotypeLookupError> {
        if has_unknown_alleles {
            return Ok(NO_RESULT.to_owned());
        }

        let keys = self
            .diplotypes
            .iter()
            .filter(|record| record.matches_key(diplotype_key))
            .filter_map(|record| record.activity_score.clone())
            .collect::<BTreeSet<_>>();

        unique_lookup("activity score", keys, NA)
    }

    /// Looks up phenotypes by activity score like Java `lookupPhenotypeByActivityScore`.
    pub fn lookup_phenotypes_by_activity_score(&self, activity_score: &str) -> BTreeSet<String> {
        let mut phenotypes = self
            .diplotypes
            .iter()
            .filter(|record| record.lookup_key == activity_score)
            .map(|record| record.gene_result.clone())
            .collect::<BTreeSet<_>>();

        if phenotypes.is_empty() {
            phenotypes.insert(INDETERMINATE.to_owned());
        }

        phenotypes
    }

    /// Looks up activity scores by phenotype like Java `lookupActivityScoresByPhenotype`.
    pub fn lookup_activity_scores_by_phenotype<'a>(
        &self,
        phenotypes: impl IntoIterator<Item = &'a str>,
    ) -> BTreeSet<String> {
        let phenotype_set = phenotypes.into_iter().collect::<BTreeSet<_>>();
        let mut scores = self
            .diplotypes
            .iter()
            .filter(|record| phenotype_set.contains(record.gene_result.as_str()))
            .filter_map(|record| record.activity_score.clone())
            .collect::<BTreeSet<_>>();

        if scores.is_empty() {
            scores.insert(NA.to_owned());
        }

        scores
    }

    fn lookup_key_by_diplotype_key(
        &self,
        key_type: &'static str,
        diplotype_key: &BTreeMap<String, i32>,
    ) -> Result<String, PhenotypeLookupError> {
        let keys = self
            .diplotypes
            .iter()
            .filter(|record| record.matches_key(diplotype_key))
            .map(|record| record.lookup_key.clone())
            .collect::<BTreeSet<_>>();

        unique_lookup(key_type, keys, NA)
    }

    /// Infers the recommendation diplotype from DPYD matcher diplotype labels like Java
    /// `LowestFunctionGeneCaller.inferFromDiplotypes`.
    pub fn infer_dpyd_lowest_function_from_diplotypes(
        &self,
        diplotypes: impl IntoIterator<Item = impl AsRef<str>>,
    ) -> Option<DiplotypeAnnotationInput> {
        let mut candidates = diplotypes
            .into_iter()
            .filter_map(|diplotype| {
                let (hap1, hap2) = split_diplotype_label(diplotype.as_ref());
                let allele1 = self.best_dpyd_haplotype(split_haplotype_label(&hap1))?;
                let allele2 = self.best_dpyd_haplotype(split_haplotype_label(hap2.as_deref()?))?;
                Some((allele1, allele2))
            })
            .collect::<Vec<_>>();

        candidates.sort_by(|left, right| self.compare_dpyd_diplotype_pairs(left, right));
        candidates.into_iter().next().map(|(allele1, allele2)| {
            DiplotypeAnnotationInput::from_alleles(&self.gene, allele1, Some(allele2))
        })
    }

    /// Infers the recommendation diplotype from unphased DPYD haplotype labels like Java
    /// `LowestFunctionGeneCaller.inferFromHaplotypeMatches`.
    pub fn infer_dpyd_lowest_function_from_haplotypes(
        &self,
        haplotypes: impl IntoIterator<Item = impl AsRef<str>>,
    ) -> Option<DiplotypeAnnotationInput> {
        let mut haplotypes = haplotypes
            .into_iter()
            .flat_map(|haplotype| split_haplotype_label(haplotype.as_ref()))
            .collect::<Vec<_>>();
        haplotypes.sort_by(|left, right| self.compare_dpyd_haplotypes(left, right));

        let allele1 = haplotypes.first()?.clone();
        let allele2 = haplotypes.get(1).cloned();
        Some(DiplotypeAnnotationInput::from_alleles(
            &self.gene, allele1, allele2,
        ))
    }

    fn best_dpyd_haplotype(&self, haplotypes: Vec<String>) -> Option<String> {
        haplotypes
            .into_iter()
            .min_by(|left, right| self.compare_dpyd_haplotypes(left, right))
    }

    fn compare_dpyd_haplotypes(&self, left: &str, right: &str) -> std::cmp::Ordering {
        match (
            self.dpyd_haplotype_activity(left),
            self.dpyd_haplotype_activity(right),
        ) {
            (None, None) => compare_haplotype_names(left, right),
            (None, Some(_)) => std::cmp::Ordering::Less,
            (Some(_), None) => std::cmp::Ordering::Greater,
            (Some(left_activity), Some(right_activity)) => left_activity
                .total_cmp(&right_activity)
                .then_with(|| compare_haplotype_names(left, right)),
        }
    }

    fn compare_dpyd_diplotype_pairs(
        &self,
        left: &(String, String),
        right: &(String, String),
    ) -> std::cmp::Ordering {
        let left_activity = self.dpyd_haplotype_activity_or_missing(&left.0)
            + self.dpyd_haplotype_activity_or_missing(&left.1);
        let right_activity = self.dpyd_haplotype_activity_or_missing(&right.0)
            + self.dpyd_haplotype_activity_or_missing(&right.1);
        left_activity
            .total_cmp(&right_activity)
            .then_with(|| compare_haplotype_names(&left.0, &right.0))
            .then_with(|| compare_haplotype_names(&left.1, &right.1))
    }

    fn dpyd_haplotype_activity_or_missing(&self, haplotype: &str) -> f32 {
        self.dpyd_haplotype_activity(haplotype).unwrap_or(2.0)
    }

    fn dpyd_haplotype_activity(&self, haplotype: &str) -> Option<f32> {
        self.haplotype_activity(haplotype)?.parse().ok()
    }

    /// Infers the recommendation diplotype from RYR1 matcher diplotype labels like Java
    /// `LowestFunctionGeneCaller.inferFromDiplotypes`.
    pub fn infer_ryr1_lowest_function_from_diplotypes(
        &self,
        diplotypes: impl IntoIterator<Item = impl AsRef<str>>,
    ) -> Option<DiplotypeAnnotationInput> {
        let mut candidates = diplotypes
            .into_iter()
            .filter_map(|diplotype| {
                let (hap1, hap2) = split_diplotype_label(diplotype.as_ref());
                let allele1 = self.best_ryr1_haplotype(split_haplotype_label(&hap1))?;
                let allele2 = self.best_ryr1_haplotype(split_haplotype_label(hap2.as_deref()?))?;
                Some(canonical_diplotype_pair(allele1, allele2))
            })
            .collect::<Vec<_>>();

        candidates.sort_by(|left, right| self.compare_ryr1_diplotype_pairs(left, right));
        candidates.into_iter().next().map(|(allele1, allele2)| {
            DiplotypeAnnotationInput::from_alleles(&self.gene, allele1, Some(allele2))
        })
    }

    /// Infers the recommendation diplotype from unphased RYR1 haplotype labels like Java
    /// `LowestFunctionGeneCaller.inferFromHaplotypeMatches`.
    pub fn infer_ryr1_lowest_function_from_haplotypes(
        &self,
        haplotypes: impl IntoIterator<Item = impl AsRef<str>>,
    ) -> Option<DiplotypeAnnotationInput> {
        let mut haplotypes = haplotypes
            .into_iter()
            .flat_map(|haplotype| split_haplotype_label(haplotype.as_ref()))
            .collect::<Vec<_>>();
        haplotypes.sort_by(|left, right| self.compare_ryr1_haplotypes(left, right));

        let allele1 = haplotypes.first()?.clone();
        let allele2 = haplotypes.get(1).cloned();
        let (allele1, allele2) = if let Some(allele2) = allele2 {
            let (allele1, allele2) = canonical_diplotype_pair(allele1, allele2);
            (allele1, Some(allele2))
        } else {
            (allele1, None)
        };

        Some(DiplotypeAnnotationInput::from_alleles(
            &self.gene, allele1, allele2,
        ))
    }

    fn best_ryr1_haplotype(&self, haplotypes: Vec<String>) -> Option<String> {
        haplotypes
            .into_iter()
            .min_by(|left, right| self.compare_ryr1_haplotypes(left, right))
    }

    fn compare_ryr1_haplotypes(&self, left: &str, right: &str) -> std::cmp::Ordering {
        let left_malignant = self.haplotype_function(left).contains("Malignant");
        let right_malignant = self.haplotype_function(right).contains("Malignant");
        right_malignant
            .cmp(&left_malignant)
            .then_with(|| compare_haplotype_names(left, right))
    }

    fn compare_ryr1_diplotype_pairs(
        &self,
        left: &(String, String),
        right: &(String, String),
    ) -> std::cmp::Ordering {
        let left_malignant =
            self.ryr1_malignant_count(&left.0) + self.ryr1_malignant_count(&left.1);
        let right_malignant =
            self.ryr1_malignant_count(&right.0) + self.ryr1_malignant_count(&right.1);
        right_malignant
            .cmp(&left_malignant)
            .then_with(|| compare_haplotype_names(&left.0, &right.0))
            .then_with(|| compare_haplotype_names(&left.1, &right.1))
    }

    fn ryr1_malignant_count(&self, haplotype: &str) -> usize {
        usize::from(self.haplotype_function(haplotype).contains("Malignant"))
    }
}

/// Input for Java `Diplotype.annotateDiplotype`-style phenotype annotation.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct DiplotypeAnnotationInput {
    /// HGNC gene symbol.
    pub gene: String,
    /// First haplotype name, if provided.
    pub allele1: Option<String>,
    /// Second haplotype name, if provided.
    pub allele2: Option<String>,
    /// Pre-existing phenotype values, usually from an outside call.
    pub phenotypes: Vec<String>,
    /// Pre-existing activity score, usually from an outside call.
    pub activity_score: Option<String>,
    /// Whether `phenotypes` came from an outside call.
    pub outside_phenotype: bool,
    /// Whether `activity_score` came from an outside call.
    pub outside_activity_score: bool,
}

impl DiplotypeAnnotationInput {
    /// Creates input for a called diplotype.
    pub fn from_alleles(
        gene: impl Into<String>,
        allele1: impl Into<String>,
        allele2: Option<impl Into<String>>,
    ) -> Self {
        Self {
            gene: gene.into(),
            allele1: Some(allele1.into()),
            allele2: allele2.map(Into::into),
            ..Self::default()
        }
    }

    /// Creates input for an outside call.
    pub fn outside_call(
        gene: impl Into<String>,
        allele1: Option<impl Into<String>>,
        allele2: Option<impl Into<String>>,
        phenotype: Option<impl Into<String>>,
        activity_score: Option<impl Into<String>>,
    ) -> Self {
        let phenotype = phenotype.map(Into::into);
        let activity_score = activity_score.map(Into::into);
        Self {
            gene: gene.into(),
            allele1: allele1.map(Into::into),
            allele2: allele2.map(Into::into),
            outside_phenotype: phenotype.is_some(),
            phenotypes: phenotype.into_iter().collect(),
            outside_activity_score: activity_score.is_some(),
            activity_score,
        }
    }

    /// Annotates this input like Java `Diplotype.annotateDiplotype`.
    pub fn annotate(
        mut self,
        gene_phenotype: Option<&GenePhenotype>,
    ) -> Result<AnnotatedDiplotype, PhenotypeLookupError> {
        if self.gene.starts_with("HLA") {
            return self.annotate_hla();
        }

        let Some(gene_phenotype) = gene_phenotype else {
            return Ok(AnnotatedDiplotype::from_input(self));
        };

        let outside_activity_score_mismatch =
            self.outside_activity_score_mismatch(gene_phenotype)?;
        let outside_phenotype_mismatch = self.outside_phenotype_mismatch(gene_phenotype)?;

        let diplotype_key = self.compute_diplotype_key();
        let diplotype = gene_phenotype.find_diplotype(&diplotype_key);
        let mut lookup_keys = Vec::new();

        if gene_phenotype.is_activity_gene() {
            if self.activity_score.is_none() {
                if self.phenotypes.is_empty() {
                    self.activity_score =
                        match diplotype.and_then(|record| record.activity_score.clone()) {
                            Some(activity_score) => Some(activity_score),
                            None => normalize_activity_score(Some(
                                &gene_phenotype.lookup_activity_by_diplotype_key(
                                    &diplotype_key,
                                    self.is_unknown_alleles(),
                                )?,
                            )),
                        };
                } else {
                    self.activity_score = Some(NA.to_owned());
                }
            }

            if self.phenotypes.is_empty() {
                if self.outside_activity_score {
                    self.phenotypes = gene_phenotype
                        .lookup_phenotypes_by_activity_score(
                            self.activity_score.as_deref().unwrap_or(NA),
                        )
                        .into_iter()
                        .collect();
                } else if let Some(diplotype) = diplotype {
                    self.phenotypes.push(diplotype.phenotype.clone());
                } else {
                    self.phenotypes
                        .push(gene_phenotype.lookup_phenotype_by_diplotype_key(
                            &diplotype_key,
                            self.is_unknown_alleles(),
                            true,
                        )?);
                }
            }

            if self.activity_score.as_deref() != Some(NA) {
                if let Some(activity_score) = &self.activity_score {
                    lookup_keys.push(activity_score.clone());
                }
            } else if self.is_unknown() {
                lookup_keys.push(NO_RESULT.to_owned());
            } else {
                lookup_keys.extend(gene_phenotype.lookup_activity_scores_by_phenotype(
                    self.phenotypes.iter().map(String::as_str),
                ));
            }
        } else {
            if self.phenotypes.is_empty() {
                self.phenotypes = vec![gene_phenotype.lookup_phenotype_by_diplotype_key(
                    &diplotype_key,
                    self.is_unknown_alleles(),
                    false,
                )?];
            }
            if self.is_unknown() {
                lookup_keys.push(NO_RESULT.to_owned());
            } else {
                lookup_keys = self.phenotypes.clone();
            }
        }

        Ok(AnnotatedDiplotype {
            gene: self.gene,
            allele1: self.allele1,
            allele2: self.allele2,
            phenotypes: self.phenotypes,
            activity_score: self.activity_score,
            lookup_keys,
            is_activity_score_type: gene_phenotype.is_activity_gene(),
            outside_phenotype: self.outside_phenotype,
            outside_phenotype_mismatch,
            outside_activity_score: self.outside_activity_score,
            outside_activity_score_mismatch,
        })
    }

    fn compute_diplotype_key(&self) -> BTreeMap<String, i32> {
        compute_lookup_map(self.allele1.as_deref(), self.allele2.as_deref(), None)
    }

    fn is_unknown_alleles(&self) -> bool {
        is_unknown_allele(self.allele1.as_deref()) && is_unknown_allele(self.allele2.as_deref())
    }

    fn is_unknown_phenotype(&self) -> bool {
        self.phenotypes.is_empty()
            || self
                .phenotypes
                .iter()
                .any(|phenotype| phenotype == NO_RESULT)
    }

    fn is_unknown(&self) -> bool {
        self.is_unknown_phenotype() && self.is_unknown_alleles()
    }

    fn has_activity_score(&self) -> bool {
        !is_unspecified(self.activity_score.as_deref())
    }

    fn annotate_hla(mut self) -> Result<AnnotatedDiplotype, PhenotypeLookupError> {
        let mut lookup_keys = Vec::new();

        if !self.is_unknown_alleles() && self.is_unknown_phenotype() {
            self.phenotypes = make_hla_phenotype(&self)?;
            lookup_keys = self.phenotypes.clone();
        } else if self.is_unknown_alleles() {
            if self.is_unknown_phenotype() {
                lookup_keys.push(NO_RESULT.to_owned());
            } else {
                lookup_keys = self.phenotypes.clone();
            }
        }

        Ok(AnnotatedDiplotype {
            gene: self.gene,
            allele1: self.allele1,
            allele2: self.allele2,
            phenotypes: self.phenotypes,
            activity_score: self.activity_score,
            lookup_keys,
            is_activity_score_type: false,
            outside_phenotype: self.outside_phenotype,
            outside_phenotype_mismatch: None,
            outside_activity_score: self.outside_activity_score,
            outside_activity_score_mismatch: None,
        })
    }

    fn outside_activity_score_mismatch(
        &self,
        gene_phenotype: &GenePhenotype,
    ) -> Result<Option<String>, PhenotypeLookupError> {
        if !self.outside_activity_score {
            return Ok(None);
        }
        let Some(activity_score) = self.activity_score.as_deref() else {
            return Ok(None);
        };

        if let Some(phenotype) = self.phenotypes.first() {
            let expected = gene_phenotype.lookup_activity_scores_by_phenotype([phenotype.as_str()]);
            if !expected.contains(activity_score) {
                return Ok(Some(expected.into_iter().collect::<Vec<_>>().join(", ")));
            }
            return Ok(None);
        }

        if !self.is_unknown_alleles() {
            let lookup_map =
                compute_lookup_map(self.allele1.as_deref(), self.allele2.as_deref(), None);
            let expected = normalize_activity_score(Some(
                &gene_phenotype.lookup_key_by_diplotype_key("activity score", &lookup_map)?,
            ))
            .unwrap_or_else(|| NA.to_owned());
            if activity_score != expected {
                return Ok(Some(expected));
            }
        }

        Ok(None)
    }

    fn outside_phenotype_mismatch(
        &self,
        gene_phenotype: &GenePhenotype,
    ) -> Result<Option<String>, PhenotypeLookupError> {
        if !self.outside_phenotype {
            return Ok(None);
        }
        let Some(phenotype) = self.phenotypes.first() else {
            return Ok(None);
        };

        if gene_phenotype.is_activity_gene() && self.has_activity_score() {
            let expected = gene_phenotype
                .lookup_phenotypes_by_activity_score(self.activity_score.as_deref().unwrap_or(NA));
            if !expected.contains(phenotype) {
                return Ok(Some(expected.into_iter().collect::<Vec<_>>().join(", ")));
            }
            return Ok(None);
        }

        if !self.is_unknown_alleles() {
            let expected = gene_phenotype.lookup_phenotype_by_diplotype_key(
                &self.compute_diplotype_key(),
                false,
                gene_phenotype.is_activity_gene(),
            )?;
            if phenotype != &expected {
                return Ok(Some(expected));
            }
        }

        Ok(None)
    }
}

/// Java `Diplotype` fields affected by phenotype annotation.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct AnnotatedDiplotype {
    /// HGNC gene symbol.
    pub gene: String,
    /// First haplotype name, if provided.
    pub allele1: Option<String>,
    /// Second haplotype name, if provided.
    pub allele2: Option<String>,
    /// Annotated phenotypes.
    pub phenotypes: Vec<String>,
    /// Annotated activity score.
    pub activity_score: Option<String>,
    /// Recommendation lookup keys.
    pub lookup_keys: Vec<String>,
    /// Whether the gene phenotype model is activity-score based.
    pub is_activity_score_type: bool,
    /// Whether phenotype came from an outside call.
    pub outside_phenotype: bool,
    /// Expected phenotype when an outside phenotype mismatches.
    pub outside_phenotype_mismatch: Option<String>,
    /// Whether activity score came from an outside call.
    pub outside_activity_score: bool,
    /// Expected activity score when an outside activity score mismatches.
    pub outside_activity_score_mismatch: Option<String>,
}

impl AnnotatedDiplotype {
    fn from_input(input: DiplotypeAnnotationInput) -> Self {
        Self {
            gene: input.gene,
            allele1: input.allele1,
            allele2: input.allele2,
            phenotypes: input.phenotypes,
            activity_score: input.activity_score,
            outside_phenotype: input.outside_phenotype,
            outside_activity_score: input.outside_activity_score,
            ..Self::default()
        }
    }
}

/// Validation data used when parsing Java outside-call TSV records.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct OutsideCallValidation {
    /// Genes known to this parser. Unknown genes are parsed without allele validation.
    pub supported_genes: BTreeSet<String>,
    /// Genes that use activity scores.
    pub activity_score_genes: BTreeSet<String>,
    /// Genes whose named alleles are variants in Java warning text.
    pub variant_genes: BTreeSet<String>,
    /// Valid named alleles keyed by gene. Missing entries skip named-allele validation.
    pub valid_named_alleles: BTreeMap<String, BTreeSet<String>>,
}

impl OutsideCallValidation {
    /// Builds validation with the provided genes marked as supported.
    pub fn for_supported_genes(genes: impl IntoIterator<Item = impl Into<String>>) -> Self {
        Self {
            supported_genes: genes.into_iter().map(Into::into).collect(),
            ..Self::default()
        }
    }

    /// Returns whether `gene` is known to the PharmCAT environment.
    pub fn supports_gene(&self, gene: &str) -> bool {
        self.supported_genes.contains(gene)
    }

    /// Returns whether `gene` uses activity scores.
    pub fn is_activity_score_gene(&self, gene: &str) -> bool {
        self.activity_score_genes.contains(gene)
    }

    /// Returns whether `gene` is represented as a variant gene in warning text.
    pub fn is_variant_gene(&self, gene: &str) -> bool {
        self.variant_genes.contains(gene)
    }

    /// Returns whether `allele` is valid for `gene`, or `true` when no allele set is loaded.
    pub fn is_valid_named_allele(&self, gene: &str, allele: &str) -> bool {
        self.valid_named_alleles
            .get(gene)
            .is_none_or(|alleles| alleles.contains(allele))
    }
}

/// One Java `OutsideCall` parsed from an outside-call TSV row.
#[derive(Clone, Debug)]
pub struct OutsideCall {
    /// HGNC gene symbol.
    pub gene: String,
    /// Normalized diplotype, if supplied.
    pub diplotype: Option<String>,
    /// Normalized phenotype, if supplied.
    pub phenotype: Option<String>,
    /// Activity score, if supplied.
    pub activity_score: Option<String>,
    /// Parsed haplotypes from `diplotype`.
    pub haplotypes: BTreeSet<String>,
    /// Java-compatible parse warnings.
    pub warnings: Vec<String>,
}

impl OutsideCall {
    /// Returns whether no diplotype, phenotype, or activity score is available.
    pub fn is_no_call(&self) -> bool {
        self.diplotype.is_none() && self.phenotype.is_none() && self.activity_score.is_none()
    }

    /// Converts this parsed outside call into Java `Diplotype(OutsideCall, Env)` annotation input.
    pub fn to_annotation_input(&self) -> Result<DiplotypeAnnotationInput, OutsideCallError> {
        let alleles = self
            .diplotype
            .as_deref()
            .map(|diplotype| split_outside_diplotype(&self.gene, diplotype))
            .transpose()?;
        let allele1 = alleles
            .as_ref()
            .and_then(|alleles| alleles.first())
            .cloned();
        let allele2 = alleles.as_ref().and_then(|alleles| alleles.get(1)).cloned();

        Ok(DiplotypeAnnotationInput::outside_call(
            self.gene.clone(),
            allele1,
            allele2,
            self.phenotype.clone(),
            self.activity_score.clone(),
        ))
    }

    fn comparison_key(&self) -> (&str, Option<&str>, Option<&str>, Option<&str>) {
        (
            &self.gene,
            self.diplotype.as_deref(),
            self.phenotype.as_deref(),
            self.activity_score.as_deref(),
        )
    }
}

impl std::fmt::Display for OutsideCall {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if let Some(diplotype) = &self.diplotype {
            write!(f, "{}:{diplotype}", self.gene)
        } else if let Some(phenotype) = &self.phenotype {
            write!(f, "{}:{phenotype}", self.gene)
        } else {
            write!(
                f,
                "{}:{}",
                self.gene,
                self.activity_score.as_deref().unwrap_or("")
            )
        }
    }
}

impl Eq for OutsideCall {}

impl PartialEq for OutsideCall {
    fn eq(&self, other: &Self) -> bool {
        self.comparison_key() == other.comparison_key()
    }
}

impl Ord for OutsideCall {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.comparison_key().cmp(&other.comparison_key())
    }
}

impl PartialOrd for OutsideCall {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

/// Error from Java outside-call parsing.
#[derive(Debug)]
pub enum OutsideCallError {
    /// Invalid outside-call data.
    Invalid(String),
    /// File read failure.
    Io(io::Error),
}

impl std::fmt::Display for OutsideCallError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Invalid(message) => f.write_str(message),
            Self::Io(error) => error.fmt(f),
        }
    }
}

impl std::error::Error for OutsideCallError {}

impl From<io::Error> for OutsideCallError {
    fn from(error: io::Error) -> Self {
        Self::Io(error)
    }
}

fn split_outside_diplotype(
    gene: &str,
    diplotype_text: &str,
) -> Result<Vec<String>, OutsideCallError> {
    if diplotype_text.contains('/') {
        if is_single_ploidy_gene(gene) && !is_x_chromosome_gene(gene) {
            return Err(OutsideCallError::Invalid(format!(
                "Cannot specify two genotypes [{diplotype_text}] for single chromosome gene {gene}"
            )));
        }
        let mut alleles = diplotype_text
            .split('/')
            .map(str::to_owned)
            .collect::<Vec<_>>();
        if alleles.len() != 2 {
            return Err(OutsideCallError::Invalid(format!(
                "Diplotype for {gene} has {} alleles: ({diplotype_text})",
                alleles.len()
            )));
        }
        // Java OutsideCall stores haplotypes in a TreeSet ordered by HaplotypeNameComparator, so an
        // unordered outside diplotype (e.g. *2/*1) normalizes to its ordered form (*1/*2).
        alleles.sort_by(|left, right| crate::matcher::compare_haplotype_names(left, right));
        Ok(alleles)
    } else if !is_single_ploidy_gene(gene) && !is_phenotype_only_gene(gene) {
        Err(OutsideCallError::Invalid(format!(
            "Expected two genotypes separated by a '/' but saw [{diplotype_text}] for {gene}"
        )))
    } else {
        Ok(vec![diplotype_text.to_owned()])
    }
}

fn is_single_ploidy_gene(gene: &str) -> bool {
    matches!(gene, "G6PD" | "MT-RNR1")
}

fn is_x_chromosome_gene(gene: &str) -> bool {
    gene == "G6PD"
}

fn is_phenotype_only_gene(gene: &str) -> bool {
    matches!(gene, "HLA-A" | "HLA-B")
}

/// Parses one Java outside-call TSV line.
pub fn parse_outside_call_line(
    validation: &OutsideCallValidation,
    line: &str,
    line_number: usize,
) -> Result<OutsideCall, OutsideCallError> {
    let fields = line.split('\t').map(str::trim).collect::<Vec<_>>();
    if fields.len() < 2 {
        return Err(OutsideCallError::Invalid(format!(
            "Line {line_number}: Expected at least 2 TSV fields, got {}",
            fields.len()
        )));
    }

    let gene = strip_to_null(fields[0])
        .ok_or_else(|| OutsideCallError::Invalid(format!("Line {line_number}: No gene specified")))?
        .to_owned();
    let mut diplotype = strip_to_null(fields[1]).map(str::to_owned);
    if diplotype.as_deref() == Some("/") {
        diplotype = None;
    }
    let mut phenotype = fields
        .get(2)
        .and_then(|value| normalize_phenotype(Some(&value.replace(&gene, ""))));
    let mut activity_score = fields
        .get(3)
        .and_then(|value| strip_to_null(value))
        .map(str::to_owned);
    let mut haplotypes = BTreeSet::new();
    let mut warnings = Vec::new();

    if !validation.supports_gene(&gene) {
        return Ok(OutsideCall {
            gene,
            diplotype,
            phenotype,
            activity_score,
            haplotypes,
            warnings,
        });
    }

    if phenotype.is_none() && activity_score.is_none() && diplotype.is_none() {
        return Err(OutsideCallError::Invalid(format!(
            "Specify a diplotype, phenotype, or activity score for {gene}"
        )));
    }

    if fields.len() == 2 && diplotype.is_none() {
        let message = if fields[1].trim().is_empty() {
            "No diplotype specified"
        } else {
            "Invalid diplotype specified"
        };
        return Err(OutsideCallError::Invalid(format!(
            "Line {line_number}: {message}"
        )));
    }

    if let Some(current_diplotype) = diplotype.as_deref() {
        if current_diplotype == "." || current_diplotype == "./." {
            diplotype = None;
        } else {
            let mut alleles = current_diplotype
                .split('/')
                .map(str::trim)
                .map(|allele| strip_gene_prefix(&gene, allele))
                .collect::<Vec<_>>();
            if alleles.len() > 2 {
                return Err(OutsideCallError::Invalid(format!(
                    "Line {line_number}: Too many alleles specified in {current_diplotype}"
                )));
            }

            if gene == "CYP2D6" {
                for allele in &mut alleles {
                    if let Some(base) = cyp2d6_suballele_base(allele) {
                        warnings.push(format!(
                            "PharmCAT does not support sub-alleles for {gene}. Using '{base}' instead of '{allele}'."
                        ));
                        *allele = base;
                    }
                }
            } else if gene == "HLA-A" || gene == "HLA-B" {
                normalize_hla_alleles(&gene, &mut alleles, &mut warnings)?;
            }

            for allele in &alleles {
                if !validation.is_valid_named_allele(&gene, allele) {
                    let named_type = if validation.is_variant_gene(&gene) {
                        "variant"
                    } else {
                        "allele"
                    };
                    warnings.push(format!(
                        "Undocumented {gene} named {named_type} in outside call: {allele}"
                    ));
                }
            }

            diplotype = Some(alleles.join("/"));
            if let Some(first) = alleles.first() {
                haplotypes.insert(first.clone());
            }
            if let Some(second) = alleles.get(1) {
                haplotypes.insert(second.clone());
            }
        }
    } else if validation.is_activity_score_gene(&gene) {
        warnings.push(format!(
            "{gene} is not an activity score gene but has outside call with only an activity score.  PharmCAT will not be able to provide any recommendations based on this gene."
        ));
    }

    if let Some(value) = fields.get(2) {
        phenotype = normalize_phenotype(Some(&value.replace(&gene, "")));
    }
    if let Some(value) = fields.get(3) {
        activity_score = strip_to_null(value).map(str::to_owned);
    }

    Ok(OutsideCall {
        gene,
        diplotype,
        phenotype,
        activity_score,
        haplotypes,
        warnings,
    })
}

/// Parses outside-call TSV file data like Java `OutsideCallParser.parse(Env, Path)`.
pub fn parse_outside_calls_file(
    validation: &OutsideCallValidation,
    path: &Path,
) -> Result<Vec<OutsideCall>, OutsideCallError> {
    let data = fs::read_to_string(path)?;
    let mut calls = Vec::new();
    for (index, line) in data.lines().enumerate() {
        if is_non_comment_line(line) {
            calls.push(parse_outside_call_line(validation, line, index + 1)?);
        }
    }
    Ok(calls)
}

/// Parses outside-call TSV string data like Java `OutsideCallParser.parse(Env, String)`.
pub fn parse_outside_calls_str(
    validation: &OutsideCallValidation,
    data: &str,
) -> Result<BTreeSet<OutsideCall>, OutsideCallError> {
    let stripped = data.trim();
    let mut calls = BTreeSet::new();
    for (index, line) in stripped.split('\n').enumerate() {
        if is_non_comment_line(line) {
            calls.insert(parse_outside_call_line(validation, line, index + 1)?);
        }
    }
    Ok(calls)
}

fn strip_to_null(value: &str) -> Option<&str> {
    let trimmed = value.trim();
    (!trimmed.is_empty()).then_some(trimmed)
}

fn is_non_comment_line(line: &str) -> bool {
    !line.trim().is_empty() && !line.starts_with('#')
}

fn strip_gene_prefix(gene: &str, allele: &str) -> String {
    allele
        .strip_prefix(gene)
        .map(str::trim_start)
        .unwrap_or(allele)
        .to_owned()
}

fn cyp2d6_suballele_base(allele: &str) -> Option<String> {
    let (base, suballele) = allele.split_once('.')?;
    if base.len() > 1
        && base.starts_with('*')
        && base[1..].chars().all(|char| char.is_ascii_digit())
        && !suballele.is_empty()
        && suballele.chars().all(|char| char.is_ascii_digit())
    {
        return Some(base.to_owned());
    }
    None
}

fn normalize_hla_alleles(
    gene: &str,
    alleles: &mut [String],
    warnings: &mut Vec<String>,
) -> Result<(), OutsideCallError> {
    let prefix = format!("{}*", &gene[gene.len() - 1..]);
    for allele in alleles {
        let original = allele.clone();
        if !allele.starts_with('*') {
            if !allele.starts_with(&prefix) {
                return Err(OutsideCallError::Invalid(format!(
                    "Invalid {gene} allele: '{original}'."
                )));
            }
            *allele = allele[1..].to_owned();
        }

        let mut is_suballele = false;
        if let Some(base) = hla_suballele_base(allele) {
            *allele = base;
            is_suballele = true;
        }

        if original != *allele {
            if is_suballele {
                warnings.push(format!(
                    "PharmCAT does not support sub-alleles for {gene}. Using '{allele}' instead of '{original}'."
                ));
            } else {
                warnings.push(format!(
                    "Converting outside call for {gene} from '{original}', to '{allele}'."
                ));
            }
        }
    }
    Ok(())
}

fn hla_suballele_base(allele: &str) -> Option<String> {
    let parts = allele.split(':').collect::<Vec<_>>();
    if parts.len() < 3 {
        return None;
    }
    let first = parts[0];
    if first.len() > 1
        && first.starts_with('*')
        && first[1..].chars().all(|char| char.is_ascii_digit())
        && parts[1].chars().all(|char| char.is_ascii_digit())
        && parts[2..]
            .iter()
            .all(|part| !part.is_empty() && part.chars().all(|char| char.is_ascii_digit()))
    {
        return Some(format!("{first}:{}", parts[1]));
    }
    None
}

fn compute_lookup_map(
    allele1: Option<&str>,
    allele2: Option<&str>,
    phenotype: Option<&str>,
) -> BTreeMap<String, i32> {
    let mut lookup_map = BTreeMap::new();

    match (allele1, allele2) {
        (None, None) => {
            if let Some(phenotype) = phenotype {
                lookup_map.insert(phenotype.to_owned(), 1);
            }
        }
        (Some(allele1), Some(allele2)) if allele1 == allele2 => {
            lookup_map.insert(allele1.to_owned(), 2);
        }
        (Some(allele1), Some(allele2)) => {
            lookup_map.insert(allele1.to_owned(), 1);
            lookup_map.insert(allele2.to_owned(), 1);
        }
        (Some(allele1), None) => {
            lookup_map.insert(allele1.to_owned(), 1);
        }
        (None, Some(_)) => {}
    }

    lookup_map
}

fn is_unknown_allele(allele: Option<&str>) -> bool {
    allele.is_none_or(|allele| allele == "Unknown")
}

fn is_unspecified(value: Option<&str>) -> bool {
    value.is_none_or(|value| value.trim().is_empty() || value.trim().eq_ignore_ascii_case(NA))
}

fn make_hla_phenotype(
    input: &DiplotypeAnnotationInput,
) -> Result<Vec<String>, PhenotypeLookupError> {
    let alleles = match input.gene.as_str() {
        "HLA-A" => ["*31:01"].as_slice(),
        "HLA-B" => ["*15:02", "*57:01", "*58:01"].as_slice(),
        _ => {
            return Err(PhenotypeLookupError::UnsupportedHlaGene(input.gene.clone()));
        }
    };

    Ok(alleles
        .iter()
        .map(|allele| {
            if hla_name_contains_allele(input.allele1.as_deref(), allele)
                || hla_name_contains_allele(input.allele2.as_deref(), allele)
            {
                format!("{allele} positive")
            } else {
                format!("{allele} negative")
            }
        })
        .collect())
}

fn hla_name_contains_allele(haplotype: Option<&str>, allele: &str) -> bool {
    haplotype.is_some_and(|haplotype| haplotype.contains(allele))
}

fn split_diplotype_label(diplotype: &str) -> (String, Option<String>) {
    let mut depth = 0;
    for (index, character) in diplotype.char_indices() {
        match character {
            '[' => depth += 1,
            ']' => depth -= 1,
            '/' if depth == 0 => {
                return (
                    diplotype[..index].to_owned(),
                    Some(diplotype[index + 1..].to_owned()),
                );
            }
            _ => {}
        }
    }
    (diplotype.to_owned(), None)
}

fn split_haplotype_label(haplotype: &str) -> Vec<String> {
    if haplotype.starts_with('[') && haplotype.ends_with(']') && haplotype.contains(" + ") {
        haplotype
            .trim_start_matches('[')
            .trim_end_matches(']')
            .split(" + ")
            .map(str::to_owned)
            .collect()
    } else {
        vec![haplotype.to_owned()]
    }
}

fn canonical_diplotype_pair(left: String, right: String) -> (String, String) {
    if compare_haplotype_names(&right, &left).is_lt() {
        (right, left)
    } else {
        (left, right)
    }
}

/// One Java `DiplotypeRecord`.
#[derive(Clone, Debug, Deserialize, Eq, PartialEq)]
pub struct DiplotypeRecord {
    /// Gene-level result.
    #[serde(rename = "generesult")]
    pub gene_result: String,
    /// Diplotype display name.
    pub diplotype: String,
    /// Optional description.
    #[serde(default)]
    pub description: Option<String>,
    /// Recommendation lookup key.
    #[serde(rename = "lookupkey")]
    pub lookup_key: String,
    /// Diplotype key.
    #[serde(rename = "diplotypekey", default)]
    pub diplotype_key: BTreeMap<String, i32>,
    /// Activity score.
    #[serde(rename = "activityScore", default)]
    pub activity_score: Option<String>,
    /// Phenotype.
    pub phenotype: String,
}

impl DiplotypeRecord {
    /// Returns whether this record matches a diplotype key like Java `matchesKey`.
    pub fn matches_key(&self, other_key: &BTreeMap<String, i32>) -> bool {
        !other_key.is_empty()
            && !self.diplotype_key.is_empty()
            && other_key.len() == self.diplotype_key.len()
            && self
                .diplotype_key
                .iter()
                .all(|(key, value)| other_key.get(key) == Some(value))
    }
}

impl Ord for DiplotypeRecord {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.diplotype
            .cmp(&other.diplotype)
            .then_with(|| self.activity_score.cmp(&other.activity_score))
            .then_with(|| self.gene_result.cmp(&other.gene_result))
            .then_with(|| self.phenotype.cmp(&other.phenotype))
    }
}

impl PartialOrd for DiplotypeRecord {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

/// One Java `HaplotypeRecord`.
#[derive(Clone, Debug, Deserialize, Eq, PartialEq)]
#[serde(rename_all = "camelCase")]
pub struct HaplotypeRecord {
    /// Haplotype name.
    pub name: String,
    /// Activity value.
    #[serde(default)]
    pub activity_value: Option<String>,
    /// Function value.
    #[serde(default)]
    pub function_value: Option<String>,
    /// Lookup key.
    pub lookup_key: String,
}

/// Phenotype lookup error.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum PhenotypeLookupError {
    /// More than one unique Java lookup key matched a single diplotype key.
    MultipleMatches {
        /// Lookup key type.
        key_type: &'static str,
        /// Matched values.
        values: BTreeSet<String>,
    },
    /// HLA phenotype calling was requested for an unsupported HLA gene.
    UnsupportedHlaGene(String),
}

impl std::fmt::Display for PhenotypeLookupError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::MultipleMatches { key_type, values } => {
                write!(
                    f,
                    "More than one {key_type} matched: {}",
                    values.iter().cloned().collect::<Vec<_>>().join("; ")
                )
            }
            Self::UnsupportedHlaGene(gene) => {
                write!(f, "Gene not supported for HLA phenotype calling: {gene}")
            }
        }
    }
}

impl std::error::Error for PhenotypeLookupError {}

fn unique_lookup(
    key_type: &'static str,
    keys: BTreeSet<String>,
    empty_value: &str,
) -> Result<String, PhenotypeLookupError> {
    if keys.len() > 1 {
        Err(PhenotypeLookupError::MultipleMatches {
            key_type,
            values: keys,
        })
    } else {
        Ok(keys
            .into_iter()
            .next()
            .unwrap_or_else(|| empty_value.to_owned()))
    }
}

/// Reads one phenotype JSON file.
pub fn read_gene_phenotype_file(path: &Path) -> Result<GenePhenotype, PhenotypeLoadError> {
    // Phenotype resources are up to 20 MB. Avoid `from_reader` on an unbuffered `File`, whose
    // byte-at-a-time adapter otherwise turns resource initialization into millions of syscalls.
    let data = fs::read(path)?;
    Ok(serde_json::from_slice(&data)?)
}

/// Phenotype loading error.
#[derive(Debug)]
pub enum PhenotypeLoadError {
    /// Path was not a directory.
    NotDirectory(std::path::PathBuf),
    /// Directory had no phenotype JSON files.
    NoPhenotypeFiles(std::path::PathBuf),
    /// Duplicate gene file.
    DuplicateGene(String),
    /// I/O error.
    Io(io::Error),
    /// JSON parse error.
    Json(serde_json::Error),
}

impl std::fmt::Display for PhenotypeLoadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NotDirectory(path) => write!(f, "{} is not a directory", path.display()),
            Self::NoPhenotypeFiles(path) => {
                write!(f, "Cannot find phenotype files in {}", path.display())
            }
            Self::DuplicateGene(gene) => write!(f, "Multiple GenePhenotypes for {gene}"),
            Self::Io(error) => error.fmt(f),
            Self::Json(error) => error.fmt(f),
        }
    }
}

impl std::error::Error for PhenotypeLoadError {}

impl From<io::Error> for PhenotypeLoadError {
    fn from(error: io::Error) -> Self {
        Self::Io(error)
    }
}

impl From<serde_json::Error> for PhenotypeLoadError {
    fn from(error: serde_json::Error) -> Self {
        Self::Json(error)
    }
}

/// Normalizes an activity value or score like Java `ActivityUtils.normalize`.
pub fn normalize_activity_score(value: Option<&str>) -> Option<String> {
    let trimmed = value?.trim();
    if trimmed.is_empty() {
        return None;
    }

    if let Some(rest) = trimmed.strip_prefix('>')
        && let Some(normalized) = normalize_decimal(rest)
    {
        return Some(format!(">{normalized}"));
    }
    if let Some(rest) = trimmed.strip_prefix(GTE)
        && let Some(normalized) = normalize_decimal(rest)
    {
        return Some(format!("{GTE}{normalized}"));
    }
    if let Some(normalized) = normalize_decimal(trimmed) {
        return Some(normalized);
    }

    Some(trimmed.to_owned())
}

/// Normalizes common phenotype strings like Java `PhenotypeUtils.normalize`.
pub fn normalize_phenotype(value: Option<&str>) -> Option<String> {
    let value = value?;
    if value.trim().is_empty() {
        return None;
    }

    let collapsed = value.split_whitespace().collect::<Vec<_>>().join(" ");
    let trimmed = normalize_metabolizer_spelling(&collapsed);
    let mut lowered = trimmed.to_lowercase();

    if let Some(expanded) = phenotype_abbreviation(&lowered) {
        return Some(expanded.to_owned());
    }

    if lowered.contains(&INDETERMINATE.to_lowercase()) {
        return Some(INDETERMINATE.to_owned());
    }

    if lowered.contains("extensive") {
        lowered = lowered.replace("extensive", "normal");
    }

    for pheno_name in PHENOTYPE_NAMES {
        if lowered.contains(&pheno_name.to_lowercase()) {
            if lowered.starts_with("likely") {
                return Some(format!("Likely {pheno_name}"));
            }
            if lowered.starts_with("possible") {
                return Some(format!("Possible {pheno_name}"));
            }
            return Some(pheno_name.to_owned());
        }
    }

    Some(trimmed)
}

fn normalize_decimal(value: &str) -> Option<String> {
    if value.is_empty() || !value.chars().all(|c| c.is_ascii_digit() || c == '.') {
        return None;
    }

    let mut decimal_points = value.chars().filter(|c| *c == '.');
    let has_decimal = decimal_points.next().is_some();
    if decimal_points.next().is_some() {
        return None;
    }

    let (whole, fraction) = value.split_once('.').unwrap_or((value, ""));
    if whole.is_empty() || !whole.chars().all(|c| c.is_ascii_digit()) {
        return None;
    }
    if has_decimal && (fraction.is_empty() || !fraction.chars().all(|c| c.is_ascii_digit())) {
        return None;
    }

    if has_decimal {
        Some(value.to_owned())
    } else {
        Some(format!("{value}.0"))
    }
}

const PHENOTYPE_NAMES: [&str; 5] = [
    "Poor Metabolizer",
    "Intermediate Metabolizer",
    "Normal Metabolizer",
    "Normal Metabolizer",
    "Ultrarapid Metabolizer",
];

fn phenotype_abbreviation(value: &str) -> Option<&'static str> {
    match value {
        "pm" => Some("Poor Metabolizer"),
        "im" => Some("Intermediate Metabolizer"),
        "nm" | "em" => Some("Normal Metabolizer"),
        "um" => Some("Ultrarapid Metabolizer"),
        _ => None,
    }
}

fn normalize_metabolizer_spelling(value: &str) -> String {
    let variants = ["metabolizers", "metabolisers", "metabolizer", "metaboliser"];
    let mut normalized = String::new();
    let mut remaining = value;

    while !remaining.is_empty() {
        let lowered = remaining.to_lowercase();
        let Some((start, variant)) = variants
            .iter()
            .filter_map(|variant| lowered.find(variant).map(|start| (start, *variant)))
            .min_by_key(|(start, _)| *start)
        else {
            normalized.push_str(remaining);
            break;
        };

        normalized.push_str(&remaining[..start]);
        normalized.push_str("Metabolizer");
        remaining = &remaining[start + variant.len()..];
    }

    normalized
}

#[cfg(test)]
mod tests {
    use std::{
        collections::{BTreeMap, BTreeSet},
        fs,
        path::{Path, PathBuf},
        time::{SystemTime, UNIX_EPOCH},
    };

    use super::{
        DiplotypeAnnotationInput, GTE, GenePhenotype, HaplotypeRecord, INDETERMINATE, NA,
        NO_RESULT, OutsideCall, OutsideCallValidation, PhenotypeMap, normalize_activity_score,
        normalize_phenotype, parse_outside_call_line, parse_outside_calls_file,
        parse_outside_calls_str, read_gene_phenotype_file,
    };

    #[test]
    fn normalizes_activity_scores_like_java_activity_utils() {
        assert_eq!(normalize_activity_score(Some(" ")), None);
        assert_eq!(normalize_activity_score(None), None);
        assert_eq!(
            normalize_activity_score(Some("1.0")).as_deref(),
            Some("1.0")
        );
        assert_eq!(normalize_activity_score(Some("1")).as_deref(), Some("1.0"));
        assert_eq!(
            normalize_activity_score(Some(&format!("{GTE}1"))).as_deref(),
            Some("\u{2265}1.0")
        );
        assert_eq!(
            normalize_activity_score(Some(">3")).as_deref(),
            Some(">3.0")
        );
        assert_eq!(normalize_activity_score(Some(NA)).as_deref(), Some(NA));
    }

    #[test]
    fn keeps_non_decimal_activity_scores_trimmed_like_java() {
        assert_eq!(
            normalize_activity_score(Some(" activity unknown ")).as_deref(),
            Some("activity unknown")
        );
        assert_eq!(normalize_activity_score(Some("1.")).as_deref(), Some("1."));
        assert_eq!(
            normalize_activity_score(Some(">=1")).as_deref(),
            Some(">=1")
        );
    }

    #[test]
    fn normalizes_phenotypes_like_java_phenotype_utils() {
        assert_eq!(
            normalize_phenotype(Some("Normal Metaboliser")).as_deref(),
            Some("Normal Metabolizer")
        );
        assert_eq!(normalize_phenotype(Some("  ")), None);
        assert_eq!(
            normalize_phenotype(Some("  Normal    Metabolizer")).as_deref(),
            Some("Normal Metabolizer")
        );
        assert_eq!(
            normalize_phenotype(Some("  Normal    Metabolizers")).as_deref(),
            Some("Normal Metabolizer")
        );
        assert_eq!(
            normalize_phenotype(Some("NM")).as_deref(),
            Some("Normal Metabolizer")
        );
        assert_eq!(
            normalize_phenotype(Some("EM")).as_deref(),
            Some("Normal Metabolizer")
        );
        assert_eq!(
            normalize_phenotype(Some("PM")).as_deref(),
            Some("Poor Metabolizer")
        );
        assert_eq!(
            normalize_phenotype(Some("Extensive Metabolisers")).as_deref(),
            Some("Normal Metabolizer")
        );
        assert_eq!(
            normalize_phenotype(Some("Normal MetaboliSer")).as_deref(),
            Some("Normal Metabolizer")
        );
        assert_eq!(
            normalize_phenotype(Some("normal metabolizer (NM)    ")).as_deref(),
            Some("Normal Metabolizer")
        );
        assert_eq!(
            normalize_phenotype(Some("normal   metabolizer (NM)\n")).as_deref(),
            Some("Normal Metabolizer")
        );
        assert_eq!(
            normalize_phenotype(Some("foo\t\t\t bar")).as_deref(),
            Some("foo bar")
        );
        assert_eq!(
            normalize_phenotype(Some("Likely Normal Metaboliser")).as_deref(),
            Some("Likely Normal Metabolizer")
        );
        assert_eq!(
            normalize_phenotype(Some("Likely Normal Metaboliser (NM)")).as_deref(),
            Some("Likely Normal Metabolizer")
        );
    }

    #[test]
    fn normalizes_indeterminate_and_possible_phenotypes_like_java() {
        assert_eq!(
            normalize_phenotype(Some("Possible poor metabolizer")).as_deref(),
            Some("Possible Poor Metabolizer")
        );
        assert_eq!(
            normalize_phenotype(Some("likely intermediate metabolizer")).as_deref(),
            Some("Likely Intermediate Metabolizer")
        );
        assert_eq!(
            normalize_phenotype(Some("DPYD indeterminate result")).as_deref(),
            Some("Indeterminate")
        );
    }

    #[test]
    fn parses_minimal_outside_call_like_java() {
        let validation = outside_call_validation();
        let call = parse_outside_call_line(&validation, "CYP2C9\t*1/*2", 1).expect("parse");

        assert_eq!(call.gene, "CYP2C9");
        assert_eq!(call.diplotype.as_deref(), Some("*1/*2"));
        assert_eq!(call.phenotype, None);
        assert_eq!(
            call.haplotypes,
            BTreeSet::from(["*1".to_owned(), "*2".to_owned()])
        );
    }

    #[test]
    fn strips_gene_prefix_and_normalizes_outside_phenotype_like_java() {
        let validation = outside_call_validation();
        let call = parse_outside_call_line(
            &validation,
            "CYP2C9\tCYP2C9*1/CYP2C9      *2\tCYP2C9 Normal Metabolizer",
            1,
        )
        .expect("parse");

        assert_eq!(call.diplotype.as_deref(), Some("*1/*2"));
        assert_eq!(call.phenotype.as_deref(), Some("Normal Metabolizer"));
    }

    #[test]
    fn parses_outside_call_file_with_comments_like_java() {
        let validation = outside_call_validation();
        let path = write_temp_file(
            "pharmcat-outside-calls",
            "# comment\nCYP2C9\t*1/*2\n\n## another comment\n\nCYP2C9\t*3/*4\n",
        );

        let calls = parse_outside_calls_file(&validation, &path).expect("parse file");

        assert_eq!(calls.len(), 2);
        assert_eq!(calls[0].diplotype.as_deref(), Some("*1/*2"));
        assert_eq!(calls[1].diplotype.as_deref(), Some("*3/*4"));
    }

    #[test]
    fn rejects_bad_outside_call_diplotype_like_java() {
        let validation = outside_call_validation();
        let error =
            parse_outside_call_line(&validation, "CYP2C19\t*3/*4/*2", 2).expect_err("bad call");

        assert_eq!(
            error.to_string(),
            "Line 2: Too many alleles specified in *3/*4/*2"
        );
    }

    #[test]
    fn normalizes_hla_and_cyp2d6_outside_call_suballeles_like_java() {
        let validation = outside_call_validation();
        let cyp2d6 =
            parse_outside_call_line(&validation, "CYP2D6\t*1.001/*4", 1).expect("parse cyp2d6");
        let hla =
            parse_outside_call_line(&validation, "HLA-A\tA*31:01:02/*01:01", 1).expect("parse hla");

        assert_eq!(cyp2d6.diplotype.as_deref(), Some("*1/*4"));
        assert_eq!(
            cyp2d6.warnings,
            ["PharmCAT does not support sub-alleles for CYP2D6. Using '*1' instead of '*1.001'."]
        );
        assert_eq!(hla.diplotype.as_deref(), Some("*31:01/*01:01"));
        assert_eq!(
            hla.warnings,
            [
                "PharmCAT does not support sub-alleles for HLA-A. Using '*31:01' instead of 'A*31:01:02'.",
                "Undocumented HLA-A named allele in outside call: *01:01"
            ]
        );
    }

    #[test]
    fn parses_outside_call_string_as_set_like_java() {
        let validation = outside_call_validation();
        let calls = parse_outside_calls_str(
            &validation,
            "\nCYP2C9\t*1/*2\n# comment\nCYP2C9\t*1/*2\nCYP2C19\t*3/*4\n",
        )
        .expect("parse string");

        assert_eq!(calls.len(), 2);
        assert!(calls.contains(&OutsideCall {
            gene: "CYP2C9".to_owned(),
            diplotype: Some("*1/*2".to_owned()),
            phenotype: None,
            activity_score: None,
            haplotypes: BTreeSet::new(),
            warnings: Vec::new(),
        }));
    }

    #[test]
    fn converts_outside_call_to_annotation_input_like_java_diplotype_constructor() {
        let validation = outside_call_validation();
        let call = parse_outside_call_line(
            &validation,
            "CYP2D6\t*1/*3\tIntermediate Metabolizer\t1.0",
            1,
        )
        .expect("outside call");
        let cyp2d6 = read_test_phenotype("CYP2D6");

        let annotated = call
            .to_annotation_input()
            .expect("annotation input")
            .annotate(Some(&cyp2d6))
            .expect("annotated outside call");

        assert_eq!(annotated.allele1.as_deref(), Some("*1"));
        assert_eq!(annotated.allele2.as_deref(), Some("*3"));
        assert_eq!(annotated.phenotypes, ["Intermediate Metabolizer"]);
        assert_eq!(annotated.activity_score.as_deref(), Some("1.0"));
        assert_eq!(annotated.lookup_keys, ["1.0"]);
        assert!(annotated.is_activity_score_type);
        assert!(annotated.outside_phenotype);
        assert!(annotated.outside_activity_score);
    }

    #[test]
    fn rejects_outside_call_annotation_input_bad_diplotype_shapes_like_java_factory() {
        let single = OutsideCall {
            gene: "CYP2C9".to_owned(),
            diplotype: Some("*1".to_owned()),
            phenotype: None,
            activity_score: None,
            haplotypes: BTreeSet::new(),
            warnings: Vec::new(),
        };
        let haploid_pair = OutsideCall {
            gene: "MT-RNR1".to_owned(),
            diplotype: Some("m.1555A>G/Reference".to_owned()),
            phenotype: None,
            activity_score: None,
            haplotypes: BTreeSet::new(),
            warnings: Vec::new(),
        };

        assert_eq!(
            single
                .to_annotation_input()
                .expect_err("bad single")
                .to_string(),
            "Expected two genotypes separated by a '/' but saw [*1] for CYP2C9"
        );
        assert_eq!(
            haploid_pair
                .to_annotation_input()
                .expect_err("bad haploid pair")
                .to_string(),
            "Cannot specify two genotypes [m.1555A>G/Reference] for single chromosome gene MT-RNR1"
        );
    }

    #[test]
    fn reads_gene_phenotype_json_like_java_model() {
        let phenotype = read_gene_phenotype_file(Path::new(
            "../../repos/PharmCAT/src/main/resources/org/pharmgkb/pharmcat/phenotype/DPYD.json",
        ))
        .expect("DPYD phenotype JSON");

        assert_eq!(phenotype.gene, "DPYD");
        assert!(phenotype.is_activity_gene());
        assert_eq!(phenotype.haplotype_function("Reference"), "Normal function");
        assert_eq!(phenotype.haplotype_activity("c.61C>T"), Some("0.0"));
        assert_eq!(
            phenotype.haplotype_function("missing"),
            GenePhenotype::UNASSIGNED_FUNCTION
        );
    }

    #[test]
    fn infers_dpyd_lowest_function_diplotypes_like_java_caller() {
        let dpyd = read_dpyd_phenotype();

        let inferred = dpyd
            .infer_dpyd_lowest_function_from_diplotypes([
                "c.2657G>A (*9B)/c.1679T>G (*13)",
                "c.1774C>T/c.1679T>G (*13)",
            ])
            .expect("inferred DPYD");
        assert_eq!(inferred.allele1.as_deref(), Some("c.1774C>T"));
        assert_eq!(inferred.allele2.as_deref(), Some("c.1679T>G (*13)"));

        let inferred = dpyd
            .infer_dpyd_lowest_function_from_diplotypes([
                "c.2657G>A (*9B)/[c.1774C>T + c.1679T>G (*13)]",
            ])
            .expect("inferred DPYD");
        assert_eq!(inferred.allele1.as_deref(), Some("c.2657G>A (*9B)"));
        assert_eq!(inferred.allele2.as_deref(), Some("c.1679T>G (*13)"));

        let inferred = dpyd
            .infer_dpyd_lowest_function_from_diplotypes([
                "c.498G>A/[c.2933A>G + c.1905+1G>A (*2A)]",
            ])
            .expect("inferred DPYD");
        assert_eq!(inferred.allele1.as_deref(), Some("c.498G>A"));
        assert_eq!(inferred.allele2.as_deref(), Some("c.1905+1G>A (*2A)"));
    }

    #[test]
    fn infers_dpyd_lowest_function_haplotypes_like_java_caller() {
        let dpyd = read_dpyd_phenotype();

        let inferred = dpyd
            .infer_dpyd_lowest_function_from_haplotypes([
                "c.498G>A",
                "c.2582A>G",
                "c.2846A>T",
                "c.2933A>G",
            ])
            .expect("inferred DPYD");

        assert_eq!(inferred.allele1.as_deref(), Some("c.2933A>G"));
        assert_eq!(inferred.allele2.as_deref(), Some("c.2846A>T"));
    }

    #[test]
    fn infers_ryr1_lowest_function_diplotypes_like_java_caller() {
        let ryr1 = read_ryr1_phenotype();

        let inferred = ryr1
            .infer_ryr1_lowest_function_from_diplotypes(["c.38T>G/c.51_53del"])
            .expect("inferred RYR1");
        assert_eq!(inferred.allele1.as_deref(), Some("c.38T>G"));
        assert_eq!(inferred.allele2.as_deref(), Some("c.51_53del"));

        let inferred = ryr1
            .infer_ryr1_lowest_function_from_diplotypes([
                "c.38T>G/c.51_53del",
                "c.51_53del/c.51_53del",
            ])
            .expect("inferred RYR1");
        assert_eq!(inferred.allele1.as_deref(), Some("c.38T>G"));
        assert_eq!(inferred.allele2.as_deref(), Some("c.51_53del"));

        let inferred = ryr1
            .infer_ryr1_lowest_function_from_diplotypes([
                "c.38T>G/c.38T>G",
                "c.51_53del/c.51_53del",
            ])
            .expect("inferred RYR1");
        assert_eq!(inferred.allele1.as_deref(), Some("c.38T>G"));
        assert_eq!(inferred.allele2.as_deref(), Some("c.38T>G"));

        let inferred = ryr1
            .infer_ryr1_lowest_function_from_diplotypes(["c.38T>G/c.38T>G", "c.38T>G/c.51_53del"])
            .expect("inferred RYR1");
        assert_eq!(inferred.allele1.as_deref(), Some("c.38T>G"));
        assert_eq!(inferred.allele2.as_deref(), Some("c.38T>G"));
    }

    #[test]
    fn infers_ryr1_lowest_function_haplotypes_like_java_caller() {
        let ryr1 = read_ryr1_phenotype();

        let inferred = ryr1
            .infer_ryr1_lowest_function_from_haplotypes(["c.152C>A", "c.97A>G", "c.418G>A"])
            .expect("inferred RYR1");

        assert_eq!(inferred.allele1.as_deref(), Some("c.97A>G"));
        assert_eq!(inferred.allele2.as_deref(), Some("c.152C>A"));

        let inferred = ryr1
            .infer_ryr1_lowest_function_from_haplotypes(["[c.418G>A + c.14918C>T]", "c.152C>A"])
            .expect("inferred RYR1");

        assert_eq!(inferred.allele1.as_deref(), Some("c.152C>A"));
        assert_eq!(inferred.allele2.as_deref(), Some("c.14918C>T"));
    }

    #[test]
    fn loads_phenotype_map_directory_like_java_phenotype_map() {
        let phenotype_map = PhenotypeMap::from_dir(Path::new(
            "../../repos/PharmCAT/src/main/resources/org/pharmgkb/pharmcat/phenotype",
        ))
        .expect("phenotype map");

        let genes = phenotype_map.genes().collect::<Vec<_>>();
        assert_eq!(
            genes,
            [
                "ABCG2", "CACNA1S", "CFTR", "CYP2B6", "CYP2C19", "CYP2C9", "CYP2D6", "CYP3A4",
                "CYP3A5", "DPYD", "F2", "F5", "G6PD", "IFNL3", "MT-RNR1", "NAT2", "NUDT15", "RYR1",
                "SLCO1B1", "TPMT", "UGT1A1", "VKORC1",
            ]
        );
        assert!(genes.iter().all(|gene| !gene.starts_with("HLA")));
        assert_eq!(
            phenotype_map
                .phenotype("CYP2C9")
                .expect("CYP2C9")
                .haplotypes
                .get("*6")
                .map(String::as_str),
            Some("0.0")
        );
        assert_eq!(
            phenotype_map
                .phenotype("DPYD")
                .expect("DPYD")
                .haplotype_function("Reference"),
            "Normal function"
        );
    }

    #[test]
    fn finds_diplotype_records_by_key_like_java_gene_phenotype() {
        let cyp2c19 = read_test_phenotype("CYP2C19");
        let key = [("*1".to_owned(), 2)]
            .into_iter()
            .collect::<BTreeMap<_, _>>();

        let diplotype = cyp2c19.find_diplotype(&key).expect("*1/*1");

        assert_eq!(diplotype.diplotype, "*1/*1");
        assert_eq!(diplotype.lookup_key, "Normal Metabolizer");
        assert_eq!(diplotype.phenotype, "Normal Metabolizer");
    }

    #[test]
    fn looks_up_non_activity_phenotype_by_diplotype_key_like_java() {
        let cyp2c19 = read_test_phenotype("CYP2C19");
        let key = [("*1".to_owned(), 2)]
            .into_iter()
            .collect::<BTreeMap<_, _>>();

        assert_eq!(
            cyp2c19
                .lookup_phenotype_by_diplotype_key(&key, false, false)
                .expect("phenotype"),
            "Normal Metabolizer"
        );
        assert_eq!(
            cyp2c19
                .lookup_activity_by_diplotype_key(&key, false)
                .expect("activity score"),
            NA
        );
    }

    #[test]
    fn looks_up_activity_score_gene_by_diplotype_and_activity_score_like_java() {
        let cyp2d6 = read_test_phenotype("CYP2D6");
        let key = [("*1".to_owned(), 1), ("*3".to_owned(), 1)]
            .into_iter()
            .collect::<BTreeMap<_, _>>();

        assert_eq!(
            cyp2d6
                .lookup_activity_by_diplotype_key(&key, false)
                .expect("activity score"),
            "1.0"
        );
        assert_eq!(
            cyp2d6
                .lookup_phenotype_by_diplotype_key(&key, false, true)
                .expect("phenotype"),
            "Intermediate Metabolizer"
        );
        assert_eq!(
            cyp2d6.lookup_phenotypes_by_activity_score("1.0"),
            ["Intermediate Metabolizer".to_owned()]
                .into_iter()
                .collect()
        );
    }

    #[test]
    fn looks_up_activity_scores_by_phenotype_like_java() {
        let cyp2d6 = read_test_phenotype("CYP2D6");

        assert_eq!(
            cyp2d6.lookup_activity_scores_by_phenotype(["Intermediate Metabolizer"]),
            ["0.25", "0.5", "0.75", "1.0"]
                .into_iter()
                .map(str::to_owned)
                .collect()
        );
        assert_eq!(
            cyp2d6.lookup_activity_scores_by_phenotype(["Missing phenotype"]),
            [NA.to_owned()].into_iter().collect()
        );
    }

    #[test]
    fn handles_unknown_and_missing_diplotype_lookup_like_java() {
        let cyp2c19 = read_test_phenotype("CYP2C19");
        let cyp2d6 = read_test_phenotype("CYP2D6");
        let missing_key = [("*missing".to_owned(), 2)]
            .into_iter()
            .collect::<BTreeMap<_, _>>();

        assert_eq!(
            cyp2c19
                .lookup_phenotype_by_diplotype_key(&missing_key, false, false)
                .expect("missing non-activity phenotype"),
            NA
        );
        assert_eq!(
            cyp2d6
                .lookup_phenotype_by_diplotype_key(&missing_key, false, true)
                .expect("missing activity phenotype"),
            INDETERMINATE
        );
        assert_eq!(
            cyp2d6
                .lookup_activity_by_diplotype_key(&missing_key, false)
                .expect("missing activity score"),
            NA
        );
        assert_eq!(
            cyp2d6
                .lookup_phenotype_by_diplotype_key(&missing_key, true, true)
                .expect("unknown alleles"),
            NO_RESULT
        );
        assert_eq!(
            cyp2d6.lookup_phenotypes_by_activity_score("missing"),
            [INDETERMINATE.to_owned()].into_iter().collect()
        );
    }

    #[test]
    fn annotates_called_non_activity_diplotype_like_java() {
        let cyp2c19 = read_test_phenotype("CYP2C19");

        let diplotype = DiplotypeAnnotationInput::from_alleles("CYP2C19", "*1", Some("*1"))
            .annotate(Some(&cyp2c19))
            .expect("annotated diplotype");

        assert_eq!(diplotype.phenotypes, ["Normal Metabolizer"]);
        assert_eq!(diplotype.lookup_keys, ["Normal Metabolizer"]);
        assert_eq!(diplotype.activity_score, None);
    }

    #[test]
    fn annotates_called_activity_score_diplotypes_like_java() {
        let cyp2d6 = read_test_phenotype("CYP2D6");
        let dpyd = read_test_phenotype("DPYD");

        let cyp2d6_diplotype = DiplotypeAnnotationInput::from_alleles("CYP2D6", "*1", Some("*3"))
            .annotate(Some(&cyp2d6))
            .expect("CYP2D6 annotation");
        assert_eq!(cyp2d6_diplotype.activity_score.as_deref(), Some("1.0"));
        assert_eq!(cyp2d6_diplotype.phenotypes, ["Intermediate Metabolizer"]);
        assert_eq!(cyp2d6_diplotype.lookup_keys, ["1.0"]);

        let dpyd_diplotype =
            DiplotypeAnnotationInput::from_alleles("DPYD", "Reference", Some("c.2846A>T"))
                .annotate(Some(&dpyd))
                .expect("DPYD annotation");
        assert_eq!(dpyd_diplotype.activity_score.as_deref(), Some("1.5"));
        assert_eq!(dpyd_diplotype.lookup_keys, ["1.5"]);

        let missing_dpyd = DiplotypeAnnotationInput::from_alleles("DPYD", "foo", Some("bar"))
            .annotate(Some(&dpyd))
            .expect("missing DPYD annotation");
        assert_eq!(missing_dpyd.activity_score.as_deref(), Some(NA));
        assert_eq!(missing_dpyd.lookup_keys, [NA]);
    }

    #[test]
    fn annotates_outside_activity_score_calls_like_java_diplotype_test() {
        let cyp2d6 = read_test_phenotype("CYP2D6");

        let phenotype_and_score = DiplotypeAnnotationInput::outside_call(
            "CYP2D6",
            None::<String>,
            None::<String>,
            Some("Normal Metabolizer"),
            Some("4.0"),
        )
        .annotate(Some(&cyp2d6))
        .expect("outside phenotype and score");
        assert_eq!(phenotype_and_score.phenotypes, ["Normal Metabolizer"]);
        assert_eq!(phenotype_and_score.activity_score.as_deref(), Some("4.0"));
        assert_eq!(phenotype_and_score.lookup_keys, ["4.0"]);
        assert!(phenotype_and_score.outside_phenotype_mismatch.is_some());
        assert!(
            phenotype_and_score
                .outside_activity_score_mismatch
                .is_some()
        );

        let phenotype_only = DiplotypeAnnotationInput::outside_call(
            "CYP2D6",
            None::<String>,
            None::<String>,
            Some("Normal Metabolizer"),
            None::<String>,
        )
        .annotate(Some(&cyp2d6))
        .expect("outside phenotype only");
        assert_eq!(phenotype_only.phenotypes, ["Normal Metabolizer"]);
        assert_eq!(phenotype_only.activity_score.as_deref(), Some(NA));
        assert!(phenotype_only.lookup_keys.len() > 1);
        assert_eq!(phenotype_only.outside_phenotype_mismatch, None);
        assert_eq!(phenotype_only.outside_activity_score_mismatch, None);

        let score_only = DiplotypeAnnotationInput::outside_call(
            "CYP2D6",
            None::<String>,
            None::<String>,
            None::<String>,
            Some("4.0"),
        )
        .annotate(Some(&cyp2d6))
        .expect("outside score only");
        assert_eq!(score_only.phenotypes, ["Ultrarapid Metabolizer"]);
        assert_eq!(score_only.activity_score.as_deref(), Some("4.0"));
        assert_eq!(score_only.lookup_keys, ["4.0"]);
        assert_eq!(score_only.outside_phenotype_mismatch, None);
        assert_eq!(score_only.outside_activity_score_mismatch, None);
    }

    #[test]
    fn annotates_outside_allele_calls_with_mismatches_like_java_diplotype_test() {
        let cyp2d6 = read_test_phenotype("CYP2D6");

        let phenotype_and_score = DiplotypeAnnotationInput::outside_call(
            "CYP2D6",
            Some("*1"),
            Some("*1"),
            Some("Normal Metabolizer"),
            Some("4.0"),
        )
        .annotate(Some(&cyp2d6))
        .expect("outside diplotype phenotype and score");
        assert_eq!(phenotype_and_score.phenotypes, ["Normal Metabolizer"]);
        assert_eq!(phenotype_and_score.activity_score.as_deref(), Some("4.0"));
        assert_eq!(phenotype_and_score.lookup_keys, ["4.0"]);
        assert!(phenotype_and_score.outside_phenotype_mismatch.is_some());
        assert!(
            phenotype_and_score
                .outside_activity_score_mismatch
                .is_some()
        );

        let phenotype_only = DiplotypeAnnotationInput::outside_call(
            "CYP2D6",
            Some("*1"),
            Some("*1"),
            Some("Normal Metabolizer"),
            None::<String>,
        )
        .annotate(Some(&cyp2d6))
        .expect("outside diplotype phenotype only");
        assert_eq!(phenotype_only.phenotypes, ["Normal Metabolizer"]);
        assert_eq!(phenotype_only.activity_score.as_deref(), Some(NA));
        assert!(phenotype_only.lookup_keys.len() > 1);
        assert_eq!(phenotype_only.outside_phenotype_mismatch, None);
        assert_eq!(phenotype_only.outside_activity_score_mismatch, None);

        let score_only = DiplotypeAnnotationInput::outside_call(
            "CYP2D6",
            Some("*1"),
            Some("*1"),
            None::<String>,
            Some("4.0"),
        )
        .annotate(Some(&cyp2d6))
        .expect("outside diplotype score only");
        assert_eq!(score_only.phenotypes, ["Ultrarapid Metabolizer"]);
        assert_eq!(score_only.activity_score.as_deref(), Some("4.0"));
        assert_eq!(score_only.lookup_keys, ["4.0"]);
        assert_eq!(score_only.outside_phenotype_mismatch, None);
        assert!(score_only.outside_activity_score_mismatch.is_some());
    }

    #[test]
    fn annotates_hla_a_and_hla_b_allele_presence_like_java() {
        let hla_a = DiplotypeAnnotationInput::from_alleles("HLA-A", "*31:01", Some("*02:01"))
            .annotate(None)
            .expect("HLA-A annotation");
        assert_eq!(hla_a.phenotypes, ["*31:01 positive"]);
        assert_eq!(hla_a.lookup_keys, hla_a.phenotypes);

        let hla_b = DiplotypeAnnotationInput::from_alleles("HLA-B", "*15:02", Some("*57:01"))
            .annotate(None)
            .expect("HLA-B annotation");
        assert_eq!(
            hla_b.phenotypes,
            ["*15:02 positive", "*57:01 positive", "*58:01 negative"]
        );
        assert_eq!(hla_b.lookup_keys, hla_b.phenotypes);
    }

    #[test]
    fn annotates_hla_unknown_and_outside_phenotype_like_java() {
        let unknown = DiplotypeAnnotationInput::from_alleles("HLA-B", "Unknown", Some("Unknown"))
            .annotate(None)
            .expect("unknown HLA-B annotation");
        assert!(unknown.phenotypes.is_empty());
        assert_eq!(unknown.lookup_keys, [NO_RESULT]);

        let phenotype_only = DiplotypeAnnotationInput::outside_call(
            "HLA-B",
            None::<String>,
            None::<String>,
            Some("*57:01 positive"),
            None::<String>,
        )
        .annotate(None)
        .expect("phenotype-only HLA-B annotation");
        assert_eq!(phenotype_only.phenotypes, ["*57:01 positive"]);
        assert_eq!(phenotype_only.lookup_keys, ["*57:01 positive"]);

        let known_with_outside_phenotype = DiplotypeAnnotationInput::outside_call(
            "HLA-B",
            Some("*15:02"),
            Some("*57:01"),
            Some("*57:01 positive"),
            None::<String>,
        )
        .annotate(None)
        .expect("known HLA-B with outside phenotype");
        assert_eq!(known_with_outside_phenotype.phenotypes, ["*57:01 positive"]);
        assert!(known_with_outside_phenotype.lookup_keys.is_empty());
    }

    #[test]
    fn rejects_unsupported_hla_gene_like_java_factory() {
        let error = DiplotypeAnnotationInput::from_alleles("HLA-C", "*01:01", Some("*02:01"))
            .annotate(None)
            .expect_err("unsupported HLA gene");

        assert_eq!(
            error.to_string(),
            "Gene not supported for HLA phenotype calling: HLA-C"
        );
    }

    fn outside_call_validation() -> OutsideCallValidation {
        let mut validation =
            OutsideCallValidation::for_supported_genes(["CYP2C9", "CYP2C19", "CYP2D6", "HLA-A"]);
        validation.activity_score_genes.insert("CYP2D6".to_owned());
        validation.valid_named_alleles.insert(
            "CYP2C9".to_owned(),
            BTreeSet::from([
                "*1".to_owned(),
                "*2".to_owned(),
                "*3".to_owned(),
                "*4".to_owned(),
            ]),
        );
        validation.valid_named_alleles.insert(
            "CYP2C19".to_owned(),
            BTreeSet::from(["*3".to_owned(), "*4".to_owned()]),
        );
        validation.valid_named_alleles.insert(
            "CYP2D6".to_owned(),
            BTreeSet::from(["*1".to_owned(), "*4".to_owned()]),
        );
        validation
            .valid_named_alleles
            .insert("HLA-A".to_owned(), BTreeSet::from(["*31:01".to_owned()]));
        validation
    }

    fn read_test_phenotype(gene: &str) -> GenePhenotype {
        let path = format!(
            "../../repos/PharmCAT/src/main/resources/org/pharmgkb/pharmcat/phenotype/{gene}.json"
        );
        read_gene_phenotype_file(Path::new(&path)).expect("gene phenotype JSON")
    }

    fn read_ryr1_phenotype() -> GenePhenotype {
        let functions = [
            ("c.38T>G", "Malignant Hyperthermia associated"),
            ("c.51_53del", "Uncertain function"),
            ("c.97A>G", "Malignant Hyperthermia associated"),
            ("c.152C>A", "Uncertain function"),
            ("c.418G>A", "Normal function"),
            ("c.14918C>T", "Malignant Hyperthermia associated"),
        ];

        GenePhenotype {
            gene: "RYR1".to_owned(),
            haplotypes: functions
                .iter()
                .map(|(name, function)| ((*name).to_owned(), (*function).to_owned()))
                .collect(),
            activity_values: BTreeMap::new(),
            diplotypes: Vec::new(),
            named_alleles: functions
                .into_iter()
                .map(|(name, function)| HaplotypeRecord {
                    name: name.to_owned(),
                    activity_value: None,
                    function_value: Some(function.to_owned()),
                    lookup_key: name.to_owned(),
                })
                .collect(),
            version: None,
        }
    }

    fn read_dpyd_phenotype() -> GenePhenotype {
        read_test_phenotype("DPYD")
    }

    fn write_temp_file(prefix: &str, contents: &str) -> PathBuf {
        let path = temp_path(prefix);
        fs::write(&path, contents).expect("write temp file");
        path
    }

    fn temp_path(prefix: &str) -> PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("time")
            .as_nanos();
        std::env::temp_dir().join(format!("{prefix}-{nanos}.tsv"))
    }
}
