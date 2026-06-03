//! Named allele matching primitives.

use std::collections::{BTreeMap, BTreeSet};

use crate::{
    common::chromosome,
    definition::{DefinitionExemption, DefinitionFile, NamedAllele, VariantLocus},
    vcf::SampleAlleleSummary,
};

/// Sample alleles prepared for matching one gene.
#[derive(Clone, Debug, PartialEq)]
pub struct MatchData {
    /// Selected sample ID.
    pub sample_id: String,
    /// Gene being matched.
    pub gene: String,
    /// Definition positions that have sample data.
    pub positions: Vec<VariantLocus>,
    /// Definition positions missing from the sample.
    pub missing_positions: BTreeSet<VariantLocus>,
    /// Definition positions whose sample calls contain undocumented variations.
    pub positions_with_undocumented_variations: BTreeSet<VariantLocus>,
    /// Whether undocumented variation calls were treated as reference during matching.
    pub treat_undocumented_variations_as_reference: bool,
    /// Missing required positions from definition exemptions.
    pub missing_required_positions: Vec<String>,
    /// Missing AMP Tier 1 positions from definition exemptions.
    pub missing_amp1_positions: Vec<String>,
    sample_alleles: Vec<SampleAllele>,
    haplotypes: Vec<NamedAllele>,
    permutations: BTreeSet<String>,
    /// Whether every sample allele was phased in the VCF.
    pub phased: bool,
    /// Whether Java would treat this data as effectively phased.
    pub effectively_phased: bool,
    /// Whether the available sample alleles are haploid.
    pub haploid: bool,
    /// Whether all sample alleles are homozygous or haploid.
    pub homozygous: bool,
}

impl MatchData {
    /// Builds match data from definition positions and VCF sample allele calls.
    pub fn new(
        sample_id: impl Into<String>,
        gene: impl Into<String>,
        definition: &DefinitionFile,
        allele_map: &BTreeMap<String, SampleAlleleSummary>,
    ) -> Self {
        Self::new_with_exemption(sample_id, gene, definition, None, allele_map)
    }

    /// Builds match data and applies exemption-driven missing-position tracking.
    pub fn new_with_exemption(
        sample_id: impl Into<String>,
        gene: impl Into<String>,
        definition: &DefinitionFile,
        exemption: Option<&DefinitionExemption>,
        allele_map: &BTreeMap<String, SampleAlleleSummary>,
    ) -> Self {
        let mut positions = Vec::new();
        let mut missing_positions = BTreeSet::new();
        let mut positions_with_undocumented_variations = BTreeSet::new();
        let mut treat_undocumented_variations_as_reference = false;
        let mut sample_alleles = Vec::new();
        let mut phased = true;

        for variant in &definition.variants {
            let chr_position = variant.vcf_chr_position();
            let Some(call) = allele_map.get(&chr_position) else {
                missing_positions.insert(variant.clone());
                continue;
            };

            positions.push(variant.clone());
            if !call.undocumented_variations.is_empty() {
                positions_with_undocumented_variations.insert(variant.clone());
                if call.treat_undocumented_variations_as_reference {
                    treat_undocumented_variations_as_reference = true;
                }
            }
            phased &= call.phased;
            sample_alleles.push(SampleAllele::from_summary(call));
        }

        sample_alleles.sort();
        let haploid = are_sample_alleles_haploid(&sample_alleles);
        let homozygous = haploid || sample_alleles.iter().all(SampleAllele::is_homozygous);
        let missing_required_positions = exemption
            .map(|exemption| {
                missing_positions
                    .iter()
                    .filter(|missing| exemption.required_positions.contains(&missing.position))
                    .map(VariantLocus::vcf_chr_position)
                    .collect::<Vec<_>>()
            })
            .unwrap_or_default();
        let missing_amp1_positions = exemption
            .map(|exemption| {
                missing_positions
                    .iter()
                    .filter(|missing| exemption.amp1_positions.contains(&missing.position))
                    .map(VariantLocus::vcf_chr_position)
                    .collect::<Vec<_>>()
            })
            .unwrap_or_default();

        Self {
            sample_id: sample_id.into(),
            gene: gene.into(),
            positions,
            missing_positions,
            positions_with_undocumented_variations,
            treat_undocumented_variations_as_reference,
            missing_required_positions,
            missing_amp1_positions,
            sample_alleles,
            haplotypes: Vec::new(),
            permutations: BTreeSet::new(),
            phased,
            effectively_phased: false,
            haploid,
            homozygous,
        }
    }

    /// Organizes named alleles against available sample positions.
    pub fn marshall_haplotypes(&mut self, definition: &DefinitionFile) {
        if self.missing_positions.is_empty() {
            self.haplotypes = definition.named_alleles.clone();
            return;
        }

        let mut haplotypes = Vec::new();
        for haplotype in &definition.named_alleles {
            let mut alleles = Vec::with_capacity(self.positions.len());
            let mut cpic_alleles = Vec::with_capacity(self.positions.len());
            let mut score = 0;

            for variant in &self.positions {
                let Some(index) = definition.index_for_position(variant.position) else {
                    continue;
                };
                let allele = haplotype.alleles.get(index).cloned().unwrap_or_default();
                if allele.is_some() {
                    score += 1;
                }
                alleles.push(allele);
                cpic_alleles.push(
                    haplotype
                        .cpic_alleles
                        .get(index)
                        .cloned()
                        .unwrap_or_default(),
                );
            }

            if score > 0 {
                let mut available = haplotype.clone();
                available.alleles = alleles;
                available.cpic_alleles = cpic_alleles;
                available.core_positions = self
                    .positions
                    .iter()
                    .filter_map(|position| {
                        haplotype
                            .core_positions
                            .contains(&position.position)
                            .then_some(position.position)
                    })
                    .collect();
                available.missing_positions = self
                    .missing_positions
                    .iter()
                    .filter(|missing| {
                        definition
                            .index_for_position(missing.position)
                            .and_then(|index| haplotype.alleles.get(index))
                            .is_some_and(Option::is_some)
                    })
                    .cloned()
                    .collect();
                available.score_override = Some(score - haplotype.num_partials);
                haplotypes.push(available);
            }
        }

        self.haplotypes = haplotypes;
    }

    /// Fills missing named-allele positions with reference alleles like Java `defaultMissingAllelesToReference`.
    pub fn default_missing_alleles_to_reference(&mut self) -> Result<(), MatchError> {
        let reference = self
            .haplotypes
            .iter()
            .find(|haplotype| haplotype.reference)
            .cloned()
            .ok_or_else(|| MatchError::NoReference(self.gene.clone()))?;

        for haplotype in &mut self.haplotypes {
            if haplotype.reference {
                continue;
            }

            let original_score = named_allele_score(haplotype);
            for index in 0..reference.alleles.len() {
                if haplotype.alleles.get(index).is_some_and(Option::is_none) {
                    let reference_allele =
                        reference.alleles.get(index).cloned().unwrap_or_default();
                    haplotype.alleles[index] =
                        if reference_allele.as_deref().is_some_and(is_iupac_wobble) {
                            self.positions
                                .get(index)
                                .map(|position| Some(position.reference.clone()))
                                .unwrap_or(reference_allele)
                        } else {
                            reference_allele
                        };
                    haplotype.cpic_alleles[index] = reference
                        .cpic_alleles
                        .get(index)
                        .cloned()
                        .unwrap_or_default();
                }
            }
            haplotype.score_override = Some(original_score);
        }

        Ok(())
    }

    /// Generates sample allele permutations with Java's phasing rules.
    pub fn generate_sample_permutations(&mut self) -> Result<(), MatchError> {
        if self.sample_alleles.is_empty() {
            return Err(MatchError::NoSampleAlleles);
        }

        self.permutations = generate_permutations(&self.sample_alleles)?;
        self.effectively_phased = self.permutations.len() <= 2;
        Ok(())
    }

    /// Returns generated sample permutations.
    pub fn permutations(&self) -> &BTreeSet<String> {
        &self.permutations
    }

    /// Returns the named alleles that remain matchable after marshalling missing positions.
    pub fn haplotypes(&self) -> &[NamedAllele] {
        &self.haplotypes
    }

    /// Returns whether any sample allele has a Java-retained phase set.
    pub fn is_using_phase_sets(&self) -> bool {
        self.sample_alleles
            .iter()
            .any(|sample_allele| sample_allele.phase_set.is_some())
    }

    /// Returns the sample allele retained for a definition position.
    pub fn sample_allele_at_position(&self, position: u64) -> Option<&SampleAllele> {
        self.sample_alleles
            .iter()
            .find(|sample_allele| sample_allele.position as u64 == position)
    }

    /// Returns whether any sample allele has a missing side in the VCF call.
    pub fn has_partial_missing_alleles(&self) -> bool {
        self.sample_alleles.iter().any(|sample_allele| {
            sample_allele.allele1.as_deref() == Some(".")
                || sample_allele.allele2.as_deref() == Some(".")
        })
    }

    /// Returns the Java-retained phase set for a definition position.
    pub fn phase_set(&self, position: u64) -> Option<i32> {
        self.sample_alleles
            .iter()
            .find(|sample_allele| sample_allele.position as u64 == position)
            .and_then(|sample_allele| sample_allele.phase_set)
    }

    /// Compares generated sample permutations against marshalled haplotypes.
    pub fn compare_permutations(&self) -> Vec<HaplotypeMatch> {
        let mut matches = Vec::new();

        for haplotype in &self.haplotypes {
            let sequences = self
                .permutations
                .iter()
                .filter(|sequence| haplotype_matches_sequence(haplotype, &self.positions, sequence))
                .cloned()
                .collect::<BTreeSet<_>>();

            if !sequences.is_empty() {
                matches.push(HaplotypeMatch {
                    name: haplotype.name.clone(),
                    haplotype: haplotype.clone(),
                    positions: self.positions.clone(),
                    sequences,
                });
            }
        }

        matches.sort_by(compare_haplotype_matches);
        matches
    }

    /// Computes diplotype matches from currently matched haplotypes.
    pub fn compute_diplotypes(&self, top_candidate_only: bool) -> Vec<DiplotypeMatch> {
        let haplotype_matches = self.compare_permutations();
        self.compute_diplotypes_from_matches(haplotype_matches, top_candidate_only)
    }

    /// Computes diplotypes using Java `CombinationMatcher` strand matching.
    pub fn compute_combination_diplotypes(
        &self,
        definition: &DefinitionFile,
        find_partials: bool,
    ) -> Vec<DiplotypeMatch> {
        let combination_matches = self.compute_combination_matches(definition, find_partials);
        self.compute_diplotypes_from_matches(combination_matches, false)
    }

    /// Computes Java `CombinationMatcher` strand matches for current sample permutations.
    pub fn compute_combination_matches(
        &self,
        definition: &DefinitionFile,
        find_partials: bool,
    ) -> Vec<HaplotypeMatch> {
        let mut matches = Vec::new();

        for sequence in &self.permutations {
            let allele_map = parse_sequence(sequence);
            let var_positions = self
                .positions
                .iter()
                .filter_map(|variant| {
                    let allele = allele_map.get(&variant.position)?;
                    if *allele != variant.reference {
                        Some(variant.position)
                    } else {
                        None
                    }
                })
                .collect::<BTreeSet<_>>();

            let covered_haps = self
                .haplotypes
                .iter()
                .filter(|haplotype| {
                    !haplotype.reference
                        && sample_has_named_allele(&allele_map, haplotype, &self.positions)
                })
                .cloned()
                .collect::<Vec<_>>();

            if covered_haps.len() <= 1 {
                let haplotypes = if covered_haps.is_empty() {
                    self.haplotypes
                        .iter()
                        .find(|haplotype| haplotype.reference)
                        .cloned()
                        .into_iter()
                        .collect::<Vec<_>>()
                } else {
                    covered_haps
                };
                let partials = calculate_partial_names(
                    definition,
                    &allele_map,
                    &var_positions,
                    &haplotypes,
                    find_partials,
                );
                if partials.is_empty() {
                    if let Some(haplotype) = haplotypes.first() {
                        matches.push(HaplotypeMatch::from_haplotype(
                            haplotype.clone(),
                            self.positions.clone(),
                            sequence.clone(),
                        ));
                    }
                } else {
                    matches.push(build_combination_match(
                        &self.positions,
                        sequence,
                        &haplotypes,
                        partials,
                    ));
                }
            } else {
                for combo in compute_viable_combinations(covered_haps) {
                    let partials = calculate_partial_names(
                        definition,
                        &allele_map,
                        &var_positions,
                        &combo,
                        find_partials,
                    );
                    matches.push(build_combination_match(
                        &self.positions,
                        sequence,
                        &combo,
                        partials,
                    ));
                }
            }
        }

        matches.sort_by(compare_haplotype_matches);
        matches.dedup();
        matches
    }

    fn compute_diplotypes_from_matches(
        &self,
        haplotype_matches: Vec<HaplotypeMatch>,
        top_candidate_only: bool,
    ) -> Vec<DiplotypeMatch> {
        if haplotype_matches.is_empty() {
            return Vec::new();
        }

        let mut diplotypes = if self.permutations.len() == 1 {
            self.determine_homozygous_pairs(&haplotype_matches)
        } else {
            self.determine_heterozygous_pairs(&haplotype_matches)
        };

        diplotypes.sort();
        diplotypes.dedup();

        if top_candidate_only && diplotypes.len() > 1 {
            let top_score = diplotypes[0].score;
            diplotypes.retain(|diplotype| diplotype.score == top_score);
        }

        diplotypes
    }

    fn determine_homozygous_pairs(
        &self,
        haplotype_matches: &[HaplotypeMatch],
    ) -> Vec<DiplotypeMatch> {
        let seq = self
            .permutations
            .iter()
            .next()
            .expect("permutation exists")
            .clone();

        if haplotype_matches.len() == 1 {
            let haplotype1 = haplotype_matches[0].clone();
            let haplotype2 = if self.haploid {
                None
            } else {
                Some(haplotype1.clone())
            };
            let sequence_pair = if self.haploid {
                vec![seq]
            } else {
                vec![seq.clone(), seq]
            };
            return vec![DiplotypeMatch::new(haplotype1, haplotype2, sequence_pair)];
        }

        perfect_pairs(haplotype_matches)
            .into_iter()
            .map(|(haplotype1, haplotype2)| {
                DiplotypeMatch::new(
                    haplotype1.clone(),
                    Some(haplotype2),
                    vec![seq.clone(), seq.clone()],
                )
            })
            .collect()
    }

    fn determine_heterozygous_pairs(
        &self,
        haplotype_matches: &[HaplotypeMatch],
    ) -> Vec<DiplotypeMatch> {
        let mut by_name: BTreeMap<String, Vec<HaplotypeMatch>> = BTreeMap::new();
        for haplotype_match in haplotype_matches {
            by_name
                .entry(haplotype_match.name.clone())
                .or_default()
                .push(haplotype_match.clone());
        }

        let mut names = by_name.keys().cloned().collect::<Vec<_>>();
        names.sort_by(|left, right| compare_haplotype_names(left, right));

        let mut diplotypes = Vec::new();
        for (name1, name2) in perfect_pairs(&names) {
            let hm1s = by_name.get(name1).expect("haplotype matches");
            let hm2s = by_name.get(&name2).expect("haplotype matches");

            if *name1 == name2 && hm1s.len() == 1 && hm1s[0].sequences.len() == 1 {
                continue;
            }

            for haplotype1 in hm1s {
                for haplotype2 in hm2s {
                    let sequence_pairs = self.find_sequence_pairs(haplotype1, haplotype2);
                    if !sequence_pairs.is_empty() {
                        for sequence_pair in sequence_pairs {
                            diplotypes.push(DiplotypeMatch::new(
                                haplotype1.clone(),
                                Some(haplotype2.clone()),
                                sequence_pair,
                            ));
                        }
                    }
                }
            }
        }

        diplotypes
    }

    fn find_sequence_pairs(
        &self,
        haplotype1: &HaplotypeMatch,
        haplotype2: &HaplotypeMatch,
    ) -> Vec<Vec<String>> {
        let mut sequence_pairs = Vec::new();

        for seq1 in &haplotype1.sequences {
            for seq2 in &haplotype2.sequences {
                if self.is_viable_complement(seq1, seq2) {
                    sequence_pairs.push(vec![seq1.clone(), seq2.clone()]);
                }
            }
        }

        sequence_pairs
    }

    fn is_viable_complement(&self, sequence1: &str, sequence2: &str) -> bool {
        let sequence1 = parse_sequence(sequence1);
        let sequence2 = parse_sequence(sequence2);

        for sample_allele in &self.sample_alleles {
            let position = sample_allele.position as u64;
            let Some(allele1) = sequence1.get(&position) else {
                return false;
            };
            let Some(allele2) = sequence2.get(&position) else {
                return false;
            };

            if sample_allele.is_homozygous() {
                if allele1 != allele2 {
                    return false;
                }
            } else if allele1 == allele2 {
                return false;
            }
        }

        true
    }
}

/// Calls a standard, non-lowest-function gene using Java `NamedAlleleMatcher.callAssumingReference`.
pub fn call_standard_gene(
    sample_id: impl Into<String>,
    definition: &DefinitionFile,
    allele_map: &BTreeMap<String, SampleAlleleSummary>,
    top_candidate_only: bool,
    find_combinations: bool,
) -> Result<GeneCallResult, MatchError> {
    call_standard_gene_with_exemption(
        sample_id,
        definition,
        None,
        allele_map,
        top_candidate_only,
        find_combinations,
    )
}

/// Calls a standard gene with Java `ResultBuilder.build` finalization data.
pub fn call_standard_gene_with_exemption(
    sample_id: impl Into<String>,
    definition: &DefinitionFile,
    exemption: Option<&DefinitionExemption>,
    allele_map: &BTreeMap<String, SampleAlleleSummary>,
    top_candidate_only: bool,
    find_combinations: bool,
) -> Result<GeneCallResult, MatchError> {
    let sample_id = sample_id.into();
    let gene = definition.gene_symbol.as_str();
    let data = initialize_match_data_with_exemption(
        &sample_id, gene, definition, exemption, allele_map, true,
    )?;
    if data.permutations.is_empty() {
        return finalize_gene_call_result(
            GeneCallResult::no_call(gene, data),
            definition,
            exemption,
        );
    }
    if !data.missing_required_positions.is_empty() {
        let missing = data.missing_required_positions.clone();
        let mut result = GeneCallResult::no_call(gene, data);
        result
            .warnings
            .insert(GeneCallWarning::MissingRequiredPosition(missing));
        return finalize_gene_call_result(result, definition, exemption);
    }
    if !data.positions_with_undocumented_variations.is_empty()
        && !data.treat_undocumented_variations_as_reference
        && !find_combinations
    {
        return finalize_gene_call_result(
            GeneCallResult::no_call(gene, data),
            definition,
            exemption,
        );
    }

    if data.has_partial_missing_alleles() {
        if find_combinations {
            let result = call_standard_gene_combination_with_exemption(
                &sample_id, definition, exemption, allele_map,
            )?;
            return finalize_gene_call_result(result, definition, exemption);
        }
        return finalize_gene_call_result(
            GeneCallResult::no_call(gene, data),
            definition,
            exemption,
        );
    }

    let unphased_priority_mode = !data.effectively_phased
        && exemption.is_some_and(|exemption| !exemption.unphased_diplotype_priorities.is_empty());
    let matches = data
        .compute_diplotypes(top_candidate_only && !find_combinations && !unphased_priority_mode);
    if matches.is_empty() && find_combinations {
        let result = call_standard_gene_combination_with_exemption(
            &sample_id, definition, exemption, allele_map,
        )?;
        return finalize_gene_call_result(result, definition, exemption);
    }
    if matches.is_empty() {
        return finalize_gene_call_result(
            GeneCallResult::no_call(gene, data),
            definition,
            exemption,
        );
    }

    let result = GeneCallResult::diplotypes(gene, data, matches, BTreeSet::new());
    finalize_gene_call_result(result, definition, exemption)
}

/// Applies the first Java `ResultBuilder.build` finalization behaviors.
pub fn finalize_gene_call_result(
    mut result: GeneCallResult,
    definition: &DefinitionFile,
    exemption: Option<&DefinitionExemption>,
) -> Result<GeneCallResult, MatchError> {
    if !result.match_data.missing_amp1_positions.is_empty() {
        result.warnings.insert(GeneCallWarning::MissingAmp1Position(
            result.match_data.missing_amp1_positions.clone(),
        ));
    }

    if let Some(exemption) = exemption
        && let GeneCallKind::Diplotypes(diplotypes) = &mut result.kind
        && diplotypes.len() > 1
        && !result.match_data.effectively_phased
        && !exemption.unphased_diplotype_priorities.is_empty()
    {
        let diplotype_names = diplotypes
            .iter()
            .map(|diplotype| diplotype.name.clone())
            .collect::<BTreeSet<_>>();
        for priority in &exemption.unphased_diplotype_priorities {
            if priority
                .list
                .iter()
                .all(|candidate| diplotype_names.contains(candidate))
            {
                let Some(priority_diplotype) = diplotypes
                    .iter()
                    .find(|diplotype| diplotype.name == priority.pick)
                    .cloned()
                else {
                    return Err(MatchError::MissingPriorityDiplotype(priority.pick.clone()));
                };
                *diplotypes = vec![priority_diplotype];
                result.warnings.insert(GeneCallWarning::UnphasedPriority);
                break;
            }
        }
    }

    if !definition.suballeles_map.is_empty()
        && let GeneCallKind::Diplotypes(diplotypes) = &mut result.kind
    {
        for diplotype in diplotypes {
            diplotype.handle_suballele_conversion(definition);
        }
    }

    Ok(result)
}

fn call_standard_gene_combination_with_exemption(
    sample_id: &str,
    definition: &DefinitionFile,
    exemption: Option<&DefinitionExemption>,
    allele_map: &BTreeMap<String, SampleAlleleSummary>,
) -> Result<GeneCallResult, MatchError> {
    let gene = definition.gene_symbol.as_str();
    let data = initialize_match_data_with_exemption(
        sample_id, gene, definition, exemption, allele_map, false,
    )?;
    if data.permutations.is_empty() {
        return finalize_gene_call_result(
            GeneCallResult::no_call(gene, data),
            definition,
            exemption,
        );
    }

    let matches = data.compute_combination_diplotypes(definition, true);
    if matches.is_empty() {
        return finalize_gene_call_result(
            GeneCallResult::no_call(gene, data),
            definition,
            exemption,
        );
    }
    let result = GeneCallResult::diplotypes(gene, data, matches, BTreeSet::new());
    finalize_gene_call_result(result, definition, exemption)
}

/// Computes the initial Java lowest-function DPYD diplotype path.
///
/// This covers the phased/effectively-phased branch in `NamedAlleleMatcher.callLowestFunctionGene`:
/// exact matching without HapB3 alleles, followed by Java `DpydHapB3Matcher` merge behavior.
pub fn compute_dpyd_lowest_function_diplotypes(
    sample_id: impl Into<String>,
    definition: &DefinitionFile,
    allele_map: &BTreeMap<String, SampleAlleleSummary>,
) -> Result<Vec<DiplotypeMatch>, MatchError> {
    let sample_id = sample_id.into();
    let orig_data = initialize_match_data(&sample_id, "DPYD", definition, allele_map, true)?;
    if orig_data.permutations.is_empty() {
        return Ok(Vec::new());
    }

    let mut hap_b3_matcher = DpydHapB3Matcher::new(definition, allele_map, &orig_data)?;
    let working_definition = if hap_b3_matcher.has_hap_b3_variants() {
        definition_without_dpyd_hap_b3(definition)
    } else {
        definition.clone()
    };
    let working_data =
        initialize_match_data(&sample_id, "DPYD", &working_definition, allele_map, true)?;

    if orig_data.effectively_phased {
        let diplotype_matches = working_data.compute_diplotypes(false);
        if hap_b3_matcher.has_hap_b3_variants() {
            let merger_data =
                initialize_match_data(&sample_id, "DPYD", definition, allele_map, false)?;
            if !diplotype_matches.is_empty() {
                return hap_b3_matcher.merge_phased_hap_b3_call(&merger_data, &diplotype_matches);
            }
            if !hap_b3_matcher.has_non_hap_b3_variants() {
                return hap_b3_matcher.add_phased_hap_b3_call_to_ref(&merger_data);
            }
        }
        if !diplotype_matches.is_empty() {
            return Ok(diplotype_matches);
        }
    }

    if orig_data.effectively_phased || orig_data.is_using_phase_sets() {
        let combo_definition = if hap_b3_matcher.has_hap_b3_variants() {
            definition_without_dpyd_hap_b3(definition)
        } else {
            definition.clone()
        };
        let combo_data =
            initialize_match_data(&sample_id, "DPYD", &combo_definition, allele_map, false)?;
        let mut diplotype_matches =
            combo_data.compute_combination_diplotypes(&combo_definition, true);

        if hap_b3_matcher.has_hap_b3_variants() {
            let merger_data =
                initialize_match_data(&sample_id, "DPYD", definition, allele_map, false)?;
            if !diplotype_matches.is_empty() {
                diplotype_matches =
                    hap_b3_matcher.fix_partials(&merger_data, &diplotype_matches)?;
                return hap_b3_matcher.merge_phased_hap_b3_call(&merger_data, &diplotype_matches);
            }
            if !hap_b3_matcher.has_non_hap_b3_variants() {
                return hap_b3_matcher.add_phased_hap_b3_call_to_ref(&merger_data);
            }
        }
        if !diplotype_matches.is_empty() {
            return Ok(diplotype_matches);
        }
    }

    Ok(Vec::new())
}

/// Calls DPYD through Java's lowest-function gene flow.
///
/// This is the first Rust result-level wrapper for `NamedAlleleMatcher.callLowestFunctionGene`.
/// It returns the same broad outcomes Java stores through `ResultBuilder`: no call, diplotype
/// matches, or fallback haplotype matches.
pub fn call_dpyd_lowest_function_gene(
    sample_id: impl Into<String>,
    definition: &DefinitionFile,
    allele_map: &BTreeMap<String, SampleAlleleSummary>,
) -> Result<GeneCallResult, MatchError> {
    let sample_id = sample_id.into();
    let orig_data = initialize_match_data(&sample_id, "DPYD", definition, allele_map, true)?;
    if orig_data.permutations.is_empty() {
        return Ok(GeneCallResult::no_call("DPYD", orig_data));
    }

    let mut hap_b3_matcher = DpydHapB3Matcher::new(definition, allele_map, &orig_data)?;
    let working_definition = if hap_b3_matcher.has_hap_b3_variants() {
        definition_without_dpyd_hap_b3(definition)
    } else {
        definition.clone()
    };
    let working_data =
        initialize_match_data(&sample_id, "DPYD", &working_definition, allele_map, true)?;

    if orig_data.effectively_phased {
        let diplotype_matches = working_data.compute_diplotypes(false);
        if hap_b3_matcher.has_hap_b3_variants() {
            let merger_data =
                initialize_match_data(&sample_id, "DPYD", definition, allele_map, false)?;
            if !diplotype_matches.is_empty() {
                let matches =
                    hap_b3_matcher.merge_phased_hap_b3_call(&merger_data, &diplotype_matches)?;
                let warnings = hap_b3_matcher.warnings().clone();
                return Ok(GeneCallResult::diplotypes(
                    "DPYD",
                    merger_data,
                    matches,
                    warnings,
                ));
            }
            if !hap_b3_matcher.has_non_hap_b3_variants() {
                let matches = hap_b3_matcher.add_phased_hap_b3_call_to_ref(&merger_data)?;
                let warnings = hap_b3_matcher.warnings().clone();
                return Ok(GeneCallResult::diplotypes(
                    "DPYD",
                    merger_data,
                    matches,
                    warnings,
                ));
            }
        }
        if !diplotype_matches.is_empty() {
            return Ok(GeneCallResult::diplotypes(
                "DPYD",
                working_data,
                diplotype_matches,
                BTreeSet::new(),
            ));
        }
    }

    if orig_data.effectively_phased || orig_data.is_using_phase_sets() {
        let combo_definition = if hap_b3_matcher.has_hap_b3_variants() {
            definition_without_dpyd_hap_b3(definition)
        } else {
            definition.clone()
        };
        let combo_data =
            initialize_match_data(&sample_id, "DPYD", &combo_definition, allele_map, false)?;
        let mut diplotype_matches =
            combo_data.compute_combination_diplotypes(&combo_definition, true);

        if hap_b3_matcher.has_hap_b3_variants() {
            let merger_data =
                initialize_match_data(&sample_id, "DPYD", definition, allele_map, false)?;
            if !diplotype_matches.is_empty() {
                diplotype_matches =
                    hap_b3_matcher.fix_partials(&merger_data, &diplotype_matches)?;
                let matches =
                    hap_b3_matcher.merge_phased_hap_b3_call(&merger_data, &diplotype_matches)?;
                let warnings = hap_b3_matcher.warnings().clone();
                return Ok(GeneCallResult::diplotypes(
                    "DPYD",
                    merger_data,
                    matches,
                    warnings,
                ));
            }
            if !hap_b3_matcher.has_non_hap_b3_variants() {
                let matches = hap_b3_matcher.add_phased_hap_b3_call_to_ref(&merger_data)?;
                let warnings = hap_b3_matcher.warnings().clone();
                return Ok(GeneCallResult::diplotypes(
                    "DPYD",
                    merger_data,
                    matches,
                    warnings,
                ));
            }
        }
        if !diplotype_matches.is_empty() {
            return Ok(GeneCallResult::diplotypes(
                "DPYD",
                combo_data,
                diplotype_matches,
                BTreeSet::new(),
            ));
        }
    }

    let combo_definition = if hap_b3_matcher.has_hap_b3_variants() {
        definition_without_dpyd_hap_b3(definition)
    } else {
        definition.clone()
    };
    let combo_data =
        initialize_match_data(&sample_id, "DPYD", &combo_definition, allele_map, false)?;
    let combo_dip_matches = combo_data.compute_combination_diplotypes(&combo_definition, false);
    if !combo_dip_matches.is_empty()
        && (orig_data.effectively_phased || combo_data.is_using_phase_sets())
    {
        if hap_b3_matcher.has_hap_b3_variants() {
            let merger_data =
                initialize_match_data(&sample_id, "DPYD", definition, allele_map, false)?;
            let matches =
                hap_b3_matcher.merge_phased_hap_b3_call(&merger_data, &combo_dip_matches)?;
            let warnings = hap_b3_matcher.warnings().clone();
            return Ok(GeneCallResult::diplotypes(
                "DPYD",
                merger_data,
                matches,
                warnings,
            ));
        }
        return Ok(GeneCallResult::diplotypes(
            "DPYD",
            combo_data,
            combo_dip_matches,
            BTreeSet::new(),
        ));
    }

    let has_hap_b3_variants = hap_b3_matcher.has_hap_b3_variants();
    let haplotype_matches = if has_hap_b3_variants {
        call_haplotypes_for_lowest_function_gene(
            &combo_data,
            &combo_dip_matches,
            Some(&mut hap_b3_matcher),
        )?
    } else {
        call_haplotypes_for_lowest_function_gene(&combo_data, &combo_dip_matches, None)?
    };
    let warnings = if has_hap_b3_variants {
        hap_b3_matcher.warnings().clone()
    } else {
        BTreeSet::new()
    };
    Ok(GeneCallResult::haplotypes(
        "DPYD",
        orig_data,
        haplotype_matches,
        warnings,
    ))
}

/// Calls RYR1 through Java's lowest-function gene flow.
///
/// RYR1 shares the generic `NamedAlleleMatcher.callLowestFunctionGene` branch with DPYD, but
/// without DPYD's HapB3 preprocessing and merge rules.
pub fn call_ryr1_lowest_function_gene(
    sample_id: impl Into<String>,
    definition: &DefinitionFile,
    allele_map: &BTreeMap<String, SampleAlleleSummary>,
) -> Result<GeneCallResult, MatchError> {
    call_non_dpyd_lowest_function_gene(sample_id, "RYR1", definition, None, allele_map)
}

fn call_non_dpyd_lowest_function_gene(
    sample_id: impl Into<String>,
    gene: &str,
    definition: &DefinitionFile,
    exemption: Option<&DefinitionExemption>,
    allele_map: &BTreeMap<String, SampleAlleleSummary>,
) -> Result<GeneCallResult, MatchError> {
    let sample_id = sample_id.into();
    let orig_data = initialize_match_data_with_exemption(
        &sample_id, gene, definition, exemption, allele_map, true,
    )?;
    if orig_data.permutations.is_empty() {
        return finalize_gene_call_result(
            GeneCallResult::no_call(gene, orig_data),
            definition,
            exemption,
        );
    }
    if !orig_data.missing_required_positions.is_empty() {
        let missing = orig_data.missing_required_positions.clone();
        let mut result = GeneCallResult::no_call(gene, orig_data);
        result
            .warnings
            .insert(GeneCallWarning::MissingRequiredPosition(missing));
        return finalize_gene_call_result(result, definition, exemption);
    }

    if orig_data.effectively_phased {
        let diplotype_matches = orig_data.compute_diplotypes(false);
        if !diplotype_matches.is_empty() {
            return finalize_gene_call_result(
                GeneCallResult::diplotypes(gene, orig_data, diplotype_matches, BTreeSet::new()),
                definition,
                exemption,
            );
        }
    }

    if orig_data.effectively_phased || orig_data.is_using_phase_sets() {
        let combo_data = initialize_match_data_with_exemption(
            &sample_id, gene, definition, exemption, allele_map, false,
        )?;
        let diplotype_matches = combo_data.compute_combination_diplotypes(definition, true);
        if !diplotype_matches.is_empty() {
            return finalize_gene_call_result(
                GeneCallResult::diplotypes(gene, combo_data, diplotype_matches, BTreeSet::new()),
                definition,
                exemption,
            );
        }
    }

    let combo_data = initialize_match_data_with_exemption(
        &sample_id, gene, definition, exemption, allele_map, false,
    )?;
    let combo_dip_matches = combo_data.compute_combination_diplotypes(definition, false);
    if !combo_dip_matches.is_empty()
        && (orig_data.effectively_phased || combo_data.is_using_phase_sets())
    {
        return finalize_gene_call_result(
            GeneCallResult::diplotypes(gene, combo_data, combo_dip_matches, BTreeSet::new()),
            definition,
            exemption,
        );
    }

    let haplotype_matches =
        call_haplotypes_for_lowest_function_gene(&combo_data, &combo_dip_matches, None)?;
    finalize_gene_call_result(
        GeneCallResult::haplotypes(gene, orig_data, haplotype_matches, BTreeSet::new()),
        definition,
        exemption,
    )
}

/// Calls the Java lowest-function fallback haplotype list after diplotype matching fails.
///
/// This mirrors `NamedAlleleMatcher.callHaplotypesForLowestFunctionGene`: preserve homozygous
/// components discovered through combination diplotypes, strip Reference when more than two
/// candidate haplotypes are present, then append DPYD HapB3 haplotype calls when applicable.
pub fn call_haplotypes_for_lowest_function_gene(
    combo_data: &MatchData,
    matches: &[DiplotypeMatch],
    mut dpyd_hap_b3_matcher: Option<&mut DpydHapB3Matcher<'_>>,
) -> Result<Vec<HaplotypeMatch>, MatchError> {
    if let Some(matcher) = dpyd_hap_b3_matcher.as_deref_mut() {
        matcher.call_hap_b3_haplotype_matches()?;
    }

    let mut homozygous = BTreeSet::new();
    for diplotype in matches {
        let mut haps = BTreeMap::<String, usize>::new();
        add_homozygous_candidates(&diplotype.haplotype1, combo_data, &mut haps);
        if let Some(haplotype2) = &diplotype.haplotype2 {
            add_homozygous_candidates(haplotype2, combo_data, &mut haps);
        }

        for (name, count) in haps {
            if count > 1 {
                homozygous.insert(name);
            }
        }
    }

    if dpyd_hap_b3_matcher
        .as_deref()
        .is_some_and(DpydHapB3Matcher::is_hap_b3_present)
    {
        homozygous.remove("Reference");
    }

    let mut hap_matches = combo_data.compare_permutations();
    let mut num_matches = hap_matches.len();
    if let Some(matcher) = dpyd_hap_b3_matcher.as_deref() {
        num_matches += matcher.num_hap_b3_called();
    }
    if num_matches > 2 {
        hap_matches.retain(|haplotype_match| haplotype_match.name != "Reference");
    }

    let mut final_haps = Vec::new();
    for haplotype_match in hap_matches {
        final_haps.push(haplotype_match.clone());
        if homozygous.remove(&haplotype_match.name) {
            final_haps.push(haplotype_match);
        }
    }

    if let Some(matcher) = dpyd_hap_b3_matcher {
        final_haps.extend(matcher.build_hap_b3_haplotype_matches()?);
    }

    if !homozygous.is_empty() {
        return Err(MatchError::MissingLowestFunctionHaplotype(
            homozygous.into_iter().collect(),
        ));
    }

    Ok(final_haps)
}

fn initialize_match_data(
    sample_id: &str,
    gene: &str,
    definition: &DefinitionFile,
    allele_map: &BTreeMap<String, SampleAlleleSummary>,
    assume_reference: bool,
) -> Result<MatchData, MatchError> {
    initialize_match_data_with_exemption(
        sample_id,
        gene,
        definition,
        None,
        allele_map,
        assume_reference,
    )
}

fn initialize_match_data_with_exemption(
    sample_id: &str,
    gene: &str,
    definition: &DefinitionFile,
    exemption: Option<&DefinitionExemption>,
    allele_map: &BTreeMap<String, SampleAlleleSummary>,
    assume_reference: bool,
) -> Result<MatchData, MatchError> {
    let mut data =
        MatchData::new_with_exemption(sample_id, gene, definition, exemption, allele_map);
    if data.sample_alleles.is_empty() {
        return Ok(data);
    }
    data.marshall_haplotypes(definition);
    if assume_reference {
        data.default_missing_alleles_to_reference()?;
    }
    data.generate_sample_permutations()?;
    Ok(data)
}

fn definition_without_dpyd_hap_b3(definition: &DefinitionFile) -> DefinitionFile {
    let mut definition = definition.clone();
    definition.named_alleles.retain(|allele| {
        allele.name != DpydHapB3Matcher::HAP_B3_ALLELE
            && allele.name != DpydHapB3Matcher::HAP_B3_INTRONIC_ALLELE
    });
    definition
}

/// Java DPYD HapB3 warning categories.
#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub enum DpydHapB3Warning {
    /// The intronic and exonic HapB3 calls disagree.
    IntronicMismatchExonic,
    /// HapB3 was assigned from the exonic SNP because the intronic SNP is missing.
    ExonicOnly,
    /// Only one allele was available for a HapB3 variant.
    AlleleCount {
        /// Number of observed alleles.
        count: usize,
        /// Variant rsid.
        rsid: Option<String>,
    },
}

/// DPYD-specific HapB3 matching state.
#[derive(Clone, Debug)]
pub struct DpydHapB3Matcher<'a> {
    definition: &'a DefinitionFile,
    orig_data: &'a MatchData,
    allele_map: &'a BTreeMap<String, SampleAlleleSummary>,
    hap_b3_exon_locus: VariantLocus,
    hap_b3_intron_locus: VariantLocus,
    has_hap_b3_variants: bool,
    has_non_hap_b3_variants: bool,
    hap_b3_intron_call: Option<Vec<String>>,
    hap_b3_call: Option<Vec<String>>,
    num_hap_b3_called: usize,
    warnings: BTreeSet<DpydHapB3Warning>,
}

impl<'a> DpydHapB3Matcher<'a> {
    /// Full HapB3 named allele.
    pub const HAP_B3_ALLELE: &'static str = "c.1129-5923C>G, c.1236G>A (HapB3)";
    /// Intronic HapB3 named allele.
    pub const HAP_B3_INTRONIC_ALLELE: &'static str = "c.1129-5923C>G";
    /// HapB3 exonic rsID.
    pub const HAP_B3_EXONIC_RSID: &'static str = "rs56038477";
    /// HapB3 intronic rsID.
    pub const HAP_B3_INTRONIC_RSID: &'static str = "rs75017182";
    /// HapB3 exonic partial name emitted by combination matching.
    pub const HAP_B3_EXONIC_PARTIAL: &'static str = "g.97573863C>T";
    /// HapB3 intronic partial name emitted by combination matching.
    pub const HAP_B3_INTRONIC_PARTIAL: &'static str = "g.97579893G>C";
    /// HapB3 exonic VCF position.
    pub const HAP_B3_EXONIC_POSITION: u64 = 97573863;
    /// HapB3 intronic VCF position.
    pub const HAP_B3_INTRONIC_POSITION: u64 = 97579893;

    /// Creates a DPYD HapB3 matcher from full DPYD match data.
    pub fn new(
        definition: &'a DefinitionFile,
        allele_map: &'a BTreeMap<String, SampleAlleleSummary>,
        orig_data: &'a MatchData,
    ) -> Result<Self, MatchError> {
        let hap_b3_exon_locus = definition
            .variants
            .iter()
            .find(|variant| variant.rsid.as_deref() == Some(Self::HAP_B3_EXONIC_RSID))
            .cloned()
            .ok_or_else(|| {
                MatchError::MissingDpydHapB3Variant(Self::HAP_B3_EXONIC_RSID.to_owned())
            })?;
        let hap_b3_intron_locus = definition
            .variants
            .iter()
            .find(|variant| variant.rsid.as_deref() == Some(Self::HAP_B3_INTRONIC_RSID))
            .cloned()
            .ok_or_else(|| {
                MatchError::MissingDpydHapB3Variant(Self::HAP_B3_INTRONIC_RSID.to_owned())
            })?;

        let has_hap_b3_variants = [&hap_b3_exon_locus, &hap_b3_intron_locus]
            .into_iter()
            .filter_map(|variant| {
                allele_map
                    .get(&variant.vcf_chr_position())
                    .map(|call| (variant, call))
            })
            .any(|(variant, call)| sample_call_has_non_reference(call, &variant.reference));

        let hap_b3_positions = [
            hap_b3_exon_locus.vcf_chr_position(),
            hap_b3_intron_locus.vcf_chr_position(),
        ];
        let has_non_hap_b3_variants = allele_map.iter().any(|(key, call)| {
            !hap_b3_positions.contains(key)
                && definition
                    .variant_for_position(call.position as u64)
                    .is_some_and(|variant| sample_call_has_non_reference(call, &variant.reference))
        });

        Ok(Self {
            definition,
            orig_data,
            allele_map,
            hap_b3_exon_locus,
            hap_b3_intron_locus,
            has_hap_b3_variants,
            has_non_hap_b3_variants,
            hap_b3_intron_call: None,
            hap_b3_call: None,
            num_hap_b3_called: 0,
            warnings: BTreeSet::new(),
        })
    }

    /// Gets whether HapB3 variants are present in the sample.
    pub fn has_hap_b3_variants(&self) -> bool {
        self.has_hap_b3_variants
    }

    /// Gets whether non-HapB3 variants are present in the sample.
    pub fn has_non_hap_b3_variants(&self) -> bool {
        self.has_non_hap_b3_variants
    }

    /// Gets the phased/inferred intronic HapB3 calls.
    pub fn hap_b3_intron_call(&self) -> Option<&[String]> {
        self.hap_b3_intron_call.as_deref()
    }

    /// Gets the phased/inferred HapB3 calls.
    pub fn hap_b3_call(&self) -> Option<&[String]> {
        self.hap_b3_call.as_deref()
    }

    /// Gets Java-equivalent warning categories emitted by this primitive.
    pub fn warnings(&self) -> &BTreeSet<DpydHapB3Warning> {
        &self.warnings
    }

    /// Gets how many HapB3 or HapB3-intronic strands were called.
    pub fn num_hap_b3_called(&self) -> usize {
        self.num_hap_b3_called
    }

    /// Gets whether HapB3 or HapB3-intronic named allele is called.
    pub fn is_hap_b3_present(&self) -> bool {
        self.num_hap_b3_called > 0
    }

    /// Calls HapB3 strand state for the non-diplotype path.
    pub fn call_hap_b3_haplotype_matches(&mut self) -> Result<(), MatchError> {
        let mut is_effectively_phased = self.orig_data.effectively_phased;
        if !is_effectively_phased && self.orig_data.is_using_phase_sets() {
            let exon_ps = self.orig_data.phase_set(self.hap_b3_exon_locus.position);
            let intron_ps = self.orig_data.phase_set(self.hap_b3_intron_locus.position);
            if exon_ps.is_some() && exon_ps == intron_ps {
                is_effectively_phased = true;
            }
        }

        let hap_b3_allele = self.find_hap_b3_allele(Self::HAP_B3_ALLELE)?.clone();
        let hap_b3_exon_sample = self
            .allele_map
            .get(&self.hap_b3_exon_locus.vcf_chr_position());
        let hap_b3_intron_sample = self
            .allele_map
            .get(&self.hap_b3_intron_locus.vcf_chr_position());
        let hap_b3_exon_locus = self.hap_b3_exon_locus.clone();
        let hap_b3_intron_locus = self.hap_b3_intron_locus.clone();

        if hap_b3_exon_sample.is_some() || hap_b3_intron_sample.is_some() {
            if let Some(intron_sample) = hap_b3_intron_sample {
                let intronic_hap_b3 =
                    self.call_hap_b3(&hap_b3_allele, &hap_b3_intron_locus, intron_sample)?;
                let num_intronic_strands = count_called_hap_b3_strands(&intronic_hap_b3);

                if let Some(exon_sample) = hap_b3_exon_sample {
                    let exonic_hap_b3 =
                        self.call_hap_b3(&hap_b3_allele, &hap_b3_exon_locus, exon_sample)?;
                    let num_exonic_strands = count_called_hap_b3_strands(&exonic_hap_b3);

                    if num_intronic_strands == 0 {
                        self.hap_b3_call = Some(exonic_hap_b3);
                        if num_exonic_strands != 0 {
                            self.warnings
                                .insert(DpydHapB3Warning::IntronicMismatchExonic);
                        }
                    } else if num_intronic_strands == 1 {
                        if num_exonic_strands == 0 {
                            self.hap_b3_intron_call = Some(intronic_hap_b3);
                        } else if is_effectively_phased {
                            self.handle_phased_call(&intronic_hap_b3, &exonic_hap_b3, 0);
                            self.handle_phased_call(&intronic_hap_b3, &exonic_hap_b3, 1);
                        }
                    } else if num_intronic_strands == 2 {
                        if num_exonic_strands == 0 {
                            self.hap_b3_intron_call = Some(intronic_hap_b3);
                        } else if is_effectively_phased {
                            self.handle_phased_call(&intronic_hap_b3, &exonic_hap_b3, 0);
                            self.handle_phased_call(&intronic_hap_b3, &exonic_hap_b3, 1);
                        } else if num_exonic_strands == 2 {
                            self.handle_unphased_double_intronic_double_exonic(
                                &intronic_hap_b3,
                                &exonic_hap_b3,
                            );
                        }
                    }
                } else {
                    self.hap_b3_intron_call = Some(intronic_hap_b3);
                }
            } else if let Some(exon_sample) = hap_b3_exon_sample {
                self.hap_b3_call =
                    Some(self.call_hap_b3(&hap_b3_allele, &hap_b3_exon_locus, exon_sample)?);
                if self
                    .hap_b3_call
                    .as_ref()
                    .is_some_and(|call| count_called_hap_b3_strands(call) == 2)
                {
                    self.warnings.insert(DpydHapB3Warning::ExonicOnly);
                }
            }
        }

        self.num_hap_b3_called = self
            .hap_b3_intron_call
            .iter()
            .chain(self.hap_b3_call.iter())
            .flat_map(|call| call.iter())
            .filter(|call| call.as_str() == "1")
            .count();

        Ok(())
    }

    /// Builds HapB3 haplotype matches for the unphased lowest-function fallback path.
    pub fn build_hap_b3_haplotype_matches(&self) -> Result<Vec<HaplotypeMatch>, MatchError> {
        let mut matches = Vec::new();
        if let Some(hap_b3_call) = &self.hap_b3_call {
            let hap_b3_allele = self.find_hap_b3_allele(Self::HAP_B3_ALLELE)?;
            matches.extend(
                hap_b3_call
                    .iter()
                    .filter(|call| call.as_str() == "1")
                    .map(|_| HaplotypeMatch::from_named_allele(hap_b3_allele.clone())),
            );
        }
        if let Some(hap_b3_intron_call) = &self.hap_b3_intron_call {
            let hap_b3_intron_allele = self.find_hap_b3_allele(Self::HAP_B3_INTRONIC_ALLELE)?;
            matches.extend(
                hap_b3_intron_call
                    .iter()
                    .filter(|call| call.as_str() == "1")
                    .map(|_| HaplotypeMatch::from_named_allele(hap_b3_intron_allele.clone())),
            );
        }
        Ok(matches)
    }

    /// Adds the HapB3 call to phased sequence data that is reference everywhere else.
    pub fn add_phased_hap_b3_call_to_ref(
        &mut self,
        match_data: &MatchData,
    ) -> Result<Vec<DiplotypeMatch>, MatchError> {
        if !self.has_hap_b3_variants || self.has_non_hap_b3_variants {
            return Err(MatchError::InvalidDpydHapB3State(
                "add_phased_hap_b3_call_to_ref requires HapB3-only sample data".to_owned(),
            ));
        }

        let haps = match_data
            .permutations()
            .iter()
            .map(|sequence| self.call_phased_hap_b3(sequence, match_data, None))
            .collect::<Result<Vec<_>, _>>()?;

        let diplotype = if haps.len() == 2 {
            DiplotypeMatch::new(
                haps[0].clone(),
                Some(haps[1].clone()),
                match_data.permutations().iter().cloned().collect(),
            )
        } else {
            let sequence = match_data
                .permutations()
                .iter()
                .next()
                .cloned()
                .ok_or(MatchError::NoSampleAlleles)?;
            DiplotypeMatch::new(
                haps[0].clone(),
                Some(haps[0].clone()),
                vec![sequence.clone(), sequence],
            )
        };

        let mut diplotypes = vec![diplotype];
        diplotypes.sort();
        diplotypes.dedup();
        Ok(diplotypes)
    }

    /// Merges HapB3 strand calls into existing phased diplotype matches.
    pub fn merge_phased_hap_b3_call(
        &mut self,
        match_data: &MatchData,
        diplotype_matches: &[DiplotypeMatch],
    ) -> Result<Vec<DiplotypeMatch>, MatchError> {
        if !self.has_hap_b3_variants {
            return Ok(diplotype_matches.to_vec());
        }

        let mut final_matches = Vec::new();
        for diplotype in diplotype_matches {
            let seq1 = diplotype
                .haplotype1
                .sequences
                .iter()
                .next()
                .ok_or(MatchError::NoSampleAlleles)?;
            let haplotype1 =
                self.call_phased_hap_b3(seq1, match_data, Some(&diplotype.haplotype1))?;
            let haplotype2 = diplotype
                .haplotype2
                .as_ref()
                .ok_or_else(|| {
                    MatchError::InvalidDpydHapB3State(
                        "HapB3 merge requires diploid matches".to_owned(),
                    )
                })
                .and_then(|haplotype2| {
                    let seq2 = haplotype2
                        .sequences
                        .iter()
                        .next()
                        .ok_or(MatchError::NoSampleAlleles)?;
                    self.call_phased_hap_b3(seq2, match_data, Some(haplotype2))
                })?;

            final_matches.push(DiplotypeMatch::new(
                haplotype1,
                Some(haplotype2),
                diplotype
                    .sequence_pairs
                    .first()
                    .cloned()
                    .unwrap_or_else(|| vec![seq1.clone()]),
            ));
        }

        final_matches.sort();
        final_matches.dedup();
        Ok(final_matches)
    }

    /// Removes HapB3 partial calls from combination matches before HapB3 is merged as a named allele.
    pub fn fix_partials(
        &self,
        match_data: &MatchData,
        diplotype_matches: &[DiplotypeMatch],
    ) -> Result<Vec<DiplotypeMatch>, MatchError> {
        let mut updated_matches = Vec::new();
        for diplotype in diplotype_matches {
            let (haplotype1, modified1) =
                self.remove_hap_b3_partials(match_data, &diplotype.haplotype1)?;
            let haplotype2 = diplotype.haplotype2.as_ref().ok_or_else(|| {
                MatchError::InvalidDpydHapB3State(
                    "HapB3 partial fixing requires diploid matches".to_owned(),
                )
            })?;
            let (haplotype2, modified2) = self.remove_hap_b3_partials(match_data, haplotype2)?;

            if modified1 || modified2 {
                updated_matches.push(DiplotypeMatch::new(
                    haplotype1,
                    Some(haplotype2),
                    diplotype
                        .sequence_pairs
                        .first()
                        .cloned()
                        .unwrap_or_default(),
                ));
            } else {
                updated_matches.push(diplotype.clone());
            }
        }

        updated_matches.sort();
        updated_matches.dedup();
        Ok(updated_matches)
    }

    fn remove_hap_b3_partials(
        &self,
        match_data: &MatchData,
        haplotype_match: &HaplotypeMatch,
    ) -> Result<(HaplotypeMatch, bool), MatchError> {
        if !has_hap_b3_partial(&haplotype_match.name) {
            return Ok((haplotype_match.clone(), false));
        }

        let sequence = haplotype_match
            .sequences
            .iter()
            .next()
            .cloned()
            .ok_or(MatchError::NoSampleAlleles)?;
        let component_identifiers = component_names_for_match(haplotype_match);
        let components = component_identifiers
            .iter()
            .map(|identifier| {
                match_data
                    .haplotypes
                    .iter()
                    .find(|haplotype| haplotype.name == *identifier || haplotype.id == *identifier)
                    .cloned()
                    .ok_or_else(|| MatchError::MissingDpydComponentAllele(identifier.clone()))
            })
            .collect::<Result<Vec<_>, _>>()?;

        let remaining_partials = partials_for_match(haplotype_match)
            .into_iter()
            .filter(|(_, name)| {
                name != Self::HAP_B3_EXONIC_PARTIAL && name != Self::HAP_B3_INTRONIC_PARTIAL
            })
            .collect::<BTreeMap<_, _>>();

        if remaining_partials.is_empty()
            && components.len() == 1
            && components[0].name == "Reference"
        {
            return Ok((
                HaplotypeMatch::from_haplotype(
                    components[0].clone(),
                    match_data.positions.clone(),
                    sequence,
                ),
                true,
            ));
        }

        Ok((
            build_combination_match(
                &match_data.positions,
                &sequence,
                &components,
                remaining_partials,
            ),
            true,
        ))
    }

    fn call_phased_hap_b3(
        &mut self,
        sequence: &str,
        match_data: &MatchData,
        base_match: Option<&HaplotypeMatch>,
    ) -> Result<HaplotypeMatch, MatchError> {
        let allele_map = parse_sequence(sequence);
        let intron_allele = allele_map.get(&self.hap_b3_intron_locus.position);
        let mut has_intron_locus = false;
        let mut has_intron = false;
        if let Some(allele) = intron_allele
            && *allele != "."
        {
            has_intron_locus = true;
            has_intron = allele != &self.hap_b3_intron_locus.reference;
        }

        let has_exon = allele_map
            .get(&self.hap_b3_exon_locus.position)
            .is_some_and(|allele| *allele != "." && allele != &self.hap_b3_exon_locus.reference);

        let hap_b3_allele = if has_intron {
            if has_exon {
                self.find_hap_b3_allele(Self::HAP_B3_ALLELE)?
            } else {
                self.find_hap_b3_allele(Self::HAP_B3_INTRONIC_ALLELE)?
            }
        } else {
            if has_intron_locus {
                if has_exon {
                    self.warnings
                        .insert(DpydHapB3Warning::IntronicMismatchExonic);
                }
            } else if has_exon {
                self.warnings.insert(DpydHapB3Warning::ExonicOnly);
                return self.build_match(
                    match_data,
                    self.find_hap_b3_allele(Self::HAP_B3_ALLELE)?,
                    base_match,
                    sequence,
                );
            }
            self.find_hap_b3_allele("Reference")?
        };

        self.build_match(match_data, hap_b3_allele, base_match, sequence)
    }

    fn build_match(
        &self,
        match_data: &MatchData,
        hap_b3_allele: &NamedAllele,
        base_match: Option<&HaplotypeMatch>,
        sequence: &str,
    ) -> Result<HaplotypeMatch, MatchError> {
        let Some(base_match) = base_match else {
            return Ok(HaplotypeMatch::from_haplotype(
                hap_b3_allele.clone(),
                match_data.positions.clone(),
                sequence.to_owned(),
            ));
        };
        if hap_b3_allele.name == "Reference" {
            return Ok(base_match.clone());
        }
        if base_match.name == "Reference" {
            return Ok(HaplotypeMatch::from_haplotype(
                hap_b3_allele.clone(),
                match_data.positions.clone(),
                sequence.to_owned(),
            ));
        }

        let mut components = vec![hap_b3_allele.clone()];
        for identifier in component_names_for_match(base_match) {
            let component = match_data
                .haplotypes
                .iter()
                .find(|haplotype| haplotype.name == identifier || haplotype.id == identifier)
                .cloned()
                .ok_or_else(|| MatchError::MissingDpydComponentAllele(identifier.clone()))?;
            components.push(component);
        }

        Ok(build_combination_match(
            &match_data.positions,
            sequence,
            &components,
            BTreeMap::new(),
        ))
    }

    fn handle_phased_call(
        &mut self,
        intronic_calls: &[String],
        exonic_calls: &[String],
        index: usize,
    ) {
        let intronic_call = intronic_calls.get(index).map(String::as_str).unwrap_or(".");
        let exonic_call = exonic_calls.get(index).map(String::as_str).unwrap_or(".");

        let hap_b3_intron_call = self.hap_b3_intron_call.get_or_insert_with(Vec::new);
        let hap_b3_call = self.hap_b3_call.get_or_insert_with(Vec::new);

        if intronic_call == "." && exonic_call == "." {
            hap_b3_intron_call.push(".".to_owned());
            hap_b3_call.push(".".to_owned());
        }

        if intronic_call == "0" {
            hap_b3_intron_call.push("0".to_owned());
            hap_b3_call.push("0".to_owned());
            if exonic_call != "0" {
                self.warnings
                    .insert(DpydHapB3Warning::IntronicMismatchExonic);
            }
        } else if intronic_call == "1" {
            if exonic_call == "0" {
                hap_b3_intron_call.push("1".to_owned());
                hap_b3_call.push("0".to_owned());
            } else if exonic_call == "." {
                hap_b3_intron_call.push("1".to_owned());
                hap_b3_call.push(".".to_owned());
            } else {
                hap_b3_intron_call.push("0".to_owned());
                hap_b3_call.push("1".to_owned());
            }
        } else {
            hap_b3_intron_call.push(".".to_owned());
            hap_b3_call.push(".".to_owned());
        }
    }

    fn handle_unphased_double_intronic_double_exonic(
        &mut self,
        intronic_hap_b3: &[String],
        exonic_hap_b3: &[String],
    ) {
        self.hap_b3_intron_call = Some(Vec::new());
        self.hap_b3_call = Some(Vec::new());

        let num_intronic = intronic_hap_b3
            .iter()
            .filter(|call| call.as_str() == "1")
            .count();
        let num_exonic = exonic_hap_b3
            .iter()
            .filter(|call| call.as_str() == "1")
            .count();

        if num_intronic == num_exonic {
            if let Some(call) = &mut self.hap_b3_call {
                call.extend((0..num_intronic).map(|_| "1".to_owned()));
            }
        } else if num_intronic > num_exonic {
            if let Some(call) = &mut self.hap_b3_intron_call {
                call.push("1".to_owned());
                if num_exonic == 0 && num_intronic == 2 {
                    call.push("1".to_owned());
                }
            }
            if num_exonic == 1
                && let Some(call) = &mut self.hap_b3_call
            {
                call.push("1".to_owned());
            }
        } else if num_intronic == 1 {
            if let Some(call) = &mut self.hap_b3_call {
                call.push("1".to_owned());
            }
            self.warnings
                .insert(DpydHapB3Warning::IntronicMismatchExonic);
        } else {
            self.warnings
                .insert(DpydHapB3Warning::IntronicMismatchExonic);
        }
    }

    fn call_hap_b3(
        &mut self,
        hap_b3_allele: &NamedAllele,
        locus: &VariantLocus,
        sample_allele: &SampleAlleleSummary,
    ) -> Result<Vec<String>, MatchError> {
        let exon_allele = self.named_allele_vcf_allele(hap_b3_allele, locus)?;
        let mut calls = Vec::new();
        let mut count = 0;

        if let Some(allele1) = &sample_allele.allele1 {
            calls.push(if exon_allele == allele1 { "1" } else { "0" }.to_owned());
            count += 1;
        } else if sample_allele.phased {
            calls.push(".".to_owned());
        }

        if let Some(allele2) = &sample_allele.allele2 {
            calls.push(if exon_allele == allele2 { "1" } else { "0" }.to_owned());
            count += 1;
        }

        if count < 2 {
            self.warnings.insert(DpydHapB3Warning::AlleleCount {
                count,
                rsid: locus.rsid.clone(),
            });
        }

        Ok(calls)
    }

    fn find_hap_b3_allele(&self, name: &str) -> Result<&NamedAllele, MatchError> {
        self.definition
            .named_allele(name)
            .ok_or_else(|| MatchError::MissingDpydHapB3Allele(name.to_owned()))
    }

    fn named_allele_vcf_allele<'b>(
        &self,
        named_allele: &'b NamedAllele,
        locus: &VariantLocus,
    ) -> Result<&'b str, MatchError> {
        self.definition
            .index_for_position(locus.position)
            .and_then(|index| named_allele.alleles.get(index))
            .and_then(Option::as_deref)
            .ok_or_else(|| MatchError::MissingDpydHapB3Allele(named_allele.name.clone()))
    }
}

/// A named allele and the sample permutations that matched it.
#[derive(Clone, Debug, PartialEq)]
pub struct HaplotypeMatch {
    /// Named allele name.
    pub name: String,
    /// Matched named allele.
    pub haplotype: NamedAllele,
    /// Positions aligned to the matched haplotype allele arrays.
    pub positions: Vec<VariantLocus>,
    /// Matching sample permutation strings.
    pub sequences: BTreeSet<String>,
}

impl Eq for HaplotypeMatch {}

impl HaplotypeMatch {
    fn from_named_allele(haplotype: NamedAllele) -> Self {
        Self {
            name: haplotype.name.clone(),
            haplotype,
            positions: Vec::new(),
            sequences: BTreeSet::new(),
        }
    }

    fn from_haplotype(
        haplotype: NamedAllele,
        positions: Vec<VariantLocus>,
        sequence: String,
    ) -> Self {
        Self {
            name: haplotype.name.clone(),
            haplotype,
            positions,
            sequences: [sequence].into_iter().collect(),
        }
    }

    fn score(&self) -> i32 {
        named_allele_score(&self.haplotype)
    }
}

/// A matched diplotype.
#[derive(Clone, Debug)]
pub struct DiplotypeMatch {
    /// Rendered diplotype name.
    pub name: String,
    /// First haplotype match.
    pub haplotype1: HaplotypeMatch,
    /// Second haplotype match, absent for haploid calls.
    pub haplotype2: Option<HaplotypeMatch>,
    /// Java diplotype score.
    pub score: i32,
    /// Sequence pairs that support this diplotype.
    pub sequence_pairs: Vec<Vec<String>>,
}

impl PartialEq for DiplotypeMatch {
    fn eq(&self, other: &Self) -> bool {
        self.score == other.score
            && self.haplotype1 == other.haplotype1
            && self.haplotype2 == other.haplotype2
    }
}

impl Eq for DiplotypeMatch {}

impl DiplotypeMatch {
    fn new(
        haplotype1: HaplotypeMatch,
        haplotype2: Option<HaplotypeMatch>,
        sequence_pair: Vec<String>,
    ) -> Self {
        let (haplotype1, haplotype2) = if let Some(haplotype2) = haplotype2 {
            if compare_haplotype_matches(&haplotype2, &haplotype1).is_lt() {
                (haplotype2, Some(haplotype1))
            } else {
                (haplotype1, Some(haplotype2))
            }
        } else {
            (haplotype1, None)
        };

        let name = if let Some(haplotype2) = &haplotype2 {
            format!("{}/{}", haplotype1.name, haplotype2.name)
        } else {
            haplotype1.name.clone()
        };
        let score = diplotype_score(&haplotype1, haplotype2.as_ref(), &sequence_pair);

        Self {
            name,
            haplotype1,
            haplotype2,
            score,
            sequence_pairs: vec![sequence_pair],
        }
    }

    fn handle_suballele_conversion(&mut self, definition: &DefinitionFile) {
        for (suballele, core_allele) in &definition.suballeles_map {
            if !self.name.contains(suballele) {
                continue;
            }

            self.name = self.name.replace(suballele, core_allele);
            if self.haplotype1.name.contains(suballele) {
                self.haplotype1.name = self.haplotype1.name.replace(suballele, core_allele);
            }
            if let Some(haplotype2) = &mut self.haplotype2
                && haplotype2.name.contains(suballele)
            {
                haplotype2.name = haplotype2.name.replace(suballele, core_allele);
            }
        }
    }
}

impl Ord for DiplotypeMatch {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        other
            .score
            .cmp(&self.score)
            .then_with(|| compare_haplotype_matches(&self.haplotype1, &other.haplotype1))
            .then_with(|| match (&self.haplotype2, &other.haplotype2) {
                (None, None) => std::cmp::Ordering::Equal,
                (None, Some(_)) => std::cmp::Ordering::Less,
                (Some(_), None) => std::cmp::Ordering::Greater,
                (Some(left), Some(right)) => compare_haplotype_matches(left, right),
            })
    }
}

impl PartialOrd for DiplotypeMatch {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

/// Result-level gene call outcome matching Java `ResultBuilder` branches.
#[derive(Clone, Debug, PartialEq)]
pub enum GeneCallKind {
    /// No callable sample alleles were available, or a required-position gate failed.
    NoCall,
    /// One or more diplotype matches were called.
    Diplotypes(Vec<DiplotypeMatch>),
    /// Potential haplotype matches were returned because no diplotype call could be made.
    Haplotypes(Vec<HaplotypeMatch>),
}

/// Java result-builder warning categories ported so far.
#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub enum GeneCallWarning {
    /// Java `unphased-priority` note emitted when a priority diplotype is selected.
    UnphasedPriority,
    /// Java `missing-required-position` note emitted when a required variant is absent.
    MissingRequiredPosition(Vec<String>),
    /// Java `missing-amp1-position` note emitted when AMP Tier 1 variants are absent.
    MissingAmp1Position(Vec<String>),
}

/// Minimal Rust gene-call result used while porting `NamedAlleleMatcher`.
#[derive(Clone, Debug, PartialEq)]
pub struct GeneCallResult {
    /// Gene symbol.
    pub gene: String,
    /// Match data used for the Java-equivalent result branch.
    pub match_data: MatchData,
    /// Called result branch.
    pub kind: GeneCallKind,
    /// DPYD HapB3 warning categories emitted during matching.
    pub dpyd_hap_b3_warnings: BTreeSet<DpydHapB3Warning>,
    /// Generic gene-call warnings emitted by Java `ResultBuilder`.
    pub warnings: BTreeSet<GeneCallWarning>,
}

impl GeneCallResult {
    fn no_call(gene: impl Into<String>, match_data: MatchData) -> Self {
        Self {
            gene: gene.into(),
            match_data,
            kind: GeneCallKind::NoCall,
            dpyd_hap_b3_warnings: BTreeSet::new(),
            warnings: BTreeSet::new(),
        }
    }

    fn diplotypes(
        gene: impl Into<String>,
        match_data: MatchData,
        matches: Vec<DiplotypeMatch>,
        dpyd_hap_b3_warnings: BTreeSet<DpydHapB3Warning>,
    ) -> Self {
        Self {
            gene: gene.into(),
            match_data,
            kind: GeneCallKind::Diplotypes(matches),
            dpyd_hap_b3_warnings,
            warnings: BTreeSet::new(),
        }
    }

    fn haplotypes(
        gene: impl Into<String>,
        match_data: MatchData,
        matches: Vec<HaplotypeMatch>,
        dpyd_hap_b3_warnings: BTreeSet<DpydHapB3Warning>,
    ) -> Self {
        Self {
            gene: gene.into(),
            match_data,
            kind: GeneCallKind::Haplotypes(matches),
            dpyd_hap_b3_warnings,
            warnings: BTreeSet::new(),
        }
    }
}

/// One sample allele in Java `SampleAllele` shape.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SampleAllele {
    chromosome: String,
    position: usize,
    allele1: Option<String>,
    allele2: Option<String>,
    computed_allele1: Option<String>,
    computed_allele2: Option<String>,
    vcf_alleles: Vec<String>,
    vcf_call: String,
    undocumented_variations: BTreeSet<String>,
    treat_undocumented_variations_as_reference: bool,
    phased: bool,
    effectively_phased: bool,
    phase_set: Option<i32>,
}

impl SampleAllele {
    /// Builds a matcher sample allele from a VCF call summary.
    pub fn from_summary(summary: &SampleAlleleSummary) -> Self {
        let computed_allele1 = compute_sample_allele(
            summary.allele1.as_deref(),
            summary,
            summary.vcf_call.split(['|', '/']).next(),
        );
        let computed_allele2 = compute_sample_allele(
            summary.allele2.as_deref(),
            summary,
            summary.vcf_call.split(['|', '/']).nth(1),
        );
        Self {
            chromosome: summary.chromosome.clone(),
            position: summary.position,
            allele1: summary.allele1.clone(),
            allele2: summary.allele2.clone(),
            computed_allele1,
            computed_allele2,
            vcf_alleles: summary.vcf_alleles.clone(),
            vcf_call: summary.vcf_call.clone(),
            undocumented_variations: summary.undocumented_variations.clone(),
            treat_undocumented_variations_as_reference: summary
                .treat_undocumented_variations_as_reference,
            phased: summary.phased,
            effectively_phased: summary.effectively_phased,
            phase_set: summary.phase_set,
        }
    }

    /// Java `SampleAllele.getVcfCall`.
    pub fn vcf_call(&self) -> &str {
        &self.vcf_call
    }

    /// Java `SampleAllele.getVcfAlleles`.
    pub fn vcf_alleles(&self) -> &[String] {
        &self.vcf_alleles
    }

    /// Whether the original VCF genotype was phased.
    pub fn phased(&self) -> bool {
        self.phased
    }

    /// Java-retained phase set.
    pub fn phase_set(&self) -> Option<i32> {
        self.phase_set
    }

    /// Java `SampleAllele.getUndocumentedVariations`.
    pub fn undocumented_variations(&self) -> &BTreeSet<String> {
        &self.undocumented_variations
    }

    /// Java `SampleAllele.isTreatUndocumentedVariationsAsReference`.
    pub fn treat_undocumented_variations_as_reference(&self) -> bool {
        self.treat_undocumented_variations_as_reference
    }

    #[cfg(test)]
    fn new(
        chromosome: &str,
        position: usize,
        allele1: Option<&str>,
        allele2: Option<&str>,
        phased: bool,
        effectively_phased: bool,
        phase_set: Option<i32>,
    ) -> Self {
        Self {
            chromosome: chromosome.to_owned(),
            position,
            allele1: allele1.map(str::to_owned),
            allele2: allele2.map(str::to_owned),
            computed_allele1: allele1.map(str::to_owned).or_else(|| Some(".".to_owned())),
            computed_allele2: allele2.map(str::to_owned),
            vcf_alleles: allele1
                .into_iter()
                .chain(allele2)
                .map(str::to_owned)
                .collect(),
            vcf_call: format!(
                "{}{}{}",
                allele1.unwrap_or("."),
                if phased { "|" } else { "/" },
                allele2.unwrap_or(".")
            ),
            undocumented_variations: BTreeSet::new(),
            treat_undocumented_variations_as_reference: false,
            phased,
            effectively_phased,
            phase_set,
        }
    }

    fn append_allele(&self, first_allele: bool, out: &mut String) {
        out.push_str(&self.position.to_string());
        out.push(':');
        if first_allele {
            out.push_str(self.computed_allele1.as_deref().unwrap_or("."));
        } else {
            out.push_str(self.computed_allele2.as_deref().unwrap_or("."));
        }
        out.push(';');
    }

    fn is_homozygous(&self) -> bool {
        self.computed_allele1.is_some() && self.computed_allele1 == self.computed_allele2
    }
}

fn compute_sample_allele(
    allele: Option<&str>,
    summary: &SampleAlleleSummary,
    rendered_part: Option<&str>,
) -> Option<String> {
    let Some(allele) = allele else {
        return (rendered_part == Some(".")).then(|| ".".to_owned());
    };
    if summary.treat_undocumented_variations_as_reference
        && summary.undocumented_variations.contains(allele)
    {
        return summary.vcf_alleles.first().cloned();
    }
    Some(allele.to_owned())
}

impl Ord for SampleAllele {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        chromosome::compare_names(Some(&self.chromosome), Some(&other.chromosome))
            .then_with(|| self.position.cmp(&other.position))
    }
}

impl PartialOrd for SampleAllele {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

/// Matching error.
#[derive(Debug, Eq, PartialEq)]
pub enum MatchError {
    /// Cannot generate permutations with no sample alleles.
    NoSampleAlleles,
    /// The definition lacks a reference named allele.
    NoReference(String),
    /// The DPYD definition lacks a HapB3 variant.
    MissingDpydHapB3Variant(String),
    /// The DPYD definition lacks a HapB3 named allele.
    MissingDpydHapB3Allele(String),
    /// The DPYD definition lacks a component allele needed for HapB3 merging.
    MissingDpydComponentAllele(String),
    /// DPYD HapB3 merge was requested for an invalid state.
    InvalidDpydHapB3State(String),
    /// Combination matching found a homozygous component that haplotype matching did not.
    MissingLowestFunctionHaplotype(Vec<String>),
    /// Java unphased priority listed a diplotype that is not in the current result.
    MissingPriorityDiplotype(String),
}

impl std::fmt::Display for MatchError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NoSampleAlleles => write!(f, "No alleles to generate permutations for"),
            Self::NoReference(gene) => write!(f, "{gene} does not have a reference"),
            Self::MissingDpydHapB3Variant(rsid) => {
                write!(f, "DPYD definition is missing HapB3 variant {rsid}")
            }
            Self::MissingDpydHapB3Allele(name) => {
                write!(f, "DPYD definition is missing HapB3 allele ({name})")
            }
            Self::MissingDpydComponentAllele(name) => {
                write!(f, "Cannot find DPYD allele '{name}'")
            }
            Self::InvalidDpydHapB3State(message) => write!(f, "{message}"),
            Self::MissingLowestFunctionHaplotype(names) => write!(
                f,
                "Combination matching found {names:?} but haplotype matching didn't"
            ),
            Self::MissingPriorityDiplotype(name) => {
                write!(f, "Cannot find priority diplotype {name}")
            }
        }
    }
}

impl std::error::Error for MatchError {}

fn sample_call_has_non_reference(call: &SampleAlleleSummary, reference: &str) -> bool {
    call.allele1
        .iter()
        .chain(call.allele2.iter())
        .any(|allele| allele != "." && allele != reference)
}

fn count_called_hap_b3_strands(calls: &[String]) -> usize {
    calls.iter().filter(|call| call.as_str() != ".").count()
}

fn add_homozygous_candidates(
    haplotype_match: &HaplotypeMatch,
    match_data: &MatchData,
    haps: &mut BTreeMap<String, usize>,
) {
    for name in component_display_names_for_match(haplotype_match, match_data) {
        if !(haplotype_match.haplotype.is_combination_or_partial
            && haplotype_match.haplotype.num_partials > 0
            && name == "Reference")
        {
            *haps.entry(name).or_default() += 1;
        }
    }
}

fn component_display_names_for_match(
    haplotype_match: &HaplotypeMatch,
    match_data: &MatchData,
) -> Vec<String> {
    component_names_for_match(haplotype_match)
        .into_iter()
        .map(|identifier| {
            match_data
                .haplotypes
                .iter()
                .find(|haplotype| haplotype.name == identifier || haplotype.id == identifier)
                .map(|haplotype| haplotype.name.clone())
                .unwrap_or(identifier)
        })
        .collect()
}

fn component_names_for_match(haplotype_match: &HaplotypeMatch) -> Vec<String> {
    if haplotype_match.haplotype.is_combination_or_partial {
        let components = haplotype_match
            .haplotype
            .id
            .split(" + ")
            .map(str::to_owned)
            .filter(|name| !is_partial_name(name))
            .collect::<Vec<_>>();
        if !components.is_empty() {
            return components;
        }
    }

    if is_combination_name(&haplotype_match.name) {
        split_combination_name(&haplotype_match.name)
            .into_iter()
            .filter(|name| !is_partial_name(name))
            .collect()
    } else {
        vec![haplotype_match.name.clone()]
    }
}

fn has_hap_b3_partial(name: &str) -> bool {
    name.contains(DpydHapB3Matcher::HAP_B3_EXONIC_PARTIAL)
        || name.contains(DpydHapB3Matcher::HAP_B3_INTRONIC_PARTIAL)
}

fn partials_for_match(haplotype_match: &HaplotypeMatch) -> BTreeMap<u64, String> {
    partial_names_for_match(&haplotype_match.name)
        .filter_map(|name| {
            let position = if name.as_str() == DpydHapB3Matcher::HAP_B3_EXONIC_PARTIAL {
                DpydHapB3Matcher::HAP_B3_EXONIC_POSITION
            } else if name.as_str() == DpydHapB3Matcher::HAP_B3_INTRONIC_PARTIAL {
                DpydHapB3Matcher::HAP_B3_INTRONIC_POSITION
            } else {
                parse_partial_position(&name)?
            };
            Some((position, name))
        })
        .collect()
}

fn partial_names_for_match(name: &str) -> Box<dyn Iterator<Item = String> + '_> {
    if is_combination_name(name) {
        Box::new(
            split_combination_name(name)
                .into_iter()
                .filter(|name| is_partial_name(name)),
        )
    } else if is_partial_name(name) {
        Box::new(std::iter::once(name.to_owned()))
    } else {
        Box::new(std::iter::empty())
    }
}

fn parse_partial_position(name: &str) -> Option<u64> {
    let rest = name.strip_prefix("g.")?;
    let digits = rest
        .chars()
        .take_while(char::is_ascii_digit)
        .collect::<String>();
    digits.parse().ok()
}

fn generate_permutations(sample_alleles: &[SampleAllele]) -> Result<BTreeSet<String>, MatchError> {
    let is_s1_blank = sample_alleles
        .iter()
        .all(|sample_allele| sample_allele.allele1.is_none());
    let is_s2_blank = sample_alleles
        .iter()
        .all(|sample_allele| sample_allele.allele2.is_none());
    let has_phase_sets = sample_alleles
        .iter()
        .any(|sample_allele| sample_allele.phase_set.is_some());

    let mut permutations = BTreeSet::new();
    if !is_s1_blank {
        let phase_sets = if has_phase_sets {
            Some(BTreeMap::new())
        } else {
            None
        };
        permutations.extend(generate_permutation_branch(
            sample_alleles,
            0,
            is_s2_blank,
            true,
            String::new(),
            phase_sets,
        ));
    }

    if !is_s2_blank {
        let phase_sets = if has_phase_sets {
            Some(BTreeMap::new())
        } else {
            None
        };
        permutations.extend(generate_permutation_branch(
            sample_alleles,
            0,
            is_s1_blank,
            false,
            String::new(),
            phase_sets,
        ));
    }

    if permutations.is_empty() {
        Err(MatchError::NoSampleAlleles)
    } else {
        Ok(permutations)
    }
}

fn generate_permutation_branch(
    sample_alleles: &[SampleAllele],
    position: usize,
    haploid: bool,
    first_allele: bool,
    allele_so_far: String,
    phase_sets: Option<BTreeMap<i32, bool>>,
) -> BTreeSet<String> {
    if position >= sample_alleles.len() {
        return [allele_so_far].into_iter().collect();
    }

    let allele = &sample_alleles[position];
    let mut out = BTreeSet::new();

    if allele.effectively_phased || haploid {
        let mut next = allele_so_far;
        allele.append_allele(first_allele, &mut next);
        out.extend(generate_permutation_branch(
            sample_alleles,
            position + 1,
            haploid,
            first_allele,
            next,
            phase_sets,
        ));
    } else if let (Some(phase_set), Some(phase_sets)) = (allele.phase_set, phase_sets.as_ref()) {
        if let Some(in_phase) = phase_sets.get(&phase_set).copied() {
            let mut next = allele_so_far;
            allele.append_allele(in_phase, &mut next);
            out.extend(generate_permutation_branch(
                sample_alleles,
                position + 1,
                false,
                first_allele,
                next,
                Some(phase_sets.clone()),
            ));
        } else {
            for in_phase in [true, false] {
                let mut next_phase_sets = phase_sets.clone();
                next_phase_sets.insert(phase_set, in_phase);
                let mut next = allele_so_far.clone();
                allele.append_allele(in_phase, &mut next);
                out.extend(generate_permutation_branch(
                    sample_alleles,
                    position + 1,
                    false,
                    first_allele,
                    next,
                    Some(next_phase_sets),
                ));
            }
        }
    } else {
        for first in [true, false] {
            let mut next = allele_so_far.clone();
            allele.append_allele(first, &mut next);
            out.extend(generate_permutation_branch(
                sample_alleles,
                position + 1,
                false,
                first_allele,
                next,
                phase_sets.clone(),
            ));
        }
    }

    out
}

fn haplotype_matches_sequence(
    haplotype: &NamedAllele,
    positions: &[VariantLocus],
    sequence: &str,
) -> bool {
    let alleles_by_position = parse_sequence(sequence);

    for (index, variant) in positions.iter().enumerate() {
        let expected = haplotype.alleles.get(index).and_then(Option::as_deref);
        let Some(expected) = expected else {
            continue;
        };
        let Some(actual) = alleles_by_position.get(&variant.position) else {
            return false;
        };
        if !allele_matches(expected, actual) {
            return false;
        }
    }

    true
}

fn parse_sequence(sequence: &str) -> BTreeMap<u64, &str> {
    sequence
        .split(';')
        .filter(|part| !part.is_empty())
        .filter_map(|part| {
            let (position, allele) = part.split_once(':')?;
            Some((position.parse().ok()?, allele))
        })
        .collect()
}

fn allele_matches(expected: &str, actual: &str) -> bool {
    expected == actual || iupac_bases(expected).is_some_and(|bases| bases.contains(&actual))
}

fn sample_has_named_allele(
    allele_map: &BTreeMap<u64, &str>,
    named_allele: &NamedAllele,
    positions: &[VariantLocus],
) -> bool {
    for (index, allele) in named_allele.alleles.iter().enumerate() {
        let Some(expected) = allele.as_deref() else {
            continue;
        };
        let Some(position) = positions.get(index) else {
            return false;
        };
        let Some(sample_allele) = allele_map.get(&position.position) else {
            return false;
        };
        if !allele_matches(expected, sample_allele) {
            return false;
        }
    }

    true
}

fn calculate_partial_names(
    definition: &DefinitionFile,
    allele_map: &BTreeMap<u64, &str>,
    var_positions: &BTreeSet<u64>,
    covered_haps: &[NamedAllele],
    find_partials: bool,
) -> BTreeMap<u64, String> {
    if !find_partials {
        return BTreeMap::new();
    }

    let mut partial_positions = var_positions.clone();
    for haplotype in covered_haps {
        for position in &haplotype.core_positions {
            partial_positions.remove(position);
        }
    }

    partial_positions
        .into_iter()
        .filter_map(|position| {
            let variant = definition.variant_for_position(position)?;
            let allele = allele_map.get(&position)?;
            Some((position, variant.hgvs_for_vcf_allele(allele)))
        })
        .collect()
}

fn compute_viable_combinations(mut covered_haps: Vec<NamedAllele>) -> Vec<Vec<NamedAllele>> {
    covered_haps.sort_by(|left, right| {
        right
            .core_positions
            .len()
            .cmp(&left.core_positions.len())
            .then_with(|| compare_haplotype_names(&left.name, &right.name))
    });

    let mut combos: Vec<Vec<NamedAllele>> = Vec::new();
    for allele in covered_haps {
        let mut added = false;
        let mut overlap_count = 0;
        for combo in &mut combos {
            let mut overlaps_combo = false;
            for existing in combo.iter() {
                if overlaps(&existing.core_positions, &allele.core_positions) {
                    overlaps_combo = true;
                    if existing.core_positions.len() != allele.core_positions.len() {
                        overlap_count += 1;
                        break;
                    }
                }
            }
            if !overlaps_combo {
                combo.push(allele.clone());
                combo.sort_by(|left, right| compare_haplotype_names(&left.name, &right.name));
                added = true;
            }
        }
        if !added && overlap_count == 0 {
            combos.push(vec![allele]);
        }
    }

    combos
}

fn overlaps(left: &BTreeSet<u64>, right: &BTreeSet<u64>) -> bool {
    right.iter().any(|position| left.contains(position))
}

fn build_combination_match(
    positions: &[VariantLocus],
    sequence: &str,
    components: &[NamedAllele],
    partials: BTreeMap<u64, String>,
) -> HaplotypeMatch {
    let name = build_combination_name(components, &partials);
    let haplotype = build_combination_haplotype(positions, &name, components, &partials);
    HaplotypeMatch::from_haplotype(haplotype, positions.to_vec(), sequence.to_owned())
}

fn build_combination_name(components: &[NamedAllele], partials: &BTreeMap<u64, String>) -> String {
    let mut parts = Vec::new();
    let mut count = 0;

    let mut components = components.to_vec();
    components.sort_by(|left, right| compare_haplotype_names(&left.name, &right.name));

    for component in &components {
        if component.reference {
            continue;
        }
        if is_combination_name(&component.name) {
            parts.extend(split_combination_name(&component.name));
            count += 2;
        } else {
            parts.push(component.name.clone());
            count += 1;
        }
    }

    for partial in partials.values() {
        parts.push(partial.clone());
        count += 1;
    }

    let name = parts.join(" + ");
    if count > 1 { format!("[{name}]") } else { name }
}

fn build_combination_haplotype(
    positions: &[VariantLocus],
    name: &str,
    components: &[NamedAllele],
    partials: &BTreeMap<u64, String>,
) -> NamedAllele {
    let mut components = components.to_vec();
    components.sort_by(|left, right| compare_haplotype_names(&left.name, &right.name));

    let mut alleles = vec![None; positions.len()];
    let mut cpic_alleles = vec![None; positions.len()];
    let mut missing_positions = BTreeSet::new();
    let id = components
        .iter()
        .map(|component| {
            missing_positions.extend(component.missing_positions.iter().cloned());
            component.id.clone()
        })
        .collect::<Vec<_>>()
        .join(" + ");

    for component in &components {
        for index in 0..positions.len() {
            merge_allele(
                &mut alleles[index],
                component.alleles.get(index).cloned().unwrap_or_default(),
            );
            merge_allele(
                &mut cpic_alleles[index],
                component
                    .cpic_alleles
                    .get(index)
                    .cloned()
                    .unwrap_or_default(),
            );
        }
    }

    let mut core_positions = BTreeSet::new();
    for (index, allele) in alleles.iter().enumerate() {
        if allele.is_some() {
            core_positions.insert(positions[index].position);
        }
    }

    NamedAllele {
        id,
        name: name.to_owned(),
        alleles,
        cpic_alleles,
        population_frequency: None,
        reference: false,
        is_combination_or_partial: !partials.is_empty() || components.len() > 1,
        structural_variant: false,
        score: None,
        core_positions,
        missing_positions,
        score_override: None,
        num_combinations: components.len() as i32,
        num_partials: partials.len() as i32,
    }
}

fn merge_allele(target: &mut Option<String>, candidate: Option<String>) {
    if target.is_none() {
        *target = candidate;
    } else if let (Some(existing), Some(candidate)) = (target.as_ref(), candidate)
        && existing != &candidate
        && !alleles_are_compatible(existing, &candidate)
    {
        panic!("conflicting combination alleles: {existing} vs {candidate}");
    }
}

fn alleles_are_compatible(existing: &str, candidate: &str) -> bool {
    iupac_bases(existing).is_some_and(|bases| bases.contains(&candidate))
        || iupac_bases(candidate).is_some_and(|bases| bases.contains(&existing))
}

fn is_combination_name(name: &str) -> bool {
    name.starts_with('[') && name.contains(" + ") && name.ends_with(']')
}

fn split_combination_name(name: &str) -> Vec<String> {
    name.trim_start_matches('[')
        .trim_end_matches(']')
        .split(" + ")
        .map(str::to_owned)
        .collect()
}

fn iupac_bases(allele: &str) -> Option<&'static [&'static str]> {
    match allele {
        "A" | "C" | "G" | "T" => None,
        "R" => Some(&["A", "G"]),
        "Y" => Some(&["C", "T"]),
        "S" => Some(&["G", "C"]),
        "W" => Some(&["A", "T"]),
        "K" => Some(&["G", "T"]),
        "M" => Some(&["A", "C"]),
        "B" => Some(&["C", "G", "T"]),
        "D" => Some(&["A", "G", "T"]),
        "H" => Some(&["A", "C", "T"]),
        "V" => Some(&["A", "C", "G"]),
        "N" => Some(&["A", "C", "G", "T"]),
        _ => None,
    }
}

fn compare_haplotype_matches(left: &HaplotypeMatch, right: &HaplotypeMatch) -> std::cmp::Ordering {
    compare_haplotype_names(&left.name, &right.name)
        .then_with(|| left.score().cmp(&right.score()))
        .then_with(|| compare_sequences(&left.sequences, &right.sequences))
        .then_with(|| left.haplotype.id.cmp(&right.haplotype.id))
}

pub(crate) fn compare_haplotype_names(left: &str, right: &str) -> std::cmp::Ordering {
    if left == right {
        return std::cmp::Ordering::Equal;
    }

    let left_combination = is_combination_name(left);
    let right_combination = is_combination_name(right);
    match (left_combination, right_combination) {
        (true, false) => return std::cmp::Ordering::Greater,
        (false, true) => return std::cmp::Ordering::Less,
        (true, true) => return compare_combination_names(left, right),
        (false, false) => {}
    }

    const TOP: &[&str] = &["Any", "All", "Reference"];
    const BOTTOM: &[&str] = &["Other", "Unknown"];

    if TOP.contains(&left) || BOTTOM.contains(&right) {
        return std::cmp::Ordering::Less;
    }
    if TOP.contains(&right) || BOTTOM.contains(&left) {
        return std::cmp::Ordering::Greater;
    }

    natural_compare(left, right)
}

fn compare_combination_names(left: &str, right: &str) -> std::cmp::Ordering {
    let left_parts = split_combination_name(left);
    let right_parts = split_combination_name(right);

    left_parts.len().cmp(&right_parts.len()).then_with(|| {
        left_parts
            .iter()
            .zip(&right_parts)
            .map(
                |(left, right)| match (is_partial_name(left), is_partial_name(right)) {
                    (true, false) => std::cmp::Ordering::Greater,
                    (false, true) => std::cmp::Ordering::Less,
                    _ => natural_compare(left, right),
                },
            )
            .find(|ordering| !ordering.is_eq())
            .unwrap_or(std::cmp::Ordering::Equal)
    })
}

fn is_partial_name(name: &str) -> bool {
    let Some(rest) = name.strip_prefix("g.") else {
        return false;
    };
    let mut chars = rest.chars().peekable();
    if chars.next_if(char::is_ascii_digit).is_none() {
        return false;
    }
    while chars.next_if(char::is_ascii_digit).is_some() {}

    let tail = chars.collect::<String>();
    !tail.is_empty()
        && tail
            .chars()
            .all(|ch| matches!(ch, '?' | '=' | '>' | 'A' | 'C' | 'G' | 'T') || ch.is_ascii_digit())
}

fn natural_compare(left: &str, right: &str) -> std::cmp::Ordering {
    let left = natural_parts(left);
    let right = natural_parts(right);

    for (left, right) in left.iter().zip(&right) {
        let ordering = match (left, right) {
            (NaturalPart::Number(left), NaturalPart::Number(right)) => left.cmp(right),
            (NaturalPart::Text(left), NaturalPart::Text(right)) => left.cmp(right),
            (NaturalPart::Text(_), NaturalPart::Number(_)) => std::cmp::Ordering::Greater,
            (NaturalPart::Number(_), NaturalPart::Text(_)) => std::cmp::Ordering::Less,
        };
        if !ordering.is_eq() {
            return ordering;
        }
    }

    left.len().cmp(&right.len())
}

#[derive(Debug)]
enum NaturalPart {
    Number(i32),
    Text(String),
}

fn natural_parts(s: &str) -> Vec<NaturalPart> {
    let mut parts = Vec::new();
    let mut start = 0;
    let mut was_digit = false;
    let chars = s.char_indices().collect::<Vec<_>>();

    for (idx, (pos, ch)) in chars.iter().enumerate() {
        let is_digit = ch.is_ascii_digit();
        if idx == 0 {
            was_digit = is_digit;
        }
        if is_digit != was_digit {
            push_natural_part(&mut parts, &s[start..*pos], was_digit);
            start = *pos;
            was_digit = is_digit;
        }
    }

    if start < s.len() {
        push_natural_part(&mut parts, &s[start..], was_digit);
    }

    parts
}

fn push_natural_part(parts: &mut Vec<NaturalPart>, part: &str, digit: bool) {
    if part.is_empty() {
        return;
    }
    if digit {
        parts.push(NaturalPart::Number(part.parse().unwrap_or(i32::MAX)));
    } else {
        parts.push(NaturalPart::Text(part.to_owned()));
    }
}

fn compare_sequences(left: &BTreeSet<String>, right: &BTreeSet<String>) -> std::cmp::Ordering {
    left.len()
        .cmp(&right.len())
        .then_with(|| left.iter().cmp(right.iter()))
}

fn named_allele_score(named_allele: &NamedAllele) -> i32 {
    if let Some(score) = named_allele.score_override {
        return score;
    }
    if let Some(score) = named_allele.score {
        return score;
    }
    named_allele
        .alleles
        .iter()
        .filter(|allele| allele.is_some())
        .count() as i32
        - named_allele.num_partials
}

fn diplotype_score(
    haplotype1: &HaplotypeMatch,
    haplotype2: Option<&HaplotypeMatch>,
    sequence_pair: &[String],
) -> i32 {
    let is_homozygous = haplotype2.is_some_and(|haplotype2| haplotype1.name == haplotype2.name);
    let score1_sequences = if is_homozygous {
        sequence_pair[..1].to_vec()
    } else {
        sequences_for_haplotype(haplotype1, sequence_pair)
    };
    let mut score = score_for_sample(haplotype1, &score1_sequences);

    if let Some(haplotype2) = haplotype2 {
        let score2_sequences = if is_homozygous && sequence_pair.len() > 1 {
            sequence_pair[1..].to_vec()
        } else {
            sequences_for_haplotype(haplotype2, sequence_pair)
        };
        score += score_for_sample(haplotype2, &score2_sequences);
    }

    score
}

fn sequences_for_haplotype(
    haplotype_match: &HaplotypeMatch,
    sequence_pair: &[String],
) -> Vec<String> {
    if haplotype_match.sequences.len() == 1 {
        return haplotype_match.sequences.iter().cloned().collect();
    }

    let sequences = sequence_pair
        .iter()
        .filter(|sequence| haplotype_match.sequences.contains(*sequence))
        .cloned()
        .collect::<Vec<_>>();
    if sequences.is_empty() {
        sequence_pair.to_vec()
    } else {
        sequences
    }
}

fn score_for_sample(haplotype_match: &HaplotypeMatch, sequences: &[String]) -> i32 {
    let mut score = haplotype_match.score();
    for (index, allele) in haplotype_match.haplotype.alleles.iter().enumerate() {
        let Some(allele) = allele.as_deref() else {
            continue;
        };
        if !is_iupac_wobble(allele) {
            continue;
        }

        let Some(position) = haplotype_match.positions.get(index) else {
            continue;
        };

        let all_reference = sequences.iter().all(|sequence| {
            parse_sequence(sequence)
                .get(&position.position)
                .is_some_and(|actual| *actual == position.reference)
        });
        if all_reference {
            score -= 1;
        }
    }
    score
}

fn is_iupac_wobble(allele: &str) -> bool {
    matches!(
        allele,
        "R" | "Y" | "S" | "W" | "K" | "M" | "B" | "D" | "H" | "V" | "N"
    )
}

fn perfect_pairs<T: Clone>(values: &[T]) -> Vec<(&T, T)> {
    let mut pairs = Vec::new();
    for x in 0..values.len() {
        for y in x..values.len() {
            pairs.push((&values[x], values[y].clone()));
        }
    }
    pairs
}

fn are_sample_alleles_haploid(sample_alleles: &[SampleAllele]) -> bool {
    let s1 = sample_alleles
        .iter()
        .filter(|sample_allele| sample_allele.allele1.is_none())
        .count();
    let s2 = sample_alleles
        .iter()
        .filter(|sample_allele| sample_allele.allele2.is_none())
        .count();

    s1 == sample_alleles.len() || s2 == sample_alleles.len()
}

#[cfg(test)]
mod tests {
    use std::collections::{BTreeMap, BTreeSet};

    use super::{
        DiplotypeMatch, DpydHapB3Matcher, DpydHapB3Warning, GeneCallKind, GeneCallResult,
        GeneCallWarning, HaplotypeMatch, MatchData, SampleAllele, build_combination_match,
        call_dpyd_lowest_function_gene, call_haplotypes_for_lowest_function_gene,
        call_ryr1_lowest_function_gene, call_standard_gene, call_standard_gene_with_exemption,
        compute_dpyd_lowest_function_diplotypes, definition_without_dpyd_hap_b3,
        finalize_gene_call_result, generate_permutations,
    };
    use crate::{
        definition::{
            DefinitionExemption, DefinitionFile, NamedAllele, UnphasedDiplotypePriority,
            VariantLocus, read_definition_file, read_exemptions_file,
        },
        vcf::{VcfRecords, VcfWarnings, allele_calls_for_locations, read_record_summaries},
    };

    #[test]
    fn generates_sample_permutations_like_java_combination_util() {
        let sample_alleles = vec![
            SampleAllele::new("chr1", 1, Some("T"), Some("T"), true, true, None),
            SampleAllele::new("chr1", 2, Some("A"), Some("T"), false, false, None),
            SampleAllele::new("chr1", 3, Some("C"), Some("C"), false, true, None),
            SampleAllele::new("chr1", 4, Some("C"), Some("G"), false, false, None),
        ];

        assert_eq!(
            generate_permutations(&sample_alleles).expect("permutations"),
            [
                "1:T;2:A;3:C;4:C;".to_owned(),
                "1:T;2:A;3:C;4:G;".to_owned(),
                "1:T;2:T;3:C;4:C;".to_owned(),
                "1:T;2:T;3:C;4:G;".to_owned()
            ]
            .into_iter()
            .collect()
        );
    }

    #[test]
    fn compares_permutations_to_named_alleles_like_java_match_data() {
        let definition = synthetic_definition();
        let allele_map = [
            sample_call("chr1", 1, Some("T"), Some("T"), true, true, None),
            sample_call("chr1", 2, Some("A"), Some("T"), false, false, None),
            sample_call("chr1", 3, Some("C"), Some("C"), false, true, None),
            sample_call("chr1", 4, Some("C"), Some("G"), false, false, None),
        ]
        .into_iter()
        .map(|call| (format!("{}:{}", call.chromosome, call.position), call))
        .collect::<BTreeMap<_, _>>();

        let mut data = MatchData::new("Sample_1", "GENE", &definition, &allele_map);
        data.marshall_haplotypes(&definition);
        data.generate_sample_permutations().expect("permutations");

        let matches = data.compare_permutations();
        assert_eq!(
            matches
                .iter()
                .map(|haplotype_match| haplotype_match.name.as_str())
                .collect::<Vec<_>>(),
            ["*1", "*2"]
        );
    }

    #[test]
    fn builds_match_data_for_java_haplotyper_fixture() {
        let definition = read_definition_file(std::path::Path::new(
            "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/haplotyper.json",
        ))
        .expect("definition");
        let records = read_record_summaries(
            "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/haplotyper.vcf",
            Some("NA12878"),
        )
        .expect("vcf");
        let allele_map = records
            .records
            .iter()
            .filter_map(|record| {
                record
                    .allele_call
                    .clone()
                    .map(|call| (format!("{}:{}", call.chromosome, call.position), call))
            })
            .collect::<BTreeMap<_, _>>();

        let mut data = MatchData::new("NA12878", "CYP3A5", &definition, &allele_map);
        assert_eq!(data.positions.len(), 3);
        assert_eq!(data.missing_positions.len(), 0);

        data.marshall_haplotypes(&definition);
        assert_eq!(data.haplotypes.len(), 2);

        data.generate_sample_permutations().expect("permutations");
        assert_eq!(
            data.permutations(),
            &[
                "1:C;2:C;4:TG;".to_owned(),
                "1:C;2:CT;4:TG;".to_owned(),
                "1:T;2:C;4:T;".to_owned(),
                "1:T;2:CT;4:T;".to_owned()
            ]
            .into_iter()
            .collect()
        );

        assert_eq!(
            data.compare_permutations()
                .iter()
                .map(|haplotype_match| haplotype_match.name.as_str())
                .collect::<Vec<_>>(),
            ["*1", "*2"]
        );
        assert_eq!(
            data.compute_diplotypes(false)
                .iter()
                .map(|diplotype| diplotype.name.as_str())
                .collect::<Vec<_>>(),
            ["*1/*2"]
        );
    }

    #[test]
    fn computes_diplotype_pairs_like_java_diplotype_matcher_tests() {
        let definition = synthetic_diplotype_definition();

        assert_eq!(
            compute_diplotype_names(
                &definition,
                [
                    sample_call("chr1", 1, Some("A"), Some("G"), false, false, None),
                    sample_call("chr1", 2, Some("C"), Some("C"), false, true, None),
                    sample_call("chr1", 3, Some("C"), Some("C"), false, true, None),
                ],
                false,
            ),
            ["*1/*4a"]
        );

        assert_eq!(
            compute_diplotype_names(
                &definition,
                [
                    sample_call("chr1", 1, Some("A"), Some("G"), false, false, None),
                    sample_call("chr1", 2, Some("C"), Some("T"), false, false, None),
                    sample_call("chr1", 3, Some("C"), Some("T"), false, false, None),
                ],
                false,
            ),
            ["*1/*4b", "*1/*17", "*1/*4a", "*4a/*17"]
        );

        assert_eq!(
            compute_diplotype_names(
                &definition,
                [
                    sample_call("chr1", 1, Some("A"), Some("A"), false, true, None),
                    sample_call("chr1", 2, Some("C"), Some("T"), false, false, None),
                    sample_call("chr1", 3, Some("C"), Some("T"), false, false, None),
                ],
                false,
            ),
            ["*1/*17"]
        );

        assert_eq!(
            compute_diplotype_names(
                &definition,
                [
                    sample_call("chr1", 1, Some("G"), Some("G"), false, true, None),
                    sample_call("chr1", 2, Some("T"), Some("T"), false, true, None),
                    sample_call("chr1", 3, Some("C"), Some("T"), false, false, None),
                ],
                false,
            ),
            ["*4a/*4b", "*4a/*17", "*4a/*4a"]
        );

        assert_eq!(
            compute_diplotype_names(
                &definition,
                [
                    sample_call("chr1", 1, Some("G"), Some("G"), false, true, None),
                    sample_call("chr1", 2, Some("C"), Some("T"), false, false, None),
                    sample_call("chr1", 3, Some("C"), Some("T"), false, false, None),
                ],
                false,
            ),
            ["*4a/*4b", "*4a/*17", "*4a/*4a"]
        );
    }

    #[test]
    fn filters_to_top_candidate_diplotype_scores_like_java() {
        let definition = synthetic_diplotype_definition();

        assert_eq!(
            compute_diplotype_names(
                &definition,
                [
                    sample_call("chr1", 1, Some("A"), Some("G"), false, false, None),
                    sample_call("chr1", 2, Some("C"), Some("T"), false, false, None),
                    sample_call("chr1", 3, Some("C"), Some("T"), false, false, None),
                ],
                true,
            ),
            ["*1/*4b"]
        );
    }

    #[test]
    fn missing_position_marshalling_tracks_missing_haplotype_positions_and_keeps_score() {
        let definition = synthetic_diplotype_definition();
        let allele_map = [
            sample_call("chr1", 1, Some("A"), Some("G"), false, false, None),
            sample_call("chr1", 2, Some("C"), Some("T"), false, false, None),
        ]
        .into_iter()
        .map(|call| (format!("{}:{}", call.chromosome, call.position), call))
        .collect::<BTreeMap<_, _>>();

        let mut data = MatchData::new("Sample_1", "CYP2B6", &definition, &allele_map);
        assert_eq!(
            data.missing_positions
                .iter()
                .map(|position| position.position)
                .collect::<Vec<_>>(),
            [3]
        );

        data.marshall_haplotypes(&definition);
        let hap4b = data
            .haplotypes
            .iter()
            .find(|haplotype| haplotype.name == "*4b")
            .expect("*4b");
        assert_eq!(hap4b.alleles, [Some("G".to_owned()), Some("T".to_owned())]);
        assert_eq!(
            hap4b
                .missing_positions
                .iter()
                .map(|position| position.position)
                .collect::<Vec<_>>(),
            [3]
        );
        assert_eq!(hap4b.score_override, Some(2));

        data.default_missing_alleles_to_reference()
            .expect("reference haplotype");
        let hap4a = data
            .haplotypes
            .iter()
            .find(|haplotype| haplotype.name == "*4a")
            .expect("*4a");
        assert_eq!(hap4a.alleles, [Some("G".to_owned()), Some("C".to_owned())]);
        assert_eq!(hap4a.score_override, Some(1));
    }

    #[test]
    fn wobble_alleles_do_not_score_when_all_supporting_sequences_are_reference() {
        let definition = synthetic_wobble_definition();
        let names_and_scores = compute_diplotype_names_and_scores(
            &definition,
            [sample_call(
                "chr1",
                1,
                Some("A"),
                Some("A"),
                false,
                true,
                None,
            )],
            false,
        );

        assert_eq!(
            names_and_scores,
            [
                ("*1/*1".to_owned(), 2),
                ("*1/*wobble".to_owned(), 1),
                ("*wobble/*wobble".to_owned(), 0),
            ]
        );
    }

    #[test]
    fn computes_combination_matches_and_partials_like_java_combination_matcher() {
        let definition = synthetic_combination_definition();
        let names = compute_combination_diplotype_names(
            &definition,
            [
                sample_call("chr1", 1, Some("T"), Some("T"), false, true, None),
                sample_call("chr1", 2, Some("T"), Some("T"), false, true, None),
                sample_call("chr1", 3, Some("T"), Some("T"), false, true, None),
            ],
            true,
        );

        assert_eq!(names, ["[*2 + *5 + g.3C>T]/[*2 + *5 + g.3C>T]"]);
    }

    #[test]
    fn computes_combination_baseline_from_java_fixture() {
        assert_eq!(
            compute_combination_fixture_names("NamedAlleleMatcher-combinationBaseline.vcf"),
            ["*1/*1"]
        );
    }

    #[test]
    fn computes_combination_phased_from_java_fixture() {
        assert_eq!(
            compute_combination_fixture_names("NamedAlleleMatcher-combinationPhased.vcf"),
            ["*1/[*6 + *27 + *28 + *80]"]
        );
    }

    #[test]
    fn computes_combination_unphased_from_java_fixture() {
        assert_eq!(
            compute_combination_fixture_names("NamedAlleleMatcher-combinationUnphased.vcf"),
            [
                "*1/[*6 + *27 + *28 + *80]",
                "*6/[*27 + *28 + *80]",
                "*27/[*6 + *28 + *80]",
                "*28/[*6 + *27 + *80]",
                "*80/[*6 + *27 + *28]",
                "[*6 + *27]/[*28 + *80]",
                "[*6 + *28]/[*27 + *80]",
                "[*6 + *80]/[*27 + *28]",
            ]
        );
    }

    #[test]
    fn computes_partial_with_combination_from_java_fixture() {
        assert_eq!(
            compute_combination_fixture_names("NamedAlleleMatcher-partialWithCombination.vcf"),
            [
                "*1/[*6 + *28 + g.233760973C>T]",
                "g.233760973C>T/[*6 + *28]",
                "*6/[*28 + g.233760973C>T]",
                "*28/[*6 + g.233760973C>T]",
            ]
        );
    }

    #[test]
    fn computes_partial_from_java_fixture() {
        assert_eq!(
            compute_combination_fixture_names("NamedAlleleMatcher-partial.vcf"),
            ["*6/[*6 + g.233760973C>T]"]
        );
    }

    #[test]
    fn dpyd_hap_b3_matcher_detects_hap_b3_and_non_hap_b3_variants() {
        let definition = read_dpyd_definition();
        let non_hap_variant = definition
            .variants
            .iter()
            .find(|variant| {
                variant.rsid.as_deref() != Some(DpydHapB3Matcher::HAP_B3_EXONIC_RSID)
                    && variant.rsid.as_deref() != Some(DpydHapB3Matcher::HAP_B3_INTRONIC_RSID)
                    && !variant.alts.is_empty()
            })
            .expect("non-HapB3 variant");
        let allele_map = allele_map([
            sample_call("chr1", 97573863, Some("C"), Some("T"), false, true, None),
            sample_call(
                "chr1",
                non_hap_variant.position as usize,
                Some(&non_hap_variant.reference),
                Some(&non_hap_variant.alts[0]),
                false,
                true,
                None,
            ),
        ]);
        let data = dpyd_match_data(&definition, &allele_map);

        let matcher =
            DpydHapB3Matcher::new(&definition, &allele_map, &data).expect("HapB3 matcher");
        assert!(matcher.has_hap_b3_variants());
        assert!(matcher.has_non_hap_b3_variants());
    }

    #[test]
    fn dpyd_hap_b3_matcher_calls_phased_full_hap_b3_like_java() {
        let definition = read_dpyd_definition();
        let allele_map = allele_map([
            sample_call("chr1", 97573863, Some("C"), Some("T"), true, true, None),
            sample_call("chr1", 97579893, Some("G"), Some("C"), true, true, None),
        ]);
        let data = dpyd_match_data(&definition, &allele_map);
        let mut matcher =
            DpydHapB3Matcher::new(&definition, &allele_map, &data).expect("HapB3 matcher");

        matcher
            .call_hap_b3_haplotype_matches()
            .expect("HapB3 calls");

        assert_eq!(
            matcher.hap_b3_intron_call().expect("intronic calls"),
            ["0", "0"]
        );
        assert_eq!(matcher.hap_b3_call().expect("HapB3 calls"), ["0", "1"]);
        assert_eq!(matcher.num_hap_b3_called(), 1);
        assert!(matcher.warnings().is_empty());
    }

    #[test]
    fn dpyd_hap_b3_matcher_warns_when_phased_exon_conflicts_with_reference_intron() {
        let definition = read_dpyd_definition();
        let allele_map = allele_map([
            sample_call("chr1", 97573863, Some("C"), Some("T"), true, true, None),
            sample_call("chr1", 97579893, Some("G"), Some("G"), true, true, None),
        ]);
        let data = dpyd_match_data(&definition, &allele_map);
        let mut matcher =
            DpydHapB3Matcher::new(&definition, &allele_map, &data).expect("HapB3 matcher");

        matcher
            .call_hap_b3_haplotype_matches()
            .expect("HapB3 calls");

        assert_eq!(
            matcher.hap_b3_intron_call().expect("intronic calls"),
            ["0", "0"]
        );
        assert_eq!(matcher.hap_b3_call().expect("HapB3 calls"), ["0", "0"]);
        assert_eq!(matcher.num_hap_b3_called(), 0);
        assert!(
            matcher
                .warnings()
                .contains(&DpydHapB3Warning::IntronicMismatchExonic)
        );
    }

    #[test]
    fn dpyd_hap_b3_matcher_calls_intronic_only_and_exonic_only_like_java() {
        let definition = read_dpyd_definition();

        let intron_only = allele_map([sample_call(
            "chr1",
            97579893,
            Some("G"),
            Some("C"),
            false,
            true,
            None,
        )]);
        let intron_data = dpyd_match_data(&definition, &intron_only);
        let mut intron_matcher =
            DpydHapB3Matcher::new(&definition, &intron_only, &intron_data).expect("HapB3 matcher");
        intron_matcher
            .call_hap_b3_haplotype_matches()
            .expect("HapB3 calls");
        assert_eq!(
            intron_matcher.hap_b3_intron_call().expect("intronic calls"),
            ["0", "1"]
        );
        assert_eq!(intron_matcher.num_hap_b3_called(), 1);
        assert!(intron_matcher.hap_b3_call().is_none());

        let exon_only = allele_map([sample_call(
            "chr1",
            97573863,
            Some("C"),
            Some("T"),
            false,
            true,
            None,
        )]);
        let exon_data = dpyd_match_data(&definition, &exon_only);
        let mut exon_matcher =
            DpydHapB3Matcher::new(&definition, &exon_only, &exon_data).expect("HapB3 matcher");
        exon_matcher
            .call_hap_b3_haplotype_matches()
            .expect("HapB3 calls");
        assert_eq!(exon_matcher.hap_b3_call().expect("HapB3 calls"), ["0", "1"]);
        assert_eq!(exon_matcher.num_hap_b3_called(), 1);
        assert!(
            exon_matcher
                .warnings()
                .contains(&DpydHapB3Warning::ExonicOnly)
        );
    }

    #[test]
    fn dpyd_hap_b3_matcher_uses_matching_phase_sets_like_java() {
        let definition = read_dpyd_definition();
        let allele_map = allele_map([
            sample_call(
                "chr1",
                97573863,
                Some("C"),
                Some("T"),
                false,
                false,
                Some(97571276),
            ),
            sample_call(
                "chr1",
                97579893,
                Some("G"),
                Some("C"),
                false,
                false,
                Some(97571276),
            ),
        ]);
        let data = dpyd_match_data(&definition, &allele_map);
        let mut matcher =
            DpydHapB3Matcher::new(&definition, &allele_map, &data).expect("HapB3 matcher");

        matcher
            .call_hap_b3_haplotype_matches()
            .expect("HapB3 calls");

        assert_eq!(matcher.hap_b3_call().expect("HapB3 calls"), ["0", "1"]);
        assert_eq!(matcher.num_hap_b3_called(), 1);
    }

    #[test]
    fn dpyd_hap_b3_matcher_adds_phased_hap_b3_to_reference_like_java() {
        let definition = read_dpyd_definition();
        let allele_map = allele_map([
            sample_call("chr1", 97573863, Some("C"), Some("T"), true, true, None),
            sample_call("chr1", 97579893, Some("G"), Some("C"), true, true, None),
        ]);
        let data = dpyd_match_data(&definition, &allele_map);
        let mut matcher =
            DpydHapB3Matcher::new(&definition, &allele_map, &data).expect("HapB3 matcher");

        let names = matcher
            .add_phased_hap_b3_call_to_ref(&data)
            .expect("HapB3 ref merge")
            .into_iter()
            .map(|diplotype| diplotype.name)
            .collect::<Vec<_>>();

        assert_eq!(names, ["Reference/c.1129-5923C>G, c.1236G>A (HapB3)"]);
    }

    #[test]
    fn dpyd_hap_b3_matcher_merges_phased_hap_b3_into_existing_matches_like_java() {
        let definition = read_dpyd_definition();
        let c61_position = definition
            .named_allele("c.61C>T")
            .expect("c.61C>T")
            .core_positions
            .iter()
            .next()
            .copied()
            .expect("c.61 position");
        let c61_variant = definition
            .variant_for_position(c61_position)
            .expect("c.61 variant");
        let allele_map = allele_map([
            sample_call(
                "chr1",
                c61_variant.position as usize,
                Some(&c61_variant.reference),
                Some(&c61_variant.alts[0]),
                true,
                true,
                None,
            ),
            sample_call("chr1", 97573863, Some("C"), Some("T"), true, true, None),
            sample_call("chr1", 97579893, Some("G"), Some("C"), true, true, None),
        ]);
        let data = dpyd_match_data(&definition, &allele_map);
        let sequences = data.permutations().iter().cloned().collect::<Vec<_>>();
        let reference = data
            .haplotypes
            .iter()
            .find(|haplotype| haplotype.name == "Reference")
            .cloned()
            .expect("Reference");
        let c61 = data
            .haplotypes
            .iter()
            .find(|haplotype| haplotype.name == "c.61C>T")
            .cloned()
            .expect("c.61C>T");
        let base_diplotype = DiplotypeMatch::new(
            HaplotypeMatch::from_haplotype(reference, data.positions.clone(), sequences[0].clone()),
            Some(HaplotypeMatch::from_haplotype(
                c61,
                data.positions.clone(),
                sequences[1].clone(),
            )),
            sequences.clone(),
        );
        let mut matcher =
            DpydHapB3Matcher::new(&definition, &allele_map, &data).expect("HapB3 matcher");

        let names = matcher
            .merge_phased_hap_b3_call(&data, &[base_diplotype])
            .expect("HapB3 merge")
            .into_iter()
            .map(|diplotype| diplotype.name)
            .collect::<Vec<_>>();

        assert_eq!(
            names,
            ["Reference/[c.61C>T + c.1129-5923C>G, c.1236G>A (HapB3)]"]
        );
    }

    #[test]
    fn dpyd_hap_b3_matcher_fix_partials_collapses_hap_b3_partial_reference_like_java() {
        let definition = read_dpyd_definition();
        let allele_map = allele_map([sample_call(
            "chr1",
            97573863,
            Some("C"),
            Some("T"),
            true,
            true,
            None,
        )]);
        let data = dpyd_match_data(&definition, &allele_map);
        let sequence = data
            .permutations()
            .iter()
            .next()
            .cloned()
            .expect("sequence");
        let reference = data
            .haplotypes
            .iter()
            .find(|haplotype| haplotype.name == "Reference")
            .cloned()
            .expect("Reference");
        let partial = build_combination_match(
            &data.positions,
            &sequence,
            std::slice::from_ref(&reference),
            [(
                DpydHapB3Matcher::HAP_B3_EXONIC_POSITION,
                DpydHapB3Matcher::HAP_B3_EXONIC_PARTIAL.to_owned(),
            )]
            .into_iter()
            .collect(),
        );
        let base_diplotype = DiplotypeMatch::new(
            partial,
            Some(HaplotypeMatch::from_haplotype(
                reference,
                data.positions.clone(),
                sequence.clone(),
            )),
            vec![sequence.clone(), sequence],
        );
        let matcher =
            DpydHapB3Matcher::new(&definition, &allele_map, &data).expect("HapB3 matcher");

        let names = matcher
            .fix_partials(&data, &[base_diplotype])
            .expect("fixed partials")
            .into_iter()
            .map(|diplotype| diplotype.name)
            .collect::<Vec<_>>();

        assert_eq!(names, ["Reference/Reference"]);
    }

    #[test]
    fn dpyd_hap_b3_matcher_fix_partials_keeps_non_hap_b3_components_like_java() {
        let definition = read_dpyd_definition();
        let c61_position = definition
            .named_allele("c.61C>T")
            .expect("c.61C>T")
            .core_positions
            .iter()
            .next()
            .copied()
            .expect("c.61 position");
        let c61_variant = definition
            .variant_for_position(c61_position)
            .expect("c.61 variant");
        let allele_map = allele_map([
            sample_call(
                "chr1",
                c61_variant.position as usize,
                Some(&c61_variant.reference),
                Some(&c61_variant.alts[0]),
                true,
                true,
                None,
            ),
            sample_call("chr1", 97573863, Some("C"), Some("T"), true, true, None),
        ]);
        let data = dpyd_match_data(&definition, &allele_map);
        let sequence = data
            .permutations()
            .iter()
            .next()
            .cloned()
            .expect("sequence");
        let reference = data
            .haplotypes
            .iter()
            .find(|haplotype| haplotype.name == "Reference")
            .cloned()
            .expect("Reference");
        let c61 = data
            .haplotypes
            .iter()
            .find(|haplotype| haplotype.name == "c.61C>T")
            .cloned()
            .expect("c.61C>T");
        let partial = build_combination_match(
            &data.positions,
            &sequence,
            &[c61],
            [(
                DpydHapB3Matcher::HAP_B3_EXONIC_POSITION,
                DpydHapB3Matcher::HAP_B3_EXONIC_PARTIAL.to_owned(),
            )]
            .into_iter()
            .collect(),
        );
        let base_diplotype = DiplotypeMatch::new(
            HaplotypeMatch::from_haplotype(reference, data.positions.clone(), sequence.clone()),
            Some(partial),
            vec![sequence.clone(), sequence],
        );
        let matcher =
            DpydHapB3Matcher::new(&definition, &allele_map, &data).expect("HapB3 matcher");

        let names = matcher
            .fix_partials(&data, &[base_diplotype])
            .expect("fixed partials")
            .into_iter()
            .map(|diplotype| diplotype.name)
            .collect::<Vec<_>>();

        assert_eq!(names, ["Reference/c.61C>T"]);
    }

    #[test]
    fn computes_initial_dpyd_lowest_function_hap_b3_only_path_like_java() {
        let definition = read_dpyd_definition();
        let allele_map = allele_map([
            sample_call("chr1", 97573863, Some("C"), Some("T"), true, true, None),
            sample_call("chr1", 97579893, Some("G"), Some("C"), true, true, None),
        ]);

        let names = compute_dpyd_lowest_function_diplotypes("Sample_1", &definition, &allele_map)
            .expect("DPYD diplotypes")
            .into_iter()
            .map(|diplotype| diplotype.name)
            .collect::<Vec<_>>();

        assert_eq!(names, ["Reference/c.1129-5923C>G, c.1236G>A (HapB3)"]);
    }

    #[test]
    fn computes_initial_dpyd_lowest_function_phased_hap_b3_merge_like_java() {
        let definition = read_dpyd_definition();
        let c61_position = definition
            .named_allele("c.61C>T")
            .expect("c.61C>T")
            .core_positions
            .iter()
            .next()
            .copied()
            .expect("c.61 position");
        let c61_variant = definition
            .variant_for_position(c61_position)
            .expect("c.61 variant");
        let allele_map = allele_map([
            sample_call(
                "chr1",
                c61_variant.position as usize,
                Some(&c61_variant.reference),
                Some(&c61_variant.alts[0]),
                true,
                true,
                None,
            ),
            sample_call("chr1", 97573863, Some("C"), Some("T"), true, true, None),
            sample_call("chr1", 97579893, Some("G"), Some("C"), true, true, None),
        ]);

        let names = compute_dpyd_lowest_function_diplotypes("Sample_1", &definition, &allele_map)
            .expect("DPYD diplotypes")
            .into_iter()
            .map(|diplotype| diplotype.name)
            .collect::<Vec<_>>();

        assert_eq!(
            names,
            ["Reference/[c.61C>T + c.1129-5923C>G, c.1236G>A (HapB3)]"]
        );
    }

    #[test]
    fn lowest_function_haplotype_fallback_strips_reference_when_more_than_two_matches_like_java() {
        let definition = synthetic_lowest_function_definition();
        let allele_map = allele_map([
            sample_call("chr1", 1, Some("C"), Some("T"), false, false, None),
            sample_call("chr1", 2, Some("C"), Some("T"), false, false, None),
            sample_call("chr1", 3, Some("C"), Some("C"), false, true, None),
        ]);
        let mut data = MatchData::new("Sample_1", "COMBO", &definition, &allele_map);
        data.marshall_haplotypes(&definition);
        data.generate_sample_permutations().expect("permutations");

        let names = call_haplotypes_for_lowest_function_gene(&data, &[], None)
            .expect("lowest-function haplotypes")
            .into_iter()
            .map(|haplotype| haplotype.name)
            .collect::<Vec<_>>();

        assert_eq!(names, ["*2", "*5"]);
    }

    #[test]
    fn lowest_function_haplotype_fallback_expands_homozygous_components_like_java() {
        let definition = synthetic_lowest_function_definition();
        let allele_map = allele_map([
            sample_call("chr1", 1, Some("C"), Some("T"), false, false, None),
            sample_call("chr1", 2, Some("C"), Some("T"), false, false, None),
            sample_call("chr1", 3, Some("C"), Some("C"), false, true, None),
        ]);
        let mut data = MatchData::new("Sample_1", "COMBO", &definition, &allele_map);
        data.marshall_haplotypes(&definition);
        data.generate_sample_permutations().expect("permutations");
        let star2 = data
            .compare_permutations()
            .into_iter()
            .find(|haplotype| haplotype.name == "*2")
            .expect("*2 match");
        let sequence = star2.sequences.iter().next().cloned().expect("*2 sequence");
        let diplotype =
            DiplotypeMatch::new(star2.clone(), Some(star2), vec![sequence.clone(), sequence]);

        let names = call_haplotypes_for_lowest_function_gene(&data, &[diplotype], None)
            .expect("lowest-function haplotypes")
            .into_iter()
            .map(|haplotype| haplotype.name)
            .collect::<Vec<_>>();

        assert_eq!(names, ["*2", "*2", "*5"]);
    }

    #[test]
    fn lowest_function_haplotype_fallback_appends_dpyd_hap_b3_matches_like_java() {
        let definition = read_dpyd_definition();
        let combo_definition = definition_without_dpyd_hap_b3(&definition);
        let allele_map = allele_map([sample_call(
            "chr1",
            97573863,
            Some("C"),
            Some("T"),
            false,
            false,
            None,
        )]);
        let orig_data = dpyd_match_data(&definition, &allele_map);
        let data = dpyd_match_data(&combo_definition, &allele_map);
        let mut matcher =
            DpydHapB3Matcher::new(&definition, &allele_map, &orig_data).expect("HapB3 matcher");

        let names = call_haplotypes_for_lowest_function_gene(&data, &[], Some(&mut matcher))
            .expect("lowest-function haplotypes")
            .into_iter()
            .map(|haplotype| haplotype.name)
            .collect::<Vec<_>>();

        assert_eq!(names, ["Reference", "c.1129-5923C>G, c.1236G>A (HapB3)"]);
    }

    #[test]
    fn call_dpyd_lowest_function_gene_returns_diplotypes_for_phased_hap_b3_like_java() {
        let definition = read_dpyd_definition();
        let allele_map = allele_map([
            sample_call("chr1", 97573863, Some("C"), Some("T"), true, true, None),
            sample_call("chr1", 97579893, Some("G"), Some("C"), true, true, None),
        ]);

        let result = call_dpyd_lowest_function_gene("Sample_1", &definition, &allele_map)
            .expect("DPYD call");

        let GeneCallKind::Diplotypes(diplotypes) = result.kind else {
            panic!("expected diplotype call");
        };
        assert_eq!(result.gene, "DPYD");
        assert!(result.match_data.effectively_phased);
        assert!(result.dpyd_hap_b3_warnings.is_empty());
        assert_eq!(
            diplotypes
                .into_iter()
                .map(|diplotype| diplotype.name)
                .collect::<Vec<_>>(),
            ["Reference/c.1129-5923C>G, c.1236G>A (HapB3)"]
        );
    }

    #[test]
    fn dpyd_effectively_phased_fixture_matches_java_named_allele_matcher() {
        let definition = read_dpyd_definition();
        let allele_map = allele_map([sample_call(
            "chr1",
            97883352,
            Some("C"),
            Some("T"),
            false,
            true,
            None,
        )]);

        let result = call_dpyd_lowest_function_gene("Sample_1", &definition, &allele_map)
            .expect("DPYD call");

        assert_eq!(result.gene, "DPYD");
        assert!(result.dpyd_hap_b3_warnings.is_empty());
        let GeneCallKind::Diplotypes(diplotypes) = result.kind else {
            panic!("expected diplotype call");
        };
        assert_eq!(diplotypes.len(), 1);
        let diplotype = diplotypes.first().expect("diplotype");
        assert_eq!(diplotype.name, "Reference/c.62G>A");
        assert_eq!(
            [
                diplotype.haplotype1.name.as_str(),
                diplotype
                    .haplotype2
                    .as_ref()
                    .expect("second haplotype")
                    .name
                    .as_str(),
            ],
            ["Reference", "c.62G>A"]
        );
    }

    #[test]
    fn dpyd_effectively_phased_combination_fixture_matches_java_named_allele_matcher() {
        let definition = read_dpyd_definition();
        let allele_map = allele_map([
            sample_call("chr1", 97573881, Some("T"), Some("T"), false, true, None),
            sample_call("chr1", 97515839, Some("C"), Some("C"), false, true, None),
            sample_call("chr1", 97883329, Some("A"), Some("G"), false, true, None),
        ]);

        let result = call_dpyd_lowest_function_gene("Sample_1", &definition, &allele_map)
            .expect("DPYD call");

        assert_eq!(result.gene, "DPYD");
        assert!(result.dpyd_hap_b3_warnings.is_empty());
        let GeneCallKind::Diplotypes(diplotypes) = result.kind else {
            panic!("expected diplotype call");
        };
        assert_eq!(diplotypes.len(), 1);
        let diplotype = diplotypes.first().expect("diplotype");
        assert_eq!(
            diplotype.name,
            "[c.1218G>A + c.1627A>G (*5)]/[c.85T>C (*9A) + c.1218G>A + c.1627A>G (*5)]"
        );
        let component_names = [&diplotype.haplotype1.name]
            .into_iter()
            .chain(
                diplotype
                    .haplotype2
                    .as_ref()
                    .map(|haplotype| &haplotype.name),
            )
            .flat_map(|name| {
                if super::is_combination_name(name) {
                    super::split_combination_name(name)
                } else {
                    vec![name.clone()]
                }
            })
            .collect::<BTreeSet<_>>();
        assert_eq!(
            component_names.into_iter().collect::<Vec<_>>(),
            ["c.1218G>A", "c.1627A>G (*5)", "c.85T>C (*9A)"]
        );
    }

    #[test]
    fn dpyd_phased_fixture_matches_java_named_allele_matcher() {
        let definition = read_dpyd_definition();
        let allele_map = allele_map([
            sample_call("chr1", 97883352, Some("C"), Some("T"), true, true, None),
            sample_call("chr1", 97078987, Some("G"), Some("T"), true, true, None),
        ]);

        let result = call_dpyd_lowest_function_gene("Sample_1", &definition, &allele_map)
            .expect("DPYD call");

        assert_eq!(result.gene, "DPYD");
        assert!(result.dpyd_hap_b3_warnings.is_empty());
        let GeneCallKind::Diplotypes(diplotypes) = result.kind else {
            panic!("expected diplotype call");
        };
        assert_eq!(diplotypes.len(), 1);
        let diplotype = diplotypes.first().expect("diplotype");
        assert_eq!(diplotype.name, "Reference/[c.62G>A + c.3067C>A]");
        assert_eq!(
            [
                diplotype.haplotype1.name.as_str(),
                diplotype
                    .haplotype2
                    .as_ref()
                    .expect("second haplotype")
                    .name
                    .as_str(),
            ],
            ["Reference", "[c.62G>A + c.3067C>A]"]
        );
    }

    #[test]
    fn dpyd_unphased_fixture_matches_java_named_allele_matcher() {
        let definition = read_dpyd_definition();
        let allele_map = allele_map([
            sample_call("chr1", 97883352, Some("C"), Some("T"), false, false, None),
            sample_call("chr1", 97078987, Some("G"), Some("T"), false, false, None),
        ]);

        let result = call_dpyd_lowest_function_gene("Sample_1", &definition, &allele_map)
            .expect("DPYD call");

        assert_eq!(result.gene, "DPYD");
        assert!(result.dpyd_hap_b3_warnings.is_empty());
        let GeneCallKind::Haplotypes(haplotypes) = result.kind else {
            panic!("expected haplotype fallback");
        };
        assert_eq!(
            haplotypes
                .into_iter()
                .map(|haplotype| haplotype.name)
                .collect::<Vec<_>>(),
            ["c.62G>A", "c.3067C>A"]
        );
    }

    #[test]
    fn dpyd_phased_double_fixture_matches_java_named_allele_matcher() {
        let definition = read_dpyd_definition();
        let allele_map = allele_map([
            sample_call("chr1", 97579893, Some("G"), Some("C"), true, true, None),
            sample_call("chr1", 97573863, Some("C"), Some("T"), true, true, None),
            sample_call("chr1", 97883352, Some("T"), Some("T"), true, true, None),
            sample_call("chr1", 97078987, Some("T"), Some("T"), true, true, None),
        ]);

        let result = call_dpyd_lowest_function_gene("Sample_1", &definition, &allele_map)
            .expect("DPYD call");

        assert_eq!(result.gene, "DPYD");
        assert!(result.dpyd_hap_b3_warnings.is_empty());
        let GeneCallKind::Diplotypes(diplotypes) = result.kind else {
            panic!("expected diplotype call");
        };
        assert_eq!(diplotypes.len(), 1);
        let diplotype = diplotypes.first().expect("diplotype");
        assert_eq!(
            diplotype.name,
            "[c.62G>A + c.3067C>A]/[c.62G>A + c.1129-5923C>G, c.1236G>A (HapB3) + c.3067C>A]"
        );
        assert_eq!(
            [
                diplotype.haplotype1.name.as_str(),
                diplotype
                    .haplotype2
                    .as_ref()
                    .expect("second haplotype")
                    .name
                    .as_str(),
            ],
            [
                "[c.62G>A + c.3067C>A]",
                "[c.62G>A + c.1129-5923C>G, c.1236G>A (HapB3) + c.3067C>A]",
            ]
        );
    }

    #[test]
    fn dpyd_homozygous_fixture_matches_java_named_allele_matcher() {
        let definition = read_dpyd_definition();
        let allele_map = allele_map([
            sample_call("chr1", 97883352, Some("T"), Some("T"), true, true, None),
            sample_call("chr1", 97078987, Some("T"), Some("T"), true, true, None),
        ]);

        let result = call_dpyd_lowest_function_gene("Sample_1", &definition, &allele_map)
            .expect("DPYD call");

        assert_eq!(result.gene, "DPYD");
        assert!(result.dpyd_hap_b3_warnings.is_empty());
        let GeneCallKind::Diplotypes(diplotypes) = result.kind else {
            panic!("expected diplotype call");
        };
        assert_eq!(diplotypes.len(), 1);
        let diplotype = diplotypes.first().expect("diplotype");
        assert_eq!(
            diplotype.name,
            "[c.62G>A + c.3067C>A]/[c.62G>A + c.3067C>A]"
        );
        assert_eq!(
            [
                diplotype.haplotype1.name.as_str(),
                diplotype
                    .haplotype2
                    .as_ref()
                    .expect("second haplotype")
                    .name
                    .as_str(),
            ],
            ["[c.62G>A + c.3067C>A]", "[c.62G>A + c.3067C>A]"]
        );
    }

    #[test]
    fn dpyd_effectively_phased_combination_missing_fixture_matches_java_named_allele_matcher() {
        let definition = read_dpyd_definition();
        let mut allele_map = reference_allele_map(&definition);
        for rsid in [
            "rs148799944",
            "rs140114515",
            "rs1801268",
            "rs72547601",
            "rs72547602",
            "rs141044036",
            "rs147545709",
            "rs55674432",
            "rs146529561",
            "rs137999090",
            "rs138545885",
            "rs55971861",
            "rs72549303",
            "rs147601618",
            "rs145773863",
            "rs138616379",
            "rs148994843",
            "rs138391898",
            "rs111858276",
            "rs72549304",
            "rs142512579",
            "rs764666241",
            "rs140602333",
            "rs78060119",
            "rs143154602",
            "rs72549306",
            "rs145112791",
            "rs150437414",
            "rs1801266",
            "rs72549307",
            "rs72549308",
            "rs139834141",
            "rs141462178",
            "rs150385342",
            "rs72549309",
            "rs80081766",
            "rs150036960",
        ] {
            allele_map.remove(&variant_by_rsid(&definition, rsid).vcf_chr_position());
        }
        allele_map.insert(
            "chr1:97515839".to_owned(),
            sample_call("chr1", 97515839, Some("C"), Some("C"), false, true, None),
        );
        allele_map.insert(
            "chr1:97883329".to_owned(),
            sample_call("chr1", 97883329, Some("A"), Some("G"), false, true, None),
        );

        let result = call_dpyd_lowest_function_gene("Sample_1", &definition, &allele_map)
            .expect("DPYD call");

        assert_eq!(result.gene, "DPYD");
        assert!(result.dpyd_hap_b3_warnings.is_empty());
        let GeneCallKind::Diplotypes(diplotypes) = result.kind else {
            panic!("expected diplotype call");
        };
        assert_eq!(diplotypes.len(), 1);
        let diplotype = diplotypes.first().expect("diplotype");
        assert_eq!(
            diplotype.name,
            "c.1627A>G (*5)/[c.85T>C (*9A) + c.1627A>G (*5)]"
        );
        assert_eq!(
            [
                diplotype.haplotype1.name.as_str(),
                diplotype
                    .haplotype2
                    .as_ref()
                    .expect("second haplotype")
                    .name
                    .as_str(),
            ],
            ["c.1627A>G (*5)", "[c.85T>C (*9A) + c.1627A>G (*5)]",]
        );
    }

    #[test]
    fn dpyd_unphased_homozygous_no_function_fixture_matches_java_named_allele_matcher() {
        let definition = read_dpyd_definition();
        let mut allele_map = reference_allele_map(&definition);
        for (rsid, allele1, allele2) in [
            ("rs72547601", "C", "C"),
            ("rs67376798", "A", "T"),
            ("rs60139309", "T", "C"),
        ] {
            let variant = variant_by_rsid(&definition, rsid);
            allele_map.insert(
                variant.vcf_chr_position(),
                sample_call(
                    &variant.chromosome,
                    variant.position as usize,
                    Some(allele1),
                    Some(allele2),
                    false,
                    false,
                    None,
                ),
            );
        }

        let result = call_dpyd_lowest_function_gene("Sample_1", &definition, &allele_map)
            .expect("DPYD call");

        assert_eq!(result.gene, "DPYD");
        assert!(result.dpyd_hap_b3_warnings.is_empty());
        let GeneCallKind::Haplotypes(haplotypes) = result.kind else {
            panic!("expected haplotype fallback");
        };
        assert_eq!(
            haplotypes
                .into_iter()
                .map(|haplotype| haplotype.name)
                .collect::<Vec<_>>(),
            ["c.2582A>G", "c.2846A>T", "c.2933A>G", "c.2933A>G"]
        );
    }

    #[test]
    fn cyp2c19_permutation_generation_fixture_matches_java_named_allele_matcher() {
        let definition = read_cyp2c19_definition();

        for (case, phased_overrides) in [
            ("phased", [true, true, true, true]),
            ("unphased", [false, false, false, false]),
            ("mixed", [false, true, true, false]),
        ] {
            let allele_map = cyp2c19_permutation_allele_map(&definition, phased_overrides);
            let result = call_standard_gene("Sample_1", &definition, &allele_map, true, false)
                .unwrap_or_else(|err| panic!("{case}: CYP2C19 call failed: {err}"));
            let GeneCallKind::Diplotypes(diplotypes) = result.kind else {
                panic!("{case}: expected diplotype call");
            };

            assert_eq!(
                diplotypes
                    .into_iter()
                    .map(|diplotype| diplotype.name)
                    .collect::<Vec<_>>(),
                ["*2/*17"],
                "{case}"
            );
        }
    }

    #[test]
    fn cyp2b6_partial_missing_allele_combination_fixture_matches_java_named_allele_matcher() {
        let definition = read_cyp2b6_definition();
        let variant = variant_by_rsid(&definition, "rs8192709");
        let mut allele_map = reference_allele_map(&definition);
        allele_map.insert(
            variant.vcf_chr_position(),
            sample_call(
                &variant.chromosome,
                variant.position as usize,
                Some("T"),
                Some("."),
                false,
                false,
                None,
            ),
        );

        let result = call_standard_gene("Sample_1", &definition, &allele_map, true, true)
            .expect("CYP2B6 call");

        assert!(result.match_data.has_partial_missing_alleles());
        let GeneCallKind::Diplotypes(diplotypes) = result.kind else {
            panic!("expected diplotype call");
        };
        assert_eq!(
            diplotypes
                .into_iter()
                .map(|diplotype| diplotype.name)
                .collect::<Vec<_>>(),
            ["*2/g.40991369?"]
        );
    }

    #[test]
    fn cyp2b6_phased_partial_missing_allele_combination_fixture_matches_java_named_allele_matcher()
    {
        let definition = read_cyp2b6_definition();
        let variant = variant_by_rsid(&definition, "rs8192709");
        let mut allele_map = reference_allele_map(&definition);
        for call in allele_map.values_mut() {
            call.phased = true;
            call.genotype = "0|0".to_owned();
            call.vcf_call = format!(
                "{}|{}",
                call.allele1.as_deref().unwrap_or("."),
                call.allele2.as_deref().unwrap_or(".")
            );
        }
        allele_map.insert(
            variant.vcf_chr_position(),
            sample_call(
                &variant.chromosome,
                variant.position as usize,
                Some("."),
                Some("T"),
                true,
                true,
                None,
            ),
        );

        let result = call_standard_gene("Sample_1", &definition, &allele_map, true, true)
            .expect("CYP2B6 call");

        assert!(result.match_data.has_partial_missing_alleles());
        let GeneCallKind::Diplotypes(diplotypes) = result.kind else {
            panic!("expected diplotype call");
        };
        assert_eq!(
            diplotypes
                .into_iter()
                .map(|diplotype| diplotype.name)
                .collect::<Vec<_>>(),
            ["*2/g.40991369?"]
        );
    }

    #[test]
    fn cyp2b6_wobble_scoring_multiple_sequence_fixture_matches_java_named_allele_matcher() {
        let definition = read_cyp2b6_definition();
        let allele_map = cyp2b6_wobble_scoring_allele_map(&definition);

        let result = call_standard_gene("Sample_1", &definition, &allele_map, true, false)
            .expect("CYP2B6 call");

        assert!(result.warnings.is_empty());
        let GeneCallKind::Diplotypes(diplotypes) = result.kind else {
            panic!("expected diplotype call");
        };
        assert_eq!(
            diplotypes
                .into_iter()
                .map(|diplotype| diplotype.name)
                .collect::<Vec<_>>(),
            ["*6/*18", "*9/*18"]
        );
    }

    #[test]
    fn nat2_unphased_priority_same_score_fixture_matches_java_named_allele_matcher() {
        let definition = read_nat2_definition();
        let allele_map = nat2_unphased_priority_allele_map(
            &definition,
            [
                ("rs1799930", "G", "A"),
                ("rs1208", "G", "A"),
                ("rs1799931", "G", "A"),
            ],
        );

        let without_priority = call_standard_gene("Sample_1", &definition, &allele_map, true, true)
            .expect("NAT2 call without priority");
        assert!(without_priority.warnings.is_empty());
        let GeneCallKind::Diplotypes(diplotypes) = without_priority.kind else {
            panic!("expected diplotype call without priority");
        };
        assert_eq!(
            diplotypes
                .into_iter()
                .map(|diplotype| diplotype.name)
                .collect::<Vec<_>>(),
            ["*6/*40", "*7/*34"]
        );

        let exemptions = read_default_exemptions();
        let with_priority = call_standard_gene_with_exemption(
            "Sample_1",
            &definition,
            exemptions.get("NAT2"),
            &allele_map,
            true,
            true,
        )
        .expect("NAT2 call with priority");
        assert!(
            with_priority
                .warnings
                .contains(&GeneCallWarning::UnphasedPriority)
        );
        let GeneCallKind::Diplotypes(diplotypes) = with_priority.kind else {
            panic!("expected diplotype call with priority");
        };
        assert_eq!(
            diplotypes
                .into_iter()
                .map(|diplotype| diplotype.name)
                .collect::<Vec<_>>(),
            ["*7/*34"]
        );
    }

    #[test]
    fn nat2_unphased_priority_different_score_fixture_matches_java_named_allele_matcher() {
        let definition = read_nat2_definition();
        let allele_map = nat2_unphased_priority_allele_map(
            &definition,
            [("rs1799930", "G", "A"), ("rs1208", "G", "A")],
        );

        let without_priority = call_standard_gene("Sample_1", &definition, &allele_map, true, true)
            .expect("NAT2 call without priority");
        assert!(without_priority.warnings.is_empty());
        let GeneCallKind::Diplotypes(diplotypes) = without_priority.kind else {
            panic!("expected diplotype call without priority");
        };
        assert_eq!(
            diplotypes
                .into_iter()
                .map(|diplotype| diplotype.name)
                .collect::<Vec<_>>(),
            ["*1/*6", "*4/*34"]
        );

        let exemptions = read_default_exemptions();
        let with_priority = call_standard_gene_with_exemption(
            "Sample_1",
            &definition,
            exemptions.get("NAT2"),
            &allele_map,
            true,
            true,
        )
        .expect("NAT2 call with priority");
        assert_eq!(
            with_priority.warnings.iter().collect::<Vec<_>>(),
            [&GeneCallWarning::UnphasedPriority]
        );
        let GeneCallKind::Diplotypes(diplotypes) = with_priority.kind else {
            panic!("expected diplotype call with priority");
        };
        assert_eq!(
            diplotypes
                .into_iter()
                .map(|diplotype| diplotype.name)
                .collect::<Vec<_>>(),
            ["*1/*6"]
        );
    }

    #[test]
    fn default_reference_filled_alleles_keep_original_scores_for_slco1b1_like_java() {
        let definition = read_slco1b1_definition();
        let mut allele_map = reference_allele_map(&definition);
        for (rsid, allele1, allele2) in [("rs2306283", "A", "G"), ("rs4149056", "T", "C")] {
            let variant = variant_by_rsid(&definition, rsid);
            allele_map.insert(
                variant.vcf_chr_position(),
                sample_call(
                    &variant.chromosome,
                    variant.position as usize,
                    Some(allele1),
                    Some(allele2),
                    false,
                    false,
                    None,
                ),
            );
        }

        let result = call_standard_gene("Sample_1", &definition, &allele_map, true, false)
            .expect("SLCO1B1 call");
        let GeneCallKind::Diplotypes(diplotypes) = result.kind else {
            panic!("expected diplotype call");
        };

        assert_eq!(
            diplotypes
                .into_iter()
                .map(|diplotype| (diplotype.name, diplotype.score))
                .collect::<Vec<_>>(),
            [("*1/*15".to_owned(), 37)]
        );
    }

    #[test]
    fn nat2_required_position_fixture_matches_java_named_allele_matcher() {
        let definition = read_nat2_definition();
        let missing_variant = variant_by_rsid(&definition, "rs1801279");
        let mut allele_map = reference_allele_map(&definition);
        allele_map.remove(&missing_variant.vcf_chr_position());

        let without_required_positions =
            call_standard_gene("Sample_1", &definition, &allele_map, true, true)
                .expect("NAT2 call without required positions");
        assert!(without_required_positions.warnings.is_empty());
        let GeneCallKind::Diplotypes(diplotypes) = without_required_positions.kind else {
            panic!("expected diplotype call without required positions");
        };
        assert_eq!(
            diplotypes
                .into_iter()
                .map(|diplotype| diplotype.name)
                .collect::<Vec<_>>(),
            ["*1/*1"]
        );

        let exemptions = read_default_exemptions();
        let with_required_positions = call_standard_gene_with_exemption(
            "Sample_1",
            &definition,
            exemptions.get("NAT2"),
            &allele_map,
            true,
            true,
        )
        .expect("NAT2 call with required positions");
        assert!(matches!(with_required_positions.kind, GeneCallKind::NoCall));
        assert_eq!(
            with_required_positions
                .match_data
                .missing_required_positions,
            ["chr8:18400194"]
        );
        assert_eq!(
            with_required_positions.warnings.iter().collect::<Vec<_>>(),
            [&GeneCallWarning::MissingRequiredPosition(vec![
                "chr8:18400194".to_owned()
            ])]
        );
    }

    #[test]
    fn nat2_combination_fixture_matches_java_named_allele_matcher() {
        let definition = read_nat2_definition();
        let allele_map = nat2_phased_combination_allele_map(&definition);

        let result = call_standard_gene("Sample_1", &definition, &allele_map, true, true)
            .expect("NAT2 combination call");

        assert!(result.warnings.is_empty());
        let GeneCallKind::Diplotypes(diplotypes) = result.kind else {
            panic!("expected diplotype call");
        };
        assert_eq!(
            diplotypes
                .into_iter()
                .map(|diplotype| diplotype.name)
                .collect::<Vec<_>>(),
            ["*1/[*15 + *44]", "*1/[*36 + *46]"]
        );
    }

    #[test]
    fn call_dpyd_lowest_function_gene_returns_haplotype_fallback_for_unphased_hap_b3_like_java() {
        let definition = read_dpyd_definition();
        let c61_position = definition
            .named_allele("c.61C>T")
            .expect("c.61C>T")
            .core_positions
            .iter()
            .next()
            .copied()
            .expect("c.61 position");
        let c61_variant = definition
            .variant_for_position(c61_position)
            .expect("c.61 variant");
        let allele_map = allele_map([
            sample_call(
                "chr1",
                c61_variant.position as usize,
                Some(&c61_variant.reference),
                Some(&c61_variant.alts[0]),
                false,
                false,
                None,
            ),
            sample_call("chr1", 97573863, Some("C"), Some("T"), false, false, None),
        ]);

        let result = call_dpyd_lowest_function_gene("Sample_1", &definition, &allele_map)
            .expect("DPYD call");

        let GeneCallKind::Haplotypes(haplotypes) = result.kind else {
            panic!("expected haplotype fallback");
        };
        assert_eq!(result.gene, "DPYD");
        assert!(!result.match_data.effectively_phased);
        assert!(
            result
                .dpyd_hap_b3_warnings
                .contains(&DpydHapB3Warning::ExonicOnly)
        );
        assert_eq!(
            haplotypes
                .into_iter()
                .map(|haplotype| haplotype.name)
                .collect::<Vec<_>>(),
            ["c.61C>T", "c.1129-5923C>G, c.1236G>A (HapB3)"]
        );
    }

    #[test]
    fn call_dpyd_lowest_function_gene_returns_no_call_without_sample_alleles_like_java() {
        let definition = read_dpyd_definition();
        let allele_map = BTreeMap::new();

        let result = call_dpyd_lowest_function_gene("Sample_1", &definition, &allele_map)
            .expect("DPYD call");

        assert_eq!(result.gene, "DPYD");
        assert!(matches!(result.kind, GeneCallKind::NoCall));
        assert!(result.match_data.permutations().is_empty());
    }

    #[test]
    fn call_ryr1_lowest_function_gene_returns_reference_diplotype_like_java() {
        let definition = read_ryr1_definition();
        let allele_map = reference_allele_map(&definition);

        let result = call_ryr1_lowest_function_gene("Sample_1", &definition, &allele_map)
            .expect("RYR1 call");

        let GeneCallKind::Diplotypes(diplotypes) = result.kind else {
            panic!("expected diplotype call");
        };
        assert_eq!(result.gene, "RYR1");
        assert!(result.match_data.effectively_phased);
        assert_eq!(result.match_data.missing_positions.len(), 0);
        assert_eq!(
            diplotypes
                .into_iter()
                .map(|diplotype| diplotype.name)
                .collect::<Vec<_>>(),
            ["Reference/Reference"]
        );
    }

    #[test]
    fn call_ryr1_lowest_function_gene_returns_heterozygous_diplotype_like_java() {
        let definition = read_ryr1_definition();
        let mut allele_map = reference_allele_map(&definition);
        let variant = variant_by_rsid(&definition, "rs137933390");
        allele_map.insert(
            variant.vcf_chr_position(),
            sample_call(
                &variant.chromosome,
                variant.position as usize,
                Some(&variant.reference),
                Some("G"),
                false,
                false,
                None,
            ),
        );

        let result = call_ryr1_lowest_function_gene("Sample_1", &definition, &allele_map)
            .expect("RYR1 call");

        let GeneCallKind::Diplotypes(diplotypes) = result.kind else {
            panic!("expected diplotype call");
        };
        assert_eq!(result.match_data.missing_positions.len(), 0);
        assert_eq!(
            diplotypes
                .into_iter()
                .map(|diplotype| diplotype.name)
                .collect::<Vec<_>>(),
            ["Reference/c.4178A>G"]
        );
    }

    #[test]
    fn call_ryr1_lowest_function_gene_preserves_missing_position_like_java() {
        let definition = read_ryr1_definition();
        let mut allele_map = reference_allele_map(&definition);
        let called_variant = variant_by_rsid(&definition, "rs137933390");
        allele_map.insert(
            called_variant.vcf_chr_position(),
            sample_call(
                &called_variant.chromosome,
                called_variant.position as usize,
                Some(&called_variant.reference),
                Some("G"),
                false,
                false,
                None,
            ),
        );
        let missing_variant = variant_by_rsid(&definition, "rs193922753");
        allele_map.remove(&missing_variant.vcf_chr_position());

        let result = call_ryr1_lowest_function_gene("Sample_1", &definition, &allele_map)
            .expect("RYR1 call");

        let GeneCallKind::Diplotypes(diplotypes) = result.kind else {
            panic!("expected diplotype call");
        };
        assert_eq!(
            result
                .match_data
                .missing_positions
                .iter()
                .map(VariantLocus::vcf_chr_position)
                .collect::<Vec<_>>(),
            ["chr19:38444212"]
        );
        assert_eq!(
            diplotypes
                .into_iter()
                .map(|diplotype| diplotype.name)
                .collect::<Vec<_>>(),
            ["Reference/c.4178A>G"]
        );
    }

    #[test]
    fn definition_aware_vcf_ref_mismatch_warnings_match_java_named_allele_matcher() {
        let definition = read_definition_file(std::path::Path::new(
            "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-mismatchedRefAllele.json",
        ))
        .expect("definition");
        let reader = crate::definition::DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition)]
                .into_iter()
                .collect(),
        );
        let records = read_record_summaries(
            "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-mismatchedRefAllele.vcf",
            Some("NA12878"),
        )
        .expect("vcf");

        let mut warnings = records.warnings.clone();
        let allele_calls =
            allele_calls_for_locations(&records, reader.locations_of_interest(), &mut warnings);

        assert_eq!(warnings.len(), 2);
        assert_eq!(
            warnings
                .get("chr10:94942205")
                .expect("warning at chr10:94942205")
                .iter()
                .collect::<Vec<_>>(),
            [&"Discarded genotype at this position because REF in VCF (C) does not match expected reference (CAATGGAAAGA)".to_owned()]
        );
        assert_eq!(
            warnings
                .get("chr10:94949281")
                .expect("warning at chr10:94949281")
                .iter()
                .collect::<Vec<_>>(),
            [&"Discarded genotype at this position because REF in VCF (G) does not match expected reference (GA)".to_owned()]
        );
        assert!(!allele_calls.iter().any(|call| call.position == 94942205));
        assert!(!allele_calls.iter().any(|call| call.position == 94949281));
    }

    #[test]
    fn cyp2d6_sort_exception_fixture_produces_many_diplotypes_like_java_named_allele_matcher() {
        let definitions_dir = std::path::Path::new(
            "../../repos/PharmCAT/src/main/resources/org/pharmgkb/pharmcat/definition/alleles",
        );
        let mut paths = std::fs::read_dir(definitions_dir)
            .expect("definition directory")
            .map(|entry| entry.map(|entry| entry.path()))
            .collect::<Result<Vec<_>, _>>()
            .expect("definition paths");
        paths.retain(|path| {
            path.file_name()
                .and_then(|name| name.to_str())
                .is_some_and(|name| name.ends_with("_translation.json"))
        });
        paths.sort();
        let reader = crate::definition::DefinitionReader::from_paths_with_exemptions(
            paths,
            definitions_dir.join("exemptions.json"),
        )
        .expect("definition reader");
        let definition = reader.definition_file("CYP2D6").expect("CYP2D6 definition");
        let records = read_record_summaries(
            "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-sortError.vcf",
            None,
        )
        .expect("vcf");

        let mut warnings = records.warnings.clone();
        let allele_calls =
            allele_calls_for_locations(&records, reader.locations_of_interest(), &mut warnings);
        let allele_map = allele_calls
            .into_iter()
            .map(|call| (format!("{}:{}", call.chromosome, call.position), call))
            .collect::<BTreeMap<_, _>>();
        let result = call_standard_gene_with_exemption(
            "2266231",
            definition,
            reader.exemption("CYP2D6"),
            &allele_map,
            true,
            true,
        )
        .expect("CYP2D6 call");

        assert_eq!(
            warnings.values().map(BTreeSet::len).sum::<usize>(),
            34,
            "expected Java SortedSetMultimap warning size"
        );
        let GeneCallKind::Diplotypes(diplotypes) = result.kind else {
            panic!("CYP2D6 result was not a diplotype call: {:?}", result.kind);
        };
        assert!(
            diplotypes.len() >= 10,
            "CYP2D6 diplotypes size was {}, expected at least 10",
            diplotypes.len()
        );
    }

    #[test]
    fn partial_reference_unphased_fixture_matches_java_named_allele_matcher() {
        let (result, warnings) = call_cyp2b6_partial_reference_fixture(
            "NamedAlleleMatcher-partialReferenceUnphased.vcf",
        );
        let GeneCallKind::Diplotypes(diplotypes) = result.kind else {
            panic!("CYP2B6 result was not a diplotype call: {:?}", result.kind);
        };

        assert_eq!(warnings.values().map(BTreeSet::len).sum::<usize>(), 1);
        assert_eq!(
            diplotypes
                .iter()
                .map(|diplotype| diplotype.name.as_str())
                .collect::<Vec<_>>(),
            ["*1/[*9 + g.41012316T>G]", "*9/g.41012316T>G"]
        );
        let first = diplotypes.first().expect("top diplotype");
        assert!(!first.haplotype1.haplotype.is_combination());
        assert!(!first.haplotype1.haplotype.is_partial());
        let second = first.haplotype2.as_ref().expect("second haplotype");
        assert!(!second.haplotype.is_combination());
        assert!(second.haplotype.is_partial());
    }

    #[test]
    fn partial_reference_phased_fixture_matches_java_named_allele_matcher() {
        let (result, warnings) =
            call_cyp2b6_partial_reference_fixture("NamedAlleleMatcher-partialReferencePhased.vcf");
        let GeneCallKind::Diplotypes(diplotypes) = result.kind else {
            panic!("CYP2B6 result was not a diplotype call: {:?}", result.kind);
        };

        assert_eq!(warnings.values().map(BTreeSet::len).sum::<usize>(), 2);
        assert_eq!(diplotypes.len(), 1);
        let diplotype = diplotypes.first().expect("diplotype");
        assert_eq!(diplotype.name, "*9/[g.41006968T>A + g.41012316T>G]");
        assert!(!diplotype.haplotype1.haplotype.is_combination());
        assert!(!diplotype.haplotype1.haplotype.is_partial());
        let second = diplotype.haplotype2.as_ref().expect("second haplotype");
        assert!(!second.haplotype.is_combination());
        assert!(second.haplotype.is_partial());
    }

    #[test]
    fn partial_reference_double_fixture_matches_java_named_allele_matcher() {
        let (result, warnings) =
            call_cyp2b6_partial_reference_fixture("NamedAlleleMatcher-partialReferenceDouble.vcf");
        let GeneCallKind::Diplotypes(diplotypes) = result.kind else {
            panic!("CYP2B6 result was not a diplotype call: {:?}", result.kind);
        };

        assert_eq!(warnings.values().map(BTreeSet::len).sum::<usize>(), 1);
        assert_eq!(diplotypes.len(), 1);
        let diplotype = diplotypes.first().expect("diplotype");
        assert_eq!(diplotype.name, "g.41006968T>A/g.41006968T>A");
        assert!(!diplotype.haplotype1.haplotype.is_combination());
        assert!(diplotype.haplotype1.haplotype.is_partial());
        let second = diplotype.haplotype2.as_ref().expect("second haplotype");
        assert!(!second.haplotype.is_combination());
        assert!(second.haplotype.is_partial());
    }

    #[test]
    fn unknown_alt_multisample_selected_sample_has_no_warnings_like_java_named_allele_matcher() {
        let definition = read_definition_file(std::path::Path::new(
            "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-ryr1.json",
        ))
        .expect("RYR1 definition");
        let reader = crate::definition::DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition.clone())]
                .into_iter()
                .collect(),
        );
        let records = read_record_summaries(
            "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-unknownAltMultisample.vcf",
            Some("Sample_2"),
        )
        .expect("vcf");
        let mut warnings = records.warnings.clone();
        let allele_map =
            allele_calls_for_locations(&records, reader.locations_of_interest(), &mut warnings)
                .into_iter()
                .map(|call| (format!("{}:{}", call.chromosome, call.position), call))
                .collect::<BTreeMap<_, _>>();

        let result =
            call_ryr1_lowest_function_gene("Sample_2", &definition, &allele_map).expect("RYR1");

        assert_eq!(warnings.values().map(BTreeSet::len).sum::<usize>(), 0);
        assert_eq!(result.gene, "RYR1");
    }

    #[test]
    fn dpyd_diplotype_matcher_fixture_repeats_like_java_named_allele_matcher() {
        let definition = read_dpyd_definition();
        let reader = crate::definition::DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition.clone())]
                .into_iter()
                .collect(),
        );
        let records = read_record_summaries(
            "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-diplotypeMatcher.vcf",
            None,
        )
        .expect("vcf");
        let mut warnings = records.warnings.clone();
        let allele_map =
            allele_calls_for_locations(&records, reader.locations_of_interest(), &mut warnings)
                .into_iter()
                .map(|call| (format!("{}:{}", call.chromosome, call.position), call))
                .collect::<BTreeMap<_, _>>();

        assert_eq!(records.sample_name, "123456");
        assert_eq!(warnings.values().map(BTreeSet::len).sum::<usize>(), 0);
        for _ in 0..10 {
            let result = call_dpyd_lowest_function_gene("123456", &definition, &allele_map)
                .expect("DPYD call");
            assert_eq!(result.gene, "DPYD");
            match result.kind {
                GeneCallKind::Diplotypes(diplotypes) => assert!(!diplotypes.is_empty()),
                GeneCallKind::Haplotypes(haplotypes) => assert!(!haplotypes.is_empty()),
                GeneCallKind::NoCall => panic!("DPYD fixture produced no call"),
            }
        }
    }

    #[test]
    fn call_standard_gene_returns_exact_diplotype_for_java_haplotyper_fixture() {
        let definition = read_definition_file(std::path::Path::new(
            "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/haplotyper.json",
        ))
        .expect("definition");
        let records = read_record_summaries(
            "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/haplotyper.vcf",
            Some("NA12878"),
        )
        .expect("vcf");
        let allele_map = allele_map_from_records(&records);

        let result = call_standard_gene("NA12878", &definition, &allele_map, false, false)
            .expect("standard gene call");

        let GeneCallKind::Diplotypes(diplotypes) = result.kind else {
            panic!("expected diplotype call");
        };
        assert_eq!(result.gene, "CYP3A5");
        assert_eq!(
            diplotypes
                .into_iter()
                .map(|diplotype| diplotype.name)
                .collect::<Vec<_>>(),
            ["*1/*2"]
        );
    }

    #[test]
    fn call_standard_gene_uses_combination_fallback_when_enabled_like_java() {
        let (definition, allele_map) =
            read_combination_fixture("NamedAlleleMatcher-partialWithCombination.vcf");

        let result = call_standard_gene("PharmCAT", &definition, &allele_map, false, true)
            .expect("standard gene call");

        let GeneCallKind::Diplotypes(diplotypes) = result.kind else {
            panic!("expected diplotype call");
        };
        assert_eq!(result.gene, "UGT1A1");
        assert_eq!(
            diplotypes
                .into_iter()
                .map(|diplotype| diplotype.name)
                .collect::<Vec<_>>(),
            [
                "*1/[*6 + *28 + g.233760973C>T]",
                "g.233760973C>T/[*6 + *28]",
                "*6/[*28 + g.233760973C>T]",
                "*28/[*6 + g.233760973C>T]",
            ]
        );
    }

    #[test]
    fn call_standard_gene_returns_no_call_for_partial_missing_without_combinations_like_java() {
        let definition = synthetic_definition();
        let allele_map = allele_map([sample_call(
            "chr1",
            2,
            Some("A"),
            Some("."),
            false,
            false,
            None,
        )]);

        let result = call_standard_gene("Sample_1", &definition, &allele_map, false, false)
            .expect("standard gene call");

        assert!(matches!(result.kind, GeneCallKind::NoCall));
        assert!(result.match_data.has_partial_missing_alleles());
    }

    #[test]
    fn call_standard_gene_returns_no_call_for_cyp2c19_empty_match_set_like_java() {
        let definition = read_cyp2c19_definition();
        let mut allele_map = reference_allele_map(&definition);
        for (rsid, allele1, allele2) in [("rs12769205", "A", "G"), ("rs4244285", "A", "A")] {
            let variant = variant_by_rsid(&definition, rsid);
            allele_map.insert(
                variant.vcf_chr_position(),
                sample_call(
                    &variant.chromosome,
                    variant.position as usize,
                    Some(allele1),
                    Some(allele2),
                    false,
                    allele1 == allele2,
                    None,
                ),
            );
        }

        let result = call_standard_gene("PharmCAT", &definition, &allele_map, false, false)
            .expect("CYP2C19 call");

        assert!(matches!(result.kind, GeneCallKind::NoCall));
        assert!(!result.match_data.has_partial_missing_alleles());
    }

    #[test]
    fn call_standard_gene_returns_no_call_without_sample_alleles_like_java() {
        let definition = synthetic_definition();
        let allele_map = BTreeMap::new();

        let result = call_standard_gene("Sample_1", &definition, &allele_map, false, false)
            .expect("standard gene call");

        assert!(matches!(result.kind, GeneCallKind::NoCall));
        assert!(result.match_data.permutations().is_empty());
    }

    #[test]
    fn call_standard_gene_with_exemption_no_calls_when_required_position_missing_like_java() {
        let definition = synthetic_definition();
        let exemption = synthetic_position_exemption("GENE", [4], []);
        let allele_map = allele_map([sample_call(
            "chr1",
            2,
            Some("A"),
            Some("T"),
            false,
            false,
            None,
        )]);

        let result = call_standard_gene_with_exemption(
            "Sample_1",
            &definition,
            Some(&exemption),
            &allele_map,
            false,
            false,
        )
        .expect("standard gene call");

        assert!(matches!(result.kind, GeneCallKind::NoCall));
        assert_eq!(result.match_data.missing_required_positions, ["chr1:4"]);
        assert!(
            result
                .warnings
                .contains(&GeneCallWarning::MissingRequiredPosition(vec![
                    "chr1:4".to_owned()
                ]))
        );
    }

    #[test]
    fn call_standard_gene_with_exemption_warns_when_amp1_position_missing_like_java() {
        let definition = synthetic_definition();
        let exemption = synthetic_position_exemption("GENE", [], [4]);
        let allele_map = allele_map([sample_call(
            "chr1",
            2,
            Some("A"),
            Some("T"),
            false,
            false,
            None,
        )]);

        let result = call_standard_gene_with_exemption(
            "Sample_1",
            &definition,
            Some(&exemption),
            &allele_map,
            false,
            false,
        )
        .expect("standard gene call");

        assert!(matches!(result.kind, GeneCallKind::Diplotypes(_)));
        assert_eq!(result.match_data.missing_amp1_positions, ["chr1:4"]);
        assert!(
            result
                .warnings
                .contains(&GeneCallWarning::MissingAmp1Position(vec![
                    "chr1:4".to_owned()
                ]))
        );
    }

    #[test]
    fn call_standard_gene_converts_suballeles_like_java_result_builder() {
        let definition = synthetic_suballele_definition();
        let allele_map = allele_map([
            sample_call("chr1", 1, Some("A"), Some("G"), false, false, None),
            sample_call("chr1", 2, Some("C"), Some("C"), false, true, None),
            sample_call("chr1", 3, Some("C"), Some("C"), false, true, None),
        ]);

        let result = call_standard_gene("Sample_1", &definition, &allele_map, false, false)
            .expect("standard gene call");

        let GeneCallKind::Diplotypes(diplotypes) = result.kind else {
            panic!("expected diplotype call");
        };
        assert_eq!(
            diplotypes
                .iter()
                .map(|diplotype| diplotype.name.as_str())
                .collect::<Vec<_>>(),
            ["*1/*4"]
        );
        assert_eq!(diplotypes[0].haplotype2.as_ref().expect("hap2").name, "*4");
    }

    #[test]
    fn finalize_gene_call_result_applies_unphased_diplotype_priority_like_java() {
        let (definition, allele_map) =
            read_combination_fixture("NamedAlleleMatcher-combinationUnphased.vcf");
        let mut result = call_standard_gene("PharmCAT", &definition, &allele_map, false, true)
            .expect("standard gene call");
        let GeneCallKind::Diplotypes(diplotypes) = &result.kind else {
            panic!("expected diplotypes");
        };
        assert!(diplotypes.len() > 1);
        assert!(!result.match_data.effectively_phased);

        let exemption = synthetic_priority_exemption(
            "UGT1A1",
            ["*1/[*6 + *27 + *28 + *80]", "*6/[*27 + *28 + *80]"],
            "*6/[*27 + *28 + *80]",
        );

        result = finalize_gene_call_result(result, &definition, Some(&exemption))
            .expect("finalized result");

        let GeneCallKind::Diplotypes(diplotypes) = result.kind else {
            panic!("expected diplotypes");
        };
        assert_eq!(
            diplotypes
                .into_iter()
                .map(|diplotype| diplotype.name)
                .collect::<Vec<_>>(),
            ["*6/[*27 + *28 + *80]"]
        );
        assert!(result.warnings.contains(&GeneCallWarning::UnphasedPriority));
    }

    #[test]
    fn call_standard_gene_with_exemption_applies_unphased_priority_like_java_result_builder() {
        let (definition, allele_map) =
            read_combination_fixture("NamedAlleleMatcher-combinationUnphased.vcf");
        let exemption = synthetic_priority_exemption(
            "UGT1A1",
            ["*1/[*6 + *27 + *28 + *80]", "*27/[*6 + *28 + *80]"],
            "*27/[*6 + *28 + *80]",
        );

        let result = call_standard_gene_with_exemption(
            "PharmCAT",
            &definition,
            Some(&exemption),
            &allele_map,
            false,
            true,
        )
        .expect("standard gene call");

        let GeneCallKind::Diplotypes(diplotypes) = result.kind else {
            panic!("expected diplotypes");
        };
        assert_eq!(
            diplotypes
                .into_iter()
                .map(|diplotype| diplotype.name)
                .collect::<Vec<_>>(),
            ["*27/[*6 + *28 + *80]"]
        );
        assert!(result.warnings.contains(&GeneCallWarning::UnphasedPriority));
    }

    fn synthetic_definition() -> DefinitionFile {
        let variants = vec![
            variant("chr1", 1),
            variant("chr1", 2),
            variant("chr1", 3),
            variant("chr1", 4),
        ];
        let mut definition = DefinitionFile {
            format_version: "2".to_owned(),
            data_version: None,
            source: None,
            version: None,
            modification_date: None,
            gene_symbol: "GENE".to_owned(),
            orientation: None,
            chromosome: "chr1".to_owned(),
            genome_build: "GRCh38.p13".to_owned(),
            ref_seq_chromosome_id: "NC_000001.11".to_owned(),
            variants,
            named_alleles: vec![
                named_allele("*1", true, [Some("T"), Some("A"), Some("C"), Some("C")]),
                named_allele("*2", false, [None, Some("T"), Some("C"), None]),
                named_allele("*3", false, [None, None, Some("GG"), None]),
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

    fn synthetic_diplotype_definition() -> DefinitionFile {
        let variants = vec![variant("chr1", 1), variant("chr1", 2), variant("chr1", 3)];
        let mut definition = DefinitionFile {
            format_version: "2".to_owned(),
            data_version: None,
            source: None,
            version: None,
            modification_date: None,
            gene_symbol: "CYP2B6".to_owned(),
            orientation: None,
            chromosome: "chr1".to_owned(),
            genome_build: "GRCh38.p13".to_owned(),
            ref_seq_chromosome_id: "NC_000001.11".to_owned(),
            variants,
            named_alleles: vec![
                named_allele3("*1", true, [Some("A"), Some("C"), Some("C")]),
                named_allele3("*4a", false, [Some("G"), None, None]),
                named_allele3("*4b", false, [Some("G"), Some("T"), Some("T")]),
                named_allele3("*17", false, [None, Some("T"), Some("T")]),
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

    fn synthetic_suballele_definition() -> DefinitionFile {
        let mut definition = synthetic_diplotype_definition();
        definition
            .suballeles_map
            .insert("*4a".to_owned(), "*4".to_owned());
        definition
    }

    fn synthetic_wobble_definition() -> DefinitionFile {
        let mut variants = vec![variant("chr1", 1)];
        variants[0].reference = "A".to_owned();
        let mut definition = DefinitionFile {
            format_version: "2".to_owned(),
            data_version: None,
            source: None,
            version: None,
            modification_date: None,
            gene_symbol: "WOBBLE".to_owned(),
            orientation: None,
            chromosome: "chr1".to_owned(),
            genome_build: "GRCh38.p13".to_owned(),
            ref_seq_chromosome_id: "NC_000001.11".to_owned(),
            variants,
            named_alleles: vec![
                named_allele1("*1", true, Some("A")),
                named_allele1("*wobble", false, Some("R")),
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

    fn synthetic_combination_definition() -> DefinitionFile {
        let mut variants = vec![variant("chr1", 1), variant("chr1", 2), variant("chr1", 3)];
        for (index, variant) in variants.iter_mut().enumerate() {
            variant.reference = "C".to_owned();
            variant.alts = vec!["T".to_owned()];
            variant.chromosome_hgvs_name = format!("g.{}C>T", index + 1);
        }
        let mut definition = DefinitionFile {
            format_version: "2".to_owned(),
            data_version: None,
            source: None,
            version: None,
            modification_date: None,
            gene_symbol: "COMBO".to_owned(),
            orientation: None,
            chromosome: "chr1".to_owned(),
            genome_build: "GRCh38.p13".to_owned(),
            ref_seq_chromosome_id: "NC_000001.11".to_owned(),
            variants,
            named_alleles: vec![
                named_allele3("*1", true, [Some("C"), Some("C"), Some("C")]),
                named_allele3("*2", false, [Some("T"), None, None]),
                named_allele3("*5", false, [None, Some("T"), None]),
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

    fn synthetic_lowest_function_definition() -> DefinitionFile {
        let mut definition = synthetic_combination_definition();
        let reference = definition
            .named_alleles
            .iter_mut()
            .find(|allele| allele.reference)
            .expect("reference allele");
        reference.id = "Reference".to_owned();
        reference.name = "Reference".to_owned();
        definition.initialize_derived_fields();
        definition
    }

    fn synthetic_priority_exemption<const N: usize>(
        gene: &str,
        diplotypes: [&str; N],
        pick: &str,
    ) -> DefinitionExemption {
        let list = diplotypes
            .into_iter()
            .map(str::to_owned)
            .collect::<BTreeSet<_>>();
        DefinitionExemption {
            gene: gene.to_owned(),
            required_positions: Default::default(),
            ignored_positions: Default::default(),
            extra_positions: Default::default(),
            ignored_alleles: Default::default(),
            ignored_alleles_lc: Default::default(),
            unphased_diplotype_priorities: [UnphasedDiplotypePriority {
                id: list.iter().cloned().collect::<Vec<_>>().join("|"),
                list,
                pick: pick.to_owned(),
            }]
            .into_iter()
            .collect(),
            amp1_alleles: Vec::new(),
            amp1_positions: Default::default(),
        }
    }

    fn synthetic_position_exemption<const R: usize, const A: usize>(
        gene: &str,
        required_positions: [u64; R],
        amp1_positions: [u64; A],
    ) -> DefinitionExemption {
        DefinitionExemption {
            gene: gene.to_owned(),
            required_positions: required_positions.into_iter().collect(),
            ignored_positions: Default::default(),
            extra_positions: Default::default(),
            ignored_alleles: Default::default(),
            ignored_alleles_lc: Default::default(),
            unphased_diplotype_priorities: Default::default(),
            amp1_alleles: Vec::new(),
            amp1_positions: amp1_positions.into_iter().collect(),
        }
    }

    fn variant(chromosome: &str, position: u64) -> VariantLocus {
        VariantLocus {
            chromosome: chromosome.to_owned(),
            position,
            cpic_position: position,
            rsid: None,
            chromosome_hgvs_name: format!("g.{position}T>A"),
            cpic_alleles: Default::default(),
            cpic_to_vcf_allele_map: Default::default(),
            reference: "T".to_owned(),
            alts: Vec::new(),
        }
    }

    fn named_allele(name: &str, reference: bool, alleles: [Option<&str>; 4]) -> NamedAllele {
        let alleles = alleles
            .into_iter()
            .map(|allele| allele.map(str::to_owned))
            .collect::<Vec<_>>();
        NamedAllele {
            id: name.to_owned(),
            name: name.to_owned(),
            alleles: alleles.clone(),
            cpic_alleles: alleles,
            population_frequency: None,
            reference,
            is_combination_or_partial: false,
            structural_variant: false,
            score: None,
            core_positions: Default::default(),
            missing_positions: Default::default(),
            score_override: None,
            num_combinations: 0,
            num_partials: 0,
        }
    }

    fn named_allele3(name: &str, reference: bool, alleles: [Option<&str>; 3]) -> NamedAllele {
        let alleles = alleles
            .into_iter()
            .map(|allele| allele.map(str::to_owned))
            .collect::<Vec<_>>();
        NamedAllele {
            id: name.to_owned(),
            name: name.to_owned(),
            alleles: alleles.clone(),
            cpic_alleles: alleles,
            population_frequency: None,
            reference,
            is_combination_or_partial: false,
            structural_variant: false,
            score: None,
            core_positions: Default::default(),
            missing_positions: Default::default(),
            score_override: None,
            num_combinations: 0,
            num_partials: 0,
        }
    }

    fn named_allele1(name: &str, reference: bool, allele: Option<&str>) -> NamedAllele {
        let alleles = vec![allele.map(str::to_owned)];
        NamedAllele {
            id: name.to_owned(),
            name: name.to_owned(),
            alleles: alleles.clone(),
            cpic_alleles: alleles,
            population_frequency: None,
            reference,
            is_combination_or_partial: false,
            structural_variant: false,
            score: None,
            core_positions: Default::default(),
            missing_positions: Default::default(),
            score_override: None,
            num_combinations: 0,
            num_partials: 0,
        }
    }

    fn compute_diplotype_names<const N: usize>(
        definition: &DefinitionFile,
        calls: [crate::vcf::SampleAlleleSummary; N],
        top_candidate_only: bool,
    ) -> Vec<String> {
        let allele_map = calls
            .into_iter()
            .map(|call| (format!("{}:{}", call.chromosome, call.position), call))
            .collect::<BTreeMap<_, _>>();
        let mut data = MatchData::new("Sample_1", "CYP2B6", definition, &allele_map);
        data.marshall_haplotypes(definition);
        data.generate_sample_permutations().expect("permutations");
        data.compute_diplotypes(top_candidate_only)
            .into_iter()
            .map(|diplotype| diplotype.name)
            .collect()
    }

    fn compute_diplotype_names_and_scores<const N: usize>(
        definition: &DefinitionFile,
        calls: [crate::vcf::SampleAlleleSummary; N],
        top_candidate_only: bool,
    ) -> Vec<(String, i32)> {
        let allele_map = calls
            .into_iter()
            .map(|call| (format!("{}:{}", call.chromosome, call.position), call))
            .collect::<BTreeMap<_, _>>();
        let mut data = MatchData::new("Sample_1", &definition.gene_symbol, definition, &allele_map);
        data.marshall_haplotypes(definition);
        data.generate_sample_permutations().expect("permutations");
        data.compute_diplotypes(top_candidate_only)
            .into_iter()
            .map(|diplotype| (diplotype.name, diplotype.score))
            .collect()
    }

    fn compute_combination_diplotype_names<const N: usize>(
        definition: &DefinitionFile,
        calls: [crate::vcf::SampleAlleleSummary; N],
        find_partials: bool,
    ) -> Vec<String> {
        let allele_map = calls
            .into_iter()
            .map(|call| (format!("{}:{}", call.chromosome, call.position), call))
            .collect::<BTreeMap<_, _>>();
        let mut data = MatchData::new("Sample_1", &definition.gene_symbol, definition, &allele_map);
        data.marshall_haplotypes(definition);
        data.generate_sample_permutations().expect("permutations");
        data.compute_combination_diplotypes(definition, find_partials)
            .into_iter()
            .map(|diplotype| diplotype.name)
            .collect()
    }

    fn compute_combination_fixture_names(vcf_fixture: &str) -> Vec<String> {
        let (definition, allele_map) = read_combination_fixture(vcf_fixture);

        let mut data = MatchData::new("PharmCAT", "UGT1A1", &definition, &allele_map);
        data.marshall_haplotypes(&definition);
        data.generate_sample_permutations().expect("permutations");
        data.compute_combination_diplotypes(&definition, true)
            .into_iter()
            .map(|diplotype| diplotype.name)
            .collect()
    }

    fn read_combination_fixture(
        vcf_fixture: &str,
    ) -> (
        DefinitionFile,
        BTreeMap<String, crate::vcf::SampleAlleleSummary>,
    ) {
        let definition = read_definition_file(std::path::Path::new(
            "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-combination.json",
        ))
        .expect("definition");
        let records = read_record_summaries(
            format!(
                "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/{vcf_fixture}"
            ),
            Some("PharmCAT"),
        )
        .expect("vcf");
        (definition, allele_map_from_records(&records))
    }

    fn call_cyp2b6_partial_reference_fixture(vcf_fixture: &str) -> (GeneCallResult, VcfWarnings) {
        let definition = read_definition_file(std::path::Path::new(
            "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-cyp2b6.json",
        ))
        .expect("CYP2B6 definition");
        let reader = crate::definition::DefinitionReader::from_definitions(
            [(definition.gene_symbol.clone(), definition.clone())]
                .into_iter()
                .collect(),
        );
        let records = read_record_summaries(
            format!(
                "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/{vcf_fixture}"
            ),
            Some("PharmCAT"),
        )
        .expect("vcf");
        let mut warnings = records.warnings.clone();
        let allele_map =
            allele_calls_for_locations(&records, reader.locations_of_interest(), &mut warnings)
                .into_iter()
                .map(|call| (format!("{}:{}", call.chromosome, call.position), call))
                .collect::<BTreeMap<_, _>>();
        let result = call_standard_gene("PharmCAT", &definition, &allele_map, true, true)
            .expect("CYP2B6 call");
        (result, warnings)
    }

    fn allele_map_from_records(
        records: &VcfRecords,
    ) -> BTreeMap<String, crate::vcf::SampleAlleleSummary> {
        records
            .records
            .iter()
            .filter_map(|record| {
                record
                    .allele_call
                    .clone()
                    .map(|call| (format!("{}:{}", call.chromosome, call.position), call))
            })
            .collect()
    }

    fn read_dpyd_definition() -> DefinitionFile {
        read_definition_file(std::path::Path::new(
            "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-dpyd.json",
        ))
        .expect("DPYD definition")
    }

    fn read_cyp2c19_definition() -> DefinitionFile {
        read_definition_file(std::path::Path::new(
            "../../repos/PharmCAT/src/main/resources/org/pharmgkb/pharmcat/definition/alleles/CYP2C19_translation.json",
        ))
        .expect("CYP2C19 definition")
    }

    fn read_cyp2b6_definition() -> DefinitionFile {
        read_definition_file(std::path::Path::new(
            "../../repos/PharmCAT/src/main/resources/org/pharmgkb/pharmcat/definition/alleles/CYP2B6_translation.json",
        ))
        .expect("CYP2B6 definition")
    }

    fn read_nat2_definition() -> DefinitionFile {
        read_definition_file(std::path::Path::new(
            "../../repos/PharmCAT/src/main/resources/org/pharmgkb/pharmcat/definition/alleles/NAT2_translation.json",
        ))
        .expect("NAT2 definition")
    }

    fn read_slco1b1_definition() -> DefinitionFile {
        read_definition_file(std::path::Path::new(
            "../../repos/PharmCAT/src/main/resources/org/pharmgkb/pharmcat/definition/alleles/SLCO1B1_translation.json",
        ))
        .expect("SLCO1B1 definition")
    }

    fn read_default_exemptions() -> BTreeMap<String, DefinitionExemption> {
        read_exemptions_file(std::path::Path::new(
            "../../repos/PharmCAT/src/main/resources/org/pharmgkb/pharmcat/definition/alleles/exemptions.json",
        ))
        .expect("default exemptions")
    }

    fn read_ryr1_definition() -> DefinitionFile {
        read_definition_file(std::path::Path::new(
            "../../repos/PharmCAT/src/main/resources/org/pharmgkb/pharmcat/definition/alleles/RYR1_translation.json",
        ))
        .expect("RYR1 definition")
    }

    fn reference_allele_map(
        definition: &DefinitionFile,
    ) -> BTreeMap<String, crate::vcf::SampleAlleleSummary> {
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
                    None,
                );
                (variant.vcf_chr_position(), call)
            })
            .collect()
    }

    fn variant_by_rsid<'a>(definition: &'a DefinitionFile, rsid: &str) -> &'a VariantLocus {
        definition
            .variants
            .iter()
            .find(|variant| variant.rsid.as_deref() == Some(rsid))
            .unwrap_or_else(|| panic!("variant {rsid}"))
    }

    fn cyp2c19_permutation_allele_map(
        definition: &DefinitionFile,
        phased_overrides: [bool; 4],
    ) -> BTreeMap<String, crate::vcf::SampleAlleleSummary> {
        let mut allele_map = reference_allele_map(definition);
        for ((rsid, allele1, allele2), phased) in [
            ("rs12248560", "C", "T"),
            ("rs12769205", "G", "A"),
            ("rs4244285", "A", "G"),
            ("rs3758581", "G", "G"),
        ]
        .into_iter()
        .zip(phased_overrides)
        {
            let variant = variant_by_rsid(definition, rsid);
            allele_map.insert(
                variant.vcf_chr_position(),
                sample_call(
                    &variant.chromosome,
                    variant.position as usize,
                    Some(allele1),
                    Some(allele2),
                    phased,
                    phased || allele1 == allele2,
                    None,
                ),
            );
        }
        allele_map
    }

    fn nat2_unphased_priority_allele_map<const N: usize>(
        definition: &DefinitionFile,
        variants: [(&str, &str, &str); N],
    ) -> BTreeMap<String, crate::vcf::SampleAlleleSummary> {
        let mut allele_map = reference_allele_map(definition);
        for (rsid, allele1, allele2) in variants {
            let variant = variant_by_rsid(definition, rsid);
            allele_map.insert(
                variant.vcf_chr_position(),
                sample_call(
                    &variant.chromosome,
                    variant.position as usize,
                    Some(allele1),
                    Some(allele2),
                    false,
                    false,
                    None,
                ),
            );
        }
        allele_map
    }

    fn nat2_phased_combination_allele_map(
        definition: &DefinitionFile,
    ) -> BTreeMap<String, crate::vcf::SampleAlleleSummary> {
        let mut allele_map = reference_allele_map(definition);
        for (rsid, allele1, allele2) in [
            ("rs1801279", "A", "G"),
            ("rs12720065", "G", "C"),
            ("rs1799930", "A", "G"),
            ("rs1208", "A", "G"),
        ] {
            let variant = variant_by_rsid(definition, rsid);
            allele_map.insert(
                variant.vcf_chr_position(),
                sample_call(
                    &variant.chromosome,
                    variant.position as usize,
                    Some(allele1),
                    Some(allele2),
                    true,
                    true,
                    None,
                ),
            );
        }
        allele_map
    }

    fn cyp2b6_wobble_scoring_allele_map(
        definition: &DefinitionFile,
    ) -> BTreeMap<String, crate::vcf::SampleAlleleSummary> {
        let mut allele_map = reference_allele_map(definition);
        for (rsid, allele1, allele2) in [
            ("rs3745274", "G", "T"),
            ("rs2279343", "A", "G"),
            ("rs28399499", "T", "C"),
        ] {
            let variant = variant_by_rsid(definition, rsid);
            allele_map.insert(
                variant.vcf_chr_position(),
                sample_call(
                    &variant.chromosome,
                    variant.position as usize,
                    Some(allele1),
                    Some(allele2),
                    false,
                    false,
                    None,
                ),
            );
        }
        allele_map
    }

    fn allele_map<const N: usize>(
        calls: [crate::vcf::SampleAlleleSummary; N],
    ) -> BTreeMap<String, crate::vcf::SampleAlleleSummary> {
        calls
            .into_iter()
            .map(|call| (format!("{}:{}", call.chromosome, call.position), call))
            .collect()
    }

    fn dpyd_match_data(
        definition: &DefinitionFile,
        allele_map: &BTreeMap<String, crate::vcf::SampleAlleleSummary>,
    ) -> MatchData {
        let mut data = MatchData::new("Sample_1", "DPYD", definition, allele_map);
        data.marshall_haplotypes(definition);
        data.generate_sample_permutations().expect("permutations");
        data
    }

    fn sample_call(
        chromosome: &str,
        position: usize,
        allele1: Option<&str>,
        allele2: Option<&str>,
        phased: bool,
        effectively_phased: bool,
        phase_set: Option<i32>,
    ) -> crate::vcf::SampleAlleleSummary {
        crate::vcf::SampleAlleleSummary {
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
            phase_set,
            undocumented_variations: BTreeSet::new(),
            treat_undocumented_variations_as_reference: false,
        }
    }
}
