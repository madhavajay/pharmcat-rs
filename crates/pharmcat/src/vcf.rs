//! VCF helpers for the PharmCAT Rust port.
//!
//! This starts by mirroring the Java `VcfSampleReader`: read VCF metadata,
//! collect sample IDs, and validate that contig assembly metadata is
//! consistent. Row-level PharmCAT genotype extraction will build on this.

use std::{
    collections::{BTreeMap, BTreeSet},
    fs::File,
    io::{self, BufRead, BufReader, Cursor, Read},
    path::Path,
};

use flate2::read::MultiGzDecoder;
use noodles_vcf::{
    self,
    header::record::value::{
        Map,
        map::{Format, format::Number},
    },
    variant::record::{AlternateBases as _, Filters as _, Ids as _, Samples as _},
};

use crate::definition::VariantLocus;

/// Java `VcfReader.MSG_AD_FORMAT_MISSING`.
pub const MSG_AD_FORMAT_MISSING: &str = "AD format is not defined.  Assuming AD field is valid.";

/// Metadata PharmCAT needs before reading variant rows.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct VcfHeaderSummary {
    /// Sample names in VCF header order.
    pub samples: Vec<String>,
    /// Shared contig assembly, if any contig declared one.
    pub genome_build: Option<String>,
}

/// Raw VCF records for one selected sample.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct VcfRecords {
    /// VCF header metadata.
    pub header: VcfHeaderSummary,
    /// Selected sample name.
    pub sample_name: String,
    /// FORMAT/AD handling selected from the VCF header.
    pub allelic_depth_policy: AllelicDepthPolicy,
    /// Warnings keyed by chromosome position or `VCF`, matching Java's grouping.
    pub warnings: VcfWarnings,
    /// Variant rows in source order.
    pub records: Vec<VcfRecordSummary>,
}

/// Sorted warnings keyed by chromosome position or the VCF header.
pub type VcfWarnings = BTreeMap<String, BTreeSet<String>>;

/// Whether PharmCAT should trust row-level `AD` sample values.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum AllelicDepthPolicy {
    /// `FORMAT/AD` was declared with `Number=R`.
    UseDefinedReferenceAlternate,
    /// `FORMAT/AD` was absent, but Java still trusts row-level `AD` values and warns.
    UseMissingDefinition,
    /// `FORMAT/AD` had `Number=.` and Java ignores row-level `AD`.
    IgnoreUnknownNumber,
    /// `FORMAT/AD` had an unexpected `Number` value and Java ignores row-level `AD`.
    IgnoreInvalidNumber,
}

impl AllelicDepthPolicy {
    fn should_use(self) -> bool {
        matches!(
            self,
            Self::UseDefinedReferenceAlternate | Self::UseMissingDefinition
        )
    }
}

/// PharmCAT-relevant VCF row data before allele interpretation.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct VcfRecordSummary {
    /// Chromosome or contig name.
    pub chromosome: String,
    /// 1-based variant position.
    pub position: usize,
    /// Record IDs.
    pub ids: Vec<String>,
    /// Reference allele bases.
    pub reference: String,
    /// Alternate allele values in VCF order.
    pub alternates: Vec<String>,
    /// Filter values.
    pub filters: Vec<String>,
    /// FORMAT keys for the sample columns.
    pub format_keys: Vec<String>,
    /// Selected sample values.
    pub sample: VcfSampleFields,
    /// Whether Java would skip this row because a prior valid row already exists at this position.
    pub skipped_duplicate: bool,
    /// Whether Java would discard this row because PharmCAT preprocessor REF mismatch filter is set.
    pub discarded_by_preprocessor_filter: bool,
    /// Whether Java would discard this row because `GT` and `AD` disagree.
    pub discarded_by_allelic_depth: bool,
    /// Whether Java would discard this row because a selected REF/ALT allele is invalid.
    pub discarded_by_allele_validation: bool,
    /// Java `SampleAllele`-style call data, when the row is not discarded before allele extraction.
    pub allele_call: Option<SampleAlleleSummary>,
}

/// Java `SampleAllele`-style VCF call data for one row.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SampleAlleleSummary {
    /// Chromosome or contig name.
    pub chromosome: String,
    /// 1-based variant position.
    pub position: usize,
    /// First allele selected by the first GT token.
    pub allele1: Option<String>,
    /// Second allele selected by the second GT token.
    pub allele2: Option<String>,
    /// REF followed by ALT alleles from the VCF.
    pub vcf_alleles: Vec<String>,
    /// Raw genotype string.
    pub genotype: String,
    /// Allele call rendered like Java `SampleAllele.getVcfCall()`.
    pub vcf_call: String,
    /// Whether the VCF genotype used `|`.
    pub phased: bool,
    /// Whether PharmCAT treats this row as phased.
    pub effectively_phased: bool,
    /// Raw `PS` as an integer when Java would keep it.
    pub phase_set: Option<i32>,
    /// ALT alleles observed in the VCF but not documented by the PharmCAT definition.
    pub undocumented_variations: BTreeSet<String>,
    /// Whether Java replaced undocumented variation calls with reference for matching.
    pub treat_undocumented_variations_as_reference: bool,
}

/// Raw selected sample fields split out for PharmCAT genotype logic.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct VcfSampleFields {
    /// Sample name.
    pub sample_name: String,
    /// Raw sample column value.
    pub raw: String,
    /// Raw `GT` value, when present.
    pub genotype: Option<String>,
    /// Raw `AD` value, when present.
    pub allelic_depth: Option<String>,
    /// Raw `PS` value, when present.
    pub phase_set: Option<String>,
}

/// Reads VCF metadata and sample names from `src`.
pub fn read_header_summary<P>(src: P) -> Result<VcfHeaderSummary, ReadHeaderError>
where
    P: AsRef<Path>,
{
    let src = src.as_ref();

    if !is_vcf_file(src) {
        return Err(ReadHeaderError::NotVcf(src.display().to_string()));
    }

    let input = open_vcf_input(src)?;
    let mut reader = noodles_vcf::io::Reader::new(input.reader);
    let header = reader.read_header()?;

    let samples = header.sample_names().iter().cloned().collect();
    let genome_build = shared_contig_assembly(&header)?;

    Ok(VcfHeaderSummary {
        samples,
        genome_build,
    })
}

/// Reads VCF records for `sample_name`, or the first sample when none is given.
pub fn read_record_summaries<P>(
    src: P,
    sample_name: Option<&str>,
) -> Result<VcfRecords, ReadVcfError>
where
    P: AsRef<Path>,
{
    let src = src.as_ref();

    if !is_vcf_file(src) {
        return Err(ReadHeaderError::NotVcf(src.display().to_string()).into());
    }

    let input = open_vcf_input(src)?;
    let invalid_ad_number = input.invalid_ad_number.clone();
    let mut reader = noodles_vcf::io::Reader::new(input.reader);
    let header = reader.read_header()?;

    let samples: Vec<String> = header.sample_names().iter().cloned().collect();
    let selected_sample_name = match sample_name {
        Some(name) => name.to_owned(),
        None => samples
            .first()
            .cloned()
            .ok_or(ReadVcfError::NoSamplesDeclared)?,
    };
    let selected_sample_index = samples
        .iter()
        .position(|sample| sample == &selected_sample_name)
        .ok_or_else(|| ReadVcfError::SampleNotFound(selected_sample_name.clone()))?;

    let header_summary = VcfHeaderSummary {
        samples,
        genome_build: shared_contig_assembly(&header)?,
    };
    let allelic_depth_policy = if invalid_ad_number.is_some() {
        AllelicDepthPolicy::IgnoreInvalidNumber
    } else {
        allelic_depth_policy(&header)
    };
    let mut warnings = VcfWarnings::new();

    if matches!(
        allelic_depth_policy,
        AllelicDepthPolicy::IgnoreInvalidNumber
    ) {
        let number = invalid_ad_number
            .or_else(|| header.formats().get("AD").map(format_number_string))
            .unwrap_or_default();
        add_warning(
            &mut warnings,
            "VCF",
            format!(
                "INFO header for AD has unexpected number ({number}). Expecting 'R'. Treating number as '.' and ignoring AD field."
            ),
        );
    }

    let mut records = Vec::new();
    let mut valid_positions = BTreeSet::new();
    let mut discarded_positions = BTreeSet::new();

    for result in reader.records() {
        let record = result?;
        let chromosome = record.reference_sequence_name().to_owned();
        let position = record
            .variant_start()
            .transpose()?
            .ok_or(ReadVcfError::MissingPosition)?
            .get();
        let ids = record.ids().iter().map(str::to_owned).collect();
        let reference = record.reference_bases().to_owned();
        let alternates = record
            .alternate_bases()
            .iter()
            .map(|result| result.map(str::to_owned))
            .collect::<io::Result<Vec<_>>>()?;
        let filters = record
            .filters()
            .iter(&header)
            .map(|result| result.map(str::to_owned))
            .collect::<io::Result<Vec<_>>>()?;
        let samples = record.samples();
        let format_keys = samples
            .column_names(&header)
            .map(|result| result.map(str::to_owned))
            .collect::<io::Result<Vec<_>>>()?;
        let raw_sample = samples
            .get_index(selected_sample_index)
            .ok_or(ReadVcfError::MissingSampleColumn)?
            .as_ref()
            .to_owned();
        let sample = summarize_sample(&selected_sample_name, &format_keys, &raw_sample);
        let chr_position = format!("{chromosome}:{position}");

        let skipped_duplicate = valid_positions.contains(&chr_position);
        let discarded_by_preprocessor_filter = !skipped_duplicate
            && apply_preprocessor_filter_warnings(
                &filters,
                &reference,
                &chr_position,
                &mut warnings,
            );
        let discarded_by_allelic_depth = !skipped_duplicate
            && !discarded_by_preprocessor_filter
            && check_allelic_depth(&sample, allelic_depth_policy, &chr_position, &mut warnings);
        let (discarded_by_allele_validation, allele_call) = if skipped_duplicate {
            if !has_pharmcat_preprocessor_filter(&filters) {
                add_warning(
                    &mut warnings,
                    &chr_position,
                    "Duplicate entry found in VCF; first valid entry trumps others.",
                );
            }
            (false, None)
        } else if discarded_by_preprocessor_filter || discarded_by_allelic_depth {
            (false, None)
        } else {
            sample_allele_summary(
                &chromosome,
                position,
                &reference,
                &alternates,
                &sample,
                &chr_position,
                &mut warnings,
            )?
        };

        if allele_call.is_some() {
            valid_positions.insert(chr_position.clone());
            if discarded_positions.remove(&chr_position) {
                add_warning(
                    &mut warnings,
                    &chr_position,
                    "Duplicate entry found in VCF; this entry trumps previous invalid entry.",
                );
            }
        } else if !skipped_duplicate {
            discarded_positions.insert(chr_position);
        }

        records.push(VcfRecordSummary {
            chromosome,
            position,
            ids,
            reference,
            alternates,
            filters,
            format_keys,
            sample,
            skipped_duplicate,
            discarded_by_preprocessor_filter,
            discarded_by_allelic_depth,
            discarded_by_allele_validation,
            allele_call,
        });
    }

    Ok(VcfRecords {
        header: header_summary,
        sample_name: selected_sample_name,
        allelic_depth_policy,
        warnings,
        records,
    })
}

/// Returns allele calls for PharmCAT definition locations, adding Java-style definition-aware warnings.
pub fn allele_calls_for_locations(
    records: &VcfRecords,
    locations_of_interest: &BTreeMap<String, VariantLocus>,
    warnings: &mut VcfWarnings,
) -> Vec<SampleAlleleSummary> {
    allele_calls_for_locations_with_genes(
        records,
        locations_of_interest,
        &BTreeMap::new(),
        false,
        warnings,
    )
}

/// Returns allele calls for PharmCAT definition locations with Java undocumented-as-reference rules.
pub fn allele_calls_for_locations_with_genes(
    records: &VcfRecords,
    locations_of_interest: &BTreeMap<String, VariantLocus>,
    locations_by_gene: &BTreeMap<String, String>,
    find_combinations: bool,
    warnings: &mut VcfWarnings,
) -> Vec<SampleAlleleSummary> {
    records
        .records
        .iter()
        .filter_map(|record| {
            let chr_position = format!("{}:{}", record.chromosome, record.position);
            let variant = locations_of_interest.get(&chr_position)?;
            if record.skipped_duplicate
                || record.discarded_by_preprocessor_filter
                || record.discarded_by_allelic_depth
                || record.discarded_by_allele_validation
            {
                return None;
            }
            if record.reference != variant.reference {
                add_warning(
                    warnings,
                    chr_position,
                    format!(
                        "Discarded genotype at this position because REF in VCF ({}) does not match expected reference ({})",
                        record.reference, variant.reference
                    ),
                );
                return None;
            }
            let mut allele_call = record.allele_call.clone()?;
            let undocumented_variations = selected_undocumented_variations(&allele_call, variant);
            if !undocumented_variations.is_empty() {
                let treat_as_reference = locations_by_gene.get(&chr_position).is_some_and(|gene| {
                    treat_undocumented_variations_as_reference(gene, find_combinations)
                });
                let mut warning = format!(
                    "The genetic variation at this position does not match what is in the allele definition (expected {}, found {} in VCF)",
                    variant.alts.join("/"),
                    undocumented_variations
                        .iter()
                        .cloned()
                        .collect::<Vec<_>>()
                        .join("/")
                );
                if treat_as_reference {
                    warning.push_str(".  Undocumented variations will be replaced with reference.");
                }
                add_warning(
                    warnings,
                    chr_position,
                    warning,
                );
                allele_call.treat_undocumented_variations_as_reference = treat_as_reference;
                allele_call.undocumented_variations = undocumented_variations;
            }
            Some(allele_call)
        })
        .collect()
}

fn treat_undocumented_variations_as_reference(gene: &str, find_combinations: bool) -> bool {
    matches!(
        gene,
        "CACNA1S" | "DPYD" | "G6PD" | "NUDT15" | "RYR1" | "TPMT"
    ) && (!find_combinations || matches!(gene, "DPYD" | "RYR1"))
}

fn selected_undocumented_variations(
    allele_call: &SampleAlleleSummary,
    variant: &VariantLocus,
) -> BTreeSet<String> {
    [&allele_call.allele1, &allele_call.allele2]
        .into_iter()
        .filter_map(|allele| allele.as_deref())
        .filter(|allele| *allele != variant.reference)
        .filter(|allele| !variant.alts.iter().any(|expected| expected == allele))
        .map(str::to_owned)
        .collect()
}

/// Returns whether `path` has a PharmCAT-supported VCF extension.
pub fn is_vcf_file(path: &Path) -> bool {
    let Some(filename) = path.to_str() else {
        return false;
    };

    filename.ends_with(".vcf") || filename.ends_with(".vcf.gz") || filename.ends_with(".vcf.bgz")
}

fn open_vcf(path: &Path) -> io::Result<Box<dyn BufRead>> {
    let file = File::open(path)?;

    if is_gzipped_vcf_file(path) {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

struct VcfInput {
    reader: Box<dyn BufRead>,
    invalid_ad_number: Option<String>,
}

fn open_vcf_input(path: &Path) -> io::Result<VcfInput> {
    let mut reader = open_vcf(path)?;
    let mut src = String::new();
    reader.read_to_string(&mut src)?;

    let (src, invalid_ad_number) = sanitize_invalid_ad_number(src);

    Ok(VcfInput {
        reader: Box::new(BufReader::new(Cursor::new(src))),
        invalid_ad_number,
    })
}

fn sanitize_invalid_ad_number(src: String) -> (String, Option<String>) {
    let mut invalid_ad_number = None;
    let mut out = String::with_capacity(src.len());

    for line in src.lines() {
        if line.starts_with("##FORMAT=<ID=AD,") {
            let (line, number) = sanitize_ad_format_line(line);
            invalid_ad_number = invalid_ad_number.or(number);
            out.push_str(&line);
        } else if line.starts_with("##FORMAT=<ID=FT,") {
            out.push_str(&sanitize_ft_format_line(line));
        } else {
            out.push_str(line);
        }
        out.push('\n');
    }

    (out, invalid_ad_number)
}

fn sanitize_ad_format_line(line: &str) -> (String, Option<String>) {
    let Some(number_start) = line.find("Number=") else {
        return (line.to_owned(), None);
    };
    let value_start = number_start + "Number=".len();
    let Some(value_end_offset) = line[value_start..].find(',') else {
        return (line.to_owned(), None);
    };
    let value_end = value_start + value_end_offset;
    let number = &line[value_start..value_end];

    if matches!(
        number,
        "A" | "R" | "G" | "LA" | "LR" | "LG" | "P" | "M" | "."
    ) || number.parse::<usize>().is_ok()
    {
        return (line.to_owned(), None);
    }

    let mut sanitized = String::with_capacity(line.len());
    sanitized.push_str(&line[..value_start]);
    sanitized.push('.');
    sanitized.push_str(&line[value_end..]);

    (sanitized, Some(number.to_owned()))
}

fn sanitize_ft_format_line(line: &str) -> String {
    line.replace("Number=.", "Number=1")
}

fn is_gzipped_vcf_file(path: &Path) -> bool {
    let Some(filename) = path.to_str() else {
        return false;
    };

    filename.ends_with(".vcf.gz") || filename.ends_with(".vcf.bgz")
}

fn shared_contig_assembly(header: &noodles_vcf::Header) -> Result<Option<String>, ReadHeaderError> {
    let assemblies: BTreeSet<_> = header
        .contigs()
        .values()
        .filter_map(|contig| {
            contig
                .other_fields()
                .iter()
                .find(|(key, _)| key.as_ref() == "assembly")
                .map(|(_, value)| value.clone())
        })
        .collect();

    match assemblies.len() {
        0 => Ok(None),
        1 => Ok(assemblies.into_iter().next()),
        _ => Err(ReadHeaderError::MixedAssemblies(
            assemblies.into_iter().collect(),
        )),
    }
}

fn allelic_depth_policy(header: &noodles_vcf::Header) -> AllelicDepthPolicy {
    match header.formats().get("AD") {
        None => AllelicDepthPolicy::UseMissingDefinition,
        Some(format) => match format.number() {
            Number::ReferenceAlternateBases => AllelicDepthPolicy::UseDefinedReferenceAlternate,
            Number::Unknown => AllelicDepthPolicy::IgnoreUnknownNumber,
            _ => AllelicDepthPolicy::IgnoreInvalidNumber,
        },
    }
}

fn format_number_string(format: &Map<Format>) -> String {
    match format.number() {
        Number::Count(n) => n.to_string(),
        Number::AlternateBases => "A".to_owned(),
        Number::ReferenceAlternateBases => "R".to_owned(),
        Number::Samples => "G".to_owned(),
        Number::LocalAlternateBases => "LA".to_owned(),
        Number::LocalReferenceAlternateBases => "LR".to_owned(),
        Number::LocalSamples => "LG".to_owned(),
        Number::Ploidy => "P".to_owned(),
        Number::BaseModifications => "M".to_owned(),
        Number::Unknown => ".".to_owned(),
    }
}

fn summarize_sample(sample_name: &str, format_keys: &[String], raw: &str) -> VcfSampleFields {
    let values: Vec<&str> = raw.split(':').collect();

    VcfSampleFields {
        sample_name: sample_name.to_owned(),
        raw: raw.to_owned(),
        genotype: sample_field(format_keys, &values, "GT"),
        allelic_depth: sample_field(format_keys, &values, "AD"),
        phase_set: sample_field(format_keys, &values, "PS"),
    }
}

fn sample_field(format_keys: &[String], values: &[&str], key: &str) -> Option<String> {
    format_keys
        .iter()
        .position(|format_key| format_key == key)
        .and_then(|index| values.get(index))
        .map(|value| (*value).to_owned())
}

fn apply_preprocessor_filter_warnings(
    filters: &[String],
    reference: &str,
    chr_position: &str,
    warnings: &mut VcfWarnings,
) -> bool {
    if filters.iter().any(|filter| filter == "PCATxREF") {
        add_warning(
            warnings,
            chr_position,
            format!(
                "Discarded genotype at this position because REF in VCF ({reference}) does not match expected reference"
            ),
        );
        return true;
    }

    if filters.iter().any(|filter| filter == "PCATxALT") {
        add_warning(
            warnings,
            chr_position,
            "The genetic variation at this position does not match what is in the allele definition",
        );
    }

    if filters.iter().any(|filter| filter == "PCATxINDEL") {
        add_warning(
            warnings,
            chr_position,
            "Genotype at this position uses unexpected format for INDEL",
        );
    }

    false
}

fn has_pharmcat_preprocessor_filter(filters: &[String]) -> bool {
    filters
        .iter()
        .any(|filter| matches!(filter.as_str(), "PCATxREF" | "PCATxALT" | "PCATxINDEL"))
}

fn check_allelic_depth(
    sample: &VcfSampleFields,
    policy: AllelicDepthPolicy,
    chr_position: &str,
    warnings: &mut VcfWarnings,
) -> bool {
    if !policy.should_use() {
        return false;
    }

    let Some(allelic_depth) = sample.allelic_depth.as_deref() else {
        return false;
    };

    if matches!(policy, AllelicDepthPolicy::UseMissingDefinition) {
        add_warning(warnings, "VCF", MSG_AD_FORMAT_MISSING);
    }

    if allelic_depth == "." {
        return false;
    }

    let Some(genotype) = sample.genotype.as_deref() else {
        return false;
    };

    let depths = match parse_allelic_depths(allelic_depth) {
        Ok(depths) => depths,
        Err(()) => {
            add_warning(
                warnings,
                chr_position,
                format!("Invalid allelic depth (AD) field: {allelic_depth}"),
            );
            return false;
        }
    };

    let genotype_alleles = parse_called_genotype_alleles(genotype);

    if genotype_alleles.len() > 1 && depths.iter().filter(|depth| **depth > 0).count() == 1 {
        add_warning(
            warnings,
            chr_position,
            format!(
                "Discarding genotype at this position because GT field indicates heterozygous ({genotype}) but AD field indicates homozygous ({allelic_depth})"
            ),
        );
        true
    } else {
        false
    }
}

fn sample_allele_summary(
    chromosome: &str,
    position: usize,
    reference: &str,
    alternates: &[String],
    sample: &VcfSampleFields,
    chr_position: &str,
    warnings: &mut VcfWarnings,
) -> Result<(bool, Option<SampleAlleleSummary>), ReadVcfError> {
    let Some(genotype) = sample.genotype.as_deref() else {
        add_warning(warnings, chr_position, "Ignoring: no genotype");
        return Ok((false, None));
    };

    let genotype_tokens: Vec<&str> = genotype.split(['|', '/']).collect();
    let called_genotypes = parse_called_genotype_tokens(&genotype_tokens)?;

    if called_genotypes.is_empty() {
        add_warning(
            warnings,
            chr_position,
            format!("Ignoring: no call ({genotype})"),
        );
        return Ok((false, None));
    }

    if is_haploid_chromosome(chromosome) {
        if called_genotypes.len() > 1 {
            add_warning(
                warnings,
                chr_position,
                format!(
                    "{} genotypes found (GT={genotype}) for haploid chromosome. Will only use first non-missing genotype.",
                    called_genotypes.len()
                ),
            );
        }
    } else if called_genotypes.len() > 2 {
        add_warning(
            warnings,
            chr_position,
            format!(
                "{} genotypes found (GT={genotype}). Will only use first two genotypes.",
                called_genotypes.len()
            ),
        );
    } else if called_genotypes.len() == 1 && chromosome != "chrX" {
        if called_genotypes[0] == 0 {
            add_warning(
                warnings,
                chr_position,
                format!(
                    "Ignoring: only a single genotype found (GT={genotype}).  Since it's reference, treating this as a missing position."
                ),
            );
            return Ok((false, None));
        }

        add_warning(
            warnings,
            chr_position,
            format!("1 genotype found (GT={genotype}), expecting 2."),
        );
    }

    let discarded_by_allele_validation = !validate_alleles(
        chr_position,
        reference,
        alternates,
        &called_genotypes,
        warnings,
    );
    if discarded_by_allele_validation {
        return Ok((true, None));
    }

    for allele_idx in &called_genotypes {
        if *allele_idx > alternates.len() {
            return Err(ReadVcfError::InvalidGenotypeAllele {
                allele: *allele_idx,
                chr_position: chr_position.to_owned(),
                alt_count: alternates.len(),
            });
        }
    }

    let mut vcf_alleles = Vec::with_capacity(alternates.len() + 1);
    vcf_alleles.push(reference.to_uppercase());
    vcf_alleles.extend(alternates.iter().map(|allele| allele.to_uppercase()));

    let allele1 = genotype_tokens
        .first()
        .and_then(|token| genotype_token_to_allele(token, &vcf_alleles));
    let allele2 = if !is_haploid_chromosome(chromosome) || allele1.is_none() {
        genotype_tokens
            .get(1)
            .and_then(|token| genotype_token_to_allele(token, &vcf_alleles))
    } else {
        None
    };
    let phased = !genotype.contains('/');
    let called_unique: BTreeSet<_> = called_genotypes.iter().copied().collect();
    let mut effectively_phased = true;
    let mut phase_set = None;

    if genotype.contains('/') {
        if called_unique.len() > 1 {
            effectively_phased = false;
        }
    } else if let Some(raw_phase_set) = sample.phase_set.as_deref()
        && raw_phase_set != "."
    {
        phase_set = Some(
            raw_phase_set
                .parse()
                .map_err(|_| ReadVcfError::InvalidPhaseSet {
                    phase_set: raw_phase_set.to_owned(),
                    chr_position: chr_position.to_owned(),
                })?,
        );
        if called_unique.len() > 1 {
            effectively_phased = false;
        }
    }

    let vcf_call = render_vcf_call(genotype, allele1.as_deref(), allele2.as_deref(), phased);

    Ok((
        false,
        Some(SampleAlleleSummary {
            chromosome: chromosome.to_owned(),
            position,
            allele1,
            allele2,
            vcf_alleles,
            genotype: genotype.to_owned(),
            vcf_call,
            phased,
            effectively_phased,
            phase_set,
            undocumented_variations: BTreeSet::new(),
            treat_undocumented_variations_as_reference: false,
        }),
    ))
}

fn validate_alleles(
    chr_position: &str,
    reference: &str,
    alternates: &[String],
    called_genotypes: &[usize],
    warnings: &mut VcfWarnings,
) -> bool {
    let mut is_valid = true;

    if reference.starts_with('<') {
        add_warning(
            warnings,
            chr_position,
            format!(
                "Discarded genotype at this position because REF uses structural variation '{reference}'"
            ),
        );
        is_valid = false;
    } else if reference.to_uppercase().contains('N') {
        add_warning(
            warnings,
            chr_position,
            format!(
                "Discarded genotype at this position because REF uses ambiguous allele in '{reference}'"
            ),
        );
        is_valid = false;
    } else if !is_known_dna_allele(reference) {
        add_warning(
            warnings,
            chr_position,
            format!(
                "Discarded genotype at this position because REF uses unknown base in '{reference}'"
            ),
        );
        is_valid = false;
    }

    for (index, alternate) in alternates.iter().enumerate() {
        let allele_index = index + 1;
        let is_selected = called_genotypes.contains(&allele_index);
        let prefix = if is_selected {
            "Discarded genotype at this position because "
        } else {
            ""
        };

        if alternate.starts_with('<') {
            if is_selected {
                is_valid = false;
            }
            add_warning(
                warnings,
                chr_position,
                format!("{prefix}ALT uses structural variation '{alternate}'"),
            );
        } else if alternate.to_uppercase().contains('N') {
            if is_selected {
                is_valid = false;
            }
            add_warning(
                warnings,
                chr_position,
                format!("{prefix}ALT uses ambiguous allele in '{alternate}'"),
            );
        } else if alternate.contains('*') {
            if is_selected {
                is_valid = false;
            }
            add_warning(
                warnings,
                chr_position,
                format!("{prefix}ALT uses missing allele in '{alternate}'"),
            );
        } else if !is_known_dna_allele(alternate) {
            if is_selected {
                is_valid = false;
            }
            add_warning(
                warnings,
                chr_position,
                format!("{prefix}ALT uses unknown base in '{alternate}'"),
            );
        }
    }

    is_valid
}

fn is_known_dna_allele(allele: &str) -> bool {
    !allele.is_empty()
        && allele
            .bytes()
            .all(|base| matches!(base, b'A' | b'a' | b'C' | b'c' | b'G' | b'g' | b'T' | b't'))
}

fn parse_called_genotype_tokens(genotype_tokens: &[&str]) -> Result<Vec<usize>, ReadVcfError> {
    genotype_tokens
        .iter()
        .filter(|token| **token != ".")
        .map(|token| {
            token.parse().map_err(|_| ReadVcfError::InvalidGenotype {
                genotype: (*token).to_owned(),
            })
        })
        .collect()
}

fn genotype_token_to_allele(token: &str, vcf_alleles: &[String]) -> Option<String> {
    if token == "." {
        None
    } else {
        token
            .parse::<usize>()
            .ok()
            .and_then(|index| vcf_alleles.get(index))
            .cloned()
    }
}

fn render_vcf_call(
    genotype: &str,
    allele1: Option<&str>,
    allele2: Option<&str>,
    phased: bool,
) -> String {
    let mut call = String::new();

    if let Some(allele) = allele1 {
        call.push_str(allele);
    } else {
        call.push('.');
    }

    if let Some(allele) = allele2 {
        call.push(if phased { '|' } else { '/' });
        call.push_str(allele);
    } else if genotype.split(['|', '/']).nth(1) == Some(".") {
        call.push(if phased { '|' } else { '/' });
        call.push('.');
    }

    call
}

fn is_haploid_chromosome(chromosome: &str) -> bool {
    matches!(chromosome, "chrY" | "chrM")
}

fn parse_allelic_depths(allelic_depth: &str) -> Result<Vec<i32>, ()> {
    allelic_depth
        .split(',')
        .filter(|value| *value != ".")
        .map(|value| value.parse().map_err(|_| ()))
        .collect()
}

fn parse_called_genotype_alleles(genotype: &str) -> BTreeSet<i32> {
    genotype
        .split(['|', '/'])
        .filter(|allele| *allele != ".")
        .filter_map(|allele| allele.parse().ok())
        .collect()
}

fn add_warning(
    warnings: &mut VcfWarnings,
    chr_position: impl Into<String>,
    warning: impl Into<String>,
) {
    warnings
        .entry(chr_position.into())
        .or_default()
        .insert(warning.into());
}

/// Error returned when reading a VCF header fails.
#[derive(Debug)]
pub enum ReadHeaderError {
    /// Input path is not a supported VCF filename.
    NotVcf(String),
    /// File I/O failed.
    Io(io::Error),
    /// The VCF parser rejected the header.
    Parse(io::Error),
    /// Multiple contig assemblies were declared.
    MixedAssemblies(Vec<String>),
}

impl std::fmt::Display for ReadHeaderError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NotVcf(path) => write!(f, "{path} is not a VCF file"),
            Self::Io(err) => err.fmt(f),
            Self::Parse(err) => write!(f, "failed to parse VCF header: {err}"),
            Self::MixedAssemblies(assemblies) => {
                write!(f, "VCF file uses different assemblies")?;
                for (i, assembly) in assemblies.iter().enumerate() {
                    if i == 0 {
                        write!(f, " ({assembly}")?;
                    } else {
                        write!(f, " vs. {assembly}")?;
                    }
                }
                write!(f, " for contig)")
            }
        }
    }
}

impl std::error::Error for ReadHeaderError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(err) | Self::Parse(err) => Some(err),
            Self::NotVcf(_) | Self::MixedAssemblies(_) => None,
        }
    }
}

impl From<io::Error> for ReadHeaderError {
    fn from(err: io::Error) -> Self {
        Self::Io(err)
    }
}

/// Error returned when reading VCF records fails.
#[derive(Debug)]
pub enum ReadVcfError {
    /// Header validation failed.
    Header(ReadHeaderError),
    /// File I/O failed.
    Io(io::Error),
    /// The VCF did not declare any sample columns.
    NoSamplesDeclared,
    /// The requested sample was not declared in the header.
    SampleNotFound(String),
    /// A record was missing its 1-based position.
    MissingPosition,
    /// A record did not contain the selected sample column.
    MissingSampleColumn,
    /// A genotype token was not an integer or `.`.
    InvalidGenotype {
        /// Invalid genotype token.
        genotype: String,
    },
    /// A genotype allele index exceeded the available REF/ALT allele count.
    InvalidGenotypeAllele {
        /// Invalid allele index.
        allele: usize,
        /// Chromosome position.
        chr_position: String,
        /// Number of ALT alleles in the VCF row.
        alt_count: usize,
    },
    /// A phase set was present but not an integer.
    InvalidPhaseSet {
        /// Invalid phase set value.
        phase_set: String,
        /// Chromosome position.
        chr_position: String,
    },
}

impl std::fmt::Display for ReadVcfError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Header(err) => err.fmt(f),
            Self::Io(err) => err.fmt(f),
            Self::NoSamplesDeclared => write!(f, "VCF did not declare any samples"),
            Self::SampleNotFound(sample) => write!(f, "sample {sample} was not found in VCF"),
            Self::MissingPosition => write!(f, "VCF record is missing a position"),
            Self::MissingSampleColumn => write!(f, "VCF record is missing a sample column"),
            Self::InvalidGenotype { genotype } => write!(f, "invalid genotype token {genotype}"),
            Self::InvalidGenotypeAllele {
                allele,
                chr_position,
                alt_count,
            } => write!(
                f,
                "Invalid GT allele value ({allele}) for {chr_position} (only {alt_count} ALT allele{} specified)",
                if *alt_count == 1 { "" } else { "s" }
            ),
            Self::InvalidPhaseSet {
                phase_set,
                chr_position,
            } => write!(f, "invalid phase set {phase_set} for {chr_position}"),
        }
    }
}

impl std::error::Error for ReadVcfError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Header(err) => Some(err),
            Self::Io(err) => Some(err),
            Self::NoSamplesDeclared
            | Self::SampleNotFound(_)
            | Self::MissingPosition
            | Self::MissingSampleColumn
            | Self::InvalidGenotype { .. }
            | Self::InvalidGenotypeAllele { .. }
            | Self::InvalidPhaseSet { .. } => None,
        }
    }
}

impl From<ReadHeaderError> for ReadVcfError {
    fn from(err: ReadHeaderError) -> Self {
        Self::Header(err)
    }
}

impl From<io::Error> for ReadVcfError {
    fn from(err: io::Error) -> Self {
        Self::Io(err)
    }
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use super::{
        AllelicDepthPolicy, MSG_AD_FORMAT_MISSING, ReadHeaderError, ReadVcfError, is_vcf_file,
        read_header_summary, read_record_summaries,
    };

    #[test]
    fn reads_samples_and_shared_assembly_from_multisample_fixture() {
        let summary = read_header_summary(Path::new(
            "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/VcfSampleReaderTest-multisample.vcf",
        ))
        .expect("valid VCF header");

        assert_eq!(summary.samples, ["Sample_1", "Sample_2"]);
        assert_eq!(summary.genome_build.as_deref(), Some("GRCh38.p13"));
    }

    #[test]
    fn rejects_non_vcf_extension_like_java_vcf_file_wrapper() {
        let err = read_header_summary(Path::new("../../TODO.md")).expect_err("not a VCF");

        assert!(matches!(err, ReadHeaderError::NotVcf(_)));
    }

    #[test]
    fn recognizes_pharmcat_vcf_extensions() {
        assert!(is_vcf_file(Path::new("sample.vcf")));
        assert!(is_vcf_file(Path::new("sample.vcf.gz")));
        assert!(is_vcf_file(Path::new("sample.vcf.bgz")));
        assert!(!is_vcf_file(Path::new("sample.bcf")));
    }

    #[test]
    fn reads_record_fields_for_selected_sample() {
        let records = read_record_summaries(
            Path::new(
                "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/VcfSampleReaderTest-multisample.vcf",
            ),
            Some("Sample_2"),
        )
        .expect("valid records");

        assert_eq!(records.sample_name, "Sample_2");
        assert_eq!(records.header.samples, ["Sample_1", "Sample_2"]);
        assert_eq!(records.records.len(), 18);

        let first = &records.records[0];
        assert_eq!(first.chromosome, "chr1");
        assert_eq!(first.position, 97078987);
        assert_eq!(first.ids, ["rs114096998"]);
        assert_eq!(first.reference, "G");
        assert_eq!(first.alternates, ["T"]);
        assert_eq!(first.filters, ["PASS"]);
        assert_eq!(first.format_keys, ["GT"]);
        assert_eq!(first.sample.genotype.as_deref(), Some("0/0"));
    }

    #[test]
    fn accepts_ft_format_number_unknown_like_java_vcf_reader() {
        let records = read_record_summaries(
            Path::new(
                "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/NamedAlleleMatcher-unknownAltMultisample.vcf",
            ),
            Some("Sample_2"),
        )
        .expect("valid records");

        assert_eq!(records.sample_name, "Sample_2");
        assert_eq!(records.header.samples, ["Sample_1", "Sample_2"]);
        assert_eq!(records.records.len(), 2);
    }

    #[test]
    fn preserves_multiallelic_order_for_selected_sample() {
        let records = read_record_summaries(
            Path::new(
                "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/VcfSampleReaderTest-multisample.vcf",
            ),
            Some("Sample_2"),
        )
        .expect("valid records");

        let row = records
            .records
            .iter()
            .find(|record| record.chromosome == "chr2" && record.position == 233760233)
            .expect("UGT1A1 multiallelic row");

        assert_eq!(row.reference, "CAT");
        assert_eq!(row.alternates, ["CATAT", "CATATAT", "C"]);
        assert_eq!(row.sample.genotype.as_deref(), Some("2/3"));
        assert_eq!(
            row.allele_call.as_ref().expect("allele call").vcf_call,
            "CATATAT/C"
        );
    }

    #[test]
    fn reads_compressed_bgzf_fixture_and_collapses_to_java_allele_count() {
        let records = read_record_summaries(
            Path::new(
                "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/multisample.vcf.bgz",
            ),
            None,
        )
        .expect("valid compressed VCF");

        assert_eq!(records.header.samples, ["Sample_1", "Sample_2"]);
        assert_eq!(
            records
                .records
                .iter()
                .filter(|record| record.allele_call.is_some())
                .count(),
            15
        );
    }

    #[test]
    fn skips_duplicate_positions_and_applies_preprocessor_filter_warnings() {
        let records = read_record_summaries(
            Path::new(
                "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/VcfSampleReaderTest-multisample.vcf",
            ),
            Some("Sample_1"),
        )
        .expect("valid records");

        let duplicate_indel = records
            .records
            .iter()
            .find(|record| record.position == 97450065 && record.filters == ["PCATxINDEL"])
            .expect("duplicate preprocessor row");
        assert!(duplicate_indel.skipped_duplicate);
        assert!(duplicate_indel.allele_call.is_none());
        assert!(!records.warnings.contains_key("chr1:97450065"));

        let duplicate_ref_mismatch = records
            .records
            .iter()
            .find(|record| {
                record.chromosome == "chr7"
                    && record.position == 117509035
                    && record.reference == "GA"
            })
            .expect("PCATxREF row");
        assert!(duplicate_ref_mismatch.skipped_duplicate);
        assert!(!duplicate_ref_mismatch.discarded_by_preprocessor_filter);
        assert!(duplicate_ref_mismatch.allele_call.is_none());
        assert!(!records.warnings.contains_key("chr7:117509035"));
    }

    #[test]
    fn discards_nonduplicate_pcatxref_rows_like_java_no_definition_mode() {
        let records = read_record_summaries(
            Path::new(
                "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/VcfReaderTest-filters.vcf",
            ),
            None,
        )
        .expect("valid records");

        let ref_mismatch = records
            .records
            .iter()
            .find(|record| record.position == 94938683)
            .expect("PCATxREF row");
        assert!(ref_mismatch.discarded_by_preprocessor_filter);
        assert!(ref_mismatch.allele_call.is_none());
        assert_eq!(
            records
                .warnings
                .get("chr10:94938683")
                .expect("warning")
                .first()
                .map(String::as_str),
            Some(
                "Discarded genotype at this position because REF in VCF (G) does not match expected reference"
            )
        );
    }

    #[test]
    fn extracts_phase_set_when_format_declares_ps() {
        let records = read_record_summaries(
            Path::new(
                "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/VcfReaderTest-phaseSet.vcf",
            ),
            Some("Sample_1"),
        )
        .expect("valid records");

        let phased = records
            .records
            .iter()
            .find(|record| record.position == 40991369)
            .expect("phased row");
        assert_eq!(phased.sample.genotype.as_deref(), Some("1|0"));
        assert_eq!(phased.sample.phase_set.as_deref(), Some("40991224"));

        let unphased = records
            .records
            .iter()
            .find(|record| record.position == 40991381)
            .expect("row without PS");
        assert_eq!(unphased.sample.genotype.as_deref(), Some("0|0"));
        assert_eq!(unphased.sample.phase_set, None);
    }

    #[test]
    fn builds_sample_allele_calls_and_single_genotype_warnings_like_java() {
        let records = read_record_summaries(
            Path::new(
                "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/VcfReaderTest-alleleOrder.vcf",
            ),
            None,
        )
        .expect("valid records");

        assert_eq!(
            records
                .records
                .iter()
                .filter(|record| record.allele_call.is_some())
                .count(),
            4
        );

        let first = records
            .records
            .iter()
            .find(|record| record.position == 97078987)
            .and_then(|record| record.allele_call.as_ref())
            .expect("first allele call");
        assert_eq!(first.vcf_call, "G|T");
        assert!(first.phased);
        assert!(first.effectively_phased);

        let missing_first = records
            .records
            .iter()
            .find(|record| record.position == 97079005)
            .and_then(|record| record.allele_call.as_ref())
            .expect(".|1 call");
        assert_eq!(missing_first.allele1, None);
        assert_eq!(missing_first.allele2.as_deref(), Some("T"));
        assert_eq!(missing_first.vcf_call, ".|T");
        assert_eq!(
            records
                .warnings
                .get("chr1:97079005")
                .expect("warning")
                .first()
                .map(String::as_str),
            Some("1 genotype found (GT=.|1), expecting 2.")
        );

        let missing_second = records
            .records
            .iter()
            .find(|record| record.position == 97079071)
            .and_then(|record| record.allele_call.as_ref())
            .expect("1|. call");
        assert_eq!(missing_second.allele1.as_deref(), Some("A"));
        assert_eq!(missing_second.allele2, None);
        assert_eq!(missing_second.vcf_call, "A|.");

        assert!(
            records
                .records
                .iter()
                .find(|record| record.position == 97079076)
                .expect("0|. row")
                .allele_call
                .is_none()
        );
        assert_eq!(
            records
                .warnings
                .get("chr1:97079076")
                .expect("warning")
                .first()
                .map(String::as_str),
            Some(
                "Ignoring: only a single genotype found (GT=0|.).  Since it's reference, treating this as a missing position."
            )
        );
        assert!(
            records
                .records
                .iter()
                .find(|record| record.position == 97079077)
                .expect(".|0 row")
                .allele_call
                .is_none()
        );
        assert_eq!(
            records
                .warnings
                .get("chr1:97079077")
                .expect("warning")
                .first()
                .map(String::as_str),
            Some(
                "Ignoring: only a single genotype found (GT=.|0).  Since it's reference, treating this as a missing position."
            )
        );
    }

    #[test]
    fn matches_java_effective_phasing_for_unphased_heterozygous_rows() {
        let records = read_record_summaries(
            Path::new(
                "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/VcfReaderTest-phasing.vcf",
            ),
            None,
        )
        .expect("valid records");

        let calls: Vec<_> = records
            .records
            .iter()
            .map(|record| record.allele_call.as_ref().expect("allele call"))
            .collect();

        assert_eq!(calls.len(), 6);
        for call in calls {
            if call.chromosome == "chr3" {
                assert!(call.effectively_phased);
            } else {
                assert!(!call.effectively_phased);
            }
        }
    }

    #[test]
    fn keeps_phase_set_only_for_phased_sample_rows_like_java() {
        let sample_1 = read_record_summaries(
            Path::new(
                "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/VcfReaderTest-phaseSet.vcf",
            ),
            Some("Sample_1"),
        )
        .expect("valid records");
        assert_eq!(
            sample_1
                .records
                .iter()
                .filter(|record| record
                    .allele_call
                    .as_ref()
                    .and_then(|call| call.phase_set)
                    .is_some())
                .count(),
            6
        );
        assert_eq!(
            sample_1
                .records
                .iter()
                .find(|record| record.position == 40991367)
                .and_then(|record| record.allele_call.as_ref())
                .and_then(|call| call.phase_set),
            Some(40991224)
        );
        assert_eq!(
            sample_1
                .records
                .iter()
                .find(|record| record.position == 40991381)
                .and_then(|record| record.allele_call.as_ref())
                .and_then(|call| call.phase_set),
            None
        );

        let sample_2 = read_record_summaries(
            Path::new(
                "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/VcfReaderTest-phaseSet.vcf",
            ),
            Some("Sample_2"),
        )
        .expect("valid records");
        assert_eq!(
            sample_2
                .records
                .iter()
                .filter(|record| record
                    .allele_call
                    .as_ref()
                    .and_then(|call| call.phase_set)
                    .is_some())
                .count(),
            0
        );
    }

    #[test]
    fn rejects_unknown_sample_name() {
        let err = read_record_summaries(
            Path::new(
                "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/VcfSampleReaderTest-multisample.vcf",
            ),
            Some("missing"),
        )
        .expect_err("sample should be absent");

        assert!(matches!(err, ReadVcfError::SampleNotFound(sample) if sample == "missing"));
    }

    #[test]
    fn warns_for_unselected_structural_alt_and_discards_selected_structural_alt() {
        let records = read_record_summaries(
            Path::new(
                "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/VcfReaderTest-structuralAlt.vcf",
            ),
            None,
        )
        .expect("valid records");

        assert_eq!(records.warnings.len(), 6);
        for position in 1..=3 {
            let key = format!("chr1:{position}");
            let warning = records
                .warnings
                .get(&key)
                .expect("warning")
                .first()
                .expect("warning text");
            assert_eq!(warning, "ALT uses structural variation '<*>'");
            let row = records
                .records
                .iter()
                .find(|record| record.position == position)
                .expect("row");
            assert!(!row.discarded_by_allele_validation);
            assert!(row.allele_call.is_some());
        }

        for position in 4..=6 {
            let key = format!("chr1:{position}");
            let warning = records
                .warnings
                .get(&key)
                .expect("warning")
                .first()
                .expect("warning text");
            assert_eq!(
                warning,
                "Discarded genotype at this position because ALT uses structural variation '<*>'"
            );
            let row = records
                .records
                .iter()
                .find(|record| record.position == position)
                .expect("row");
            assert!(row.discarded_by_allele_validation);
            assert!(row.allele_call.is_none());
        }
    }

    #[test]
    fn rejects_gt_allele_index_past_declared_alts_like_java() {
        let err = read_record_summaries(
            Path::new(
                "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/VcfReaderTest-badAllele.vcf",
            ),
            None,
        )
        .expect_err("bad GT allele");

        assert!(matches!(
            err,
            ReadVcfError::InvalidGenotypeAllele {
                allele: 2,
                chr_position,
                alt_count: 1,
            } if chr_position == "chr2:233760233"
        ));
    }

    #[test]
    fn haploid_chromosomes_use_first_non_missing_genotype_like_java() {
        let records = read_record_summaries(Path::new("../../fixtures/vcf/haploid.vcf"), None)
            .expect("valid records");

        let chr_y_ref_first = records
            .records
            .iter()
            .find(|record| record.chromosome == "chrY" && record.position == 1)
            .and_then(|record| record.allele_call.as_ref())
            .expect("chrY 0/1 call");
        assert_eq!(chr_y_ref_first.allele1.as_deref(), Some("A"));
        assert_eq!(chr_y_ref_first.allele2, None);
        assert_eq!(chr_y_ref_first.vcf_call, "A");
        assert_eq!(
            records
                .warnings
                .get("chrY:1")
                .expect("warning")
                .first()
                .map(String::as_str),
            Some(
                "2 genotypes found (GT=0/1) for haploid chromosome. Will only use first non-missing genotype."
            )
        );

        let chr_y_alt_first = records
            .records
            .iter()
            .find(|record| record.chromosome == "chrY" && record.position == 2)
            .and_then(|record| record.allele_call.as_ref())
            .expect("chrY 1/0 call");
        assert_eq!(chr_y_alt_first.allele1.as_deref(), Some("G"));
        assert_eq!(chr_y_alt_first.allele2, None);
        assert_eq!(chr_y_alt_first.vcf_call, "G");

        let chr_y_missing_first = records
            .records
            .iter()
            .find(|record| record.chromosome == "chrY" && record.position == 3)
            .and_then(|record| record.allele_call.as_ref())
            .expect("chrY .|1 call");
        assert_eq!(chr_y_missing_first.allele1, None);
        assert_eq!(chr_y_missing_first.allele2.as_deref(), Some("G"));
        assert_eq!(chr_y_missing_first.vcf_call, ".|G");
        assert!(!records.warnings.contains_key("chrY:3"));

        let chr_m = records
            .records
            .iter()
            .find(|record| record.chromosome == "chrM" && record.position == 4)
            .and_then(|record| record.allele_call.as_ref())
            .expect("chrM 0/0 call");
        assert_eq!(chr_m.vcf_call, "C");
        assert_eq!(
            records
                .warnings
                .get("chrM:4")
                .expect("warning")
                .first()
                .map(String::as_str),
            Some(
                "2 genotypes found (GT=0/0) for haploid chromosome. Will only use first non-missing genotype."
            )
        );
    }

    #[test]
    fn warns_when_gt_heterozygous_but_ad_has_one_nonzero_depth() {
        let records = read_record_summaries(
            Path::new(
                "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/VcfReaderTest-AD.vcf",
            ),
            None,
        )
        .expect("valid records");

        assert_eq!(
            records.allelic_depth_policy,
            AllelicDepthPolicy::UseDefinedReferenceAlternate
        );
        assert_eq!(records.warnings.len(), 2);
        assert_eq!(
            records
                .warnings
                .get("chr10:94942254")
                .expect("warning")
                .first()
                .map(String::as_str),
            Some(
                "Discarding genotype at this position because GT field indicates heterozygous (0/1) but AD field indicates homozygous (91,0)"
            )
        );
        assert_eq!(
            records
                .warnings
                .get("chr10:94942255")
                .expect("warning")
                .first()
                .map(String::as_str),
            Some(
                "Discarding genotype at this position because GT field indicates heterozygous (0/1) but AD field indicates homozygous (0,91)"
            )
        );
        assert!(
            records
                .records
                .iter()
                .find(|record| record.position == 94942254)
                .expect("AD mismatch row")
                .discarded_by_allelic_depth
        );
        assert!(
            !records
                .records
                .iter()
                .find(|record| record.position == 94942249)
                .expect("homozygous row")
                .discarded_by_allelic_depth
        );
    }

    #[test]
    fn warns_once_when_ad_column_exists_without_format_header() {
        let records = read_record_summaries(
            Path::new(
                "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/VcfReaderTest-AD-missing.vcf",
            ),
            None,
        )
        .expect("valid records");

        assert_eq!(
            records.allelic_depth_policy,
            AllelicDepthPolicy::UseMissingDefinition
        );
        assert_eq!(records.warnings.len(), 3);
        assert_eq!(
            records
                .warnings
                .get("VCF")
                .expect("VCF warning")
                .first()
                .map(String::as_str),
            Some(MSG_AD_FORMAT_MISSING)
        );
    }

    #[test]
    fn ignores_ad_when_header_number_is_unknown() {
        let records = read_record_summaries(
            Path::new(
                "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/VcfReaderTest-AD-unknown.vcf",
            ),
            None,
        )
        .expect("valid records");

        assert_eq!(
            records.allelic_depth_policy,
            AllelicDepthPolicy::IgnoreUnknownNumber
        );
        assert!(records.warnings.is_empty());
        assert!(
            records
                .records
                .iter()
                .all(|record| !record.discarded_by_allelic_depth)
        );
    }

    #[test]
    fn warns_and_ignores_ad_when_header_number_is_invalid() {
        let records = read_record_summaries(
            Path::new(
                "../../repos/PharmCAT/src/test/resources/org/pharmgkb/pharmcat/haplotype/VcfReaderTest-AD-invalid.vcf",
            ),
            None,
        )
        .expect("valid records");

        assert_eq!(
            records.allelic_depth_policy,
            AllelicDepthPolicy::IgnoreInvalidNumber
        );
        assert_eq!(records.warnings.len(), 1);
        assert_eq!(
            records
                .warnings
                .get("VCF")
                .expect("VCF warning")
                .first()
                .map(String::as_str),
            Some(
                "INFO header for AD has unexpected number (f). Expecting 'R'. Treating number as '.' and ignoring AD field."
            )
        );
        assert!(
            records
                .records
                .iter()
                .all(|record| !record.discarded_by_allelic_depth)
        );
    }
}
