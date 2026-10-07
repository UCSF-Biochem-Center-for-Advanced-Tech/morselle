use std::collections::HashMap;
use std::fs;
use std::io::BufRead;
use std::path::{Path, PathBuf};
use anyhow::{bail, Context, Result};
use tracing::info;

use super::types::{
    Annotation, ContigIndex, Gene, Interval, Strand, Transcript, BIN_SIZE,
};

pub fn add_to_bins(bins: &mut Vec<Vec<usize>>, interval: Interval, value: usize) {
    let first = (interval.start / BIN_SIZE) as usize;
    let last = ((interval.end - 1) / BIN_SIZE) as usize;
    if bins.len() <= last {
        bins.resize_with(last + 1, Vec::new);
    }
    for bin in &mut bins[first..=last] {
        bin.push(value);
    }
}

pub fn collect_bins(bins: &[Vec<usize>], interval: Interval, output: &mut std::collections::HashSet<usize>) {
    let first = (interval.start / BIN_SIZE).max(0) as usize;
    let last = ((interval.end - 1) / BIN_SIZE).max(0) as usize;
    for bin in first..=last.min(bins.len().saturating_sub(1)) {
        output.extend(bins[bin].iter().copied());
    }
}

pub fn genes_mtx_bytes(annotation: &Annotation) -> Vec<u8> {
    let mut output = Vec::new();
    for gene in &annotation.genes {
        output.extend_from_slice(gene.id.as_bytes());
        output.push(b'\t');
        output.extend_from_slice(gene.name.as_bytes());
        output.extend_from_slice(b"\tGene Expression\n");
    }
    output
}

pub fn establish_genes_mtx(path: &Path, annotation: &Annotation) -> Result<()> {
    let expected = genes_mtx_bytes(annotation);
    if path.exists() {
        let observed = fs::read(path).with_context(|| format!("reading {}", path.display()))?;
        if observed != expected {
            bail!(
                "existing gene manifest {} does not match the filtered annotation; use a matching GTF/rule set or choose a new path",
                path.display()
            );
        }
        info!(path = %path.display(), genes = annotation.genes.len(), "verified canonical gene manifest");
        return Ok(());
    }
    if let Some(parent) = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        fs::create_dir_all(parent)?;
    }
    let temporary = PathBuf::from(format!("{}.tmp.{}", path.display(), std::process::id()));
    fs::write(&temporary, &expected).with_context(|| format!("writing {}", temporary.display()))?;
    fs::rename(&temporary, path).with_context(|| format!("creating {}", path.display()))?;
    info!(path = %path.display(), genes = annotation.genes.len(), "created canonical gene manifest");
    Ok(())
}

pub fn text_reader(path: &Path) -> Result<Box<dyn BufRead>> {
    crate::io::open_text_reader(path)
}

pub fn attribute<'a>(attrs: &'a str, key: &str) -> Option<&'a str> {
    attrs.split(';').find_map(|field| {
        let field = field.trim();
        let (k, value) = field.split_once(char::is_whitespace)?;
        (k == key).then(|| value.trim().trim_matches('"'))
    })
}

pub fn has_attribute_value(attrs: &str, key: &str, wanted: &str) -> bool {
    attrs.split(';').any(|field| {
        let field = field.trim();
        let Some((k, value)) = field.split_once(char::is_whitespace) else {
            return false;
        };
        k == key && value.trim().trim_matches('"') == wanted
    })
}

pub fn allowed_biotype(value: &str) -> bool {
    matches!(
        value,
        "protein_coding"
            | "protein_coding_LoF"
            | "lncRNA"
            | "antisense"
            | "IG_C_gene"
            | "IG_D_gene"
            | "IG_J_gene"
            | "IG_LV_gene"
            | "IG_V_gene"
            | "IG_V_pseudogene"
            | "IG_J_pseudogene"
            | "IG_C_pseudogene"
            | "TR_C_gene"
            | "TR_D_gene"
            | "TR_J_gene"
            | "TR_V_gene"
            | "TR_V_pseudogene"
            | "TR_J_pseudogene"
    )
}

pub fn normalize_ensembl_id(value: &str, strip_version: bool) -> String {
    if strip_version && (value.starts_with("ENSG") || value.starts_with("ENST")) {
        if let Some((base, version)) = value.rsplit_once('.') {
            if version.chars().all(|c| c.is_ascii_digit()) {
                return base.to_owned();
            }
        }
    }
    value.to_owned()
}

impl Annotation {
    pub fn load(path: &Path, reference_filter: bool, strip_id_versions: bool) -> Result<Self> {
        let mut a = Self::default();
        let mut genes_by_id: HashMap<String, usize> = HashMap::new();
        let mut tx_by_id: HashMap<String, usize> = HashMap::new();

        for (line_no, result) in text_reader(path)?.lines().enumerate() {
            let line = result.with_context(|| format!("reading GTF line {}", line_no + 1))?;
            if line.starts_with('#') || line.is_empty() {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() != 9 || fields[2] != "exon" {
                continue;
            }
            a.exons_seen += 1;
            let Some(strand) = Strand::parse(fields[6]) else {
                continue;
            };
            let start = fields[3]
                .parse::<i64>()
                .with_context(|| format!("invalid GTF start at line {}", line_no + 1))?
                - 1;
            let end = fields[4]
                .parse::<i64>()
                .with_context(|| format!("invalid GTF end at line {}", line_no + 1))?;
            if start < 0 || end <= start {
                bail!("invalid GTF interval at line {}", line_no + 1);
            }
            let attrs = fields[8];
            if reference_filter {
                let gene_type =
                    attribute(attrs, "gene_type").or_else(|| attribute(attrs, "gene_biotype"));
                let transcript_type = attribute(attrs, "transcript_type")
                    .or_else(|| attribute(attrs, "transcript_biotype"));
                if !gene_type.is_some_and(allowed_biotype)
                    || !transcript_type.is_some_and(allowed_biotype)
                {
                    a.exons_filtered_biotype += 1;
                    continue;
                }
                if has_attribute_value(attrs, "tag", "readthrough_transcript")
                    || has_attribute_value(attrs, "tag", "PAR")
                {
                    a.exons_filtered_tag += 1;
                    continue;
                }
            }
            let gene_id = attribute(attrs, "gene_id")
                .with_context(|| format!("exon lacks gene_id at line {}", line_no + 1))?;
            let gene_name = attribute(attrs, "gene_name").unwrap_or(gene_id);
            let tx_id = attribute(attrs, "transcript_id")
                .with_context(|| format!("exon lacks transcript_id at line {}", line_no + 1))?;
            let gene_id = normalize_ensembl_id(gene_id, strip_id_versions);
            let tx_id = normalize_ensembl_id(tx_id, strip_id_versions);
            let contig = if let Some(&id) = a.contig_ids.get(fields[0]) {
                id
            } else {
                let id = a.contig_names.len();
                a.contig_names.push(fields[0].to_owned());
                a.contig_ids.insert(fields[0].to_owned(), id);
                a.indexes.push(ContigIndex::default());
                id
            };
            let exon = Interval { start, end };
            let gene = if let Some(&index) = genes_by_id.get(&gene_id) {
                let g = &mut a.genes[index];
                if g.contig != contig || g.strand != strand {
                    bail!("gene {gene_id} occurs on inconsistent contig/strand");
                }
                g.span.start = g.span.start.min(start);
                g.span.end = g.span.end.max(end);
                index
            } else {
                let index = a.genes.len();
                a.genes.push(Gene {
                    id: gene_id.clone(),
                    name: gene_name.to_owned(),
                    contig,
                    strand,
                    span: exon,
                });
                genes_by_id.insert(gene_id.clone(), index);
                index
            };
            let tx = if let Some(&index) = tx_by_id.get(&tx_id) {
                index
            } else {
                let index = a.transcripts.len();
                a.transcripts.push(Transcript {
                    id: tx_id.clone(),
                    gene,
                    contig,
                    strand,
                    exons: Vec::new(),
                });
                tx_by_id.insert(tx_id, index);
                index
            };
            a.transcripts[tx].exons.push(exon);
        }

        for tx in &mut a.transcripts {
            tx.exons.sort_unstable_by_key(|e| (e.start, e.end));
            tx.exons.dedup();
        }
        for (tx_index, tx) in a.transcripts.iter().enumerate() {
            for &exon in &tx.exons {
                add_to_bins(&mut a.indexes[tx.contig].exon_bins, exon, tx_index);
            }
        }
        for (gene_index, gene) in a.genes.iter().enumerate() {
            add_to_bins(&mut a.indexes[gene.contig].gene_bins, gene.span, gene_index);
        }
        Ok(a)
    }
}
