use rust_htslib::bam::{self, record::Cigar};
use std::collections::HashSet;

use super::index::collect_bins;
use super::types::{Annotation, Call, Interval, Strand, Transcript};

impl Annotation {
    pub fn call(
        &self,
        contig_name: &str,
        blocks: &[Interval],
        reverse: bool,
        cigar: &str,
        include_introns: bool,
        three_prime: bool,
    ) -> Call {
        let Some(&contig) = self.contig_ids.get(contig_name) else {
            return Call {
                region: b'I',
                genes: vec![],
                tx: vec![],
                antisense: vec![],
                terminal_genes: vec![],
                terminal_tx: vec![],
            };
        };
        if blocks.is_empty() {
            return Call {
                region: b'I',
                genes: vec![],
                tx: vec![],
                antisense: vec![],
                terminal_genes: vec![],
                terminal_tx: vec![],
            };
        }
        let span = Interval {
            start: blocks[0].start,
            end: blocks.last().unwrap().end,
        };
        let mut tx_candidates = HashSet::new();
        collect_bins(&self.indexes[contig].exon_bins, span, &mut tx_candidates);
        let mut sense_exon_genes = HashSet::new();
        let mut any_exon_overlap = false;
        let mut sense_tx = Vec::new();
        let mut antisense_tx = Vec::new();
        let mut terminal_tx = Vec::new();
        let mut terminal_genes = HashSet::new();

        for tx_index in tx_candidates {
            let tx = &self.transcripts[tx_index];
            let sense = tx.strand.matches_record(reverse);
            if tx
                .exons
                .iter()
                .any(|&e| blocks.iter().any(|&b| e.overlaps(b)))
            {
                any_exon_overlap = true;
            }
            if let Some(pos) = compatible_position(tx, blocks) {
                let entry = format!(
                    "{},{}{},{}",
                    tx.id,
                    if sense { '+' } else { '-' },
                    pos,
                    transcript_cigar(cigar, tx.strand)
                );
                if sense {
                    sense_exon_genes.insert(tx.gene);
                    if three_prime && overlaps_terminal_exon(tx, blocks) {
                        terminal_genes.insert(tx.gene);
                        terminal_tx.push(entry.clone());
                    }
                    sense_tx.push(entry);
                } else {
                    antisense_tx.push(entry);
                }
            }
        }

        sense_tx.sort_unstable();
        sense_tx.dedup();
        antisense_tx.sort_unstable();
        antisense_tx.dedup();
        terminal_tx.sort_unstable();
        terminal_tx.dedup();

        let mut body_candidates = HashSet::new();
        collect_bins(&self.indexes[contig].gene_bins, span, &mut body_candidates);
        let mut sense_body = Vec::new();
        let mut antisense_body = Vec::new();
        for gene_index in body_candidates {
            let gene = &self.genes[gene_index];
            if !gene.span.overlaps(span) {
                continue;
            }
            if gene.strand.matches_record(reverse) {
                sense_body.push(gene_index);
            } else {
                antisense_body.push(gene_index);
            }
        }

        let exonic = any_exon_overlap;
        let intronic = !exonic && (!sense_body.is_empty() || !antisense_body.is_empty());
        let mut genes: Vec<usize> = if !sense_exon_genes.is_empty() {
            sense_exon_genes.into_iter().collect()
        } else if intronic && include_introns && !sense_body.is_empty() {
            sense_body.clone()
        } else {
            Vec::new()
        };
        genes.sort_unstable_by(|&x, &y| self.genes[x].id.cmp(&self.genes[y].id));
        genes.dedup();
        let mut terminal_genes: Vec<usize> = terminal_genes.into_iter().collect();
        terminal_genes.sort_unstable_by(|&x, &y| self.genes[x].id.cmp(&self.genes[y].id));
        terminal_genes.dedup();

        if sense_tx.is_empty() && intronic && include_introns {
            sense_tx = genes
                .iter()
                .map(|&g| format!("{},+", self.genes[g].id))
                .collect();
        }
        if antisense_tx.is_empty() && intronic && !antisense_body.is_empty() && include_introns {
            antisense_body.sort_unstable_by(|&x, &y| self.genes[x].id.cmp(&self.genes[y].id));
            antisense_body.dedup();
            antisense_tx = antisense_body
                .iter()
                .map(|&g| format!("{},-", self.genes[g].id))
                .collect();
        }
        Call {
            region: if exonic {
                b'E'
            } else if intronic {
                b'N'
            } else {
                b'I'
            },
            genes,
            tx: sense_tx,
            antisense: antisense_tx,
            terminal_genes,
            terminal_tx,
        }
    }
}

pub fn compatible_position(tx: &Transcript, blocks: &[Interval]) -> Option<i64> {
    let mut exon_indexes = Vec::with_capacity(blocks.len());
    for &block in blocks {
        exon_indexes.push(tx.exons.iter().position(|&exon| exon.contains(block))?);
    }
    for (pair, exon_pair) in blocks.windows(2).zip(exon_indexes.windows(2)) {
        if exon_pair[1] != exon_pair[0] + 1
            || pair[0].end != tx.exons[exon_pair[0]].end
            || pair[1].start != tx.exons[exon_pair[1]].start
        {
            return None;
        }
    }
    let first = exon_indexes[0];
    let last = *exon_indexes.last()?;
    Some(match tx.strand {
        Strand::Plus => {
            tx.exons[..first]
                .iter()
                .map(|e| e.end - e.start)
                .sum::<i64>()
                + blocks[0].start
                - tx.exons[first].start
        }
        Strand::Minus => {
            tx.exons[last + 1..]
                .iter()
                .map(|e| e.end - e.start)
                .sum::<i64>()
                + tx.exons[last].end
                - blocks.last()?.end
        }
    })
}

pub fn overlaps_terminal_exon(tx: &Transcript, blocks: &[Interval]) -> bool {
    let terminal = match tx.strand {
        Strand::Plus => tx.exons.last(),
        Strand::Minus => tx.exons.first(),
    };
    terminal.is_some_and(|&exon| blocks.iter().any(|&block| exon.overlaps(block)))
}

pub fn transcript_cigar(cigar: &str, strand: Strand) -> String {
    let mut ops: Vec<(u64, char)> = Vec::new();
    let mut digits = String::new();
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            digits.push(ch);
            continue;
        }
        if ch != 'N' {
            let length = digits.parse::<u64>().unwrap_or(0);
            if let Some((previous_length, previous_op)) = ops.last_mut() {
                if *previous_op == ch {
                    *previous_length += length;
                } else {
                    ops.push((length, ch));
                }
            } else {
                ops.push((length, ch));
            }
        }
        digits.clear();
    }
    if strand == Strand::Minus {
        ops.reverse();
    }
    ops.into_iter().map(|(n, op)| format!("{n}{op}")).collect()
}

pub fn alignment_blocks(record: &bam::Record) -> Vec<Interval> {
    let mut pos = record.pos();
    let mut start = pos;
    let mut blocks = Vec::new();
    for op in record.cigar().iter() {
        match *op {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) | Cigar::Del(n) => pos += n as i64,
            Cigar::RefSkip(n) => {
                if pos > start {
                    blocks.push(Interval { start, end: pos });
                }
                pos += n as i64;
                start = pos;
            }
            Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
    if pos > start {
        blocks.push(Interval { start, end: pos });
    }
    blocks
}
