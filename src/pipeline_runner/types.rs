/// Running counts for the pipeline summary.
#[derive(Default)]
pub(crate) struct Stats {
    pub reads: u64,
    pub oriented: u64,
    pub pairs_emitted: u64,
    pub r1_eq_28: u64,
    pub r1_not_28: u64,
    pub adapter_trimmed: u64,
    pub r2_len_sum: u64,
    pub r2_len_max: u64,
    // read classes
    pub class_full: u64,
    pub class_partial: u64,
    pub class_unparseable: u64,
    // cb gate
    pub cb_on_whitelist: u64,
    pub cb_within1: u64,
    pub cb_fail: u64,
    pub gated_out: u64,
    pub realigned: u64,
    pub reframed: u64,
    pub corrected2: u64,
    /// CBs whose location is dual-anchor confirmed (R1-tail flank + polyT flank).
    pub cb_located_confident: u64,
    /// CB-fails that are nonetheless confidently located (= genuine base errors).
    pub cb_fail_confident: u64,
    // demux
    pub undetermined: u64,
    /// Reads written to the `--failed-output` side channel.
    pub failed_emitted: u64,
    pub via_dual: u64,
    pub via_i5: u64,
    pub via_i7: u64,
    pub via_conflict: u64,
    pub per_sample: Vec<u64>,
    // demux (early): the same i5/i7 assignment, computed the moment R1 is
    // located -- before the CB gate and independent of whether R2/cDNA was
    // ever found. Reachable for every read that got as far as R1, unlike
    // `per_sample`/`undetermined` above (only reachable for reads that also
    // cleared both gates). This is the closest SBX-side analog to Illumina's
    // own post-BCL-Convert index-only demultiplex count -- see
    // docs/archive/FLEX_SUPPORT_PLAN.md's read-metrics comparison.
    pub per_sample_early: Vec<u64>,
    /// Of `per_sample_early`, how many also located a real R2/cDNA 3' boundary
    /// (`Extracted::adapter_trimmed`) -- excludes reads that are structurally
    /// too short/garbled to ever locate one, regardless of matching-algorithm
    /// quality. The closest fair denominator for "of reads that had a real
    /// shot at this, how many succeeded" when scoring a matcher in isolation
    /// from upstream capture-length loss.
    pub per_sample_early_boundary_found: Vec<u64>,
    pub undetermined_early: u64,
    pub via_dual_early: u64,
    pub via_i5_only_early: u64,
    pub via_i7_only_early: u64,
    pub via_conflict_early: u64,
    /// R1 itself was never located (orientation failed, or R1 not found even
    /// though oriented) -- demux is structurally unreachable for these, same
    /// as for Illumina's own index read, which requires nothing about R1/R2
    /// content but does require the read to exist at all.
    pub index_unreachable: u64,
    // probe barcode
    pub probe_located: u64,
    pub probe_hit: u64,
    /// `probe_hit` reads that *also* landed on an assigned (non-Undetermined)
    /// sample from the outer index-based demux -- the real 4-way join (CB +
    /// UMI + Sample Barcode + Probe Assigned) requested to replace
    /// `probe_hit` alone as a "pre-dedup" proxy, since `probe_hit` on its own
    /// still includes reads whose outer sample/library demux was
    /// Undetermined.
    pub probe_hit_sample_assigned: u64,
}

impl Stats {
    /// A zeroed `Stats` sized for `n_samples` per-sample counters.
    pub(crate) fn new(n_samples: usize) -> Self {
        Stats {
            per_sample: vec![0; n_samples],
            per_sample_early: vec![0; n_samples],
            per_sample_early_boundary_found: vec![0; n_samples],
            ..Default::default()
        }
    }

    /// Accumulate a per-file `Stats` into this total.
    pub(crate) fn merge(&mut self, o: Stats) {
        self.reads += o.reads;
        self.oriented += o.oriented;
        self.pairs_emitted += o.pairs_emitted;
        self.r1_eq_28 += o.r1_eq_28;
        self.r1_not_28 += o.r1_not_28;
        self.adapter_trimmed += o.adapter_trimmed;
        self.r2_len_sum += o.r2_len_sum;
        self.r2_len_max = self.r2_len_max.max(o.r2_len_max);
        self.class_full += o.class_full;
        self.class_partial += o.class_partial;
        self.class_unparseable += o.class_unparseable;
        self.cb_on_whitelist += o.cb_on_whitelist;
        self.cb_within1 += o.cb_within1;
        self.cb_fail += o.cb_fail;
        self.gated_out += o.gated_out;
        self.realigned += o.realigned;
        self.reframed += o.reframed;
        self.corrected2 += o.corrected2;
        self.cb_located_confident += o.cb_located_confident;
        self.cb_fail_confident += o.cb_fail_confident;
        self.undetermined += o.undetermined;
        self.failed_emitted += o.failed_emitted;
        self.via_dual += o.via_dual;
        self.via_i5 += o.via_i5;
        self.via_i7 += o.via_i7;
        self.via_conflict += o.via_conflict;
        self.probe_located += o.probe_located;
        self.probe_hit += o.probe_hit;
        self.probe_hit_sample_assigned += o.probe_hit_sample_assigned;
        self.undetermined_early += o.undetermined_early;
        self.index_unreachable += o.index_unreachable;
        self.via_dual_early += o.via_dual_early;
        self.via_i5_only_early += o.via_i5_only_early;
        self.via_i7_only_early += o.via_i7_only_early;
        self.via_conflict_early += o.via_conflict_early;
        if self.per_sample.len() < o.per_sample.len() {
            self.per_sample.resize(o.per_sample.len(), 0);
        }
        for (a, b) in self.per_sample.iter_mut().zip(o.per_sample.iter()) {
            *a += b;
        }
        if self.per_sample_early.len() < o.per_sample_early.len() {
            self.per_sample_early.resize(o.per_sample_early.len(), 0);
        }
        for (a, b) in self
            .per_sample_early
            .iter_mut()
            .zip(o.per_sample_early.iter())
        {
            *a += b;
        }
        if self.per_sample_early_boundary_found.len() < o.per_sample_early_boundary_found.len() {
            self.per_sample_early_boundary_found
                .resize(o.per_sample_early_boundary_found.len(), 0);
        }
        for (a, b) in self
            .per_sample_early_boundary_found
            .iter_mut()
            .zip(o.per_sample_early_boundary_found.iter())
        {
            *a += b;
        }
    }
}

/// Wall-clock breakdown of a run, for the summary log and metrics JSON.
pub(crate) struct Timing {
    pub elapsed: f64,
    pub load_secs: f64,
    pub process_secs: f64,
    pub threads: usize,
}

/// Compact, versioned SBX provenance carried from extraction into the alignment.
#[derive(Debug, Clone, Copy)]
pub(crate) struct ProcessingTags {
    pub orientation: u8,
    pub barcode_route: u8,
    pub index_route: u8,
    pub anchor_bits: u8,
    pub extracted_cdna_len: usize,
}

/// Records read from one file, flat-packed for cheap cross-thread transfer.
pub(crate) struct Batch {
    pub orientation: Option<crate::pipeline::Orientation>,
    pub ids: Vec<u8>,
    pub seqs: Vec<u8>,
    pub quals: Vec<u8>,
    /// `(id_len, seq_len)` per record; `qual_len == seq_len` for FASTQ.
    pub lens: Vec<(u32, u32)>,
}

/// Records per batch — big enough to amortize channel/sync overhead.
pub(crate) const BATCH_READS: usize = 4096;

impl Batch {
    pub(crate) fn with_orientation(orientation: Option<crate::pipeline::Orientation>) -> Self {
        Batch {
            orientation,
            ids: Vec::new(),
            seqs: Vec::new(),
            quals: Vec::new(),
            lens: Vec::with_capacity(BATCH_READS),
        }
    }

    pub(crate) fn push(&mut self, id: &[u8], seq: &[u8], qual: &[u8]) {
        self.ids.extend_from_slice(id);
        self.seqs.extend_from_slice(seq);
        self.quals.extend_from_slice(qual);
        self.lens.push((id.len() as u32, seq.len() as u32));
    }

    pub(crate) fn len(&self) -> usize {
        self.lens.len()
    }
}
