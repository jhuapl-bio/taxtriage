//
// Stage 1 of the pre-aligned (BAM/CRAM/SAM) input path.
//
// Runs BEFORE REFERENCE_PREP because its consensus output feeds into it: when no
// --reference_fasta is supplied, the reference sequence is reconstructed from the
// alignment itself and then fuzzy-matched to the assembly summary exactly like a
// user-supplied FASTA would be.  Coverage products are computed afterwards in
// BAM_INPUT, which needs the map REFERENCE_PREP produces — hence the split.
//
// ##############################################################################################
// # Copyright 2022 The Johns Hopkins University Applied Physics Laboratory LLC
// # All rights reserved.
// # Permission is hereby granted, free of charge, to any person obtaining a copy of this
// # software and associated documentation files (the "Software"), to deal in the Software
// # without restriction, including without limitation the rights to use, copy, modify,
// # merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
// # permit persons to whom the Software is furnished to do so.
// #
// # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
// # INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
// # PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// # LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
// # TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
// # OR OTHER DEALINGS IN THE SOFTWARE.
// #

include { PREPARE_BAM } from '../../modules/local/prepare_bam'
include { BAM_CONSENSUS } from '../../modules/local/bam_consensus'

workflow BAM_PREP {
    take:
    ch_alignments      // [ meta, alignment file ]
    ch_cram_reference  // path: reference FASTA used for CRAM decoding, or NO_FILE
    val_consensus      // boolean: reconstruct reference sequence from the BAM

    main:
    ch_versions = Channel.empty()

    PREPARE_BAM(
        ch_alignments,
        ch_cram_reference
    )
    ch_versions = ch_versions.mix(PREPARE_BAM.out.versions)

    // [meta, bam, csi] — same shape as ALIGNMENT.out.bams.
    // NOTE: meta is deliberately left untouched here.  It is the join key for
    // every downstream channel (mapping, fastas, kraken placeholder …), so the
    // primary-read count is folded into meta.read_count only once all of those
    // meta-keyed joins are done — see `counts` and its use in taxtriage.nf.
    ch_bams = PREPARE_BAM.out.bam

    // ── Reference sequence recovery ────────────────────────────────────────
    // Only when the user gave us no reference to compare against.  match_paths.py
    // can read reference LENGTHS from the BAM header on its own, but sourmash
    // sketching, the shared-window conflict report and the ANI matrix all need
    // bases, so we call them off the alignment.  Empty results (no reference
    // recovered enough sequence) are dropped rather than passed on as an empty
    // FASTA that every sketch step would silently skip.
    ch_consensus = Channel.empty()
    if (val_consensus) {
        BAM_CONSENSUS(ch_bams)
        ch_versions = ch_versions.mix(BAM_CONSENSUS.out.versions)
        ch_consensus = BAM_CONSENSUS.out.fasta.filter { meta, fasta -> fasta.size() > 0 }
    }

    emit:
        bams      = ch_bams                  // [meta, bam, csi]
        counts    = PREPARE_BAM.out.count    // [meta, readcount.txt]
        consensus = ch_consensus             // [meta, consensus.fasta] (empty unless recovering)
        versions  = ch_versions
}
