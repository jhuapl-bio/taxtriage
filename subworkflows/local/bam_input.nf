//
// Pre-aligned (BAM/CRAM/SAM) input path.
//
// Samples that arrive as alignments skip everything upstream of ALIGNMENT — no
// compression, trimming, QC plots, host removal, classifier or reference
// download.  This subworkflow simply normalises the alignment file and produces
// the exact channel shapes ALIGNMENT emits, so REPORT can consume both kinds of
// sample without knowing which is which.
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
include { SAMTOOLS_COVERAGE as SAMTOOLS_COVERAGE_BAM } from '../../modules/nf-core/samtools/coverage/main'
include { SAMTOOLS_HIST_COVERAGE as SAMTOOLS_HIST_COVERAGE_BAM } from '../../modules/local/samtools_hist_coverage.nf'
include { BEDTOOLS_GENOMECOVERAGE as BEDTOOLS_GENOMECOVERAGE_BAM } from '../../modules/local/bedtools_genomcov'

workflow BAM_INPUT {
    take:
    ch_alignments      // [ meta, bamfile ]
    ch_mapping         // [ meta, merged accession->taxid map ] (from REFERENCE_PREP)
    ch_cram_reference  // path: reference FASTA used for CRAM decoding, or NO_FILE

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
    // meta-keyed joins are done — see `counts` below and its use in taxtriage.nf.
    ch_bams = PREPARE_BAM.out.bam

    SAMTOOLS_COVERAGE_BAM(
        ch_bams
    )
    ch_versions = ch_versions.mix(SAMTOOLS_COVERAGE_BAM.out.versions)
    ch_stats = SAMTOOLS_COVERAGE_BAM.out.coverage

    // Per-accession coverage histogram (published only; mirrors ALIGNMENT)
    SAMTOOLS_HIST_COVERAGE_BAM(
        ch_bams.map { meta, bam, csi -> [meta, bam] }
            .join(ch_mapping, remainder: true)
            .map { meta, bam, mapping -> [meta, bam, mapping ?: []] }
            .filter { it[1] }
    )
    ch_versions = ch_versions.mix(SAMTOOLS_HIST_COVERAGE_BAM.out.versions)

    // Bedgraph is only consumed by match_paths.py in --fast mode; otherwise the
    // sourmash shared-window comparison is used and a placeholder is enough.
    if (params.fast) {
        BEDTOOLS_GENOMECOVERAGE_BAM(
            ch_bams.map { meta, bam, csi -> [meta, bam] }
        )
        ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOVERAGE_BAM.out.versions)
        ch_bedgraphs = BEDTOOLS_GENOMECOVERAGE_BAM.out.bedgraph
    } else {
        ch_bedgraphs = ch_bams.map { meta, bam, csi ->
            [meta, file("$projectDir/assets/NO_FILE_bedgraph")]
        }
    }

    emit:
        bams      = ch_bams                 // [meta, bam, csi]
        stats     = ch_stats                // [meta, coverage.txt]
        bedgraphs = ch_bedgraphs            // [meta, bedgraph | NO_FILE_bedgraph]
        counts    = PREPARE_BAM.out.count   // [meta, readcount.txt] -> meta.read_count downstream
        versions  = ch_versions
}
