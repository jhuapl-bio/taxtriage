//
// Stage 2 of the pre-aligned (BAM/CRAM/SAM) input path.
//
// Samples that arrive as alignments skip everything upstream of ALIGNMENT — no
// compression, trimming, QC plots, host removal, classifier or reference
// download.  BAM_PREP has already normalised and indexed the alignment (and, if
// no reference FASTA was given, recovered consensus sequence from it); this
// subworkflow derives the coverage products in the exact channel shapes
// ALIGNMENT emits, so REPORT consumes both kinds of sample identically.
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

include { SAMTOOLS_COVERAGE as SAMTOOLS_COVERAGE_BAM } from '../../modules/nf-core/samtools/coverage/main'
include { SAMTOOLS_HIST_COVERAGE as SAMTOOLS_HIST_COVERAGE_BAM } from '../../modules/local/samtools_hist_coverage.nf'
include { BEDTOOLS_GENOMECOVERAGE as BEDTOOLS_GENOMECOVERAGE_BAM } from '../../modules/local/bedtools_genomcov'

workflow BAM_INPUT {
    take:
    ch_bams            // [ meta, bam, csi ] prepared alignments (from BAM_PREP)
    ch_mapping         // [ meta, merged accession->taxid map ] (from REFERENCE_PREP)

    main:
    ch_versions = Channel.empty()

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
        stats     = ch_stats                // [meta, coverage.txt]
        bedgraphs = ch_bedgraphs            // [meta, bedgraph | NO_FILE_bedgraph]
        versions  = ch_versions
}
