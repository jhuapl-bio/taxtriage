//
// Check input samplesheet and get read channels
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

include { KRONA_KTIMPORTTAXONOMY } from '../../modules/nf-core/krona/ktimporttaxonomy/main'
include { KRONA_KTIMPORTTEXT } from '../../modules/nf-core/krona/ktimporttext/main'
include { KRAKENTOOLS_KREPORT2KRONA  } from '../../modules/nf-core/krakentools/kreport2krona/main'


workflow VISUALIZE_REPORTS {
    take:
    report
    taxtab

    main:
        ch_krona = Channel.empty()

        // KRAKENTOOLS_KREPORT2KRONA(
        //     report
        // )

        // KRONA_KTIMPORTTEXT(
        //     KRAKENTOOLS_KREPORT2KRONA.out.txt
        // )
        // KRONA_KTIMPORTTEXT.out.html.set{
        //     ch_krona
        // }

        KRONA_KTIMPORTTAXONOMY(
            report,
            taxtab
        )
        KRONA_KTIMPORTTAXONOMY.out.html.set{
            ch_krona
        }

    emit:
        krona = ch_krona

}
