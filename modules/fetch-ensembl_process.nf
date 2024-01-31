process SampleSheet_SourceEnsembl {
    input:
    path samplesheet

    output:
    path "samplesheet_processed.csv"

    script:
    """
    ensembl_fetch.py -s ${samplesheet} -d ${params.out}
    """
}