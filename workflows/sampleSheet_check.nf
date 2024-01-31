workflow SAMPLESHEET_CHECK {
    take:
    samplesheet

    main:

    // Channel.from( samplesheet )
    samplesheet.splitCsv ( header:true, sep:',' )
        .map { create_metadata_channel(it) }
        .set { data_meta }

    emit:
    data_meta                              
}   

def create_metadata_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.name = row.name
    meta.id = row.id
    // meta.id = "${row.name}.${row.id}"
    meta.annotation_version = row.annotation_version
    meta.ensembl_db = row.ensembl_db
    meta.index_Folder = "${params.out}/${row.name}/${row.id}" 
    meta.gtf_path = row.gtf_path
    meta.fna_path = row.fna_path

    if (!file(row.gtf_path).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> gtf file does not exist!\n${row.gtf_path}"
    }
    if (!file(row.fna_path).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> fna file does not exist!\n${row.fna_path}"
    }
    
    return meta 
    }