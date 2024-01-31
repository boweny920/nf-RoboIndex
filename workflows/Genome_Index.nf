include { refFlat_geneData_Make; STARMake; RSEM_Make ; bedFiles_Make; genomeBPtype_Make; riboList_Make; transcriptData_Make; ens_link_Make_for_SIMR } from '../modules/Core_process.nf'

// Additional stuff to build based on requests
include { faidx_Make; bt2Index_Make; bwaIndex_Make; cellranger_singlecell_Index_Make; igv_genome_Make }  from '../modules/requested_process.nf'

workflow GENOME_INDEX_BUILD {
    take:
    meta
    
    main:
    //Core Processes: 
    refFlat_geneData_Make(meta) 
    // refFlat_geneData_Make.out.view() // One row of the sample sheet becomes one channel, view to understand :D 
    STARMake(meta, params.star_cpus)
    RSEM_Make(meta)
    bedFiles_Make(meta)
    genomeBPtype_Make(refFlat_geneData_Make.out)
    riboList_Make(refFlat_geneData_Make.out)
    transcriptData_Make(refFlat_geneData_Make.out)
    //ens_link_Make_for_SIMR(refFlat_geneData_Make.out)

    //Requested Processes: 
    faidx_Make(meta)
    bt2Index_Make(meta)
    bwaIndex_Make(meta)
    cellranger_singlecell_Index_Make(meta)
    igv_genome_Make(meta,riboList_Make.out)
}

workflow SECUNDO_INDEX_BUILD {
    take:
    meta

    main:
    // Processes required from running the Secundo pipeline: 
    refFlat_geneData_Make(meta) 
    STARMake(meta, params.star_cpus)
    RSEM_Make(meta)
    genomeBPtype_Make(refFlat_geneData_Make.out)
    riboList_Make(refFlat_geneData_Make.out)
    //ens_link_Make_for_SIMR(refFlat_geneData_Make.out)

}