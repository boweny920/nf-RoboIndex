process refFlat_geneData_Make {

    input: 
    val meta

    output:
    tuple path("${meta.id}.${meta.annotation_version}.gene_data.txt"), path("${meta.id}.${meta.annotation_version}.refFlat.txt"), val(meta)
    // This is important for downstream processes to all work on the correct channel 

    publishDir "${meta.index_Folder}/annotation/${meta.annotation_version}/tables", mode: 'copy'

    script:     

    file("${meta.index_Folder}/annotation/${meta.annotation_version}/tables").mkdirs()     // tables folder
    file("${meta.index_Folder}/annotation/${meta.annotation_version}/gtfs").mkdirs()    // gtfs folder
    gtf_folder = file(meta.gtf_path).path

    if ( !meta.ensembl_db ) 
    """
    gtfToGenePred -genePredExt -geneNameAsName2 -ignoreGroupsWithoutExons ${meta.gtf_path} genePred.txt 

    awk 'BEGIN { OFS="\t" } {print \$12, \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10}' genePred.txt > ${meta.id}.${meta.annotation_version}.refFlat.txt

    copy_tool.py ${meta.gtf_path} ${meta.index_Folder}/annotation/${meta.annotation_version}/gtfs/${meta.id}.${meta.annotation_version}.gtf

    geneData_MAKE_Version_2.py -g ${meta.gtf_path}

    rename_tool.py ${file(meta.gtf_path).name.replaceAll('.gtf', '.gene_data.txt')} ${meta.id}.${meta.annotation_version}.gene_data.txt
    """

    else 
    """
    gtfToGenePred -genePredExt -geneNameAsName2 -ignoreGroupsWithoutExons ${meta.gtf_path} genePred.txt 

    awk 'BEGIN { OFS="\t" } {print \$12, \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10}' genePred.txt > ${meta.id}.${meta.annotation_version}.refFlat.txt

    copy_tool.py ${meta.gtf_path} ${meta.index_Folder}/annotation/${meta.annotation_version}/gtfs/${meta.id}.${meta.annotation_version}.gtf

    geneData_MAKE_Version_2.py -g ${meta.gtf_path} -e ${meta.ensembl_db}

    rename_tool.py ${file(meta.gtf_path).name.replaceAll('.gtf', '.gene_data.txt')} ${meta.id}.${meta.annotation_version}.gene_data.txt
    """
}

process STARMake {
    input:
    val meta 
    val cpus

    output:
    path "STAR_76_*"

    publishDir "${meta.index_Folder}/annotation/${meta.annotation_version}", mode: 'copy'

    script:
    """
    starVer=`STAR --version`

    STAR --runThreadN ${cpus} --runMode genomeGenerate --genomeDir STAR_76_\$starVer --genomeFastaFiles ${meta.fna_path} --sjdbGTFfile ${meta.gtf_path}  --sjdbOverhang 75 --limitGenomeGenerateRAM=200000000000

    ln -s ${meta.index_Folder}/annotation/${meta.annotation_version}/STAR_76_\$starVer ${meta.index_Folder}/annotation/${meta.annotation_version}/STAR_76bp
    """
}

process RSEM_Make {
    input: 
    val  meta

    output:
    path "RSEM/*"

    publishDir "${meta.index_Folder}/annotation/${meta.annotation_version}", mode: 'copy'

    script:
    """
    mkdir RSEM 
    rsem-prepare-reference --gtf ${meta.gtf_path} ${meta.fna_path} RSEM/${meta.id}.${meta.annotation_version}.RSEM
    """
}

process bedFiles_Make {
    input:
    val meta

    output:
    path "beds/*"

    publishDir "${meta.index_Folder}/annotation/${meta.annotation_version}" , mode: 'copy'

    script:
    """
    Bed_files_generate_Version_2.py -g ${meta.gtf_path} -o beds -d ${meta.id} -v ${meta.annotation_version}
    """
}

process genomeBPtype_Make {
    input:
    tuple path(gene_data_file), path(refFlat_file), val(meta)
    // path gene_data_file
    // path refFlat_file
    // val meta

    output:
    path "GenomeBP_TYPE/${refFlat_file.name.replaceAll('.refFlat.txt', '.GenomeBpTypes.txt')}"

    // Better to use cp than publisheDir here because of how the python script output is designed 
    script:

    file("${meta.index_Folder}/annotation/${meta.annotation_version}/extras" ).mkdirs()

    """
    Bed_files_GenomeTypeBP_file_MAKE_Version_6.py -r ${refFlat_file} -n ${meta.fna_path} -g ${gene_data_file} -o GenomeBP_TYPE

    cp GenomeBP_TYPE/${refFlat_file.name.replaceAll('.refFlat.txt', '.GenomeBpTypes.txt')} ${meta.index_Folder}/annotation/${meta.annotation_version}/extras/${meta.id}.${meta.annotation_version}.GenomeBpTypes.txt
    """
}

process riboList_Make {
    input:
    tuple path(gene_data_file), path(refFlat_file), val(meta)
    // path gene_data_file
    // val meta
    
    output:
    path "*.riboList.default.txt"

    publishDir "${meta.index_Folder}/annotation/${meta.annotation_version}/extras" , mode: 'copy'

    script:

    """
    java -Xmx150g -jar ${projectDir}/bin/picard.jar CreateSequenceDictionary R=${meta.fna_path} O=${file(meta.fna_path).name}.dict 

    RiboRNA_MAKE_Version_1.py -g ${gene_data_file} -p ${file(meta.fna_path).name}.dict

    copy_tool.py ${file(meta.fna_path).name}.dict ${meta.index_Folder}

    copy_tool.py ${file(meta.fna_path)} ${meta.index_Folder}
    """
    //Include the picard.jar into the singularity in the future? 
}

process transcriptData_Make {
    input:
    tuple path(gene_data_file), path(refFlat_file), val(meta)
    // path gene_data_file
    // path refFlat_file
    // val meta

    output:
    path "*.transcript_data.txt"

    publishDir "${meta.index_Folder}/annotation/${meta.annotation_version}/tables" , mode: 'copy'

    script:
    """
    transcriptData_MAKE_Version_1.py -g ${(meta.gtf_path)} -d ${gene_data_file} -r ${refFlat_file} 
    """
}

process ens_link_Make_for_SIMR {
    input:
    tuple path(gene_data_file), path(refFlat_file), val(meta)

    output:
    val "${meta.annotation_version}"

    script:

    if (!meta.annotation_version.contains('Ens_'))
    """
    ln -s ${meta.index_Folder}/annotation/${meta.annotation_version} ${meta.index_Folder}/annotation/${params.latest_ens}
    """
    else
    """
    echo "${params.latest_ens} soft link creat skipped"
    """
}