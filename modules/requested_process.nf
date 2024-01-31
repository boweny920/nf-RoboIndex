process faidx_Make {
    input: 
    val meta

    output:
    path "*.fai"

    publishDir "${meta.index_Folder}", mode: 'copy'

    script:
    def fna = file(meta.fna_path)
    """
    samtools faidx ${fna} 
    
    mv ${meta.fna_path}.fai ./
    """
}

process bt2Index_Make {
    input:
    val meta

    output:
    path "bowtie2/*" //The naming might be different if the fasta file does not end with .fa

    publishDir "${meta.index_Folder}", mode: 'copy'

    script:
    """
    mkdir -p bowtie2

    bowtie2-build ${meta.fna_path} bowtie2/${meta.id}
    """
}

process bwaIndex_Make {
    input:
    val meta

    output:
    path "BWA/*"

    publishDir "${meta.index_Folder}", mode: 'copy'

    script:

    def fna = file(meta.fna_path)

    """
    mkdir -p BWA

    bwa index ${fna}

    mv `ls ${fna}.*|grep -v ".fai"|grep -v ".dict"` BWA
    """
}

process cellranger_singlecell_Index_Make {
    input:
    val meta

    output:
    path "10x/*"

    publishDir "${meta.index_Folder}/annotation/${meta.annotation_version}/", mode: 'copy'

    script:
    """
    cellranger_ver=`cellranger --version |awk '{print \$2}'`

    mkdir -p 10x 

    cellranger mkref --genome=\$cellranger_ver  --fasta=${meta.fna_path} --genes=${meta.gtf_path} --memgb=300 

    ln -s \$cellranger_ver cellranger

    mv \$cellranger_ver cellranger 10x
    """
}

process igv_genome_Make {
    input: 
    val meta 
    val riboList_Make_out //tells the pipeline to wait until that process is finished

    publishDir "${meta.index_Folder}/annotation/${meta.annotation_version}/extras", mode: 'copy'
    
    output:
    path "${meta.id}.${meta.annotation_version}.igv_genome.json"

    script:
    """
    igv_genome_make_v2.py -g ${meta.index_Folder}/annotation/${meta.annotation_version}/gtfs/${meta.id}.${meta.annotation_version}.gtf -n ${meta.index_Folder}/${file(meta.fna_path).getName()} -f ${meta.index_Folder}/${file(meta.fna_path).getName()}.fai -i ${meta.id} -a ${meta.name} -v ${meta.annotation_version}
    """
}
