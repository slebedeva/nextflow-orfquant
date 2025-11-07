#!/usr/bin/env nextflow

/*
 * The following pipeline parameters specify the reference genomes
 * and mapped reads and can be provided as command line options
 * Input is the directory where bam files are located and pipeline will take all bam files in that directory
 * Mapped reads for test data are downloaded from Google drive as described in Riboseqc manual
 * Reference genome and gtf are chr22 and chrM of the hg38 gencode assembly.
 */
params.input_dir = "$baseDir/test_data"
params.gtf = "$baseDir/test_data/test_human_chrM_22.gtf"
params.fasta = "$baseDir/test_data/test_human_chrM_22.fa"
params.rmd_template = "https://raw.githubusercontent.com/slebedeva/nextflow-riboseqc/refs/heads/main/riboseqc_template.Rmd"
params.outdir = "results"


workflow {
    input_ch = channel.fromPath( "${params.input_dir}/*.bam_for_ORFquant", checkIfExists: true )
    gtf = file(params.gtf)
    fasta = file(params.fasta)
    twobit_ch = UCSC_FATOTWOBIT(fasta)
    rannot_ch = ORFQUANT_ANNOTATION(gtf, twobit_ch, fasta)
    ORFQUANT(input_ch, rannot_ch, fasta)
    //ORFQUANT_REPORT(ORFQUANT.out.orfquant_results.collect(), params.rmd_template)
}

process UCSC_FATOTWOBIT {
    tag "${fasta}"

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    ? 'oras://community.wave.seqera.io/library/ucsc-fatotwobit:482--1d5005b012bd3271'
    : 'community.wave.seqera.io/library/ucsc-fatotwobit:482--f820aabce6f6870e'}"

    input:
        path fasta

    output:
        path "${twobit}"

    script:
    def extension = fasta.toString().tokenize('.')[-1]
    def name = fasta.toString() - ".${extension}"
    twobit = name + ".2bit"
    """
    faToTwoBit $fasta ${twobit}
    """

}

process ORFQUANT_ANNOTATION {
    tag "$gtf"
    publishDir params.outdir, mode: 'copy'


    // WARN: only works with given version, do not bump up!
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    ? 'https://depot.galaxyproject.org/singularity/orfquant:1.1.0--r40_1'
    : 'quay.io/biocontainers/orfquant:1.1.0--r40_1'}"

    input:
    path gtf
    path twobit
    path fasta

    output:
    path '*Rannot', emit: orfquant_annotation

    script:
    """
    #!/usr/bin/env Rscript

    library("ORFquant")
    library("magrittr")

    # prepare annotation for ORFquant
    # warning: ORFquant will not work with riboseqc anno

    # provide species and annotation name (defaults: Homo.sapiens and genc25)
    # extract them from Ensembl-style gtf file name
    gtf_file_name <- basename("${gtf}")
    ann_name <- sub(".gtf","",gtf_file_name) 
    species <- sub("_",".",sub("\\\\..+","",gtf_file_name))

    prepare_annotation_files(annotation_directory = ".",
                            twobit_file = "${twobit}"
                            ,gtf_file = "${gtf}"
                            ,scientific_name = species
                            ,annotation_name = ann_name
                            ,export_bed_tables_TxDb = T
                            ,forge_BSgenome = F
                            ,genome_seq="${fasta}")
    """
}

process ORFQUANT {
    tag "ORFQUANT on $input.simpleName"
    publishDir params.outdir, mode: 'copy'
    label 'multi'
    errorStrategy 'ignore'

    // WARN: only works with given version, do not bump up!
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    ? 'https://depot.galaxyproject.org/singularity/orfquant:1.1.0--r40_1'
    : 'quay.io/biocontainers/orfquant:1.1.0--r40_1'}"

    input:
    path input
    path Rannot
    path fasta

    output:
    path "*_Detected_ORFs.gtf"
    path "*_final_ORFquant_results", emit: orfquant_results

    script:
    """
    #!/usr/bin/env Rscript

    library("ORFquant")

    genome_dir <- "/projects/crc1678/genome_assemblies/human/hg38"
    ann_dir <- file.path(genome_dir, "orfquant-annotation")
    gtf_file_name <- "Homo_sapiens.GRCh38.112.gtf"
    suppressWarnings(run_ORFquant(
        for_ORFquant_file="${input}"
        , annotation_file="${Rannot}"
        , interactive=FALSE
        , n_cores=4
    ))
    """
}

