#!/usr/bin/env nextflow

//AI
def helpMessage() {
    log.info"""
    Usage:
      nextflow run main.nf --input_dir <path> --gtf <file> --fasta <file>
    
    Required:
      --input_dir    Input data directory
      --gtf          GTF annotation file
      --fasta        Genome fasta file
    
    Optional:
      --rannot       Pre-built annotation (skips rannot generation)
      --outdir       Output directory (default: ${params.outdir})
    """.stripIndent()
}

workflow {
    // Validate required parameters
    if (!params.input_dir || !params.gtf || !params.fasta) {
        helpMessage()
        exit 1, "Error: Missing required parameters"
    }
    input_ch = channel.fromPath( "${params.input_dir}/*.bam_for_ORFquant", checkIfExists: true )
    gtf = file(params.gtf)
    fasta = file(params.fasta)
    // make it possible to provide existing annotation
    if (params.rannot) {
        rannot_ch = channel.fromPath( "${params.rannot}").first() //convert to value to avoid consuming channel
    } else {
        twobit_ch = UCSC_FATOTWOBIT(fasta)
        rannot_ch = ORFQUANT_ANNOTATION(gtf, twobit_ch, fasta)
    } 
    ORFQUANT(input_ch, rannot_ch, fasta)
    ORFQUANT_PLOTS(input_ch, ORFQUANT.out.orfquant_results, rannot_ch)
    ORFQUANT_REPORT(ORFQUANT_PLOTS.out.plots.collect())
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
    tag "ORFQUANT on $input_fororfquant.simpleName"
    publishDir params.outdir, mode: 'copy'
    label 'multi'
    errorStrategy 'ignore'

    // WARN: only works with given version, do not bump up!
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    ? 'https://depot.galaxyproject.org/singularity/orfquant:1.1.0--r40_1'
    : 'quay.io/biocontainers/orfquant:1.1.0--r40_1'}"

    input:
    path input_fororfquant
    path Rannot
    path fasta

    output:
    path "*_Detected_ORFs.gtf"
    path "*_Protein_sequences.fasta"
    path "*_final_ORFquant_results", emit: orfquant_results

    script:
    """
    #!/usr/bin/env Rscript

    library("ORFquant")

    suppressWarnings(run_ORFquant(
        for_ORFquant_file="${input_fororfquant}"
        , annotation_file="${Rannot}"
        , interactive=FALSE
        , n_cores=4
    ))
    """
}

process ORFQUANT_PLOTS{

    publishDir params.outdir, mode: 'copy'

    // WARN: only works with given version, do not bump up!
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    ? 'https://depot.galaxyproject.org/singularity/orfquant:1.1.0--r40_1'
    : 'quay.io/biocontainers/orfquant:1.1.0--r40_1'}"


    input:
    path input_fororfquant
    path orfquant_results
    path Rannot

    output:
    path "*plots"
    path "*plots/*_plots_RData", emit: plots

    script:
    """
    #!/usr/bin/env Rscript

    library("ORFquant")
    library("magrittr")

    for_orfquant_input_files <- "${input_fororfquant}" #%>% stringr::str_split(" ") %>% unlist()
    orfquant_input_files <- "${orfquant_results}" #%>% stringr::str_split(" ") %>% unlist()
    sample_names <- orfquant_input_files %>% sub(".bam_for_ORFquant_final_ORFquant_results","",.)

    # first, create plots
    plot_ORFquant_results(
            for_ORFquant_file=for_orfquant_input_files
            , ORFquant_output_file = orfquant_input_files
            , annotation_file="${Rannot}"
            #, output_plots_path = "plots"
            , prefix = sample_names
            )

    
    """
}

process ORFQUANT_REPORT {

    publishDir params.outdir, mode: 'copy'

    // WARN: only works with given version, do not bump up!
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    ? 'https://depot.galaxyproject.org/singularity/orfquant:1.1.0--r41h9ee0642_3'
    : 'quay.io/biocontainers/orfquant:1.1.0--r41h9ee0642_3'}"

    //containerOptions { "-v /usr/share/fonts:/usr/share/fonts:ro" }

    input:

    path Rplots

    output:

    path "ORFquant_report.html"

    script:
    """
    #!/usr/bin/env Rscript

    library("ORFquant")
    library("magrittr")

    # Convert input string to R vector (input = plot files *ORFquant_plots_RData)
    input_files <- "${Rplots}" %>% stringr::str_split(" ") %>% unlist()

    # Generate sample names from file names
    input_sample_names <- input_files %>% basename() %>% sub("_ORFquant_plots_RData","",.)

    # AI: Get absolute paths to ensure files can be found during RMD rendering
    input_files <- normalizePath(input_files, mustWork = TRUE)

    # AI: Get absolute path for output file in current directory
    output_file <- file.path(getwd(), "ORFquant_report.html")

    # AI: Check that files exist before rendering
    file_check <- sapply(input_files, file.exists)
    if(!all(file_check)) {
        stop("Missing input files: ", paste(input_files[!file_check], collapse=", "))
    }

    # Path to rmd template stored in container itself
    rmd_template_source <- paste(system.file(package="ORFquant"),"/rmd/ORFquant_template.Rmd",sep="")

    # AI: SINGULARITY FIX: Copy template to writable working directory
    # This avoids permission issues when rmarkdown tries to write intermediate files
    rmd_template_local <- file.path(getwd(), "ORFquant_template.Rmd")
    file.copy(rmd_template_source, rmd_template_local, overwrite = TRUE)
    # Verify template was copied
    if(!file.exists(rmd_template_local)) {
        stop("Failed to copy RMD template to working directory")
    }

    # Finally, generate html
    # The create_ORFquant_html_report function fails - use raw commands instead
    sink(file = paste(output_file,"_ORFquant_report_output.txt",sep = ""))
    # render RMarkdown file > html report:
    suppressWarnings(rmarkdown::render(rmd_template_local, 
                            params = list(input_files = input_files,
                                          input_sample_names = input_sample_names),
                            output_file = basename(output_file),
                            output_dir = dirname(output_file),
                            knit_root_dir = getwd(),
                            intermediates_dir = getwd()  # Explicitly set intermediate files directory for singularity to avoid permission issues
                            ))
    sink()
    """

}
