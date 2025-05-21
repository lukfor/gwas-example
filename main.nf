#!/usr/bin/env nextflow

// Pipeline parameters
params.genotypes = "data/files/*.vcf.gz"
params.phenotypes = "data/phenotypes.txt"
params.pheno_name = "caffeine_consumption,purple_hair"
params.outdir = "output"

// Define processes
process QUALITY_CONTROL {

    input:
    path genotypes

    output:
    path "${genotypes.baseName}.filtered.vcf.gz", emit: filtered_genotypes

    script:
    """
    plink2 --vcf ${genotypes} \
           --maf 0.01 \
           --hwe 1e-6 \
           --autosome \
           --keep-allele-order \
           --recode vcf bgz \
           --out ${genotypes.baseName}.filtered
    """
}

process CALCULATE_PCA {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path filtered_genotypes

    output:
    path "pca.eigenvec", emit: pca_results

    script:
    """

    bcftools concat ${filtered_genotypes} -Oz -o merged.vcf.gz

    plink2 --vcf merged.vcf.gz \
           --double-id \
           --pca \
           --out pca
    """
}

process RUN_GWAS {
    input:
    tuple path(filtered_genotypes), path(phenotypes), path(pca_results), val(pheno_name)

    output:
    tuple val(pheno_name), path("*.gwas.*.glm.*"), emit: gwas_results

    script:
    """
    plink2 --vcf ${filtered_genotypes} \
           --double-id \
           --pheno ${phenotypes} \
           --pheno-name ${pheno_name} \
           --covar ${pca_results} \
           --linear hide-covar \
           --out ${filtered_genotypes.baseName}.gwas
    """
}

process PLOT_PCA {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path result

    output:
    path "*.png", emit: gwas_results

    script:
    """
    #!Rscript
    library("ggplot2")
    pca_data <- read.table("${result}", header=TRUE, comment.char="")
    png("${result.baseName}.png", width = 1200, height = 800)
    ggplot(pca_data, aes(x = PC1, y = PC2)) + geom_point()
    dev.off()
    """
}

process MERGE_GWAS_RESULTS {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(pheno_name), path(results)

    output:
    tuple val(pheno_name), path("${pheno_name}.txt"), emit: merged_results

    script:
    """
    csvtk concat -C "" ${results} > ${pheno_name}.txt
    """
}


process PLOT_RESULTS {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(pheno_name), path(result)

    output:
    path "*.png", emit: gwas_results

    script:
    """
    #!Rscript
    library("qqman")
    gwas_result <- read.table("${result}", header=TRUE, comment.char="")
    png("${result.baseName}.manhattan.png", width = 1200, height = 800)
    manhattan(x = gwas_result, chr = "X.CHROM", bp = "POS", p = "P", snp="ID")
    dev.off()
    png("${result.baseName}.qq.png", width = 1200, height = 800)
    qq(gwas_result\$P)
    dev.off()    
    """
}



// Define workflow
workflow {
    // Quality control
    genotypes_ch = channel.fromPath(params.genotypes, checkIfExists: true)
    filtered_genotypes_ch = QUALITY_CONTROL(genotypes_ch)

    // Calculate PCA
    pca_ch = CALCULATE_PCA(filtered_genotypes_ch.collect())
    PLOT_PCA(pca_ch)

    phenotypes_ch = channel.fromPath(params.phenotypes, checkIfExists: true)

    // Split the phenotype names and create a channel
    pheno_names_ch = channel.of(params.pheno_name.split(","))

    // Combine the phenotype names with the other required channels
    combined_input_ch = filtered_genotypes_ch
        .combine(phenotypes_ch)
        .combine(pca_ch)
        .combine(pheno_names_ch)

    // Run GWAS for each combination
    gwas_ch = RUN_GWAS(combined_input_ch)
    merged_gwas_ch = MERGE_GWAS_RESULTS(gwas_ch.groupTuple())
    
    PLOT_RESULTS(merged_gwas_ch)
}