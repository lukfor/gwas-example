#!/usr/bin/env nextflow

// Pipeline parameters
params.genotypes = "data/genotypes.vcf.gz"
params.phenotypes = "data/phenotypes.txt"
params.pheno_name = "caffeine_consumption"
params.output = "output"

// Define processes
process QUALITY_CONTROL {

    input:
    path genotypes

    output:
    path "genotypes.filtered.vcf.gz", emit: filtered_genotypes

    script:
    """
    plink2 --vcf ${genotypes} \
           --maf 0.01 \
           --hwe 1e-6 \
           --mind 0.03 \
           --autosome \
           --keep-allele-order \
           --recode vcf bgz \
           --out genotypes.filtered
    """
}

process CALCULATE_PCA {
    publishDir "${params.output}", mode: 'copy'

    input:
    path filtered_genotypes

    output:
    path "pca.eigenvec", emit: pca_results

    script:
    """
    plink2 --vcf ${filtered_genotypes} \
           --double-id \
           --pca \
           --out pca
    """
}

process RUN_GWAS {
    publishDir "${params.output}", mode: 'copy'

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
           --out ${pheno_name}.gwas
    """
}

process PLOT_PCA {
    publishDir "${params.output}", mode: 'copy'

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

process PLOT_RESULTS {
    publishDir "${params.output}", mode: 'copy'

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
    genotypes_ch = Channel.fromPath(params.genotypes, checkIfExists: true)
    filtered_genotypes_ch = QUALITY_CONTROL(genotypes_ch)

    // Calculate PCA
    pca_ch = CALCULATE_PCA(filtered_genotypes_ch)

    phenotypes_ch = Channel.fromPath(params.phenotypes, checkIfExists: true)

    // Split the phenotype names and create a channel
    pheno_names_ch = Channel.of(params.pheno_name)
    
    // Combine the phenotype names with the other required channels
    combined_input_ch = filtered_genotypes_ch
        .combine(phenotypes_ch)
        .combine(pca_ch)
        .combine(pheno_names_ch)

    // Run GWAS for each combination
    gwas_ch = RUN_GWAS(combined_input_ch)

    PLOT_PCA(pca_ch)
    PLOT_RESULTS(gwas_ch)
}