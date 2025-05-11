#!/usr/bin/env nextflow

// Pipeline parameters
params.genotypes = "data/genotypes.vcf.gz"
params.phenotypes = "data/phenotypes.txt"
params.pheno_name = "caffeine_consumption"
params.outdir = "output"

// Define processes
process QUALITY_CONTROL {
    conda 'environment.yml'

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
    conda 'environment.yml'
    publishDir "${params.outdir}", mode: 'copy'

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
    conda 'environment.yml'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path filtered_genotypes
    path phenotypes
    path pca_results

    output:
    path "gwas.*.linear", emit: gwas_results

    script:
    """
    plink2 --vcf ${filtered_genotypes} \
           --double-id \
           --pheno ${phenotypes} \
           --pheno-name ${params$pheno_name} \
           --covar ${pca_results} \
           --linear hide-covar \
           --out gwas
    """
}

process PLOT_PCA {
    conda 'environment.yml'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path result

    output:
    path "*.png", emit: gwas_results

    script:
    """
    #!RScript
    library("ggplot2")
    pca_data <- read.table("${result}", header=TRUE, comment.char="")
    png("${result.baseName}.png", width = 1200, height = 800)
    ggplot(pca_data, aes(x = PC1, y = PC2)) + geom_point()
    dev.off()
    """
}

process PLOT_RESULTS {
    conda 'environment.yml'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path result

    output:
    path "*.png", emit: gwas_results

    script:
    """
    #!RScript
    library("qqman")
    gwas_result <- read.table("${result}", header=TRUE, comment.char="")
    png("${result.baseName}.manhattan.png", width = 1200, height = 800)
    manhattan(x = gwas_result, chr = "X.CHROM", bp = "POS", p = "P", snp="ID")
    dev.off()
    png("${result.baseName}.qq.png", width = 1200, height = 800)
    qq(gwas_result\$P)
    dev.off()    
    #TODO: qq-plot
    """
}


// Define workflow
workflow {
    // Quality control
    genotypes_ch = channel.fromPath(params.genotypes, checkIfExists: true)
    filtered_genotypes_ch = QUALITY_CONTROL(genotypes_ch)

    // Calculate PCA
    pca_ch = CALCULATE_PCA(filtered_genotypes_ch)

    // Run GWAS
    phenotypes_ch = channel.fromPath(params.phenotypes, checkIfExists: true)
    gwas_ch = RUN_GWAS(filtered_genotypes_ch, phenotypes_ch, pca_ch)

    PLOT_PCA(pca_ch)
    PLOT_RESULTS(gwas_ch)
}