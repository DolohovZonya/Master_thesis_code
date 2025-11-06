# Install everything (optional)

To work with the following pipeline, you will need some specific tools -
**bcftools**, **plink1.9**, **VEP**, **R version &gt;=3.1.0 and &lt;=
4.5.0**.

If you already have everything needed or you know how to install it -
you can skip this step. If you require some assistance - please use the
commands below.

## Install dependencies

    sudo apt install bcftools
    sudo apt install plink1.9

    git clone https://github.com/Ensembl/ensembl-vep.git
    cd ensembl-vep
    perl INSTALL.pl

### Install R

    wget https://cran.r-project.org/src/base/R-4/R-4.2.0.tar.gz
    tar -xzf R-4.2.0.tar.gz
    cd R-4.2.0
    sudo apt install gfortran
    sudo apt install libreadline-dev
    sudo apt install xorg-dev libxt-dev
    sudo apt install liblzma-dev
    sudo apt install libpcre2-dev
    ./configure --prefix=/usr/local --enable-R-shlib
    make
    sudo make install

On step 9 you might receive this error:

    configure: error: libcurl >= 7.28.0 library and headers are required with support for https

To fix it, install a compatible version of curl:

    sudo apt remove --purge curl libcurl4 libcurl4-openssl-dev

    wget https://curl.se/download/curl-7.88.1.tar.gz
    tar -xzf curl-7.88.1.tar.gz
    cd curl-7.88.1
    sudo apt install libssl-dev
    sudo apt install default-jdk
    sudo ldconfig /usr/local/lib
    ./configure --prefix=/usr/local
    make -j$(nproc)
    sudo make install

After that, repeat the R installation from step 9.

When successful, you should see something like:

    R version 4.2.0 (2022-04-22) -- "Vigorous Calisthenics"

Then install the SVDFunctions package:

    devtools::install_github("https://github.com/alexloboda/SVDFunctions")

# The code below has been performed for GBR sampleset

## Prepare the data

Data from the 1k Genomes consortium can be downloaded from:

<https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/vcf_with_sample_level_annotation/>

Remember to download also the `.tbi` index files.

Merge chromosome files (avoid mito, X, Y) and keep only exonic regions:

    for i in {1..22}; do
      bcftools sort ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5_extra_anno.20130502.genotypes.vcf.gz -Oz -o chr${i}_sorted.vcf.gz
    done

    cat chr{1..22}_sorted.vcf.gz > genomes_merged.vcf.gz

    bedtools intersect -a genomes_merged.vcf.gz -b 20130108.exome.targets.bed -header -wa | bgzip > exomes_merged.vcf.gz

## Relatedness analysis

### MAF filtering

    bcftools view -i 'MAF > 0.01' exomes_merged.vcf.gz -Oz -o exomes_maf_001.vcf.gz

### LD pruning

    plink1.9 --vcf exomes_maf_001.vcf.gz --indep-pairwise 50 5 0.2 --out exomes_maf_pruned
    plink --vcf exomes_maf_001.vcf.gz --extract exomes_maf_pruned.prune.in --recode vcf bgz --out exomes_maf_pruned

### IBD calculation and relatives exclusion

    plink1.9 --vcf exomes_maf_pruned.vcf.gz --genome --out ibd_results
    awk '$10 > 0.2 {print $1, $2}' ibd_results.genome | uniq > samples_to_remove.txt
    plink1.9 --vcf exomes_maf_pruned.vcf.gz --remove samples_to_remove.txt --recode vcf bgz --out exomes_maf_pruned_no_related

As a result, you will get the file
`exomes_maf_pruned_no_related.vcf.gz`. You can use it as a source for
generating the genotype matrix through the pipeline on
[DNAScore](https://dnascore.net/).

# Variant annotation

Split multiallelic variants and convert into “sites-only” format:

    bcftools norm -m -any exomes_merged.vcf.gz -Oz -o exomes_biallelic.vcf.gz
    zcat exomes_biallelic.vcf.gz | awk 'BEGIN {OFS="\t"} !/^#/ {print $1, $2, $3, $4}' > exomes_biallelic_sites_only.vcf

Annotation with Ensembl VEP (GRCh37). Activate: - gnomAD (exomes) allele
frequencies - canonical transcripts

## Parsing VEP output

Use this parser:
<https://github.com/alexloboda/vep-parser/blob/main/vep_parser.py>

Modify to include allele frequency column if needed.

# Generate binary file

    scanVCF(
      exomes_biallelic.vcf.gz,
      DP = 10L,
      GQ = 20L,
      samples = "your_samples.txt",
      bannedPositions = NULL,
      variants = NULL,
      returnGenotypeMatrix = TRUE,
      predictMissing = FALSE,
      missingRateThreshold = 0.1,
      regions = NULL,
      binaryPathPrefix = "binary_file",
      verbose = FALSE
    )

# Generate genotype matrix

    variants <- read.table("your_vep_output.tsv", header=TRUE, sep="\t") 
    vars <- paste(paste(paste("chr", variants$chr, sep=""), variants$pos, sep=":"), variants$ref, variants$alt, sep="\t")
    samples <- scan("your_gbr_sample_ids.txt", what=character())
    result <- SVDFunctions::scanBinaryFile("bin_file_name", "meta_file_name", samples=samples, vars=vars, DP=10, GQ=20)

# Merge annotation and genotype data

    variants <- as.data.frame(cbind(vars, variants))
    master_table <- inner_join(result, variants)
    master_table$af <- (master_table$het + master_table$hom_alt*2) / (master_table$hom_ref*2 + master_table$het*2 + master_table$hom_alt*2)
    master_table$var <- master_table$af - master_table$freq

# Allele frequency filtering

    master_syns <- master_table %>% filter(reason == "synonymous_variant")
    master_mis_ptv <- master_table %>% filter(reason %in% c("missense_variant", "frameshift_variant", "stop_gained", "stop_lost", "splice_acceptor_variant", "splice_donor_variant", "start_lost"))

# Then, we need to perform filtering and get variants we actually want to work with. It will include:

calculator variance between gnomAD allele frequency and internal allele
frequency filtering by |variance| &lt;= 2\*std of this variance
filtering by gnomAD allele frequency filtering by internal allele
frequency

# This should help in filtering discordant variants

# synonymous variants:

    master_syns_1p_filt <- master_syns %>% filter((freq <= 0.01 & freq !=0) & af <= 0.01 & (var <= 2*sd(master_table$var) & var >= -2*sd(master_table$var)))

# missense+ptv variants:

    master_mis_ptv_1p_filt <- master_mis_ptv %>% filter((freq <= 0.01 & freq !=0) & internal_af <= 0.01 & (var <= 2*sd(master_table$var) & var >= -2*sd(master_table$var))

# Platform parameters: DP=10, GQ=20, gnomAD\_MAF = 0.01, divide to syn and mis+ptv separate analyses

# Group variants by gene

# synonymous variants:

    master_1p_syns_gr <- master_syns_1p_filt %>%                                                                                                   group_by(gene) %>%
      summarise(
        hom_ref = sum(hom_ref, na.rm = TRUE),
        het = sum(het, na.rm = TRUE),
        hom_alt = sum(hom_alt, na.rm = TRUE),
        n_variants = sum(n_variants, na.rm = TRUE),
        call_rate = mean(call_rate, na.rm = TRUE)
      )

# missense+ptv variants:

    master_1p_mis_ptv_gr <- master_mis_ptv_1p_filt %>%                                                                                                   group_by(gene) %>%
      summarise(
        hom_ref = sum(hom_ref, na.rm = TRUE),
        het = sum(het, na.rm = TRUE),
        hom_alt = sum(hom_alt, na.rm = TRUE),
        n_variants = sum(n_variants, na.rm = TRUE),
        call_rate = mean(call_rate, na.rm = TRUE)
      )

# Import files from platform

# synonymous variants:

    syns_controls_1p <- read.table("syns_genes_controls_1p.tsv", header=TRUE, sep="\t")

# missense+ptv variants:

    mis_ptv_controls_1p <- read.table("mis_ptv_genes_controls_1p.tsv", header=TRUE, sep="\t")

# Merge cases matrix and controls matrix

# Cases samples get the flag ‘\_test’

# synonymous variants:

    merged_data_syns_1p <- merge(master_1p_syns_gr, syns_controls_1p, all = TRUE, suffix = c('_test', '_control'), by = 'gene')

    merged_data_syns_1p$hom_ref_control[is.na(merged_data_syns_1p$hom_ref_control)] <- merged_data_syns_1p$n_variants_test[is.na(merged_data_syns_1p$hom_ref_control)] * 496

    merged_data_syns_1p$hom_ref_test[is.na(merged_data_syns_1p$hom_ref_test)] <- 91
    merged_data_syns_1p[is.na(merged_data_syns_1p)] <- 0

    merged_data_syns_1p$hom_ref_control <- merged_data_syns_1p$hom_ref_control + merged_data_syns_1p$hom_alt_control

    merged_data_syns_1p$hom_alt_control <- 0

# missense+ptv variants:

    merged_data_mis_ptv_1p <- merge(master_1p_mis_ptv_gr, mis_ptv_controls_1p, all = TRUE, suffix = c('_test', '_control'), by = 'gene')

    merged_data_mis_ptv_1p$n_variants_test[is.na(merged_data_mis_ptv_1p$hom_ref_control)] * 496

    merged_data_mis_ptv_1p$hom_ref_test[is.na(merged_data_mis_ptv_1p$hom_ref_test)] <- 91

    merged_data_mis_ptv_1p[is.na(merged_data_mis_ptv_1p)] = 0

    merged_data_mis_ptv_1p$hom_ref_control <- merged_data_mis_ptv_1p$hom_ref_control + merged_data_mis_ptv_1p$hom_alt_control

    merged_data_mis_ptv_1p$hom_alt_control <- 0

    merged_data_mis_ptv_1p[is.na(merged_data_mis_ptv_1p)] = 0

# Allele based fisher test

# synonymous variants:

    merged_data_syns_fisher <- merged_data_syns_1p %>%
      rowwise() %>%
      mutate(
        fisher_p_value = fisher.test(matrix(c(
          (het_control + hom_alt_control * 2),
          (hom_ref_control * 2 + het_control),
          (het_test + hom_alt_test * 2),
          (hom_ref_test * 2 + het_test)
        ), nrow = 2, byrow = TRUE))$p.value,
        
        odds_ratio = fisher.test(matrix(c(
          (het_control + hom_alt_control * 2),
          (hom_ref_control * 2 + het_control),
          (het_test + hom_alt_test * 2),
          (hom_ref_test * 2 + het_test)
        ), nrow = 2, byrow = TRUE))$estimate)

# missense+ptv variants:

    merged_data_mis_ptv_fisher <- merged_data_mis_ptv_1p %>%
      rowwise() %>%
      mutate(
        fisher_p_value = fisher.test(matrix(c(
          (het_control + hom_alt_control * 2),
          (hom_ref_control * 2 + het_control),
          (het_test + hom_alt_test * 2),
          (hom_ref_test * 2 + het_test)
        ), nrow = 2, byrow = TRUE))$p.value,
        
        odds_ratio = fisher.test(matrix(c(
          (het_control + hom_alt_control * 2),
          (hom_ref_control * 2 + het_control),
          (het_test + hom_alt_test * 2),
          (hom_ref_test * 2 + het_test)
        ), nrow = 2, byrow = TRUE))$estimate
      )

# Calculate p value distributions for visualization

# synonymous variants:

    merged_data_syns_fisher$fisher_p_value[merged_data_syns_fisher$fisher_p_value == 0] <- .Machine$double.xmin

    merged_data_syns_fisher$fisher_p_value <- as.numeric(merged_data_syns_fisher$fisher_p_value)

    merged_data_syns_fisher$neg_log_p <- -log10(merged_data_syns_fisher$fisher_p_value)

    merged_data_syns_fisher$theoretical <- -log10((rank(merged_data_syns_fisher$fisher_p_value) - 0.5) / length(merged_data_syns_fisher$fisher_p_value))

# missense+ptv variants:

    merged_data_mis_ptv_fisher$fisher_p_value[merged_data_mis_ptv_fisher$fisher_p_value == 0] <- .Machine$double.xmin

    merged_data_mis_ptv_fisher$fisher_p_value <- as.numeric(merged_data_mis_ptv_fisher$fisher_p_value)

    merged_data_mis_ptv_fisher$neg_log_p <- -log10(merged_data_mis_ptv_fisher$fisher_p_value)
    merged_data_mis_ptv_fisher$theoretical <- -log10((rank(merged_data_mis_ptv_fisher$fisher_p_value) - 0.5) / length(merged_data_mis_ptv_fisher$fisher_p_value))

# Save the result table

# synonymous variants:

    write.table(merged_data_syns_fisher, "for_qpplot_carriers_syns_1p_new_p_value.tsv", sep="\t")

# missense+ptv variants:

    write.table(merged_data_mis_ptv_fisher, "for_qpplot_carriers_mis_ptv_1p_new_p_value.tsv", sep="\t")

# Genotype based Fisher test

# synonymous variants:

    fisher_p_values_syns <- numeric(nrow(merged_data_syns_1p))

    odds_ratios_syns <- numeric(nrow(merged_data_syns_1p))

    total_control <- 496

    total_gbr <- 91

    for (i in 1:nrow(merged_data_syns_1p)) {
      het_control <- merged_data_syns_1p$het_control[i]
      het_test <- merged_data_syns_1p$het_test[i]
      
      remaining_control <- total_control - het_control
      remaining_gbr <- total_gbr - het_test
        if (any(c(het_control, het_test, remaining_control, remaining_gbr) < 0)) {
        next
      }
        matrix_data <- matrix(c(het_control, remaining_control, het_test, remaining_gbr), nrow = 2)
        if (any(rowSums(matrix_data) == 0) || any(colSums(matrix_data) == 0)) {
        next
      }
      fisher_test <- fisher.test(matrix_data)
      fisher_p_values_syns[i] <- fisher_test$p.value
      odds_ratios_syns[i] <- ifelse(is.finite(fisher_test$estimate), fisher_test$estimate, NA)
    }

    results_syns_genes_fisher <- data.frame(
      gene = merged_data_syns_1p$gene,  fisher_p_value = fisher_p_values_syns,
      odds_ratio = odds_ratios_syns)

# missense+ptv variants:

    fisher_p_values_mis_ptv <- numeric(nrow(merged_data_mis_ptv_1p))

    odds_ratios_mis_ptv <- numeric(nrow(merged_data_mis_ptv_1p))

    total_control <- 496

    total_gbr <- 91

    for (i in 1:nrow(merged_data_syns_1p)) {
      het_control <- merged_data_syns_1p$het_control[i]
      het_test <- merged_data_syns_1p$het_test[i]
      
      remaining_control <- total_control - het_control
      remaining_gbr <- total_gbr - het_test
        if (any(c(het_control, het_test, remaining_control, remaining_gbr) < 0)) {
        next
      }
        matrix_data <- matrix(c(het_control, remaining_control, het_test, remaining_gbr), nrow = 2)
        if (any(rowSums(matrix_data) == 0) || any(colSums(matrix_data) == 0)) {
        next
      }
      fisher_test <- fisher.test(matrix_data)
      fisher_p_values_syns[i] <- fisher_test$p.value
      odds_ratios_syns[i] <- ifelse(is.finite(fisher_test$estimate), fisher_test$estimate, NA)
    }

    results_mis_ptv_genes_fisher <- data.frame(
      gene = merged_data_mis_ptv_1p$gene,  fisher_p_value = fisher_p_values_mis_ptv,
      odds_ratio = odds_ratios_mis_ptv)

# Calculate p value distributions for visualization

# synonymous variants:

    results_syns_genes_fisher$fisher_p_value[results_syns_genes_fisher$fisher_p_value == 0] <- .Machine$double.xmin

    results_syns_genes_fisher$fisher_p_value <- as.numeric(results_syns_genes_fisher$fisher_p_value)

    results_syns_genes_fisher$neg_log_p <- -log10(results_syns_genes_fisher$fisher_p_value)

    results_syns_genes_fisher$theoretical <- -log10((rank(results_syns_genes_fisher$fisher_p_value) - 0.5) / length(results_syns_genes_fisher$fisher_p_value))

# missense+ptv variants:

    results_mis_ptv_genes_fisher$fisher_p_value[results_mis_ptv_genes_fisher$fisher_p_value == 0] <- .Machine$double.xmin

    results_mis_ptv_genes_fisher$fisher_p_value <- as.numeric(results_mis_ptv_genes_fisher$fisher_p_value)

    results_mis_ptv_genes_fisher$neg_log_p <- -log10(results_mis_ptv_genes_fisher$fisher_p_value)

    results_mis_ptv_genes_fisher$theoretical <- -log10((rank(results_mis_ptv_genes_fisher$fisher_p_value) - 0.5) / length(results_mis_ptv_genes_fisher$fisher_p_value))

# Save the result table

# synonymous variants:

    write.table(results_syns_genes_fisher, "for_qpplot_genes_syns_1p_new_p_value.tsv", sep="\t")

# missense+ptv variants:

    write.table(results_mis_ptv_genes_fisher, "for_qpplot_genes_mis_ptv_1p_new_p_value.tsv", sep="\t")

# QQPlots

    import pandas as pd 
    import numpy as np 
    import scipy.stats as stats 
    import matplotlib.pyplot as plt 

    results_df = pd.read_csv("input table" ,sep = "\t") 
    print(results_df) 
    expected = results_df['theoretical'] 
    observed = results_df['neg_log_p']
    plt.figure(figsize=(6, 6)) 
    plt.scatter(expected, observed, edgecolor='black') plt.plot([expected.min(), expected.max()], [expected.min(), expected.max()], color='red', linestyle='--') 
    plt.xlabel("Expected -log10(p-value)") 
    plt.ylabel("Observed -log10(p-value)") 
    plt.title("QQ-plot of p-values") 
    plt.show()

# Calculate genomic inflation according to alpha

alpha is a value which is multiplying standard deviation to get a good
filtering threshold. This parameter is introduced to balance the data in
conditions of hidden controls and different exome coverage qualities.

in the script below you can set alpha and estimate the lambda for your
resulted test

    for (alpha in seq(0.5, 1, by = 0.1)) {
      master_syns_1p_filt <- master_syns %>% 
        filter((freq <= 0.01 & freq != 0) & 
               af <= 0.01 & (var <= alpha*sd(master_table$var) & var >= -alpha*sd(master_table$var)))
      
      master_1p_syns_gr <- master_syns_1p_filt %>%
        group_by(gene) %>%
        summarise(
          hom_ref = sum(hom_ref, na.rm = TRUE),
          het = sum(het, na.rm = TRUE),
          hom_alt = sum(hom_alt, na.rm = TRUE),
          n_variants = sum(n_variants, na.rm = TRUE),
          call_rate = mean(call_rate, na.rm = TRUE)
        )

      master_1p_syns_gr <- as.data.frame(master_1p_syns_gr)
      
      merged_data_syns_1p <- merge(master_1p_syns_gr, syns_controls_1p, all=TRUE, suffix = c('_test', '_control'), by = 'gene')

    merged_data_syns_1p$hom_ref_control[is.na(merged_data_syns_1p$hom_ref_control)] <-merged_data_syns_1p$n_variants_test[is.na(merged_data_syns_1p$hom_ref_control)] * 1696
      merged_data_syns_1p$hom_ref_test[is.na(merged_data_syns_1p$hom_ref_test)] <- 87
      merged_data_syns_1p$hom_ref_control <- merged_data_syns_1p$hom_ref_control + merged_data_syns_1p$hom_alt_control
      merged_data_syns_1p$hom_alt_control <- 0

      merged_data_syns_fisher <- merged_data_syns_1p %>%
        rowwise() %>%
        mutate(
          fisher_p_value = fisher.test(matrix(c(
            (het_control + hom_alt_control * 2),
            (hom_ref_control * 2 + het_control),
            (het_test + hom_alt_test * 2),
            (hom_ref_test * 2 + het_test)
          ), nrow = 2, byrow = TRUE))$p.value,
          odds_ratio = fisher.test(matrix(c(
            (het_control + hom_alt_control * 2),
            (hom_ref_control * 2 + het_control),
            (het_test + hom_alt_test * 2),
            (hom_ref_test * 2 + het_test)
          ), nrow = 2, byrow = TRUE))$estimate)
      
    merged_data_syns_fisher <- as.data.frame(merged_data_syns_fisher) 
    merged_data_syns_fisher$fisher_p_value[merged_data_syns_fisher$fisher_p_value == 0] <- .Machine$double.xmin
    merged_data_syns_fisher$fisher_p_value <- as.numeric(merged_data_syns_fisher$fisher_p_value)

    merged_data_syns_fisher$neg_log_p <- -log10(merged_data_syns_fisher$fisher_p_value)

    merged_data_syns_fisher$theoretical <- -log10((rank(merged_data_syns_fisher$fisher_p_value) - 0.5) / length(merged_data_syns_fisher$fisher_p_value))
     p_values <- merged_data_syns_fisher$fisher_p_value
    safe_p <- pmin(pmax(p_values, .Machine$double.eps), 1 - .Machine$double.eps)
    chisq <- qchisq(1 - safe_p, df = 1)
    current_lambda <- median(chisq,) / qchisq(0.5, df = 1)
    filename <- paste0("1kg_gbr_syns_fisher_all_1p1p_", alpha, ".tsv")
    write.table(merged_data_syns_fisher, filename, sep = "\t")
    cat("Alpha:", alpha, "Lambda:", current_lambda, "\n")
    }

    best_alpha <- results$alpha[which.min(results$lambda)]
    cat("\nBest alpha:", best_alpha, "with lambda:", min(results$lambda), "\n")
    print(results)
