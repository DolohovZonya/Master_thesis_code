``` {r}
library(SVDFunctions)
library(dplyr)
```

# The code below has been performed for GBR sampleset

# Generate genotype matrix

``` {r}
variants <- read.table("one_consequence.tsv", header=TRUE, sep="\t")
```

1 476 344 rows annotation file from vep

``` {r}
vars <- paste(paste(paste("chr", variants$chr, sep = ""), variants$pos, sep = ":"), variants$ref, variants$alt, sep = "\t")
samples <- scan("intersected_samples.txt", what = character())
```

read 91 samples

``` {r}
result <- SVDFunctions::scanBinaryFile(“new_bin_bin”, “new_bin_meta”, samples = samples, vars = vars, DP = 0, GQ = 0)
```

# Generate merged file of annotaion and genotype distribution of samples

``` {r}
variants <- as.data.frame(cbind(vars, variants))
master_table <- inner_join(result, variants)
```

161558 rows

# Calculate the internal allele frequency

``` {r}
master_table$af <- (master_table$het + master_table$hom_alt*2) / 182 (n_samples*2)
```

# Calculate the variance between gnomAD and internal frequencies

``` {r}
master_table$var <- master_table$af - master_table$freq
```

# Allele frequency filtering

Firstly, we can apply basic filters, dividing our data by variant
annotation

## synonymous variants:

``` {r}
master_syns <- master_table %>% filter(reason == "synonymous_variant")
```

## missense+ptv variants:

``` {r}
master_mis_ptv <- master_table %>% filter(reason == "missense_variant" | reason == "frameshift_variant" | reason == "stop_gained" | reason == "stop_lost" | reason == "splice_acceptor_variant" | reason == "splice_donor_variant" | reason == "start_lost")
```

Then, we need to perform filtering and get variants we actually want to
work with. It will include:

1.  calculator variance between gnomAD allele frequency and internal
    allele frequency
2.  filtering by \|variance\| \<= 2\*std of this variance
3.  filtering by gnomAD allele frequency
4.  filtering by internal allele frequency

This should help in filtering discordant variants

## synonymous variants:

``` {r}
master_syns_1p_filt <- master_syns %>% filter((freq <= 0.01 & freq !=0) & af <= 0.01 & (var <= 2*sd(master_table$var) & var >= -2*sd(master_table$var)))
```

## missense+ptv variants:

``` {r}
master_mis_ptv_1p_filt <- master_mis_ptv %>% filter((freq <= 0.01 & freq !=0) & internal_af <= 0.01 & (var <= 2*sd(master_table$var) & var >= -2*sd(master_table$var))
```

Platform parameters: DP=0, GQ=0, gnomAD_MAF = 0.01, отдельно syn,
отдельно mis+ptv

# Group variants by gene

## synonymous variants:

``` {r}
master_1p_syns_gr <- master_syns_1p_filt %>%                                                                                                   group_by(gene) %>%
  summarise(
    hom_ref = sum(hom_ref, na.rm = TRUE),
    het = sum(het, na.rm = TRUE),
    hom_alt = sum(hom_alt, na.rm = TRUE),
    n_variants = sum(n_variants, na.rm = TRUE),
    call_rate = mean(call_rate, na.rm = TRUE)
  )
```

## missense+ptv variants:

``` {r}
master_1p_mis_ptv_gr <- master_mis_ptv_1p_filt %>%                                                                                                   group_by(gene) %>%
  summarise(
    hom_ref = sum(hom_ref, na.rm = TRUE),
    het = sum(het, na.rm = TRUE),
    hom_alt = sum(hom_alt, na.rm = TRUE),
    n_variants = sum(n_variants, na.rm = TRUE),
    call_rate = mean(call_rate, na.rm = TRUE)
  )
```

# Import files from platform

## synonymous variants:

``` {r}
syns_controls_1p <- read.table("syns_genes_controls_1p.tsv", header=TRUE, sep="\t")
```

## missense+ptv variants:

``` {r}
mis_ptv_controls_1p <- read.table("mis_ptv_genes_controls_1p.tsv", header=TRUE, sep="\t")
```

# Merge cases matrix and controls matrix

Cases samples get the flag '\_test'

## synonymous variants:

``` {r}
merged_data_syns_1p <- merge(master_1p_syns_gr, syns_controls_1p, all = TRUE, suffix = c('_test', '_control'), by = 'gene')

merged_data_syns_1p$hom_ref_control[is.na(merged_data_syns_1p$hom_ref_control)] <- merged_data_syns_1p$n_variants_test[is.na(merged_data_syns_1p$hom_ref_control)] * 496

merged_data_syns_1p$hom_ref_test[is.na(merged_data_syns_1p$hom_ref_test)] <- 91

merged_data_syns_1p[is.na(merged_data_syns_1p)] <- 0

merged_data_syns_1p$hom_ref_control <- merged_data_syns_1p$hom_ref_control + merged_data_syns_1p$hom_alt_control

merged_data_syns_1p$hom_alt_control <- 0
```

## missense+ptv variants:

``` {r}
merged_data_mis_ptv_1p <- merge(master_1p_mis_ptv_gr, mis_ptv_controls_1p, all = TRUE, suffix = c('_test', '_control'), by = 'gene')

merged_data_mis_ptv_1p$n_variants_test[is.na(merged_data_mis_ptv_1p$hom_ref_control)] * 496

merged_data_mis_ptv_1p$hom_ref_test[is.na(merged_data_mis_ptv_1p$hom_ref_test)] <- 91

merged_data_mis_ptv_1p[is.na(merged_data_mis_ptv_1p)] = 0

merged_data_mis_ptv_1p$hom_ref_control <- merged_data_mis_ptv_1p$hom_ref_control + merged_data_mis_ptv_1p$hom_alt_control

merged_data_mis_ptv_1p$hom_alt_control <- 0

merged_data_mis_ptv_1p[is.na(merged_data_mis_ptv_1p)] = 0
```

# Allele based fisher test

## synonymous variants:

``` {r}
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
```

## missense+ptv variants:

``` {r}
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
```

# Calculate p value distributions for visualization

## synonymous variants:

``` {r}
merged_data_syns_fisher$fisher_p_value[merged_data_syns_fisher$fisher_p_value == 0] <- .Machine$double.xmin

merged_data_syns_fisher$fisher_p_value <- as.numeric(merged_data_syns_fisher$fisher_p_value)

merged_data_syns_fisher$neg_log_p <- -log10(merged_data_syns_fisher$fisher_p_value)

merged_data_syns_fisher$theoretical <- -log10((rank(merged_data_syns_fisher$fisher_p_value) - 0.5) / length(merged_data_syns_fisher$fisher_p_value))
```

## missense+ptv variants:

``` {r}
merged_data_mis_ptv_fisher$fisher_p_value[merged_data_mis_ptv_fisher$fisher_p_value == 0] <- .Machine$double.xmin

merged_data_mis_ptv_fisher$fisher_p_value <- as.numeric(merged_data_mis_ptv_fisher$fisher_p_value)

merged_data_mis_ptv_fisher$neg_log_p <- -log10(merged_data_mis_ptv_fisher$fisher_p_value)

merged_data_mis_ptv_fisher$theoretical <- -log10((rank(merged_data_mis_ptv_fisher$fisher_p_value) - 0.5) / length(merged_data_mis_ptv_fisher$fisher_p_value))
```

# Save the result table

## synonymous variants:

``` {r}
write.table(merged_data_syns_fisher, "for_qpplot_carriers_syns_1p_new_p_value.tsv", sep="\t")
```

## missense+ptv variants:

``` {r}
write.table(merged_data_mis_ptv_fisher, "for_qpplot_carriers_mis_ptv_1p_new_p_value.tsv", sep="\t")
```

# Genotype based Fisher test

## synonymous variants:

``` {r}
fisher_p_values_syns <- numeric(nrow(merged_data_syns_1p))

odds_ratios_syns <- numeric(nrow(merged_data_syns_1p))

total_control <- 496

total_gbr <- 91

for (i in 1:nrow(merged_data_syns_1p)) {
  het_control <- merged_data_syns_1p$het_control[i]  het_test <- merged_data_syns_1p$het_test[i]
  remaining_control <- total_control - het_control remaining_gbr <- total_gbr - het_test if (any(c(het_control, het_test, remaining_control, remaining_gbr) < 0)) { next } 
  matrix_data <- matrix(c(het_control, remaining_control, het_test, remaining_gbr), nrow = 2) if (any(rowSums(matrix_data) == 0) || any(colSums(matrix_data) == 0)) { next } 
  fisher_test <- fisher.test(matrix_data) fisher_p_values_syns[i] <- fisher_test$p.value  odds_ratios_syns[i] <- ifelse(is.finite(fisher_test$estimate), fisher_test$estimate, NA)
  }

results_syns_genes_fisher <- data.frame( gene = merged_data_syns_1p$gene, fisher_p_value = fisher_p_values_syns, odds_ratio = odds_ratios_syns)
```

## missense+ptv variants:

``` {r}
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
```

# Calculate p value distributions for visualization

## synonymous variants:

``` {r}
results_syns_genes_fisher$fisher_p_value[results_syns_genes_fisher$fisher_p_value == 0] <- .Machine$double.xmin

results_syns_genes_fisher$fisher_p_value <- as.numeric(results_syns_genes_fisher$fisher_p_value)

results_syns_genes_fisher$neg_log_p <- -log10(results_syns_genes_fisher$fisher_p_value)

results_syns_genes_fisher$theoretical <- -log10((rank(results_syns_genes_fisher$fisher_p_value) - 0.5) / length(results_syns_genes_fisher$fisher_p_value))
```

## missense+ptv variants:

``` {r}
results_mis_ptv_genes_fisher$fisher_p_value[results_mis_ptv_genes_fisher$fisher_p_value == 0] <- .Machine$double.xmin

results_mis_ptv_genes_fisher$fisher_p_value <- as.numeric(results_mis_ptv_genes_fisher$fisher_p_value)

results_mis_ptv_genes_fisher$neg_log_p <- -log10(results_mis_ptv_genes_fisher$fisher_p_value)

results_mis_ptv_genes_fisher$theoretical <- -log10((rank(results_mis_ptv_genes_fisher$fisher_p_value) - 0.5) / length(results_mis_ptv_genes_fisher$fisher_p_value))
```

# Save the result table

## synonymous variants:

``` {r}
write.table(results_syns_genes_fisher, "for_qpplot_genes_syns_1p_new_p_value.tsv", sep="\t")
```

## missense+ptv variants:

``` {r}
write.table(results_mis_ptv_genes_fisher, "for_qpplot_genes_mis_ptv_1p_new_p_value.tsv", sep="\t")
```

# QQPlots

``` {r}
import pandas as pd 
import numpy as np 
import scipy.stats as stats 
import matplotlib.pyplot as plt 
```

``` {r}
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
```

# Calculate genomic inflation according to alpha

alpha is a value which is multiplying standard deviation to get a good
filtering threshold.

``` {r}
results <- data.frame(alpha = numeric(), lambda = numeric())

for (alpha in seq(0, 5, by = 0.1)) {
  master_syns_1p_filt <- master_syns %>% 
    filter((freq_afr <= 0.01 & freq_afr != 0) & 
           af <= 0.01 & 
           (var <= alpha*sd(master_table$var) & var >= -alpha*sd(master_table$var)))
  
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
  
  merged_data_syns_1p <- merge(master_1p_syns_gr, syns_controls_1p, all = TRUE, 
                              suffix = c('_test', '_control'), by = 'gene')
  
  merged_data_syns_1p$hom_ref_control[is.na(merged_data_syns_1p$hom_ref_control)] <-merged_data_syns_1p$n_variants_test[is.na(merged_data_syns_1p$hom_ref_control)] * 896
  merged_data_syns_1p$hom_ref_test[is.na(merged_data_syns_1p$hom_ref_test)] <- 61
  merged_data_syns_1p[is.na(merged_data_syns_1p)] <- 0
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
 p_values <- merged_data_syns_fisher$fisher_p_value
safe_p <- pmin(pmax(p_values, .Machine$double.eps), 1 - .Machine$double.eps)
chisq <- qchisq(1 - safe_p, df = 1)
current_lambda <- median(chisq,) / qchisq(0.5, df = 1)
results <- rbind(results, data.frame(alpha = alpha, lambda = current_lambda))
  
cat("Alpha:", alpha, "Lambda:", current_lambda, "\n")
}


best_alpha <- results$alpha[which.min(results$lambda)]
cat("\nBest alpha:", best_alpha, "with lambda:", min(results$lambda), "\n")


print(results)
```
