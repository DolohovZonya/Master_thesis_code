\# library(SVDFunctions)

\# library(dplyr)

\# The code below has been performed for GBR sampleset

### Generate genotype matrix

variants <- read.table("one\_consequence.tsv", header=TRUE, sep="\\t")

1 476 344 rows #annotation file from vep

vars <- paste(paste(paste("chr", variants$chr, sep = ""), variants$pos, sep = ":"), variants$ref, variants$alt, sep = "\\t")

samples <- scan("intersected\_samples.txt", what = character())

read 91 samples #read samples

result <- SVDFunctions::scanBinaryFile(“new\_bin\_bin”, “new\_bin\_meta”, samples = samples, vars = vars, DP = 0, GQ = 0)

### Generate merged file of annotaion and genotype distribution of samples

variants <- as.data.frame(cbind(vars, variants))

master\_table <- inner\_join(result, variants)

\# 161558 rows

### Calculate the internal allele frequency

master\_table$internal\_af <- (master\_table$het + master\_table$hom\_alt\*2) / 182 (n\_samples\*2)

Calculate the variance between gnomAD and internal frequencies

master\_table$var <- master\_table$internal\_af - master\_table$freq

### Allele frequency filtering

\# Firstly, we can apply basic filters, dividing our data by variant annotation

\# synonymous variants:

master\_syns <- master\_table %>% filter(reason == “synonymous\_variant”)

\# missense+ptv variants:

master\_mis\_ptv <- master\_table %>% filter(reason == "missense\_variant" | reason == "frameshift\_variant" | reason == "stop\_gained" | reason == "stop\_lost" | reason == "splice\_acceptor\_variant" | reason == "splice\_donor\_variant" | reason == "start\_lost")

Then, we need to perform filtering and get variants we actually want to work with. It will include:

1.  calculator variance between gnomAD allele frequency and internal allele frequency
2.  filtering by |variance| <= 2\*std of this variance
3.  filtering by gnomAD allele frequency
4.  filtering by internal allele frequency

\# This should help in filtering discordant variants

\# synonymous variants:

master\_syns\_1p\_filt <- master\_syns %>% filter((freq <= 0.01 & freq !=0) & internal\_af <= 0.01 & (var <= 2\*sd(master\_table$var) & var >= -2\*sd(master\_table$var))

\# missense+ptv variants:

master\_mis\_ptv\_1p\_filt <- master\_mis\_ptv %>% filter((freq <= 0.01 & freq !=0) & internal\_af <= 0.01 & (var <= 2\*sd(master\_table$var) & var >= -2\*sd(master\_table$var))

\№ Platform parameters: DP=0, GQ=0, gnomAD\_MAF = 0.01, отдельно syn, отдельно mis+ptv

### Group variants by gene

\# synonymous variants:

master\_1p\_syns\_gr <- master\_syns\_1p\_filt %>%                                                                                                   group\_by(gene) %>%

  summarise(

    hom\_ref = sum(hom\_ref, na.rm = TRUE),

    het = sum(het, na.rm = TRUE),

    hom\_alt = sum(hom\_alt, na.rm = TRUE),

    n\_variants = sum(n\_variants, na.rm = TRUE),

    call\_rate = mean(call\_rate, na.rm = TRUE)

  )

\# missense+ptv variants:

master\_1p\_mis\_ptv\_gr <- master\_mis\_ptv\_1p\_filt %>%                                                                                                   group\_by(gene) %>%

  summarise(

    hom\_ref = sum(hom\_ref, na.rm = TRUE),

    het = sum(het, na.rm = TRUE),

    hom\_alt = sum(hom\_alt, na.rm = TRUE),

    n\_variants = sum(n\_variants, na.rm = TRUE),

    call\_rate = mean(call\_rate, na.rm = TRUE)

  )

### Import files from platform

\# synonymous variants:

syns\_controls\_1p <- read.table("syns\_genes\_controls\_1p.tsv", header=TRUE, sep="\\t")

\# missense+ptv variants:

mis\_ptv\_controls\_1p <- read.table("mis\_ptv\_genes\_controls\_1p.tsv", header=TRUE, sep="\\t")

### Merge cases matrix and controls matrix

\# Cases samples get the flag ‘\_test’

\# synonymous variants:

merged\_data\_syns\_1p <- merge(master\_1p\_syns\_gr, syns\_controls\_1p, all = TRUE, suffix = c('\_test', '\_control'), by = 'gene')

merged\_data\_syns\_1p$hom\_ref\_control\[is.na(merged\_data\_syns\_1p$hom\_ref\_control)\] <- merged\_data\_syns\_1p$n\_variants\_test\[is.na(merged\_data\_syns\_1p$hom\_ref\_control)\] \* 496

fisher.test(matrix(hom\_ref\*2 + het, het + 2\* hom\_alt,)

\# missense+ptv variants:

merged\_data\_mis\_ptv\_1p <- merge(master\_1p\_mis\_ptv\_gr, mis\_ptv\_controls\_1p, all = TRUE, suffix = c('\_test', '\_control'), by = 'gene')

merged\_data\_mis\_ptv\_1p\[is.na(merged\_data\_mis\_ptv\_1p)\] = 0

### Allele based fisher test

\# synonymous variants:

merged\_data\_syns\_fisher <- merged\_data\_syns\_1p %>%

  rowwise() %>%

  mutate(

    fisher\_p\_value = fisher.test(matrix(c(

      (het\_control + hom\_alt\_control \* 2),

      (hom\_ref\_control \* 2 + het\_control),

      (het\_test + hom\_alt\_test \* 2),

      (hom\_ref\_test \* 2 + het\_test)

    ), nrow = 2, byrow = TRUE))$p.value,

    odds\_ratio = fisher.test(matrix(c(

      (het\_control + hom\_alt\_control \* 2),

      (hom\_ref\_control \* 2 + het\_control),

      (het\_test + hom\_alt\_test \* 2),

      (hom\_ref\_test \* 2 + het\_test)

    ), nrow = 2, byrow = TRUE))$estimate)

\# missense+ptv variants:

merged\_data\_mis\_ptv\_fisher <- merged\_data\_mis\_ptv\_1p %>%

  rowwise() %>%

  mutate(

    fisher\_p\_value = fisher.test(matrix(c(

      (het\_control + hom\_alt\_control \* 2),

      (hom\_ref\_control \* 2 + het\_control),

      (het\_test + hom\_alt\_test \* 2),

      (hom\_ref\_test \* 2 + het\_test)

    ), nrow = 2, byrow = TRUE))$p.value,

    odds\_ratio = fisher.test(matrix(c(

      (het\_control + hom\_alt\_control \* 2),

      (hom\_ref\_control \* 2 + het\_control),

      (het\_test + hom\_alt\_test \* 2),

      (hom\_ref\_test \* 2 + het\_test)

    ), nrow = 2, byrow = TRUE))$estimate

  )

### Calculate p value distributions for visualization

\# synonymous variants:

merged\_data\_syns\_fisher$fisher\_p\_value\[merged\_data\_syns\_fisher$fisher\_p\_value == 0\] <- .Machine$double.xmin

merged\_data\_syns\_fisher$fisher\_p\_value <- as.numeric(merged\_data\_syns\_fisher$fisher\_p\_value)

merged\_data\_syns\_fisher$neg\_log\_p <- -log10(merged\_data\_syns\_fisher$fisher\_p\_value)

merged\_data\_syns\_fisher$theoretical <- -log10((rank(merged\_data\_syns\_fisher$fisher\_p\_value) - 0.5) / length(merged\_data\_syns\_fisher$fisher\_p\_value))

\# missense+ptv variants:

merged\_data\_mis\_ptv\_fisher$fisher\_p\_value\[merged\_data\_mis\_ptv\_fisher$fisher\_p\_value == 0\] <- .Machine$double.xmin

merged\_data\_mis\_ptv\_fisher$fisher\_p\_value <- as.numeric(merged\_data\_mis\_ptv\_fisher$fisher\_p\_value)

merged\_data\_mis\_ptv\_fisher$neg\_log\_p <- -log10(merged\_data\_mis\_ptv\_fisher$fisher\_p\_value)

merged\_data\_mis\_ptv\_fisher$theoretical <- -log10((rank(merged\_data\_mis\_ptv\_fisher$fisher\_p\_value) - 0.5) / length(merged\_data\_mis\_ptv\_fisher$fisher\_p\_value))

### Save the result table

\# synonymous variants:

write.table(merged\_data\_syns\_fisher, "for\_qpplot\_carriers\_syns\_1p\_new\_p\_value.tsv", sep="\\t")

\# missense+ptv variants:

write.table(merged\_data\_mis\_ptv\_fisher, "for\_qpplot\_carriers\_mis\_ptv\_1p\_new\_p\_value.tsv", sep="\\t")

### Genotype based Fisher test

\# synonymous variants:

fisher\_p\_values\_syns <- numeric(nrow(merged\_data\_syns\_1p))

odds\_ratios\_syns <- numeric(nrow(merged\_data\_syns\_1p))

total\_control <- 496

total\_gbr <- 91

for (i in 1:nrow(merged\_data\_syns\_1p)) {

  het\_control <- merged\_data\_syns\_1p$het\_control\[i\]

  het\_test <- merged\_data\_syns\_1p$het\_test\[i\]

  remaining\_control <- total\_control - het\_control

  remaining\_gbr <- total\_gbr - het\_test

    if (any(c(het\_control, het\_test, remaining\_control, remaining\_gbr) < 0)) {

    next

  }

    matrix\_data <- matrix(c(het\_control, remaining\_control, het\_test, remaining\_gbr), nrow = 2)

    if (any(rowSums(matrix\_data) == 0) || any(colSums(matrix\_data) == 0)) {

    next

  }

  fisher\_test <- fisher.test(matrix\_data)

  fisher\_p\_values\_syns\[i\] <- fisher\_test$p.value

  odds\_ratios\_syns\[i\] <- ifelse(is.finite(fisher\_test$estimate), fisher\_test$estimate, NA)

}

results\_syns\_genes\_fisher <- data.frame(

  gene = merged\_data\_syns\_1p$gene,  fisher\_p\_value = fisher\_p\_values\_syns,

  odds\_ratio = odds\_ratios\_syns)

\# missense+ptv variants:

fisher\_p\_values\_mis\_ptv <- numeric(nrow(merged\_data\_mis\_ptv\_1p))

odds\_ratios\_mis\_ptv <- numeric(nrow(merged\_data\_mis\_ptv\_1p))

total\_control <- 496

total\_gbr <- 91

for (i in 1:nrow(merged\_data\_syns\_1p)) {

  het\_control <- merged\_data\_syns\_1p$het\_control\[i\]

  het\_test <- merged\_data\_syns\_1p$het\_test\[i\]

  remaining\_control <- total\_control - het\_control

  remaining\_gbr <- total\_gbr - het\_test

    if (any(c(het\_control, het\_test, remaining\_control, remaining\_gbr) < 0)) {

    next

  }

    matrix\_data <- matrix(c(het\_control, remaining\_control, het\_test, remaining\_gbr), nrow = 2)

    if (any(rowSums(matrix\_data) == 0) || any(colSums(matrix\_data) == 0)) {

    next

  }

  fisher\_test <- fisher.test(matrix\_data)

  fisher\_p\_values\_syns\[i\] <- fisher\_test$p.value

  odds\_ratios\_syns\[i\] <- ifelse(is.finite(fisher\_test$estimate), fisher\_test$estimate, NA)

}

results\_mis\_ptv\_genes\_fisher <- data.frame(

  gene = merged\_data\_mis\_ptv\_1p$gene,  fisher\_p\_value = fisher\_p\_values\_mis\_ptv,

  odds\_ratio = odds\_ratios\_mis\_ptv)

### Calculate p value distributions for visualization

\# synonymous variants:

results\_syns\_genes\_fisher$fisher\_p\_value\[results\_syns\_genes\_fisher$fisher\_p\_value == 0\] <- .Machine$double.xmin

results\_syns\_genes\_fisher$fisher\_p\_value <- as.numeric(results\_syns\_genes\_fisher$fisher\_p\_value)

results\_syns\_genes\_fisher$neg\_log\_p <- -log10(results\_syns\_genes\_fisher$fisher\_p\_value)

results\_syns\_genes\_fisher$theoretical <- -log10((rank(results\_syns\_genes\_fisher$fisher\_p\_value) - 0.5) / length(results\_syns\_genes\_fisher$fisher\_p\_value))

\# missense+ptv variants:

results\_mis\_ptv\_genes\_fisher$fisher\_p\_value\[results\_mis\_ptv\_genes\_fisher$fisher\_p\_value == 0\] <- .Machine$double.xmin

results\_mis\_ptv\_genes\_fisher$fisher\_p\_value <- as.numeric(results\_mis\_ptv\_genes\_fisher$fisher\_p\_value)

results\_mis\_ptv\_genes\_fisher$neg\_log\_p <- -log10(results\_mis\_ptv\_genes\_fisher$fisher\_p\_value)

results\_mis\_ptv\_genes\_fisher$theoretical <- -log10((rank(results\_mis\_ptv\_genes\_fisher$fisher\_p\_value) - 0.5) / length(results\_mis\_ptv\_genes\_fisher$fisher\_p\_value))

### Save the result table

\# synonymous variants:

write.table(results\_syns\_genes\_fisher, "for\_qpplot\_genes\_syns\_1p\_new\_p\_value.tsv", sep="\\t")

\# missense+ptv variants:

write.table(results\_mis\_ptv\_genes\_fisher, "for\_qpplot\_genes\_mis\_ptv\_1p\_new\_p\_value.tsv", sep="\\t")

### QQPlots

import pandas as pd

import numpy as np

import scipy.stats as stats

import matplotlib.pyplot as plt

results\_df = pd.read\_csv("input table" ,sep = "\\t")

print(results\_df)

expected = results\_df\['theoretical'\]

observed = results\_df\['neg\_log\_p'\]

plt.figure(figsize=(6, 6))

plt.scatter(expected, observed, edgecolor='black') plt.plot(\[expected.min(), expected.max()\], \[expected.min(), expected.max()\], color='red', linestyle='--')

plt.xlabel("Expected -log10(p-value)")

plt.ylabel("Observed -log10(p-value)")

plt.title("QQ-plot of p-values")

plt.show()
