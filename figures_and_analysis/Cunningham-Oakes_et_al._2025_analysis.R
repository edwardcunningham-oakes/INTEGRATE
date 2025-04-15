#----Generating Figure 1 for Cunningham-Oakes et al. 2025----
#Convert biom to taxonomy for INTEGRATE manuscript

# Ensure BiocManager is installed (used for Bioconductor packages)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Define required CRAN packages
cran_packages <- c(
  "ggplot2", "dplyr", "tidyverse", "tidyr", "data.table", "patchwork",
  "firatheme", "svglite", "corrplot", "ggcorrplot", "viridis",
  "hrbrthemes", "ggpubr"
)

# Define required Bioconductor packages
bioc_packages <- c(
  "phyloseq", "DESeq2", "Maaslin2", "MicrobiotaProcess"
)

# Install missing CRAN packages
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Install missing Bioconductor packages
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

# Load all packages into session
all_packages <- c(cran_packages, bioc_packages)
invisible(lapply(all_packages, library, character.only = TRUE))

# Import biom file and convert using Greengenes taxonomy
INTEGRATE_biom <- "INTEGRATE_taxonomy_table.biom"
INTEGRATE_phyloseq <- import_biom(INTEGRATE_biom, parseFunction=parse_taxonomy_greengenes)

# Filter: retain rows with both Genus and Species annotations
no.na = !is.na(phyloseq::tax_table(INTEGRATE_phyloseq)[,"Genus"]) & !is.na(phyloseq::tax_table(INTEGRATE_phyloseq)[,"Species"])

# Concatenate Genus and Species to create full species names
phyloseq::tax_table(INTEGRATE_phyloseq)[no.na][,"Species"] = paste(
  phyloseq::tax_table(INTEGRATE_phyloseq)[no.na][,"Genus"],
  phyloseq::tax_table(INTEGRATE_phyloseq)[no.na][,"Species"]
)

# Create complete taxa abundance data frame
alltaxatab <- get_alltaxadf(INTEGRATE_phyloseq, method = "count")
INTEGRATE_df = alltaxatab

# Load sample metadata and merge with abundance data
INTEGRATE_meta_df <- read.table(file = "INTEGRATE_sample_metadata.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
INTEGRATE_df <- merge(INTEGRATE_df, INTEGRATE_meta_df, by=0, all=TRUE, row.names = TRUE)
INTEGRATE_df <- INTEGRATE_df %>% column_to_rownames(var = "Row.names")

# Separate into DNA and RNA datasets based on sample identifiers
INTEGRATE_DNA_df <- INTEGRATE_df %>% filter(row.names(INTEGRATE_df) %like% "DNA")
INTEGRATE_RNA_df <- INTEGRATE_df %>% filter(row.names(INTEGRATE_df) %like% "RNA")

#----Generate DNA vs RNA comparison plots (Figure 1 panels)----#

# Each block below processes and plots a specific taxon
# Steps:
# 1. Create data frame of DNA and RNA read counts
# 2. Log-transform values
# 3. Remove non-finite values
# 4. Create scatterplot with identity line

# Adenovirus
adenovirus_frame <- data.frame(
  DNA_reads = INTEGRATE_DNA_df$f__Adenoviridae,
  RNA_reads = INTEGRATE_RNA_df$f__Adenoviridae
)
adenovirus_frame$log_RNA_reads = log(adenovirus_frame$RNA_reads)
adenovirus_frame$log_DNA_reads = log(adenovirus_frame$DNA_reads)
adenovirus_frame <- adenovirus_frame %>% mutate(newCol = 100)
adenovirus_frame[is.na(adenovirus_frame) | adenovirus_frame == "Inf"] = NA
adenovirus_frame[is.na(adenovirus_frame) | adenovirus_frame == "-Inf"] = NA

adenovirus_corr_plot <- ggplot(adenovirus_frame, aes(log_DNA_reads, log_RNA_reads)) +
  geom_point(alpha = 0.5, color = 'black', aes(color = newCol)) +
  geom_abline(intercept = 0, slope = 1, size = 0.5) +
  xlim(0, 20) + ylim(0, 20) + theme_fira() +
  theme(legend.position = "none", plot.title = element_text(size = 14, hjust = 0.5, face = "bold")) +
  ggtitle("Adenovirus") + xlab("log DNA kmers") + ylab("log RNA kmers") +
  scale_fill_fira()

adenovirus_fit <- lm(log_DNA_reads ~ log_RNA_reads, data = adenovirus_frame)

# This structure is repeated for each organism below
# (Only comment shown once per block to avoid redundancy)

campylobacter_frame <- data.frame(DNA_reads=(INTEGRATE_DNA_df$g__Campylobacter),
                   RNA_reads=(INTEGRATE_RNA_df$g__Campylobacter))

campylobacter_frame$log_RNA_reads = log(campylobacter_frame$RNA_reads)
campylobacter_frame$log_DNA_reads = log(campylobacter_frame$DNA_reads)
campylobacter_frame <- campylobacter_frame %>% mutate(newCol = 100)
campylobacter_frame[is.na(campylobacter_frame) | campylobacter_frame=="Inf"] = NA
campylobacter_frame[is.na(campylobacter_frame) | campylobacter_frame=="-Inf"] = NA

campylobacter_corr_plot <- ggplot(campylobacter_frame,aes(log_DNA_reads, log_RNA_reads)) +
 geom_point(alpha=0.5, color='black', aes(color=newCol)) +
 #geom_smooth(method='lm', color='#2e2ea1') + 
 geom_abline(intercept = 0, slope = 1, size = 0.5) +
 xlim(0, 20) + ylim(0, 20) + theme_fira() +
 theme(legend.position="none", plot.title = element_text(size=14, hjust=0.5, face = "bold.italic")) +
 ggtitle("Campylobacter") + xlab("log DNA kmers") + ylab ("log RNA kmers") + 
 scale_fill_fira()

campylobacter_fit <- lm(log_DNA_reads ~ log_RNA_reads, data = campylobacter_frame)

clostridioides_difficile_frame <- data.frame(DNA_reads=(INTEGRATE_DNA_df$`s__Clostridioides difficile`),
                   RNA_reads=(INTEGRATE_RNA_df$`s__Clostridioides difficile`))

clostridioides_difficile_frame$log_RNA_reads = log(clostridioides_difficile_frame$RNA_reads)
clostridioides_difficile_frame$log_DNA_reads = log(clostridioides_difficile_frame$DNA_reads)
clostridioides_difficile_frame <- clostridioides_difficile_frame %>% mutate(newCol = 100)
clostridioides_difficile_frame[is.na(clostridioides_difficile_frame) | clostridioides_difficile_frame=="Inf"] = NA
clostridioides_difficile_frame[is.na(clostridioides_difficile_frame) | clostridioides_difficile_frame=="-Inf"] = NA

clostridioides_difficile_corr_plot <- ggplot(clostridioides_difficile_frame,aes(log_DNA_reads, log_RNA_reads)) +
 geom_point(alpha=0.5, color='black', aes(color=newCol)) +
 #geom_smooth(method='lm', color='#2e2ea1') + 
 geom_abline(intercept = 0, slope = 1, size = 0.5) +
 xlim(0, 20) + ylim(0, 20) + theme_fira() +
 theme(legend.position="none", plot.title = element_text(size=14, hjust=0.5, face = "bold.italic")) +
 ggtitle("Clostridioides difficile") + xlab("log DNA kmers") + ylab ("log RNA kmers") + 
 scale_fill_fira()

clostridioides_difficile_fit <- lm(log_DNA_reads ~ log_RNA_reads, data = clostridioides_difficile_frame)

Yersinia_enterocolitica_frame <- data.frame(DNA_reads=(INTEGRATE_DNA_df$`s__Escherichia coli`),
                   RNA_reads=(INTEGRATE_RNA_df$`s__Escherichia coli`))

Yersinia_enterocolitica_frame$log_RNA_reads = log(Yersinia_enterocolitica_frame$RNA_reads)
Yersinia_enterocolitica_frame$log_DNA_reads = log(Yersinia_enterocolitica_frame$DNA_reads)
Yersinia_enterocolitica_frame <- Yersinia_enterocolitica_frame %>% mutate(newCol = 100)
Yersinia_enterocolitica_frame[is.na(Yersinia_enterocolitica_frame) | Yersinia_enterocolitica_frame=="Inf"] = NA
Yersinia_enterocolitica_frame[is.na(Yersinia_enterocolitica_frame) | Yersinia_enterocolitica_frame=="-Inf"] = NA

Yersinia_enterocolitica_corr_plot <- ggplot(Yersinia_enterocolitica_frame,aes(log_DNA_reads, log_RNA_reads)) +
 geom_point(alpha=0.5, color='black', aes(color=newCol)) +
 #geom_smooth(method='lm', color='#2e2ea1') + 
 geom_abline(intercept = 0, slope = 1, size = 0.5) +
 xlim(0, 20) + ylim(0, 20) + theme_fira() +
 theme(legend.position="none", plot.title = element_text(size=14, hjust=0.5, face = "bold.italic")) +
 ggtitle("Escherichia coli") + xlab("log DNA kmers") + ylab ("log RNA kmers") + 
 scale_fill_fira()

Yersinia_enterocolitica_fit <- lm(log_DNA_reads ~ log_RNA_reads, data = Yersinia_enterocolitica_frame)

salmonella_frame <- data.frame(DNA_reads=(INTEGRATE_DNA_df$g__Salmonella),
                   RNA_reads=(INTEGRATE_RNA_df$g__Salmonella))

salmonella_frame$log_RNA_reads = log(salmonella_frame$RNA_reads)
salmonella_frame$log_DNA_reads = log(salmonella_frame$DNA_reads)
salmonella_frame <- salmonella_frame %>% mutate(newCol = 100)
salmonella_frame[is.na(salmonella_frame) | salmonella_frame=="Inf"] = NA
salmonella_frame[is.na(salmonella_frame) | salmonella_frame=="-Inf"] = NA

salmonella_corr_plot <- ggplot(salmonella_frame,aes(log_DNA_reads, log_RNA_reads)) +
 geom_point(alpha=0.5, color='black', aes(color=newCol)) +
 #geom_smooth(method='lm', color='#2e2ea1') + 
 geom_abline(intercept = 0, slope = 1, size = 0.5) +
 xlim(0, 20) + ylim(0, 20) + theme_fira() +
 theme(legend.position="none", plot.title = element_text(size=14, hjust=0.5, face = "bold.italic")) +
 ggtitle("Salmonella") + xlab("log DNA kmers") + ylab ("log RNA kmers") + 
 scale_fill_fira()

salmonella_fit <- lm(log_DNA_reads ~ log_RNA_reads, data = salmonella_frame)

shigella_frame <- data.frame(DNA_reads=(INTEGRATE_DNA_df$g__Shigella),
                   RNA_reads=(INTEGRATE_RNA_df$g__Shigella))

shigella_frame$log_RNA_reads = log(shigella_frame$RNA_reads)
shigella_frame$log_DNA_reads = log(shigella_frame$DNA_reads)
shigella_frame <- shigella_frame %>% mutate(newCol = 100)
shigella_frame[is.na(shigella_frame) | shigella_frame=="Inf"] = NA
shigella_frame[is.na(shigella_frame) | shigella_frame=="-Inf"] = NA

shigella_corr_plot <- ggplot(shigella_frame,aes(log_DNA_reads, log_RNA_reads)) +
 geom_point(alpha=0.5, color='black', aes(color=newCol)) +
 #geom_smooth(method='lm', color='#2e2ea1') + 
 geom_abline(intercept = 0, slope = 1, size = 0.5) +
 xlim(0, 20) + ylim(0, 20) + theme_fira() +
 theme(legend.position="none", plot.title = element_text(size=14, hjust=0.5, face = "bold.italic")) +
 ggtitle("Shigella") + xlab("log DNA kmers") + ylab ("log RNA kmers") + 
 scale_fill_fira()

shigella_fit <- lm(log_DNA_reads ~ log_RNA_reads, data = shigella_frame)

yersinia_enterocolitica_frame <- data.frame(DNA_reads=(INTEGRATE_DNA_df$`s__Yersinia enterocolitica`),
                   RNA_reads=(INTEGRATE_RNA_df$`s__Yersinia enterocolitica`))

yersinia_enterocolitica_frame$log_RNA_reads = log(yersinia_enterocolitica_frame$RNA_reads)
yersinia_enterocolitica_frame$log_DNA_reads = log(yersinia_enterocolitica_frame$DNA_reads)
yersinia_enterocolitica_frame <- yersinia_enterocolitica_frame %>% mutate(newCol = 100)
yersinia_enterocolitica_frame[is.na(yersinia_enterocolitica_frame) | yersinia_enterocolitica_frame=="Inf"] = NA
yersinia_enterocolitica_frame[is.na(yersinia_enterocolitica_frame) | yersinia_enterocolitica_frame=="-Inf"] = NA

yersinia_enterocolitica_corr_plot <- ggplot(yersinia_enterocolitica_frame,aes(log_DNA_reads, log_RNA_reads)) +
 geom_point(alpha=0.5, color='black', aes(color=newCol)) +
 #geom_smooth(method='lm', color='#2e2ea1') + 
 geom_abline(intercept = 0, slope = 1, size = 0.5) +
 xlim(0, 20) + ylim(0, 20) + theme_fira() +
 theme(legend.position="none", plot.title = element_text(size=14, hjust=0.5, face = "bold.italic")) +
 ggtitle("Yersinia enterocolitica") + xlab("log DNA kmers") + ylab ("log RNA kmers") + 
 scale_fill_fira()

yersinia_enterocolitica_fit <- lm(log_DNA_reads ~ log_RNA_reads, data = yersinia_enterocolitica_frame)

vibrio_cholerae_frame <- data.frame(DNA_reads=(INTEGRATE_DNA_df$`s__Vibrio cholerae`),
                   RNA_reads=(INTEGRATE_RNA_df$`s__Vibrio cholerae`))

vibrio_cholerae_frame$log_RNA_reads = log(vibrio_cholerae_frame$RNA_reads)
vibrio_cholerae_frame$log_DNA_reads = log(vibrio_cholerae_frame$DNA_reads)
vibrio_cholerae_frame <- vibrio_cholerae_frame %>% mutate(newCol = 100)
vibrio_cholerae_frame[is.na(vibrio_cholerae_frame) | vibrio_cholerae_frame=="Inf"] = NA
vibrio_cholerae_frame[is.na(vibrio_cholerae_frame) | vibrio_cholerae_frame=="-Inf"] = NA

vibrio_cholerae_corr_plot <- ggplot(vibrio_cholerae_frame,aes(log_DNA_reads, log_RNA_reads)) +
 geom_point(alpha=0.5, color='black', aes(color=newCol)) +
 #geom_smooth(method='lm', color='#2e2ea1') + 
 geom_abline(intercept = 0, slope = 1, size = 0.5) +
 xlim(0, 20) + ylim(0, 20) + theme_fira() +
 theme(legend.position="none", plot.title = element_text(size=14, hjust=0.5, face = "bold.italic")) +
 ggtitle("Vibrio cholerae") + xlab("log DNA kmers") + ylab ("log RNA kmers") + 
 scale_fill_fira()

vibrio_cholerae_fit <- lm(log_DNA_reads ~ log_RNA_reads, data = vibrio_cholerae_frame)

cryptosporidium_frame <- data.frame(DNA_reads=(INTEGRATE_DNA_df$g__Cryptosporidium),
                   RNA_reads=(INTEGRATE_RNA_df$g__Cryptosporidium))

cryptosporidium_frame$log_RNA_reads = log(cryptosporidium_frame$RNA_reads)
cryptosporidium_frame$log_DNA_reads = log(cryptosporidium_frame$DNA_reads)
cryptosporidium_frame <- cryptosporidium_frame %>% mutate(newCol = 100)
cryptosporidium_frame[is.na(cryptosporidium_frame) | cryptosporidium_frame=="Inf"] = NA
cryptosporidium_frame[is.na(cryptosporidium_frame) | cryptosporidium_frame=="-Inf"] = NA

cryptosporidium_corr_plot <- ggplot(cryptosporidium_frame,aes(log_DNA_reads, log_RNA_reads)) +
 geom_point(alpha=0.5, color='black', aes(color=newCol)) +
 #geom_smooth(method='lm', color='#2e2ea1') + 
 geom_abline(intercept = 0, slope = 1, size = 0.5) +
 xlim(0, 20) + ylim(0, 20) + theme_fira() +
 theme(legend.position="none", plot.title = element_text(size=14, hjust=0.5, face = "bold.italic")) +
 ggtitle("Cryptosporidium") + xlab("log DNA kmers") + ylab ("log RNA kmers") + 
 scale_fill_fira()

cryptosporidium_fit <- lm(log_DNA_reads ~ log_RNA_reads, data = cryptosporidium_frame)

giardia_frame <- data.frame(DNA_reads=(INTEGRATE_DNA_df$g__Giardia),
                   RNA_reads=(INTEGRATE_RNA_df$g__Giardia))

giardia_frame$log_RNA_reads = log(giardia_frame$RNA_reads)
giardia_frame$log_DNA_reads = log(giardia_frame$DNA_reads)
giardia_frame <- giardia_frame %>% mutate(newCol = 100)
giardia_frame[is.na(giardia_frame) | giardia_frame=="Inf"] = NA
giardia_frame[is.na(giardia_frame) | giardia_frame=="-Inf"] = NA

giardia_corr_plot <- ggplot(giardia_frame,aes(log_DNA_reads, log_RNA_reads)) +
 geom_point(alpha=0.5, color='black', aes(color=newCol)) +
 #geom_smooth(method='lm', color='#2e2ea1') + 
 geom_abline(intercept = 0, slope = 1, size = 0.5) +
 xlim(0, 20) + ylim(0, 20) + theme_fira() +
 theme(legend.position="none", plot.title = element_text(size=14, hjust=0.5, face = "bold.italic")) +
 ggtitle("Giardia") + xlab("log DNA kmers") + ylab ("log RNA kmers") + 
 scale_fill_fira()

giardia_fit <- lm(log_DNA_reads ~ log_RNA_reads, data = giardia_frame)

Figure_1 <- adenovirus_corr_plot + campylobacter_corr_plot + clostridioides_difficile_corr_plot +
  Yersinia_enterocolitica_corr_plot + salmonella_corr_plot + shigella_corr_plot +
  vibrio_cholerae_corr_plot + yersinia_enterocolitica_corr_plot + cryptosporidium_corr_plot +
  giardia_corr_plot + plot_layout(ncol = 2, widths = 2, heights = 2)

svglite("Figure_1_Cunningham-Oakes_et_al.svg", width = 10, height = 10)
plot(Figure_1)
dev.off()

#----Generate Figure 4 for Cunningham-Oakes et al. 2025----
ratios <- read.csv("INTEGRATE_DNA_RNA_ratios.tab",sep = "\t", header = T, row.names = 1)
ratios <- log(ratios)
metadata <- read.csv("INTEGRATE_ratio_metadata.tab",sep = "\t", header = T, row.names = 1)

merged_table <- merge(ratios,metadata, by = 0, all = TRUE)
merged_table2 <- merged_table[,-1]
rownames(merged_table2) <- merged_table[,1]

merged_table2$campylobacter_traditional_result <- as.factor(merged_table2$campylobacter_traditional_result)
merged_table2$campylobacter_luminex_result <- as.factor(merged_table2$campylobacter_luminex_result)
merged_table2$clostridioides_difficile_traditional_result <- as.factor(merged_table2$clostridioides_difficile_traditional_result)
merged_table2$clostridioides_difficile_luminex_result <- as.factor(merged_table2$clostridioides_difficile_luminex_result)
merged_table2$vibrio_cholerae_traditional_result <- as.factor(merged_table2$vibrio_cholerae_traditional_result)
merged_table2$vibrio_cholerae_luminex_result <- as.factor(merged_table2$vibrio_cholerae_luminex_result)

Campylobacter_traditional_ratios <- merged_table2 %>% drop_na(Campylobacter,campylobacter_traditional_result) %>% 
  ggplot( aes(x=campylobacter_traditional_result, y=Campylobacter, fill=campylobacter_traditional_result)) +
    geom_violin() +
    scale_fill_manual(values = c('#bcbddc', '#fdae6b')) +
    stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               colour = "white") +
    theme_ipsum() +
    theme(
      legend.position="none",
      plot.title = element_text(size=11), axis.title.x = element_text(angle=0, hjust = 0.5), axis.title.y = element_text(angle=0, vjust = 0.5)) +
    ggtitle("Campylobacter Traditional") + geom_jitter(alpha = 0.5) +
    xlab("Diagnostic result") + ylab("log DNA/RNA ratio") + stat_compare_means(method="wilcox.test", size = 3.5, position = position_nudge(x= 0.25)) + ylim(-5.0,7.0)

Campylobacter_luminex_ratios <- merged_table2 %>% drop_na(Campylobacter,campylobacter_luminex_result) %>% 
  ggplot( aes(x=campylobacter_luminex_result, y=Campylobacter, fill=campylobacter_luminex_result)) +
    geom_violin() +
    scale_fill_manual(values = c('#bcbddc', '#fdae6b')) +
    stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               colour = "white") +
    theme_ipsum() +
    theme(
      legend.position="none",
      plot.title = element_text(size=11), axis.title.x = element_text(angle=0, hjust = 0.5), axis.title.y = element_text(angle=0, vjust = 0.5)) +
    ggtitle("Campylobacter Luminex") + geom_jitter(alpha = 0.5) +
    xlab("Diagnostic result") + ylab("log DNA/RNA ratio") + stat_compare_means(method="wilcox.test", size = 3.5, position = position_nudge(x= 0.25)) + ylim(-5.0,7.0)

clostridioides_difficile_traditional_ratios <- merged_table2 %>% drop_na(Clostridioides_difficile,clostridioides_difficile_traditional_result) %>% 
  ggplot( aes(x=clostridioides_difficile_traditional_result, y=Clostridioides_difficile, fill=clostridioides_difficile_traditional_result)) +
    geom_violin() +
    scale_fill_manual(values = c('#bcbddc', '#fdae6b')) +
    stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               colour = "white") +
    theme_ipsum() +
    theme(
      legend.position="none",
      plot.title = element_text(size=11), axis.title.x = element_text(angle=0, hjust = 0.5), axis.title.y = element_text(angle=0, vjust = 0.5)) +
    ggtitle("Clostridioides difficile Traditional") + geom_jitter(alpha = 0.5) +
    xlab("Diagnostic result") + ylab("log DNA/RNA ratio") + stat_compare_means(method="wilcox.test", size = 3.5, position = position_nudge(x= 0.25)) + ylim(-5.0,7.0)

clostridioides_difficile_luminex_ratios <- merged_table2 %>% drop_na(Clostridioides_difficile,clostridioides_difficile_luminex_result) %>% 
  ggplot( aes(x=clostridioides_difficile_luminex_result, y=Clostridioides_difficile, fill=clostridioides_difficile_luminex_result)) +
    geom_violin() +
    scale_fill_manual(values = c('#bcbddc', '#fdae6b')) +
    stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               colour = "white") +
    theme_ipsum() +
   theme(
      legend.position="none",
      plot.title = element_text(size=11), axis.title.x = element_text(angle=0, hjust = 0.5), axis.title.y = element_text(angle=0, vjust = 0.5)) +
    ggtitle("Clostridioides difficile Luminex") + geom_jitter(alpha = 0.5) +
   xlab("Diagnostic result") + ylab("log DNA/RNA ratio") + stat_compare_means(method="wilcox.test", size = 3.5, position = position_nudge(x= 0.25)) + ylim(-5.0,7.0)


Figure_4 <- Campylobacter_traditional_ratios + Campylobacter_luminex_ratios +
clostridioides_difficile_traditional_ratios + clostridioides_difficile_luminex_ratios +
plot_layout(ncol = 2, nrow = 2, widths = 12, heights = 12)

Figure_4 <- Figure_4 + 
  plot_annotation(tag_levels = 'A')

ggsave(file="Figure_4_Cunningham-Oakes_et_al.svg", plot=Figure_4, width=12, height=12)

#----Generate Supplementary Figure 1 (correlation between pathogen kmers and all diagnostic test results)----

plot_size <- as.matrix(read.table("INTEGRATE_kmer_diagnostics_coefficients.txt", sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE))

plot_significance <- as.matrix(read.table("INTEGRATE_kmer_diagnostics_significance.txt", sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE))

#Define colours
col <- colorRampPalette(c("#67001f","#b2182b","#d6604d","#f4a582","#92c5de","#4393c3","#2166ac","#053061"))

#Correlation plot of significance
svglite("Supplementary_Figure_1_Cunningham-Oakes_et_al.svg", width = 10, height = 10)
corrplot(plot_size, method="color", is.corr = FALSE, p.mat = plot_significance, sig.level = c(0.001, 0.01, 0.05, 0.25), insig = "label_sig", pch.cex = 0.8, pch.col = "white",
        na.label = "square", na.label.col = "black", col = col(200), addgrid.col = "#636363", tl.cex = .8, tl.srt = 30,
        tl.col = ("black"))
dev.off()
