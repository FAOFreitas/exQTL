library(dplyr)
library(ggplot2)

# Helper function to calculate start and stop positions
calculate_positions <- function(gene, range, Genes.df, gbuild) {
  rangebp <- range * 1000
  gene_info <- Genes.df %>% filter(Gene == gene & Build == gbuild)
  startpos <- min(gene_info$Start) - rangebp
  stoppos <- max(gene_info$Stop) + rangebp
  chromosome <- gene_info$CHR[1]
  list(startpos = startpos, stoppos = stoppos, chromosome = chromosome)
}

# Helper function to subset data based on gene, trait, and range
subset_gwas_data <- function(GWAS.df, chromosome, startpos, stoppos, trait, sigpvalue_GWAS) {
  gwas_data <- GWAS.df %>%
    filter(CHR == chromosome & PHE == trait & BP >= startpos & BP <= stoppos & !is.na(P) & !is.na(BETA))
  if (nrow(gwas_data) == 0) stop('No SNPs found in the specified range for the gene and trait.')
  if (nrow(gwas_data %>% filter(P <= sigpvalue_GWAS)) == 0) warning('No significant SNPs found within the range.')
  gwas_data
}

# Helper function to filter eQTL data
filter_eqtl_data <- function(eQTL.df, gene, sigpvalue_eQTL, tissue) {
  eqtl_data <- eQTL.df %>%
    filter(Gene.Symbol == gene & P.Value <= sigpvalue_eQTL & !is.na(NES) & !is.na(P.Value))
  if (any(tissue == "all") || length(tissue) >= 2) eqtl_data else eqtl_data %>% filter(Tissue == tissue)
}

# Function to compute co-localization
compute_colocalization <- function(GWAS.df, eQTL.df, Genes.df, gene, trait, sigpvalue_GWAS, sigpvalue_eQTL, tissue, range, gbuild) {
  pos <- calculate_positions(gene, range, Genes.df, gbuild)
  gwas_data <- subset_gwas_data(GWAS.df, pos$chromosome, pos$startpos, pos$stoppos, trait, sigpvalue_GWAS)
  eqtl_data <- filter_eqtl_data(eQTL.df, gene, sigpvalue_eQTL, tissue)
  
  combinedSNPS <- union(gwas_data$SNP, eqtl_data$SNP.Id)
  Combined.eQTL.GWAS.Data <- left_join(
    gwas_data %>% mutate(SNP = factor(SNP, levels = combinedSNPS)),
    eqtl_data %>% mutate(SNP.Id = factor(SNP.Id, levels = combinedSNPS)) %>% rename(SNP = SNP.Id),
    by = "SNP"
  )
  
  if (nrow(Combined.eQTL.GWAS.Data) == 0) stop('No overlapping SNPs between GWAS and eQTL data.')
  
  Combined.eQTL.GWAS.Data <- Combined.eQTL.GWAS.Data %>%
    mutate(
      Congruence = ifelse(BETA * NES < 0, "Incongruent", "Congruent"),
      NeglogeQTLpValue = -log10(P.Value),
      Neglog10pvalue_GWAS = -log10(P),
      Congruence = ifelse(is.na(Congruence), "Non-eQTL", Congruence)
    )
  
  Combined.eQTL.GWAS.Data
}

# Example usage
GWAS.df.example <- data.frame(CHR = c(1, 1), BP = c(100000, 200000), SNP = c("rs1", "rs2"), P = c(1e-7, 1e-6), BETA = c(0.5, -0.4), PHE = "LDL")
eQTL.df.example <- data.frame(SNP.Id = c("rs1", "rs2"), Gene.Symbol = c("ACTN3", "ACTN3"), P.Value = c(1e-5, 1e-4), NES = c(0.2, -0.3), Tissue = c("Liver", "Liver"))
Genes.df.example <- data.frame(Gene = "ACTN3", CHR = 1, Start = 90000, Stop = 210000, Build = "hg19")

results <- compute_colocalization(
  GWAS.df = GWAS.df.example, 
  eQTL.df = eQTL.df.example, 
  Genes.df = Genes.df.example, 
  gene = "ACTN3", 
  trait = "LDL", 
  sigpvalue_GWAS = 5e-8, 
  sigpvalue_eQTL = 0.05, 
  tissue = "all", 
  range = 200, 
  gbuild = "hg19"
)

print(results)
