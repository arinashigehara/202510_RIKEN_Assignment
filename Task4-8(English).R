##202510 RIKEN

##Task4
# 1. Square root of 10
sqrt(10)

# 2. Logarithm base 2 of 32
log2(32)

# 3. Sum of numbers from 1 to 1000
sum(1:1000)

# 4. Sum of even numbers from 2 to 1000
sum(seq(2, 1000, by = 2))

# 5. Number of pairwise comparisons among 100 genes (combinations)
choose(100, 2)

# 6. Number of permutations choosing 3 genes from 100
factorial(100) / factorial(97)

##Task5
# 1. Load built-in CO2 dataset
data(CO2)
head(CO2)

# 2. Check documentation of CO2 dataset
?CO2
# → Dataset of CO2 uptake in plants. Includes Plant, Type, Treatment, conc, uptake, etc.

# 3. Compare average and median CO2 uptake between Quebec and Mississippi
library(dplyr)
CO2 %>%
  group_by(Type) %>%
  summarise(
    mean_uptake = mean(uptake),
    median_uptake = median(uptake)
  )


##Task6
# 1. Function to compute ratio of mean to median
mean_median_ratio <- function(x) {
  mean(x) / median(x)
}

# 2. Function to compute mean after removing min and max
mean_trimmed <- function(x) {
  x_sorted <- sort(x)
  mean(x_sorted[-c(1, length(x_sorted))])
}

# 3. Pipe usage explanation (within 300 characters, no spaces)

# 4. Usefulness of apply-family functions (within 300 characters, no spaces)

##Task7
# Load magic_guys.csv and inspect
library(ggplot2)
magic <- read.csv("/Users/arinashigehara/Documents/202510RIKEN/magic_guys.csv")
glimpse(magic)

# Base R histogram for each species
hist(
  magic$length[magic$species == "jedi"],
  main = "Height Distribution - Jedi",
  xlab = "Height (cm)", col = "skyblue", breaks = 20
)

hist(
  magic$length[magic$species == "sith"],
  main = "Height Distribution - Sith",
  xlab = "Height (cm)", col = "lightpink", breaks = 20
)

# Boxplot of height by species using ggplot2
ggplot(magic, aes(x = species, y = length, fill = species)) +
  geom_boxplot() +
  labs(title = "Height Boxplot by Species", x = "Species", y = "Height (cm)") +
  theme_minimal()

# Save histogram plot in PNG/PDF/SVG formats
p <- ggplot(magic, aes(x = length, fill = species)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
  labs(title = "Histogram of Height", x = "Height (cm)", y = "Count") +
  theme_minimal()

ggsave("hist_plot.png", plot = p, width = 6, height = 4, dpi = 300)
ggsave("hist_plot.pdf", plot = p, width = 6, height = 4)
ggsave("hist_plot.svg", plot = p, width = 6, height = 4)

# Usage of each format:
# - png: suitable for slides or web
# - pdf: high-resolution printing, manuscripts
# - svg: editable vector format (e.g., Illustrator)

# Load and process microarray data
expr <- read.table(
  "/Users/arinashigehara/Downloads/microarray_data.tab",
  header = TRUE,
  sep = "\t",
  row.names = NULL,
  check.names = FALSE
)

# Count and plot missing values per gene
na_per_gene <- rowSums(is.na(expr))
hist(na_per_gene, main = "Missing Values per Gene", xlab = "Number of NA", col = "tomato")

# Identify genes with >10% missing values and impute with row mean
threshold <- 0.1
high_na <- expr[rowSums(is.na(expr)) > ncol(expr) * threshold, ]

expr_imputed <- expr
for (i in 1:nrow(expr_imputed)) {
  row_mean <- mean(expr_imputed[i, ], na.rm = TRUE)
  expr_imputed[i, is.na(expr_imputed[i, ])] <- row_mean
}

# Visualize CO2 uptake by concentration and treatment
ggplot(CO2, aes(x = conc, y = uptake, color = Type)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~Treatment) +
  labs(title = "CO2 Uptake by Concentration", x = "CO2 Concentration", y = "Uptake") +
  theme_minimal()


##Task8
# Load tidybiology and datasets
library(tidybiology)
library(ggplot2)
library(dplyr)

data("chromosome")
data("proteins")

# a. Summary stats for variations, coding genes, and miRNAs
chromosome %>%
  summarise(
    mean_var   = mean(variations),
    median_var = median(variations),
    max_var    = max(variations),
    
    mean_gene   = mean(protein_codinggenes),
    median_gene = median(protein_codinggenes),
    max_gene    = max(protein_codinggenes),
    
    mean_mirna   = mean(mi_rna),
    median_mirna = median(mi_rna),
    max_mirna    = max(mi_rna)
  )

# b. Distribution of chromosome sizes (mm and base pairs)
ggplot(chromosome, aes(x = length_mm)) +
  geom_histogram(bins = 15, fill = "steelblue", color = "black") +
  labs(title = "Distribution of Chromosome Lengths (mm)", x = "Length (mm)", y = "Count") +
  theme_minimal()

ggplot(chromosome, aes(x = basepairs)) +
  geom_histogram(bins = 15, fill = "darkorange", color = "black") +
  scale_x_log10() +
  labs(title = "Distribution of Chromosome Size (Base Pairs)", x = "Base pairs (log10)", y = "Count") +
  theme_minimal()

# c. Correlation between chromosome size and gene/miRNA count
ggplot(chromosome, aes(x = basepairs, y = protein_codinggenes)) +
  geom_point(color = "forestgreen", size = 3) +
  geom_smooth(method = "lm", color = "black") +
  labs(title = "Protein Coding Genes vs Chromosome Size", x = "Base pairs", y = "Protein coding genes") +
  theme_light()

ggplot(chromosome, aes(x = basepairs, y = mi_rna)) +
  geom_point(color = "purple", size = 3) +
  geom_smooth(method = "lm", color = "black") +
  labs(title = "miRNAs vs Chromosome Size", x = "Base pairs", y = "miRNAs") +
  theme_light()

# d. Summary and visualization of protein length and mass
proteins %>%
  summarise(
    mean_length  = mean(length),
    median_length = median(length),
    max_length    = max(length),
    
    mean_mass  = mean(mass),
    median_mass = median(mass),
    max_mass    = max(mass)
  )

ggplot(proteins, aes(x = length, y = mass)) +
  geom_point(alpha = 0.6, color = "dodgerblue", size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(title = "Protein Length vs Mass", x = "Protein Length", y = "Protein Mass") +
  theme_bw(base_size = 13)
