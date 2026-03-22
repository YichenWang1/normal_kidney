library(ggplot2)
library(nlme)
library(dplyr)
library(tidyr)

# -----------Normal SBS----------------
metadata_path <- "../data/bulk_normal_kidney_metadata.csv"
weights_path  <- "../data/SBS96_kidney_signature_weights_table.csv"

metadata <- read.csv(metadata_path, header = TRUE, stringsAsFactors = FALSE)
metadata <- metadata %>%
  mutate(
    Country = recode(Country, "United Kingdom" = "UK", "Czech Republic" = "Czechia"),
    Country = factor(Country)
  )
df_normal <- metadata %>%
  transmute(
    normal = Normal_Kidney,
    cancer = Tumor,
    type = "Normal",
    Country = Country,
    Age = Age,
    Sex = Sex,
    Alcohol = Alcohol,
    Tobacco = Tobacco,
    Burden_SBS = Normal_SBS_burden,
    Burden_ID = Normal_ID_burden
  )
df_normal$Country <- factor(df_normal$Country)

weights_raw <- read.csv(path.expand(weights_path), check.names = FALSE, stringsAsFactors = FALSE)
# keep only samples present in metadata
weights_raw <- weights_raw[weights_raw$Sample %in% metadata$Normal_Kidney, , drop = FALSE]
rownames(weights_raw) <- weights_raw$Sample

# align and subset weights to the same order as df_normal
weights <- weights_raw[ df_normal$normal, , drop = FALSE ]

# ensure numeric signature columns (exclude "Sample" if present)
sig_cols <- setdiff(colnames(weights), "Sample")
weights[sig_cols] <- lapply(weights[sig_cols], function(x) as.numeric(as.character(x)))

# zero small contributions, but preserve SBS2/SBS13 when their sum > 0.05 (based on original raw values)
weights[sig_cols] <- lapply(weights[sig_cols], function(x) { x[x < 0.05] <- 0; x })
rows_restore <- rownames(weights)[ (weights_raw[rownames(weights), "SBS2"] + weights_raw[rownames(weights), "SBS13"]) > 0.05 ]
if (length(rows_restore) > 0) {
  weights[rows_restore, "SBS2"]  <- weights_raw[rows_restore, "SBS2"]
  weights[rows_restore, "SBS13"] <- weights_raw[rows_restore, "SBS13"]
}

# compute Unassigned and clamp negatives to zero
weights$Unassigned <- 1 - rowSums(weights[, sig_cols, drop = FALSE], na.rm = TRUE)
weights[sig_cols] <- lapply(weights[sig_cols], function(x) pmax(x, 0))
weights$Unassigned <- pmax(weights$Unassigned, 0)

# add signature burdens to df_normal
df_normal <- df_normal %>%
  mutate(
    SBS1 = weights$SBS1 * Burden_SBS,
    SBS5 = weights$SBS5 * Burden_SBS,
    SBS12 = weights$SBS12 * Burden_SBS,
    SBS2 = weights$SBS2 * Burden_SBS,
    SBS13 = weights$SBS13 * Burden_SBS,
    APOBEC = (weights$SBS2+weights$SBS13) * Burden_SBS,
    SBS18 = weights$SBS18 * Burden_SBS,
    SBS22a = weights$SBS22a * Burden_SBS,
    SBS22b = weights$SBS22b * Burden_SBS,
    SBS22c = weights$SBS22c * Burden_SBS,
    SBS40a = weights$SBS40a * Burden_SBS,
    SBS40b = weights$SBS40b * Burden_SBS,
    SBS40c = weights$SBS40c * Burden_SBS,
    SBSB = weights$SBSB * Burden_SBS,
    SBSC = weights$SBSC * Burden_SBS,
    SBSD = weights$SBSD * Burden_SBS,
    MSI = (weights$SBS21+weights$SBS44) * Burden_SBS,
    SBS21 = weights$SBS21 * Burden_SBS,
    SBS44 = weights$SBS44 * Burden_SBS,
    # Unassigned = weights$Unassigned * Burden_SBS
  )
df_normal <- df_normal[df_normal$normal!='PD47592c_ds0003',]

signature_names <- c(
  "SBS40a","SBS40b","SBS40c","SBS5","SBS22a",
  "SBS22b","SBS22c","SBS12","SBS18","SBS1",
  # "APOBEC","SBSD","MSI",
  "SBSB","SBSC","Others", "Unassigned"
)

cols <- c(
  "#8DD3C7", "#CFECBB", "#F4F3B9", "#BD98A2",  "#1f78b4","#8AB1C9",
   "#759696", "#F5847A", "#D3B387", "#17BEBB", 
   # "#BFD767", "#CECBD0", "purple",
  "#cab2d6","#FCCDE5", "#E5E4E6","#8A8A8A"
)
names(cols) <- signature_names

sig_levels <- signature_names

# Load Indel data
weights_path  <- "../data/ID83_kidney_signature_ID_weights_table.csv"

weights_raw <- read.csv(path.expand(weights_path), check.names = FALSE, stringsAsFactors = FALSE)
# keep only samples present in metadata
weights_raw <- weights_raw[weights_raw$Sample %in% metadata$Normal_Kidney, , drop = FALSE]
rownames(weights_raw) <- weights_raw$Sample
weights <- weights_raw[df_normal$normal, ,drop = FALSE ]

# ensure numeric signature columns (exclude "Sample" if present)
sig_cols <- setdiff(colnames(weights), "Sample")
weights[sig_cols] <- lapply(weights[sig_cols], function(x) as.numeric(as.character(x)))

# zero small contributions, but preserve SBS2/SBS13 when their sum > 0.05 (based on original raw values)
weights[sig_cols] <- lapply(weights[sig_cols], function(x) { x[x < 0.05] <- 0; x })

# compute Unassigned and clamp negatives to zero
weights$Unassigned <- 1 - rowSums(weights[, sig_cols, drop = FALSE], na.rm = TRUE)
weights[sig_cols] <- lapply(weights[sig_cols], function(x) pmax(x, 0))
weights$Unassigned <- pmax(weights$Unassigned, 0)

df_normal<- df_normal %>%
  mutate(
    ID5 = weights$ID5 * Burden_ID,
    ID8 = weights$ID8 * Burden_ID,
    ID9 = weights$ID9 * Burden_ID,
    ID11 = weights$ID11 * Burden_ID,
    ID21 = weights$ID21 * Burden_ID,
    ID23 = weights$ID23 * Burden_ID,
    ID1 = weights$ID1 * Burden_ID,
    ID2 = weights$ID2 * Burden_ID,
    ID3 = weights$ID3 * Burden_ID,
    ID12 = weights$ID12 * Burden_ID,
    # Unassigned = weights$Unassigned * Burden_ID
  )


# Tobacco associated signatures
pdf("../../figure/mutation rate/country/bulk/SBSB_bulk_rates.pdf",width = 3, height = 3)
lmm <- lme(
  SBSB ~ 0+Age:Country,
  random = ~ 0 + Age | Tobacco/Patient,
  data   = df_normal,
  method = "ML"
)


ci <- intervals(lmm, which = "fixed")$fixed
rates <- data.frame(
  Country = gsub("Age:Country", "", names(fixef(lmm))),
  Rate = as.numeric(fixef(lmm))
)
rates$lower <- ci[, "lower"]
rates$upper <- ci[, "upper"]
rates <- rates[rates$Country!="(Intercept)",]
rates$Std.Error <- summary(lmm)$tTable[, "Std.Error"]
rates$p_value   <- summary(lmm)$tTable[, "p-value"]
# write.csv(rates, "../../data_S1/S7_SBSB_bulk_rates_by_country.csv", row.names = FALSE)
country_order <- c("Japan", "Thailand", "Serbia", "Brazil", "Canada",
"Romania", "UK", "Russia", "Lithuania", "Czechia") #ordered by cancer incidence
rates$Country <- factor(rates$Country, levels = country_order)
avg_rate <- mean(rates$Rate, na.rm = TRUE)
ggplot(rates, aes(x = reorder(Country, Rate), y = Rate)) +
  geom_point(shape = 21, size = 3, fill=cols['SBSB'],color = 'black') +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  coord_flip() + theme_bw()+
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(),
        legend.position = "none",
        axis.line = element_line(size = 1, colour = "black"),
         title=element_text(size=6),
              axis.text.y=element_text(size=6,color="black"),
              axis.text.x=element_text(size=6,color="black"),
              legend.text=element_text(size=6))+
  geom_hline(yintercept = avg_rate,
  linetype = "dashed",
  linewidth = 0.5)+
  labs(x = "Country",
    y = "SBSB mutations per year",
    title = "Bulk kidney cortex")
dev.off()

# Sex associated signatures
pdf("../../figure/mutation rate/country/bulk/SBS40b_bulk_rates.pdf",width = 3, height = 3)
lmm <- lme(
  SBS40b ~ 0+Age:Country,
  random = ~ 0 + Age | Sex,
  data   = df_normal,
  method = "ML"
)
summary(lmm)
ci <- intervals(lmm, which = "fixed")$fixed
rates <- data.frame(
  Country = gsub("Age:Country", "", names(fixef(lmm))),
  Rate = as.numeric(fixef(lmm))
)
rates$lower <- ci[, "lower"]
rates$upper <- ci[, "upper"]
rates$Std.Error <- summary(lmm)$tTable[, "Std.Error"]
rates$p_value   <- summary(lmm)$tTable[, "p-value"]
write.csv(rates, "../../data_S1/S7_SBS40b_bulk_rates_by_country.csv", row.names = FALSE)
country_order <- c("Japan", "Thailand", "Serbia", "Brazil", "Canada",
"Romania", "UK", "Russia", "Lithuania", "Czechia") #ordered by cancer incidence
rates$Country <- factor(rates$Country, levels = country_order)
avg_rate <- mean(rates$Rate, na.rm = TRUE)
# rates <- rates[rates$Country!="(Intercept)",]
ggplot(rates, aes(x = reorder(Country, Rate), y = Rate)) +
  geom_point(shape = 21, size = 3, fill=cols['SBS40b'],color = 'black') +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  coord_flip() + theme_bw()+
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(),
        legend.position = "none",
        axis.line = element_line(size = 1, colour = "black"),
         title=element_text(size=6),
              axis.text.y=element_text(size=6,color="black"),
              axis.text.x=element_text(size=6,color="black"),
              legend.text=element_text(size=6))+
  geom_hline(yintercept = avg_rate,
  linetype = "dashed",
  linewidth = 0.5)+
  labs(x = "Country",
    y = "SBS40b mutations per year",
    title = "Bulk kidney cortex")
dev.off()

# Rest of signatures
pdf("../../figure/mutation rate/country/bulk/SBS40c_bulk_rates.pdf",width = 3, height = 3)
lmm <- lm(
  SBS40c ~  0+Age:Country,
  data   = df_normal
)
summary(lmm)
rates <- data.frame(
  Country = gsub("Age:Country", "", names(coef(lmm))),
  Rate = as.numeric(coef(lmm))
)
ci <- confint(lmm, level = 0.95)
rates$lower <- ci[, "2.5 %"]
rates$upper <- ci[, "97.5 %"]
rates$Std.Error <- summary(lmm)$coefficients[, "Std. Error"]
rates$p_value   <- summary(lmm)$coefficients[, "Pr(>|t|)"]
write.csv(rates, "../../data_S1/S7_SBS40c_bulk_rates_by_country.csv", row.names = FALSE)

country_order <- c("Japan", "Thailand", "Serbia", "Brazil", "Canada",
"Romania", "UK", "Russia", "Lithuania", "Czechia") #ordered by cancer incidence
rates$Country <- factor(rates$Country, levels = country_order)
avg_rate <- mean(rates$Rate, na.rm = TRUE)
# rates <- rates[rates$Country!="(Intercept)",]
ggplot(rates, aes(x = reorder(Country, Rate), y = Rate)) +
  geom_point(shape = 21, size = 3, fill=cols['SBS40c'],color = 'black') +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  coord_flip() + theme_bw()+
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(),
        legend.position = "none",
        axis.line = element_line(size = 1, colour = "black"),
         title=element_text(size=6),
              axis.text.y=element_text(size=6,color="black"),
              axis.text.x=element_text(size=6,color="black"),
              legend.text=element_text(size=6))+
  geom_hline(yintercept = avg_rate,
  linetype = "dashed",
  linewidth = 0.5)+
  labs(x = "Country",
    y = "SBS40c mutations per year",
    title = "Bulk kidney cortex")
dev.off()

# -----------Normal LCM----------------
metadata_path <- "../data/LCM_normal_kidney_metadata.csv"
weights_path  <- "../data/SBS96_kidney_signature_weights_table.csv"
LCM_metadata <- read.csv(metadata_path, header = TRUE, stringsAsFactors = FALSE)
LCM_metadata  <- LCM_metadata  %>%
  mutate(Country = recode(Country, "United Kingdom" = "UK", "Czech Republic" = "Czechia"))

df = data.frame(Patient = LCM_metadata$Patient, Sample = LCM_metadata$Sample, Structure = LCM_metadata$Structure, Country = LCM_metadata$Country, Sex = LCM_metadata$Sex, Age = LCM_metadata$Age, Tobacco = LCM_metadata$Tobacco, Burden_SBS= LCM_metadata$Normal_SBS_burden,Burden_ID= LCM_metadata$Normal_ID_burden)
df$Structure[df$Structure=='Distal_tubules'] ='Distal tubules'
df$Structure[df$Structure=='Proximal_tubules'] ='Proximal tubules'
df$Structure = factor(df$Structure, levels = c('Proximal tubules','Medulla','Glomeruli','Distal tubules'))

weights_raw <- read.csv(path.expand(weights_path), check.names = FALSE, stringsAsFactors = FALSE)
weights_raw$Sample <- sub(" .*", "", weights_raw$Sample)
# keep only samples present in metadata
weights_raw <- weights_raw[weights_raw$Sample %in% df$Sample, , drop = FALSE]
rownames(weights_raw) <- weights_raw$Sample

# align and subset weights to the same order as df_normal
weights <- weights_raw[df$Sample, , drop = FALSE ]

# ensure numeric signature columns (exclude "Sample" if present)
sig_cols <- setdiff(colnames(weights), "Sample")
weights[sig_cols] <- lapply(weights[sig_cols], function(x) as.numeric(as.character(x)))

# zero small contributions, but preserve SBS2/SBS13 when their sum > 0.05 (based on original raw values)
weights[sig_cols] <- lapply(weights[sig_cols], function(x) { x[x < 0.05] <- 0; x })

# compute Unassigned and clamp negatives to zero
weights$Unassigned <- 1 - rowSums(weights[, sig_cols, drop = FALSE], na.rm = TRUE)
weights[sig_cols] <- lapply(weights[sig_cols], function(x) pmax(x, 0))
weights$Unassigned <- pmax(weights$Unassigned, 0)

df<- df %>%
  mutate(
    SBS1 = weights$SBS1 * Burden_SBS,
    SBS5 = weights$SBS5 * Burden_SBS,
    SBS12 = weights$SBS12 * Burden_SBS,
    SBS2 = weights$SBS2 * Burden_SBS,
    SBS13 = weights$SBS13 * Burden_SBS,
    SBS18 = weights$SBS18 * Burden_SBS,
    SBS22a = weights$SBS22a * Burden_SBS,
    SBS22b = weights$SBS22b * Burden_SBS,
    SBS22c = weights$SBS22c * Burden_SBS,
    SBS40a = weights$SBS40a * Burden_SBS,
    SBS40b = weights$SBS40b * Burden_SBS,
    SBS40c = weights$SBS40c * Burden_SBS,
    SBSB = weights$SBSB * Burden_SBS,
    SBSC = weights$SBSC * Burden_SBS,
    SBSD = weights$SBSD * Burden_SBS,
    SBS21 = weights$SBS21 * Burden_SBS,
    SBS44 = weights$SBS44 * Burden_SBS,
    Unassigned = weights$Unassigned * Burden_SBS
  )

head(df)

# Load Indel data
weights_path  <- "../data/ID83_kidney_signature_ID_weights_table.csv"

weights_raw <- read.csv(path.expand(weights_path), check.names = FALSE, stringsAsFactors = FALSE)
weights_raw$Sample <- sub(" .*", "", weights_raw$Sample)
# keep only samples present in metadata
weights_raw <- weights_raw[weights_raw$Sample %in% df$Sample, , drop = FALSE]
rownames(weights_raw) <- weights_raw$Sample

# align and subset weights to the same order as df_normal
weights <- weights_raw[df$Sample, , drop = FALSE ]

# ensure numeric signature columns (exclude "Sample" if present)
sig_cols <- setdiff(colnames(weights), "Sample")
weights[sig_cols] <- lapply(weights[sig_cols], function(x) as.numeric(as.character(x)))

# zero small contributions, but preserve SBS2/SBS13 when their sum > 0.05 (based on original raw values)
weights[sig_cols] <- lapply(weights[sig_cols], function(x) { x[x < 0.05] <- 0; x })

# compute Unassigned and clamp negatives to zero
weights$Unassigned <- 1 - rowSums(weights[, sig_cols, drop = FALSE], na.rm = TRUE)
weights[sig_cols] <- lapply(weights[sig_cols], function(x) pmax(x, 0))
weights$Unassigned <- pmax(weights$Unassigned, 0)

df<- df %>%
  mutate(
    ID5 = weights$ID5 * Burden_ID,
    ID8 = weights$ID8 * Burden_ID,
    ID9 = weights$ID9 * Burden_ID,
    ID11 = weights$ID11 * Burden_ID,
    ID21 = weights$ID21 * Burden_ID,
    ID23 = weights$ID23 * Burden_ID,
    ID1 = weights$ID1 * Burden_ID,
    ID2 = weights$ID2 * Burden_ID,
    ID3 = weights$ID3 * Burden_ID,
    ID12 = weights$ID12 * Burden_ID,
    # Unassigned = weights$Unassigned * Burden_ID
  )



# Sex associated signatures
pdf("../../figure/mutation rate/country/PT/SBS40b_PT_rates.pdf",width = 3, height = 3)
lmm <- lme(
  SBS40b ~ 0+Age:Country,
  random = ~ 0 + Age | Sex/Patient,
  data   = df[df$Structure=='Proximal tubules',],
  method = "ML"
)
summary(lmm)

ci <- intervals(lmm, which = "fixed")$fixed
rates <- data.frame(
  Country = gsub("Age:Country", "", names(fixef(lmm))),
  Rate = as.numeric(fixef(lmm))
)
rates$lower <- ci[, "lower"]
rates$upper <- ci[, "upper"]
rates$Std.Error <- summary(lmm)$tTable[, "Std.Error"]
rates$p_value   <- summary(lmm)$tTable[, "p-value"]
write.csv(rates, "../../data_S1/S7_SBS40b_PT_rates_by_country.csv", row.names = FALSE)
country_order <- c("Japan", "Thailand", "Serbia", "Brazil", "Canada",
"Romania", "UK", "Russia", "Lithuania", "Czechia") #ordered by cancer incidence
rates$Country <- factor(rates$Country, levels = country_order)
avg_rate <- mean(rates$Rate, na.rm = TRUE)
ggplot(rates, aes(x = Country, y = Rate)) +
  geom_point(shape = 21, size = 3, fill=cols['SBS40b'],color = 'black') +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  coord_flip() + theme_bw()+
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(),
        legend.position = "none",
        axis.line = element_line(size = 1, colour = "black"),
         title=element_text(size=6),
              axis.text.y=element_text(size=6,color="black"),
              axis.text.x=element_text(size=6,color="black"),
              legend.text=element_text(size=6))+
  geom_hline(yintercept = avg_rate,
  linetype = "dashed",
  linewidth = 0.5)+
  labs( x = "Country",
    y = "SBS40b mutations per year",
    title = "Proximal tubules")
dev.off()

# Tobacco associated signatures
pdf("../../figure/mutation rate/country/PT/SBSB_PT_rates.pdf",width = 3, height = 3)
lmm <- lme(
  SBSB ~ 0+Age:Country,
  random = ~ 0 + Age | Tobacco/Patient,
  data   = df[df$Structure=='Proximal tubules',],
  method = "ML"
)
summary(lmm)

ci <- intervals(lmm, which = "fixed")$fixed
rates <- data.frame(
  Country = gsub("Age:Country", "", names(fixef(lmm))),
  Rate = as.numeric(fixef(lmm))
)
rates$lower <- ci[, "lower"]
rates$upper <- ci[, "upper"]
rates$Std.Error <- summary(lmm)$tTable[, "Std.Error"]
rates$p_value   <- summary(lmm)$tTable[, "p-value"]
# write.csv(rates, "../../data_S1/S7_SBSB_PT_rates_by_country.csv", row.names = FALSE)

country_order <- c("Japan", "Thailand", "Serbia", "Brazil", "Canada",
"Romania", "UK", "Russia", "Lithuania", "Czechia") #ordered by cancer incidence
rates$Country <- factor(rates$Country, levels = country_order)
avg_rate <- mean(rates$Rate, na.rm = TRUE)
ggplot(rates, aes(x = Country, y = Rate)) +
  geom_point(shape = 21, size = 3, fill=cols['SBSB'],color = 'black') +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  coord_flip() + theme_bw()+
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(),
        legend.position = "none",
        axis.line = element_line(size = 1, colour = "black"),
         title=element_text(size=6),
              axis.text.y=element_text(size=6,color="black"),
              axis.text.x=element_text(size=6,color="black"),
              legend.text=element_text(size=6))+
  geom_hline(yintercept = avg_rate,
  linetype = "dashed",
  linewidth = 0.5)+
  labs( x = "Country",
    y = "SBSB mutations per year",
    title = "Proximal tubules")
dev.off()

# Other signatures
pdf("../../figure/mutation rate/country/PT/SBS40c_PT_rates.pdf",width = 3, height = 3)
lmm <- lme(
  SBS40c ~ 0+Age:Country,
  random = ~ 0 + Age | Patient,
  data   = df[df$Structure=='Proximal tubules',],
  method = "ML"
)
summary(lmm)

ci <- intervals(lmm, which = "fixed")$fixed
rates <- data.frame(
  Country = gsub("Age:Country", "", names(fixef(lmm))),
  Rate = as.numeric(fixef(lmm))
)
rates$lower <- ci[, "lower"]
rates$upper <- ci[, "upper"]
rates <- rates[rates$Country!="(Intercept)",]
country_order <- c("Japan", "Thailand", "Serbia", "Brazil", "Canada",
"Romania", "UK", "Russia", "Lithuania", "Czechia") #ordered by cancer incidence
rates$Country <- factor(rates$Country, levels = country_order)
avg_rate <- mean(rates$Rate, na.rm = TRUE)
ggplot(rates, aes(x = reorder(Country, Rate), y = Rate)) +
  geom_point(shape = 21, size = 3, fill=cols['SBS40c'],color = 'black') +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  coord_flip() + theme_bw()+
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(),
        legend.position = "none",
        axis.line = element_line(size = 1, colour = "black"),
         title=element_text(size=6),
              axis.text.y=element_text(size=6,color="black"),
              axis.text.x=element_text(size=6,color="black"),
              legend.text=element_text(size=6))+
  geom_hline(yintercept = avg_rate,
  linetype = "dashed",
  linewidth = 0.5)+
  labs(x = "Country",
    y = "SBS40c mutations per year",
    title = "Proximal tubules")
dev.off()

#SBS40b vs ASR in PT
lmm <- lme(
  SBS40b ~ 0+Age:Country,
  random = ~ 0 + Age | Sex/Patient,
  data   = df[df$Structure=='Proximal tubules',],
  method = "ML"
)
rates <- data.frame(
  Country = gsub("Age:Country", "", names(fixef(lmm))),
  Rate = as.numeric(fixef(lmm))
)
country_df <- rates %>%
  left_join(
    df %>%
      distinct(Country, ASR),
    by = "Country"
  )

cor.test(country_df$Rate, country_df$ASR, method = "spearman")


# -----------Cancer SBS----------------
metadata_path <- "../data/bulk_normal_kidney_metadata.csv"
weights_path  <- "../data/SBS96_kidney_signature_weights_table.csv"

metadata <- read.csv(metadata_path, header = TRUE, stringsAsFactors = FALSE)
metadata <- metadata %>%
  filter(!is.na(Tumor_SBS_burden)) %>%
  mutate(
    Country = recode(Country, "United Kingdom" = "UK", "Czech Republic" = "Czechia"),
    Country = factor(Country)
  )

df_cancer_SBS <- metadata %>%
  transmute(
    normal = Normal_Kidney,
    cancer = Tumor,
    type = "Cancer",
    Country = Country,
    Age = Age,
    Sex = Sex,
    Alcohol = Alcohol,
    Tobacco = Tobacco,
    Burden_SBS = Tumor_SBS_burden
  )
df_cancer_SBS$Country <- factor(df_cancer_SBS$Country)
df_cancer_SBS <- df_cancer_SBS[df_cancer_SBS$cancer!='PD47592a',]
weights_raw <- read.csv(path.expand(weights_path), check.names = FALSE, stringsAsFactors = FALSE)
# keep only samples present in metadata
weights_raw <- weights_raw[weights_raw$Sample %in% df_cancer_SBS$cancer, , drop = FALSE]
rownames(weights_raw) <- weights_raw$Sample

# align and subset weights to the same order as df_cancer_SBS
weights <- weights_raw[df_cancer_SBS$cancer, , drop = FALSE]

# ensure numeric signature columns (exclude "Sample" if present)
sig_cols <- setdiff(colnames(weights), "Sample")
weights[sig_cols] <- lapply(weights[sig_cols], function(x) as.numeric(as.character(x)))

# zero small contributions, but preserve SBS2/SBS13 when their sum > 0.05
weights[sig_cols] <- lapply(weights[sig_cols], function(x) { x[x < 0.05] <- 0; x })
rows_restore <- rownames(weights)[(weights_raw[rownames(weights), "SBS2"] + weights_raw[rownames(weights), "SBS13"]) > 0.05]
if (length(rows_restore) > 0) {
  weights[rows_restore, "SBS2"]  <- weights_raw[rows_restore, "SBS2"]
  weights[rows_restore, "SBS13"] <- weights_raw[rows_restore, "SBS13"]
}

# compute Unassigned and clamp negatives to zero
weights$Unassigned <- 1 - rowSums(weights[, sig_cols, drop = FALSE], na.rm = TRUE)
weights[sig_cols] <- lapply(weights[sig_cols], function(x) pmax(x, 0))
weights$Unassigned <- pmax(weights$Unassigned, 0)

# add signature Burden_SBSs to df_cancer_SBS
df_cancer_SBS <- df_cancer_SBS %>%
  mutate(
    SBS1 = weights$SBS1 * Burden_SBS,
    SBS5 = weights$SBS5 * Burden_SBS,
    SBS12 = weights$SBS12 * Burden_SBS,
    APOBEC = (weights$SBS2+weights$SBS13) * Burden_SBS,
    SBS18 = weights$SBS18 * Burden_SBS,
    SBS22a = weights$SBS22a * Burden_SBS,
    SBS22b = weights$SBS22b * Burden_SBS,
    SBS22c = weights$SBS22c * Burden_SBS,
    SBS40a = weights$SBS40a * Burden_SBS,
    SBS40b = weights$SBS40b * Burden_SBS,
    SBS40c = weights$SBS40c * Burden_SBS,
    SBSB = weights$SBSB * Burden_SBS,
    SBSC = weights$SBSC * Burden_SBS,
    SBSD = weights$SBSD * Burden_SBS,
    MSI = (weights$SBS21+weights$SBS44) * Burden_SBS,
    Unassigned = weights$Unassigned * Burden_SBS
  )
# write.csv(df_cancer_SBS, "../../data_S1/S7_cancer_SBS_samples.csv", row.names = FALSE)

# Sex associated signatures
pdf("../../figure/mutation rate/country/cancer/SBS40b_cancer_rates.pdf",width = 3, height = 3)
lmm <- lme(
  SBS40b ~  0+Age:Country,
  random = ~ 0 + Age | Sex,
  data   = df_cancer_SBS,
  method = "ML"
)


ci <- intervals(lmm, which = "fixed")$fixed
rates <- data.frame(
  Country = gsub("Age:Country", "", names(fixef(lmm))),
  Rate = as.numeric(fixef(lmm))
)
rates$lower <- ci[, "lower"]
rates$upper <- ci[, "upper"]
rates <- rates[rates$Country!="(Intercept)",]
rates$Std.Error <- summary(lmm)$tTable[, "Std.Error"]
rates$p_value   <- summary(lmm)$tTable[, "p-value"]
write.csv(rates, "../../data_S1/S7_SBS40b_cancer_rates_by_country.csv", row.names = FALSE)
country_order <- c("Japan", "Thailand", "Serbia", "Brazil", "Canada",
"Romania", "UK", "Russia", "Lithuania", "Czechia") #ordered by cancer incidence
rates$Country <- factor(rates$Country, levels = country_order)
avg_rate <- mean(rates$Rate, na.rm = TRUE)
ggplot(rates, aes(x = Country, y = Rate)) +
  geom_point(shape = 21, size = 3, fill=cols['SBS40b'],color = 'black') +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  coord_flip() + theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position = "none",
        axis.line = element_line(size = 1, colour = "black"),
         title=element_text(size=6),
              axis.text.y=element_text(size=6,color="black"),
              axis.text.x=element_text(size=6,color="black"),
              legend.text=element_text(size=6))+
  geom_hline(yintercept = avg_rate,
  linetype = "dashed",
  linewidth = 0.5)+
  labs(x = "Country",
    y = "SBS40b mutations per year",
    title = "Kidney cancer")
dev.off()

# Tobacco associated signatures
pdf("../../figure/mutation rate/country/cancer/SBSB_cancer_rates.pdf",width = 3, height = 3)
lmm <- lme(
  SBSB ~  0+Age:Country,
  random = ~ 0 + Age | Tobacco,
  data   = df_cancer_SBS,
  method = "ML"
)


ci <- intervals(lmm, which = "fixed")$fixed
rates <- data.frame(
  Country = gsub("Age:Country", "", names(fixef(lmm))),
  Rate = as.numeric(fixef(lmm))
)
rates$lower <- ci[, "lower"]
rates$upper <- ci[, "upper"]
rates <- rates[rates$Country!="(Intercept)",]
rates$Std.Error <- summary(lmm)$tTable[, "Std.Error"]
rates$p_value   <- summary(lmm)$tTable[, "p-value"]
country_order <- c("Japan", "Thailand", "Serbia", "Brazil", "Canada",
"Romania", "UK", "Russia", "Lithuania", "Czechia") #ordered by cancer incidence
rates$Country <- factor(rates$Country, levels = country_order)
avg_rate <- mean(rates$Rate, na.rm = TRUE)
ggplot(rates, aes(x = reorder(Country, Rate), y = Rate)) +
  geom_point(shape = 21, size = 3, fill=cols['SBSB'],color = 'black') +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  coord_flip() + theme_bw()+
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(),
        legend.position = "none",
        axis.line = element_line(size = 1, colour = "black"),
         title=element_text(size=6),
              axis.text.y=element_text(size=6,color="black"),
              axis.text.x=element_text(size=6,color="black"),
              legend.text=element_text(size=6))+
  geom_hline(yintercept = avg_rate,
  linetype = "dashed",
  linewidth = 0.5)+
  labs(x = "Country",
    y = "SBSB mutations per year",
    title = "Kidney cancer")
dev.off()

#Other signatures
pdf("../../figure/mutation rate/country/cancer/SBS40c_cancer_rates.pdf",width = 3, height = 3)
lmm <- lm(
  SBS40c ~  0+Age:Country,
  data   = df_cancer_SBS
)
summary(lmm)

rates <- data.frame(
  Country = gsub("Age:Country", "", names(coef(lmm))),
  Rate = as.numeric(coef(lmm))
)
ci <- confint(lmm, level = 0.95)
rates$lower <- ci[, "2.5 %"]
rates$upper <- ci[, "97.5 %"]
rates$Std.Error <- summary(lmm)$coefficients[, "Std. Error"]
rates$p_value   <- summary(lmm)$coefficients[, "Pr(>|t|)"]
# write.csv(rates, "../../data_S1/S7_SBS40c_cancer_rates_by_country.csv", row.names = FALSE)
country_order <- c("Japan", "Thailand", "Serbia", "Brazil", "Canada",
"Romania", "UK", "Russia", "Lithuania", "Czechia") #ordered by cancer incidence
rates$Country <- factor(rates$Country, levels = country_order)
avg_rate <- mean(rates$Rate, na.rm = TRUE)
ggplot(rates, aes(x = reorder(Country, Rate), y = Rate)) +
  geom_point(shape = 21, size = 3, fill=cols['SBS40c'],color = 'black') +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  coord_flip() + theme_bw()+
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(),
        legend.position = "none",
        axis.line = element_line(size = 1, colour = "black"),
         title=element_text(size=6),
              axis.text.y=element_text(size=6,color="black"),
              axis.text.x=element_text(size=6,color="black"),
              legend.text=element_text(size=6))+
  geom_hline(yintercept = avg_rate,
  linetype = "dashed",
  linewidth = 0.5)+
  labs(x = "Country",
    y = "SBS40c mutations per year",
    title = "Kidney cancer")
dev.off()
