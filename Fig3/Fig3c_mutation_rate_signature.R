library(ggplot2)
library(gridExtra)
library(ggpubr)
library(nlme)
library(tidyverse)
library(RColorBrewer)
# library("cowplot")
library(dplyr)

# load SBS for LCM
metadata_path = "../data/LCM_normal_kidney_metadata.csv"
weights_path  <- "../data/SBS96_kidney_signature_weights_table.csv"

LCM_metadata <- read.csv(metadata_path, header = TRUE, stringsAsFactors = FALSE)

df = data.frame(Patient = LCM_metadata$Patient, Sample = LCM_metadata$Sample, Structure = LCM_metadata$Structure, Country = LCM_metadata$Country, Sex = LCM_metadata$Sex, Age = LCM_metadata$Age, Burden= LCM_metadata$Normal_SBS_burden)
df$Structure[df$Structure=='Distal_tubules'] ='Distal tubules'
df$Structure[df$Structure=='Proximal_tubules'] ='Proximal tubules'
df$Structure = factor(df$Structure, levels = c('Proximal tubules','Medulla','Glomeruli','Distal tubules'))

weights_raw <- read.csv(path.expand(weights_path), check.names = FALSE, stringsAsFactors = FALSE)
weights_raw$Sample <- sub(" .*", "", weights_raw$Sample)
# keep only samples present in metadata
weights_raw <- weights_raw[weights_raw$Sample %in% df$Sample, , drop = FALSE]
rownames(weights_raw) <- weights_raw$Sample

# align and subset weights to the same order as df_normal_SBS
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
    SBS1 = weights$SBS1 * Burden,
    SBS5 = weights$SBS5 * Burden,
    SBS12 = weights$SBS12 * Burden,
    SBS2 = weights$SBS2 * Burden,
    SBS13 = weights$SBS13 * Burden,
    SBS18 = weights$SBS18 * Burden,
    SBS22a = weights$SBS22a * Burden,
    SBS22b = weights$SBS22b * Burden,
    SBS22c = weights$SBS22c * Burden,
    SBS40a = weights$SBS40a * Burden,
    SBS40b = weights$SBS40b * Burden,
    SBS40c = weights$SBS40c * Burden,
    SBSB = weights$SBSB * Burden,
    SBSC = weights$SBSC * Burden,
    SBSD = weights$SBSD * Burden,
    SBS21 = weights$SBS21 * Burden,
    SBS44 = weights$SBS44 * Burden,
    Unassigned = weights$Unassigned * Burden
  )

head(df)


# merge samples from the same patient and structure
sbs_cols <- grep("^SBS", colnames(df), value = TRUE)

df_avg <- df %>%
  group_by(Patient, Structure) %>%
  summarise(
    Country  = first(Country),
    Sex      = first(Sex),
    Age      = first(Age),
    # n_samples = n(),
    across(all_of(sbs_cols), mean, na.rm = TRUE)
  )

# ---------save regression summary for all signatures in one table------
all_reg_summaries <- lapply(sbs_cols, function(sig) {
  fit <- tryCatch(
    lme(
      as.formula(paste0(sig, " ~ Age:Structure")),
      random = ~ 1 | Country/Patient,
      data = df_avg,
      method = "REML"
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  data.frame(
    Signature = sig,
    Term      = rownames(summary(fit)$tTable),
    Estimate  = summary(fit)$tTable[, "Value"],
    Std.Error = summary(fit)$tTable[, "Std.Error"],
    CI_lower  = intervals(fit, which = "fixed")$fixed[, "lower"],
    CI_upper  = intervals(fit, which = "fixed")$fixed[, "upper"],
    p_value   = summary(fit)$tTable[, "p-value"]
  )
})
all_reg_summaries <- do.call(rbind, all_reg_summaries)
write.csv(all_reg_summaries, "../../data_S1/Fig3c_LCM_all_SBS_regression.csv", row.names = FALSE)


# ---------repeat the following for all signatures------
# r modelling SBS burden
lmm1 <- lme(
  SBS1 ~  Age:Structure,
  random = ~ 1 | Country/Patient,
  data = df_avg,
  method = "REML" 
)

summary(lmm1)

fixed.m1 <- data.frame(fixef(lmm1))
intervals(lmm1, which = "fixed")

# plotting
structure_cols <- c(
  '#C77CFF','#00BFC4','#7CAE00','#F8766D'
)
names(structure_cols) <- c('Distal tubules','Glomeruli','Medulla','Proximal tubules')
pdf("../../figure/mutation rate/merged/SBS1.pdf",width = 2.5, height = 2.5)
p_sbs=ggplot(data = df_avg, mapping = aes(x = Age, y = SBS1 ))+geom_point(data=df_avg, aes(colour = Structure,fill=Structure),alpha = 0.8) +theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.25,colour="black"),axis.ticks = element_line(linewidth = 0.25))+

geom_abline(intercept = fixed.m1[1,], slope =  fixed.m1['Age:StructureDistal tubules',],colour='#C77CFF')+
geom_ribbon(aes(ymin = fixed.m1[1,]+Age*intervals(lmm1, which = "fixed")[["fixed"]]['Age:StructureDistal tubules','lower'], ymax = fixed.m1[1,]+Age*intervals(lmm1, which = "fixed")[["fixed"]]['Age:StructureDistal tubules','upper']),fill='#C77CFF', alpha = 0.1)+

geom_abline(intercept = fixed.m1[1,], slope =  fixed.m1['Age:StructureGlomeruli',],colour='#00BFC4')+
geom_ribbon(aes(ymin = fixed.m1[1,]+
                  Age*intervals(lmm1, which = "fixed")[["fixed"]]['Age:StructureGlomeruli','lower'], ymax = fixed.m1[1,]+Age*intervals(lmm1, which = "fixed")[["fixed"]]['Age:StructureGlomeruli','upper']), fill='#00BFC4',alpha = 0.1)+

geom_abline(intercept = fixed.m1[1,], slope =  fixed.m1['Age:StructureMedulla',],colour='#7CAE00')+
geom_ribbon(aes(ymin = fixed.m1[1,]+Age*intervals(lmm1, which = "fixed")[["fixed"]]['Age:StructureMedulla','lower'], ymax = fixed.m1[1,]+Age*intervals(lmm1, which = "fixed")[["fixed"]]['Age:StructureMedulla','upper']), fill='#7CAE00',alpha = 0.1)+

# geom_abline(intercept = fixed.m1[1,], slope =  fixed.m1['Age:StructureProximal tubules',],colour='#F8766D')+
# geom_ribbon(aes(ymin = fixed.m1[1,]+Age*intervals(lmm1, which = "fixed")[["fixed"]]['Age:StructureProximal tubules','lower'], ymax = fixed.m1[1,]+Age*intervals(lmm1, which = "fixed")[["fixed"]]['Age:StructureProximal tubules','upper']), fill='#F8766D',alpha = 0.1)+
        labs(y='Substitutions / genome',x="Age (yrs)",  fill = "Structure", color = "Structure", title='SBS1')+
        scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05), add = c(0, 0)))+
        theme(legend.key.size = unit(0.25, "cm"),
              title=element_text(size=12),
              axis.text.y=element_text(size=12,color="black"),
              axis.text.x=element_text(size=12,color="black"),
              legend.text=element_text(size=12),
              axis.title.x = element_blank(),
              axis.title.y  = element_blank(),
              legend.position = "none")

p_sbs
dev.off()
