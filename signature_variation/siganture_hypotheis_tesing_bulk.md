    metadata_path <- "../data/bulk_normal_kidney_metadata.csv"
    weights_path  <- "../data/SBS96_kidney_signature_weights_table.csv"

    metadata <- read.csv(metadata_path, header = TRUE, stringsAsFactors = FALSE)
    metadata$Country[metadata$Country=='United Kingdom']='UK'
    metadata$Country[metadata$Country=='Czech Republic']='Czechia'

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
    df_normal <- df_normal[df_normal$normal!='PD47592c_ds0003',] #remove one repetitive sample for PD47592

    weights_path <- "../data/ID83_kidney_signature_ID_weights_table.csv"

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

    ASR_df <- data.frame(
      Country = c('Brazil','Czechia', "Japan", 'Lithuania', "Romania", "Russia", "Serbia", "Thailand", "UK", "Canada"),  
      ASR = c(4.5,14.4,7.6,14.5,7.7,10.3,7.4,1.8,10.3,10.4) #age-standardized kidney cancer incidence per 100k people from GLOBOCAN 2020
    )

    df_normal <- merge(df_normal, ASR_df, by = "Country", all.x = TRUE)
    df_normal$Country = as.factor(df_normal$Country)

    metadata_path <- "../data/bulk_normal_kidney_metadata.csv"
    weights_path  <- "../data/SBS96_kidney_signature_weights_table.csv"

    metadata <- read.csv(metadata_path, header = TRUE, stringsAsFactors = FALSE)
    metadata <- dplyr::filter(metadata, !is.na(Tumor_SBS_burden))
    metadata$Country[metadata$Country=='United Kingdom']='UK'
    metadata$Country[metadata$Country=='Czech Republic']='Czechia'

    df_cancer <- metadata %>%
      transmute(
        normal = Normal_Kidney,
        cancer = Tumor,
        type = "Cancer",
        Country = Country,
        Age = Age,
        Sex = Sex,
        Alcohol = Alcohol,
        Tobacco = Tobacco,
        Burden_SBS = Tumor_SBS_burden,
        Burden_ID = Tumor_ID_burden
      )
    df_cancer$Country <- factor(df_cancer$Country)
    df_cancer <- df_cancer[df_cancer$cancer!='PD47592a',]
    weights_raw <- read.csv(path.expand(weights_path), check.names = FALSE, stringsAsFactors = FALSE)
    # keep only samples present in metadata
    weights_raw <- weights_raw[weights_raw$Sample %in% df_cancer$cancer, , drop = FALSE]
    rownames(weights_raw) <- weights_raw$Sample

    # align and subset weights to the same order as df_cancer
    weights <- weights_raw[df_cancer$cancer, , drop = FALSE]

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

    # add signature Burden_SBSs to df_cancer
    df_cancer <- df_cancer %>%
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

    weights_path  <- "../data/ID83_kidney_signature_ID_weights_table.csv"

    weights_raw <- read.csv(path.expand(weights_path), check.names = FALSE, stringsAsFactors = FALSE)
    # keep only samples present in metadata
    weights_raw <- weights_raw[weights_raw$Sample %in% metadata$Tumor, , drop = FALSE]
    rownames(weights_raw) <- weights_raw$Sample
    weights <- weights_raw[df_cancer$cancer, ,drop = FALSE ]

    # ensure numeric signature columns (exclude "Sample" if present)
    sig_cols <- setdiff(colnames(weights), "Sample")
    weights[sig_cols] <- lapply(weights[sig_cols], function(x) as.numeric(as.character(x)))

    # zero small contributions, but preserve SBS2/SBS13 when their sum > 0.05 (based on original raw values)
    weights[sig_cols] <- lapply(weights[sig_cols], function(x) { x[x < 0.05] <- 0; x })

    # compute Unassigned and clamp negatives to zero
    weights$Unassigned <- 1 - rowSums(weights[, sig_cols, drop = FALSE], na.rm = TRUE)
    weights[sig_cols] <- lapply(weights[sig_cols], function(x) pmax(x, 0))
    weights$Unassigned <- pmax(weights$Unassigned, 0)

    df_cancer<- df_cancer %>%
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

    ASR_df <- data.frame(
      Country = c('Brazil','Czechia', "Japan", 'Lithuania', "Romania", "Russia", "Serbia", "Thailand", "UK", "Canada"),  
      ASR = c(4.5,14.4,7.6,14.5,7.7,10.3,7.4,1.8,10.3,10.4) #age-standardized kidney cancer incidence per 100k people from GLOBOCAN 2020
    )

    df_cancer <- merge(df_cancer, ASR_df, by = "Country", all.x = TRUE)
    df_cancer$Country = as.factor(df_cancer$Country)

## Smoking

    df_normal$Tobacco[df_normal$Tobacco=='Current smoker'] = 'Smoker'
    df_normal$Tobacco[df_normal$Tobacco=='Ex-smoker'] = 'Smoker'
    df_normal$Tobacco <- factor(df_normal$Tobacco)
    df_normal$Tobacco <- relevel(df_normal$Tobacco, ref = "Never")
    lmm <- lme(SBSB ~ Age + Tobacco,
               random = ~1|Country, weights = varIdent(form = ~1), data = df_normal[df_normal$Tobacco!='Missing',], method = "REML")
    summary(lmm)

    ## Linear mixed-effects model fit by REML
    ##   Data: df_normal[df_normal$Tobacco != "Missing", ] 
    ##        AIC     BIC    logLik
    ##   4163.855 4182.29 -2076.928
    ## 
    ## Random effects:
    ##  Formula: ~1 | Country
    ##         (Intercept) Residual
    ## StdDev:    97.28809  261.084
    ## 
    ## Fixed effects:  SBSB ~ Age + Tobacco 
    ##                    Value Std.Error  DF   t-value p-value
    ## (Intercept)   -129.91539  86.86789 286 -1.495551  0.1359
    ## Age              6.65172   1.24324 286  5.350306  0.0000
    ## TobaccoSmoker  121.21035  31.22565 286  3.881756  0.0001
    ##  Correlation: 
    ##               (Intr) Age   
    ## Age           -0.893       
    ## TobaccoSmoker -0.299  0.150
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -2.23704404 -0.77446871  0.09270948  0.68718443  2.62524285 
    ## 
    ## Number of Observations: 298
    ## Number of Groups: 10

    anova(lmm)["Tobacco", "p-value"]

    ## [1] 0.0001288199

    df_cancer$Tobacco[df_cancer$Tobacco=='Current smoker'] = 'Smoker'
    df_cancer$Tobacco[df_cancer$Tobacco=='Ex-smoker'] = 'Smoker'

    df_cancer$Tobacco <- factor(df_cancer$Tobacco)
    df_cancer$Tobacco <- relevel(df_cancer$Tobacco, ref = "Never")
    lmm <- lme(SBSB ~ Age + Tobacco,
               random = ~1|Country, weights = varIdent(form = ~1), data = df_cancer[df_cancer$Tobacco!='Missing',], method = "REML")
    summary(lmm)

    ## Linear mixed-effects model fit by REML
    ##   Data: df_cancer[df_cancer$Tobacco != "Missing", ] 
    ##        AIC      BIC    logLik
    ##   4695.534 4713.831 -2342.767
    ## 
    ## Random effects:
    ##  Formula: ~1 | Country
    ##         (Intercept) Residual
    ## StdDev:    304.7027 801.0045
    ## 
    ## Fixed effects:  SBSB ~ Age + Tobacco 
    ##                   Value Std.Error  DF   t-value p-value
    ## (Intercept)   -295.4845 271.79399 278 -1.087164  0.2779
    ## Age             23.9786   3.87373 278  6.190058  0.0000
    ## TobaccoSmoker  402.9495  97.23395 278  4.144124  0.0000
    ##  Correlation: 
    ##               (Intr) Age   
    ## Age           -0.893       
    ## TobaccoSmoker -0.303  0.154
    ## 
    ## Standardized Within-Group Residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -2.8953160 -0.6163033 -0.1246945  0.4356544  4.5383880 
    ## 
    ## Number of Observations: 290
    ## Number of Groups: 10

    anova(lmm)["Tobacco", "p-value"]

    ## [1] 4.532184e-05

## Sex

    signatures <- c("SBS40b", "ID5", "ID8", "SBS40a", "SBS40c", "SBS1", "SBS5", "SBSB", "ID3","SBSC", "ID9", "ID11","ID21", "SBS22a", "SBS22b", "SBS22c", "ID23","SBS12")

    with_tobacco <- c("SBSB", "ID3")

    pvals <- sapply(signatures, function(sig) {
      print(sig)
      if (sig %in% with_tobacco){
        formula <- as.formula(paste(sig, "~ Age + Tobacco + Sex"))
      } else {
        formula <- as.formula(paste(sig, "~ Age + Sex"))}
      fit <- lme(
        formula,
        random = ~ 1 | Country,
        data   = df_normal,
        method = "ML"
      )
      anova(fit)["Sex", "p-value"]
    })

    ## [1] "SBS40b"
    ## [1] "ID5"
    ## [1] "ID8"
    ## [1] "SBS40a"
    ## [1] "SBS40c"
    ## [1] "SBS1"
    ## [1] "SBS5"
    ## [1] "SBSB"
    ## [1] "ID3"
    ## [1] "SBSC"
    ## [1] "ID9"
    ## [1] "ID11"
    ## [1] "ID21"
    ## [1] "SBS22a"
    ## [1] "SBS22b"
    ## [1] "SBS22c"
    ## [1] "ID23"
    ## [1] "SBS12"

    pvals_adj <- p.adjust(pvals, method = "BH")

    results <- data.frame(
      Signature = names(pvals),
      p_value   = unname(pvals),
      p_adj     = unname(pvals_adj)
    )

    results <- results[order(results$p_adj), ]
    results

    ##    Signature      p_value       p_adj
    ## 2        ID5 0.0001499034 0.002698262
    ## 1     SBS40b 0.0010024136 0.009021722
    ## 3        ID8 0.0089224970 0.045745157
    ## 6       SBS1 0.0103395151 0.045745157
    ## 12      ID11 0.0127069879 0.045745157
    ## 17      ID23 0.0555441479 0.166632444
    ## 13      ID21 0.0904199587 0.232508465
    ## 5     SBS40c 0.1689658149 0.357313679
    ## 14    SBS22a 0.1786568396 0.357313679
    ## 15    SBS22b 0.2046912214 0.368444198
    ## 7       SBS5 0.2571441784 0.420781383
    ## 11       ID9 0.3250002534 0.487500380
    ## 4     SBS40a 0.3699533016 0.512243033
    ## 18     SBS12 0.4243891929 0.545643248
    ## 8       SBSB 0.5461559071 0.602696924
    ## 9        ID3 0.6026969242 0.602696924
    ## 10      SBSC 0.5997381052 0.602696924
    ## 16    SBS22c 0.5256473293 0.602696924

## Age

    signatures <- c("SBS40b", "ID5", "ID8", "SBS40a", "SBS40c", "SBS1", "SBS5", "SBSB", "ID3","SBSC", "ID9", "ID11","ID21", "SBS22a", "SBS22b", "SBS22c", "ID23","SBS12")

    baseline        <- c("SBS40a","SBS40c","SBS1","SBS5","SBSC","ID8","ID9","ID11","ID21")
    with_tobacco    <- c("SBSB","ID3")
    with_sex    <- c("SBS40b","ID5")
    region_specific <- c("SBS22a","SBS22b","SBS22c","ID23","SBS12")

    pvals <- sapply(signatures, function(sig) {
      formula <- if (sig %in% with_tobacco) {
        as.formula(paste0(sig, " ~ Age + Tobacco"))
      } else if (sig %in% with_sex) {
        as.formula(paste0(sig, " ~ Age + Sex"))
      } else {
        as.formula(paste0(sig, " ~ Age"))
      }
      fit <- lme(
        formula,
        random = ~ 1 | Country,
        data   = df_normal,
        method = "ML"
      )
      anova(fit)["Age", "p-value"]
    })

    pvals_adj <- p.adjust(pvals, method = "BH")

    results <- data.frame(
      Signature = names(pvals),
      p_value   = unname(pvals),
      p_adj     = unname(pvals_adj)
    )

    results <- results[order(results$p_adj), ]
    results

    ##    Signature      p_value        p_adj
    ## 2        ID5 1.365574e-14 2.458034e-13
    ## 3        ID8 2.358402e-11 2.122562e-10
    ## 1     SBS40b 1.230908e-08 5.907684e-08
    ## 5     SBS40c 1.312819e-08 5.907684e-08
    ## 11       ID9 1.609448e-07 5.794013e-07
    ## 8       SBSB 2.398235e-06 7.194704e-06
    ## 10      SBSC 3.559598e-05 9.153251e-05
    ## 15    SBS22b 6.589913e-05 1.482731e-04
    ## 14    SBS22a 1.475877e-04 2.951754e-04
    ## 13      ID21 3.729866e-04 6.713759e-04
    ## 7       SBS5 1.920707e-03 3.142974e-03
    ## 9        ID3 3.145264e-03 4.397306e-03
    ## 18     SBS12 3.175832e-03 4.397306e-03
    ## 17      ID23 5.801512e-03 7.459087e-03
    ## 16    SBS22c 1.846426e-02 2.215712e-02
    ## 12      ID11 2.147818e-02 2.416295e-02
    ## 6       SBS1 3.425549e-02 3.627052e-02
    ## 4     SBS40a 2.395899e-01 2.395899e-01

## ASR

    df_normal <- df_normal %>% mutate(Age_c = Age - median(Age))
    signatures <- c("SBS40b", "ID5", "ID8", "SBS40a")

    extract_country_burden <- function(sig) {
        formula <- if (sig %in% c('SBS40b','ID5')) {
        as.formula(paste(sig, "~  Age_c + Sex"))
      } else {
        as.formula(paste0(sig, " ~  Age_c"))
      }
        
      model <- lme(
      fixed  = formula,
      random = ~ Age_c | Country,
      data   = df_normal,
      method = "REML",
      control = lmeControl(opt = "optim")
    )

      fe <- fixef(model)
      re <- ranef(model)

      tibble(
        Country = rownames(re),
        Signature = sig,
        burden_med_age =
          fe["(Intercept)"] + re$`(Intercept)`
      )
    }

    asr_tbl <- df_normal %>%
      group_by(Country) %>%
      summarise(ASR = first(ASR), .groups = "drop")

    country_burden <- map_dfr(signatures, extract_country_burden) %>%
      left_join(asr_tbl, by = "Country")


    results <- country_burden %>%
      group_by(Signature) %>%
      do({
        ct <- cor.test(
          .$burden_med_age,
          .$ASR,
          method = "spearman",
          exact = FALSE
        )
        tibble(
          rho = unname(ct$estimate),
          p_value = ct$p.value
        )
      }) %>%
      ungroup() %>%
      mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
      arrange(p_adj)

    results

    ## # A tibble: 4 × 4
    ##   Signature    rho p_value   p_adj
    ##   <chr>      <dbl>   <dbl>   <dbl>
    ## 1 SBS40b     0.839 0.00241 0.00966
    ## 2 ID8        0.754 0.0118  0.0236 
    ## 3 ID5        0.462 0.179   0.238  
    ## 4 SBS40a    -0.219 0.544   0.544

## Country

    signatures <- c("SBS40b", "ID5", "ID8", "SBS40a", "SBS40c", "SBS1", "SBS5", "SBSB", "ID3","SBSC", "ID9", "ID11","ID21", "SBS22a", "SBS22b", "SBS22c", "ID23","SBS12")
    with_tobacco    <- c("SBSB","ID3")
    with_sex    <- c("SBS40b","ID5")
    region_specific <- c("SBS22a","SBS22b","SBS22c","ID23","SBS12")

    fit_one_signature <- function(sig, df) {

      # build formula based on your logic
      form <- if (sig %in% region_specific) {
        as.formula(paste0(sig, " ~ Country"))
      } else if (sig %in% with_tobacco) {
        as.formula(paste0(sig, " ~ Age + Tobacco + Country"))
      } else if (sig %in% with_sex) {
        as.formula(paste0(sig, " ~ Age + Sex+ Country"))
      } else {
        as.formula(paste0(sig, " ~ Age + Country"))
      }

      model <- lm(form, data = df)

      # adjusted country means (for non-region-specific sigs, this is adjusted for Age/Tobacco)
      emm   <- emmeans(model, ~ Country)

      # Country vs grand mean (one test per country)
      contr <- contrast(emm, method = "eff")

      res <- as.data.frame(summary(contr, infer = c(TRUE, TRUE), adjust = "none")) |>
        mutate(
          Signature = sig,
          Country   = sub(" effect$", "", contrast)  # "France effect" -> "France"
        )

      res
    }

    # 1) run all signatures
    all_tests <- map_dfr(signatures, ~fit_one_signature(.x, df_normal))

    # 2) Global BH across ALL (Country × Signature) tests
    all_tests_normal <- all_tests |>
      mutate(p_adj_global = p.adjust(p.value, method = "BH")) |>
      arrange(p_adj_global, p.value)

    all_tests_normal

    ##             contrast     estimate         SE  df      lower.CL     upper.CL
    ## 1       Japan effect  230.0954788  21.857465 289   187.0754764  273.1154813
    ## 2     Romania effect  985.2574350 128.489622 289   732.3633336 1238.1515363
    ## 3       Japan effect -108.2791364  14.146399 287  -136.1229868  -80.4352861
    ## 4     Romania effect  319.2339455  43.739296 289   233.1459824  405.3219086
    ## 5     Czechia effect   93.5771590  13.683053 287    66.6452961  120.5090219
    ## 6     Romania effect 1446.1921667 226.087761 289  1001.2047801 1891.1795533
    ## 7     Czechia effect  266.5057885  45.208101 287   177.5243056  355.4872714
    ## 8      Serbia effect 1940.2254266 346.479866 289  1258.2815251 2622.1693280
    ## 9     Czechia effect   15.8891987   2.950712 288    10.0815032   21.6968941
    ## 10    Romania effect  -87.1065904  16.339136 288  -119.2658529  -54.9473278
    ## 11     Serbia effect 1000.6213169 196.910557 289   613.0606909 1388.1819429
    ## 12      Japan effect -235.8503739  46.738972 287  -327.8450161 -143.8557317
    ## 13    Romania effect   58.4023420  11.886980 289    35.0063112   81.7983728
    ## 14    Czechia effect  171.7860179  36.231518 286   100.4717664  243.1002695
    ## 15     Serbia effect   76.8718080  18.216817 289    41.0173524  112.7262636
    ## 16    Czechia effect   16.4356709   4.232728 288     8.1046669   24.7666748
    ## 17      Japan effect  -11.1863024   3.003842 288   -17.0985690   -5.2740359
    ## 18      Japan effect -225.9699369  60.842624 288  -345.7225287 -106.2173451
    ## 19    Romania effect  261.3313757  76.608082 288   110.5486540  412.1140973
    ## 20      Japan effect   12.9599626   3.902232 288     5.2794526   20.6404726
    ## 21   Thailand effect  -42.6933809  13.943377 288   -70.1372265  -15.2495354
    ## 22    Czechia effect -531.7852323 176.427564 289  -879.0310989 -184.5393657
    ## 23    Czechia effect -295.4257720 100.266865 289  -492.7716610  -98.0798831
    ## 24    Romania effect -263.4079410  91.361574 288  -443.2290054  -83.5868767
    ## 25      Japan effect -508.6451765 178.986343 289  -860.9272480 -156.3631049
    ## 26      Japan effect -281.8136495 101.721064 289  -482.0217017  -81.6055974
    ## 27     Serbia effect  -68.5083851  25.096838 288  -117.9048639  -19.1119062
    ## 28    Czechia effect  -25.0494360   9.276004 289   -43.3065265   -6.7923455
    ## 29     Russia effect -111.8874592  41.877748 289  -194.3115114  -29.4634069
    ## 30     Russia effect   45.9347734  17.250186 287    11.9818514   79.8876954
    ## 31      Japan effect  105.1539625  39.954563 288    26.5139859  183.7939392
    ## 32    Czechia effect    6.6693376   2.538417 288     1.6731354   11.6655397
    ## 33    Romania effect -119.1207735  46.437393 286  -210.5231801  -27.7183669
    ## 34     Serbia effect  171.6938466  67.030542 289    39.7639016  303.6237916
    ## 35    Czechia effect  -86.6729651  34.131955 289  -153.8516993  -19.4942309
    ## 36         UK effect  263.6215331 104.010724 288    58.9039679  468.3390982
    ## 37      Japan effect  -86.8244955  34.626981 289  -154.9775418  -18.6714491
    ## 38     Russia effect -535.1744697 216.465447 289  -961.2231562 -109.1257833
    ## 39     Russia effect -291.4897677 123.021093 289  -533.6206750  -49.3588603
    ## 40     Russia effect  107.8826477  45.830996 286    17.6738069  198.0914884
    ## 41         UK effect  120.5509439  52.974082 286    16.2824152  224.8194726
    ## 42  Lithuania effect   57.1774826  26.082060 288     5.8418537  108.5131114
    ## 43    Romania effect   37.2018885  17.509419 287     2.7387278   71.6650493
    ## 44     Russia effect    7.8328894   3.713655 288     0.5235423   15.1422365
    ## 45         UK effect -538.7594002 256.684538 289 -1043.9675589  -33.5512415
    ## 46     Brazil effect  -83.3912893  40.683279 289  -163.4643815   -3.3181971
    ## 47         UK effect -298.4474223 145.878304 289  -585.5660399  -11.3288047
    ## 48     Serbia effect -143.8230702  71.563276 286  -284.6805830   -2.9655575
    ## 49     Russia effect  -22.6907488  11.381069 289   -45.0910426   -0.2904550
    ## 50     Russia effect  177.8155517  89.706002 288     1.2530439  354.3780595
    ## 51         UK effect -109.9428251  57.272691 288  -222.6689495    2.7832993
    ## 52     Brazil effect   29.1601048  15.271191 288    -0.8971901   59.2173996
    ## 53    Czechia effect  132.7494221  71.276563 288    -7.5396154  273.0384596
    ## 54  Lithuania effect -126.0325118  69.725608 289  -263.2669024   11.2018788
    ## 55     Brazil effect    8.9313808   5.070853 288    -1.0492497   18.9120113
    ## 56     Brazil effect -361.1084100 210.291259 289  -775.0050193   52.7881993
    ## 57      Japan effect   -4.3934356   2.584123 288    -9.4795975    0.6927262
    ## 58  Lithuania effect  -12.8481083   7.568830 286   -27.7457861    2.0495695
    ## 59         UK effect   10.3964217   6.176632 288    -1.7606426   22.5534859
    ## 60     Serbia effect -235.6003346 140.330957 288  -511.8046595   40.6039904
    ## 61     Brazil effect -200.3555284 119.512194 289  -435.5801981   34.8691414
    ## 62  Lithuania effect   14.2953628   8.660640 288    -2.7508134   31.3415390
    ## 63     Canada effect  179.8520509 110.447386 288   -37.5343781  397.2384800
    ## 64    Romania effect   81.6044908  50.307536 288   -17.4125704  180.6215519
    ## 65         UK effect  -79.9121160  49.658596 289  -177.6504849   17.8262529
    ## 66      Japan effect  -14.5703361   9.410536 289   -33.0922143    3.9515422
    ## 67      Japan effect  -58.2642057  37.657092 286  -132.3844067   15.8559953
    ## 68   Thailand effect  -15.0375348   9.720184 288   -34.1691422    4.0940726
    ## 69         UK effect  -20.7519757  13.495662 289   -47.3142244    5.8102730
    ## 70     Canada effect   14.3122424   9.332086 286    -4.0560398   32.6805245
    ## 71  Lithuania effect -542.9566789 360.410620 289 -1252.3191798  166.4058219
    ## 72     Brazil effect   78.9002973  52.490562 288   -24.4134739  182.2140684
    ## 73    Romania effect   -4.7849107   3.253718 288   -11.1889920    1.6191705
    ## 74  Lithuania effect -298.4474223 204.827648 289  -701.5905171  104.6956725
    ## 75     Brazil effect   -6.4501251   4.429858 286   -15.1693844    2.2691342
    ## 76     Russia effect    4.6505118   3.194757 288    -1.6375204   10.9385440
    ## 77     Serbia effect  112.1908536  77.272143 288   -39.8988954  264.2806026
    ## 78  Lithuania effect  -27.4447649  18.949252 289   -64.7408039    9.8512740
    ## 79  Lithuania effect  132.6126963  92.357601 287   -49.1714548  314.3968473
    ## 80     Russia effect   79.5321752  56.993723 287   -32.6465241  191.7108745
    ## 81         UK effect -121.4130141  87.214588 288  -293.0718337   50.2458055
    ## 82         UK effect    5.9510776   4.305843 288    -2.5238345   14.4259898
    ## 83     Brazil effect  -95.8497943  71.600887 288  -236.7771777   45.0775892
    ## 84    Czechia effect  -28.6157567  21.544991 289   -71.0207476   13.7892342
    ## 85    Czechia effect   -5.0498279   3.833213 288   -12.5944916    2.4948359
    ## 86     Canada effect  -44.7093501  34.426554 287  -112.4699006   23.0512003
    ## 87     Brazil effect  -33.1589175  25.680360 289   -83.7031677   17.3853326
    ## 88     Serbia effect   -9.1671798   7.304996 286   -23.5455539    5.2111944
    ## 89    Czechia effect   15.9953305  12.747126 288    -9.0940105   41.0846714
    ## 90     Russia effect  -33.1589175  26.434340 289   -85.1871557   18.8693207
    ## 91   Thailand effect -183.7103888 148.995867 287  -476.9736074  109.5528298
    ## 92     Canada effect -542.9566789 444.321559 289 -1417.4732314  331.5598736
    ## 93   Thailand effect  136.0329383 112.256685 289   -84.9113917  356.9772683
    ## 94      Japan effect  -54.0099781  44.603687 288  -141.8005236   33.7805674
    ## 95     Canada effect -117.8360460  98.935253 288  -312.5638906   76.8917987
    ## 96     Russia effect    6.2235150   5.327152 288    -4.2615733   16.7086033
    ## 97  Lithuania effect   32.3182545  27.953706 287   -22.7020225   87.3385314
    ## 98     Brazil effect   97.8756171  85.390073 288   -70.1921284  265.9433627
    ## 99     Russia effect  -56.3308512  49.395907 288  -153.5536137   40.8919113
    ## 100        UK effect  -72.2283673  63.936956 288  -198.0713338   53.6145992
    ## 101    Brazil effect    5.0947761   4.592229 288    -3.9438106   14.1333627
    ## 102   Czechia effect  -48.4789363  43.814775 288  -134.7167178   37.7588452
    ## 103   Romania effect   -5.3799278   4.913373 288   -15.0506012    4.2907456
    ## 104    Brazil effect  -11.9976015  11.056450 289   -33.7589774    9.7637745
    ## 105    Canada effect -273.5710174 252.515700 289  -770.5740434  223.4320086
    ## 106        UK effect  -33.1589175  31.345817 289   -94.8539568   28.5361217
    ## 107  Thailand effect  207.8455488 196.881717 288  -179.6639710  595.3550687
    ## 108    Serbia effect  -28.2015469  26.976327 287   -81.2980836   24.8949898
    ## 109   Romania effect   -5.6234246   5.425468 288   -16.3020209    5.0551716
    ## 110  Thailand effect -147.4342138 144.333856 288  -431.5171852  136.6487576
    ## 111    Canada effect  143.9128236 150.658144 288  -152.6178354  440.4434827
    ## 112     Japan effect   66.0960713  72.559941 288   -76.7189561  208.9110987
    ## 113    Serbia effect   77.3788938  86.263550 288   -92.4080593  247.1658468
    ## 114    Canada effect   -9.5350139  10.669774 288   -30.5356376   11.4656098
    ## 115    Canada effect   -5.5488520   6.398790 288   -18.1431759    7.0454719
    ## 116    Serbia effect -101.5233105 117.669660 288  -333.1248751  130.0782541
    ## 117    Russia effect  -46.7957838  55.143629 288  -155.3314131   61.7398455
    ## 118        UK effect    3.0946997   3.704200 288    -4.1960364   10.3854358
    ## 119 Lithuania effect   71.0887538  89.649985 288  -105.3634996  247.5410073
    ## 120  Thailand effect    9.5032203  12.192131 286   -14.4944692   33.5009097
    ## 121  Thailand effect   31.9046872  41.991355 288   -50.7441750  114.5535494
    ## 122 Lithuania effect  -33.1589175  44.012645 289  -119.7848886   53.4670535
    ## 123  Thailand effect  -32.5292503  45.096307 287  -121.2906929   56.2321924
    ## 124    Canada effect  -16.7447649  23.361024 289   -62.7240821   29.2345522
    ## 125     Japan effect    9.1248154  12.976645 288   -16.4162745   34.6659052
    ## 126        UK effect   13.0235477  18.601315 288   -23.5882138   49.6353092
    ## 127  Thailand effect  -90.4694201 129.289674 288  -344.9419048  164.0030646
    ## 128    Canada effect -125.6939713 179.672496 288  -479.3316967  227.9437542
    ## 129   Romania effect  -38.2727164  56.161334 288  -148.8114294   72.2659966
    ## 130     Japan effect   -2.8048025   4.308941 288   -11.2858112    5.6762063
    ## 131   Romania effect  -17.6286181  27.609398 289   -71.9696136   36.7123774
    ## 132    Brazil effect   -2.2318891   3.534984 288    -9.1895695    4.7257912
    ## 133        UK effect   12.4699031  19.941296 287   -26.7798346   51.7196407
    ## 134    Serbia effect   -3.1016551   4.997695 288   -12.9382943    6.7349842
    ## 135    Canada effect  -52.2398936  85.959151 289  -221.4252481  116.9454608
    ## 136 Lithuania effect   43.4343162  80.305603 288  -114.6259968  201.4946291
    ## 137  Thailand effect  -63.8462684 119.440018 286  -298.9392497  171.2467129
    ## 138    Serbia effect    4.3742708   8.333493 288   -12.0280038   20.7765454
    ## 139        UK effect   34.4131421  65.885012 287   -95.2659623  164.0922465
    ## 140  Thailand effect   -6.5335362  12.627301 288   -31.3870337   18.3199614
    ## 141    Brazil effect    1.5463335   3.041051 288    -4.4391694    7.5318363
    ## 142    Russia effect  -38.1431853  75.219859 288  -186.1935573  109.9071867
    ## 143 Lithuania effect   -3.8920744   7.843186 288   -19.3293096   11.5451608
    ## 144    Serbia effect   -2.8314650   5.809431 288   -14.2657911    8.6028611
    ## 145    Brazil effect   -7.7826947  16.475332 287   -40.2104994   24.6451099
    ## 146  Thailand effect  -33.1589175  70.859384 289  -172.6248127  106.3069777
    ## 147    Russia effect   -2.2489679   4.824337 288   -11.7443980    7.2464622
    ## 148  Thailand effect    3.8911103   8.362009 288   -12.5672904   20.3495110
    ## 149    Serbia effect    3.5013056   7.546918 288   -11.3528041   18.3554153
    ## 150        UK effect   -2.5571546   5.593637 288   -13.5667475    8.4524382
    ## 151 Lithuania effect  -66.3969269 145.839908 288  -353.4441627  220.6503089
    ## 152   Czechia effect    1.5978793   3.698421 286    -5.6816969    8.8774555
    ## 153    Canada effect  -48.7853277 113.743554 287  -272.6626840  175.0920286
    ## 154    Canada effect    4.1054445   9.662683 288   -14.9129889   23.1238780
    ## 155    Brazil effect  -21.8645866  54.433645 287  -129.0043769   85.2752038
    ## 156 Lithuania effect   -2.0231393   5.193889 288   -12.2459340    8.1996553
    ## 157    Brazil effect   18.1114931  47.019375 288   -74.4336951  110.6566813
    ## 158   Romania effect  -21.4756904  57.850215 287  -135.3401922   92.3888113
    ## 159   Czechia effect   14.0840262  39.247881 288   -63.1650339   91.3330862
    ## 160    Serbia effect  -13.0350779  42.311448 289   -96.3127428   70.2425871
    ## 161   Czechia effect  -18.2151697  59.766491 288  -135.8496802   99.4193409
    ## 162  Thailand effect  174.9684532 580.253117 289  -967.0894652 1317.0263717
    ## 163    Russia effect    1.2500276   4.678311 286    -7.9582598   10.4583151
    ## 164 Lithuania effect    1.4503004   6.037491 288   -10.4329013   13.3335021
    ## 165   Romania effect    0.7923309   3.782193 288    -6.6519141    8.2365760
    ## 166  Thailand effect  -47.0590215 234.797989 288  -509.1966825  415.0786395
    ## 167   Romania effect    0.9435751   4.740210 286    -8.3865485   10.2736987
    ## 168 Lithuania effect  -14.1010906  74.147926 286  -160.0459526  131.8437714
    ## 169     Japan effect    0.7243457   3.843940 286    -6.8416546    8.2903460
    ## 170    Russia effect   -2.9381082  16.043053 288   -34.5146094   28.6383929
    ## 171  Thailand effect  -46.3281722 329.767978 289  -695.3796278  602.7232833
    ## 172  Thailand effect    3.9754779  30.507875 289   -56.0703189   64.0212747
    ## 173 Lithuania effect  -11.9753374 122.289000 288  -252.6688491  228.7181743
    ## 174    Canada effect   -5.0214387  54.259686 289  -111.8157018  101.7728244
    ## 175    Canada effect   -0.6286057   7.438095 288   -15.2685251   14.0113137
    ## 176    Canada effect    2.1671156  32.132692 288   -61.0775799   65.4118111
    ## 177    Brazil effect   -2.7353940  43.397032 286   -88.1534786   82.6826905
    ## 178        UK effect    0.1341228   5.407459 286   -10.5093416   10.7775872
    ## 179    Canada effect    1.6711930  91.421632 286  -178.2733893  181.6157753
    ## 180    Serbia effect   -1.3774347  89.128391 287  -176.8056475  174.0507781
    ##         t.ratio      p.value Signature   Country p_adj_global
    ## 1   10.52708918 3.769709e-22     SBS12     Japan 6.785475e-20
    ## 2    7.66799231 2.679184e-13    SBS22b   Romania 1.787426e-11
    ## 3   -7.65418345 2.979043e-13       ID5     Japan 1.787426e-11
    ## 4    7.29856165 2.824450e-12    SBS22c   Romania 1.271003e-10
    ## 5    6.83890916 4.806362e-11       ID5   Czechia 1.730290e-09
    ## 6    6.39659643 6.379148e-10    SBS22a   Romania 1.913744e-08
    ## 7    5.89508921 1.049530e-08    SBS40b   Czechia 2.698790e-07
    ## 8    5.59982157 4.992811e-08    SBS22a    Serbia 1.123382e-06
    ## 9    5.38486909 1.506780e-07       ID8   Czechia 3.013559e-06
    ## 10  -5.33116252 1.974210e-07      SBS1   Romania 3.553577e-06
    ## 11   5.08160319 6.729733e-07    SBS22b    Serbia 1.101229e-05
    ## 12  -5.04611807 8.015088e-07    SBS40b     Japan 1.202263e-05
    ## 13   4.91313533 1.503254e-06      ID23   Romania 2.081429e-05
    ## 14   4.74134206 3.353292e-06      SBSB   Czechia 4.311376e-05
    ## 15   4.21982665 3.276692e-05      ID23    Serbia 3.932030e-04
    ## 16   3.88299726 1.280101e-04       ID9   Czechia 1.440114e-03
    ## 17  -3.72399885 2.358820e-04       ID8     Japan 2.449512e-03
    ## 18  -3.71400711 2.449512e-04    SBS40c     Japan 2.449512e-03
    ## 19   3.41127684 7.391459e-04    SBS40c   Romania 7.002435e-03
    ## 20   3.32116671 1.012016e-03      ID21     Japan 9.108146e-03
    ## 21  -3.06191106 2.407001e-03       ID9  Thailand 2.063144e-02
    ## 22  -3.01418452 2.805241e-03    SBS22a   Czechia 2.295197e-02
    ## 23  -2.94639482 3.476816e-03    SBS22b   Czechia 2.720986e-02
    ## 24  -2.88313706 4.234082e-03      SBS5   Romania 3.175561e-02
    ## 25  -2.84180999 4.805102e-03    SBS22a     Japan 3.459673e-02
    ## 26  -2.77045520 5.960514e-03    SBS22b     Japan 4.126510e-02
    ## 27  -2.72976164 6.728336e-03      SBS1    Serbia 4.485557e-02
    ## 28  -2.70045554 7.332942e-03      ID23   Czechia 4.714034e-02
    ## 29  -2.67176402 7.973462e-03    SBS22c    Russia 4.911228e-02
    ## 30   2.66285670 8.185379e-03       ID5    Russia 4.911228e-02
    ## 31   2.63183862 8.950171e-03      SBSC     Japan 5.099609e-02
    ## 32   2.62736060 9.065972e-03      ID11   Czechia 5.099609e-02
    ## 33  -2.56519083 1.082206e-02      SBSB   Romania 5.787488e-02
    ## 34   2.56142708 1.093192e-02    SBS22c    Serbia 5.787488e-02
    ## 35  -2.53934954 1.163009e-02    SBS22c   Czechia 5.894257e-02
    ## 36   2.53456108 1.178851e-02      SBS5        UK 5.894257e-02
    ## 37  -2.50742322 1.271005e-02    SBS22c     Japan 6.183268e-02
    ## 38  -2.47233208 1.399919e-02    SBS22a    Russia 6.631194e-02
    ## 39  -2.36942918 1.847263e-02    SBS22b    Russia 8.525830e-02
    ## 40   2.35392325 1.925296e-02      SBSB    Russia 8.663831e-02
    ## 41   2.27565894 2.360606e-02      SBSB        UK 1.036364e-01
    ## 42   2.19221494 2.916421e-02      SBS1 Lithuania 1.249895e-01
    ## 43   2.12467866 3.446798e-02       ID5   Romania 1.442846e-01
    ## 44   2.10921279 3.578922e-02       ID8    Russia 1.464104e-01
    ## 45  -2.09891646 3.669177e-02    SBS22a        UK 1.467671e-01
    ## 46  -2.04976813 4.128937e-02    SBS22c    Brazil 1.596048e-01
    ## 47  -2.04586572 4.167458e-02    SBS22b        UK 1.596048e-01
    ## 48  -2.00973291 4.539850e-02      SBSB    Serbia 1.702444e-01
    ## 49  -1.99372733 4.712152e-02      ID23    Russia 1.730995e-01
    ## 50   1.98220351 4.840675e-02      SBS5    Russia 1.742643e-01
    ## 51  -1.91963785 5.589147e-02      SBSC        UK 1.972640e-01
    ## 52   1.90948467 5.719308e-02      SBS1    Brazil 1.979761e-01
    ## 53   1.86245545 6.355708e-02      SBS5   Czechia 2.158542e-01
    ## 54  -1.80754984 7.171627e-02    SBS22c Lithuania 2.390542e-01
    ## 55   1.76131736 7.924580e-02       ID9    Brazil 2.593499e-01
    ## 56  -1.71718222 8.701720e-02    SBS22a    Brazil 2.795387e-01
    ## 57  -1.70016504 9.017947e-02      ID11     Japan 2.795387e-01
    ## 58  -1.69750248 9.068950e-02       ID3 Lithuania 2.795387e-01
    ## 59   1.68318621 9.342268e-02       ID9        UK 2.795387e-01
    ## 60  -1.67889067 9.425793e-02      SBS5    Serbia 2.795387e-01
    ## 61  -1.67644423 9.473256e-02    SBS22b    Brazil 2.795387e-01
    ## 62   1.65061275 9.990772e-02       ID9 Lithuania 2.900547e-01
    ## 63   1.62839573 1.045344e-01    SBS40a    Canada 2.977689e-01
    ## 64   1.62211267 1.058734e-01      SBSC   Romania 2.977689e-01
    ## 65  -1.60923026 1.086578e-01    SBS22c        UK 3.008984e-01
    ## 66  -1.54830028 1.226444e-01      ID23     Japan 3.245628e-01
    ## 67  -1.54723060 1.229134e-01      SBSB     Japan 3.245628e-01
    ## 68  -1.54704225 1.229512e-01       ID8  Thailand 3.245628e-01
    ## 69  -1.53767750 1.252216e-01      ID23        UK 3.245628e-01
    ## 70   1.53365948 1.262189e-01       ID3    Canada 3.245628e-01
    ## 71  -1.50649467 1.330325e-01    SBS22a Lithuania 3.347512e-01
    ## 72   1.50313301 1.339005e-01    SBS40a    Brazil 3.347512e-01
    ## 73  -1.47059808 1.424915e-01      ID11   Romania 3.429436e-01
    ## 74  -1.45706610 1.461838e-01    SBS22b Lithuania 3.429436e-01
    ## 75  -1.45605690 1.464738e-01       ID3    Brazil 3.429436e-01
    ## 76   1.45567017 1.465730e-01      ID11    Russia 3.429436e-01
    ## 77   1.45189262 1.476199e-01      SBSC    Serbia 3.429436e-01
    ## 78  -1.44832975 1.486089e-01      ID23 Lithuania 3.429436e-01
    ## 79   1.43586121 1.521309e-01    SBS40b Lithuania 3.466273e-01
    ## 80   1.39545499 1.639574e-01    SBS40b    Russia 3.665792e-01
    ## 81  -1.39211819 1.649606e-01    SBS40c        UK 3.665792e-01
    ## 82   1.38209342 1.680140e-01       ID8        UK 3.688112e-01
    ## 83  -1.33866770 1.817346e-01    SBS40c    Brazil 3.941233e-01
    ## 84  -1.32818603 1.851646e-01     SBS12   Czechia 3.967813e-01
    ## 85  -1.31738792 1.887557e-01      ID21   Czechia 3.997180e-01
    ## 86  -1.29868794 1.950935e-01       ID5    Canada 4.083353e-01
    ## 87  -1.29121700 1.976603e-01     SBS12    Brazil 4.089523e-01
    ## 88  -1.25491922 2.105326e-01       ID3    Serbia 4.214292e-01
    ## 89   1.25481861 2.105619e-01      SBS1   Czechia 4.214292e-01
    ## 90  -1.25438795 2.107146e-01     SBS12    Russia 4.214292e-01
    ## 91  -1.23298983 2.185885e-01    SBS40b  Thailand 4.323729e-01
    ## 92  -1.22199040 2.227070e-01    SBS22a    Canada 4.345503e-01
    ## 93   1.21180256 2.265778e-01    SBS22c  Thailand 4.345503e-01
    ## 94  -1.21088596 2.269318e-01    SBS40a     Japan 4.345503e-01
    ## 95  -1.19104204 2.346178e-01      SBSC    Canada 4.445389e-01
    ## 96   1.16826304 2.436668e-01       ID9    Russia 4.568753e-01
    ## 97   1.15613486 2.485875e-01       ID5 Lithuania 4.612964e-01
    ## 98   1.14621775 2.526565e-01      SBS5    Brazil 4.637622e-01
    ## 99  -1.14039511 2.550692e-01      SBSC    Russia 4.637622e-01
    ## 100 -1.12968104 2.595507e-01    SBS40a        UK 4.671912e-01
    ## 101  1.10943423 2.681684e-01      ID21    Brazil 4.755077e-01
    ## 102 -1.10645179 2.694544e-01    SBS40a   Czechia 4.755077e-01
    ## 103 -1.09495614 2.744507e-01      ID21   Romania 4.792146e-01
    ## 104 -1.08512240 2.787718e-01      ID23    Brazil 4.792146e-01
    ## 105 -1.08338221 2.795418e-01    SBS22b    Canada 4.792146e-01
    ## 106 -1.05784186 2.910109e-01     SBS12        UK 4.912079e-01
    ## 107  1.05568741 2.919958e-01    SBS40c  Thailand 4.912079e-01
    ## 108 -1.04541833 2.967094e-01       ID5    Serbia 4.945156e-01
    ## 109 -1.03648662 3.008447e-01       ID9   Romania 4.968077e-01
    ## 110 -1.02148046 3.078841e-01    SBS40a  Thailand 5.038104e-01
    ## 111  0.95522764 3.402637e-01    SBS40c    Canada 5.517789e-01
    ## 112  0.91091683 3.631010e-01      SBS5     Japan 5.835551e-01
    ## 113  0.89700567 3.704651e-01    SBS40a    Serbia 5.877739e-01
    ## 114 -0.89364726 3.722568e-01       ID9    Canada 5.877739e-01
    ## 115 -0.86717205 3.865697e-01      ID11    Canada 6.035819e-01
    ## 116 -0.86278239 3.889750e-01    SBS40c    Serbia 6.035819e-01
    ## 117 -0.84861633 3.967996e-01    SBS40a    Russia 6.104609e-01
    ## 118  0.83545704 4.041530e-01      ID11        UK 6.165045e-01
    ## 119  0.79295891 4.284545e-01    SBS40a Lithuania 6.480825e-01
    ## 120  0.77945522 4.363561e-01       ID3  Thailand 6.545342e-01
    ## 121  0.75979180 4.480005e-01      SBS1  Thailand 6.664471e-01
    ## 122 -0.75339525 4.518258e-01     SBS12 Lithuania 6.666282e-01
    ## 123 -0.72132847 4.712946e-01       ID5  Thailand 6.816921e-01
    ## 124 -0.71678216 4.740874e-01      ID23    Canada 6.816921e-01
    ## 125  0.70317213 4.825165e-01      SBS1     Japan 6.816921e-01
    ## 126  0.70014126 4.844043e-01      SBS1        UK 6.816921e-01
    ## 127 -0.69974204 4.846532e-01      SBSC  Thailand 6.816921e-01
    ## 128 -0.69957269 4.847588e-01      SBS5    Canada 6.816921e-01
    ## 129 -0.68147805 4.961167e-01    SBS40a   Romania 6.922558e-01
    ## 130 -0.65092620 5.156131e-01       ID9     Japan 7.139259e-01
    ## 131 -0.63850062 5.236535e-01     SBS12   Romania 7.191166e-01
    ## 132 -0.63137173 5.282983e-01       ID8    Brazil 7.191166e-01
    ## 133  0.62533062 5.322510e-01       ID5        UK 7.191166e-01
    ## 134 -0.62061711 5.353424e-01      ID11    Serbia 7.191166e-01
    ## 135 -0.60772929 5.438443e-01    SBS22c    Canada 7.251257e-01
    ## 136  0.54086283 5.890202e-01      SBSC Lithuania 7.739477e-01
    ## 137 -0.53454671 5.933786e-01      SBSB  Thailand 7.739477e-01
    ## 138  0.52490242 6.000548e-01       ID9    Serbia 7.739477e-01
    ## 139  0.52232126 6.018496e-01    SBS40b        UK 7.739477e-01
    ## 140 -0.51741352 6.052647e-01      ID21  Thailand 7.739477e-01
    ## 141  0.50848658 6.115014e-01      ID11    Brazil 7.739477e-01
    ## 142 -0.50708930 6.124802e-01    SBS40c    Russia 7.739477e-01
    ## 143 -0.49623638 6.201062e-01      ID21 Lithuania 7.739477e-01
    ## 144 -0.48739112 6.263520e-01       ID8    Serbia 7.739477e-01
    ## 145 -0.47238470 6.370112e-01       ID5    Brazil 7.739477e-01
    ## 146 -0.46795379 6.401703e-01     SBS12  Thailand 7.739477e-01
    ## 147 -0.46617137 6.414452e-01      ID21    Russia 7.739477e-01
    ## 148  0.46533198 6.420454e-01      ID11  Thailand 7.739477e-01
    ## 149  0.46393846 6.430423e-01      ID21    Serbia 7.739477e-01
    ## 150 -0.45715421 6.479049e-01      ID21        UK 7.739477e-01
    ## 151 -0.45527269 6.492561e-01      SBS5 Lithuania 7.739477e-01
    ## 152  0.43204370 6.660352e-01       ID3   Czechia 7.845711e-01
    ## 153 -0.42890631 6.683132e-01    SBS40b    Canada 7.845711e-01
    ## 154  0.42487623 6.712442e-01      ID21    Canada 7.845711e-01
    ## 155 -0.40167412 6.882228e-01    SBS40b    Brazil 7.992265e-01
    ## 156 -0.38952304 6.971772e-01      ID11 Lithuania 8.029828e-01
    ## 157  0.38519213 7.003794e-01      SBSC    Brazil 8.029828e-01
    ## 158 -0.37122923 7.107407e-01    SBS40b   Romania 8.097046e-01
    ## 159  0.35884806 7.199717e-01      SBSC   Czechia 8.150623e-01
    ## 160 -0.30807449 7.582476e-01     SBS12    Serbia 8.480234e-01
    ## 161 -0.30477228 7.607597e-01    SBS40c   Czechia 8.480234e-01
    ## 162  0.30153815 7.632210e-01    SBS22a  Thailand 8.480234e-01
    ## 163  0.26719636 7.895105e-01       ID3    Russia 8.718521e-01
    ## 164  0.24021576 8.103339e-01       ID8 Lithuania 8.893909e-01
    ## 165  0.20948983 8.342140e-01       ID8   Romania 9.051012e-01
    ## 166 -0.20042344 8.412909e-01      SBS5  Thailand 9.051012e-01
    ## 167  0.19905766 8.423591e-01       ID3   Romania 9.051012e-01
    ## 168 -0.19017512 8.493068e-01      SBSB Lithuania 9.051012e-01
    ## 169  0.18843838 8.506667e-01       ID3     Japan 9.051012e-01
    ## 170 -0.18313897 8.548178e-01      SBS1    Russia 9.051012e-01
    ## 171 -0.14048718 8.883729e-01    SBS22b  Thailand 9.351294e-01
    ## 172  0.13030989 8.964119e-01      ID23  Thailand 9.381055e-01
    ## 173 -0.09792653 9.220588e-01    SBS40c Lithuania 9.582719e-01
    ## 174 -0.09254456 9.263295e-01     SBS12    Canada 9.582719e-01
    ## 175 -0.08451166 9.327084e-01       ID8    Canada 9.593572e-01
    ## 176  0.06744270 9.462761e-01      SBS1    Canada 9.658832e-01
    ## 177 -0.06303182 9.497852e-01      SBSB    Brazil 9.658832e-01
    ## 178  0.02480329 9.802292e-01       ID3        UK 9.876803e-01
    ## 179  0.01828006 9.854282e-01      SBSB    Canada 9.876803e-01
    ## 180 -0.01545450 9.876803e-01    SBS40b    Serbia 9.876803e-01
