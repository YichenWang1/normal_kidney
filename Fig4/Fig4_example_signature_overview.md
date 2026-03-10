    LCM_metadata_path <- "../data/LCM_normal_kidney_metadata.csv"
    kidney_weights_path  <- "../data/SBS96_kidney_signature_weights_table.csv"
    blood_weights_path  <- "../data/SBS96_blood_signature_weights_table.csv"
    bulk_metadata_path <- "../data/bulk_normal_kidney_metadata.csv"
    blood_path <- "../data/SBS96_blood_signature_weights_table.csv"

    blood_matrix <- read.table(blood_path, header = TRUE)
    LCM_metadata <- read.csv(LCM_metadata_path, header = TRUE, stringsAsFactors = FALSE)
    bulk_metadata <- read.csv(bulk_metadata_path, header = TRUE, stringsAsFactors = FALSE)

    df <- data.frame(Patient = LCM_metadata$Patient, Sample = LCM_metadata$Sample, Structure = LCM_metadata$Structure, Country = LCM_metadata$Country, Sex = LCM_metadata$Sex, Age = LCM_metadata$Age, Burden= LCM_metadata$Normal_SBS_burden)

    df2 <- na.omit(unique(data.frame(Patient = LCM_metadata$Patient, Sample = LCM_metadata$Tumor, Structure = "Tumor", Country = LCM_metadata$Country, Sex = LCM_metadata$Sex, Age = LCM_metadata$Age, Burden= LCM_metadata$Tumor_SBS_burden)))
    df2 <- df2[!df2$Sample%in% c('PD47592a'),]

    df3 <- unique(data.frame(Patient = LCM_metadata$Patient, Sample = NA, Structure = "Blood", Country = LCM_metadata$Country, Sex = LCM_metadata$Sex, Age = LCM_metadata$Age, Burden=NA))
    df3 <- df3[df3$Patient %in% substr(colnames(blood_matrix),1,7),]
    col_idx <- match(df3$Patient, substr(colnames(blood_matrix),1,7))
    df3$Sample <- colnames(blood_matrix)[col_idx]
    df3$Burden <- colSums(blood_matrix[col_idx])

    new_rows <- data.frame(
      Patient  = "PD55516",
      Sample   = c("PD55516b_ds0001", "PD55516a", "PD55516c"),
      Structure = c("Bulk cortex", "Tumor", "Blood"),
      Country  = "Japan",
      Sex      = "Male",
      Age      = 62,
      Burden   = c(4749, 8275, 1834)
    )
    df <- rbind(df,df2,df3,new_rows)
    df$Structure[df$Structure=='Distal_tubules'] <- 'Distal tubules'
    df$Structure[df$Structure=='Proximal_tubules'] <- 'Proximal tubules'
    df$Structure <- as.factor(df$Structure)

    # Load kidney weights
    kidney_weights_raw <- read.csv(path.expand(kidney_weights_path), check.names = FALSE, stringsAsFactors = FALSE)
    kidney_weights_raw$Sample <- sub(" .*", "", kidney_weights_raw$Sample)
    kidney_weights_raw <- kidney_weights_raw[kidney_weights_raw$Sample %in% df$Sample, , drop = FALSE]
    rownames(kidney_weights_raw) <- kidney_weights_raw$Sample

    # Load blood weights
    blood_weights_raw <- read.csv(path.expand(blood_weights_path), check.names = FALSE, stringsAsFactors = FALSE)
    blood_weights_raw$Sample <- sub(" .*", "", blood_weights_raw$Sample)
    blood_weights_raw <- blood_weights_raw[blood_weights_raw$Sample %in% df$Sample, , drop = FALSE]
    rownames(blood_weights_raw) <- blood_weights_raw$Sample

    # Combine kidney and blood weights, using all columns from both
    all_cols <- unique(c(colnames(kidney_weights_raw), colnames(blood_weights_raw)))
    weights <- data.frame(matrix(NA, nrow = 0, ncol = length(all_cols)))
    colnames(weights) <- all_cols

    for (sample in df$Sample) {
      row <- data.frame(matrix(0, nrow = 1, ncol = length(all_cols)))
      colnames(row) <- all_cols
      row$Sample <- sample
      
      if (sample %in% rownames(blood_weights_raw)) {
        for (col in colnames(blood_weights_raw)) {
          row[[col]] <- blood_weights_raw[sample, col]
        }
      } else if (sample %in% rownames(kidney_weights_raw)) {
        for (col in colnames(kidney_weights_raw)) {
          row[[col]] <- kidney_weights_raw[sample, col]
        }
      }
      weights <- rbind(weights, row)
    }
    rownames(weights) <- weights$Sample

    # Ensure numeric signature columns 
    sig_cols <- setdiff(colnames(weights), "Sample")
    weights[sig_cols] <- lapply(weights[sig_cols], function(x) as.numeric(as.character(x)))

    # Zero small contributions
    weights[sig_cols] <- lapply(weights[sig_cols], function(x) { x[x < 0.05] <- 0; x })

    # Compute Unassigned and clamp negatives to zero
    weights$Unassigned <- 1 - rowSums(weights[, sig_cols, drop = FALSE], na.rm = TRUE)
    weights[sig_cols] <- lapply(weights[sig_cols], function(x) pmax(x, 0))
    weights$Unassigned <- pmax(weights$Unassigned, 0)

    df <- df %>%
      mutate(
        SBS1  = weights$SBS1 * Burden,
        SBS5  = weights$SBS5 * Burden,
        # SBS7a = weights$SBS7a * Burden,
        # SBS7d = weights$SBS7d * Burden,
        # SBS8  = weights$SBS8 * Burden,
        # SBS9  = weights$SBS9 * Burden,
        SBS12 = weights$SBS12 * Burden,
        # SBS17a = weights$SBS17a * Burden,
        # SBS17b = weights$SBS17b * Burden,
        SBS18 = weights$SBS18 * Burden,
        # SBS19 = weights$SBS19 * Burden,
        SBS22a = weights$SBS22a * Burden,
        SBS22b = weights$SBS22b * Burden,
        SBS22c = weights$SBS22c * Burden,
        SBS40a = weights$SBS40a * Burden,
        SBS40b = weights$SBS40b * Burden,
        SBS40c = weights$SBS40c * Burden,
        SBSB = weights$SBSB * Burden,
        SBSC = weights$SBSC * Burden,
        # SBSD = weights$SBSD * Burden,
        Others = (weights$SBS21 + weights$SBS44 + weights$SBS2 + weights$SBS13 + weights$SBSD) * Burden,
        SBSblood = weights$SBSblood * Burden,
        Unassigned = weights$Unassigned * Burden
      )


    head(df)

    ##   Patient          Sample        Structure        Country  Sex Age Burden
    ## 1 PD30346 PD30346b_ds0001        Glomeruli United Kingdom Male  40   1709
    ## 2 PD30346 PD30346b_ds0002        Glomeruli United Kingdom Male  40   1489
    ## 3 PD30346 PD30346b_ds0003 Proximal tubules United Kingdom Male  40   1777
    ## 4 PD30354 PD30354c_ds0001        Glomeruli United Kingdom Male  80   1722
    ## 5 PD30354 PD30354c_ds0002        Glomeruli United Kingdom Male  80   2151
    ## 6 PD30354 PD30354c_ds0003 Proximal tubules United Kingdom Male  80   4013
    ##       SBS1      SBS5 SBS12 SBS18 SBS22a   SBS22b   SBS22c   SBS40a   SBS40b
    ## 1 228.5586 1055.0134     0     0      0   0.0000   0.0000   0.0000   0.0000
    ## 2 188.1606  895.7687     0     0      0   0.0000   0.0000   0.0000   0.0000
    ## 3   0.0000    0.0000     0     0      0   0.0000 238.2843 468.6515 334.7914
    ## 4 202.8919    0.0000     0     0      0   0.0000 527.2313 991.8768   0.0000
    ## 5 253.2902    0.0000     0     0      0   0.0000   0.0000 489.7692   0.0000
    ## 6   0.0000  856.6914     0     0      0 248.3237 412.8552   0.0000 740.0137
    ##      SBS40c     SBSB     SBSC   Others SBSblood   Unassigned
    ## 1    0.0000 367.5850   0.0000   0.0000        0 5.784306e+01
    ## 2    0.0000 405.0708   0.0000   0.0000        0 0.000000e+00
    ## 3  671.2453   0.0000   0.0000   0.0000        0 6.402748e+01
    ## 4    0.0000   0.0000   0.0000   0.0000        0 4.596965e-07
    ## 5  544.5098 344.1134 379.3227 139.9947        0 7.656311e-07
    ## 6 1200.0817 458.4374   0.0000   0.0000        0 9.659698e+01

    write.csv(df, "../../data_S1/Fig4_all_samples.csv", row.names = FALSE)
    # SBS40c "#D1D69C"
    signature_names <- c(
      "SBS40a","SBS40b","SBS40c","SBS5","SBS22a",
      "SBS22b","SBS22c","SBS12","SBS18","SBS1",
      "SBSblood",
      # "APOBEC","SBSD","MSI",
      "SBSB","SBSC","Others", "Unassigned"
    )

    cols <- c(
      "#8DD3C7", "#CFECBB", "#F4F3B9", "#BD98A2",  "#1f78b4","#8AB1C9",
       "#759696", "#F5847A", "#D3B387", "#17BEBB", 
        "orange", #"#CECBD0", "purple",
      "#cab2d6","#FCCDE5", "#E5E4E6","#8A8A8A"
    )
    names(cols) <- signature_names
    sig_levels <- signature_names

    summary_AA_case <-df[df$Patient == 'PD42887'&df$Structure!='Bulk cortex',] %>% 
      group_by(Structure) %>% 
      summarize(SBS1 = mean(SBS1),
                SBS5 = mean(SBS5),
                # SBS12 = mean(SBS12),
                SBS22a = mean(SBS22a),
                SBS22b = mean(SBS22b),
                SBS22c = mean(SBS22c),
                SBS40a = mean(SBS40a),
                SBS40b = mean(SBS40b),
                SBS40c = mean(SBS40c),
                SBSB = mean(SBSB),
                SBSC = mean(SBSC),
                SBSblood = mean(SBSblood),
                Others = mean(Others+SBS18+SBS12),
                Unassigned = mean(Unassigned))

    write.csv(summary_AA_case, "../../data_S1/Fig4a_PD42887_signature_summary.csv", row.names = FALSE)

    test <- gather(summary_AA_case, E1, E2, -Structure)
    test$E1<-factor(test$E1,levels = sig_levels)
    test$Structure <-  factor(test$Structure,
      levels = factor(c('Blood','Distal tubules', 'Medulla','Glomeruli', 'Proximal tubules', 'Tumor')))

    p_sbs1 <- ggplot(test, aes(x = Structure, y = E2, fill = E1)) +
      geom_bar(stat = "identity",show.legend = TRUE) +scale_fill_manual(values=cols)+ theme_bw()+theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust=1.2),axis.ticks.x = element_blank(), legend.key.size = unit(0.2, "cm"), panel.grid=element_blank(),panel.border=element_blank(),text = element_text(size = 6))+
      labs(x = '', y = 'SBS average burden', fill = "Signatures", title = 'PD42887')

    p_sbs1

![](Fig4_example_signature_overview_files/figure-markdown_strict/plot%20PD42887,%20AA%20high-1.png)

    summary_SBS40b_case <-df[df$Patient == 'PD46868'&df$Structure!='Bulk cortex',] %>% 
      group_by(Structure) %>% 
      summarize(SBS1 = mean(SBS1),
                SBS5 = mean(SBS5),
                # SBS12 = mean(SBS12),
                # SBS22a = mean(SBS22a),
                # SBS22b = mean(SBS22b),
                # SBS22c = mean(SBS22c),
                SBS40a = mean(SBS40a),
                SBS40b = mean(SBS40b),
                SBS40c = mean(SBS40c),
                SBSB = mean(SBSB),
                SBSC = mean(SBSC),
                SBSblood = mean(SBSblood),
                Others = mean(Others+SBS18+SBS12+SBS22a+SBS22b+SBS22c),
                Unassigned = mean(Unassigned))

    write.csv(summary_SBS40b_case, "../../data_S1/Fig4b_PD46868_signature_summary.csv", row.names = FALSE)

    test <- gather(summary_SBS40b_case, E1, E2, -Structure)
    test$E1<-factor(test$E1,levels = sig_levels)
    test$Structure <-  factor(test$Structure,
      levels = factor(c('Blood','Distal tubules', 'Medulla','Glomeruli', 'Proximal tubules', 'Tumor')))

    p_sbs2 <- ggplot(test, aes(x = Structure, y = E2, fill = E1)) +
      geom_bar(stat = "identity",show.legend = TRUE) +scale_fill_manual(values=cols)+ theme_bw()+theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust=1.2),axis.ticks.x = element_blank(),legend.key.size = unit(0.2, "cm"), panel.grid=element_blank(),panel.border=element_blank(),text = element_text(size = 6))+
      labs(x = "", y = 'SBS average burden', fill = "Signatures", title = 'PD46868')

    p_sbs2

![](Fig4_example_signature_overview_files/figure-markdown_strict/plot%20PD46868,%20SBS40b%20high-1.png)

    summary_SBS12_case <-df[df$Patient == 'PD67546'&df$Structure!='Bulk cortex',] %>% 
      group_by(Structure) %>% 
      summarize(SBS1 = mean(SBS1),
                SBS5 = mean(SBS5),
                SBS12 = mean(SBS12),
                # SBS22a = mean(SBS22a),
                # SBS22b = mean(SBS22b),
                # SBS22c = mean(SBS22c),
                SBS40a = mean(SBS40a),
                SBS40b = mean(SBS40b),
                SBS40c = mean(SBS40c),
                SBSB = mean(SBSB),
                SBSC = mean(SBSC),
                SBSblood = mean(SBSblood),
                Others = mean(Others+SBS18+SBS22a+SBS22b+SBS22c),
                Unassigned = mean(Unassigned))

    write.csv(summary_SBS12_case, "../../data_S1/Fig4c_PD67546_signature_summary.csv", row.names = FALSE)

    test <- gather(summary_SBS12_case, E1, E2, -Structure)
    test$E1<-factor(test$E1,levels = sig_levels)
    test$Structure <-  factor(test$Structure,
      levels = factor(c('Distal tubules', 'Medulla','Glomeruli', 'Proximal tubules', 'Tumor')))

    p_sbs3 <- ggplot(test, aes(x = Structure, y = E2, fill = E1)) +
      geom_bar(stat = "identity",show.legend = TRUE) +scale_fill_manual(values=cols)+ theme_bw()+theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust=1.2),axis.ticks.x = element_blank(),legend.key.size = unit(0.2, "cm"), panel.grid=element_blank(),panel.border=element_blank(),text = element_text(size = 6))+
      labs(x = "", y = 'SBS average burden', fill = "Signatures", title = 'PD67546')

    p_sbs3

![](Fig4_example_signature_overview_files/figure-markdown_strict/plot%20PD67546,%20SBS12%20high-1.png)

    summary_SBS12_case <-df[df$Patient == 'PD55516',] %>% 
      group_by(Structure) %>% 
      summarize(SBS1 = mean(SBS1),
                SBS5 = mean(SBS5),
                SBS12 = mean(SBS12),
                # SBS22a = mean(SBS22a),
                # SBS22b = mean(SBS22b),
                # SBS22c = mean(SBS22c),
                SBS40a = mean(SBS40a),
                SBS40b = mean(SBS40b),
                SBS40c = mean(SBS40c),
                SBSB = mean(SBSB),
                SBSC = mean(SBSC),
                SBSblood = mean(SBSblood),
                Others = mean(Others+SBS18+SBS22a+SBS22b+SBS22c),
                Unassigned = mean(Unassigned))

    write.csv(summary_SBS12_case, "../../data_S1/Fig4d_PD55516_signature_summary.csv", row.names = FALSE)

    test <- gather(summary_SBS12_case, E1, E2, -Structure)
    test$E1<-factor(test$E1,levels = sig_levels)
    test$Structure <-  factor(test$Structure,
      levels = factor(c('Blood', 'Bulk cortex', 'Tumor')))

    p_sbs4 <- ggplot(test, aes(x = Structure, y = E2, fill = E1)) +
      geom_bar(stat = "identity",show.legend = TRUE) +scale_fill_manual(values=cols)+ theme_bw()+theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust=1.2),axis.ticks.x = element_blank(),legend.key.size = unit(0.2, "cm"), panel.grid=element_blank(),panel.border=element_blank(),text = element_text(size = 6))+
      labs(x = "", y = 'SBS average burden', fill = "Signatures", title = 'PD55516')

    p_sbs4

![](Fig4_example_signature_overview_files/figure-markdown_strict/plot%20PD55516,%20SBS12%20high-1.png)

    pdf("../../figure/Fig4/4a.pdf",width = 2.2, height = 2.5)
    p_sbs1
    dev.off()

    ## quartz_off_screen 
    ##                 2

    pdf("../../figure/Fig4/4b.pdf",width = 2.2, height = 2.5)
    p_sbs2
    dev.off()

    ## quartz_off_screen 
    ##                 2

    pdf("../../figure/Fig4/4c.pdf",width = 2.2, height = 2.5)
    p_sbs3
    dev.off()

    ## quartz_off_screen 
    ##                 2

    pdf("../../figure/Fig4/4d.pdf",width = 1.85, height = 2.3)
    p_sbs4
    dev.off()

    ## quartz_off_screen 
    ##                 2
