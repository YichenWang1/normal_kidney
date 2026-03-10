    metadata_path <- "../data/LCM_normal_kidney_metadata.csv"
    weights_path  <- "../data/SBS96_kidney_signature_weights_table.csv"

    LCM_metadata <- read.csv(metadata_path, header = TRUE, stringsAsFactors = FALSE)

    df = data.frame(Patient = LCM_metadata$Patient, Sample = LCM_metadata$Sample, Structure = LCM_metadata$Structure, Country = LCM_metadata$Country, Sex = LCM_metadata$Sex, Age = LCM_metadata$Age, Burden= LCM_metadata$Normal_SBS_burden)
    df2 = na.omit(unique(data.frame(Patient = LCM_metadata$Patient, Sample = LCM_metadata$Tumor, Structure = "Tumor", Country = LCM_metadata$Country, Sex = LCM_metadata$Sex, Age = LCM_metadata$Age, Burden= LCM_metadata$Tumor_SBS_burden)))
    df2 <- df2[!df2$Sample%in% c('PD47592a'),]
    df = rbind(df,df2)
    df$Structure[df$Structure=='Distal_tubules'] <- 'Distal tubules'
    df$Structure[df$Structure=='Proximal_tubules'] <- 'Proximal tubules'
    df$Structure = as.factor(df$Structure)

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

    # weights[,2:ncol(weights)] = weights[,2:ncol(weights)]/apply(weights[,2:ncol(weights)],1,sum)
    df<- df %>%
      mutate(
        SBS1 = weights$SBS1 * Burden,
        SBS5 = weights$SBS5 * Burden,
        SBS12 = weights$SBS12 * Burden,
        APOBEC = (weights$SBS2+weights$SBS13) * Burden,
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
        MSI = (weights$SBS21+weights$SBS44) * Burden,
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
    ##       SBS1      SBS5 SBS12 APOBEC SBS18 SBS22a   SBS22b   SBS22c   SBS40a
    ## 1 228.5586 1055.0134     0      0     0      0   0.0000   0.0000   0.0000
    ## 2 188.1606  895.7687     0      0     0      0   0.0000   0.0000   0.0000
    ## 3   0.0000    0.0000     0      0     0      0   0.0000 238.2843 468.6515
    ## 4 202.8919    0.0000     0      0     0      0   0.0000 527.2313 991.8768
    ## 5 253.2902    0.0000     0      0     0      0   0.0000   0.0000 489.7692
    ## 6   0.0000  856.6914     0      0     0      0 248.3237 412.8552   0.0000
    ##     SBS40b    SBS40c     SBSB     SBSC SBSD      MSI   Unassigned
    ## 1   0.0000    0.0000 367.5850   0.0000    0   0.0000 5.784306e+01
    ## 2   0.0000    0.0000 405.0708   0.0000    0   0.0000 0.000000e+00
    ## 3 334.7914  671.2453   0.0000   0.0000    0   0.0000 6.402748e+01
    ## 4   0.0000    0.0000   0.0000   0.0000    0   0.0000 4.596965e-07
    ## 5   0.0000  544.5098 344.1134 379.3227    0 139.9947 7.656311e-07
    ## 6 740.0137 1200.0817 458.4374   0.0000    0   0.0000 9.659698e+01

    # write.csv(df, "../../data_S1/Fig3d_LCM_SBS_all_samples.csv", row.names = FALSE)
    # SBS40c "#D1D69C"
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

    SBS22_postive_cancer = unique(df$Patient[(df$SBS22a+df$SBS22b+df$SBS22c)>0.05*df$Burden & df$Structure=='Tumor'])
    SBS22_postive_normal = unique(df$Patient[(df$SBS22a+df$SBS22b+df$SBS22c)>0.05*df$Burden & df$Structure!='Tumor'])
    SBS22_postive = SBS22_postive_normal[SBS22_postive_normal%in% SBS22_postive_cancer]

    summary_lcm_AApos <-df[df$Patient %in% SBS22_postive,] %>% 
      group_by(Structure) %>% 
      summarize(SBS1 = mean(SBS1),
                SBS5 = mean(SBS5),
                SBS12 = mean(SBS12),
                # APOBEC = mean(APOBEC),
                # SBS18 = mean(SBS18),
                SBS22a = mean(SBS22a),
                SBS22b = mean(SBS22b),
                SBS22c = mean(SBS22c),
                SBS40a = mean(SBS40a),
                SBS40b = mean(SBS40b),
                SBS40c = mean(SBS40c),
                SBSB = mean(SBSB),
                SBSC = mean(SBSC),
                # SBSD = mean(SBSD),
                # MSI = mean(MSI),
                Others = mean(APOBEC+SBS18+SBSD+MSI),
                Unassigned = mean(Unassigned))

    write.csv(summary_lcm_AApos, "../../data_S1/Fig3d_LCM_SBS_AA_exposed.csv", row.names = FALSE)

    test <- gather(summary_lcm_AApos, E1, E2, -Structure)
    test$E1<-factor(test$E1,levels = sig_levels)
    test$Structure <-  factor(test$Structure,
      levels = factor(c('Distal tubules', 'Medulla','Glomeruli', 'Proximal tubules', 'Tumor')))

    p_sbs1 <- ggplot(test, aes(x = Structure, y = E2, fill = E1)) +
      geom_bar(stat = "identity",show.legend = TRUE) +scale_fill_manual(values=cols)+ theme_bw()+
      theme(axis.ticks.x = element_blank(),
        panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.title.x = element_text(size = 8, margin = margin(t = -4)),
        axis.title.y = element_text(size = 8),
        axis.text.x  = element_text(size = 8, angle = 60, hjust = 1, vjust = 1.1),
        axis.text.y  = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text  = element_text(size = 8),
        plot.title   = element_text(size = 8))+
      labs(x = 'Structure', y = 'SBS average burden', fill = "Signatures", title = 'Aristolochic acid exposed')

    p_sbs1

![](Fig3d_LCM_SBS_overview_files/figure-markdown_strict/plot%20LCM%20data,%20AA%20exposed%20cases-1.png)

    SBS12_postive_cancer = unique(df$Patient[(df$SBS12)>0.05*df$Burden & df$Structure=='Tumor'])
    SBS12_postive_normal = unique(df$Patient[(df$SBS12)>0.05*df$Burden & df$Structure!='Tumor'])
    SBS12_postive = SBS12_postive_normal[SBS12_postive_normal%in% SBS12_postive_cancer]
    SBS12_postive = SBS12_postive[SBS12_postive%in%df$Patient[df$Country=='Japan']]
    summary_lcm <-df[df$Patient %in% SBS12_postive,] %>% 
      group_by(Structure) %>% 
      summarize(SBS1 = mean(SBS1),
                SBS5 = mean(SBS5),
                SBS12 = mean(SBS12),
                # APOBEC = mean(APOBEC),
                # SBS18 = mean(SBS18),
                SBS22a = mean(SBS22a),
                SBS22b = mean(SBS22b),
                SBS22c = mean(SBS22c),
                SBS40a = mean(SBS40a),
                SBS40b = mean(SBS40b),
                SBS40c = mean(SBS40c),
                SBSB = mean(SBSB),
                SBSC = mean(SBSC),
                # SBSD = mean(SBSD),
                # MSI = mean(MSI),
                Others = mean(APOBEC+SBS18+SBSD+MSI),
                Unassigned = mean(Unassigned))

    write.csv(summary_lcm, "../../data_S1/Fig3d_LCM_SBS_SBS12_positive.csv", row.names = FALSE)

    test <- gather(summary_lcm, E1, E2, -Structure)
    test$E1<-factor(test$E1,levels = sig_levels)
    test$Structure <-  factor(test$Structure,
      levels = factor(c('Distal tubules', 'Medulla','Glomeruli', 'Proximal tubules', 'Tumor')))

    p_sbs2 <- ggplot(test, aes(x = Structure, y = E2, fill = E1)) +
      geom_bar(stat = "identity",show.legend = TRUE) +scale_fill_manual(values=cols)+ theme_bw()+
      theme(axis.ticks.x = element_blank(),
        panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.title.x = element_text(size = 8, margin = margin(t = -4)),
        axis.title.y = element_text(size = 8),
        axis.text.x  = element_text(size = 8, angle = 60, hjust = 1, vjust = 1.1),
        axis.text.y  = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text  = element_text(size = 8),
        plot.title   = element_text(size = 8))+
      labs(x = "Structure", y = 'SBS average burden', fill = "Signatures", title = 'SBS12 positive')

    p_sbs2

![](Fig3d_LCM_SBS_overview_files/figure-markdown_strict/plot%20LCM%20data,%20SBS12%20positive%20cases-1.png)

    summary_lcm <-df[!df$Patient %in% c(SBS12_postive,SBS22_postive,df$Patient[df$Country == 'UK']),] %>% 
      group_by(Structure) %>% 
      summarize(SBS1 = mean(SBS1),
                SBS5 = mean(SBS5),
                SBS12 = mean(SBS12),
                # APOBEC = mean(APOBEC),
                # SBS18 = mean(SBS18),
                SBS22a = mean(SBS22a),
                SBS22b = mean(SBS22b),
                SBS22c = mean(SBS22c),
                SBS40a = mean(SBS40a),
                SBS40b = mean(SBS40b),
                SBS40c = mean(SBS40c),
                SBSB = mean(SBSB),
                SBSC = mean(SBSC),
                # SBSD = mean(SBSD),
                # MSI = mean(MSI),
                Others = mean(APOBEC+SBS18+SBSD+MSI),
                Unassigned = mean(Unassigned))

    write.csv(summary_lcm, "../../data_S1/Fig3d_LCM_SBS_remaining.csv", row.names = FALSE)

    test <- gather(summary_lcm, E1, E2, -Structure)
    test$E1<-factor(test$E1,levels = sig_levels)
    test$Structure <-  factor(test$Structure,
      levels = factor(c('Distal tubules', 'Medulla','Glomeruli', 'Proximal tubules', 'Tumor')))

    p_sbs3 <- ggplot(test, aes(x = Structure, y = E2, fill = E1)) +
      geom_bar(stat = "identity",show.legend = TRUE) +scale_fill_manual(values=cols)+ theme_bw()+
      theme(axis.ticks.x = element_blank(),
        panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.title.x = element_text(size = 8, margin = margin(t = -4)),
        axis.title.y = element_text(size = 8),
        axis.text.x  = element_text(size = 8, angle = 60, hjust = 1, vjust = 1.1),
        axis.text.y  = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text  = element_text(size = 8),
        plot.title   = element_text(size = 8))+
      labs(x = "Structure", y = 'SBS average burden', fill = "Signatures", title = 'Remaining cases')

    p_sbs3

![](Fig3d_LCM_SBS_overview_files/figure-markdown_strict/plot%20LCM%20data,%20remaining%20cases-1.png)

    library(ggpubr)
    pdf("../../figure/LCM/sigs_compare.pdf",width = 7, height = 3.5)

    ggarrange(p_sbs1, p_sbs2,p_sbs3, ncol = 3,common.legend = TRUE, legend = "right")

    dev.off()

    ## quartz_off_screen 
    ##                 2
