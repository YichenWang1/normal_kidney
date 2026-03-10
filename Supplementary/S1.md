    metadata_path <- "../data/bulk_normal_kidney_metadata.csv"
    metadata <- read.csv(metadata_path, header = TRUE, stringsAsFactors = FALSE)
    metadata <- metadata %>%
      mutate(
        Country = recode(Country, "United Kingdom" = "UK", "Czech Republic" = "Czechia"),
        Country = factor(Country)
      )
    metadata <- metadata[metadata$Normal_Kidney!='PD47592c_ds0003',]
    metadata <- metadata[,c("Normal_Kidney","Country","Normal_SBS_burden","Normal_ID_burden")]
    write.csv(metadata, "../../data_S1/S1_bulk_normal_mutation_burden.csv", row.names = FALSE)

    pdf("../../figure/Extended Data Fig/EF1/1a_sbs_burden.pdf",width = 3, height = 3)

    kw <- kruskal.test(Normal_SBS_burden ~ Country, data = metadata)

    p_sbs <- ggplot(metadata, aes(x = Country, y = Normal_SBS_burden, fill = Country)) +
      geom_boxplot(
        outlier.shape = NA,
        width = 0.7,
        alpha = 0.9
      ) +
    geom_jitter(
      width = 0.15,
      size = 1.2,
      alpha = 0.2,
      shape = 16,
      # fill = "grey30",
      colour = "grey30"
    )+
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 7),
        plot.title = element_text(size = 7),
        legend.position = "none"
      ) +
      labs(
        x = NULL,
        y = "SBS mutational burden",
        title = paste0(
          "Kruskal-Wallis, p = ",
          formatC(kw$p.value, format = "e", digits = 2)
        )
      )

    p_sbs

    dev.off()

    ## quartz_off_screen 
    ##                 2

    p_sbs

![](S1_files/figure-markdown_strict/plot%20SBS-1.png)

    pdf("../../figure/Extended Data Fig/EF1/1b_indel_burden.pdf",width = 3, height =3)

    kw <- kruskal.test(Normal_ID_burden ~ Country, data = metadata)


    p_id <- ggplot(metadata, aes(x = Country, y = Normal_ID_burden, fill = Country)) +
      geom_boxplot(
        outlier.shape = NA,
        width = 0.7,
        alpha = 0.9
      ) +
    geom_jitter(
      width = 0.15,
      size = 1.2,
      alpha = 0.2,
      shape = 16,
      # fill = "grey30",
      colour = "grey30"
    )+
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 7),
        plot.title = element_text(size = 7),
        legend.position = "none"
      ) +
      labs(
        x = NULL,
        y = "Indel mutational burden",
        title = paste0(
          "Kruskal-Wallis, p = ",
          formatC(kw$p.value, format = "e", digits = 2)
        )
      )

    p_id

    dev.off()

    ## quartz_off_screen 
    ##                 2

    p_id

![](S1_files/figure-markdown_strict/plot%20indel-1.png)
