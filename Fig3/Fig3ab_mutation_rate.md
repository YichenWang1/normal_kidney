    metadata = read.csv('../data/LCM_normal_kidney_metadata.csv') 
    metadata$Country[metadata$Country=='United Kingdom']='UK'
    metadata = metadata[metadata$Country=='UK',]

    metadata$Structure[metadata$Structure=='Distal_tubules'] ='Distal tubules'
    metadata$Structure[metadata$Structure=='Proximal_tubules'] ='Proximal tubules'
    metadata$Structure = factor(metadata$Structure, levels = c('Proximal tubules','Medulla','Glomeruli','Distal tubules'))
    metadata=metadata[!metadata$Patient %in% c('PD30369'),] #excluded by pathology review due to sever chronic changes
    table(metadata$Structure)

    ## 
    ## Proximal tubules          Medulla        Glomeruli   Distal tubules 
    ##               34                9               36                4

    lmm1 <- lme( Normal_SBS_burden  ~ Age:Structure , random = list(Patient = pdDiag(form = ~ Age:Structure-1)), weights = varIdent(form = ~1 ), data = metadata, method = "ML")
    summary(lmm1)

    ## Linear mixed-effects model fit by maximum likelihood
    ##   Data: metadata 
    ##        AIC      BIC    logLik
    ##   1239.424 1263.612 -609.7119
    ## 
    ## Random effects:
    ##  Formula: ~Age:Structure - 1 | Patient
    ##  Structure: Diagonal
    ##         Age:StructureProximal tubules Age:StructureMedulla
    ## StdDev:                      13.70238             4.847571
    ##         Age:StructureGlomeruli Age:StructureDistal tubules Residual
    ## StdDev:               5.828601                 0.001049374 261.7025
    ## 
    ## Fixed effects:  Normal_SBS_burden ~ Age:Structure 
    ##                                   Value Std.Error DF   t-value p-value
    ## (Intercept)                   -35.99451  90.29636 64 -0.398626  0.6915
    ## Age:StructureProximal tubules  60.15736   4.50885 64 13.342071  0.0000
    ## Age:StructureMedulla           43.53637   4.58002 64  9.505717  0.0000
    ## Age:StructureGlomeruli         38.93419   2.74030 64 14.208020  0.0000
    ## Age:StructureDistal tubules    23.66945   3.12799 64  7.566980  0.0000
    ##  Correlation: 
    ##                               (Intr) Ag:SPt Ag:StM Ag:StG
    ## Age:StructureProximal tubules -0.500                     
    ## Age:StructureMedulla          -0.529  0.265              
    ## Age:StructureGlomeruli        -0.717  0.359  0.380       
    ## Age:StructureDistal tubules   -0.476  0.238  0.252  0.342
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -2.91350075 -0.32079906 -0.03153375  0.49313965  3.10478400 
    ## 
    ## Number of Observations: 83
    ## Number of Groups: 15

When variance of the residuals is allowed to vary with respect to the
location. The model fitness was also improved.

    lmm1.2 <- lme( Normal_SBS_burden  ~ Age:Structure , random = list(Patient = pdDiag(form = ~ Age:Structure-1)), weights = varIdent(form = ~1 |Structure ), data = metadata, method = "ML")
    anova(lmm1, lmm1.2)

    ##        Model df      AIC      BIC    logLik   Test  L.Ratio p-value
    ## lmm1       1 10 1239.424 1263.612 -609.7119                        
    ## lmm1.2     2 13 1215.084 1246.529 -594.5420 1 vs 2 30.33968  <.0001

    lmm1.2 <- lme( Normal_SBS_burden  ~ Age:Structure , random = list(Patient = pdDiag(form = ~ Age:Structure-1)), weights = varIdent(form = ~1 |Structure ), data = metadata, method = "REML")
    summary(lmm1.2)

    ## Linear mixed-effects model fit by REML
    ##   Data: metadata 
    ##        AIC      BIC    logLik
    ##   1190.816 1221.453 -582.4081
    ## 
    ## Random effects:
    ##  Formula: ~Age:Structure - 1 | Patient
    ##  Structure: Diagonal
    ##         Age:StructureProximal tubules Age:StructureMedulla
    ## StdDev:                      13.47211             9.072564
    ##         Age:StructureGlomeruli Age:StructureDistal tubules Residual
    ## StdDev:               6.447087                    1.900139 131.3272
    ## 
    ## Variance function:
    ##  Structure: Different standard deviations per stratum
    ##  Formula: ~1 | Structure 
    ##  Parameter estimates:
    ##        Glomeruli Proximal tubules   Distal tubules          Medulla 
    ##        1.0000000        2.9816838        0.5754943        1.2408318 
    ## Fixed effects:  Normal_SBS_burden ~ Age:Structure 
    ##                                  Value Std.Error DF   t-value p-value
    ## (Intercept)                   -8.74088  53.71606 64 -0.162724  0.8712
    ## Age:StructureProximal tubules 59.91852   4.12786 64 14.515621  0.0000
    ## Age:StructureMedulla          41.98165   5.48218 64  7.657838  0.0000
    ## Age:StructureGlomeruli        38.43299   2.23615 64 17.187104  0.0000
    ## Age:StructureDistal tubules   23.52450   1.92233 64 12.237489  0.0000
    ##  Correlation: 
    ##                               (Intr) Ag:SPt Ag:StM Ag:StG
    ## Age:StructureProximal tubules -0.302                     
    ## Age:StructureMedulla          -0.351  0.106              
    ## Age:StructureGlomeruli        -0.605  0.182  0.212       
    ## Age:StructureDistal tubules   -0.636  0.192  0.223  0.385
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -1.90335387 -0.53268962 -0.07504407  0.54019933  2.15318633 
    ## 
    ## Number of Observations: 83
    ## Number of Groups: 15

    fixed.m1 <- data.frame(fixef(lmm1.2))
    intervals(lmm1.2, which = "fixed")

    ## Approximate 95% confidence intervals
    ## 
    ##  Fixed effects:
    ##                                    lower      est.    upper
    ## (Intercept)                   -116.05104 -8.740879 98.56928
    ## Age:StructureProximal tubules   51.67216 59.918522 68.16488
    ## Age:StructureMedulla            31.02973 41.981649 52.93356
    ## Age:StructureGlomeruli          33.96576 38.432990 42.90022
    ## Age:StructureDistal tubules     19.68420 23.524500 27.36480

    reg_summary_SBS <- data.frame(
      Term      = rownames(summary(lmm1.2)$tTable),
      Estimate  = summary(lmm1.2)$tTable[, "Value"],
      Std.Error = summary(lmm1.2)$tTable[, "Std.Error"],
      CI_lower  = intervals(lmm1.2, which = "fixed")$fixed[, "lower"],
      CI_upper  = intervals(lmm1.2, which = "fixed")$fixed[, "upper"],
      p_value   = summary(lmm1.2)$tTable[, "p-value"]
    )
    write.csv(reg_summary_SBS, "../../data_S1/Fig3a_LCM_SBS_mutation_rate_regression.csv", row.names = FALSE)

    p_sbs=ggplot(data = metadata, mapping = aes(x = Age, y = Normal_SBS_burden))+geom_point(aes(colour = Structure,fill=Structure),alpha = 0.8) +theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.25,colour="black"),axis.ticks = element_line(linewidth = 0.25))+
      
            geom_abline(intercept = fixed.m1[1,], slope =  fixed.m1['Age:StructureDistal tubules',],colour='#C77CFF')+
            geom_ribbon(aes(ymin = fixed.m1[1,]+Age*intervals(lmm1.2, which = "fixed")[["fixed"]]['Age:StructureDistal tubules','lower'], ymax = fixed.m1[1,]+Age*intervals(lmm1.2, which = "fixed")[["fixed"]]['Age:StructureDistal tubules','upper']),fill='#C77CFF', alpha = 0.1)+
      
            geom_abline(intercept = fixed.m1[1,], slope =  fixed.m1['Age:StructureGlomeruli',],colour='#00BFC4')+
            geom_ribbon(aes(ymin = fixed.m1[1,]+
                              Age*intervals(lmm1.2, which = "fixed")[["fixed"]]['Age:StructureGlomeruli','lower'], ymax = fixed.m1[1,]+Age*intervals(lmm1.2, which = "fixed")[["fixed"]]['Age:StructureGlomeruli','upper']), fill='#00BFC4',alpha = 0.1)+
      
            geom_abline(intercept = fixed.m1[1,], slope =  fixed.m1['Age:StructureMedulla',],colour='#7CAE00')+
            geom_ribbon(aes(ymin = fixed.m1[1,]+Age*intervals(lmm1.2, which = "fixed")[["fixed"]]['Age:StructureMedulla','lower'], ymax = fixed.m1[1,]+Age*intervals(lmm1.2, which = "fixed")[["fixed"]]['Age:StructureMedulla','upper']), fill='#7CAE00',alpha = 0.1)+
      
            geom_abline(intercept = fixed.m1[1,], slope =  fixed.m1['Age:StructureProximal tubules',],colour='#F8766D')+
            geom_ribbon(aes(ymin = fixed.m1[1,]+Age*intervals(lmm1.2, which = "fixed")[["fixed"]]['Age:StructureProximal tubules','lower'], ymax = fixed.m1[1,]+Age*intervals(lmm1.2, which = "fixed")[["fixed"]]['Age:StructureProximal tubules','upper']), fill='#F8766D',alpha = 0.1)+
            labs(y='Substitutions / genome',x="Age (yrs)",  fill = "Structure", color = "Structure")+
            scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05), add = c(0, 0)))+
            theme(legend.key.size = unit(0.25, "cm"),
                  title=element_text(size=7),
                  axis.text.y=element_text(size=7,color="black"),
                  axis.text.x=element_text(size=7,color="black"),
                  legend.text=element_text(size=7))

    ## Warning: The `size` argument of `element_line()` is deprecated as of ggplot2 3.4.0.
    ## ℹ Please use the `linewidth` argument instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    p_sbs

![](Fig3ab_mutation_rate_files/figure-markdown_strict/plotting%20SBS%20mutation%20rate-1.png)

    lmm2 <- lme( Normal_ID_burden  ~ Age:Structure , random = list(Patient = pdDiag(form = ~ Age:Structure-1)), weights = varIdent(form = ~1 ), data = metadata, method = "ML")
    summary(lmm2)

    ## Linear mixed-effects model fit by maximum likelihood
    ##   Data: metadata 
    ##        AIC      BIC    logLik
    ##   918.5243 942.7127 -449.2621
    ## 
    ## Random effects:
    ##  Formula: ~Age:Structure - 1 | Patient
    ##  Structure: Diagonal
    ##         Age:StructureProximal tubules Age:StructureMedulla
    ## StdDev:                      3.259111             1.668108
    ##         Age:StructureGlomeruli Age:StructureDistal tubules Residual
    ## StdDev:           0.0003215077                2.070156e-06 38.24657
    ## 
    ## Fixed effects:  Normal_ID_burden ~ Age:Structure 
    ##                                  Value Std.Error DF  t-value p-value
    ## (Intercept)                   2.460292 11.981147 64 0.205347  0.8380
    ## Age:StructureProximal tubules 7.699966  0.951258 64 8.094511  0.0000
    ## Age:StructureMedulla          3.944491  1.094368 64 3.604355  0.0006
    ## Age:StructureGlomeruli        1.950147  0.245990 64 7.927752  0.0000
    ## Age:StructureDistal tubules   1.966665  0.447935 64 4.390515  0.0000
    ##  Correlation: 
    ##                               (Intr) Ag:SPt Ag:StM Ag:StG
    ## Age:StructureProximal tubules -0.337                     
    ## Age:StructureMedulla          -0.365  0.123              
    ## Age:StructureGlomeruli        -0.852  0.287  0.311       
    ## Age:StructureDistal tubules   -0.441  0.149  0.161  0.376
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -2.73569155 -0.33152434 -0.03304409  0.49900281  2.91187391 
    ## 
    ## Number of Observations: 83
    ## Number of Groups: 15

    fixed.m2 <- data.frame(fixef(lmm2))
    intervals(lmm2, which = "fixed")

    ## Approximate 95% confidence intervals
    ## 
    ##  Fixed effects:
    ##                                    lower     est.     upper
    ## (Intercept)                   -20.742666 2.460292 25.663251
    ## Age:StructureProximal tubules   5.857739 7.699966  9.542192
    ## Age:StructureMedulla            1.825113 3.944491  6.063868
    ## Age:StructureGlomeruli          1.473757 1.950147  2.426537
    ## Age:StructureDistal tubules     1.099184 1.966665  2.834145

    lmm2.2 <- lme( Normal_ID_burden  ~ Age:Structure , random = list(Patient = pdDiag(form = ~ Age:Structure-1)), weights = varIdent(form = ~1 |Structure ), data = metadata, method = "ML")
    anova(lmm2, lmm2.2)

    ##        Model df      AIC      BIC    logLik   Test  L.Ratio p-value
    ## lmm2       1 10 918.5243 942.7127 -449.2621                        
    ## lmm2.2     2 13 874.9723 906.4173 -424.4862 1 vs 2 49.55195  <.0001

    lmm2.2 <- lme( Normal_ID_burden  ~ Age:Structure , random = list(Patient = pdDiag(form = ~ Age:Structure-1)), weights = varIdent(form = ~1 |Structure ), data = metadata, method = "REML")
    summary(lmm2.2)

    ## Linear mixed-effects model fit by REML
    ##   Data: metadata 
    ##        AIC     BIC    logLik
    ##   872.8308 903.468 -423.4154
    ## 
    ## Random effects:
    ##  Formula: ~Age:Structure - 1 | Patient
    ##  Structure: Diagonal
    ##         Age:StructureProximal tubules Age:StructureMedulla
    ## StdDev:                       3.30868             1.949913
    ##         Age:StructureGlomeruli Age:StructureDistal tubules Residual
    ## StdDev:              0.3199542                5.332935e-05 14.34734
    ## 
    ## Variance function:
    ##  Structure: Different standard deviations per stratum
    ##  Formula: ~1 | Structure 
    ##  Parameter estimates:
    ##        Glomeruli Proximal tubules   Distal tubules          Medulla 
    ##         1.000000         4.032902         0.547227         3.248083 
    ## Fixed effects:  Normal_ID_burden ~ Age:Structure 
    ##                                  Value Std.Error DF   t-value p-value
    ## (Intercept)                   2.889993  4.818021 64  0.599830  0.5507
    ## Age:StructureProximal tubules 7.760436  0.914434 64  8.486603  0.0000
    ## Age:StructureMedulla          3.944413  1.175511 64  3.355488  0.0013
    ## Age:StructureGlomeruli        1.998663  0.146019 64 13.687646  0.0000
    ## Age:StructureDistal tubules   1.959574  0.112780 64 17.375172  0.0000
    ##  Correlation: 
    ##                               (Intr) Ag:SPt Ag:StM Ag:StG
    ## Age:StructureProximal tubules -0.133                     
    ## Age:StructureMedulla          -0.135  0.018              
    ## Age:StructureGlomeruli        -0.718  0.096  0.097       
    ## Age:StructureDistal tubules   -0.705  0.094  0.095  0.506
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -1.81731878 -0.46366171 -0.06479164  0.50860815  1.99224667 
    ## 
    ## Number of Observations: 83
    ## Number of Groups: 15

    fixed.m2 <- data.frame(fixef(lmm2.2))
    intervals(lmm2.2, which = "fixed")

    ## Approximate 95% confidence intervals
    ## 
    ##  Fixed effects:
    ##                                   lower     est.     upper
    ## (Intercept)                   -6.735111 2.889993 12.515097
    ## Age:StructureProximal tubules  5.933645 7.760436  9.587228
    ## Age:StructureMedulla           1.596060 3.944413  6.292766
    ## Age:StructureGlomeruli         1.706955 1.998663  2.290370
    ## Age:StructureDistal tubules    1.734269 1.959574  2.184878

    reg_summary_ID <- data.frame(
      Term      = rownames(summary(lmm2.2)$tTable),
      Estimate  = summary(lmm2.2)$tTable[, "Value"],
      Std.Error = summary(lmm2.2)$tTable[, "Std.Error"],
      CI_lower  = intervals(lmm2.2, which = "fixed")$fixed[, "lower"],
      CI_upper  = intervals(lmm2.2, which = "fixed")$fixed[, "upper"],
      p_value   = summary(lmm2.2)$tTable[, "p-value"]
    )
    write.csv(reg_summary_ID, "../../data_S1/Fig3b_LCM_ID_mutation_rate_regression.csv", row.names = FALSE)

    p_id=ggplot(data = metadata, mapping = aes(x = Age, y = Normal_ID_burden))+geom_point(aes(colour = Structure,fill=Structure),alpha = 0.8) +theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.25,colour="black"),axis.ticks = element_line(linewidth = 0.25))+
      
           geom_abline(intercept = fixed.m2[1,], slope =  fixed.m2['Age:StructureDistal tubules',],colour='#C77CFF')+
            geom_ribbon(aes(ymin = fixed.m2[1,]+Age*intervals(lmm2.2, which = "fixed")[["fixed"]]['Age:StructureDistal tubules','lower'], ymax = fixed.m2[1,]+Age*intervals(lmm2.2, which = "fixed")[["fixed"]]['Age:StructureDistal tubules','upper']),fill='#C77CFF', alpha = 0.1)+
      
            geom_abline(intercept = fixed.m2[1,], slope =  fixed.m2['Age:StructureGlomeruli',],colour='#00BFC4')+
            geom_ribbon(aes(ymin = fixed.m2[1,]+
                              Age*intervals(lmm2.2, which = "fixed")[["fixed"]]['Age:StructureGlomeruli','lower'], ymax = fixed.m2[1,]+Age*intervals(lmm2.2, which = "fixed")[["fixed"]]['Age:StructureGlomeruli','upper']), fill='#00BFC4',alpha = 0.1)+
      
            geom_abline(intercept = fixed.m2[1,], slope =  fixed.m2['Age:StructureMedulla',],colour='#7CAE00')+
            geom_ribbon(aes(ymin = fixed.m2[1,]+Age*intervals(lmm2.2, which = "fixed")[["fixed"]]['Age:StructureMedulla','lower'], ymax = fixed.m2[1,]+Age*intervals(lmm2.2, which = "fixed")[["fixed"]]['Age:StructureMedulla','upper']), fill='#7CAE00',alpha = 0.1)+
      
            geom_abline(intercept = fixed.m2[1,], slope =  fixed.m2['Age:StructureProximal tubules',],colour='#F8766D')+
            geom_ribbon(aes(ymin = fixed.m2[1,]+Age*intervals(lmm2.2, which = "fixed")[["fixed"]]['Age:StructureProximal tubules','lower'], ymax = fixed.m2[1,]+Age*intervals(lmm2.2, which = "fixed")[["fixed"]]['Age:StructureProximal tubules','upper']), fill='#F8766D',alpha = 0.1)+
            labs(y='Indels / genome',x="Age (yrs)",  fill = "Structure", color = "Structure")+
            scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05), add = c(0, 0)))+
            theme(legend.key.size = unit(0.25, "cm"),
                  title=element_text(size=7),
                  axis.text.y=element_text(size=7,color="black"),
                  axis.text.x=element_text(size=7,color="black"),
                  legend.text=element_text(size=7))

    p_id

![](Fig3ab_mutation_rate_files/figure-markdown_strict/plotting%20indel%20mutation%20rate-1.png)

    pdf("../../figure/mutation rate/mutation_rate.pdf",width = 6.04, height = 2)

    ggarrange(p_sbs, p_id, ncol = 2,common.legend = FALSE, legend = "right")#labels = c('a','b'),


    dev.off()

    ## quartz_off_screen 
    ##                 2
