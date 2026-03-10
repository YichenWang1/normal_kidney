    # Load reference signatures
    ref=read.csv("../data/COSMIC_v3.4_SBS_GRCh37_with_liver_platinum_sig.txt", header=T, stringsAsFactors = F, sep='\t')

    features<-ref$Type
    rownames(ref)<-features
    ref=ref[,-1]

    mut.cols = rep(c("dodgerblue","black","red","grey70","olivedrab3","plum2"),each=16)

    # Load HDP signatures
    hdp_sigs=read.csv("../data/bulk_merged_LCM_HDP_sigs.csv",check.names = F)
    rownames(hdp_sigs)=hdp_sigs[,1]
    hdp_sigs=hdp_sigs[,-1]

    colnames(hdp_sigs) = paste0('Component',c(0:(ncol(hdp_sigs)-1)))
    ref=ref[rownames(hdp_sigs),]

![](S2_files/figure-markdown_strict/Assess%20cosine%20similarities%20for%20RCC%20reference%20signatures-1.png)
