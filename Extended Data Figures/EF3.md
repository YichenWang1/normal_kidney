    # Load reference signatures
    ref=read.csv("../data/COSMIC_v3.4_SBS_GRCh37_with_liver_platinum_sig.txt", header=T, stringsAsFactors = F, sep='\t')
    features<-ref$Type
    rownames(ref)<-features
    ref=ref[,-1]

    mut.cols = rep(c("dodgerblue","black","red","grey70","olivedrab3","plum2"),each=16)

    # Load Sigprofiler signatures
    sigprofiler_sigs=read.table("../data/S3_Sigprofiler_Signatures.txt",check.names = F, header=T)
    rownames(sigprofiler_sigs)=sigprofiler_sigs[,1]
    sigprofiler_sigs=sigprofiler_sigs[,-1]

![](S3_files/figure-markdown_strict/Assess%20cosine%20similarities%20for%20RCC%20reference%20signatures-1.png)
