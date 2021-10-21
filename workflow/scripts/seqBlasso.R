train_blasso <- function(fasta_file, wl_file, pos_file = NULL, amino_feature_file, base_wl_file = NULL, n_sample = 1000, dim_prescreening = 100)
{
    debug_mode = F
    ## debug_mode = T
    set.seed(1);
    seed_blr = 1000
    n_sample_blr = n_sample + 1
    ## n_sample_blr = 10001
    ## n_sample_blr = 11
    
    use_aa_feat = T # T: aa feature, F: ohv
    standardize = F
    # dim_prescreening = 100
    nfold <- 10

    if (!is.null(pos_file)) {
        position <- as.matrix(read.table(pos_file))
    } else {
        position <- NULL
    }

    cat("Loading FASTA file ...\n")
    fasta <- load_fasta_data(fasta_file, position)    
    sample_names <- unlist(getAnnot(fasta$fasta))
    cat("Transforming dataset ...\n")    
    X <- fasta$seq.bin
    wavelength <- as.matrix(read.table(wl_file))
    Y = wavelength
    n_all <- length(Y)

    if(n_all != nrow(X)) {
        stop("Length of wavelength and FASTA file is not consistent")
    }
    if (!is.null(base_wl_file)) {
        base_wavelen <- as.matrix(read.table(base_wl_file))
    } else {
        base_wavelen <- max(Y,na.rm=T)
    }    

    symbols_to_rm <- c("-", "*", "x")
    remove_feat_idx <- which(is.element(colnames(X), symbols_to_rm))

    te_idx <- which(is.na(Y))    
    tr_idx <- setdiff(1:n_all, te_idx)
    n_tr <- length(tr_idx)
    n_te <- length(te_idx)
    # aseq_full <- fasta$fasta_full_seq
    aseq_sub <- fasta$fasta.seq
    aa_names <- fasta$amino.names[-na.omit(match(symbols_to_rm,fasta$amino.names))]
    seq_len <- fasta$seq.len

    ## amino-acid feature
    aa_feat <- read.table(amino_feature_file,row.names=1,header=T,sep=",")
    n_aa_feat <- ncol(aa_feat)
    aa_feat <- aa_feat / outer(rep(1,nrow(aa_feat)),apply(abs(aa_feat),2,max))

    ## Create design matrix with amino features
    if(use_aa_feat) {
        X = do.call("cbind", lapply(1:ncol(aa_feat), function(i) {
            matrix(aa_feat[aseq_sub,i],n_all,seq_len)
        }))
        for(i in 1:nrow(X)) {
            X[i,] = c(t(matrix(X[i,],ncol(aseq_sub),ncol(aa_feat))))
        }
    }

    ## Find indices having non-zero variances
    std_X = apply(X[tr_idx,],2,std)
    x_idx = which(std_X != 0)

    ## Use high correlation dimsion only
    if(length(x_idx) > dim_prescreening) {
        corXY <- cor(X[tr_idx,x_idx], Y[tr_idx])
        cor_idx <- order(abs(corXY),decreasing=T)
        dim_tmp <- min(length(cor_idx),dim_prescreening)
        x_idx <- x_idx[cor_idx[1:dim_tmp]]
    }
    
    cat(" Num. training samples =", n_tr, "\n")
    cat(" Num. test samples =", n_te, "\n")
    cat(" feature type =", ifelse(use_aa_feat,"Amino-acid descriptor","One-hot-feature"), "\n")
    cat(" Length of sequence =", seq_len, "\n")
    cat(" Num. features =", ncol(X), "\n")
    cat(" Num. non-constant features =", length(x_idx), "\n")
    cat(" Num. posterior-samples in BLASSO =", n_sample_blr-1, "\n")
    
    ## -------------------------
    cat("Estimating BLASSO ...\n")
    if(debug_mode) {
        stop()
    } else {
        set.seed(seed_blr)
        ## Traning BLASSO
        blr <- blasso(X[tr_idx,x_idx],Y[tr_idx],normalize=F,T=n_sample_blr)
        beta <- matrix(0,ncol(X)+1,n_sample_blr-1)
        beta[c(1,1+x_idx),] <- t(cbind(blr$mu[2:n_sample_blr], blr$beta[2:n_sample_blr,]))
    }

    beta_mean <- rowMeans(beta)
    
    pred_wavelen_sample <- cbind(1,X) %*% beta
    pred_wavelen_mean <- cbind(1,X) %*% beta_mean

    base_wl_mat <- matrix(rep(base_wavelen,n_sample),length(base_wavelen),n_sample)
    
    e_gain <- rowMeans(pmax(pred_wavelen_sample - base_wl_mat, 0))
    e_gain[tr_idx] = NA
    
    list(pred_wavelen_mean = pred_wavelen_mean,
         pred_wavelen_sample = pred_wavelen_sample,
         e_gain = e_gain,
         base_wavelen = base_wavelen,
         fasta_info = fasta,
         sample_names = sample_names,
         aa_names = aa_names,
         position = position,
         tr_idx = tr_idx,
         te_idx = te_idx,
         X = X,
         Y = Y,
         beta_sample = beta,
         beta_mean = beta_mean,
         aa_feat = aa_feat, 
         wavelength = wavelength,
         remove_feat_idx = remove_feat_idx)    
}


load_fasta_data <- function(fasta_file, position = NULL)
{
    fasta_raw_data <- read.fasta(fasta_file)
    fasta.seq <- t(sapply(fasta_raw_data, function(s) {
        s
    }))
    position_in_aligned_seq <- NULL
    fasta_full_seq <- fasta.seq
    if (!is.null(position)) {
        fasta.seq <- fasta.seq[, position]
    }
    binary.data.concat <- model.matrix(~. + 0, data = data.frame(factor(fasta.seq)))
    colnames(binary.data.concat) <- sub("factor.fasta.seq.",
        "", colnames(binary.data.concat))
    amino.names <- colnames(binary.data.concat)
    poswise.amino.list <- (lapply(1:ncol(fasta.seq), function(i) {
        (binary.data.concat[((i - 1) * nrow(fasta.seq) + 1):(i * nrow(fasta.seq)), ])
    }))
    seq.bin <- do.call("cbind", poswise.amino.list)
    seq.len <- ncol(fasta.seq)
    list(fasta = fasta_raw_data,
         fasta_full_seq = fasta_full_seq,
         position_in_aligned_seq = position_in_aligned_seq,
         fasta.seq = fasta.seq,
         amino.names = amino.names,
         seq.bin = seq.bin,
         seq.len = seq.len)
}

