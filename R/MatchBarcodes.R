#' Match given cell barcodes to closest whitelist barcode.
#' @param x dataframe containing at least the column "CB" with cell barcodes
#' @param WhiteList vector of whitelist cell barcodes
#' @param max_dist maximum hamming distance (mismatches) allowed during barcode matching against whitelist. Defaults to 1 for consistency with the cellranger pipeline for atac-seq
#' @param ncores number of cores to use for parallelisation
#' @return dataframe x with added column "WhiteListMatch"
#'
#'
#'
MatchBarcodes <- function(x, WhiteList, max_dist = 1L, ncores = 16L) {
    stopifnot(is.data.frame(x))
    stopifnot(all(c("CB", "R1", "R2") %in% colnames(x)))
    stopifnot(is.character(x$CB))
    stopifnot(is.character(WhiteList))
    
    CB <- x$CB
    N  <- length(CB)
    
    if (!length(CB) || !length(WhiteList)) {
        x$WhiteListMatch <- NA_character_
        return(x)
    }
    
    ## ---- sanity: Hamming requires fixed length -------------------------------
    w1 <- nchar(CB)
    w2 <- nchar(WhiteList)
    if (length(unique(w1)) != 1L || length(unique(w2)) != 1L || w1[1L] != w2[1L]) {
        stop("All barcodes in CB and WhiteList must have the same length for Hamming distance.")
    }
    L <- w1[1L]
    
    ## ---- chunk boundaries (k = max_dist + 1) ---------------------------------
    split_bounds <- function(L, k = 2L) {
        part_lengths <- rep(floor(L / k), k)
        if (L %% k > 0L) {
            part_lengths[seq_len(L %% k)] <- part_lengths[seq_len(L %% k)] + 1L
        }
        start <- cumsum(c(1L, head(part_lengths, -1L)))
        end   <- cumsum(part_lengths)
        list(start = start, end = end)
    }
    
    k <- max_dist + 1L
    b <- split_bounds(L, k)
    
    ## ---- build WhiteList index: for each part, segment -> y indices ----------
    WL_index <- vector("list", k)
    nY <- length(WhiteList)
    
    for (p in seq_len(k)) {
        segs <- substr(WhiteList, b$start[p], b$end[p])
        WL_index[[p]] <- split(seq_len(nY), segs)
    }
    
    ## ---- helper: match one CB index --------------------------
    match_one_cb <- function(i) {
        cand_list <- vector("list", k)
        for (p in seq_len(k)) {
            seg <- substr(CB[i], b$start[p], b$end[p])
            idx_vec <- WL_index[[p]][[seg]]
            if (!is.null(idx_vec)) cand_list[[p]] <- idx_vec
        }
        
        cand_idx <- unique(unlist(cand_list, use.names = FALSE))
        if (length(cand_idx) == 0L) return(NA_character_)
        
        d <- stringdist(CB[i], WhiteList[cand_idx], method = "hamming")
        
        ## old nested-if logic
        if (min(d) > max_dist) return(NA_character_)
        
        chosen <- NA_integer_
        if (sum(d == 0L) == 1L)      chosen <- which(d == 0L)
        else if (sum(d == 1L) == 1L) chosen <- which(d == 1L)
        else return(NA_character_)
        
        WhiteList[cand_idx[chosen]]
    }
    
    ## ---- MAIN LOOP + progress reporting every 100M reads ----
    WhiteListMatch <- rep(NA_character_, N)
    report_every <- 1000000L
    
    if (ncores <= 1L) {
        
        ## ---- SERIAL path: progress reporting works perfectly ----
        for (i in seq_len(N)) {
            WhiteListMatch[i] <- match_one_cb(i)
            
            if (i %% report_every == 0L) {
                message("Processed ", format(i, big.mark=","), " barcodes...")
            }
        }
        
    } else {
        
        ## ---- PARALLEL path: chunk-level progress (optional) ----
        # Currently no per-read progress in parallel (stdout not safe)
        ncores_use <- min(ncores, N)
        grp_id <- ((seq_len(N) - 1L) %% ncores_use) + 1L
        groups <- split(seq_len(N), grp_id)
        
        worker <- function(idxs) vapply(idxs, match_one_cb, character(1))
        
        res_groups <- mclapply(groups, worker, mc.cores = ncores_use)
        
        # stitch back into correct order
        for (g in seq_along(groups)) {
            WhiteListMatch[groups[[g]]] <- res_groups[[g]]
        }
    }
    
    x$WhiteListMatch <- WhiteListMatch
    x
}
