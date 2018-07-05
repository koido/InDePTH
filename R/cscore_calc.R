## LINCS CMap method Analysis
cscore_LINCS <- function( up.sig,
                          dn.sig,
                          input.dir = "init",
                          output.dir = "output",
                          write.name = NULL,
                          gene.symbol = FALSE){

  #----- Reference
  cat("=====")
  cat("Script of Cscore calculation was quoted from:\n")
  cat("Ushijima, Masaru, et al. 'Development of a gene expression database and related analysis programs for evaluation of anticancer compounds.' Cancer science 104.3 (2013): 360-368.")
  cat("\nAlso see\n")
  cat("http://scads.jfcr.or.jp/db/cs/index.html#Rcode\n")
  cat("http://www.broadinstitute.org/cmap/\n")
  cat("=====\n\n")

  #----- Load packages
  cat("Loading Required packages ...\n")
  require( readr )
  require( dplyr )
  require( rhdf5 )
  require( xlsx )
  require( igraph )
  require( compiler )

  #----- functions

  # direction check
  calc.pre.score <- function(up.score, dn.score) {
    ifelse(up.score * dn.score > 0.0, 0, up.score - dn.score)
  }
  calc.pre.score <- cmpfun(calc.pre.score)

  # Normalize ks value
  calc.normalized.score <- function(up.score, dn.score) {

    pre.score <- calc.pre.score(up.score, dn.score)
    pre.max <- max(pre.score)
    pre.min <- min(pre.score)
    out <- numeric(length(pre.score))
    out[pre.score > 0] <- pre.score[pre.score > 0] / pre.max
    out[pre.score < 0] <- pre.score[pre.score < 0] / (-pre.min)
    return( out )
  }
  calc.normalized.score <- cmpfun(calc.normalized.score)

  # sort by ranking
  sort.rank <- function(rank.vec) {
    rank.vec[order(rank.vec)]
  }
  sort.rank <- cmpfun(sort.rank)

  # v = rankmat.*.sort
  calc.score <- function(rank.sort, nprobe) {
    # get UP or DOWN score
    a.vec <- apply(rank.sort,
                   2,
                   function(v){
                     nv <- length(v)
                     return( max((1:nv) / nv - v[1:nv] / nprobe))
                   })
    b.vec <- apply(rank.sort,
                   2,
                   function(v) {
                     nv <- length(v)
                     return( max(v[1:nv] / nprobe - (1:nv - 1) / nv) )
                   })

    # set ksi = a, if a > b. Set ksi = -b if b > a
    ifelse(a.vec > b.vec,
           a.vec,
           -b.vec)
  }
  calc.score <- cmpfun(calc.score)

  # make rank
  makerank <- function(s1, s2) {
    d1 <- calc.pre.score(s1, s2)

    count.lo <- sum(d1 == 0)
    count.dn <- sum(d1 < 0)
    d2 <- s1[d1==0] - s2[d1==0]

    # d1 rank
    r <- rank(d1, ties.method="min")
    # d2 rank
    rsub <- rank(d2, ties.method="random") + count.dn

    # 0 rank
    if(count.lo != 0){
      r[r == (count.dn + 1)] <- rsub
    }

    # redefine
    res <- length(s1) + 1 - r
    return( res )
  }
  makerank <- cmpfun(makerank)

  ## calculation of connectivity scores
  calc.cscore <- function(up.sig, dn.sig, rank.mat, N.ROUND = 5 ){

    # make rank matrix
    rankmat.up <- rank.mat[up.sig,]
    rankmat.dn <- rank.mat[dn.sig,]

    # get order
    rankmat.up.sort <- apply(rankmat.up, 2, sort)
    rankmat.dn.sort <- apply(rankmat.dn, 2, sort)

    # get ks value
    up.score <- calc.score(rank.sort = rankmat.up.sort, nprobe = 978)
    dn.score <- calc.score(rank.sort = rankmat.dn.sort, nprobe = 978)

    # get normalized ks value
    norm.score <- calc.normalized.score(up.score, dn.score)

    # make rank
    rank.vec <- makerank(up.score, dn.score)

    # rounding
    res <- cbind(rank.vec, round(norm.score, N.ROUND), round(up.score, N.ROUND), round(dn.score, N.ROUND))
    colnames(res) <- c("rank", "score", "up_score", "down_score")
    return( res )
  }
  calc.cscore <- cmpfun(calc.cscore)

  ## main function
  cscore_analysis <- function(up.sig,
                              dn.sig,
                              rank.mat,
                              meta.CMAP){

    # calculation of connectivity scores
    res <- as.data.frame(calc.cscore(up.sig, dn.sig, rank.mat))
    res.out <- cbind(res[,1], meta.CMAP, res[,-1])[order(res[,1]),]
    colnames(res.out)[1] <- "RANK.SCORE"

    return( res.out )
  }
  cscore_analysis <- cmpfun(cscore_analysis)

  #----- Check input dataset
  cat(sprintf( "Checking required files from %s ...\n", input.dir))
  if( any( c( "RANK_MATRIX_sig978.RData",
              "inst.info.with.HPinfo.RData" ) %in% list.files( input.dir ) == FALSE ) ){
    cat( "Following files do NOT existed:")
    c( "RANK_MATRIX_sig978.RData",
       "inst.info.with.HPinfo.RData" )[ !( c( "RANK_MATRIX_sig978.RData",
                                              "inst.info.with.HPinfo.RData" ) %in% list.files( input.dir ) ) ] %>% print
    cat( "Please select the appropreate directory or re-run init functions\n")
    return( invisible() )
  }


  #----- Load dataset
  cat(sprintf( "Loading required files from %s ...\n", input.dir))
  # load rank matrix
  load( file.path(input.dir, "RANK_MATRIX_sig978.RData" ) )

  # load inst info
  load( file.path(input.dir, "inst.info.with.HPinfo.RData") )
  inst.info <- inst.info.wt.hp[ c( "distil_id", "INDEX",
                                   "pert_id", "pert_desc.md", "pert_type",
                                   "cell_id",
                                   "pert_time", "pert_time_unit", "pert_dose", "pert_dose_unit" ) ]

  # get probe ids
  PROBE.ID.L1000_sig978 <- rownames(rank.matrix.z)
  
  # gene symbol2probe id
  if(gene.symbol == TRUE){
      require( "hgu133plus2.db" )
      # Get the probe identifiers that are mapped to a gene symbol, hgu133
      mapped_probes <- mappedkeys( hgu133plus2SYMBOL )
      # Convert to a data.frame
      hgu133plus2SYMBOL.data <- as.data.frame(hgu133plus2SYMBOL[mapped_probes])
      
      # Limit: qc'd probe set (please see MATERIALS AND METHODS in InDePTH paper)
      hgu133plus2SYMBOL.data <- hgu133plus2SYMBOL.data %>% dplyr::filter(probe_id %in% qcd_probe)

      # update
      up.sig <- hgu133plus2SYMBOL.data %>% dplyr::filter(symbol %in% up.sig) %>% dplyr::select(probe_id) %>% unlist %>% as.vector 
      dn.sig <- hgu133plus2SYMBOL.data %>% dplyr::filter(symbol %in% dn.sig) %>% dplyr::select(probe_id) %>% unlist %>% as.vector 
  }

  # Limit probe ids for landmark genes
  up.sig.L1000_sig978 <- up.sig[up.sig %in% PROBE.ID.L1000_sig978]
  dn.sig.L1000_sig978 <- dn.sig[dn.sig %in% PROBE.ID.L1000_sig978]

  if(any(length(up.sig.L1000_sig978) <= 1, length(dn.sig.L1000_sig978) <= 1) == TRUE){
    stop("This analysis requires more than 2 UP or DOWN genes of LINCS 978 landmark genes, respectively.")
  }

  #----- RUN
  cat("Calculating CMap Score of 1.3 million LINCS perturbations. Please wait for about more than ten minutes... \n")
  res <- cscore_analysis(up.sig = up.sig.L1000_sig978,
                         dn.sig = dn.sig.L1000_sig978,
                         rank.mat = rank.matrix.z,
                         meta.CMAP = inst.info )
  if( length( write.name ) == 1 ){
    dir.create( path = output.dir, recursive = TRUE )
    write.table( res, file.path( output.dir, write.name ), sep = "\t", row.names = F, col.names = T )
  }

  cat("Done.\n")
  return( res )
}
