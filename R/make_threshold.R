make_threshold <- function( pert_id_vec = c( "TRCN0000010389", "TRCN0000010390", "TRCN0000010391"),
                            cscore = res1,
                            cell = "HT29",
                            input.dir = "init",
                            output.dir = "output" ){

  #----- Load packages
  cat("Loading Required packages ...\n")
  require( dplyr )
  require( pROC )

  #----- Check input dataset
  cat(sprintf( "Checking required files from %s ...\n", input.dir))
  if( any( c( "inst.info.with.HPinfo.RData" ) %in% list.files( input.dir ) == FALSE ) ){
    cat( "Following files do NOT existed:")
    c( "inst.info.with.HPinfo.RData" )[ !( c( "inst.info.with.HPinfo.RData" ) %in% list.files( input.dir ) ) ] %>% print
    cat( "Please select the appropreate directory or re-run init functions\n")
    return( invisible() )
  }

  #----- Load dataset
  cat(sprintf( "Loading required files from %s ...\n", input.dir))

  # load inst info
  load( file.path(input.dir, "inst.info.with.HPinfo.RData") )
  inst.info <- inst.info.wt.hp[ c( "distil_id", "pert_id", "cell_id" ) ]

  #----- RUN
  cat("Processing files... \n")

  # limit cell line
  inst.info <- inst.info %>% dplyr::filter( cell_id %in% cell ) %>% dplyr::select( -(cell_id))

  # add true info
  inst.info <- inst.info %>% dplyr::mutate( outcome = as.factor(as.numeric( pert_id %in% pert_id_vec )))

  # merge cscore
  inst.info <- inst.info %>% inner_join( x = ., y = cscore, by = "distil_id" )

  # ROC analysis
  cat("ROC analysis... \n")
  ROC_res <- roc( formula = outcome ~ score, data = inst.info, direction = "<" )
  c_index <- ROC_res$auc
  cutoff <- coords( ROC_res,
                    x = "best",　# the point with the best sum of sensitivity and specificity
                    ret = c( "threshold", "sensitivity", "specificity", "ppv", "npv"),
                    best.method = "youden" )


  cat("Done.\n")
  return( list( c_index = as.vector(c_index), cutoff = as.vector(cutoff[1] ) ) )
}
