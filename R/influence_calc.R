influence_calc <- function( up.sig,
                            dn.sig,
                            up.ratio,
                            dn.ratio,
                            tot.thr,
                            input.dir = "init",
                            output.dir = "output",
                            cscore,
                            write.name = NULL ){

  #----- Load packages
  cat("Loading Required packages ...\n")
  require( readr )
  require( dplyr )
  require( rhdf5 )
  require( xlsx )
  require( igraph )

  #----- init
  options(warn=-1)

  #----- check files
  cat(sprintf( "Loading required files from %s ...\n", input.dir))
  if( any( c( sprintf( "HG-U133_Plus_2.na%s.annot_for_mk_gmt.RData", 34 ),
              "inst.info.with.HPinfo.RData",
              "DEG.L1000.zspc.978.RData",
              sprintf( 'DEG.L1000.zspc.22268.%s.RData', 2:134 ) ) %in% list.files( input.dir ) == FALSE ) ){
    cat( "Following files do NOT existed:")
    c( sprintf( "HG-U133_Plus_2.na%s.annot_for_mk_gmt.RData", 34 ),
       "inst.info.with.HPinfo.RData",
       "DEG.L1000.zspc.978.RData",
       sprintf( 'DEG.L1000.zspc.22268.%s.RData', 2:134 ) )[ !( c( sprintf( "HG-U133_Plus_2.na%s.annot_for_mk_gmt.RData", 34 ),
                                                                  "inst.info.with.HPinfo.RData",
                                                                  "DEG.L1000.zspc.978.RData",
                                                                  sprintf( 'DEG.L1000.zspc.22268.%s.RData', 2:134 ) ) %in% list.files( input.dir ) ) ] %>% print
    cat( "Please select the appropreate directory or re-run init functions\n")
    return( invisible() )
  }

  #----- basic preprocess
  cat("Processing input files ...\n")

  # check args
  if( any(duplicated( up.sig )) ){ stop( "Up signature contains duplicated IDs.")}
  if( any(duplicated( dn.sig )) ){ stop( "Down signature contains duplicated IDs.")}

  # get
  load( file = file.path( input.dir,
                          sprintf( "HG-U133_Plus_2.na%s.annot_for_mk_gmt.RData", 34 ) ) )
  # removing multiple genes detecting probes
  dupli.probe <- anno$Probe.ID[ duplicated( anno$Probe.ID ) ] %>% unique
  anno <- anno[ !( anno$Probe.ID %in% dupli.probe ), ]

  #-------- Init process for UP and DOWN ratio
  up.probe.ratio <- data.frame( probe.id = up.sig,
                                ratio = up.ratio,
                                stringsAsFactors = FALSE ) %>%
    inner_join( x = ., y = anno[ c( "Probe.ID", "NAME" ) ], by = c( "probe.id" = "Probe.ID" ) ) %>%
    dplyr::filter( NAME != "---" ) %>% # remove non annotated probe
    group_by( NAME ) %>%
    summarise( mean.ratio   = mean( ratio   ),
               median.ratio = median( ratio ) )
  row.names( up.probe.ratio ) <- as.vector( up.probe.ratio$NAME )
  dn.probe.ratio <- data.frame( probe.id = dn.sig,
                                ratio = dn.ratio,
                                stringsAsFactors = FALSE ) %>%
    inner_join( x = ., y = anno[ c( "Probe.ID", "NAME" ) ], by = c( "probe.id" = "Probe.ID" ) ) %>%
    dplyr::filter( NAME != "---" ) %>% # remove non annotated probe
    group_by( NAME ) %>%
    summarise( mean.ratio   = mean( ratio   ),
               median.ratio = median( ratio ) )
  row.names( dn.probe.ratio ) <- as.vector( dn.probe.ratio$NAME )

  # Unify
  updn.probe.ratio <- rbind( up.probe.ratio, dn.probe.ratio )

  # Remove duplicated genes (due to probe level deg judgement)
  dup.gs <- as.vector( updn.probe.ratio$NAME )[duplicated( as.vector( updn.probe.ratio$NAME ))]
  updn.probe.ratio <- updn.probe.ratio[ !(updn.probe.ratio$NAME %in% dup.gs), ];row.names( updn.probe.ratio ) <- as.vector( updn.probe.ratio$NAME )

  # Separate up and down
  up.gene <- updn.probe.ratio %>%
    dplyr::filter( updn.probe.ratio$NAME %in% up.probe.ratio$NAME ) %>% dplyr::select( NAME ) %>% unlist %>% as.vector %>% unique
  dn.gene <- updn.probe.ratio %>%
    dplyr::filter( updn.probe.ratio$NAME %in% dn.probe.ratio$NAME ) %>% dplyr::select( NAME ) %>% unlist %>% as.vector %>% unique

  # inst info
  load( file = file.path( input.dir, "inst.info.with.HPinfo.RData" ) )
  inst.info <- inst.info.wt.hp

  # zspc, DEG
  load( file = file.path( input.dir, "DEG.L1000.zspc.978.RData" ))
  row.names( DEG.L1000.zspc.978 ) <- DEG.L1000.zspc.978 %>% dplyr::select(row_names) %>% unlist %>% as.vector

  #-------- Construct Gene Regulatory Network
  cat("Making directed graph of DEGs ...\n")

  # Extract from total score and add info
  cutoff.tot <- cscore %>%
    dplyr::filter( score > tot.thr ) %>%
    dplyr::select( distil_id, RANK.SCORE, score, up_score, down_score ) %>%
    inner_join( x = .,
                y = inst.info,
                by = "distil_id" ) %>%
    mutate( cutoff.tot = tot.thr,
            DEG.UP = ifelse( pert_desc.md %in% up.gene, TRUE, FALSE ),
            DEG.DN = ifelse( pert_desc.md %in% dn.gene, TRUE, FALSE ) )

  # Extract KD, OE (including mut and ligand), compound trt
  cutoff.tot.kd <- cutoff.tot[ cutoff.tot$pert_type == "trt_sh", ]
  cutoff.tot.oe <- cutoff.tot[ cutoff.tot$pert_type %in% c( "trt_lig", "trt_oe", "trt_oe.mut"), ]
  cutoff.tot.cp <- cutoff.tot[ cutoff.tot$pert_type == "trt_cp", ]

  # Extract data with the same direction (pert and DEG)
  cutoff.tot.kd.cor <- cutoff.tot.kd[ cutoff.tot.kd$DEG.DN == TRUE, ]
  cutoff.tot.oe.cor <- cutoff.tot.oe[ cutoff.tot.oe$DEG.UP == TRUE, ]
  cutoff.tot.kd.cor.uni.pert_desc <- unique( as.vector( cutoff.tot.kd.cor$pert_desc ) )
  cutoff.tot.oe.cor.uni.pert_desc <- unique( as.vector( cutoff.tot.oe.cor$pert_desc ) )

  # "first by" process (like SAS) against pert_id
  cutoff.tot.kd.cor<- cutoff.tot.kd.cor[ , c( "INDEX", "pert_desc.md", "score",
                                              "pert_id", "cell_id", "pert_time", "pert_time_unit", "pert_type",
                                              "DEG.UP", "DEG.DN",
                                              "distil_id", "pert_desc",
                                              "HP.Reagent.Name", "HP.Reagent.Type") ] %>%
    arrange( pert_id ) %>%
    arrange( desc( score ) )
  first.pert_id <- duplicated( cutoff.tot.kd.cor$pert_id, fromLast = FALSE )
  cutoff.tot.kd.cor <- cutoff.tot.kd.cor[ !(first.pert_id), ]
  cutoff.tot.oe.cor <- cutoff.tot.oe.cor[ c( "INDEX", "pert_desc.md", "score",
                                             "pert_id", "cell_id", "pert_time", "pert_time_unit", "pert_type",
                                             "DEG.UP", "DEG.DN",
                                             "distil_id", "pert_desc",
                                             "HP.Reagent.Name", "HP.Reagent.Type") ] %>%
    arrange( pert_id ) %>%
    arrange( desc( score ) )
  first.pert_id <- duplicated( cutoff.tot.oe.cor$pert_id, fromLast = FALSE )
  cutoff.tot.oe.cor <- cutoff.tot.oe.cor[ !(first.pert_id), ]

  # check
  if( all( nrow( cutoff.tot.kd.cor ) == 0, nrow(cutoff.tot.oe.cor) == 0 ) ){
    cat( "\nNo gene for making directed graph, due to no correspondance expression change between input and LINCS\n" )
    return( invisible() )
  }

  # Get upstream genes
  # KD: more than 2 pert_ids are required.
  kd.upstream.genes <- names( sort( table( cutoff.tot.kd.cor$pert_desc.md ), decreasing = TRUE ) )[ ( sort( table( cutoff.tot.kd.cor$pert_desc.md ), decreasing = TRUE ) ) >= 2 ]
  oe.upstream.genes <- names( sort( table( cutoff.tot.oe.cor$pert_desc.md ), decreasing = TRUE ) )

  # Get ZSC of upstream genes
  kd.upstream.genes.id <- cutoff.tot.kd.cor[ cutoff.tot.kd.cor$pert_desc.md %in% kd.upstream.genes, c( "pert_desc.md", "distil_id" ) ]
  oe.upstream.genes.id <- cutoff.tot.oe.cor[ cutoff.tot.oe.cor$pert_desc.md %in% oe.upstream.genes, c( "pert_desc.md", "distil_id" ) ]
  oe.upstream.deg.978 <- DEG.L1000.zspc.978[ as.vector( oe.upstream.genes.id$distil_id ) ]
  kd.upstream.deg.978 <- DEG.L1000.zspc.978[ as.vector( kd.upstream.genes.id$distil_id ) ]

  # count DEGs in LINCS
  n.deg.up.fun <- function( x ){ return( sum( x == 1  ) ) }
  n.deg.dn.fun <- function( x ){ return( sum( x == -1 ) ) }
  N.DEG.UP.oe.upstream.deg.978 <- apply( oe.upstream.deg.978, 2, n.deg.up.fun )
  N.DEG.DN.oe.upstream.deg.978 <- apply( oe.upstream.deg.978, 2, n.deg.dn.fun )
  N.DEG.UP.kd.upstream.deg.978 <- apply( kd.upstream.deg.978, 2, n.deg.up.fun )
  N.DEG.DN.kd.upstream.deg.978 <- apply( kd.upstream.deg.978, 2, n.deg.dn.fun )

  # percentage of input degs per LINCS978 DEGs
  pct.deg.up.fun <- function( x, names.x, cmpr.probe ){
    y <- names.x[ x == 1 ]
    y.cor <- sum( y %in% cmpr.probe )
    return( y.cor / length( y ) ) }
  pct.deg.dn.fun <- function( x, names.x, cmpr.probe ){
    y <- names.x[ x == -1 ]
    y.cor <- sum( y %in% cmpr.probe )
    return( y.cor / length( y ) ) }
  up.sig.978 <- up.sig[ up.sig %in% as.vector( unlist( DEG.L1000.zspc.978$row_names ) ) ]
  dn.sig.978 <- dn.sig[ dn.sig %in% as.vector( unlist( DEG.L1000.zspc.978$row_names ) ) ]
  N.COR.DEG.UP.oe.upstream.deg.978 <- apply( oe.upstream.deg.978,
                                             2,
                                             pct.deg.up.fun,
                                             names.x = row.names( oe.upstream.deg.978 ),
                                             cmpr.probe = up.sig.978 )
  N.COR.DEG.DN.oe.upstream.deg.978 <- apply( oe.upstream.deg.978,
                                             2,
                                             pct.deg.dn.fun,
                                             names.x = row.names( oe.upstream.deg.978 ),
                                             cmpr.probe = dn.sig.978 )
  N.COR.DEG.UP.kd.upstream.deg.978 <- apply( kd.upstream.deg.978,
                                             2,
                                             pct.deg.up.fun,
                                             names.x = row.names( kd.upstream.deg.978 ),
                                             cmpr.probe = up.sig.978 )
  N.COR.DEG.DN.kd.upstream.deg.978 <- apply( kd.upstream.deg.978,
                                             2,
                                             pct.deg.dn.fun,
                                             names.x = row.names( kd.upstream.deg.978 ),
                                             cmpr.probe = dn.sig.978 )

  # prepare for summary
  df.n.deg.up.oe <- data.frame( distil_id = names( N.DEG.UP.oe.upstream.deg.978 ),
                                LINCS978.UP = N.DEG.UP.oe.upstream.deg.978,
                                stringsAsFactors = FALSE )
  df.n.deg.dn.oe <- data.frame( distil_id = names( N.DEG.DN.oe.upstream.deg.978 ),
                                LINCS978.DOWN = N.DEG.DN.oe.upstream.deg.978,
                                stringsAsFactors = FALSE )
  df.n.deg.up.kd <- data.frame( distil_id = names( N.DEG.UP.kd.upstream.deg.978 ),
                                LINCS978.UP = N.DEG.UP.kd.upstream.deg.978,
                                stringsAsFactors = FALSE )
  df.n.deg.dn.kd <- data.frame( distil_id = names( N.DEG.DN.kd.upstream.deg.978 ),
                                LINCS978.DOWN = N.DEG.DN.kd.upstream.deg.978,
                                stringsAsFactors = FALSE )
  df.pct.deg.up.oe <- data.frame( distil_id = names( N.COR.DEG.UP.oe.upstream.deg.978 ),
                                  COMMON.CP.LINCS978.UP = N.COR.DEG.UP.oe.upstream.deg.978,
                                  stringsAsFactors = FALSE )
  df.pct.deg.dn.oe <- data.frame( distil_id = names( N.COR.DEG.DN.oe.upstream.deg.978 ),
                                  COMMON.CP.LINCS978.DOWN = N.COR.DEG.DN.oe.upstream.deg.978,
                                  stringsAsFactors = FALSE )
  df.pct.deg.up.kd <- data.frame( distil_id = names( N.COR.DEG.UP.kd.upstream.deg.978 ),
                                  COMMON.CP.LINCS978.UP = N.COR.DEG.UP.kd.upstream.deg.978,
                                  stringsAsFactors = FALSE )
  df.pct.deg.dn.kd <- data.frame( distil_id = names( N.COR.DEG.DN.kd.upstream.deg.978 ),
                                  COMMON.CP.LINCS978.DOWN = N.COR.DEG.DN.kd.upstream.deg.978,
                                  stringsAsFactors = FALSE )

  # merge
  if( nrow( kd.upstream.genes.id ) >= 1 ){
    kd.upstream.genes.id <- kd.upstream.genes.id %>%
      left_join( x = ., y = df.n.deg.up.kd, by = "distil_id" ) %>%
      left_join( x = ., y = df.n.deg.dn.kd, by = "distil_id" ) %>%
      left_join( x = ., y = df.pct.deg.up.kd, by = "distil_id" ) %>%
      left_join( x = ., y = df.pct.deg.dn.kd, by = "distil_id" )
  }
  if( nrow( oe.upstream.genes.id ) >= 1 ){
    oe.upstream.genes.id <- oe.upstream.genes.id %>%
      left_join( x = ., y = df.n.deg.up.oe, by = "distil_id" ) %>%
      left_join( x = ., y = df.n.deg.dn.oe, by = "distil_id" ) %>%
      left_join( x = ., y = df.pct.deg.up.oe, by = "distil_id" ) %>%
      left_join( x = ., y = df.pct.deg.dn.oe, by = "distil_id" )
  }

  # check
  if( nrow( kd.upstream.genes.id ) == 0 ){
    cat( "\nNo gene knockdowned in LINCS dataset remained.\n")
  }
  if( nrow( oe.upstream.genes.id ) == 0 ){
    cat( "\nNo gene overexpressed in LINCS dataset remained\n")
  }
  if( all( nrow( kd.upstream.genes.id ) == 0, nrow( oe.upstream.genes.id ) == 0 ) ){
    cat( "\nNo gene for making directed graph!!\n" )
    return( invisible() )
  }

  #=== summary of downstream genes
  kd.downstream.genes.df <- data.frame()
  oe.downstream.genes.df <- data.frame()

  #- Get downstream genes
  if( length( kd.upstream.genes ) >= 1 ){
    for( kd.upstream.gene in kd.upstream.genes ){
      # debug kd.upstream.gene <- kd.upstream.genes[1]

      # down gene (sh) -including pert_id
      kd.upstream.index.s <- cutoff.tot.kd.cor %>%
        dplyr::filter( pert_desc.md == kd.upstream.gene ) %>%
        dplyr::select( INDEX ) %>%
        unlist %>%
        as.vector

      # iterative compute
      TMP.COUNT <- c( 0, seq( 10000, 1320000, by = 10000 ), 1328098 )
      tmp.index.s <- as.numeric( gsub("^L", "", kd.upstream.index.s ) )
      tmp.df.s <- ceiling( tmp.index.s / 10000 ) + 1
      uni.df.numbers <- unique( tmp.df.s )

      for( uni.df.number in uni.df.numbers ){

        # get deg deg object
        eval( parse( text = sprintf( "load( file = file.path( input.dir, 'DEG.L1000.zspc.22268.%s.RData' ))",
                                     uni.df.number )))
        # rename
        eval( parse( text = sprintf( "DEG.L1000.zspc.22268 <- DEG.L1000.zspc.22268.%s; rm( DEG.L1000.zspc.22268.%s );gc <- gc(F,T)",
                                     uni.df.number, uni.df.number)))
        row.names( DEG.L1000.zspc.22268 ) <- DEG.L1000.zspc.22268 %>%
          dplyr::select( row_names ) %>%
          unlist() %>%
          as.vector()
        DEG.L1000.zspc.22268  <- DEG.L1000.zspc.22268 %>%
          dplyr::select( -( row_names ) )
        tmp.col.names <- colnames( DEG.L1000.zspc.22268 )

        # make INDEX
        tmp.inst.info <- inst.info.wt.hp %>%
          dplyr::filter( distil_id %in% tmp.col.names)
        tmp.col   <- as.vector( tmp.inst.info$distil_id )
        tmp.index <- as.vector( tmp.inst.info$INDEX )
        DEG.L1000.zspc.22268             <- DEG.L1000.zspc.22268[ tmp.col ]
        colnames( DEG.L1000.zspc.22268 ) <- tmp.index

        # extract index
        tmp.kd.upstream.index.s <- kd.upstream.index.s[ tmp.df.s == uni.df.number ]
        for( tmp.kd.upstream.index in tmp.kd.upstream.index.s ){

          # up/down-kd.downstream.genes
          tmp.deg <- DEG.L1000.zspc.22268[ tmp.kd.upstream.index ]
          colnames( tmp.deg ) <- "UPSTREAM"
          up_downstream.genes <- unique( as.vector( anno[ anno$Probe.ID %in% row.names( tmp.deg )[ tmp.deg$UPSTREAM == 1 ], "NAME" ] ) )
          dn_downstream.genes <- unique( as.vector( anno[ anno$Probe.ID %in% row.names( tmp.deg )[ tmp.deg$UPSTREAM == -1 ], "NAME" ] ) )
          up_downstream.genes <- up_downstream.genes[ up_downstream.genes != "---" ]
          dn_downstream.genes <- dn_downstream.genes[ dn_downstream.genes != "---" ]

          # summary
          one.kd.downstream.genes.df <- rbind( data.frame( upstream.gene = rep( kd.upstream.gene, length( up_downstream.genes )),
                                                           dnstream.gene = up_downstream.genes,
                                                           dn.direction = "UP" ),
                                               data.frame( upstream.gene = rep( kd.upstream.gene, length( dn_downstream.genes )),
                                                           dnstream.gene = dn_downstream.genes,
                                                           dn.direction = "DOWN" ) )

          # binding
          kd.downstream.genes.df <- rbind( kd.downstream.genes.df,
                                           one.kd.downstream.genes.df )
        }
      }

      rm( DEG.L1000.zspc.22268 );gc <- gc( F, T)

    }
  }

  # extract down stream genes of oe-upstream
  if( length( oe.upstream.genes ) >= 1 ){
    for( oe.upstream.gene in oe.upstream.genes ){

      # get pert id including up gene (oe/mut/ligand)
      oe.upstream.index.s <- cutoff.tot.oe.cor %>%
        dplyr::filter( pert_desc.md == oe.upstream.gene ) %>%
        dplyr::select( INDEX ) %>%
        unlist %>%
        as.vector

      # iterative compute
      TMP.COUNT <- c( 0, seq( 10000, 1320000, by = 10000 ), 1328098 )
      tmp.index.s <- as.numeric( gsub("^L", "", oe.upstream.index.s ) )
      tmp.df.s <- ceiling( tmp.index.s / 10000 ) + 1
      uni.df.numbers <- unique( tmp.df.s )

      for( uni.df.number in uni.df.numbers ){

        # load deg deg object
        eval( parse( text = sprintf( "load( file = file.path( input.dir, 'DEG.L1000.zspc.22268.%s.RData' ))",
                                     uni.df.number )))
        # rename
        eval( parse( text = sprintf( "DEG.L1000.zspc.22268 <- DEG.L1000.zspc.22268.%s; rm( DEG.L1000.zspc.22268.%s );gc <- gc(F,T)",
                                     uni.df.number, uni.df.number)))
        row.names( DEG.L1000.zspc.22268 ) <- DEG.L1000.zspc.22268 %>%
          dplyr::select( row_names ) %>%
          unlist() %>%
          as.vector()
        DEG.L1000.zspc.22268  <- DEG.L1000.zspc.22268 %>%
          dplyr::select( -( row_names ) )
        tmp.col.names <- colnames( DEG.L1000.zspc.22268 )

        # make index
        tmp.inst.info <- inst.info.wt.hp %>%
          dplyr::filter( distil_id %in% tmp.col.names)
        tmp.col   <- as.vector( tmp.inst.info$distil_id )
        tmp.index <- as.vector( tmp.inst.info$INDEX )
        DEG.L1000.zspc.22268             <- DEG.L1000.zspc.22268[ tmp.col ]
        colnames( DEG.L1000.zspc.22268 ) <- tmp.index

        # extract index
        tmp.oe.upstream.index.s <- oe.upstream.index.s[ tmp.df.s == uni.df.number ]
        for( tmp.oe.upstream.index in tmp.oe.upstream.index.s ){

          # extract up/down-downstream.genes
          tmp.deg <- DEG.L1000.zspc.22268[ tmp.oe.upstream.index ]
          colnames( tmp.deg ) <- "UPSTREAM"
          up_downstream.genes <- unique( as.vector( anno[ anno$Probe.ID %in% row.names( tmp.deg )[ tmp.deg$UPSTREAM == 1 ], "NAME" ] ) )
          dn_downstream.genes <- unique( as.vector( anno[ anno$Probe.ID %in% row.names( tmp.deg )[ tmp.deg$UPSTREAM == -1 ], "NAME" ] ) )
          up_downstream.genes <- up_downstream.genes[ up_downstream.genes != "---" ]
          dn_downstream.genes <- dn_downstream.genes[ dn_downstream.genes != "---" ]

          # summary
          one.oe.downstream.genes.df <- rbind( data.frame( upstream.gene = rep( oe.upstream.gene, length( up_downstream.genes )),
                                                           dnstream.gene = up_downstream.genes,
                                                           dn.direction = "UP" ),
                                               data.frame( upstream.gene = rep( oe.upstream.gene, length( dn_downstream.genes )),
                                                           dnstream.gene = dn_downstream.genes,
                                                           dn.direction = "DOWN" ) )

          # binding
          oe.downstream.genes.df <- rbind( oe.downstream.genes.df,
                                           one.oe.downstream.genes.df )
        }
      }
      rm( DEG.L1000.zspc.22268 );gc <- gc( F, T)

    }
  }

  # KD: more than 2 pert_ids are required.
  kd.downstream.genes.df.2 <- kd.downstream.genes.df[ duplicated( kd.downstream.genes.df ), ] %>% unique

  # UP: unique
  oe.downstream.genes.df.2 <- oe.downstream.genes.df %>% unique

  # get the same direction between lincs and in house degs.
  kd.downstream.genes.df.3 <- kd.downstream.genes.df.2 %>% dplyr::filter( (dnstream.gene %in% up.gene & dn.direction == "UP" ) | ( dnstream.gene %in% dn.gene & dn.direction == "DOWN" ) )
  oe.downstream.genes.df.3 <- oe.downstream.genes.df.2 %>% dplyr::filter( (dnstream.gene %in% up.gene & dn.direction == "UP" ) | ( dnstream.gene %in% dn.gene & dn.direction == "DOWN" ) )

  if( nrow( kd.downstream.genes.df.3 ) >= 1 ){

    # add upstream gene's FC
    kd.downstream.genes.df.3 <- kd.downstream.genes.df.3 %>%
      left_join( x = ., y = updn.probe.ratio, by = c( "upstream.gene" = "NAME")) %>%
      left_join( x = ., y = kd.upstream.genes.id, by = c( "upstream.gene" = "pert_desc.md" ) )

    # make Weight
    Weight <- ifelse( kd.downstream.genes.df.3$dn.direction == "UP",
                      ( 1 / kd.downstream.genes.df.3$median.ratio ) * kd.downstream.genes.df.3$COMMON.CP.LINCS978.UP,
                      ( 1 / kd.downstream.genes.df.3$median.ratio ) * kd.downstream.genes.df.3$COMMON.CP.LINCS978.DOWN )
    kd.downstream.genes.df.3 <- kd.downstream.genes.df.3 %>% mutate( weight = Weight ) %>%
      dplyr::select( upstream.gene, dnstream.gene, weight, median.ratio, dn.direction,
              LINCS978.UP, LINCS978.DOWN, COMMON.CP.LINCS978.UP, COMMON.CP.LINCS978.DOWN )

    # deconstract for 4 network making
    kd.downstream.up.genes.df.3 <- kd.downstream.genes.df.3 %>% dplyr::filter( dn.direction == "UP" )
    kd.downstream.dn.genes.df.3 <- kd.downstream.genes.df.3 %>% dplyr::filter( dn.direction == "DOWN" )

    # remove multiple edges and self-loops
    if( nrow( kd.downstream.up.genes.df.3 ) >= 1 ){
      igraph.kd.up <- graph.data.frame( kd.downstream.up.genes.df.3[1:3], directed = TRUE )
      igraph.kd.up <- simplify( igraph.kd.up,
                                remove.loops = TRUE,
                                remove.multiple = TRUE )
    }

    if( nrow( kd.downstream.dn.genes.df.3 ) >= 1 ){
      igraph.kd.dn <- graph.data.frame( kd.downstream.dn.genes.df.3[1:3], directed = TRUE )
      igraph.kd.dn <- simplify( igraph.kd.dn,
                                remove.loops = TRUE,
                                remove.multiple = TRUE )
    }
  }

  if( nrow( oe.downstream.genes.df.3 ) >= 1 ){

    oe.downstream.genes.df.3 <- oe.downstream.genes.df.3 %>%
      left_join( x = ., y = updn.probe.ratio, by = c( "upstream.gene" = "NAME")) %>%
      left_join( x = ., y = oe.upstream.genes.id, by = c( "upstream.gene" = "pert_desc.md" ) )

    # make Weight
    Weight <- ifelse( oe.downstream.genes.df.3$dn.direction == "UP",
                      oe.downstream.genes.df.3$median.ratio * oe.downstream.genes.df.3$COMMON.CP.LINCS978.UP,
                      oe.downstream.genes.df.3$median.ratio * oe.downstream.genes.df.3$COMMON.CP.LINCS978.DOWN )
    oe.downstream.genes.df.3 <- oe.downstream.genes.df.3 %>% mutate( weight = Weight ) %>%
      dplyr::select( upstream.gene, dnstream.gene, weight, median.ratio, dn.direction,
              LINCS978.UP, LINCS978.DOWN, COMMON.CP.LINCS978.UP, COMMON.CP.LINCS978.DOWN )

    # collapse
    oe.downstream.up.genes.df.3 <- oe.downstream.genes.df.3 %>% dplyr::filter( dn.direction == "UP" )
    oe.downstream.dn.genes.df.3 <- oe.downstream.genes.df.3 %>% dplyr::filter( dn.direction == "DOWN" )

  }

  # Join network
  if( nrow( kd.downstream.genes.df.2 ) >= 1 ){
    kd.upstream.all.genes <- kd.downstream.genes.df.2 %>% dplyr::select( upstream.gene ) %>% unlist %>% as.vector %>% unique
  }else{
    kd.upstream.all.genes <- NULL
  }
  if( nrow( oe.downstream.genes.df.2 ) >= 1 ){
    oe.upstream.all.genes <- oe.downstream.genes.df.2 %>% dplyr::select( upstream.gene ) %>% unlist %>% as.vector %>% unique
  }else{
    oe.upstream.all.genes <- NULL
  }

  # remove multiple edge and self loops
  if( length( grep( "^kd.downstream.up.genes.df.3$", objects() ) ) >= 1 ){
    if( nrow( kd.downstream.up.genes.df.3 ) >= 1 ){
      tmp.kd.up <- kd.downstream.up.genes.df.3[1:3]
      tmp.kd.up <- dplyr::filter( tmp.kd.up, weight != 0 )
      igraph.tmp.kd.up <- graph.data.frame( tmp.kd.up, directed = TRUE )
      igraph.tmp.kd.up <- simplify( igraph.tmp.kd.up,
                                    remove.loops = TRUE,
                                    remove.multiple = TRUE )
    }else{
      tmp.kd.up <- NULL
      igraph.tmp.kd.up <- NULL
    }
  }else{ tmp.kd.up <- NULL;igraph.tmp.kd.up <- NULL }
  if( length( grep( "^kd.downstream.dn.genes.df.3$", objects() ) ) >= 1 ){
    if( nrow( kd.downstream.dn.genes.df.3 ) >= 1 ){
      tmp.kd.dn <- kd.downstream.dn.genes.df.3[1:3]
      tmp.kd.dn <- dplyr::filter( tmp.kd.dn, weight != 0 )
      igraph.tmp.kd.dn <- graph.data.frame( tmp.kd.dn, directed = TRUE )
      igraph.tmp.kd.dn <- simplify( igraph.tmp.kd.dn,
                                    remove.loops = TRUE,
                                    remove.multiple = TRUE )
    }else{
      tmp.kd.dn <- NULL
      igraph.tmp.kd.dn <- NULL
    }
  }else{ tmp.kd.dn <- NULL;igraph.tmp.kd.dn <- NULL }
  if( length( grep( "^oe.downstream.up.genes.df.3$", objects() ) ) >= 1 ){
    if( nrow( oe.downstream.up.genes.df.3 ) >= 1 ){
      tmp.oe.up <- oe.downstream.up.genes.df.3[1:3]
      tmp.oe.up <- dplyr::filter( tmp.oe.up, weight != 0 )
      igraph.tmp.oe.up <- graph.data.frame( tmp.oe.up, directed = TRUE )
      igraph.tmp.oe.up <- simplify( igraph.tmp.oe.up,
                                    remove.loops = TRUE,
                                    remove.multiple = TRUE )
    }else{
      tmp.oe.up <- NULL
      igraph.tmp.oe.up <- NULL
    }
  }else{ tmp.oe.up <- NULL;igraph.tmp.oe.up <- NULL }
  if( length( grep( "^oe.downstream.dn.genes.df.3$", objects() ) ) >= 1 ){
    if( nrow( oe.downstream.dn.genes.df.3 ) >= 1 ){
      tmp.oe.dn <- oe.downstream.dn.genes.df.3[1:3]
      tmp.oe.dn <- dplyr::filter( tmp.oe.dn, weight != 0 )
      igraph.tmp.oe.dn <- graph.data.frame( tmp.oe.dn, directed = TRUE )
      igraph.tmp.oe.dn <- simplify( igraph.tmp.oe.dn,
                                    remove.loops = TRUE,
                                    remove.multiple = TRUE )

    }else{
      tmp.oe.dn <- NULL
      igraph.tmp.oe.dn <- NULL
    }
  }else{ tmp.oe.dn <- NULL;igraph.tmp.oe.dn <- NULL }

  # Unify Network
  igraph.all <-
    rbind( tmp.kd.up,
           tmp.kd.dn,
           tmp.oe.up,
           tmp.oe.dn ) %>%
    graph.data.frame( ., directed = TRUE ) %>%
    simplify( .,
              remove.loops = TRUE,
              remove.multiple = TRUE )

  if( nrow( rbind( tmp.kd.up,
                   tmp.kd.dn,
                   tmp.oe.up,
                   tmp.oe.dn ) ) == 0 ){ stop("No Network was developed by the input conditions!")}

  if( length( write.name ) == 1 ){
    # Save for igraph
    saveRDS(igraph.all, file = file.path( output.dir , paste0(write.name, ".rds") ))

    # Save for cytoscape
    require( NetPathMiner )
    plotCytoscapeGML(igraph.all, file = file.path( output.dir , paste0(write.name, ".gml") ))
  }

  # Kleinberg's hub scores.
  cat("Calculating Kleinberg's hub score from directed graph of DEGs ...\n")
  hub.igraph.all <- hub.score( igraph.all )$vector %>% sort(., decreasing = TRUE )
  res <- data.frame( UPSTREAM.GENE = names( hub.igraph.all ),
                     HUB.SCORE = hub.igraph.all,
                     stringsAsFactors = FALSE ) %>%
    left_join( x = ., y = updn.probe.ratio, by = c( "UPSTREAM.GENE" = "NAME" ) )
  return( res )

  cat("Done.\n")
}

