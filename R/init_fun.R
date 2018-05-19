init <- function( input.dir = "input",
                  out.dir  = "init" ) {

  #-------- Init
  ANS <- readline( "This function is ONLY requried for the first time of InDePTH. \nIt may use large memory and a lot of time (maybe more than a day). \nDo you run this function? [Y/N]")
  if( ANS == "Y" ){
    cat( "Running init function..." )
  }else if( ANS == "N" ){
    cat("Execution was stopped")
    return( invisible() )
  }else{
    cat( "Please enter Y/N." )
    return( invisible() )
  }

  #----- Load packages
  cat("Loading Required packages ...\n")
  require( readr )
  require( dplyr )
  require( rhdf5 )
  require( xlsx )
  require( igraph )

  #----- check files
  cat(sprintf( "Loading required files from %s ...\n", input.dir))
  if( any( c( "inst.info",
              "l000-chem-pert-20413 (1).xlsx",
              #"l000-genetic-perturbagens.xlsx",
              "zspc_n1328098x22268.gctx",
              "q2norm_n1328098x22268.gctx",
              "HG-U133_Plus_2.na34.annot.csv" ) %in% list.files( input.dir ) == FALSE ) ){
    cat( "Following files do NOT existed:")
    c( "inst.info",
       "l000-chem-pert-20413 (1).xlsx",
       #"l000-genetic-perturbagens.xlsx",
       "zspc_n1328098x22268.gctx",
       "q2norm_n1328098x22268.gctx",
       "HG-U133_Plus_2.na34.annot.csv" )[ !( c( "inst.info",
                                          "l000-chem-pert-20413 (1).xlsx",
                                          #"l000-genetic-perturbagens.xlsx",
                                          "zspc_n1328098x22268.gctx",
                                          "q2norm_n1328098x22268.gctx",
                                          "HG-U133_Plus_2.na34.annot.csv" ) %in% list.files( input.dir ) ) ] %>% print
    return( invisible() )
  }

  #----- Process annotation files
  cat( "Processing HG-U133_Plus_2.na34.annot.csv annotation files ...\n")
  process.anno <- function(NUMBER){
    NAME <- sprintf("HG-U133_Plus_2.na%s.annot.csv", NUMBER)
    anno <- read.table(file = file.path(input.dir, NAME),
                       sep = ",",
                       skip = 25,
                       header = TRUE)
    anno <- anno[, c("Probe.Set.ID", "Gene.Symbol")]
    colnames(anno) <- c("Probe.ID", "NAME")

    # collapse multiple Gene Symbol
    nrow.anno <- nrow(anno)
    anno2 <- data.frame(matrix(nrow = 0, ncol = 2))
    colnames(anno2) <- c("Probe.ID", "NAME")
    for(I in 1:nrow.anno){
      TMP.PROBE <- as.vector(anno[I, "Probe.ID"])
      TMP.NAME <- as.vector(anno[I, "NAME"])
      TMP.NAME.S <- strsplit(TMP.NAME, " /// ")[[1]]
      LEN.TMP.NAME.S <- length(TMP.NAME.S)
      tmp.anno <- data.frame(Probe.ID = rep(TMP.PROBE, LEN.TMP.NAME.S), NAME = TMP.NAME.S)
      anno2 <- rbind(anno2, tmp.anno)
    }

    anno <- anno2
    save(anno,
         file = file.path(out.dir, sprintf( "HG-U133_Plus_2.na%s.annot_for_mk_gmt.RData", NUMBER) ) )

  }
  process.anno(34)

  #---------- inst.info with HP information
  cat( "Processing inst.info and HP information files ...\n")
  inst.info <- read_delim( file = file.path(input.dir, "inst.info"),
                           delim = "\t",
                           col_types = "ccccdcdccccddcccccccccccccccccccc",
                           col_names = TRUE )
  l000.chem.pert.anno <- read.xlsx2( file.path( input.dir, "l000-chem-pert-20413 (1).xlsx" ),
                                     sheetIndex = 1)
  # l000.genetic.pert.anno <- read.xlsx2( file = file.path(input.dir, "l000-genetic-perturbagens.xlsx"),
  #                                       sheetIndex = 1)

  # combine HP annotation file
  hp.pert.info <- l000.chem.pert.anno %>%
    mutate( Reagent.Type = rep( "TRT_CP", nrow( l000.chem.pert.anno ) ) ) # %>%
    #bind_rows( l000.genetic.pert.anno )

  # join HP annotation file
  inst.info$pert_id <- as.vector( inst.info$pert_id )
  hp.pert.info$Reagent.ID <- as.vector( hp.pert.info$Reagent.ID )
  colnames( hp.pert.info )[ 1 ] <-c( "HP.Reagent.Name" )
  colnames( hp.pert.info )[ 3 ] <-c( "HP.Reagent.Type" )
  inst.info.2 <- inst.info %>%
    left_join( x = ., y = hp.pert.info,
               by = c( "pert_id" = "Reagent.ID" ) )

  # Basically, we believe inst.info
  inst.info.wt.hp <- inst.info.2 %>%
    mutate( pert_desc.md = ifelse( inst.info.2$pert_desc == "-666" &
                                     inst.info.2$pert_type == "trt_cp",
                                   as.vector( inst.info.2$HP.Reagent.Name ),
                                   as.vector( inst.info.2$pert_desc ) ) )

  # make RObject
  save( inst.info.wt.hp,
        file = file.path( out.dir, "inst.info.with.HPinfo.RData" ) )

  #---------- LINCS datasets
  cat( "Processing LINCS-providing datasets ...\n")

  #-------- Load ZSPC

  # landmark genes
  L1000.zspc.978 <- h5read(file = file.path(input.dir, "zspc_n1328098x22268.gctx"),
                           name = "/0/DATA/0/matrix",
                           start = c( 1,1 ),
                           count = c( 978, 1328098 ) )
  L1000.zspc.978 <- as.data.frame( L1000.zspc.978 )

  # get rownames
  row.978 <- h5read(file = file.path(input.dir, "zspc_n1328098x22268.gctx"),
                    name = "/0/META/ROW/id",
                    start = 1,
                    count = 978 )
  row.978 <- gsub(" ", "", row.978 )
  row.names( L1000.zspc.978 ) <- row.978

  # separate colnames
  TMP.COUNT <- c( 0, seq( 100000, 1000000, by = 100000 ), 1100000, 1200000, 1300000, 1328098 )
  for( I in 2:length( TMP.COUNT ) ){

    START = TMP.COUNT[ I - 1 ] + 1
    if( I %in% 2:( length( TMP.COUNT ) - 1 ) ){
      COUNT <- 100000
    }else{
      COUNT <- 28098
    }

    tmp.col.1328098 <- h5read(file = file.path(input.dir, "zspc_n1328098x22268.gctx"),
                              name = "/0/META/COL/id",
                              start = START,
                              count = COUNT)
    tmp.col.1328098 <- gsub(" ", "", tmp.col.1328098 )
    colnames( L1000.zspc.978 )[ START:( START + COUNT - 1 ) ] <- tmp.col.1328098
  }

  # make RObject
  save( L1000.zspc.978, file = file.path( out.dir, "L1000.zspc.978.RData"))

  #-------- Load Q2NORM

  # landmark genes
  L1000.q2norm.978 <- h5read(file = file.path(input.dir, "q2norm_n1328098x22268.gctx"),
                             name = "/0/DATA/0/matrix",
                             start = c( 1,1 ),
                             count = c( 978, 1328098 ) )
  L1000.q2norm.978 <- as.data.frame( L1000.q2norm.978 )

  # get rownames
  row.978 <- h5read(file = file.path(input.dir, "q2norm_n1328098x22268.gctx"),
                    name = "/0/META/ROW/id",
                    start = 1,
                    count = 978 )
  row.978 <- gsub(" ", "", row.978 )
  row.names( L1000.q2norm.978 ) <- row.978

  # get colnames
  TMP.COUNT <- c( 0, seq( 100000, 1000000, by = 100000 ), 1100000, 1200000, 1300000, 1328098 )
  for( I in 2:length( TMP.COUNT ) ){

    START = TMP.COUNT[ I - 1 ] + 1
    if( I %in% 2:( length( TMP.COUNT ) - 1 ) ){
      COUNT <- 100000
    }else{
      COUNT <- 28098
    }

    tmp.col.1328098 <- h5read(file = file.path(input.dir, "q2norm_n1328098x22268.gctx"),
                              name = "/0/META/COL/id",
                              start = START,
                              count = COUNT)# 1328098 )
    tmp.col.1328098 <- gsub(" ", "", tmp.col.1328098 )
    colnames( L1000.q2norm.978 )[ START:( START + COUNT - 1 ) ] <- tmp.col.1328098
  }

  # make RObject
  save( L1000.q2norm.978, file = file.path( out.dir, "L1000.q2norm.978.RData"))

  #-------- Make Rank Matrix
  cat( "Making rank matrix from LINCS datasets ...\n")

  # separate row_names
  row_names.df <- row.names( L1000.zspc.978 )

  # ranking function
  rank.z <-
    function( z, q ){
      # z : column of zspc
      # q : column of q2norm

      # z <- L1000.zspc.978[, 1]
      # q <- L1000.q2norm.978[, 1]

      # make matrix by z and q, (tie: min)
      set.seed( 123 )
      tmp.matrix <- cbind( z, q )
      tmp.matrix <- transform( tmp.matrix,
                               rank.z = rank( -z, ties.method = "min" ),
                               rank.q = rank( -q, ties.method = "random" ),
                               index = 1:length( z ) )

      # if rank.z was duplicated, consider rank.q
      dup.rank.z <- tmp.matrix[ duplicated( tmp.matrix$rank.z ), "rank.z" ]
      mstr <- tmp.matrix[ !( tmp.matrix$rank.z %in% dup.rank.z ),  ]

      # unique ranking
      for( i in unique( dup.rank.z ) ){
        tmp.tmp.matrix <- tmp.matrix[ tmp.matrix$rank.z == i, ]
        # ranking by rank.q ascending order
        tmp.tmp.matrix <- tmp.tmp.matrix[ order( tmp.tmp.matrix$rank.q, decreasing = FALSE ), ]

        # make unique rank
        tmp.tmp.matrix$rank.z <- tmp.tmp.matrix$rank.z + ( 0:( nrow( tmp.tmp.matrix ) - 1 ) )
        mstr <- rbind( mstr, tmp.tmp.matrix )
      }

      # sort by original index
      mstr <- mstr[ order( mstr$index ), ]
      return( mstr$rank.z )
    }

  # ranking
  rank.matrix.z <- matrix( 0, nrow = nrow( L1000.zspc.978 ), ncol = ncol( L1000.zspc.978 ))
  row.names( rank.matrix.z ) <- row.names( L1000.zspc.978 )
  colnames( rank.matrix.z ) <- colnames( L1000.zspc.978 )
  for( tmp.col in 1:ncol( L1000.zspc.978 ) ){
    rank.matrix.z[ ,tmp.col ] <- rank.z( z = L1000.zspc.978[   ,tmp.col ],
                                         q = L1000.q2norm.978[ ,tmp.col ] )
  }

  # write
  save( rank.matrix.z, file = file.path( out.dir, "RANK_MATRIX_sig978.RData") )

  #-------- Make DEG Matrix
  cat( "Making DEG matrix from LINCS datasets ...\n")

  deg.judge.fun <- function( x ){
    y <- rep( 0, length( x ) )
    y[ x >= 2  ] <- 1
    y[ x <= -2 ] <- -1
    return( y )
  }

  #- judgement
  row_names.df <- row.names( L1000.zspc.978 )
  DEG.L1000.zspc.978 <-
    L1000.zspc.978 %>%
    apply( ., 2, deg.judge.fun ) %>%
    as.data.frame() %>%
    mutate( row_names = row_names.df )

  # write
  save( DEG.L1000.zspc.978, file = file.path( out.dir, "DEG.L1000.zspc.978.RData") )

  #- 22268 genes, for save memories
  TMP.COUNT <- c( 0, seq( 10000, 1320000, by = 10000 ), 1328098 )
  meta.row.id <- h5read(file = file.path(input.dir,"zspc_n1328098x22268.gctx"),
                        name = "/0/META/ROW/id",
                        start = 1,
                        count = 22268)
  meta.row.id <- gsub(" ", "", meta.row.id)
  for( I in 2:length( TMP.COUNT ) ){

    START = TMP.COUNT[ I - 1 ] + 1
    if( I %in% 2:( length( TMP.COUNT ) - 1 ) ){
      COUNT <- 10000
    }else{
      COUNT <- 8098
    }

    # extract data
    tmp.L1000.zspc.22268 <- h5read(file = file.path(input.dir, "zspc_n1328098x22268.gctx"),
                               name = "/0/DATA/0/matrix",
                               start = c( 1, START ),
                               count = c( 22268, COUNT ) ) #1328098))
    tmp.L1000.zspc.22268 <- as.data.frame( tmp.L1000.zspc.22268 )

    # add row.names
    row.names( tmp.L1000.zspc.22268 ) <- meta.row.id

    # get colnames
    tmp.col.1328098 <- h5read(file = file.path(input.dir,"zspc_n1328098x22268.gctx"),
                              name = "/0/META/COL/id",
                              start = START,
                              count = COUNT)# 1328098 )
    tmp.col.1328098 <- gsub(" ", "", tmp.col.1328098 )
    colnames( tmp.L1000.zspc.22268 ) <- tmp.col.1328098

    #- deg judgement
    tmp.DEG.L1000.zspc.22268 <- apply( tmp.L1000.zspc.22268,
                                       2,
                                       deg.judge.fun )
    rm( tmp.L1000.zspc.22268 )
    tmp.DEG.L1000.zspc.22268 <- as.data.frame( tmp.DEG.L1000.zspc.22268 )
    tmp.DEG.L1000.zspc.22268 <- mutate( tmp.DEG.L1000.zspc.22268,
                                        row_names = meta.row.id )

    # rename
    eval( parse( text = sprintf( "DEG.L1000.zspc.22268.%s <- tmp.DEG.L1000.zspc.22268", I ) ) );rm(tmp.DEG.L1000.zspc.22268 )
    eval( parse( text = sprintf( "save( DEG.L1000.zspc.22268.%s, file = file.path( out.dir, 'DEG.L1000.zspc.22268.%s.RData') );rm( DEG.L1000.zspc.22268.%s );gc( T, T )",
                                 I, I, I )))
  }
  cat("Done.\n")
}
