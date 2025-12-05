#'                                 Estimating the species distribution range with alpha-hulls
#'
#' Make an alpha-hull from the spatial distribution of species occurrence records
#'
#' @param loc_data A two-column matrix or data.frame or data.table or a character string specifying the path to your occurrence files directory (i.e. where species occurrence data files are stored) OR path to your species occurrence file (must be a csv file)
#' @param output_dir An optional character string specifying the path to the directory where to save the shapefiles generated.
#' @param binom_col A numeric value indicating which column number contains the binomials for the species (i.e. Genus species).Ignored if loc_data is not a csv file.
#' @param fraction A numeric between 0 and 1 specifying the minimum fraction of occurrences that must be included in the polygon.
#' @param partCount A numeric integer specifying the maximum number of disjunct polygons that are allowed.
#' @param initialAlpha A numeric specifying the starting value for the parameter alpha.
#' @param coordHeaders A character string vector of length 2 specifying the names of the columns indicating the longitude and latitude respectively.
#' @param clipToCoast A logical. Should the terrestrial \code{clipToCoast='terrestrial'}, the aquatic part \code{clipToCoast='aquatic'} or none part \code{clipToCoast='no'} part of the range be kept ?
#' @param land_file  An optional sf object or character string specifying a path to a land basemap (shapefile or rds file) for clipping species polygon ranges by lands.
#' @param proj A character string specifying the projection of the coordinates
#' @param alphaIncrement A numeric specifying the amount to increase alpha at each iteration
#' @param alphaDecrement A numeric specifying the amount to decrease alpha if the function fails after the first iterations
#' @param maxIter A numeric integer specifying the maximum number of iterations before trying to build circular buffers.
#' @param default_buffer A numeric specifying the buffering distance in meters around each point in case the buffer method equals 0. Default is 200 km.
#' @param other_buffers A numeric value or a vector of values of length \code{maxIter} specifying alternative buffer sizes (in meters) if the function fails after the first iterations.
#' @return None
#' @export
make_alpha_hulls <- function(loc_data, output_dir=NULL,
                             binom_col = NULL,
                             coordHeaders = NULL,
                             fraction=0.95,
                             partCount=10,
                             initialAlpha=2,
                             alphaIncrement=1,
                             alphaDecrement=1,
                             maxIter=2,
                             default_buffer = 200000,
                             other_buffers=c(200000),
                             tolerance=5,
                             clipToCoast="terrestrial",
                             land_file = NULL,
                             proj='+proj=longlat +datum=WGS84',
                             save.outputs = FALSE,
                             do.parallel = FALSE,
                             ncores = 1,
                             verbose = TRUE) {
  if(verbose){
    message('#=================')
    message('#= 1. Check inputs')
    message('#=================\n')
  }
  if(missing(loc_data))
    stop("Please a csv file or a data.frame or a matrix with longitude and latitude coordinates of species occurrence records.")

  file_flag = tryCatch(file.exists(loc_data), error=function(err) FALSE) && !tryCatch(dir.exists(loc_dat), error=function(err) FALSE)
  dir_flag  = tryCatch(dir.exists(loc_data), error=function(err) FALSE) && !file_flag
  data_flag = !file_flag & any(inherits(loc_data, c("data.frame","matrix")))

  if(!file_flag & !dir_flag & !data_flag)
    stop('Unable to read input data. Please provide valid input data')

  if(file_flag && length(readLines(loc_data))==0L)
    stop(paste("The file provided:",loc_data," is empty!"))

  # ensure that the output directory provided exists
  if(!is.null(output_dir) && !dir.exists(output_dir)){
    warning(paste("Output directory:", output_dir,"does not exist or could not be found."));flush.console()
    if(save.outputs)
      stop("Unable to save the outputs if an output directory is not specified.")
  }
  if(is.null(output_dir) & save.outputs){
    stop("Unable to save the outputs since the output directory is not specified.")
  }

  coastline = NULL
  if(!is.null(land_file)){
    coastline = suppressMessages(
      if(tryCatch(file.exists(land_file) && !dir.exists(land_file),error=function(err) FALSE))  switch(stringr::str_extract(basename(land_file),'(?>.{3}$)'), 'shp' = sf::st_read(land_file), 'rds' = readRDS(land_file)) else if(inherits(land_file,"sf")) land_file else NULL
    )
  }
  if(is.null(coastline)){
    # The default is world land map at 50m resolution prodived by @naturalearth via the package 'rangeBuilder')
    # Tip: high land basemap resolution may help to avoid crashes of the script
    coastline = rangeBuilder:::loadWorldMap()
  }

  if(do.parallel){
    if(!is.null(ncores) & !is.numeric(ncores))
      stop("Argument 'ncores' must be a numeric integer")

  }

  if(verbose){
    message('#===================')
    message('#= 2. Set parameters')
    message('#===================\n')
  }


  if(verbose){
    message('#--------------------------------')
    message('#= 2.a Make a list of all species')
    message('#--------------------------------')
  }
  if(dir_flag){
    list.species <- list.files(loc_data, pattern="\\.csv$", full.names=TRUE, recursive=TRUE) # select species
  }else if(file_flag){
    occ_data 	<- read.csv(loc_data, fill=TRUE)
    if(!is.null(binom_col) & !is.numeric(binom_col))
      stop("Argument 'binom_col' is not set or must be an integer.")
    if(is.null(binom_col)){
      warning("Argument 'binom_col' null. Data will be treated as one single unknown species")
      list.species = "Unknown"
      occ_data <- as.data.frame(cbind(Species=rep(list.species,nrow(occ_data)),occ_data))
      binom_col=1
    }
    if(binom_col<1 || binom_col>ncol(occ_data))
      stop("Argument 'binom_col' is invalid")
    is.duplicated <- duplicated(occ_data)
    if(sum(is.duplicated)>0L)
      warning("The dataset contains duplicated records!")
    list.species 	<- c(as.character(unique(occ_data[, binom_col]))) # select species
  }else if(data_flag){
    occ_data <- loc_data
    if(!is.null(binom_col) & !is.numeric(binom_col))
      stop("Argument 'binom_col' is not set or must be an integer.")
    if(is.null(binom_col)){
      warning("Argument 'binom_col' null. Data will be treated as one single unknown species")
      list.species = "Unknown"
      occ_data <- as.data.frame(cbind(Species=rep(list.species,nrow(occ_data)),occ_data))
      binom_col=1
    }
    if(binom_col<1 || binom_col>ncol(occ_data))
      stop("Argument 'binom_col' is invalid")
    is.duplicated <- duplicated(occ_data)
    if(sum(is.duplicated)>0L)
      warning("The dataset contains duplicated records!")
    list.species 	<- c(as.character(unique(occ_data[, binom_col]))) # select species
  }
  else{
    stop(paste0("Unable to find directory or file",loc_data,"Please provide a valid path"))
  }

  if(verbose){
    message('#-----------------------------')
    message('#= 2.b Check parameters values')
    message('#-----------------------------')
  }
  if(fraction>1 | fraction<=0)
    stop("Parameter 'fraction' must be between > 0 and <= 1")
  if(partCount<=0)
    stop("Parameter 'partCount' must be > 0")
  if(initialAlpha<=0)
    stop("Parameter 'initialAlpha' must be > 0")
  if(alphaIncrement<=0)
    stop("Parameter 'alphaIncrement' must be > 0")

  if(verbose & do.parallel){
    message('#============================')
    message('#= 3. Set parallel processing')
    message('#============================\n')
  }
  if(do.parallel & length(list.species)>1){
    toExport <- c("fraction","partCount","initialAlpha","clipToCoast","alphaIncrement","alphaDecrement","default_buffer","save.outputs")# send objects to cluster nodes
    if(exists('occ_data', envir = environment()))
      toExport <- append(toExport,c("occ_data", "binom_col"))
    if(!is.null(coastline)){
      if(inherits(coastline,c("sf","SpatialPolygonsDataFrame","SpatialPolygons")))
        toExport <- append(toExport,"coastline")
    }

    if(is.null(ncores))
      ncores = min(6,parallel::detectCores()-1) # number of cores to use
    else
      ncores = min(ncores,parallel::detectCores())
    doParallel::registerDoParallel(ncores)
    `%execute%` <- foreach::`%dopar%`
  }else{
    toExport <- NULL
    # Sequential processing
    foreach::registerDoSEQ()
    `%execute%` <- foreach::`%do%`
  }

  if(verbose){
    cat(paste0('#',paste(rep('-',times=100),collapse="")))
    cat('\n')
    cat(paste0('> Run started on: ',Sys.time()));
    cat('\n')
  }
  # Time starts
  started.at <- Sys.time()
  list.species.polygons <- foreach::foreach(k = list.species,
                                            .packages = c("AlphaHullRangeModeller","sf"),
                                            .export=toExport) %execute% {

    out = tryCatch({

      # Get data for just that species
      if(file.exists(k)){
        # reads as data.table
        species.data <- tryCatch(data.table::fread(k) ,	error = function(err) return(NULL))
        k = gsub("\\.csv$","",basename(k))
      }else{
        species.data 	<- tryCatch({
          subset(occ_data, as.character(occ_data[, binom_col]) == k)
        },
        error = function(err) {
          return(NULL)
        })
      }
      # Print species name
      if(verbose) print(sprintf("Species name: %s", k))
      flush.console()

      # if an error occurred during the loading/reading returns an empty polygon list
      if (is.null(species.data) || ncol(species.data) < 2) {
        species.polygons <- list(NA)
        # Assign species name to the list
        names(species.polygons)[[1]] <- k
        if(save.outputs){
          if(!dir.exists(file.path(output_dir,"problems"))) dir.create(file.path(output_dir,"problems"), recursive =TRUE)
          write.table(data.frame(Species=k), file=file.path(output_dir,"problems","Cannot_read_file.csv"), row.names=FALSE,
                      col.names = !file.exists(file.path(output_dir,"problems","Cannot_read_file.csv")), sep=",",  append=TRUE)
        }
        return(species.polygons)
      }
      # if the species has no points
      if(nrow(species.data)< 1){
        species.polygons <- list(NA)
        # Assign species name to the list
        names(species.polygons)[[1]] <- k
        if(save.outputs){
          if(!dir.exists(file.path(output_dir,"problems"))) dir.create(file.path(output_dir,"problems"), recursive =TRUE)
          write.table(data.frame(Species=k), file=file.path(output_dir,"problems","Has_no_points.csv"),row.names=FALSE,
                      col.names = !file.exists(file.path(output_dir,"problems","Has_no_points.csv")), sep=",",  append=TRUE)
        }
        return(species.polygons)
      }

      # convert to data.frame
      if(inherits(species.data,"matrix"))
        species.data <- as.data.frame(species.data)
      else
        data.table::setDF(species.data)
      id_x_lon 	<- grep(pattern = "[Ll][Oo][Nn]|^[xX]$",x = names(species.data))[1]
      id_y_lat 	<- grep(pattern = "[Ll][Aa][Tt]|^[yY]$",x = names(species.data))[1]
      colnames(species.data)[id_x_lon] = "long"
      colnames(species.data)[id_y_lat] = "lat"

      # make sure coordinates are numeric
      if(any(!sapply(species.data[, c("long", "lat")], is.numeric))){
        species.data[, c("long", "lat")][, !sapply(species.data[, c("long", "lat")], is.numeric)] <-  
          sapply(species.data[, c("long", "lat")][, !sapply(species.data[, c("long", "lat")], is.numeric)], 
                 function(f)if(is.factor(f)) as.numeric(levels(f))[f] else as.numeric(f))        
      }

      # clean data
      NA_to_remove <- !complete.cases(species.data[, c("long", "lat")]) # find NA's
      Dup_to_remove <- duplicated(round(species.data[, c("long", "lat")], tolerance)) # find duplicates
      Rows_to_remove <- NA_to_remove | Dup_to_remove
      species.data.cleaned	<- species.data[!Rows_to_remove, ] # remove Na's and duplicates

      # compute distance to meridians +/- 180°
      distTo180	<- geosphere::distGeo(p1 = species.data.cleaned[, c("long", "lat")],
                                      p2 = cbind(long = 180, lat = species.data.cleaned[, "lat"]))

      # distTo180 <- sf::st_distance(sf::st_as_sf(species.data.cleaned[, c("long", "lat")], coords=1:2, crs=4326) %>%
      #                                sf::st_transform("+proj=eqearth +wktext"),
      #                               sf::st_as_sf(data.frame(long = 180, lat = species.data.cleaned[, "lat"]), coords=1:2, crs=4326) %>%
      #                                sf::st_transform("+proj=eqearth +wktext"), by_element = TRUE) %>%
      #   units::drop_units()
      
      # compute the "1/10th max" buffer
      one_tenth 	<- get_OneTenth_distmax(species.data.cleaned, default_buffer)
      id.edges	<- which(distTo180 < one_tenth)

      # Get alpha hull polygon(s)
      # If we're close to 180 then let's do this as west, centre, east
      if (length(id.edges) > 0L) {
        species.polygons.west <- species.polygons.east <- NA

        #central
        species.data.central <- species.data.cleaned[-id.edges, ]
        onetenth = get_OneTenth_distmax(species.data.central,default_buffer)
        species.polygons.central = suppressMessages(PolygonMaker(species.data.central,
                                                fraction = fraction,
                                                partCount = partCount,
                                                buffer = onetenth,
                                                initialAlpha = initialAlpha,
                                                alphaIncrement = alphaIncrement,
                                                alphaDecrement = alphaDecrement,
                                                clipToCoast = clipToCoast,
                                                coastline = coastline,
                                                proj=proj))

        # West
        if (nrow(subset(species.data.cleaned[id.edges, ], long <0)) > 0L) {
          species.data.west <- subset(species.data.cleaned[id.edges, ], long < 0)
          onetenth = get_OneTenth_distmax(species.data.west,default_buffer)

          west.buff <- floor(min(distTo180[which(distTo180 < one_tenth & species.data.cleaned$long < 0)])) / 2.

          species.polygons.west=suppressMessages(PolygonMaker(species.data.west,
                                             fraction = 1,
                                             partCount = partCount,
                                             buffer = onetenth,
                                             initialAlpha = initialAlpha,
                                             alphaIncrement = alphaIncrement,
                                             alphaDecrement = alphaDecrement,
                                             clipToCoast = clipToCoast,
                                             coastline = coastline,
                                             proj=proj,
                                             other_buffers = west.buff))
        }

        # East
        if (nrow(subset(species.data.cleaned[id.edges, ], long >0)) > 0L) {
          species.data.east <- subset(species.data.cleaned[id.edges, ], long > 0)
          onetenth = get_OneTenth_distmax(species.data.east,default_buffer)
          east.buff <- floor(min(distTo180[which(distTo180 < one_tenth & species.data.cleaned$long > 0)])) / 2.
          species.polygons.east=suppressMessages(PolygonMaker(species.data.east,
                                             fraction = 1,
                                             partCount = partCount,
                                             buffer = onetenth,
                                             initialAlpha = initialAlpha,
                                             alphaIncrement = alphaIncrement,
                                             alphaDecrement = alphaDecrement,
                                             clipToCoast = clipToCoast,
                                             coastline = coastline,
                                             proj=proj,
                                             other_buffers = east.buff))
        }

        #bind polygons
        edged.polygons	<- list(
          if (inherits(species.polygons.east,"list"))
            unlist(species.polygons.east)
          else
            species.polygons.east,
          if (inherits(species.polygons.west,"list"))
            unlist(species.polygons.west)
          else
            species.polygons.west
        )

        centered.polygons 	<-  if (inherits(species.polygons.central,"list"))  species.polygons.central else list(species.polygons.central)
        all.polygons 		<-  if (length(edged.polygons[!is.na(edged.polygons)]) > 0L) c(centered.polygons, edged.polygons[!is.na(edged.polygons)]) else centered.polygons
        species.polygons	<- 	if (length(all.polygons[!is.na(all.polygons)]) > 0L) Reduce("c",all.polygons[!is.na(all.polygons)])	else list()

        #put them in a list to add points
        if (length(species.polygons) == 0L){
          species.polygons <- list(NA)
        }else{
          # Assign species name to the polygon(s)
          species.polygons <- list(
            sf::st_sf(data.frame(species=rep(k, length(species.polygons)),
                                                   geom=species.polygons)) %>%
              sf::st_union()
          )
          }
        
        }else{
        # Get polygons
        species.polygons	<- suppressMessages(PolygonMaker(species.data.cleaned,
                                         fraction = fraction,
                                         partCount = partCount,
                                         buffer = one_tenth,
                                         initialAlpha = initialAlpha,
                                         alphaIncrement = alphaIncrement,
                                         alphaDecrement  = alphaDecrement,
                                         other_buffers = one_tenth,
                                         clipToCoast = clipToCoast,
                                         coastline=coastline,
                                         proj=proj))
        species.polygons <- list(
          sf::st_sf(data.frame(species=rep(k, length(species.polygons)),
                               geom=sf::st_geometry(species.polygons)))
        )
      }

      species.polygons$occ_points <-	species.data.cleaned[, c("long", "lat")]
      species.polygons$to_remove	<- 	which(Rows_to_remove)

      # Get potential error in coordinates
      # cooErrs <- rangeBuilder::coordError(species.data.cleaned[, c("long", "lat")],
      #                                     nthreads = 1)

      #save polygons and points
      if(save.outputs){
        dir.create(file.path(output_dir,"Outputs","polygons", names(species.polygons)[1]), recursive =TRUE)
        dir.create(file.path(output_dir,"Outputs","points", names(species.polygons)[1]), recursive = TRUE)

        if (!inherits(species.polygons[[1]],"sf")) return(NULL)
        sf::st_write(
            species.polygons[[1]],
            file.path(output_dir,
                      "Outputs",
                      "polygons", names(species.polygons)[1]),
            names(species.polygons)[1],
            driver = "ESRI Shapefile",
            delete_layer = TRUE,
            quiet=TRUE
        )
        sf::st_write(
          sf::st_as_sf(species.polygons$occ_points,
                       coords=c(1,2),
                       crs=sf::st_crs(species.polygons[[1]])),
          file.path(output_dir,
                    "Outputs",
                    "points", names(species.polygons)[1]),
          names(species.polygons)[1],
          driver = "ESRI Shapefile",
          delete_layer = TRUE,
          quiet=TRUE
        )
      }

      species.polygons[[1]]

    }, 
    error=function(err){
        return(err)
    },
    finally ={
      ## Stop the cluster
      doParallel::stopImplicitCluster()
    })

    if(inherits(out,"error")){
      message(out)
      return(NULL)
    }
    return(out)
  }

  # Time stops
  finished.at <- Sys.time()
  doParallel::stopImplicitCluster()
  time_taken <- finished.at - started.at # calculate amount of time elapsed
  if(verbose){
    message(paste0('> End of Run on: ',Sys.time()),'\n');
    message(paste0('> Running time: ',as.numeric(time_taken),' ', attr(time_taken,'units')))
    message(paste0('#',paste(rep('-',times=100),collapse="")))
    message('\n')
  }

  if(length(list.species.polygons)==1)
    return(list.species.polygons[[1]])

  return(list.species.polygons)
}



