#' From Lon/Lat to UTM zone
#'
#' Converts coordinates from longitude/latitude to UTM zone
#'
#' @param longitude A numeric vector of longitude coordinates in decimal degrees
#' @param latitude A numeric vector of latitude coordinates in decimal degrees
#' @return A numeric vector specifying the UTM zone the coordinates belong to.
LL2UTMZone <- function(longitude, latitude, add_latitude_band=FALSE) {
  
  if(add_latitude_band){
    
    #  band letters
    band <- LETTERS[3:24]
    band <- band[!band %in% c("I","O")] # letters "I" and "0" are skipped
    
    # horizontal bands spanning 8 degrees of latitude
    brks <- seq(from=-80, to=84, by=8)
    brks[length(brks)] <- 84 # band 'X' spans 12 degree
    
    latitude_band <- as.character(cut(latitude, breaks=brks, labels=band))
  }
  
  # Special zones for Svalbard and Norway
  if(latitude >= 72.0 && latitude < 84.0 )
    if (longitude >= 0.0  && longitude <  9.0)
      return(if(add_latitude_band) paste0(31,latitude_band) else 31);
  if (longitude >= 9.0  && longitude < 21.0)
    return(if(add_latitude_band) paste0(33,latitude_band) else 33)
  if (longitude >= 21.0 && longitude < 33.0)
    return(if(add_latitude_band) paste0(35,latitude_band) else 35)
  if (longitude >= 33.0 && longitude < 42.0)
    return(if(add_latitude_band) paste0(37,latitude_band) else 37)
  
  if(add_latitude_band) paste0( (floor((longitude + 180) / 6) %% 60) + 1, latitude_band) else (floor((longitude + 180) / 6) %% 60) + 1
}


#'              Most frequent value(s)
#'
#' Find the most frequent n value(s) from a vector
#'
#' @param x A vector
#' @param n A numeric integer specifying the number of values to return. Default is 1, i.e. the most frequent value.
#' @param ... Additional parameters to be passed to the function \code{table}.
#' @seealso \code{table}
modal <- function(x, n=1, ...){
  sort(table(x,...),decreasing=TRUE)[1:n]
}

#' utility function: find the minimum position for each row of a matrix, breaking ties at random.
min.col <- function(m, ...) max.col(-m,...)


#'                                 Kmeans clustering with near equal group size
#'
#' Modified from https://rviews.rstudio.com/2019/06/13/equal-size-kmeans/
#'
#'
#' @param kdat A two-column matrix or data.frame or data.table with longitude and latitude coordinates.
#' @param k A numeric integer specifying the number of cluster required.
#' @param iter_max A numeric integer specifying the maximum number of iterations to run.
#' @return A list with two elements:
#'         \code{Data} A two-column matrix or data.frame or data.table with longitude and latitude coordinates and a column 'assigned' specifying
#'         which of the k cluster each point belongs to.
#'         \code{Centers} A two-column data.frame with longitude and latitude coordinates of the cluster centers.
kmeansEqual <- function(kdat, k = 3, iter_max = 30, tolerance.iter=1e-08, random.seed=1234, method=c("iqr.outliers","top.outliers","utmzone.outliers"), nstart=20, n.outliers=10, verbose=TRUE, plot=TRUE){
  
  if(verbose){
    cat('#---------------------------------------------------------------------------\n')
    cat('#=step 0: Check input data\n')
    cat('#---------------------------------------------------------------------------\n')
  }
  
  if(!any(inherits(kdat, c("matrix","data.frame","data.table"))))
    stop("Argument kdat must be a matrix, a data.frame or a data.table")
  if(nrow(kdat) < 2)
    stop("Input data must have at least 2 coordinates.")
  
  if(!is.numeric(k))
    stop("Argument 'k' must be numeric.")
  if(k<0)
    stop("The number of cluster must be > 0.")
  if(k<2)
    warning("The number of group requested is < 2.")
  
  if(!is.numeric(iter_max))
    stop("Argument 'iter_max' must be numeric")
  if(iter_max<0)
    stop("The number of iterations must be > 0.")
  
  if(verbose){
    cat('#---------------------------------------------------------------------------\n')
    cat('#=step 1: Initialize "naive" clusters with a basic kmeans\n')
    cat('#---------------------------------------------------------------------------\n')
  }
  
  set.seed(random.seed)
  
  kclust <- kdat %>%
    stats::kmeans(k, nstart = nstart)
  
  # initial clusters
  centers = kclust$centers
  
  # transform to cartesian coordinates
  centers %<>%
    as.data.frame(row.names=NULL) %>%
    sf::st_as_sf(coords=1:2, crs=sf::st_crs(4326)) %>%
    sf::st_transform(crs=sf::st_crs("+proj=eqearth")) %>%
    sf::st_coordinates()
  
  if(is.matrix(kdat))
    kdat %<>%
    as.data.frame(row.names=NULL)
  
  # retrieve coordinates if needed
  if(!all(c("lon","lat") %in% colnames(kdat))){
    id_x_lon 	<- grep(pattern = "[Ll][Oo][Nn]|[Xx]",x = names(kdat))[1]
    id_y_lat 	<- grep(pattern = "[Ll][Aa][Tt]|[Yy]",x = names(kdat))[1]
    coordHeaders <- c(id_x_lon, id_y_lat)
    kdat %<>%
      dplyr::rename(lon = coordHeaders[1], lat = coordHeaders[2])
  }
  
  if(verbose){
    cat('#---------------------------------------------------------------------------\n')
    cat('#=step 2: Re-assign point to clusters\n')
    cat('#---------------------------------------------------------------------------\n')
  }
  
  # start iterations
  iter = 0
  process = TRUE
  while(iter <= iter_max & process){
    
    if(verbose) cat('iteration: ',iter+1,'\n')
    
    kclust <- kdat %>%
      dplyr::select(lon,lat) %>%
      sf::st_as_sf(coords=1:2, crs=sf::st_crs(4326)) %>%
      sf::st_transform(crs=sf::st_crs("+proj=eqearth")) %>%
      sf::st_coordinates() %>%
      stats::kmeans(centers)
    
    kdat$assigned = kclust$cluster
    #---------------------------------------------------------------------------
    #=step 2: Compute the distance matrix between each point and the centroids
    #---------------------------------------------------------------------------
    
    # compute mean distance between centers
    distMean <- mean(stats::dist(centers))
    
    # transform back to lon/lat
    centers %<>%
      as.data.frame(row.names=NULL) %>%
      sf::st_as_sf(coords=1:2, crs=sf::st_crs("+proj=eqearth")) %>%
      sf::st_transform(crs=sf::st_crs(4326)) %>%
      sf::st_coordinates()
    
    # calculate the distance to the center
    for(j in 1:nrow(centers)){
      kdat <- kdat %<>%
        dplyr::mutate(kdist = geosphere::distGeo(cbind(lon,lat), centers[j,]))
      # rename the column
      kdat <- kdat %>%
        dplyr::rename(!!paste0('kdist',j) := kdist)
    }
    if(iter > 0){
      # assign the outliers from each group to the closest cluster center
      switch(method[1],
             
             "top.outliers" = {
               j = 1
               kdat.copy <- kdat %>%
                 dplyr::filter(assigned == j) %>%
                 dplyr::arrange_at(paste0("kdist", j), dplyr::desc) %>%
                 dplyr::mutate(is_kdist_min =  dplyr::select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
                 dplyr::mutate(assigned = ifelse(dplyr::row_number() <= n.outliers, is_kdist_min, assigned))
               for(j in 2:k){
                 kdat.copy %<>%
                   rbind(kdat %>%
                           dplyr::filter(assigned == j) %>%
                           dplyr::arrange_at(paste0("kdist", j), dplyr::desc) %>%
                           dplyr::mutate(is_kdist_min =  dplyr::select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
                           dplyr::mutate(assigned = ifelse(dplyr::row_number() <= n.outliers, is_kdist_min, assigned))
                   )
               }
             },
             
             "iqr.outliers" = {
               j = 1
               kdat.copy <- kdat %>%
                 dplyr::filter(assigned == j) %>%
                 dplyr::mutate_at(.vars= dplyr::vars(!!paste0("kdist", j)),
                                  .funs= list(is.outliers = ~ . > stats::quantile(., probs=.75) + 1.5 * stats::IQR(.))) %>%
                 dplyr::mutate(is_kdist_min =  dplyr::select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
                 dplyr::mutate(assigned = ifelse(is.outliers, is_kdist_min, assigned))
               for(j in 2:k){
                 kdat.copy %<>%
                   rbind(kdat %>%
                           dplyr::filter(assigned == j) %>%
                           dplyr::mutate_at(.vars= dplyr::vars(!!paste0("kdist", j)),
                                            .funs= list(is.outliers = ~ . > stats::quantile(., probs=.75) + 1.5 * stats::IQR(.))) %>%
                           dplyr::mutate(is_kdist_min =  dplyr::select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
                           dplyr::mutate(assigned = ifelse(is.outliers, is_kdist_min, assigned))
                   )
               }
             },
             
             "utmzone.outliers" = {
               j = 1
               kdat.copy <- kdat %>%
                 dplyr::filter(assigned == j) %>%
                 dplyr::mutate(is.outliers = LL2UTMZone(lon,lat,add_latitude_band=TRUE) != names(modal(LL2UTMZone(lon,lat,add_latitude_band=TRUE)))) %>%
                 dplyr::mutate(is_kdist_min =  dplyr::select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
                 dplyr::mutate(assigned = ifelse(is.outliers, is_kdist_min, assigned))
               for(j in 2:k){
                 kdat.copy %<>%
                   rbind(kdat %>%
                           dplyr::filter(assigned == j) %>%
                           dplyr::mutate(is.outliers = LL2UTMZone(lon,lat,add_latitude_band=TRUE) != names(modal(LL2UTMZone(lon,lat,add_latitude_band=TRUE)))) %>%
                           dplyr::mutate(is_kdist_min =  dplyr::select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
                           dplyr::mutate(assigned = ifelse(is.outliers, is_kdist_min, assigned))
                   )
               }
             }
      )
      kdat.copy %<>% data.frame(row.names=NULL)
      kdat <- kdat.copy
      
    }
    #---------------------------------------------------------------------------
    #=step 3: Re-assign points to clusters
    #---------------------------------------------------------------------------
    kdat$index = 1:nrow(kdat)
    working = kdat
    ReAssignment = nrow(kdat) - (nrow(kdat) %% k)
    
    for(i in 1:ReAssignment){
      #cluster counts can be off by 1 due to uneven multiples of k.
      j = if(i %% k == 0) k else (i %% k)
      itemloc =  working$index[which.min(working[,(paste0("kdist", j))])][1]
      kdat$assigned[kdat$index == itemloc] = j
      working %<>%
        dplyr::filter(!index == itemloc)
    }
    # if there are some leftover points, assign to whoever's closest, without regard to k
    if(length(working$index)>0L){
      for(i in working$index){
        kdat$assigned[kdat$index == i] = which.min(working[working$index == i,grepl("kdist",colnames(working))])
      }
    }
    
    # assign the outliers from each group to the closest cluster center
    switch(method[1],
           
           "top.outliers" = {
             j = 1
             kdat.copy <- kdat %>%
               dplyr::filter(assigned == j) %>%
               dplyr::arrange_at(paste0("kdist", j), dplyr::desc) %>%
               dplyr::mutate(is_kdist_min =  dplyr::select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
               dplyr::mutate(assigned = ifelse(dplyr::row_number() <= n.outliers, is_kdist_min, assigned))
             for(j in 2:k){
               kdat.copy %<>%
                 rbind(kdat %>%
                         dplyr::filter(assigned == j) %>%
                         dplyr::arrange_at(paste0("kdist", j), dplyr::desc) %>%
                         dplyr::mutate(is_kdist_min =  dplyr::select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
                         dplyr::mutate(assigned = ifelse(dplyr::row_number() <= n.outliers, is_kdist_min, assigned))
                 )
             }
           },
           
           "iqr.outliers" = {
             j = 1
             kdat.copy <- kdat %>%
               dplyr::filter(assigned == j) %>%
               dplyr::mutate_at(.vars= dplyr::vars(!!paste0("kdist", j)),
                                .funs= list(is.outliers = ~ . > stats::quantile(., probs=.75) + 1.5 * stats::IQR(.))) %>%
               dplyr::mutate(is_kdist_min =  dplyr::select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
               dplyr::mutate(assigned = ifelse(is.outliers, is_kdist_min, assigned))
             for(j in 2:k){
               kdat.copy %<>%
                 rbind(kdat %>%
                         dplyr::filter(assigned == j) %>%
                         dplyr::mutate_at(.vars= dplyr::vars(!!paste0("kdist", j)),
                                          .funs= list(is.outliers = ~ . > stats::quantile(., probs=.75) + 1.5 * stats::IQR(.))) %>%
                         dplyr::mutate(is_kdist_min =  dplyr::select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
                         dplyr::mutate(assigned = ifelse(is.outliers, is_kdist_min, assigned))
                 )
             }
           },
           
           "utmzone.outliers" = {
             j = 1
             kdat.copy <- kdat %>%
               dplyr::filter(assigned == j) %>%
               dplyr::mutate(is.outliers = LL2UTMZone(lon,lat,add_latitude_band=TRUE) != names(modal(LL2UTMZone(lon,lat,add_latitude_band=TRUE)))) %>%
               dplyr::mutate(is_kdist_min =  dplyr::select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
               dplyr::mutate(assigned = ifelse(is.outliers, is_kdist_min, assigned))
             for(j in 2:k){
               kdat.copy %<>%
                 rbind(kdat %>%
                         dplyr::filter(assigned == j) %>%
                         dplyr::mutate(is.outliers = LL2UTMZone(lon,lat,add_latitude_band=TRUE) != names(modal(LL2UTMZone(lon,lat,add_latitude_band=TRUE)))) %>%
                         dplyr::mutate(is_kdist_min =  dplyr::select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
                         dplyr::mutate(assigned = ifelse(is.outliers, is_kdist_min, assigned))
                 )
             }
           }
    )
    kdat.copy %<>% data.frame(row.names=NULL)
    kdat <- kdat.copy
    if(verbose){
      cat('#---------------------------------------------------------------------------\n')
      cat('#=step 4: Recalculate the centroids\n')
      cat('#---------------------------------------------------------------------------\n')
    }
    j = 1
    NewCenters <- kdat %>%
      dplyr::filter(assigned == j) %>%
      dplyr::select(lon, lat) %>%
      sf::st_as_sf(coords=1:2, crs=sf::st_crs(4326)) %>%
      sf::st_transform(crs=sf::st_crs("+proj=eqearth")) %>%
      sf::st_coordinates() %>%
      stats::kmeans(1) %$% centers
    while(j < k){
      j = j+1
      NewCenters %<>%
        rbind(kdat %>%
                dplyr::filter(assigned == j) %>%
                dplyr::select(lon, lat) %>%
                sf::st_as_sf(coords=1:2, crs=sf::st_crs(4326)) %>%
                sf::st_transform(crs=sf::st_crs("+proj=eqearth")) %>%
                sf::st_coordinates() %>%
                stats::kmeans(1) %$% centers)
    }
    # New centroids
    NewCenters %<>% data.frame(row.names=NULL)
    centers <- NewCenters
    
    # keep assigments lon, lat and assigments only
    kdat <- kdat %>%
      dplyr::select(lon, lat, assigned)
    
    if(plot){
      cent <- centers %>%
        as.data.frame(row.names=NULL) %>%
        sf::st_as_sf(coords=1:2, crs=sf::st_crs("+proj=eqearth")) %>%
        sf::st_transform(crs=sf::st_crs(4326)) %>%
        sf::st_coordinates() %>%
        as.data.frame(row.names=NULL)
      
      if(!is.factor(kdat$assigned))
        kdat$assigned %<>% as.factor()
      print(
        kdat %>% ggplot2::ggplot(ggplot2::aes(x = lon, y = lat, color = assigned)) +
          ggplot2::theme_minimal() + ggplot2::geom_point()  +
          ggplot2::geom_point(data = cent, ggplot2::aes(x = X, y = Y), color = "black", size = 4)
      )
    }
    
    process = abs(mean(stats::dist(centers)) - distMean) > tolerance.iter
    iter = iter + 1
  }
  # end iterations
  
  # assignment to factor
  if(!is.factor(kdat$assigned))
    kdat$assigned %<>% as.factor()
  
  centers %<>%
    sf::st_as_sf(coords=1:2, crs=sf::st_crs("+proj=eqearth")) %>%
    sf::st_transform(crs=sf::st_crs(4326)) %>%
    sf::st_coordinates()
  
  return(list(Data=kdat, Centers=centers, converged=!process))
}

