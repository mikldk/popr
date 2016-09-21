#' #' Get distribution of number of meioses within a given radius
#' #' 
#' #' For each individual (pid) in pids, a search is performed for all individuals in the population within a certain radius. 
#' #' All the individuals within the radius are classified as either related or unrelated 
#' #' (depending on whether they belong to the same pedigree or not). If they do belong to the same pedigree,
#' #' the number of meioses are found.
#' #'  
#' #' Pedigree information is linked to the individuals in the population. All pedigrees must be trees (FIXME).
#' #' 
#' #' @param population A collection of individuals (\code{"popr_population"})
#' #' @param pids Vector of pids to get meioses_distributions for
#' #' @param radii Vector of radii (radiuses) to search within [km]
#' #' @param confirm_all_pedigrees_are_trees Set to TRUE to confirm that all pedigrees are tree (this is the caller's responsibility to ensure that as it would be too time consuming to check this requirement in this function)
#' #' 
#' #' @return A tibble with meioses distribution information
#' #' 
#' #' @export
#' meioses_distribution_for_tree_pedigrees <- function(population, pids, radii, confirm_all_pedigrees_are_trees = FALSE) {
#'   if (!is(population, "popr_population")) stop("population must be a popr_population object")
#'   
#'   if (length(confirm_all_pedigrees_are_trees) != 1L || !is.logical(confirm_all_pedigrees_are_trees) || confirm_all_pedigrees_are_trees == FALSE) {
#'     stop("You have not confirmed that all the pedigrees are trees (fulfilled e.g. for a population with only male lineages)")
#'   }
#'   
#'   radii <- sort(unique(radii), decreasing = TRUE)
#'   m <- analyse_meioses_internal_tree(population, pids, radii)
#'   
#'   if (nrow(m) == 1L && ncol(m) == 1L && is.na(m[1L, 1L])) {
#'     stop("Function was aborted") # FIXME: error class instead?
#'   }
#'   
#'   d <- tibble::tibble(pid = m[, 1L], Radius = radii[m[, 2L]], Distance = m[, 3L], count = m[, 4L])
#'   
#'   d$Type <- "Related"
#'   d$Type[d$Distance == -1L] <- "Unrelated"
#'   d$Type[d$Distance != -1L] <- "Related"
#'   d$Distance[d$Distance == -1L] <- NA
#'   
#'   return(d)
#' }
#' 

#' Get distribution of number of meioses within a given radius. 
#' The Euclidian distance is used based on ETRS89 coordinates due to speed. 
#' Hence, be sure that this is an acceptable error in comparison to e.g. great-circle distance).
#'
#' For each individual (pid) in pids, a search is performed for all individuals in the population within a certain radius.
#' All the individuals within the radius are classified as either related or unrelated
#' (depending on whether they belong to the same pedigree or not). If they do belong to the same pedigree,
#' the number of meioses are found.
#'
#' Pedigree information is linked to the individuals in the population. All pedigrees must be trees (FIXME).
#'
#' @param population A collection of individuals (\code{popr_population})
#' @param search_tree_unrelated A search kd-tree when counting unrelated individuals (typically based on \code{population}), created with e.g. \link{\code{build_kdtree_from_population}} or \link{\code{build_kdtree_from_pids}}
#' @param max_leaf_size_pedigree A parameter for building kd-trees for the pedigrees that he individuals with pid in \code{pids} belongs to
#' @param pids Vector of pids to get meioses_distributions for (must be individuals in \code{population})
#' @param radii Vector of radii (radiuses) to search within [km]
#' @param pedigree_must_be_alive specifies if individuals in pedigree must be alive to be included (typically, you want \code{search_tree_unrelated} to have same value), if FALSE it is not check if individuals are alive or not (hence, it is ignored)
#' @param pedigree_birth_year_range specifies valid birth year range, \code{c(xmin, xmax)}, such that only individuals in pedigree with \code{xmin <= birth_year <= xmax} are included (typically, you want \code{search_tree_unrelated} to have same value)
#' @param confirm_all_pedigrees_are_trees Set to TRUE to confirm that all pedigrees are tree (this is the caller's responsibility to ensure that as it would be too time consuming to check this requirement in this function)
#'
#' @return A tibble with meioses distribution information
#'
#' @export
meioses_distribution_for_tree_pedigrees <- function(population, search_tree_unrelated, max_leaf_size_pedigree = 10, pids, radii, pedigree_must_be_alive = FALSE, pedigree_birth_year_range, confirm_all_pedigrees_are_trees = FALSE) {
  if (!is(population, "popr_population")) stop("population must be a popr_population object")
  
  if (!is(search_tree_unrelated, "popr_search_data_structure")) stop("search_tree_unrelated must be a popr_search_data_structure object")
  
  if (length(confirm_all_pedigrees_are_trees) != 1L || !is.logical(confirm_all_pedigrees_are_trees) || confirm_all_pedigrees_are_trees == FALSE) {
    stop("You have not confirmed that all the pedigrees are trees (fulfilled e.g. for a population with only male lineages). Currently, only trees are supported. The plan is to change that in the future. Please notify the package maintainer and let the maintainer know that you need this functionality (in the case you actually need it).")
  }
  
  if (!missing(pedigree_must_be_alive) && (length(pedigree_must_be_alive) != 1L || !is.logical(pedigree_must_be_alive))) {
    stop("pedigree_must_be_alive must be TRUE or FALSE")
  }
  
  if (!missing(pedigree_birth_year_range) && (length(pedigree_birth_year_range) != 2L || !is.numeric(pedigree_birth_year_range))) {
    stop("if pedigree_birth_year_range is provided, it must be a numeric vector of length 2, (xmin, xmax), such that xmin <= birth_year <= xmax")
  }
  
  ##### ->   
  use_pedigree_criteria <- FALSE
  use_birth_year <- FALSE
  
  if ((!missing(pedigree_must_be_alive) && pedigree_must_be_alive == TRUE)) {
    use_pedigree_criteria <- TRUE
  }
  
  if (!missing(pedigree_birth_year_range)) {
    use_pedigree_criteria <- TRUE
    use_birth_year <- TRUE  
  }
  
  if (use_pedigree_criteria && missing(pedigree_birth_year_range)) {
    pedigree_birth_year_range <- c(0L, 0L)
  }
  ##### <- 
  
  
  radii <- sort(unique(radii), decreasing = TRUE)
  
  tmp_d1 <- radius_search_count(population, search_tree_unrelated, pids, radii)
  if (nrow(tmp_d1) == 1L && ncol(tmp_d1) == 1L && is.na(tmp_d1[1L, 1L])) {
    stop("Function was aborted") # FIXME: error class instead?
  }
  
  tmp_inds <- lapply(pids, get_individual, population = population)
  tmp_peds <- lapply(tmp_inds, get_pedigree_from_individual)
  
  #tmp_kdtrees <- lapply(tmp_peds, build_kdtree_from_pedigree, max_leaf_size = max_leaf_size_pedigree)
  tmp_kdtrees <- if (use_pedigree_criteria) {
    lapply(tmp_peds, function(ped) {
      ped_pids_criteria <- get_pids_in_pedigree_criteria(ped, 
        must_be_alive = pedigree_must_be_alive, 
        use_birth_year = use_birth_year,
        birth_year_min = pedigree_birth_year_range[1L], 
        birth_year_max = pedigree_birth_year_range[2L])
      
      print(ped_pids_criteria)
      
      build_kdtree_from_pids(pids = ped_pids_criteria, population = population, max_leaf_size = max_leaf_size_pedigree)
    })
  } else {
    lapply(tmp_peds, build_kdtree_from_pedigree, max_leaf_size = max_leaf_size_pedigree)
  }
  
  tmp_d2 <- lapply(seq_along(tmp_kdtrees), function(i) analyse_meioses_search_tree(population, tmp_kdtrees[[i]], pids[i], radii))
  tmp_d2 <- do.call(rbind, tmp_d2)

  tmp_d1_tbl <- data_frame(pid = tmp_d1[, 1L], radius_id = tmp_d1[, 2L], meiosis_d = tmp_d1[, 3L], n = tmp_d1[, 4L])
  tmp_d2_tbl <- data_frame(pid = tmp_d2[, 1L], radius_id = tmp_d2[, 2L], meiosis_d = tmp_d2[, 3L], n = tmp_d2[, 4L]) 
  
  tmp_d1_corrected <- tmp_d1_tbl %>% 
    left_join(tmp_d2_tbl %>% 
                 group_by(pid, radius_id) %>% 
                 summarise(n = sum(n)), 
               c("pid", "radius_id"))
  tmp_d1_corrected[is.na(tmp_d1_corrected)] <- 0
  
  tmp_d1_corrected <- tmp_d1_corrected %>%
    mutate(n = n.x - n.y) %>% 
    select(-n.x, -n.y)
  
  tmp_d <- rbind(tmp_d1_corrected, tmp_d2_tbl)
  
  d <- tibble::tibble(pid = tmp_d$pid, Radius = radii[tmp_d$radius_id], Meioses = tmp_d$meiosis_d, count = tmp_d$n)
  d <- d %>% arrange(pid, Radius, Meioses)

  d$Type <- "Related"
  d$Type[d$Meioses == -1L] <- "Unrelated"
  d$Type[d$Meioses != -1L] <- "Related"
  d$Meioses[d$Meioses == -1L] <- NA
  
  return(d)
}



#' Get distribution of number of meioses within a given radius. 
#' The Euclidian distance is used based on ETRS89 coordinates due to speed. 
#' Hence, be sure that this is an acceptable error in comparison to e.g. great-circle distance).
#'
#' For each individual (pid) in pids, a search is performed for all individuals in the population within a certain radius.
#' All the individuals within the radius are classified as either related or unrelated
#' (depending on whether they belong to the same pedigree or not). If they do belong to the same pedigree,
#' the number of meioses are found.
#'
#' Pedigree information is linked to the individuals in the population.
#'
#' @param population A collection of individuals (\code{popr_population})
#' @param search_tree_unrelated A search kd-tree when counting unrelated individuals (typically based on \code{population}), created with e.g. \link{\code{build_kdtree_from_population}} or \link{\code{build_kdtree_from_pids}}
#' @param max_leaf_size_pedigree A parameter for building kd-trees for the pedigrees that he individuals with pid in \code{pids} belongs to
#' @param pids Vector of pids to get meioses_distributions for (must be individuals in \code{population})
#' @param radii Vector of radii (radiuses) to search within [km]
#'
#' @return A tibble with meioses distribution information
#'
#' @export
meioses_distribution <- function(population, search_tree_unrelated, max_leaf_size_pedigree = 10, pids, radii) {
  if (!is(population, "popr_population")) stop("population must be a popr_population object")
  
  if (!is(search_tree_unrelated, "popr_search_data_structure")) stop("search_tree_unrelated must be a popr_search_data_structure object")
  
  radii <- sort(unique(radii), decreasing = TRUE)
  
  tmp_d1 <- radius_search_count(population, search_tree_unrelated, pids, radii)
  if (nrow(tmp_d1) == 1L && ncol(tmp_d1) == 1L && is.na(tmp_d1[1L, 1L])) {
    stop("Function was aborted") # FIXME: error class instead?
  }
  
  tmp_inds <- lapply(pids, get_individual, population = population)
  tmp_peds <- lapply(tmp_inds, get_pedigree_from_individual)
  tmp_igraphs <- lapply(tmp_peds, pedigree_as_igraph)
  tmp_kdtrees <- lapply(tmp_peds, build_kdtree_from_pedigree, max_leaf_size = max_leaf_size_pedigree)


  tmp_d2 <- lapply(seq_along(pids), function(pid_i) {
    ind <- tmp_inds[[pid_i]]
    g <- tmp_igraphs[[pid_i]]
    tree <- tmp_kdtrees[[pid_i]]
    pid_from <- as.character(pids[[pid_i]])
        
    ds <- do.call(rbind, (lapply(seq_along(radii), function(radius_i) {
      pids_match <- get_pids_within_radius(ind, tree, radii[radius_i])
      
      ds <- unlist(lapply(pids_match, function(pid_dest) {
        #igraph::shortest_paths(graph = g, from = pid_from, to = as.character(pid_dest), algorithm = "unweighted", mode = "all", output = "vpath")$vpath
        igraph::distances(graph = g, v = pid_from, to = as.character(pid_dest), algorithm = "unweighted", mode = "all")
      }))
      #return(ds)      
      
      ds <- table(ds)
    
      cbind(pids[[pid_i]], c(radius_i), as.integer(names(ds)), c(ds))
    })))
    
    ds
  })
  tmp_d2 <- do.call(rbind, tmp_d2)
  #tmp_d2
  
  tmp_d1_tbl <- data_frame(pid = tmp_d1[, 1L], radius_id = tmp_d1[, 2L], meiosis_d = tmp_d1[, 3L], n = tmp_d1[, 4L])
  tmp_d2_tbl <- data_frame(pid = tmp_d2[, 1L], radius_id = tmp_d2[, 2L], meiosis_d = tmp_d2[, 3L], n = tmp_d2[, 4L]) 
  
  tmp_d1_corrected <- tmp_d1_tbl %>% 
    inner_join(tmp_d2_tbl %>% 
                 group_by(pid, radius_id) %>% 
                 summarise(n = sum(n)), 
               c("pid", "radius_id")) %>% 
    mutate(n = n.x - n.y) %>% 
    select(-n.x, -n.y)
  
  tmp_d <- rbind(tmp_d1_corrected, tmp_d2_tbl)
  
  d <- tibble::tibble(pid = tmp_d$pid, Radius = radii[tmp_d$radius_id], Meioses = tmp_d$meiosis_d, count = tmp_d$n)
  d <- d %>% arrange(pid, Radius, Meioses)

  d$Type <- "Related"
  d$Type[d$Meioses == -1L] <- "Unrelated"
  d$Type[d$Meioses != -1L] <- "Related"
  d$Meioses[d$Meioses == -1L] <- NA
  
  return(d)
}





