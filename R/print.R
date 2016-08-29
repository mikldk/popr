#' @export
print.popr_population_abort <-
  function(x, ...) {
    if (!is(x, "popr_population_abort")) stop("x must be a popr_population_abort object")
    cat("Operation was cancelled, hence the assembly was not finished.\n")    
    return(invisible(NULL))
  }

#' @export
print.popr_population <-
  function(x, ...) {
    if (!is(x, "popr_population")) stop("x must be a popr_population object")
    
    cat("Population with ", formatC(pop_size(x), big.mark = ","), " individuals\n", sep = "")
    
    return(invisible(NULL))
  }

#' @export
print.popr_pedigreelist <-
  function(x, ...) {
    if (!is(x, "popr_pedigreelist")) stop("x must be a popr_pedigreelist object")
    
    sizes <- unlist(lapply(1L:pedigrees_count(x), function(i) pedigree_size(x[[i]])))
    sizes_str <- ""
    
    max_print <- 6L
    
    if (length(sizes) > 0L) {
      if (length(sizes) <= max_print) {
        sizes_str <- paste0(" (of size ", paste0(sizes, collapse = ", "), ")")
      } else {
        sizes_str <- paste0(" (of size ", paste0(sizes[1L:max_print], collapse = ", "), ", ...)")
      }
    }
    
    cat("List of ", formatC(pedigrees_count(x), big.mark = ","), " pedigrees", sizes_str, "\n", sep = "")
    
    return(invisible(NULL))
  }

  
#' @export
`[[.popr_pedigreelist` <- function(x, ...) {
  i <- ..1
  if (length(i) != 1L || !is.integer(i) || i[1L] <= 0L || i > pedigrees_count(x)) {
    stop("Wrong pedigree selected or invalid element selection criteria")
  }
  
  p <- get_pedigree(x, i - 1L) # -1 to go to 0-based indexing
  return(p)
}

#' @export
`[[.popr_population` <- function(x, ...) {
  pid <- ..1
  if (length(pid) != 1L || !is.integer(pid)) {
    stop("Wrong individual selected or invalid element selection criteria")
  }
  
  p <- get_individual(x, pid)
  return(p)
}

#' @export
print.popr_pedigree <-
  function(x, ...) {
    if (!is(x, "popr_pedigree")) stop("x must be a popr_pedigree object")
    
    print_pedigree(x)
    
    return(invisible(NULL))
  }

#' @export
pedigree_as_igraph <-
  function(x, ...) {
    if (!is(x, "popr_pedigree")) stop("x must be a popr_pedigree object")
    
    ginfo <- get_pedigree_as_graph(x)
    g <- igraph::graph_from_data_frame(ginfo$edgelist, directed = TRUE, vertices = ginfo$nodes)
    
    #co <- igraph::layout_nicely(g, dim = 2)
    co <- igraph::layout_as_tree(g, mode = "out")
    attr(g, "layout") <- co

    return(g)
  }

#' @export  
plot.popr_pedigree <-
  function(x, ...) {
    if (!is(x, "popr_pedigree")) stop("x must be a popr_pedigree object")
    
    g <- pedigree_as_igraph(x)    
    igraph::plot.igraph(g)
    
    return(invisible(NULL))
    #eturn(g)
  }
  
  
