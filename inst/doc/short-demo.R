## ------------------------------------------------------------------------
data("sample.ped", package = "kinship2")
head(sample.ped)

## ------------------------------------------------------------------------
library(dplyr)
library(tibble)

## ------------------------------------------------------------------------
library(popr)

## ------------------------------------------------------------------------
pop <- with(sample.ped, load_individuals(pid = id, 
                                         is_male = (sex == 1), 
                                         pid_mom = mother, 
                                         pid_dad = father,
                                         progress = FALSE))

## ------------------------------------------------------------------------
children_count <- as_data_frame(get_number_of_children(pop, progress = FALSE))
children_count %>% filter(boys_n > 0 | girls_n > 0)

## ------------------------------------------------------------------------
peds <- build_pedigrees(pop, progress = FALSE)
peds

## ------------------------------------------------------------------------
lapply(1L:pedigrees_count(peds), function(i) pedigree_size(peds[[i]]))

## ------------------------------------------------------------------------
ped_ids <- get_pedigree_id_from_pid(population = pop, pids = sample.ped$id)
table(ped_ids, sample.ped$ped)

## ------------------------------------------------------------------------
peds[[3L]]
pid <- get_pids_in_pedigree(peds[[3L]])
sample.ped %>% filter(id == pid | father == pid | mother == pid)

## ------------------------------------------------------------------------
peds[[2L]]
pids <- get_pids_in_pedigree(peds[[2L]])
pids
sample.ped %>% filter(id %in% pids | mother %in% pids | father %in% pids)
children_count %>% filter(pid %in% pids)

## ------------------------------------------------------------------------
plot(peds[[2L]])
g <- pedigree_as_igraph(peds[[2L]])
g
plot(g)
get_pedigree_as_graph(peds[[2L]])

## ------------------------------------------------------------------------
set.seed(1) # For reproducibility
sample.ped.coords <- as_data_frame(cbind(sample.ped, tibble(etrs89e = runif(nrow(sample.ped), 250, 650), etrs89n = runif(nrow(sample.ped), 650, 750))))
sample.ped.coords

## ------------------------------------------------------------------------
pop_coord <- with(sample.ped.coords, load_individuals(pid = id, 
                                                      is_male = (sex == 1), 
                                                      pid_mom = mother, 
                                                      pid_dad = father,
                                                      etrs89e = etrs89e,
                                                      etrs89n = etrs89n,
                                                      progress = FALSE))

## ------------------------------------------------------------------------
peds_coord <- build_pedigrees(pop_coord, progress = FALSE)

## ------------------------------------------------------------------------
kdtree <- build_kdtree_from_population(pop_coord, max_leaf_size = 10)

## ------------------------------------------------------------------------
d <- meioses_distribution(population = pop_coord, search_tree_unrelated = kdtree,
                          max_leaf_size_pedigree = 10, 
                          pids = sample.ped.coords$id[1L], 
                          radii = c(100, 50, 100))

## ------------------------------------------------------------------------
d

