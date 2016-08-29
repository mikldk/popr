// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "popr_types.hpp"
#include <Rcpp.h>

using namespace Rcpp;

// wipe_pedigrees
void wipe_pedigrees(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees);
RcppExport SEXP popr_wipe_pedigrees(SEXP pedigreesSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::XPtr< std::vector<Pedigree*> > >::type pedigrees(pedigreesSEXP);
    wipe_pedigrees(pedigrees);
    return R_NilValue;
END_RCPP
}
// build_pedigrees
Rcpp::XPtr< std::vector<Pedigree*> > build_pedigrees(Rcpp::XPtr<Population> population, bool progress);
RcppExport SEXP popr_build_pedigrees(SEXP populationSEXP, SEXP progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Population> >::type population(populationSEXP);
    Rcpp::traits::input_parameter< bool >::type progress(progressSEXP);
    __result = Rcpp::wrap(build_pedigrees(population, progress));
    return __result;
END_RCPP
}
// build_kdtree_from_population
Rcpp::XPtr<SearchDataStructure> build_kdtree_from_population(Rcpp::XPtr<Population> population, int max_leaf_size);
RcppExport SEXP popr_build_kdtree_from_population(SEXP populationSEXP, SEXP max_leaf_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Population> >::type population(populationSEXP);
    Rcpp::traits::input_parameter< int >::type max_leaf_size(max_leaf_sizeSEXP);
    __result = Rcpp::wrap(build_kdtree_from_population(population, max_leaf_size));
    return __result;
END_RCPP
}
// build_kdtree_from_pedigree
Rcpp::XPtr<SearchDataStructure> build_kdtree_from_pedigree(Rcpp::XPtr<Pedigree> pedigree, int max_leaf_size);
RcppExport SEXP popr_build_kdtree_from_pedigree(SEXP pedigreeSEXP, SEXP max_leaf_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Pedigree> >::type pedigree(pedigreeSEXP);
    Rcpp::traits::input_parameter< int >::type max_leaf_size(max_leaf_sizeSEXP);
    __result = Rcpp::wrap(build_kdtree_from_pedigree(pedigree, max_leaf_size));
    return __result;
END_RCPP
}
// build_kdtree_from_pids
Rcpp::XPtr<SearchDataStructure> build_kdtree_from_pids(IntegerVector pids, Rcpp::XPtr<Population> population, int max_leaf_size);
RcppExport SEXP popr_build_kdtree_from_pids(SEXP pidsSEXP, SEXP populationSEXP, SEXP max_leaf_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerVector >::type pids(pidsSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<Population> >::type population(populationSEXP);
    Rcpp::traits::input_parameter< int >::type max_leaf_size(max_leaf_sizeSEXP);
    __result = Rcpp::wrap(build_kdtree_from_pids(pids, population, max_leaf_size));
    return __result;
END_RCPP
}
// analyse_meioses_search_tree
IntegerMatrix analyse_meioses_search_tree(Rcpp::XPtr<Population> population, Rcpp::XPtr<SearchDataStructure> search_tree, IntegerVector pids, NumericVector radii, bool report_progress);
RcppExport SEXP popr_analyse_meioses_search_tree(SEXP populationSEXP, SEXP search_treeSEXP, SEXP pidsSEXP, SEXP radiiSEXP, SEXP report_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Population> >::type population(populationSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<SearchDataStructure> >::type search_tree(search_treeSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type pids(pidsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type radii(radiiSEXP);
    Rcpp::traits::input_parameter< bool >::type report_progress(report_progressSEXP);
    __result = Rcpp::wrap(analyse_meioses_search_tree(population, search_tree, pids, radii, report_progress));
    return __result;
END_RCPP
}
// get_pids_within_radius
IntegerVector get_pids_within_radius(Rcpp::XPtr<Individual> ind, Rcpp::XPtr<SearchDataStructure> search_tree, double radius);
RcppExport SEXP popr_get_pids_within_radius(SEXP indSEXP, SEXP search_treeSEXP, SEXP radiusSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Individual> >::type ind(indSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<SearchDataStructure> >::type search_tree(search_treeSEXP);
    Rcpp::traits::input_parameter< double >::type radius(radiusSEXP);
    __result = Rcpp::wrap(get_pids_within_radius(ind, search_tree, radius));
    return __result;
END_RCPP
}
// radius_search_count
IntegerMatrix radius_search_count(Rcpp::XPtr<Population> population, Rcpp::XPtr<SearchDataStructure> search_tree, IntegerVector pids, NumericVector radii, bool report_progress);
RcppExport SEXP popr_radius_search_count(SEXP populationSEXP, SEXP search_treeSEXP, SEXP pidsSEXP, SEXP radiiSEXP, SEXP report_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Population> >::type population(populationSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<SearchDataStructure> >::type search_tree(search_treeSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type pids(pidsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type radii(radiiSEXP);
    Rcpp::traits::input_parameter< bool >::type report_progress(report_progressSEXP);
    __result = Rcpp::wrap(radius_search_count(population, search_tree, pids, radii, report_progress));
    return __result;
END_RCPP
}
// wipe_population
void wipe_population(Rcpp::XPtr<Population> population);
RcppExport SEXP popr_wipe_population(SEXP populationSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Population> >::type population(populationSEXP);
    wipe_population(population);
    return R_NilValue;
END_RCPP
}
// load_individuals
Rcpp::XPtr<Population> load_individuals(IntegerVector pid, LogicalVector is_male, IntegerVector pid_mom, IntegerVector pid_dad, NumericVector etrs89e, NumericVector etrs89n, bool progress, bool error_on_gender_mismatch);
RcppExport SEXP popr_load_individuals(SEXP pidSEXP, SEXP is_maleSEXP, SEXP pid_momSEXP, SEXP pid_dadSEXP, SEXP etrs89eSEXP, SEXP etrs89nSEXP, SEXP progressSEXP, SEXP error_on_gender_mismatchSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerVector >::type pid(pidSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type is_male(is_maleSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type pid_mom(pid_momSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type pid_dad(pid_dadSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type etrs89e(etrs89eSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type etrs89n(etrs89nSEXP);
    Rcpp::traits::input_parameter< bool >::type progress(progressSEXP);
    Rcpp::traits::input_parameter< bool >::type error_on_gender_mismatch(error_on_gender_mismatchSEXP);
    __result = Rcpp::wrap(load_individuals(pid, is_male, pid_mom, pid_dad, etrs89e, etrs89n, progress, error_on_gender_mismatch));
    return __result;
END_RCPP
}
// pop_size
int pop_size(Rcpp::XPtr<Population> population);
RcppExport SEXP popr_pop_size(SEXP populationSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Population> >::type population(populationSEXP);
    __result = Rcpp::wrap(pop_size(population));
    return __result;
END_RCPP
}
// get_number_of_children
IntegerMatrix get_number_of_children(Rcpp::XPtr<Population> population, bool progress);
RcppExport SEXP popr_get_number_of_children(SEXP populationSEXP, SEXP progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Population> >::type population(populationSEXP);
    Rcpp::traits::input_parameter< bool >::type progress(progressSEXP);
    __result = Rcpp::wrap(get_number_of_children(population, progress));
    return __result;
END_RCPP
}
// popr_test
void popr_test();
RcppExport SEXP popr_popr_test() {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    popr_test();
    return R_NilValue;
END_RCPP
}
// pedigrees_count
int pedigrees_count(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees);
RcppExport SEXP popr_pedigrees_count(SEXP pedigreesSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::XPtr< std::vector<Pedigree*> > >::type pedigrees(pedigreesSEXP);
    __result = Rcpp::wrap(pedigrees_count(pedigrees));
    return __result;
END_RCPP
}
// pedigree_size
int pedigree_size(Rcpp::XPtr<Pedigree> ped);
RcppExport SEXP popr_pedigree_size(SEXP pedSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Pedigree> >::type ped(pedSEXP);
    __result = Rcpp::wrap(pedigree_size(ped));
    return __result;
END_RCPP
}
// pedigrees_table
std::unordered_map<int, int> pedigrees_table(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees);
RcppExport SEXP popr_pedigrees_table(SEXP pedigreesSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::XPtr< std::vector<Pedigree*> > >::type pedigrees(pedigreesSEXP);
    __result = Rcpp::wrap(pedigrees_table(pedigrees));
    return __result;
END_RCPP
}
// get_pedigree
Rcpp::XPtr<Pedigree> get_pedigree(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees, int index);
RcppExport SEXP popr_get_pedigree(SEXP pedigreesSEXP, SEXP indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::XPtr< std::vector<Pedigree*> > >::type pedigrees(pedigreesSEXP);
    Rcpp::traits::input_parameter< int >::type index(indexSEXP);
    __result = Rcpp::wrap(get_pedigree(pedigrees, index));
    return __result;
END_RCPP
}
// get_individual
Rcpp::XPtr<Individual> get_individual(Rcpp::XPtr<Population> population, int pid);
RcppExport SEXP popr_get_individual(SEXP populationSEXP, SEXP pidSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Population> >::type population(populationSEXP);
    Rcpp::traits::input_parameter< int >::type pid(pidSEXP);
    __result = Rcpp::wrap(get_individual(population, pid));
    return __result;
END_RCPP
}
// print_individual
void print_individual(Rcpp::XPtr<Individual> individual);
RcppExport SEXP popr_print_individual(SEXP individualSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Individual> >::type individual(individualSEXP);
    print_individual(individual);
    return R_NilValue;
END_RCPP
}
// print_pedigree
void print_pedigree(Rcpp::XPtr<Pedigree> ped);
RcppExport SEXP popr_print_pedigree(SEXP pedSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Pedigree> >::type ped(pedSEXP);
    print_pedigree(ped);
    return R_NilValue;
END_RCPP
}
// get_pids_in_pedigree
IntegerVector get_pids_in_pedigree(Rcpp::XPtr<Pedigree> ped);
RcppExport SEXP popr_get_pids_in_pedigree(SEXP pedSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Pedigree> >::type ped(pedSEXP);
    __result = Rcpp::wrap(get_pids_in_pedigree(ped));
    return __result;
END_RCPP
}
// get_pedigree_edgelist
CharacterMatrix get_pedigree_edgelist(Rcpp::XPtr<Pedigree> ped);
RcppExport SEXP popr_get_pedigree_edgelist(SEXP pedSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Pedigree> >::type ped(pedSEXP);
    __result = Rcpp::wrap(get_pedigree_edgelist(ped));
    return __result;
END_RCPP
}
// get_pedigree_as_graph
List get_pedigree_as_graph(Rcpp::XPtr<Pedigree> ped);
RcppExport SEXP popr_get_pedigree_as_graph(SEXP pedSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Pedigree> >::type ped(pedSEXP);
    __result = Rcpp::wrap(get_pedigree_as_graph(ped));
    return __result;
END_RCPP
}
// meiosis_dist_tree
int meiosis_dist_tree(Rcpp::XPtr<Individual> src, Rcpp::XPtr<Individual> dest);
RcppExport SEXP popr_meiosis_dist_tree(SEXP srcSEXP, SEXP destSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Individual> >::type src(srcSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<Individual> >::type dest(destSEXP);
    __result = Rcpp::wrap(meiosis_dist_tree(src, dest));
    return __result;
END_RCPP
}
// meiosis_dist_tree_matrix
IntegerMatrix meiosis_dist_tree_matrix(Rcpp::XPtr<Pedigree> ped);
RcppExport SEXP popr_meiosis_dist_tree_matrix(SEXP pedSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Pedigree> >::type ped(pedSEXP);
    __result = Rcpp::wrap(meiosis_dist_tree_matrix(ped));
    return __result;
END_RCPP
}
// get_pedigree_from_individual
Rcpp::XPtr<Pedigree> get_pedigree_from_individual(Rcpp::XPtr<Individual> individual);
RcppExport SEXP popr_get_pedigree_from_individual(SEXP individualSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Individual> >::type individual(individualSEXP);
    __result = Rcpp::wrap(get_pedigree_from_individual(individual));
    return __result;
END_RCPP
}
// get_pedigree_id_from_pid
IntegerVector get_pedigree_id_from_pid(Rcpp::XPtr<Population> population, IntegerVector pids);
RcppExport SEXP popr_get_pedigree_id_from_pid(SEXP populationSEXP, SEXP pidsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Population> >::type population(populationSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type pids(pidsSEXP);
    __result = Rcpp::wrap(get_pedigree_id_from_pid(population, pids));
    return __result;
END_RCPP
}
